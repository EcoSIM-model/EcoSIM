"""Gym-style environment wrapping an EcoSIM case for calibration.

Each step() runs one full simulation in an isolated episode directory:
  reset/prepare -> clone pristine case dir + local copy of the PFT-parameter
                   NetCDF (never touch the shared input_data file)
  step(theta)   -> sanity bounds check -> ParamEditor mutation -> run executable
                   -> h0-target extraction -> weighted reward

All external tools are invoked as documented by their skills:
  - ecosim-agent/.agents/skills/ecosim-pftpar-editor  (ParamEditor.py)
  - ecosim-agent/.agents/skills/ecosim-h0-target-extractor/.../extract_h0_targets.py
"""

import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_TOOLS = REPO_ROOT / "ecosim-agent"
PARAM_EDITOR_DIR = PYTHON_TOOLS / "applications" / "notebooks"
H0_EXTRACTOR = (
    PYTHON_TOOLS
    / ".agents/skills/ecosim-h0-target-extractor/scripts/extract_h0_targets.py"
)

RESTART_GLOBS = (
    "*.ecosim.r.*",
    "*.ecosim.rh0.*",
    "*.ecosim.h0.*",
    "rpointer.esim",
    "fort.*",
)


class BoundsError(ValueError):
    pass


class EcoSIMCalibEnv:
    def __init__(self, config_path):
        with open(config_path) as f:
            self.cfg = json.load(f)
        case_dir = Path(config_path).resolve().parent
        self.source_dir = (case_dir / self.cfg["source_run_dir"]).resolve()
        self.exe = Path(self.cfg["ecosim_exe"]).resolve()
        self.python = Path(self.cfg["python"]).absolute()
        self.pft_code = self.cfg["pft_code"]
        self.pft_file_name = self.cfg["pft_file_in"]
        self.params = {p["name"]: p for p in self.cfg["parameters"]}
        self.targets = self.cfg["targets"]
        self.spinup_years = self.cfg.get("spinup_years")
        self.keep_dirs = self.cfg.get("keep_last_n", 3)
        self.runs_dir = REPO_ROOT / "calibration_tool" / "runs"
        self.runs_dir.mkdir(parents=True, exist_ok=True)
        self._episode = 0

    # ------------------------------------------------------------------
    def reset(self):
        """Prepare a fresh pristine episode directory. Returns its path."""
        self._episode += 1
        ep = self.runs_dir / f"ep_{self._episode:04d}"
        if ep.exists():
            shutil.rmtree(ep)
        ep.mkdir(parents=True)
        self._clone_case(ep)
        self._localize_pft_file(ep)
        self._absolutize_paths(ep)
        self._prune_old_episodes()
        return ep

    def step(self, theta, episode_dir=None):
        """Run one calibration step for parameter dict theta.

        Returns (reward, info). episode_dir=None triggers reset().
        """
        ep = Path(episode_dir) if episode_dir else self.reset()
        self._check_bounds(theta)

        t0 = time.time()
        self._write_params(ep, theta)
        run_ok = self._run_model(ep)
        if not run_ok:
            return -1e6, {"episode": str(ep), "failed": "simulation"}

        sim_targets = self._extract_targets(ep)
        reward = self._reward(sim_targets)
        info = {
            "episode": str(ep),
            "sim_targets": sim_targets,
            "wall_s": round(time.time() - t0, 1),
        }
        return reward, info

    # ------------------------------------------------------------------
    def _clone_case(self, ep):
        for item in self.source_dir.iterdir():
            dst = ep / item.name
            if item.is_dir():
                shutil.copytree(item, dst)
            else:
                shutil.copy2(item, dst)
        for pattern in RESTART_GLOBS:
            for stale in ep.glob(pattern):
                stale.unlink()

    def _localize_pft_file(self, ep):
        """Copy shared pftpar NetCDF into the episode and repoint the namelist."""
        src = self._find_pft_file(ep)
        dst = ep / src.name
        shutil.copy2(src, dst)
        for nml in ep.glob("*.namelist"):
            text = nml.read_text()
            text = re.sub(
                r"pft_file_in\s*=\s*'[^']*'",
                f"pft_file_in='{dst.name}'",
                text,
            )
            nml.write_text(text)

    def _find_pft_file(self, ep):
        """Resolve pft_file_in as written in the episode namelist.

        Paths are relative to the run dir (e.g. ../../../input_data/...), so
        they must be resolved against the *original* source dir, not the clone.
        """
        for nml in ep.glob("*.namelist"):
            m = re.search(r"pft_file_in\s*=\s*'([^']+)'", nml.read_text())
            if m:
                candidate = (self.source_dir / m.group(1)).resolve()
                if candidate.exists():
                    return candidate
        raise FileNotFoundError(
            f"pft_file_in not found from namelists in {ep} (source: {self.source_dir})"
        )

    def _absolutize_paths(self, ep):
        """Rewrite namelist input paths relative to the *original* run dir.

        Clones live elsewhere in the tree, so relative input paths (grid,
        climate, GHG, management files) must be resolved against source_dir
        and replaced with absolute paths.
        """
        for nml in ep.glob("*.namelist"):
            text = nml.read_text()

            def sub(m):
                token = m.group(2)
                if (ep / Path(token).name).exists():
                    return m.group(0)  # already-localized file (pftpar)
                candidate = (self.source_dir / token).resolve()
                if candidate.is_file():
                    return f"{m.group(1)}='{candidate}'"
                if candidate.is_dir():
                    return f"{m.group(1)}='{candidate}/'"
                return m.group(0)

            text = re.sub(r"(\w*(?:_in|prefix))\s*=\s*'([^']+)'", sub, text)
            nml.write_text(text)

    def _check_bounds(self, theta):
        for name, value in theta.items():
            if name not in self.params:
                raise BoundsError(f"unknown parameter: {name}")
            lo, hi = self.params[name]["bounds"]
            if not lo <= value <= hi:
                raise BoundsError(f"{name}={value} outside [{lo}, {hi}]")

    def _write_params(self, ep, theta):
        """Apply theta via ParamEditor, imported per the pftpar-editor skill."""
        sys.path.insert(0, str(PARAM_EDITOR_DIR))
        try:
            from scripts import ParamEditor  # type: ignore
        finally:
            sys.path.remove(str(PARAM_EDITOR_DIR))
        pftpar = ep / Path(self.pft_file_name).name
        editor = ParamEditor.ParEditor(pftparfile=str(pftpar))
        editor.PlantParamModify(self.pft_code, dict(theta))

    def _run_model(self, ep):
        # drivers/ecosim/ecosim.F90 reads the namelist path from argv[1]
        namelists = sorted(ep.glob("*.namelist"))
        if not namelists:
            raise FileNotFoundError(f"no namelist in {ep}")
        proc = subprocess.run(
            [str(self.exe), namelists[0].name],
            cwd=str(ep),
            capture_output=True,
            text=True,
            timeout=self.cfg.get("timeout_s", 3600),
        )
        if proc.returncode != 0:
            (ep / "run_stderr.log").write_text(proc.stderr[-20000:])
            return False
        return True

    def _h0_file(self, ep):
        h0s = sorted(ep.glob("*.ecosim.h0.*.nc"))
        if not h0s:
            raise FileNotFoundError(f"no h0 output in {ep}")
        return h0s[-1]

    def _extract_targets(self, ep):
        out = Path(ep) / "h0_targets.json"
        subprocess.run(
            [
                str(self.python),
                str(H0_EXTRACTOR),
                str(self._h0_file(ep)),
                "--format",
                "json",
                "--output",
                str(out),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        with open(out) as f:
            raw = json.load(f)
        return {
            k: v.get("summary", {}).get("typical_value")
            for k, v in raw.get("metrics", {}).items()
            if isinstance(v, dict)
        }

    def _reward(self, sim_targets):
        """Negative weighted normalized RMSE against observation targets."""
        total_w, sq = 0.0, 0.0
        for key, spec in self.targets.items():
            sim = sim_targets.get(key)
            obs, w = spec["observed"], spec["weight"]
            if sim is None or sim != sim:  # missing or NaN
                return -1e6
            sq += w * ((sim - obs) / abs(obs)) ** 2
            total_w += w
        return -((sq / total_w) ** 0.5)

    def _prune_old_episodes(self):
        eps = sorted(
            d for d in self.runs_dir.iterdir() if d.is_dir() and d.name.startswith("ep_")
        )
        for d in eps[: max(0, len(eps) - self.keep_dirs)]:
            shutil.rmtree(d)
