"""Build timer regression checks; no model inputs or Fortran build required."""
import concurrent.futures
import csv
import os
import shutil
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

TIMER = Path(__file__).resolve().parents[1] / "cmake" / "build_timer.py"


class BuildTimerTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="ecosim timing ")
        self.addCleanup(self.temp.cleanup)
        self.log = Path(self.temp.name) / "timings.log"
        self.env = dict(os.environ, ECOSIM_BUILD_RUN_ID="test")
        self.env.pop("ECOSIM_BUILD_TIMING_LOG", None)

    def run_timer(self, code, *extra, phase="compile", env=None):
        return subprocess.run(
            [sys.executable, str(TIMER), "run", "--log", str(self.log),
             "--phase", phase, "--", sys.executable, "-c", code, *extra],
            env=env or self.env, capture_output=True, text=True)

    def test_build_script_stops_before_configure_when_cmake_is_missing(self):
        fixture = Path(self.temp.name) / "missing cmake"
        fixture.mkdir()
        empty_path = fixture / "empty-bin"
        empty_path.mkdir()
        result = subprocess.run(
            [shutil.which("bash"), str(TIMER.parents[1] / "build_EcoSIM.sh")],
            cwd=fixture, env=dict(self.env, PATH=str(empty_path)),
            capture_output=True, text=True)
        self.assertEqual(result.returncode, 127)
        self.assertIn("CMake was not found on PATH", result.stderr)
        self.assertFalse((fixture / "build").exists())

    def records(self, path=None):
        with open(path or self.log, newline="") as stream:
            return list(csv.DictReader(stream, delimiter="\t"))

    def test_output_arguments_and_failed_exit_are_preserved(self):
        result = self.run_timer(
            "import sys; print(sys.argv[1]); print('compiler error', file=sys.stderr); sys.exit(17)",
            "argument with spaces; $literal")
        self.assertEqual(result.returncode, 17)
        self.assertEqual(result.stdout, "argument with spaces; $literal\n")
        self.assertEqual(result.stderr, "compiler error\n")
        row, = self.records()
        self.assertEqual(row["exit_code"], "17")
        self.assertGreater(float(row["elapsed_seconds"]), 0)

    def test_parallel_first_writers_keep_every_record_and_one_header(self):
        def launch(i):
            result = self.run_timer(f"import time; time.sleep(0.05); print({i})")
            self.assertEqual(result.returncode, 0, result.stderr)
        with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
            list(pool.map(launch, range(16)))
        rows = self.records()
        self.assertEqual(len(rows), 16)
        self.assertEqual(self.log.read_text().count("started_utc\trun_id\t"), 1)
        self.assertTrue(all(row["exit_code"] == "0" and row["run_id"] == "test" for row in rows))

    def test_nested_launchers_share_environment_log_and_report_filter(self):
        selected = Path(self.temp.name) / "selected.log"
        env = dict(self.env, ECOSIM_BUILD_TIMING_LOG=str(selected))
        self.assertEqual(self.run_timer("pass", phase="configure", env=env).returncode, 0)
        self.assertEqual(self.run_timer("pass", phase="build", env=env).returncode, 0)
        self.assertEqual(self.run_timer("pass", phase="install", env=dict(env, ECOSIM_BUILD_RUN_ID="other")).returncode, 0)
        self.assertFalse(self.log.exists())
        result = subprocess.run([sys.executable, str(TIMER), "report", "--log", str(selected),
                                 "--run-id", "test"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Recorded commands: 2; failed: 0", result.stdout)
        self.assertIn("configure", result.stdout)
        self.assertNotIn("install", result.stdout)
        self.assertIn("their sum is not build wall time", result.stdout)

    def test_build_script_writes_summary_after_failed_configure(self):
        source = TIMER.parents[1]
        fixture = Path(self.temp.name) / "script fixture"
        (fixture / "cmake").mkdir(parents=True)
        (fixture / "tools").mkdir()
        shutil.copyfile(source / "build_EcoSIM.sh", fixture / "build_EcoSIM.sh")
        shutil.copyfile(TIMER, fixture / "cmake" / "build_timer.py")
        fake_cmake = fixture / "tools" / "cmake"
        fake_cmake.write_text("#!/bin/sh\nexit 23\n")
        fake_cmake.chmod(0o755)
        env = dict(self.env, PATH=str(fixture / "tools") + os.pathsep + self.env["PATH"])
        env.pop("ATS_ECOSIM", None)
        result = subprocess.run(["bash", "build_EcoSIM.sh", "--jobs=2"], cwd=fixture,
                                env=env, capture_output=True, text=True)
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        summary, = list((fixture / "build").glob("*/build-timings-summary.txt"))
        self.assertIn("Recorded commands: 1; failed: 1", summary.read_text())
        log, = list((fixture / "build").glob("*/build-timings.log"))
        row, = self.records(log)
        self.assertEqual((row["phase"], row["exit_code"]), ("configure", "23"))

    def test_signal_exit_status_is_not_reported_as_success(self):
        result = self.run_timer("import os,signal; os.kill(os.getpid(), signal.SIGTERM)")
        self.assertEqual(result.returncode, 143)
        self.assertEqual(self.records()[0]["exit_code"], "143")

    def test_missing_command_is_recorded(self):
        result = subprocess.run([sys.executable, str(TIMER), "run", "--log", str(self.log),
                                 "--phase", "compile", "--", "/nonexistent/ecosim-compiler"],
                                env=self.env, capture_output=True, text=True)
        self.assertEqual(result.returncode, 127)
        self.assertEqual(self.records()[0]["exit_code"], "127")


if __name__ == "__main__":
    unittest.main()
