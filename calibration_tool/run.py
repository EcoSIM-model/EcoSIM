"""CLI entry: python -m calibration_tool.run --case cases/biocrust/config.json"""

import argparse
import json
import sys
from pathlib import Path

from calibration_tool.env import EcoSIMCalibEnv
from calibration_tool.agent.random_search import RandomSearchAgent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True)
    ap.add_argument("--budget", type=int, default=10)
    ap.add_argument("--agent", choices=["random"], default="random")
    args = ap.parse_args()

    env = EcoSIMCalibEnv(args.case)
    if args.agent == "random":
        agent = RandomSearchAgent(env)
        best_r, best_theta, history = agent.run(args.budget)

    print(f"best reward: {best_r}")
    print(f"best theta: {json.dumps(best_theta, indent=2)}")
    out = Path("calibration_tool/runs") / "search_history.json"
    with open(out, "w") as f:
        json.dump(history, f, indent=2, default=str)
    print(f"history written to {out}")


if __name__ == "__main__":
    sys.exit(main())
