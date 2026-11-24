"""Baseline agent: uniform random search within parameter bounds.

Serves as the reference budget benchmark; the Bayesian/agentic loop should
beat it at equal simulation count.
"""

import random


class RandomSearchAgent:
    def __init__(self, env):
        self.env = env
        self.best_reward = -float("inf")
        self.best_theta = None
        self.history = []

    def act(self):
        return {
            name: random.uniform(spec["bounds"][0], spec["bounds"][1])
            for name, spec in self.env.params.items()
        }

    def run(self, budget):
        for _ in range(budget):
            theta = self.act()
            try:
                reward, info = self.env.step(theta)
            except Exception as exc:  # bounds/param errors: reject, don't crash
                self.history.append({"theta": theta, "reward": None, "error": str(exc)})
                continue
            self.history.append({"theta": theta, "reward": reward, **info})
            if reward > self.best_reward:
                self.best_reward, self.best_theta = reward, theta
        return self.best_reward, self.best_theta, self.history
