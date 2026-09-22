"""Compare the R LOO elpd and target with the original Python code.

Reads the data and R results written by export_for_python.R and recomputes
them with ProblemRun in simulated/general_setup.py.
Run from the repo root after the R export:
    pipeline/.venv/bin/python pipeline/checks/compare_python.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, "simulated")
from general_setup import ProblemRun  # noqa: E402

TMP = "pipeline/checks/tmp"


def load(name):
    return np.loadtxt(os.path.join(TMP, name))


with open(os.path.join(TMP, "tau2.txt")) as f:
    tau2_str = f.read().strip()
tau2 = None if tau2_str == "none" else float(tau2_str)

y = load("y.txt")
y_test = load("y_test.txt")
n_obs = 16
n_trial = y.size // n_obs
n_test = y_test.size // n_obs
y = y.reshape(n_trial, n_obs)
y_test = y_test.reshape(n_test, n_obs)
X = load("X.txt").reshape(n_trial, n_obs, -1)
X_test = load("X_test.txt").reshape(n_test, n_obs, -1)

pr = ProblemRun(n_obs=n_obs, n_obs_max=n_obs, tau2=tau2)
py = {
    "loo_a": pr.calc_loo_ti(y, X[:, :, :-1]),
    "loo_b": pr.calc_loo_ti(y, X),
    "target_a": pr.calc_elpd_tl(
        y, X[:, :, :-1], y_test, X_test[:, :, :-1]).mean(axis=1),
    "target_b": pr.calc_elpd_tl(y, X, y_test, X_test).mean(axis=1),
}

print("tau2 = {}".format(tau2_str))
for name, value in py.items():
    r = load(name + ".txt").reshape(value.shape)
    print("{:9s} max |R - Python| = {:.1e}".format(
        name, np.max(np.abs(r - value))))
