import json
import os
from pathlib import Path


RESULTS_PATH = Path(f"./tests/results/py_values.json")


def pytest_configure(config):
    if not os.path.exists(RESULTS_PATH.parent):
        os.makedirs(RESULTS_PATH.parent)
    with open(RESULTS_PATH, "w") as f:
        json.dump({}, f)


def pytest_addoption(parser):
    parser.addoption(
        "--log-results",
        action="store_true",
        default=False,
        help="Log results to a json file for comparing to Fathom",
    )


def pytest_runtest_teardown(item, nextitem):
    log_results = item.config.getoption("--log-results")

    if log_results and hasattr(item, "return_value"):
        with open(RESULTS_PATH, "r+") as f:
            data = json.load(f)
            data[item.name] = item.return_value

        with open(RESULTS_PATH, "w") as f:
            json.dump(data, f, indent=None)
