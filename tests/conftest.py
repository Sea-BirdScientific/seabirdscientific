import json
import os
from datetime import date, datetime
from pathlib import Path


now = datetime.now().isoformat(timespec='seconds').replace(':', '-')
RESULTS_PATH = Path(f"./tests/results/results_{now}.json")


def pytest_configure(config):
    if not os.path.exists(RESULTS_PATH.parent):
        os.makedirs(RESULTS_PATH.parent)
    with open(RESULTS_PATH, "w") as f:
        json.dump({}, f)


def pytest_runtest_teardown(item, nextitem):
    if hasattr(item, "return_value"):
        with open(RESULTS_PATH, "r+") as f:
            data = json.load(f)
            data[item.name] = item.return_value

        with open(RESULTS_PATH, "w") as f:
            json.dump(data, f, indent=None)
