#!/usr/bin/env python3

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


SCRIPT_PATH = Path(__file__).resolve().parent.parent / "scripts/build_runinfo_snapshot.py"
SPEC = importlib.util.spec_from_file_location("build_runinfo_snapshot", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)

RUN_LIST = [
    {"runnumber": 10, "type": "AHCAL"},
    {"runnumber": 11, "type": "Physics"},
    {"runnumber": 12, "type": "AHCAL"},
]


def detail(run_number):
    return {"runnumber": run_number, "type": "AHCAL", "configuration": {}}


class RunInfoSnapshotBuilderTest(unittest.TestCase):
    def test_fetches_only_ahcal_details_and_reuses_them(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "snapshot.json"
            calls = []

            def first_fetch(url, params=None, **_kwargs):
                calls.append((url, params))
                return RUN_LIST if params is None else detail(params["runno"])

            with patch.object(MODULE, "fetch_json", side_effect=first_fetch):
                MODULE.build_snapshot(output, timeout=1, attempts=1, retry_delay=0)

            snapshot = json.loads(output.read_text())
            self.assertEqual(snapshot["run_types"], {"10": "AHCAL", "11": "Physics", "12": "AHCAL"})
            self.assertEqual(set(snapshot["details"]), {"10", "12"})
            self.assertEqual(
                [params["runno"] for _url, params in calls if params is not None],
                [10, 12],
            )

            with patch.object(MODULE, "fetch_json", return_value=RUN_LIST) as fetch:
                MODULE.build_snapshot(output, timeout=1, attempts=1, retry_delay=0)
            fetch.assert_called_once()

    def test_resumes_after_detail_fetch_failure(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "snapshot.json"

            def failing_fetch(_url, params=None, **_kwargs):
                if params is None:
                    return RUN_LIST
                if params["runno"] == 12:
                    raise RuntimeError("temporary failure")
                return detail(params["runno"])

            with patch.object(MODULE, "fetch_json", side_effect=failing_fetch):
                with self.assertRaisesRegex(RuntimeError, "temporary failure"):
                    MODULE.build_snapshot(output, timeout=1, attempts=1, retry_delay=0)

            partial = json.loads(output.read_text())
            self.assertEqual(set(partial["details"]), {"10"})

            calls = []

            def resumed_fetch(_url, params=None, **_kwargs):
                calls.append(params)
                return RUN_LIST if params is None else detail(params["runno"])

            with patch.object(MODULE, "fetch_json", side_effect=resumed_fetch):
                MODULE.build_snapshot(output, timeout=1, attempts=1, retry_delay=0)
            self.assertEqual([params["runno"] for params in calls if params is not None], [12])


if __name__ == "__main__":
    unittest.main()
