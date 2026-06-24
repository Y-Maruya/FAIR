#!/usr/bin/env python3

import os
import re
import json
import glob
import urllib.request
import urllib.error

import ROOT


EOS_HOST = "eosproject-f.cern.ch"
EOS_DIR = "/eos/project/f/faser-upgrade/AHCAL/FAIR"

CONFIG_FILE = "/eos/user/y/ymaruya/FASER/AHCAL/FAIR/config/MIP_first_run.yaml"
BAD_RUN_OUTPUT = "bad_mipcalib_runs.txt"

TREE_NAME = "events"
RUNINFO_URL = "https://faser-runinfo.app.cern.ch/cgibin/getRunInfo.py?runno={runno}"


def list_eos_root_files(eos_dir):
    return glob.glob(f"{eos_dir}/*_MIPCalib.root")


def extract_runno(filename):
    base = os.path.basename(filename)
    m = re.match(r"(\d+)_MIPCalib\.root$", base)
    if not m:
        return None
    return int(m.group(1))


def get_runinfo_events(runno):
    """
    Get expected number of physics events from FASER runinfo.

    Expected format:
      {
        "runinfo": {
          "eventCounts": {
            "Events_sent_Physics": 128,
            ...
          }
        }
      }
    """
    url = RUNINFO_URL.format(runno=runno)

    try:
        with urllib.request.urlopen(url, timeout=20) as response:
            text = response.read().decode("utf-8", errors="replace")
    except urllib.error.HTTPError as e:
        print(f"RUNINFO_ERROR run={runno}: HTTP {e.code}")
        return None
    except urllib.error.URLError as e:
        print(f"RUNINFO_ERROR run={runno}: {e}")
        return None
    except Exception as e:
        print(f"RUNINFO_ERROR run={runno}: {e}")
        return None

    try:
        data = json.loads(text)
    except Exception as e:
        print(f"RUNINFO_ERROR run={runno}: failed to parse JSON: {e}")
        return None

    try:
        return int(data["runinfo"]["eventCounts"]["Events_sent_Physics"])
    except KeyError as e:
        print(f"RUNINFO_ERROR run={runno}: missing key {e}")
        return None
    except Exception as e:
        print(f"RUNINFO_ERROR run={runno}: {e}")
        return None


def check_root_file(filename, tree_name=TREE_NAME):
    result = {
        "filename": filename,
        "runno": extract_runno(filename),
        "open_ok": False,
        "is_zombie": True,
        "recovered": False,
        "has_tree": False,
        "entries": None,
        "first_entry_ok": False,
        "last_entry_ok": False,
        "error": None,
    }

    try:
        f = ROOT.TFile.Open(filename, "READ")
    except OSError as e:
        result["error"] = f"OSError: {e}"
        return result
    except Exception as e:
        result["error"] = f"Exception: {e}"
        return result

    if not f:
        result["error"] = "TFile.Open returned null"
        return result

    try:
        result["is_zombie"] = bool(f.IsZombie())
        result["open_ok"] = bool(f.IsOpen()) and not result["is_zombie"]

        if not result["open_ok"]:
            result["error"] = "file is zombie or not open"
            f.Close()
            return result

        result["recovered"] = bool(f.TestBit(ROOT.TFile.kRecovered))

        t = f.Get(tree_name)
        if not t:
            result["error"] = f"TTree '{tree_name}' not found"
            f.Close()
            return result

        result["has_tree"] = True
        result["entries"] = int(t.GetEntries())

        if result["entries"] <= 0:
            result["error"] = "zero entries"
            f.Close()
            return result

        result["first_entry_ok"] = bool(t.GetEntry(0) > 0)
        result["last_entry_ok"] = bool(t.GetEntry(result["entries"] - 1) > 0)

        if not result["first_entry_ok"]:
            result["error"] = "failed to read first entry"
        elif not result["last_entry_ok"]:
            result["error"] = "failed to read last entry"

        f.Close()
        return result

    except Exception as e:
        result["error"] = f"Exception while checking file: {e}"
        try:
            f.Close()
        except Exception:
            pass
        return result


def is_bad_result(res, expected_events):

    if expected_events == 0:
        return False  # If expected events is zero, we don't consider it bad

    if not res["open_ok"]:
        return True

    if not res["has_tree"]:
        return True

    if res["entries"] is None:
        return True

    if not res["first_entry_ok"] or not res["last_entry_ok"]:
        return True

    if expected_events is None:
        return True

    if res["entries"] < expected_events:
        return True
    
    if res["recovered"]:
        return True

    return False


def status_string(res, expected_events):
    if not res["open_ok"]:
        return "OPEN_FAILED"

    if not res["has_tree"]:
        return "NO_TREE"

    if res["entries"] is None:
        return "NO_ENTRIES"

    if not res["first_entry_ok"] or not res["last_entry_ok"]:
        return "READ_FAILED"

    if expected_events is None:
        return "NO_RUNINFO"

    if res["entries"] < expected_events:
        return "EVENT_MISMATCH"

    if res["recovered"]:
        return "RECOVERED"

    return "OK"


if __name__ == "__main__":
    root_files = list_eos_root_files(EOS_DIR)

    root_files = [
        f for f in root_files
        if extract_runno(f) is not None
    ]

    root_files.sort(key=lambda f: extract_runno(f))

    print(f"Found {len(root_files)} files")

    bad_runs = []

    for filename in root_files:
        runno = extract_runno(filename)

        res = check_root_file(filename)
        expected_events = get_runinfo_events(runno)

        status = status_string(res, expected_events)
        bad = is_bad_result(res, expected_events)

        print(
            f"{status:14s} "
            f"run={runno} "
            f"entries={res['entries']} "
            f"expected={expected_events} "
            f"recovered={res['recovered']} "
            f"{filename}"
        )

        if res["error"]:
            print(f"  error: {res['error']}")

        if bad:
            bad_runs.append(runno)

    with open(BAD_RUN_OUTPUT, "w") as fout:
        for runno in bad_runs:
            fout.write(f"{CONFIG_FILE} {runno}\n")

    print()
    print(f"Bad runs: {len(bad_runs)}")
    print(f"Wrote: {BAD_RUN_OUTPUT}")