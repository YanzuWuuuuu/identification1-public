"""Keep the report current and audit sources after the finite MC batch exits.

This is a companion to one already-authorized local batch, not a scheduler.
It exits when that batch exits, including abnormal exits.
"""
import argparse
import ctypes
from ctypes import wintypes
from datetime import datetime
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]


def write(path, data):
    temporary = path.with_suffix(path.suffix+".tmp")
    temporary.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf8")
    temporary.replace(path)


def audit_sources():
    original = json.loads((ROOT/"input_hashes.json").read_text(encoding="utf-8-sig"))
    changed, missing, errors = [], [], []
    checked = 0
    for item in original:
        path = Path(item["path"])
        if not path.exists():
            missing.append(str(path)); continue
        try:
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            checked += 1
            if digest.lower() != item["sha256"].lower():
                changed.append(str(path))
        except Exception as exc:
            errors.append(dict(path=str(path), error=repr(exc)))
    result = dict(checked=checked, total=len(original), changed=changed,
                  missing=missing, errors=errors,
                  note="Changes, if any, are reported; this check alone cannot attribute concurrent external edits.")
    write(ROOT/"analysis"/"source_hash_verification.json", result)
    return result


def counts():
    files = [p for p in (ROOT/"mc_qr_final/raw").glob("*.json")
             if not p.name.endswith(".started.json")]
    ok, failures, invalid = 0, 0, []
    for p in files:
        try:
            record = json.loads(p.read_text(encoding="utf8"))
            if record.get("success"): ok += 1
            else: failures += 1
        except Exception as exc:
            invalid.append(dict(file=str(p), exception=repr(exc)))
    return dict(saved=len(files), successful=ok, failed=failures,
                planned=360, remaining=360-len(files), invalid_files=invalid)


def main(pid):
    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
    kernel.OpenProcess.restype = wintypes.HANDLE
    kernel.WaitForSingleObject.argtypes = [wintypes.HANDLE, wintypes.DWORD]
    kernel.WaitForSingleObject.restype = wintypes.DWORD
    kernel.GetExitCodeProcess.argtypes = [wintypes.HANDLE, ctypes.POINTER(wintypes.DWORD)]
    kernel.CloseHandle.argtypes = [wintypes.HANDLE]
    handle = kernel.OpenProcess(0x00100000 | 0x1000, False, pid)
    last_report = -100
    final = None
    try:
        while True:
            running = bool(handle and kernel.WaitForSingleObject(handle, 0) == 258)
            state = counts()
            state.update(pid=pid, running=running, updated=datetime.now().isoformat(),
                         batch_complete_marker=(ROOT/"mc_qr_final/batch_complete.json").exists())
            if not running:
                code = wintypes.DWORD()
                if handle and kernel.GetExitCodeProcess(handle, ctypes.byref(code)):
                    state["exit_code"] = int(code.value)
                state["complete"] = bool(state["saved"] == 360 and not state["invalid_files"]
                                         and state["batch_complete_marker"])
            else:
                state["complete"] = False
            write(ROOT/"RUN_STATUS.json", state)
            if state["saved"]-last_report >= 8 or not running:
                finished = subprocess.run([sys.executable, "-B", str(ROOT/"code/summarize.py")],
                                          cwd=ROOT, capture_output=True, text=True,
                                          encoding="utf8", errors="replace")
                with (ROOT/"logs/report_updates.log").open("a", encoding="utf8") as stream:
                    stream.write(f"{state['updated']} saved={state['saved']} report_exit={finished.returncode}\n")
                    stream.write(finished.stdout+finished.stderr)
                last_report = state["saved"]
            if not running:
                state["source_audit"] = audit_sources()
                write(ROOT/"RUN_STATUS.json", state)
                final = state
                break
            time.sleep(30)
    finally:
        if handle: kernel.CloseHandle(handle)
    print(json.dumps(final, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--pid", required=True, type=int)
    main(ap.parse_args().pid)
