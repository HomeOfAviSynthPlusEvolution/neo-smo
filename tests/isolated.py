"""Run a test command with inherited Windows crash UI suppression and a timeout."""
import argparse
import ctypes
import os
import subprocess
import sys


def suppress_crash_ui():
    if os.name == "nt":
        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel32.GetErrorMode.restype = ctypes.c_uint
        kernel32.SetErrorMode.argtypes = [ctypes.c_uint]
        # Inherited by children unless CREATE_DEFAULT_ERROR_MODE is requested.
        kernel32.SetErrorMode(kernel32.GetErrorMode() | 0x0001 | 0x0002 | 0x8000)
        kernel32.WerSetFlags.argtypes = [ctypes.c_uint]
        kernel32.WerSetFlags(0x0004)  # WER_FAULT_REPORTING_NO_UI


def run(command, timeout=60):
    suppress_crash_ui()
    try:
        result = subprocess.run(command, timeout=timeout, capture_output=True, text=True,
                                creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
    except subprocess.TimeoutExpired:
        print(f"[TIMEOUT] {command!r}", file=sys.stderr, flush=True)
        return 124
    print(result.stdout, end="", flush=True)
    print(result.stderr, end="", file=sys.stderr, flush=True)
    if result.returncode:
        print(f"[CHILD FAILED] exit=0x{result.returncode & 0xffffffff:08X} {command!r}",
              file=sys.stderr, flush=True)
    return result.returncode


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timeout", type=float, default=60)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.timeout <= 0 or not args.command:
        parser.error("a positive timeout and a command are required")
    sys.exit(1 if run(args.command, timeout=args.timeout) else 0)
