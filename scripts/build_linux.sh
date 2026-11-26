#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
command -v ifx >/dev/null 2>&1 || { echo "Load the Intel oneAPI compiler and MKL environment first." >&2; exit 1; }
python3 tools/verify_sources.py
python3 tools/build.py "$@"
