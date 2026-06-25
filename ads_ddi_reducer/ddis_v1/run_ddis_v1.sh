#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
tmp_log="$(mktemp)"
status=0
wolframscript -file ddis_v1/run_ddis_v1.wl 2>&1 | tee "$tmp_log" || status=$?
mkdir -p ddi_v1_outputs
mv "$tmp_log" ddi_v1_outputs/terminal_capture.txt
echo "Output folder: ddi_v1_outputs"
exit "$status"
