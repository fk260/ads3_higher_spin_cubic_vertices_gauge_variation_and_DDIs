#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
tmp_log="$(mktemp)"
status=0
wolframscript -file ddis_vz/run_ddis_vz.wl 2>&1 | tee "$tmp_log" || status=$?
mkdir -p ddi_vz_outputs
mv "$tmp_log" ddi_vz_outputs/terminal_capture.txt
echo "Output folder: ddi_vz_outputs"
exit "$status"
