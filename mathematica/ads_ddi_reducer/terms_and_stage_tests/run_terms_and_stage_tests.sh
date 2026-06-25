#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
tmp_log="$(mktemp)"
status=0
wolframscript -file terms_and_stage_tests/test_terms.wl 2>&1 | tee "$tmp_log" || status=$?
mkdir -p term_stage_outputs
mv "$tmp_log" term_stage_outputs/terminal_capture.txt
echo "Output folder: term_stage_outputs"
exit "$status"
