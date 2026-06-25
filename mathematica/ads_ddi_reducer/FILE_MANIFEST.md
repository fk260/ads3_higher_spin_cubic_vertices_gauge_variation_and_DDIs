# File manifest for `ads_ddi_reducer`

Commit these source files:

```text
ads_ddi_reducer/README.md
ads_ddi_reducer/.gitignore
ads_ddi_reducer/FILE_MANIFEST.md
ads_ddi_reducer/package/AdSDDIReducer.wl
ads_ddi_reducer/ddis_v1/run_ddis_v1.wl
ads_ddi_reducer/ddis_v1/run_ddis_v1.sh
ads_ddi_reducer/ddis_vz/run_ddis_vz.wl
ads_ddi_reducer/ddis_vz/run_ddis_vz.sh
ads_ddi_reducer/terms_and_stage_tests/test_terms.wl
ads_ddi_reducer/terms_and_stage_tests/test_terms.nb
ads_ddi_reducer/terms_and_stage_tests/run_terms_and_stage_tests.sh
ads_ddi_reducer/helpers/*.sh
ads_ddi_reducer/run_*.sh
```

Do not commit generated folders such as:

```text
ddi_v1_outputs/
ddi_vz_outputs/
term_stage_outputs/
*_terminal_capture.txt
*_capture.txt
```

This git-ready folder deliberately excludes the paper-comparison/reference-formula layer.
