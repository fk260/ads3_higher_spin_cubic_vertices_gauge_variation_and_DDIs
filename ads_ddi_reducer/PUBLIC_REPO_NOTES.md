# Public repository notes

This folder is meant to be committed as source code, not as generated data.

Recommended policy:

- Commit `package/`, `ddis_v1/`, `ddis_vz/`, `terms_and_stage_tests/`, helper shell scripts, and documentation.
- Do not commit generated output folders by default.
- Keep paper/reference comparisons outside this folder unless they are added as a separate optional verification layer.
- Use neutral output folder names (`ddi_v1_outputs/`, `ddi_vz_outputs/`, `term_stage_outputs/`).
- Keep terminal captures as local run artifacts; they are ignored by `.gitignore`.

A useful pre-commit smoke test is:

```bash
bash run_terms_and_stage_tests.sh
bash run_ddis_v1.sh
bash run_ddis_vz.sh
```

Then inspect:

```bash
git status --ignored
```

to confirm generated outputs are ignored.
