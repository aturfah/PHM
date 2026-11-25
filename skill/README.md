# PHM agent skill

An agent skill that teaches an AI coding agent (Claude Code, Claude.ai, or any agent that reads `SKILL.md` skills) to use the PHM package and PHM Explorer: density estimation, the PHM ensemble, `PHMv2()` and the original `PHM()`, cutting and plotting the hierarchy, the command-line ensemble pipeline (serial, one machine, Slurm), and exporting fits to the interactive Explorer page.

## Layout

```
phm/                          the skill (this folder is what an agent loads)
  SKILL.md                    overview, entry points, core recipe, plots, Explorer, scale, pitfalls
  references/api.md           every exported function: signatures, defaults, return values, known bugs
  references/workflows.md     worked recipes (GMM, clusterers, weights, resolution, criteria, ensemble, localized, big data, Seurat, PHM())
  references/ensemble.md      the ensemble in R and from the command line; config.R fields; reruns; settings
  references/explorer.md      PHM Explorer export, one-command viewer, reading the display, troubleshooting
  scripts/setup_phm.R         checks for and installs PHM and its dependencies
  scripts/explorer/           copies of PHM_Explorer/r/{export_phm_json.R, view_in_explorer.R} and src/{template.html, build.py}
tests/smoke_test.R            checks the skill's claims and recipes against the installed PHM
evals/evals.json              example prompts for evaluating the skill
```

## Install the skill

Claude Code, for every project (personal skill):

```bash
ln -s "$(pwd)/skill/phm" ~/.claude/skills/phm      # run from the PHM repo root; or cp -r skill/phm ~/.claude/skills/
```

For one project, link or copy `skill/phm` into that project's `.claude/skills/`. For Claude.ai or other agents that take uploads, zip the `phm` folder (so the zip holds `phm/SKILL.md`) and upload it. Agents that do not support skills can be pointed at `skill/phm/SKILL.md` as instructions; it tells them which reference file to read for each task.

The R package must be installed where the code will run:

```bash
Rscript skill/phm/scripts/setup_phm.R --check-only         # report only
Rscript skill/phm/scripts/setup_phm.R --explorer           # install what is missing (+ jsonlite), then a short end-to-end check
```

Other flags: `--optional` (ggpubr, ppclust), `--update` (reinstall PHM from GitHub), `--lib=DIR`, `--no-verify`.

## Keeping it in sync

- `scripts/explorer/` is a copy of the Explorer tool in `PHM_Explorer/`. After changing the Explorer, refresh it from the repo root:

  ```bash
  cp PHM_Explorer/r/export_phm_json.R PHM_Explorer/r/view_in_explorer.R skill/phm/scripts/explorer/r/
  cp PHM_Explorer/src/template.html PHM_Explorer/src/build.py skill/phm/scripts/explorer/src/
  ```

- After changing the package, run `Rscript skill/tests/smoke_test.R` (about two minutes). It checks defaults, return structures, quirks, every recipe marked "tested", the ensemble functions, and the Explorer export. A failing check means that part of the skill needs revisiting; checks labelled KNOWN BUG start failing when the bug is fixed upstream.
- `skill/` is excluded from the R package build (`.Rbuildignore`).

## Credits

Builds on the earlier `xqwen/phm-skill`. Cite PHM as: Turfah A, Wen X (2025). Multiresolution clustering of genomic data. bioRxiv. doi:10.1101/2025.06.13.659529. GPL-3.0, like the package.
