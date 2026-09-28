# Paper repository workflow

This is the canonical manuscript repository for the CASF benchmark paper:
https://github.com/YerevaNN/casf-benchmark-paper

- Edit the paper in `jcim-overleaf/`; its main document is `acs-template.tex`.
- Read `jcim-overleaf/AGENTS.md` and its linked writing/evidence instructions.
- The user has requested that future paper changes be committed and pushed here
  using the authorized MenuaB account. The Overleaf-connected branch is `main`.
- Fetch and inspect incoming changes before editing/pushing; preserve edits from
  Overleaf and other authors. Do not force-push or replace repository history.
- Validate manuscript citations and cross-references in the source for text edits.
  Do not compile or update local previews unless the user explicitly requests it;
  the user reviews the compiled manuscript in Overleaf.
  Figures must remain traceable to the archived data and exact model identities.
- The benchmark code and original databases live in the separate
  `YerevaNN/casf-benchmark` repository. Paths starting with `docs/` or `src/` in
  evidence notes refer to that analysis repository unless stated otherwise.
- The old benchmark `manuscript/` directory is a historical working copy. Do not
  treat it as the source of truth or overwrite this paper repository from it.
- Keep TeX build intermediates out of commits. Leave the historical reading
  preview under `jcim-overleaf/preview/` unchanged during routine edits.
