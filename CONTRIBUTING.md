# Contributing to DockFlow-Automator

Thank you for considering a contribution!  This project is scientific
software: **reproducibility and scientific honesty outrank features**
(see the README's *Scientific limitations* section and the frozen
feature roadmap).

## Development setup

```bash
git clone https://github.com/oualid-razzok/DockFlow-Automator
cd DockFlow-Automator
python -m pip install -e ".[test,prep,engine,viz]"   # + [gui] for Qt work
python -m pip install ./bindings                      # optional C++ accelerator
pytest -m "not network and not gui and not scientific" # software tests
ruff check .
```

## Process

1. **Open an issue first** for anything larger than a typo fix,
   especially new features - the feature set is frozen behind the
   scientific-validation track (README roadmap).
2. Branch from `main`, commit with
   **conventional-commits prefixes**:
   `fix:` `feat:` `docs:` `test:` `bench:` `refactor:` `chore:`
   (e.g. `fix(prep): fall through when OBElementTable is missing`).
3. Every behaviour change needs a test that fails without it; every
   *scientific* claim change needs a scientific-validation update
   (`benchmarks/`, `docs/preparation_validation.md`, the reproducibility
   workflow).
4. `ruff check .` must pass; the full offline test suite must pass on
   all three OS (CI checks this for you).
5. **No force-pushes to `main`** - the commit history is part of the
   provenance story reviewers rely on.  Rebase your own feature branch
   freely.
6. Releases are tagged with **signed annotated tags**
   (`git tag -s v0.2.0 -m "..."`) and every tag has a CHANGELOG entry.

## Scientific contributions specifically

- New preparation/reduction behaviours must record their decisions in
  `manifest.json` (see `receptor.decisions`, `ligands[*].prep`).
- Changes that can shift scores (typing, charges, defaults) must be
  reflected in `docs/preparation_validation.md` and, ideally, exercised
  by the reproducibility workflow.
- Never describe geometric contacts as bonds, or score proxies as
  experimental quantities (see `docs/interaction_criteria.md`).
