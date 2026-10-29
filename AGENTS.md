# AGENTS.md

Guidance for AI coding agents (GitHub Copilot, Cursor, Claude Code, Codex, and others) contributing to repositories of the UMC Utrecht Genetics (UMCUGenetics) team. Read this file before making changes. If instructions in this file conflict with repo-specific instructions in a subdirectory, the more specific instruction wins.

## Context: who we are

- We are the bioinformatics team of the Genetics department of a Dutch university medical center (UMC Utrecht).
- We develop and maintain software and bioinformatics pipelines, many of which are **used in diagnostic processes** (ISO 15189 accredited laboratory, in-house IVDs under IVDR article 5.5).
- Our repositories are published on the [UMCUGenetics GitHub](https://github.com/UMCUGenetics).
- Development follows our SOP "SP.2 Software development", the coding guidelines CG.1–CG.7, and the Dutch field standard for bioinformatics (the "[Veldnorm](https://vkgl-kwaliteit.github.io/BioinformaticaVeldnorm/)"), which is based on IEC 62304 Software Safety Class A.

Because our software can affect patient diagnostics, correctness, traceability, and reviewability matter more than speed of delivery. When in doubt, be conservative and ask the developer.

## Non-negotiable rules

1. **Never commit directly to `main`/`master` or `develop`.** All work happens in feature/hotfix branches and is merged via a pull request with at least one human approval. Never merge your own PR, never force-push, never disable or work around branch protection.
2. **Never weaken or delete tests to make code pass.** If a test fails, fix the code or — if the intended behavior genuinely changed — update the test and explain why in the PR description.
3. **Never include real patient or sample data** in code, tests, fixtures, comments, commit messages, or logs. Use synthetic/anonymous test data.
4. **Do not invent requirements.** Software requirements are documented and validated per release. If a request is ambiguous, ask instead of assuming.
5. **Run linters and tests before presenting work as done.** Code that does not pass linting and the test suite is not finished.
6. **Keep changes scoped.** Do not refactor, reformat, or "improve" code outside the task at hand; mixed refactor + feature PRs are hard to review and break traceability.
7. **Flag regulatory follow-up instead of deciding it yourself** (see "Regulated environment" below).

## Regulated environment: what you must flag

Every functionally meaningful change in diagnostic software can trigger documentation obligations. As an agent you cannot perform validation or sign off documents — **list these items for the developer whenever they may apply**:

- **Software Documentation & Software Validation updates**: for tools used in diagnostics, requirements, risk analysis, and the deployment procedure are maintained in separate controlled documents. Functional changes, new features, and changed behavior usually require updating them.
- **Risk analysis**: the Veldnorm requires risk analysis to be revised for each release. New dependencies, changed inputs/outputs, or changed failure behavior are risk-relevant.
- **Validation**: each release (including hotfixes) requires a recorded validation, performed by a human (preferably an end user, independent of the developer). You cannot do this.
- **SOUP registration** (Software of Unknown Provenance): every new third-party package/dependency must be registered with name, reference/URL, intended purpose, whether it is used according to its intended purpose, whether separate validation is required, and whether it is version controlled.
- **Audit logging**: in applications that follow our layered architecture (see below), functionally meaningful actions (create/update/delete of records, status changes, starting analyses) must be audit-logged.
- **Version numbering**: releases follow Semantic Versioning (`major.minor.patch`): breaking changes → major, new features → minor, bugfixes/small changes → patch. Suggest the appropriate bump but let the developer decide.

## Git workflow

We use a Git Flow-based model:

- `main` (or `master`): released code only; each release gets a git tag / GitHub release.
- `develop`: changes ready for release.
- `feature/<name>`: created from `develop`, merged back into `develop` via reviewed PR.
- `release/<version>`: created from `main`, for release validation; fixes during validation are committed here; merges into `main` and `develop`.
- `hotfix/<version>`: created from `main` for production bugs; merges into `main` and `develop`.

Write clear commit messages and PR descriptions (what + why). PRs must pass automated checks (unit tests, linters) before merging.

## General coding style (all languages)

- **Line length: max 127 characters** (GitHub diff limit). This overrides stricter community defaults in our tooling.
- Write **clear, modular code**: short functions with one responsibility, a clear interface, no hidden side effects or state.
- **Naming**: meaningful, consistent, pronounceable names; unambiguous booleans (`is_valid`, `has_samples`); short names only for short-lived loop variables; constants named by meaning (`DAYS_PER_WEEK`, not `SEVEN`); use domain terminology.
- **Comments**: describe *what* and *why*, not *how*; never restate the code; keep comments in sync with the code; prefer rewriting confusing code over commenting it.
- **Docstrings**: every function/method/class gets a docstring describing purpose, inputs (use parameter names), outputs, and side effects, in the language's preferred docstring style.
- **Whitespace/format**: one statement per line; blank lines to group logical sections; spaces around binary operators and after commas; no spaces inside `obj.method()` chains or before argument-list parentheses; indent with spaces, never tabs.
- Don't duplicate code: if two blocks do the same thing, extract one function.

## Testing (all languages)

Automated tests are required, in a separate test folder that mirrors the application structure. Every test:

- is fast, independent, deterministic, and tests one thing;
- verifies behavior, not implementation details;
- has a descriptive name (e.g. `test_start_analysis_without_validated_metadata_raises_error`);
- covers calculations, parsing, and file name/path handling especially.

Guidelines:

- Aim for full coverage of our own code; imported/external modules are not our responsibility.
- Parameterize variations of the same behavior; write separate tests for different behaviors.
- Use temporary directories for file-based tests; choose fixture scope deliberately.
- Simple test data is created on the fly; large or structured test data goes into external files (one file per behavior, in a test-data/fixtures folder).
- Group related tests (e.g. by class); use dependency injection and mock external systems (network, databases) where necessary — with caution.
- Tests are part of the commit; CI (GitHub Actions) runs them on PRs and on `main`/`develop`.
- There is no minimum coverage threshold yet, but do not leave new logic untested.

## Architecture & design patterns

**When to apply:** our layered architecture is the standard for new standalone applications/services (see `DxSampleDatabase` as reference implementation) and for repositories that already follow it. **When not to apply:** do not retrofit this pattern onto existing/legacy code that does not use it, and do not introduce it for small scripts or one-off tools — that would be overkill. If the repo already follows the pattern, keep new code consistent with it.

Layered structure:

```text
API endpoint / Router / Controller / CLI entry
    ↓
Service layer (business logic)
    ↓
CRUD / Repository / DAO layer (database operations)
    ↓
Database
```

Rules per layer:

- **Routers/controllers** receive requests, parse parameters, call the correct service function, translate functional errors to HTTP status codes. No business logic here.
- **Services** contain business rules: existence checks, required metadata, valid status transitions, preventing invalid/duplicate actions, coordinating repositories, audit logging, calling external processes. No raw SQL/ORM query details here.
- **CRUD/repositories** are simple: get/save/update/delete records, centralized queries, transactions. No process rules here (nothing decides whether an action is *allowed*).
- **Models** define entities; simple data validation (format, required values, types) is allowed, workflow rules are not.
- **Infrastructure** (`db/`, engines, sessions, migrations) is technical configuration only.

Anti-patterns to avoid: business logic in routers, business logic in CRUD functions, the same query duplicated across the codebase, missing audit logging for functionally meaningful actions.

## Language-specific rules

### Python

- **Style**: PEP 8, except line length (127). Tooling: [Ruff](https://docs.astral.sh/ruff/) for linting and formatting — `ruff check --config "line-length=127"` and `ruff format --config "line-length=127"`.
- **Tests**: [pytest](https://docs.pytest.org/). Test files/functions start with `test_`, snake_case for functions/variables, UpperCamelCase for classes. Custom markers must be documented in `pytest.ini`; use `skip`/`xfail` only if strictly necessary. Use `tmp_path`/`tmp_path_factory` for temporary files, [pytest-datair](https://pypi.org/project/pytest-datadir/) for external test data files. `pytest --durations=5` helps find slow tests.
- **Environment**: new projects use [UV](https://docs.astral.sh/uv/) (`uv init`, `uv add <package>`, `uv run`). Legacy projects use `venv` + pip with a `requirements.txt` (always committed). Do not mix methods within a repo.
- **Preferred packages**: FastAPI (APIs), Click or Typer (CLIs), Requests or HTTPX (HTTP), SQLAlchemy or SQLModel (SQL), Flask (full-stack web), VCFPy (VCF parsing), Genologics (Clarity API client). Candidates under discussion (confirm with the team before introducing): Polars, Loguru, Seaborn, samshee.
- **Starting point**: [python_template](https://github.com/UMCUGenetics/python_template) shows the intended repo layout, tests, and formatting setup.

### R

- **Style**: tidyverse style guide, two deliberate deviations: always end functions with an explicit `return()`, and line length 127. Configure these in `.lintr`:

```text
linters: linters_with_defaults(
    line_length_linter(127),
    return_linter = NULL
  )
encoding: "UTF-8"
```

- **Lint/format**: [lintr](https://lintr.r-lib.org/) (`lintr::use_lintr(type = "tidyverse")` on first use) and [styler](https://styler.r-lib.org/).
- **Docstrings**: [roxygen2](https://roxygen2.r-lib.org/) comments above each function; `roxygen2::roxygenise()` generates package docs.
- **Tests**: [testthat](https://testthat.r-lib.org/), in `tests/testthat/` with fixtures in `tests/fixtures/`. Prefer `.txt` test data files (visible in GitHub diffs). Run with `testthat::test_dir("tests/testthat")`.
- **Structure/packages**: [devtools](https://devtools.r-lib.org/) `create_package()`; see the [R packages book](https://r-pkgs.org/). Unit test and lint workflows run via GitHub Actions ([r-lib/actions](https://github.com/r-lib/actions) examples).
- **Environment**: [renv](https://rstudio.github.io/renv/) or Docker (e.g. the DIMS pipeline uses versioned Docker images).

### Nextflow

- Follow the [nf-core guidelines](https://nf-co.re/docs/guidelines/components/overview) and the team's [nf-core template](https://github.com/UMCUGenetics/nf-core_template).

### JavaScript / HTML / CSS

- No team-specific guideline document exists yet. Follow community best practices, keep the code style consistent with the existing repository, and format with Prettier.

## Tooling notes

- Recommended editors: VS Code or JetBrains IDEs (RStudio for R). Relevant tooling the agent can rely on being present: Ruff and Python extensions (VS Code), lintr/styler addins (RStudio), markdownlint, EditorConfig.
- CI: GitHub Actions run available tests and linters on PRs and on `main`/`develop`. If a repo has no test/lint workflow and the change adds substantial functionality, suggest adding one.

## Documentation duties

- Keep the README accurate: goal of the project and how to get started. Update it when functionality changes.
- Docstrings and inline documentation are part of "done", not optional extras.

## When in doubt

- Ask the developer rather than guessing, especially for: intended behavior in diagnostics workflows, anything touching patient data handling, release/version decisions, validation scope, and new third-party dependencies.
- Explicitly report at the end of a task: what was changed, which tests were added/updated, which lint/CI checks were run, and which regulatory follow-ups (documentation, risk analysis, validation, SOUP registration) may be needed.