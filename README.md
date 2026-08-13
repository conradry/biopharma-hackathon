# biopharma-hackathon

## Setup

Requires [uv](https://docs.astral.sh/uv/).

```bash
uv sync                    # create .venv and install deps (incl. dev group)
uv run pre-commit install  # enable git hooks
```

## Common commands

```bash
uv run pytest              # run tests
uv run ruff check --fix .  # lint
uv run ruff format .       # format
uv add <pkg>               # add a runtime dependency
uv add --dev <pkg>         # add a dev dependency
uv run pre-commit run --all-files
```

Package code lives in `src/biopharma_hackathon/`, tests in `tests/`.
