"""Nox sessions for testing and linting."""
import nox

# Default sessions to run when just calling `nox`
nox.options.sessions = ["lint", "typecheck", "tests"]


@nox.session(python=["3.7", "3.8", "3.9", "3.10", "3.11"])
def tests(session):
    """Run the test suite."""
    session.install(".[test]")
    session.run("pytest", "tests/", "-v")


@nox.session
def lint(session):
    """Run ruff linting checks."""
    session.install("ruff")
    session.run("ruff", "check", "rdkit_to_params/")


@nox.session
def lint_fix(session):
    """Run ruff linting checks and auto-fix issues."""
    session.install("ruff")
    session.run("ruff", "check", "--fix", "rdkit_to_params/")


@nox.session
def format_check(session):
    """Check code formatting with ruff."""
    session.install("ruff")
    session.run("ruff", "format", "--check", "rdkit_to_params/")


@nox.session
def format(session):
    """Format code with ruff (auto-fix formatting issues)."""
    session.install("ruff")
    session.run("ruff", "format", "rdkit_to_params/")



@nox.session
def typecheck(session):
    """Run mypy type checking."""
    session.install("mypy")
    session.install(".")
    session.run("mypy", "rdkit_to_params/")

@nox.session
def complexity_check(session):
    """Check code complexity metrics using radon.

    Reports cyclomatic complexity (CC) and maintainability index (MI).
    """
    session.install("radon")
    # cyclomatic complexity (CC)
    session.run(
        "radon",
        "cc",
        "rdkit_to_params/",
        "--min",
        "D",
        "--show-complexity",
        "--average",
        *session.posargs,
    )
    # maintainability index (MI)
    session.run(
        "radon",
        "mi",
        "rdkit_to_params/",
        "--min",
        "B",
        "--show",
    )

@nox.session
def all_checks(session):
    """Run all linting and type checking."""
    session.install(".[test]")
    session.run("ruff", "check", "rdkit_to_params/")
    session.run("ruff", "format", "--check", "rdkit_to_params/")
    session.run("mypy", "rdkit_to_params/")