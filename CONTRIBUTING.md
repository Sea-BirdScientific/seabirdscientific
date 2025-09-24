# Contributing Guidelines

## Before Contributing

Thanks for your interest in contributing to the seabirdscientific repository. Before submitting any pull requests, ensure that you __have fully read and understand these guidelines__. If you're new to the process, take a look at the [Code Contribution Process Overview](#code-contribution-process-overview) section. If you've read through everything and still have any questions about contributing, please feel free to [submit a question](https://github.com/Sea-BirdScientific/seabirdscientific/issues/new?template=question.md).

## Contributing

### Contributor Agreement

We are grateful that you are considering a contribution. By being one of our contributors, you agree and confirm that:

- You created the submitted code. If it's based on someone else's work, give them due credit.
- Your work will be distributed under the [MIT license](./LICENSE) once your pull request is merged.

  - If your work is derived from code that is __already covered by the MIT License__, ask the originator if they would also like to be listed as a contributor to the toolkit before submitting your PR.
  - __If your work is derived from code that is not covered by the MIT License, get the originator's approval to release it under the MIT license. Also find out if they would also like to be listed as a contributor before submitting your PR.__

- Your submitted work is in keeping with our styles and standards below.

### Coding Style

We want your work to be readable by others; therefore, we ask you to comply with the following:

- Adhere to the [PEP 8](https://peps.python.org/pep-0008/) style guide for your code.
  - One diversion from the PEP 8 guide that we have chosen is to allow maximum line lengths of __99__ characters vs the recommended 79 characters.
- Take care when naming functions, classes, and variables. Help your reader by using __descriptive names__.
  - Avoid single letter variable names unless their life only spans a few lines.
  - Expand acronyms. Names like `gcd()` may not be obvious but `greatest_common_divisor()` is.
  - Follow the [PEP 8 Naming Conventions](https://pep8.org/#prescriptive-naming-conventions), including, but not limited to, the following:
    - Variable names and function names should be `snake_case`
    - Constants should be in `UPPER_CASE`
    - ClassNames should be `CamelCase`
- The use of [Python type hints](https://docs.python.org/3/library/typing.html) are encouraged for function parameters and return values.
- Use Python 3.9 or greater.

### Documenting

- Make use of docstrings based upon the [Sphinx docstring format](https://sphinx-rtd-tutorial.readthedocs.io/en/latest/docstrings.html) to help make your code more understandable and consistent with the rest of the toolkit.
- Ensure that your language is concise and free of grammatical errors.
- If this is your first submission to the toolkit, add your name and GitHub account to the list of [Contributors](CONTRIBUTORS.md).

### Testing

- If possible, write tests to illustrate and verify your work. We require the use of [pytest](https://docs.pytest.org). If you have doctests in your code, convert them to pytest.

### Linting

- Use [ruff](https://github.com/astral-sh/ruff) to lint your Python files before submitting your pull request. It will help identify errors and style issues in your code.

- Run `ruff check .` from the root directory of the toolkit folder to check all files in the project, excluding directories listed in the tool.ruff section of pyproject.toml.

### Formatting

- Also use [ruff](https://github.com/astral-sh/ruff) to format your Python files before submitting your pull request. It will make your code more readable and will automatically align it with much of [PEP 8](https://www.python.org/dev/peps/pep-0008/) formatting.

- Run `ruff format .` from the root directory of the toolkit folder to format all files in the project, excluding directories listed in the tool.ruff section of pyproject.toml.

- The pyproject.toml file includes a setting that sets the maximum line length used by black to 99 characters.

- If you're using VS Code, the ruff extension can be installed and set up to run on each save in Settings > search Format On Save > format a file on save.

### Type Checking

- Use [mypy](http://www.mypy-lang.org) to type check your Python files before submitting your pull request. This will help identify errors and improve code quality.

- Run `mypy .` from the root directory of the toolkit folder to check all files in the project, excluding directories listed in the tool.mypy section of pyproject.toml.

### Third Party Packages
 
- Use the [TEOS-10 GSW Toolkit](https://www.teos-10.org/software.htm) for derived physical oceanography algorithms rather than any predecessor packages.

## Code Contribution Process Overview

Below is the typical sequence for creating a submission. Feel free to [submit a question](https://github.com/Sea-BirdScientific/seabirdscientific/issues/new?template=question.md) if you need assistance.

Environment setup with native python tools:
1. Setup virtual environment `python -m venv .venv`
2. Activate environment `.venv\Scripts\activate`
3. Install all developer packages with `pip install -e .[dev]`

Environment setup with [uv](https://docs.astral.sh/uv/):
1. Setup virtual environment: `uv venv`
2. Activate environment `.venv\Scripts\activate`
3. Install all developer packages with `uv sync`

Contributions:
1. Create a fork of the seabirdscientific repository. We recommend that you only copy the main branch into your fork.
2. Clone the fork to your computer.
3. Create a new branch for your changes.
4. If this is your first submission, add your name and github username to the CONTRIBUTORS.md file.
5. Implement your changes, adhering to the guidelines above.
6. Verify that your tests pass.
7. Commit your local branch to the GitHub server.
8. Create a pull request and submit it after filling out the PR template.
9. Monitor and respond to feedback from the Sea-Bird team via comments in your pull request. We will merge your PR into the toolkit once it's passed the review process.

More details on contributing to projects in general can be found at the GitHub [contributing to a project](https://docs.github.com/en/get-started/exploring-projects-on-github/contributing-to-a-project) page.
