# Contributing Guidelines

## Before Contributing

Thanks for your interest in contributing to the seabirdscientific repository. Before submitting any pull requests, ensure that you __have fully read and understand these guidelines__. If you're new to the process, take a look at the [Code Contribution Process Overview](#code-contribution-process-overview) section. If you've read through everything and still have any questions about contributing, please feel free to [submit a question](https://github.com/Sea-BirdScientific/seabirdscientific/issues/new).

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
    - Variable_names and function_names should be lower_case
    - Constants should be in UPPERCASE
    - ClassNames should be CamelCase, etc.
- The use of [Python type hints](https://docs.python.org/3/library/typing.html) are encouraged for function parameters and return values.
- Use Python 3.9 or greater.

### Documenting

- Make use of docstrings based upon the [Sphinx docstring format](https://sphinx-rtd-tutorial.readthedocs.io/en/latest/docstrings.html) to help make your code more understandable and consistent with the rest of the toolkit.
- Ensure that your language is concise and free of grammatical errors.
- If this is your first submission to the toolkit, add your name and GitHub account to the list of [Contributors](CONTRIBUTORS.md).

### Testing

- If possible, write tests to illustrate and verify your work. We require the use of [pytest](https://docs.pytest.org). If you have doctests in your code, convert them to pytest.

### Formatting and Style Tools

- Run [black](https://github.com/python/black) on your Python file(s) before submitting your pull request. It will make your code more readable and will automatically align it with much of [PEP 8](https://www.python.org/dev/peps/pep-0008/) formatting.

  Black can be run on individual files or directories by running `black path/to/file.py`. Running `black src` from the root directory of the toolkit folder will format all the .py files within src directory and its sub-directories.

  The pyproject.toml file includes a setting that sets the maximum line length used by black to 99 characters.

- Run [pylint](https://github.com/pylint-dev/pylint) on your Python file(s) before submitting your pull request. It will help identify errors and style issues in your code.

  Pylint can be used to check an individual file by running `pylint path/to/file.py`. Running `pylint src` from the root directory of the toolkit folder will check all files in the src directory.

  __Note:__ There are currently recursion depth issues in `src\seabirdscientific\interpret_sbs_variable.py` and `src\seabirdscientific\visualization.py` that will cause pylint errors. To lint the toolkit excluding those files, run <nobr>`pylint --ignore-patterns=".*(interpret_sbs_variable|visualization)\.py" src`<nobr>

- Run [mypy](http://www.mypy-lang.org) on your Python file(s) before submitting your pull request. This will help identify errors and improve code quality.

  Mypy can be used to check an individual file by running `mypy path/to/file.py`. Running `mypy src` from the root directory of the toolkit folder will check all files in the src directory.

- Ensure that your code compiles before submitting by running `python -m compileall ./src` from the root directory of the toolkit folder.

## Package Use
 
- Use the [TEOS-10 GSW Toolkit](https://www.teos-10.org/software.htm) for derived physical oceanography algorithms rather than any predecessor packages.
- When working with 2 dimensional data, numpy is preferred over xarray.
- Xarray is recommended when working with multidimensional data.

## Code Contribution Process Overview

Below is the typical sequence for creating a submission. Feel free to [submit a question](https://github.com/Sea-BirdScientific/seabirdscientific/issues/new) if you need assistance.

1. Create a fork of the seabirdscientific repository. We recommend that you only copy the main branch into your fork.
1. Clone the fork to your computer.
1. Create a new branch for your changes.
1. If this is your first submission, add your name and github username to the CONTRIBUTORS.md file.
1. Implement your changes, adhering to the [Coding Style](#coding-style) and [Documenting](#documenting) guidelines above.
1. Run all the processes listed in the [Formatting and Style Tools](#formatting-and-style-tools) section above and ensure they are free of issues.
1. Verify that your tests pass.
1. Commit your local branch to the GitHub server.
1. Create a pull request and submit it after filling out the PR template.
1. Monitor and respond to feedback from the Sea-Bird team via comments in your pull request. We will merge your PR into the toolkit once it's passed the review process.

More details on contributing to projects in general can be found at the GitHub [contributing to a project](https://docs.github.com/en/get-started/exploring-projects-on-github/contributing-to-a-project) page.
