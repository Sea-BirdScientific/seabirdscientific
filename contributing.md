# Contributing guidelines

## Before contributing

Thanks for your interest in contributing to the seabirdscientific repository. Before submitting any pull requests, please ensure that you __have fully read and understand these guidelines__. If you have any questions about contributing, please feel free to [sumbit a question](https://github.com/Sea-BirdScientific/seabirdscientific/issues/new).

## Contributing

### Contributor Agreement

We are delighted that you are considering a contribution. By being one of our contributors, you agree and confirm that:

- You created the submitted code - no plagiarism allowed. If it's based on someone else's code, please give them due credit.
- Your work will be distributed under the [MIT License](LICENSE.md) once your pull request is merged. 
  - If your work is derived from other code, please make sure that its licensing is compatible with the MIT License
- Your submitted work is in keeping with our styles and standards (see the [Coding Style](#coding-style) section).

### Coding Style

We want your work to be readable by others; therefore, we encourage you to adhere to the following:

- Please adhere to the [PEP 8](https://peps.python.org/pep-0008/) style guide for your code.
  - One exception to the PEP 8 guide is that we have chosen to allow maximum line lengths of __99__ characters vs the recommended 79 characters.
- Please write in Python 3.9 or greater.
- Please take care when naming of functions, classes, and variables. Help your reader by using __descriptive names__.
  - Please avoid single letter variable names unless their life only spans a few lines.
  - Expand acronyms because names like `gcd()` may not be obvious but `greatest_common_divisor()` is.
  - Please follow the [PEP 8 Naming Conventions](https://pep8.org/#prescriptive-naming-conventions) 
    - Variable_names and function_names should be lower_case
    - Constants should be in UPPERCASE
    - ClassNames should be CamelCase, etc.
- The use of [Python type hints](https://docs.python.org/3/library/typing.html) is encouraged for function parameters and return values.

### Formatting and style
- Please run [black](#black) on your Python file(s) before submitting your pull request. This will make your code more readable and will automatically align it with much of [PEP 8](https://www.python.org/dev/peps/pep-0008/).
- Please run [pylint](#pylint) on your Python file(s) before submitting your pull request. This will help identify errors and style issues.
- Please run [mypy](#mypy) on your Python file(s) before submitting your pull request. This will help identify errors and improve code quality.

### Documenting
- Please make use of docstrings based upon [PEP 257](https://peps.python.org/pep-0257/) to make your code more understandable.
<!-- 
TODO: Determine how to capture author information
- At a minimum, include your github account in the header of each of your .py files. OR in  
-->
- If your work is derived from other code, please include the name(s) of the developers for the originating code and where that code can be found.

### Testing
- Write tests (see [https://docs.python.org/3/library/doctest.html](https://docs.python.org/3/library/doctest.html)) to illustrate and verify your work.

<!-- ### Other stuff to consider  
- [__List comprehensions and generators__](https://docs.python.org/3/tutorial/datastructures.html#list-comprehensions) are preferred over the use of `lambda`, `map`, `filter`, `reduce` but the important thing is to demonstrate the power of Python in code that is easy to read and maintain.

- Avoid importing external libraries for basic algorithms. Only use those libraries for complicated algorithms. 
-->

### Other Requirements for Submissions
- Strictly use snake_case (underscore_separated) in your file_name, as it will be easy to parse in future using scripts.
- If you have modified/added code work, make sure the code compiles before submitting.
- If you have modified/added documentation work, ensure your language is concise and contains no grammatical errors.

## Package Use

<!-- 
TODO: Add package suggestion information, eg. Pandas vs X-Array
 -->

## Code Contribution Process

Please see the GitHub reference for [contributing to a project](https://docs.github.com/en/get-started/exploring-projects-on-github/contributing-to-a-project) for the general steps to submit your code.

We ask that you submit new features to the relevant __community__ sub-folder within the toolkit repository. If there doesn't seem to be a suitable folder, create one as part of your PR. Only bugfixes or documentation updates should include changes to the files within the src/seabirdscientific folder.

## Tools
There are a number of tools that can be used to assist in ensuring that your code is consistent with the toolkit standards. The tools listed below have been included in the \[dev\] section of the pyproject.toml file. They can be installed by running the following:

``` pip install .[dev]```

### Black

Black [https://github.com/python/black](https://github.com/python/black) is a Python code formatting tool that helps to maintain consistency amongst the .py files in the toolkit. 

It can be run on individual files or directories by running ```black {source file or directory}```. Running ```black src``` from the root directory of the toolkit folder will format all the .py files within src directory and its sub-directories. 

The pyproject.toml file includes a setting that sets the maximum line length used by black to 99 characters.

### Pylint

Pylint [https://github.com/pylint-dev/pylint](https://github.com/pylint-dev/pylint) analyzes your code for errors and style issues.

### MyPy

MyPy [http://www.mypy-lang.org](http://www.mypy-lang.org) is a static type checker for Python that helps to improve code quality.

Instructions on how to install mypy can be found at [https://github.com/python/mypy](https://github.com/python/mypy). Please use the command `mypy --ignore-missing-imports .` to test all files or `mypy --ignore-missing-imports path/to/file.py` to test a specific file.
