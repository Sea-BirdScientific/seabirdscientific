# seabirdscientific
SeaBird Scientific Community Toolkit  

## Installation
```bash
py -m pip install seabirdscientific --extra-index-url http://wl239.sbs.ewqg.com:9000 --trusted-host wl239.sbs.ewqg.com  
```
OR from a cloned repo, in editable mode, with dev depenedencies
```bash
py -m pip install -e .[dev]
```

## Usage
```python
import sbs
```
OR
```python
from sbs import processing, visualization
```

## Included packages
- [sbs_data_processing](https://github.com/Sea-BirdScientific/python-data-processing) as `processing`  
Data conversion and processing functions for CTDs, ported from SeaSoft, as well as [GSW](https://github.com/TEOS-10/GSW-python) derivations
- [sbs_data_visualization](https://github.com/Sea-BirdScientific/sbs-data-visualization) as `visualization`  
Data visualization using Plotly

## How to add other SeaBird repos to the seabirdscientific package
1. Follow the instructions [here](https://github.com/Sea-BirdScientific/python-package-template) to build a `.tar.gz` file.

2. Remotely connect to wl239.sbs.ewqg.com using your SeaBird SSO credentials

3. Create a subfolder for the new package `c:/SBSPyPI/your-package-name`

4. Add the `.tar.gz` file to the new subfolder

5. Modify `c:/SBSPyPI/index.html` by adding a link to your package
```html
<!DOCTYPE html>
<html>
  <body>
    ...
    <a href="/your-package-name/">your-package-name</a>
  </body>
</html>
```

6. Add a dependency to the seabirdscientific pyproject.toml  
```
dependencies = [
    ...
    "your-package-name",
]
```

7. Install the package and build a new version of seabirdscientific (remember to increment the version number in pyproject.toml)
```
py -m install your-package-name --extra-index-url http://wl239.sbs.ewqg.com:9000 --trusted-host wl239.sbs.ewqg.com  
py -m build
```

8. Add the new `seabirdscientific-*.*.*.tar.gz` file to the remote `c:/SBSPyPI/seabirdscientific` folder

9. Confirm that it all works by installing seabirdscientific to a clean virtual environment and importing your package

## How to name your packages and modules
Package names (used by pip and PyPI) should be prefixed with `sbs-`, be short, and use hyphens if it improves readability.  
Module names (used when importing in python) should be short and use underscores if it improves readability. It is not necessary to prefix module names with `sbs_` because they will already be part of the main `sbs` module.  
[PEP8](https://peps.python.org/pep-0008/#package-and-module-names)

## Notes from previous commit (work in progess)
1.  Add a repo as a git submodule
```
git submodule add https://github.com/Sea-BirdScientific/sbs-data-visualization ./submodules/sbs-data-visualization 
```

2. Create a symbolic link to the git submodule's python package folder.  
```
New-Item -ItemType SymbolicLink -Path ./seabirdscientific/sbs_data_visualization -Target ./submodules/sbs-data-visualization/sbs_data_visualization/
```

3. Add dependencies to the seabirdscientific pyproject.toml  
```
dependencies = [
    "build==1.0.3",
    ...
]
```