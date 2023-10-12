# seabirdscientific
Sea-Bird Scientific Community Toolkit

## Installation
```bash
py -m pip install seabirdscientific --extra-index-url http://wl239.sbs.ewqg.com:9000 --trusted-host wl239.sbs.ewqg.com  
```
OR from a cloned repo, in editable mode, with dev depenedencies
```bash
py -m pip install -e .[dev]
```

## Examples
```python
import sbs
from sbs.processing import contour
from sbs.visualization import data_visualization as dv
```
The notebooks currently depend on data files that are not part of this repo, and can't be run as-is. For now, they are only for reference.  

## Included submodules
- [sbs_data_processing](https://github.com/Sea-BirdScientific/python-data-processing) as `processing`  
Data conversion and processing functions for CTDs, ported from SeaSoft, as well as [GSW](https://github.com/TEOS-10/GSW-python) derivations
- [sbs_data_visualization](https://github.com/Sea-BirdScientific/sbs-data-visualization) as `visualization`  
Data visualization using Plotly

## How to add other SeaBird repos to the seabirdscientific package
This project combines multiple Sea-Bird repos into one python package by adding repos as git submodules and then creating symlinks from the main python module folder to the sbs folder.  
There are 6 easy steps:
1.  From the root directory of this project, add your repo as a git submodule
```bash
git submodule add https://github.com/Sea-BirdScientific/sbs-your-project ./submodules/sbs-your-project 
```

2. Create a symbolic link to the git submodule's python module folder.  
This is an opportunity to rename the module for the community toolkit, which is imported as `sbs`. So if your module is called `sbs_your_project`, consider creating the symlink as simply `your_project`.  
The target for the symlink should be the top-level module folder within the git repo, which is not necessarily the entire git repo.
```
New-Item -ItemType SymbolicLink -Path ./sbs/your_project -Target ./submodules/sbs-your-project/sbs_your_project/
```

3. Add dependencies to the seabirdscientific pyproject.toml  
```
dependencies = [
    ...
    "your_project"
]
```

4. Build the `.whl` and `.tar.gz` files (remember to increment the version number)
```bash
py -m pip install build
py -m build
```

5. Connect remotely to `wl239.sbs.ewqg.com` with your Sea-Bird SSO

6. Add the new `seabirdscientific-*.*.*.tar.gz` file to the temporary package index `c:/SBSPyPI/seabirdscientific`