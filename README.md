# seabirdscientific
SeaBird Scientific Community Toolkit  

## Included packages
- `sbs_data_processing` Data conversion and processing functions for CTDs, ported from SeaSoft, as well as [GSW](https://github.com/TEOS-10/GSW-python) derivations
- `sbs_data_visualization` Data visualization using Plotly

## How to add other SeaBird repos to the seabirdscientific package
This repo includes other repos as git submodules and then uses symlinks to include the correct package folder in the main seabirdscientific folder.  
There are three steps:
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
