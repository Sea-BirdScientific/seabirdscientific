# seabirdscientific
Python packages by SeaBird Scientific  

## How to add other SeaBird repos to the seabirdscientific package
### Add repo as a submodule
```
git submodule add https://github.com/Sea-BirdScientific/sbs-data-visualization ./submodules/sbs-data-visualization 
```

### Create symbolic links to submodule packages
Note: the name of the symbolic link does not need to match the submodule, as shown below where dataviz maps to sbs_data_visualization.
```
New-Item -ItemType SymbolicLink -Path ./seabirdscientific/dataviz -Target ./submodules/sbs-data-visualization/sbs_data_visualization/
```



