# seabirdscientific
Python packages by SeaBird Scientific  

## How to add other SeaBird repos to the seabirdscientific package
### Add repo as a submodule
```
git submodule add https://github.com/Sea-BirdScientific/sbs-data-visualization ./submodules/sbs-data-visualization 
```

### Create symbolic links to submodule packages
```
New-Item -ItemType SymbolicLink -Path ./seabirdscientific/sbs_data_visualization -Target ./submodules/sbs-data-visualization/sbs_data_visualization/
```



