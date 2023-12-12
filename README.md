# seabirdscientific
Sea-Bird Scientific Community Toolkit

## Installation

```bash
py -m pip install seabirdscientific --extra-index-url http://winbuilder02.sbs.ewqg.com/SBSPyPi/ --trusted-host winbuilder02.sbs.ewqg.com
```
OR from a cloned repo, in editable mode, with dev depenedencies
```bash
git clone https://github.com/Sea-BirdScientific/seabirdscientific.git
py -m pip install -e .[dev]
```

## Examples
```python
import sbs
from sbs.process import contour
from sbs.visualize import visualization as viz
```
There is a ctd-processing notebook \([ctd-processing.ipynb](https://github.com/Sea-BirdScientific/seabirdscientific/blob/Install-from-winbuilder02/documentation/ctd-processing.ipynb)\) in the documentation folder of the repository. It shows examples of how to use some of the CTD related processing and visualization toolkit functions. It uses sample data from the included documentation\example_data folder.
