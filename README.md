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
from sbs.visualize import visualization as v
```
The notebooks currently depend on data files that are not part of this repo, and can't be run as-is. For now, they are only for reference.  

