# seabirdscientific
This is the repository for the Sea-Bird Scientific (SBS) Community Toolkit. It is a collection of:
- Python code to help in user developed processing of data collected with SBS instruments (see the repository src/sbs folder).
- Example SBS instrument data (see the repository documentation/example_data folder).
- A [Jupyter](https://jupyter.org/) notebook that documents the current toolkit processing options that can be applied to data collected with SBE 37 and SBE 19plus V2 CTDs (see the repository documentation folder). This notebook also serves to document the processing options available in the SBS Fathom application.

## Repostory Contents
```
seabirdscientific/
|
├── documentation/
|   ├── example_data/          [sample data files]
|   ├── ctd-processing.ipynb   [CTD Jupyter notebook]
|   ├── ctd-notebook.pdf       [static pdf notebook]
|   └── ctd-notebook.html      [static html notebook]
|
├── src/
|   └── sbs/
|       ├── process/           [processing modules]
|       └── visualize/         [visualization modules]
|
├── tests/                     [unit tests]
|
└── README.md
```
## Package Installation With pip
The seabirdscientific package uses Python 3.9 or greater. To install the package in a Python environment using pip send the command:

On Windows:
``` bash
py -m pip install seabirdscientific
```

On Unix/macOS:

``` bash
python3 -m pip install seabirdscientific
```
For additional information see the [Python.org Installing Packages](https://packaging.python.org/en/latest/tutorials/installing-packages/#installing-packages) reference.

## Example package use within python code
```python
import sbs
from sbs.process import contour
import sbs.process.conversion as conv
import sbs.process.processing as proc
```

## Required Software
You must have [Python](https://www.python.org/downloads/) version 3.9 or higher installed in order to use the toolkit.

## CTD Jupyter Notebook
The ctd-processing.ipynb notebook in the documentation folder provides examples of the methods that can be applied to SBS CTD data within both the toolkit and the SBS Fathom application. There are a number of online references available with information on Jupyter notebook setup for different platforms and environments. For those who are not interested in an interactive notebook, there are static versions available in [HTML](./documentation/ctd-notebook.html) and [PDF](./documentation/ctd-notebook.pdf).

