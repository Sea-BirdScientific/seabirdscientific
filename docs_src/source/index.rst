.. seabirdscientific documentation master file, created by
   sphinx-quickstart on Thu Dec 26 17:58:46 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

   
Sea-Bird Scientific Community Toolkit
#####################################

The toolkit is a collection of:

- Python code to help users process data collected with SBS instruments (see the repository src/seabirdscientific folder).
- Example SBS instrument data (see the repository documentation/example_data folder).
- A Jupyter notebook that demonstrates toolkit processing options that can be applied to data collected with SBE 37 and SBE 19plus V2 CTDs (see the repository documentation folder).

Installation
************

Windows:

.. code-block:: bash

   py -m pip install seabirdscientific

Mac/Linux:

.. code-block:: bash

   python3 -m pip install seabirdscientific

Example Usage:

.. code-block:: python

   import seabirdscientific
   from seabirdscientific import conversion
   import seabirdscintific.processing as p


Contents
********

Source
======

:ref:`cal_coefficients`: Classes for storing calibration coefficients used with data conversion functions.

:ref:`contour`: A class and helper functions for storing data used with TS plots.

:ref:`conversion`: A collection of functions for converting raw instrument data to SI units. Also includes some functions for converting between SI units, for exampel from depth from pressure.

:ref:`eos80_processing`: Legacy functions for legacy data.

:ref:`instrument_data`: Classes and functions for parsing instrument data files into python data types.

:ref:`interpret_sbs_variable`: A single, long function that converts various Sea-Bird unit representations into a more limited set of strings.

:ref:`processing`: Generally contains functions that modify data that has already been converted or derived from raw data. Also contains a couple functions for deriving buoyancy.

:ref:`utils`: Helper functions for example notebooks and unit tests.

:ref:`visualization`: Functions for generating Plotly charts from Sea-Bird data.

Theory
======

:ref:`ctd_processing`: In depth descriptions of some of the algorithms used in processing.py.

.. Sidebar content

.. toctree::
   :caption: Source
   :hidden:

   cal_coefficients
   contour
   conversion
   eos80_processing
   instrument_data
   interpret_sbs_variable
   processing
   utils
   visualization

.. toctree::
   :caption: Theory
   :hidden:

   ctd_processing

.. toctree::
   :caption: Links
   :hidden:

   Github <https://github.com/Sea-BirdScientific/seabirdscientific>
   PyPI <https://pypi.org/project/seabirdscientific>

- `CTD Processing (PDF) <_static/processing.pdf>`__
- :download:`Download <_static/processing.pdf>`

