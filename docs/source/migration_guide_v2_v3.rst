.. _migration_guide_v2_v3:

Migration Guide: v2 to v3
#########################

Version 3.0.0 is a major release that changes the core data model of the
toolkit and reorganizes several modules. This guide describes the changes
that affect existing code and explains how to update it.

The most significant change is that the toolkit now uses ``xarray.Dataset``
as its primary data structure in place of the ``InstrumentData`` class and
``pandas.DataFrame``. Most other changes are renames or module moves, and
the previous names are retained as deprecated aliases wherever practical so
that existing code continues to run while emitting deprecation warnings.

Summary of changes
*******************

- ``xarray.Dataset`` replaces ``InstrumentData`` and ``pandas.DataFrame`` as
  the data structure returned by the file reading functions and consumed by
  ``processing`` and ``visualization``.
- ``instrument_data.cnv_to_instrument_data`` is replaced by
  ``instrument_data.read_cnv_file``.
- The ``InstrumentData`` and ``MeasurementSeries`` classes have been removed.
- Instrument-specific hex reading functions are renamed to lower snake_case.
- EOS-80 functions have moved from ``eos80_processing`` to a new
  ``eos80_conversion`` module.
- ``buoyancy`` and ``buoyancy_frequency`` have moved from ``processing`` to
  ``conversion``.
- The ``ChartData`` class has been removed from ``visualization``; charting
  functions now take an ``xarray.Dataset`` directly.
- Several parameters that took enums now take plain strings.
- Constants have been consolidated into a new ``constants`` module.
- Dependencies and the example notebook layout have changed.

Installation and dependencies
******************************

Version 3 adds ``xarray`` as a core dependency, and the optional dependency
groups have been reorganized. The single ``dev`` extra has been replaced with
more focused groups:

- ``notebooks``: dependencies for running the example notebooks.
- ``docs``: dependencies for building the documentation.

Installing the base package is unchanged:

.. code-block:: bash

   py -m pip install seabirdscientific

Optional dependencies can be included with square brackets:

.. code-block:: bash

   py -m pip install seabirdscientific[notebooks,docs]

The data model: xarray.Dataset
*******************************

In v2, files were read into an ``InstrumentData`` object that held a
dictionary of ``MeasurementSeries`` objects and exposed a
``pandas.DataFrame``. In v3, files are read directly into an
``xarray.Dataset``. Each measurand is a variable in the dataset, and metadata
such as the sample interval and start time are stored as dataset attributes.

**v2:**

.. code-block:: python

   import seabirdscientific.instrument_data as si

   data = si.cnv_to_instrument_data("example.cnv")
   temperature = data.measurements["tv290C"].values
   df = data.to_dataframe()

**v3:**

.. code-block:: python

   import seabirdscientific.instrument_data as si

   dataset = si.read_cnv_file("example.cnv")
   temperature = dataset["tv290C"].values

Variable metadata is available on each ``DataArray`` through its ``attrs``
dictionary, and file-level metadata is available on the dataset:

.. code-block:: python

   dataset["tv290C"].attrs["long_name"]   # e.g. "Temperature"
   dataset["tv290C"].attrs["units"]       # e.g. "ITS-90, deg C"
   dataset["tv290C"].attrs["sbs_name"]    # the original Sea-Bird name
   dataset.attrs["sample_interval"]
   dataset.attrs["start_time"]

Duplicate variable names in a ``.cnv`` file are disambiguated by appending an
index. For example, a second ``depSM`` column becomes ``depSM_1``. The
original name is preserved in the ``sbs_name`` attribute.

instrument_data
***************

read_cnv_file
=============

``cnv_to_instrument_data`` has been renamed to ``read_cnv_file`` and now
returns an ``xarray.Dataset``. The original function remains as a wrapper
for ``read_cnv_file`` and emits a ``DeprecationWarning``.

.. code-block:: python

   # v2
   data = si.cnv_to_instrument_data(filepath)

   # v3
   dataset = si.read_cnv_file(filepath)

read_hex_file
=============

``read_hex_file`` keeps the same name and parameters but now returns an
``xarray.Dataset`` instead of an ``InstrumentData`` object. Each hex data type
is a variable in the returned dataset.

.. code-block:: python

   dataset = si.read_hex_file(filepath, instrument_type, enabled_sensors)
   temperature_counts = dataset["temperature"].values

Removed classes
===============

The ``InstrumentData`` and ``MeasurementSeries`` classes have been removed.
Code that referenced these types should be updated to work with
``xarray.Dataset`` and ``xarray.DataArray`` instead.

Renamed hex reading functions
==============================

The instrument-specific hex reading functions have been renamed from mixed
case to lower snake_case and no longer specify a format. The old names remain
as deprecated aliases that forward to the new functions and emit a
``DeprecationWarning``.

.. list-table::
   :header-rows: 1

   * - v2 name
     - v3 name
   * - ``read_SBE39plus_format_0``
     - ``read_sbe39plus_data``
   * - ``read_seafet_format_0``
     - ``read_seafet_data``
   * - ``read_SBE911plus_format_0``
     - ``read_sbe911plus_data``
   * - ``read_SBE19plus_format_0``
     - ``read_sbe19plus_data``
   * - ``read_SBE37SM_format_0``
     - ``read_sbe37sm_data``

The ``HexDataTypes`` enum members are now deprecated. Accessing any member
emits a ``DeprecationWarning``.

conversion
**********

The signatures of the individual ``convert_*`` functions in ``conversion``
are unchanged; they still accept and return NumPy arrays.

The seawater library has replaced implementations of EOS-80 calculations.

.. code-block:: python

   # v2
   import seabirdscientific.eos80_processing as se
   rho = se.density(salinity, temperature, pressure)

   # v3
   import seawater.eos80
   rho = seawater.eos80.dens(salinity, temperature, pressure)

buoyancy moved to conversion and was split
==========================================

The ``buoyancy`` and ``buoyancy_frequency`` functions have moved from
``processing`` to ``conversion``. The versions in ``processing`` remain as
deprecated wrappers.

``buoyancy`` now only uses TEOS-10 equations, a new EOS-80 based function
is now in ``eos80_conversion ``. 

.. code-block:: python

   # v2
   import seabirdscientific.processing as p
   result = p.buoyancy(...)

   # v3
   # TEOS-10
   import seabirdscientific.conversion as c
   result = c.buoyancy(...)

  # EOS-80
   import seabirdscientific.eos80_conversion as eos80_conversion
   result = eos80c.buoyancy_eos80(...)


Note that ``conversion.buoyancy`` returns an ``xarray.Dataset``.

processing
**********

bin_average now takes a Dataset
===============================

``bin_average`` previously accepted a ``pandas.DataFrame`` and now accepts an
``xarray.Dataset``. Pass the dataset returned by ``read_cnv_file`` or
``read_hex_file`` directly.

.. code-block:: python

   dataset = si.read_cnv_file(filepath)
   binned = p.bin_average(dataset, "prdM", bin_size=1)

String parameters replace enums
===============================

Several parameters that previously took enum values now take plain strings.

- ``loop_edit`` and ``loop_edit_pressure``/``loop_edit_depth``: the
  ``min_velocity_type`` parameter now takes ``"fixed"`` or ``"percent"``
  instead of ``MinVelocityType.FIXED`` or ``MinVelocityType.PERCENT``.
- ``window_filter``: the ``window_type`` parameter now takes one of
  ``"boxcar"``, ``"cosine"``, ``"gaussian"``, ``"median"``, or ``"triangle"``
  instead of a ``WindowFilterType`` enum value.

.. code-block:: python

   # v2
   p.window_filter(data, flags, WindowFilterType.BOXCAR, width, interval)

   # v3
   p.window_filter(data, flags, "boxcar", width, interval)

find_depth_peaks
================

``find_depth_peaks`` is now an internal function named ``_find_depth_peaks``.
The public ``find_depth_peaks`` name is retained as a deprecated alias.

visualization
*************

The ``ChartData`` class has been removed. Charting functions now take an
``xarray.Dataset`` directly in place of a ``ChartData`` object, which
reflects the move to ``xarray`` throughout the toolkit.

.. code-block:: python

   # v2
   config = vis.ChartConfig(...)
   chart_data = vis.ChartData(source, config)
   figure = vis.plot_xy_chart(chart_data, config)

   # v3
   config = vis.ChartConfig(...)
   dataset = si.read_cnv_file(filepath)
   figure = vis.plot_xy_chart(dataset, config)

The ``ChartConfig`` constructor keeps the same parameters. ``parse_instrument_data``
and ``select_subset`` now operate on and return ``xarray.Dataset`` objects.

constants
*********

Numeric constants that were previously defined within individual modules are
now consolidated in a new ``constants`` module (for example ``KELVIN_OFFSET_0C``,
``ITS90_TO_IPTS68``, ``FLAG_VALUE``, and ``COUNTS_TO_VOLTS``). If your code
imported these constants from another module, import them from ``constants``
instead.

Example notebooks and data
**************************

The example notebooks and example data have moved from the ``documentation``
directory to a new ``notebooks`` directory. If you referenced these files by
path, update the paths accordingly.

Deprecation reference
*********************

The following names still work in v3 but emit a ``DeprecationWarning`` and
should be updated:

.. list-table::
   :header-rows: 1

   * - Deprecated
     - Replacement
   * - ``instrument_data.cnv_to_instrument_data``
     - ``instrument_data.read_cnv_file``
   * - ``instrument_data.read_SBE39plus_format_0``
     - ``instrument_data.read_sbe39plus_data``
   * - ``instrument_data.read_seafet_format_0``
     - ``instrument_data.read_seafet_data``
   * - ``instrument_data.read_SBE911plus_format_0``
     - ``instrument_data.read_sbe911plus_data``
   * - ``instrument_data.read_SBE19plus_format_0``
     - ``instrument_data.read_sbe19plus_data``
   * - ``instrument_data.read_SBE37SM_format_0``
     - ``instrument_data.read_sbe37sm_data``
   * - ``instrument_data.HexDataTypes``
     - (replaced with constants)
   * - ``eos80_processing.bouyancy_frequency``
     - ``eos80_conversion.bouyancy_frequency``
   * - ``eos80_processing.density``
     - ``eos80_conversion.density``
   * - ``eos80_processing.potential_temperature``
     - ``eos80_conversion.potential_temperature``
   * - ``eos80_processing.adiabatic_temperature_gradient``
     - ``eos80_conversion.adiabatic_temperature_gradient``
   * - ``processing.buoyancy``
     - ``conversion.buoyancy``
   * - ``processing.buoyancy_frequency``
     - ``conversion.buoyancy_frequency``
   * - ``processing.find_depth_peaks``
     - ``processing._find_depth_peaks``

The following have been removed with no compatibility alias:

- ``instrument_data.InstrumentData``
- ``instrument_data.MeasurementSeries``
- ``visualization.ChartData``
