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
******************

- Python 3.11 or later is required.
- ``xarray.Dataset`` replaces ``InstrumentData`` and ``pandas.DataFrame`` as
  the data structure returned by the file reading functions and consumed by
  ``processing`` and ``visualization``.
- ``instrument_data.cnv_to_instrument_data`` is replaced by
  ``instrument_data.read_cnv_file``, which reads Fathom ``.cnv`` files by
  default.
- The ``InstrumentData`` and ``MeasurementSeries`` classes have been removed.
- Instrument-specific hex reading functions are renamed to lower snake_case.
- ``HexDataTypes`` is replaced by ``HEX_TYPE_*`` string constants.
- EOS-80 functions have moved from ``eos80_processing`` to a new
  ``eos80_conversion`` module or to the ``seawater`` library.
- ``buoyancy`` and ``buoyancy_frequency`` have moved from ``processing`` to
  ``conversion`` and now return NumPy arrays.
- ``convert_pressure`` now returns psia by default, and ``convert_conductivity``
  takes an ``instrument_type`` parameter in place of ``scalar``.
- ``loop_edit_depth`` and ``loop_edit_pressure`` are replaced by ``loop_edit``,
  which returns flags instead of a cast array.
- ``split`` now labels casts with a ``cast_type`` coordinate instead of
  returning a list of DataFrames.
- The ``ChartData`` class has been removed from ``visualization``; charting
  functions now take an ``xarray.Dataset`` directly.
- Several parameters that took enums now take plain strings.
- Several parameters with unit suffixes have been renamed.
- Constants have been consolidated into a new ``constants`` module.
- Several deprecated v2 aliases and parameters have been removed.
- Dependencies and the example notebook layout have changed.

Installation and dependencies
*****************************

Version 3 requires Python 3.11 or later. Support for Python 3.9 and 3.10 has
been removed.

Version 3 adds ``xarray`` and ``seawater`` as core dependencies, and all core
dependencies are now pinned to compatible release ranges. The optional
dependency groups have been reorganized. The ``dev`` extra has been replaced
with more focused groups:

- ``notebooks``: dependencies for running the example notebooks.
- ``docs``: dependencies for building the documentation.

Development dependencies are now only available through the ``dev``
dependency group (for use with uv).

Installing the base package is unchanged:

.. code-block:: bash

   py -m pip install seabirdscientific

Optional dependencies can be included with square brackets:

.. code-block:: bash

   py -m pip install seabirdscientific[notebooks,docs]

The data model: xarray.Dataset
******************************

In v2, ``.cnv`` files were read into an ``InstrumentData`` object that held a
dictionary of ``MeasurementSeries`` objects and exposed a
``pandas.DataFrame``, and ``.hex`` files were read into a ``pandas.DataFrame``.
In v3, both are read directly into an ``xarray.Dataset``. Each measurand is a
variable in the dataset along the ``scan`` dimension, and metadata such as the
sample interval and start time are stored as dataset attributes.

**v2:**

.. code-block:: python

   import seabirdscientific.instrument_data as si

   data = si.cnv_to_instrument_data("example.cnv")
   temperature = data.measurements["tv290C"].values
   df = data.to_dataframe()

**v3:**

.. code-block:: python

   import seabirdscientific.instrument_data as si

   dataset = si.read_cnv_file("example.cnv", software="seasoft")
   temperature = dataset["tv290C"].values
   df = dataset.to_dataframe()

Variable metadata is available on each ``DataArray`` through its ``attrs``
dictionary, and file-level metadata is available on the dataset:

.. code-block:: python

   dataset["tv290C"].attrs["long_name"]   # e.g. "Temperature"
   dataset["tv290C"].attrs["units"]       # e.g. "ITS-90, deg C"
   dataset["tv290C"].attrs["sbs_name"]    # the original Sea-Bird name
   dataset.attrs["file_name"]
   dataset.attrs["sample_interval"]
   dataset.attrs["start_time"]

Duplicate variable names in a ``.cnv`` file are disambiguated by appending an
index. For example, a second ``depSM`` column becomes ``depSM_1``. The
original name is preserved in the ``sbs_name`` attribute.

instrument_data
***************

read_cnv_file
=============

``cnv_to_instrument_data`` has been replaced by ``read_cnv_file``, which
returns an ``xarray.Dataset``. ``read_cnv_file`` takes a ``software``
parameter that selects the file format: ``"fathom"`` (default) or
``"seasoft"``. Files created by SBE Data Processing (SeaSoft) must be read
with ``software="seasoft"``.

``cnv_to_instrument_data`` remains as a wrapper for
``read_cnv_file(filepath, "seasoft")`` and emits a ``DeprecationWarning``.

.. code-block:: python

   # v2
   data = si.cnv_to_instrument_data(filepath)

   # v3
   fathom_dataset = si.read_cnv_file(filepath)
   seasoft_dataset = si.read_cnv_file(filepath, software="seasoft")

When reading a SeaSoft file, a ``scan`` column in the file is stored as
``scan_count`` because ``scan`` is used as the dataset dimension. When reading
a Fathom file, a ``scan`` column replaces the values of the ``scan``
coordinate.

read_hex_file
=============

``read_hex_file`` keeps the same name and parameters but now returns an
``xarray.Dataset`` instead of a ``pandas.DataFrame``. Each hex data type is a
variable in the returned dataset. The ``enabled_sensors`` parameter no longer
has a default and must be provided.

.. code-block:: python

   dataset = si.read_hex_file(filepath, instrument_type, enabled_sensors)
   temperature_counts = dataset["temperature"].values

Removed classes and functions
=============================

The ``InstrumentData`` and ``MeasurementSeries`` classes have been removed.
Code that referenced these types should be updated to work with
``xarray.Dataset`` and ``xarray.DataArray`` instead.

``fix_exponents`` has been removed. The ``.cnv`` readers handle exponent
formatting internally.

Renamed hex reading functions
=============================

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

The deprecated ``hex`` parameter has been removed from ``read_hex`` and the
hex reading functions. Use ``hex_segment`` instead.

Hex data types and lengths
==========================

The ``HexDataTypes`` enum is deprecated, and accessing any member emits a
``DeprecationWarning``. Use the module-level ``HEX_TYPE_*`` string constants
instead. The values are unchanged.

.. code-block:: python

   # v2
   dataset[si.HexDataTypes.temperature.value]

   # v3
   dataset[si.HEX_TYPE_TEMPERATURE]

The ``HEX_LENGTH`` dictionary has been removed and replaced with
``HEX_LEN_*`` constants, for example ``HEX_LEN_TEMPERATURE`` and
``HEX_LEN_DATE_TIME``.

cal_coefficients
****************

- ``PressureDigiquartzCoefficients``: ``AD590M`` and ``AD590B`` have been
  renamed to ``ad590m`` and ``ad590b``. New ``slope`` (default 1) and
  ``offset`` (default 0, dbar) fields have been added.
- ``TemperatureCoefficients`` and ``PressureCoefficients`` have a new
  ``offset`` field (default 0).
- ``PHSeaFETInternalCoefficients`` is now a dataclass. ``kdf0`` and ``kdf2``
  are required, and the deprecated ``k0``, ``k2``, ``int_k0``, and ``int_k2``
  parameters have been removed.
- ``PHSeaFETExternalCoefficients`` is now a dataclass. The deprecated
  ``ext_k0`` and ``ext_k2`` parameters have been removed, and
  ``k2_poly_order`` and ``fp_poly_order`` are now the last fields. Code that
  constructs this class with positional arguments must be updated; keyword
  arguments are recommended.

conversion and eos80_conversion
*******************************

The individual ``convert_*`` functions still accept and return NumPy arrays.
The following changes affect existing code.

Pressure
========

- ``convert_pressure``: the default ``units`` is now ``"psia"`` instead of
  ``"psig"``. Pass ``units="psig"`` to keep the v2 behavior. The coefficient
  ``offset`` is applied in dbar.
- ``convert_pressure_digiquartz``: ``units="psia"`` now returns absolute
  pressure. In v2 it returned gauge pressure. ``"psig"`` is now supported.
  The coefficient ``slope`` and ``offset`` are applied in dbar. The
  temperature compensation average now uses a trailing 30 second window with
  edge padding, which can produce small numerical differences.

Conductivity
============

``convert_conductivity`` now takes a required ``instrument_type`` parameter
in place of ``scalar``. The SBE 911plus scaling is applied automatically, and
the ``wbotc`` term is not applied for the SBE 16plus, SBE 19plus, and
SBE 911plus. A new ``units`` parameter selects ``"S/m"`` (default),
``"mS/cm"``, or ``"uS/cm"``.

.. code-block:: python

   # v2
   c.convert_conductivity(counts, temperature, pressure, coefs, scalar=0.1)

   # v3
   c.convert_conductivity(
       counts, temperature, pressure, coefs, si.InstrumentType.SBE911Plus
   )

SeaFET pH
=========

``convert_internal_seafet_ph`` and ``convert_external_seafet_ph`` no longer
provide defaults for the data and coefficient parameters, and the deprecated
``ph_counts`` parameter has been removed. Use ``raw_ph`` with
``ph_units="counts"``.

Other parameter changes
=======================

The following parameters are new and default to the v2 behavior:

- ``convert_sbe43_oxygen`` and ``convert_sbe63_oxygen``: ``units`` selects
  the output units (default ``"ml/l"``). ``convert_sbe63_oxygen`` also takes
  an optional ``external_temperature``.
- ``depth_from_pressure``: ``water_type`` selects ``"salt"`` (default) or
  ``"fresh"``.
- ``convert_altimeter``: ``units`` selects ``"m"`` (default) or ``"ft"``.

EOS-80 functions
================

The EOS-80 functions in ``eos80_processing`` are deprecated. Density and
potential temperature are now provided by the ``seawater`` library, and the
remaining functions have moved to ``eos80_conversion``.

.. code-block:: python

   # v2
   import seabirdscientific.eos80_processing as se
   rho = se.density(salinity, temperature, pressure)

   # v3
   import seawater
   rho = seawater.dens(salinity, temperature, pressure)

Note that the deprecated ``eos80_processing.density`` now returns density
minus 1000 kg/m^3 (sigma). In v2 it returned full density.

buoyancy moved to conversion and was split
==========================================

The ``buoyancy`` and ``bouyancy_frequency`` functions have moved from
``processing`` to ``conversion``, and ``bouyancy_frequency`` has been renamed
to ``buoyancy_frequency``. The versions in ``processing`` remain as
deprecated wrappers.

``conversion.buoyancy`` only uses TEOS-10 equations. The EOS-80 calculation
is now ``eos80_conversion.buoyancy_eos80``. Passing
``use_modern_formula=False`` to ``conversion.buoyancy`` emits a
``DeprecationWarning`` and returns arrays filled with the flag value.

.. code-block:: python

   # v2
   import seabirdscientific.processing as p
   result = p.buoyancy(..., use_modern_formula=True)
   result = p.buoyancy(..., use_modern_formula=False)

   # v3
   # TEOS-10
   import seabirdscientific.conversion as c
   n2, n, e, e_scaled = c.buoyancy(...)

   # EOS-80
   import seabirdscientific.eos80_conversion as ec
   n2, n, e, e_scaled = ec.buoyancy_eos80(...)

In v2, ``buoyancy`` returned a ``pandas.DataFrame`` with the columns ``N2``,
``N``, ``E``, and ``E10^-8``. In v3, both functions return a tuple of four
NumPy arrays in the same order. Where N^2 is negative, N is now the negative
square root of the absolute value to match SeaSoft. In v2 it was NaN.

New functions
=============

The following functions are new in v3:

- ``conversion``: ``convert_temperature_units``, ``convert_pressure_units``,
  ``convert_conductivity_units``, ``convert_oxygen_units``,
  ``convert_oxygen_to_umol_per_l``, ``derive_descent_rate``,
  ``derive_acceleration``, ``derive_oxygen_saturation_gg``,
  ``derive_oxygen_saturation_w``, and ``derive_nitrogen_saturation``.
- ``eos80_conversion``: ``derive_potential_temperature_anomaly``,
  ``derive_thermosteric_anomaly``, ``derive_sound_velocity``,
  ``derive_specific_conductance``, ``derive_average_sound_velocity``, and
  ``derive_gpa``.

contour
*******

The parameters of ``contour_from_t_s_p`` and ``contour_from_t_c_p`` have been
renamed. The old names remain as deprecated keyword arguments.

.. list-table::
   :header-rows: 1

   * - v2 name
     - v3 name
   * - ``temperature_C``
     - ``temperature``
   * - ``salinity_PSU``
     - ``salinity``
   * - ``conductivity_mScm``
     - ``conductivity``
   * - ``pressure_dbar``
     - ``pressure``

A new ``reference_pressure`` parameter (default 0) sets the reference
pressure for potential density. Samples with salinity below ``min_salinity``
are now set to NaN instead of being removed, so the output arrays have the
same length as the input arrays.

processing
**********

bin_average
===========

``bin_average`` previously accepted and returned a ``pandas.DataFrame``. It
now accepts and returns an ``xarray.Dataset``. The input must have a ``scan``
dimension, such as the dataset returned by ``read_cnv_file`` or
``read_hex_file``. The output uses a ``bin_number`` dimension and keeps the
input dataset attributes. It now raises a ValueError is CastType is invalid.

.. code-block:: python

   dataset = si.read_cnv_file(filepath, software="seasoft")
   binned = p.bin_average(dataset, "prdM", bin_size=1)

loop_edit
=========

``loop_edit_depth`` and ``loop_edit_pressure`` are deprecated and replaced by
``loop_edit``. ``loop_edit`` takes depth by default, or pressure with
``units="pressure"`` and a ``latitude``. Most parameters now have defaults.

The return value has changed. In v2, the functions modified the ``flag``
array in place and returned a cast array (-1 for downcast, 1 for upcast, 0
otherwise). In v3, ``loop_edit`` returns a new flag array and does not modify
the input.

.. code-block:: python

   # v2
   p.loop_edit_pressure(pressure, latitude, flag, interval, MinVelocityType.FIXED, ...)
   dataset["flag"] = flag

   # v3
   dataset["flag"] = p.loop_edit(
       pressure, flag, interval, "fixed", ..., latitude=latitude, units="pressure"
   )

The ``use_deck_pressure_offset`` option now adds the first depth value to the
minimum and maximum soak depths. In v2 it was subtracted.

split and cast selection
========================

``split`` now takes an ``xarray.Dataset`` and returns a copy with a
``cast_type`` coordinate that labels each scan as ``"downcast"``,
``"upcast"``, or ``None``. In v2 it returned a list of DataFrames. The
``control_name`` and ``split_mode`` parameters have been renamed to
``control_variable`` and ``cast_type``. New ``min_value`` and ``drop``
parameters have been added. Set ``drop=True`` to remove unlabeled scans.

.. code-block:: python

   # v2
   downcast, upcast = p.split(df, "prdM", CastType.BOTH)

   # v3
   ds = p.split(dataset, "prdM", CastType.BOTH)
   downcast = ds.where(ds["cast_type"] == "downcast", drop=True)
   upcast = ds.where(ds["cast_type"] == "upcast", drop=True)

``get_downcast`` and ``get_upcast`` have been removed. Use ``split``, or the
new ``find_downcast`` and ``find_upcast`` functions, which return the start
and end indices of a cast from depth and flag arrays. A new ``trim`` function
trims a DataFrame to a range of a control variable.

The ``CastType`` enum values are now strings (``"both"``, ``"downcast"``,
``"upcast"``, ``""``), and ``CastType.NA`` has been renamed to
``CastType.NONE``. ``bin_average`` raises a ``ValueError`` for an invalid
``cast_type``.

String literals replace enums
===============================

Several parameters that previously took enum values now take string literals.
The enums are deprecated and emit a ``DeprecationWarning`` when passed.

- ``loop_edit``: the ``min_velocity_type`` parameter now takes ``"fixed"`` or
  ``"percent"`` instead of ``MinVelocityType.FIXED`` or
  ``MinVelocityType.PERCENT``.
- ``window_filter``: the ``window_type`` parameter now takes one of
  ``"boxcar"``, ``"cosine"``, ``"gaussian"``, ``"median"``, or ``"triangle"``
  instead of a ``WindowFilterType`` enum value.

.. code-block:: python

   # v2
   p.window_filter(data, flags, WindowFilterType.BOXCAR, width, interval)

   # v3
   p.window_filter(data, flags, "boxcar", width, interval)


window_filter
=============

``window_filter`` now increments an even ``window_width`` by one so that the
window has a center point.

cell_thermal_mass
=================

The ``temperature_C`` and ``conductivity_Sm`` parameters have been renamed to
``temperature`` and ``conductivity``. The old names remain as deprecated
keyword arguments.

Removed aliases
===============

The deprecated v2 aliases ``min_velocity_mask``, ``mean_speed_percent_mask``,
``flag_by_minima_maxima``, and ``flag_data`` have been removed. These were
internal functions and have no public replacement.

visualization
*************

The ``ChartData`` class has been removed. Charting functions now take an
``xarray.Dataset`` directly in place of a ``ChartData`` object.

.. code-block:: python

   # v2
   config = vis.ChartConfig(...)
   chart_data = vis.ChartData(source, config)
   figure = vis.plot_xy_chart(chart_data, config)

   # v3
   config = vis.ChartConfig(...)
   dataset = si.read_cnv_file(filepath, software="seasoft")
   figure = vis.plot_xy_chart(dataset, config)

The ``ChartConfig`` constructor keeps the same parameters. ``create_single_plot``,
``create_subplots``, and ``create_overlay`` now take ``xarray.Dataset``
objects for ``x`` and ``y``.

``parse_instrument_data`` now accepts an ``xarray.Dataset`` in addition to
the v2 source types and returns an ``xarray.Dataset``. ``select_subset`` now
operates on and returns an ``xarray.Dataset``. When called with an empty name
list it returns a variable named ``"Scan Count"`` instead of
``"Sample Count"``, and it raises a ``KeyError`` for names not in the
dataset.

interpret_sbs_variable
**********************

Variable definitions are now loaded from ``resources/sbs_variables.json``.
``interpret_sbs_variable`` returns a dictionary with the additional keys
``aliases``, ``title``, and ``info``. Unit strings no longer include the
calibration standard or sensor information, which is now in ``info``. For
example, ``t090C`` has units ``"deg C"`` and info ``["ITS-90"]``. In v2 the
units were ``"ITS-90, deg C"``.

constants
*********

Numeric constants that were previously defined within individual modules are
now consolidated in a new ``constants`` module (for example
``KELVIN_OFFSET_0C``, ``ITS90_TO_IPTS68``, ``FLAG_VALUE``, and
``COUNTS_TO_VOLTS``). Constants such as ``conversion.KELVIN_OFFSET_0C`` are no
longer available from their original modules. Import them from ``constants``
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
     - ``instrument_data.HEX_TYPE_*`` constants
   * - ``eos80_processing.bouyancy_frequency``
     - ``eos80_conversion.bouyancy_frequency``
   * - ``eos80_processing.density``
     - ``seawater.dens``
   * - ``eos80_processing.potential_temperature``
     - ``seawater.ptmp``
   * - ``eos80_conversion.potential_temperature``
     - ``seawater.ptmp``
   * - ``eos80_processing.adiabatic_temperature_gradient``
     - ``eos80_conversion.adiabatic_temperature_gradient``
   * - ``processing.buoyancy``
     - ``conversion.buoyancy``
   * - ``processing.buoyancy(use_modern_formula=False)``
     - ``eos80_conversion.buoyancy_eos80``
   * - ``processing.bouyancy_frequency``
     - ``conversion.buoyancy_frequency``
   * - ``processing.loop_edit_depth``
     - ``processing.loop_edit``
   * - ``processing.loop_edit_pressure``
     - ``processing.loop_edit(..., units="pressure")``
   * - ``processing.MinVelocityType``
     - ``"fixed"`` or ``"percent"``
   * - ``processing.WindowFilterType``
     - ``"boxcar"``, ``"cosine"``, ``"gaussian"``, ``"median"``, or ``"triangle"``
   * - ``processing.find_depth_peaks``
     - ``processing._find_depth_peaks``
   * - ``processing.cell_thermal_mass`` parameters ``temperature_C``, ``conductivity_Sm``
     - ``temperature``, ``conductivity``
   * - ``contour`` parameters ``temperature_C``, ``salinity_PSU``, ``conductivity_mScm``, ``pressure_dbar``
     - ``temperature``, ``salinity``, ``conductivity``, ``pressure``

The following have been removed with no compatibility alias:

- ``instrument_data.InstrumentData``
- ``instrument_data.MeasurementSeries``
- ``instrument_data.fix_exponents``
- ``instrument_data.HEX_LENGTH``
- The ``hex`` parameter of ``instrument_data.read_hex`` and the hex reading
  functions
- ``processing.get_downcast`` and ``processing.get_upcast``
- ``processing.min_velocity_mask``, ``processing.mean_speed_percent_mask``,
  ``processing.flag_by_minima_maxima``, and ``processing.flag_data``
- ``processing.CastType.NA`` (renamed to ``CastType.NONE``)
- The ``scalar`` parameter of ``conversion.convert_conductivity``
- The ``ph_counts`` parameter of ``conversion.convert_internal_seafet_ph`` and
  ``conversion.convert_external_seafet_ph``
- The ``AD590M`` and ``AD590B`` fields of
  ``cal_coefficients.PressureDigiquartzCoefficients`` (renamed to ``ad590m``
  and ``ad590b``)
- The ``k0``, ``k2``, ``int_k0``, and ``int_k2`` parameters of
  ``cal_coefficients.PHSeaFETInternalCoefficients``
- The ``ext_k0`` and ``ext_k2`` parameters of
  ``cal_coefficients.PHSeaFETExternalCoefficients``
- Constants defined in ``conversion`` (moved to ``constants``)
- ``visualization.ChartData``
