.. _migration_guide:

Migration Guide
###############

v2.7.3 → v3.0.0
****************

Breaking Changes Summary
========================

*   **Pandas → Xarray**

    *   `read_hex_file()` now returns :class:`xarray.Dataset`.
    *   `bin_average()` now returns :class:`xarray.Dataset`.
    *   New `read_cnv_file()` replaces the old `cnv_to_instrument_data` pattern.

*   **Function Rename to Snake Case**

    Old names (`read_SBE19plus_format_0` etc.) still exist but emit a
    `DeprecationWarning`. New canonical versions:

    *   `read_sbe19plus_format_0`
    *   `read_sbe37sm_format_0`
    *   `read_sbe39plus_format_0`
    *   `read_sbe911plus_format_0`

*   **Enums Deprecated at Call Sites**

    *   `MinVelocityType` and `WindowFilterType` are deprecated at call sites.
        Functions now accept string literals instead.

*   **Hex Metadata Restructured**

    *   New `HEX_TYPE_*` and `HEX_LEN_*` constants.
    *   `HexDataTypes` Enum remains but triggers warnings when accessed.

*   **Deprecated Parameters Removed**

    *   `ph_counts`, `k0`, `k2`, `int_k0`, `ext_k0` and related variants
        removed from multiple APIs.

*   **New Utilities**

    *   `profile` decorator for local line profiling.
    *   `WarnAllMembersMeta` to warn on deprecated Enum access.

Module-level API Changes
========================

`cal_coefficients` module
-------------------------

Removed deprecated constructor parameters:

*   `PHSeaFETInternalCoefficients` no longer accepts:
    `k0`, `k2`, `int_k0`, `int_k2`.
*   `PHSeaFETExternalCoefficients` no longer accepts:
    `ext_k0`, `ext_k2`.

**Migration Example**

.. code-block:: python

   # Before
   cc = PHSeaFETInternalCoefficients(k0=1.0, k2=2.0)

   # After
   cc = PHSeaFETInternalCoefficients(kdf0=1.0, kdf2=2.0)

.. code-block:: python

   # Before
   cc = PHSeaFETExternalCoefficients(ext_k0=1.0, ext_k2=2.0)

   # After
   cc = PHSeaFETExternalCoefficients(k0=1.0, k2=2.0)


`conversion` module
-------------------

Removed `ph_counts` from:

*   `convert_internal_seafet_ph`
*   `convert_external_seafet_ph`

**Migration Example**

.. code-block:: python

   # Before
   ph = convert_internal_seafet_ph(raw_ph=None, ph_counts=arr)

   # After
   ph = convert_internal_seafet_ph(raw_ph=arr, ph_units="counts")


`eos80_processing` module
-------------------------

`bouyancy_frequency` renamed parameters:

*   New canonical names: `temperature`, `salinity`, `pressure`.
*   Old names still accepted with warnings.

**Migration Example**

.. code-block:: python

   # Before
   n2 = bouyancy_frequency(
       temp_ITS_subset=t,
       salinity_prac_subset=s,
       pressure_dbar_subset=p,
       gravity=g
   )

   # After
   n2 = bouyancy_frequency(
       temperature=t,
       salinity=s,
       pressure=p,
       gravity=g
   )

`density` similarly renamed parameters.


`instrument_data` module
------------------------

CNV Reading (New API)

New reader:

*   `read_cnv_file(filepath) -> xarray.Dataset`  
    Replaces `InstrumentData` / `MeasurementSeries` / `cnv_to_instrument_data`.

**Migration Example**

.. code-block:: python

   # Before
   inst = cnv_to_instrument_data("file.cnv")
   df = inst._to_dataframe()

   # After
   ds = read_cnv_file("file.cnv")
   df = ds.to_dataframe()

HEX Reading and Decoding

*   `read_hex_file` now returns `xarray.Dataset` instead of `pandas.DataFrame`.
*   Added `_preallocate_dataset` helper.
*   Introduced `HEX_TYPE_*` and `HEX_LEN_*` constants.
*   `HexDataTypes` still present but every member access triggers a
    `DeprecationWarning`.

**Migration Example**

.. code-block:: python

   # Before
   df = read_hex_file("cast.hex", itype, enabled_sensors=[...])

   # After
   ds = read_hex_file("cast.hex", itype, enabled_sensors=[...])

Snake-case Format Readers

Old → new names:

*   `read_SBE19plus_format_0` → `read_sbe19plus_format_0`
*   `read_SBE37SM_format_0` → `read_sbe37sm_format_0`
*   `read_SBE39plus_format_0` → `read_sbe39plus_format_0`
*   `read_SBE911plus_format_0` → `read_sbe911plus_format_0`

**Migration Example**

.. code-block:: python

   # Before
   vals = read_SBE37SM_format_0(line, enabled_sensors=...)

   # After
   vals = read_sbe37sm_format_0(line, enabled_sensors=...)


`processing` module
-------------------

Deprecated Enums at Call Sites

Use string literals instead of `MinVelocityType` or `WindowFilterType`.

**Migration Example**

.. code-block:: python

   # Before
   m = loop_edit_depth(depth, flag, dt, MinVelocityType.FIXED, 0.1)

   # After
   m = loop_edit_depth(depth, flag, dt, "fixed", 0.1)

.. code-block:: python

   # Before
   y = window_filter(x, flags, WindowFilterType.MEDIAN, 5, dt)

   # After
   y = window_filter(x, flags, "median", 5, dt)

Xarray-Based `bin_average`

*   Input: now an `xarray.Dataset`.
*   Output: `xarray.Dataset` with `bin_number` dimension.

**Migration Example**

.. code-block:: python

   # Before
   df = read_hex_file(...)
   avg = bin_average(df, "prdM", 1.0)

   # After
   ds = read_hex_file(...)
   avg = bin_average(ds, "prdM", 1.0)

Cast Splitting Now Uses Coordinates

`split` now:

*   Adds a `cast_type` coordinate.
*   Does *not* return separate datasets unless `drop=True` is used.

**Migration Example**

.. code-block:: python

   labeled = split(ds, "prdM", cast_type=CastType.BOTH, drop=False)
   down = labeled.where(labeled["cast_type"] == CastType.DOWNCAST.value,
   drop=True)
   up = labeled.where(labeled["cast_type"] == CastType.UPCAST.value,
   drop=True)

Deprecated Parameter Names

Several processing functions (`cell_thermal_mass`,
`bouyancy_frequency`, `buoyancy`) now use:

*   `temperature`, `salinity`, `pressure`

Old names are accepted but warn.

.. code-block:: python

   # Before
   ctm = cell_thermal_mass(temperature_C=t, conductivity_Sm=c,
   amplitude=0.03, time_constant=7, sample_interval=0.25)

   # After
   ctm = cell_thermal_mass(temperature=t, conductivity=c,
   amplitude=0.03, time_constant=7, sample_interval=0.25)


`utils` module
--------------

*   `get_tolerance` now accepts `flag_value` explicitly.
*   New decorator `profile` for line profiling.
*   New `WarnAllMembersMeta` metaclass (used for deprecation warnings).


End-to-End Migration Example
============================

.. code-block:: python

   from seabirdscientific.instrument_data import read_hex_file, Sensors, InstrumentType
   from seabirdscientific.processing import bin_average, split, CastType

   ds = read_hex_file(
   filepath="cast.hex",
   instrument_type=InstrumentType.SBE37SM,
   enabled_sensors=[Sensors.Temperature, Sensors.Conductivity, Sensors.Pressure],
   moored_mode=False,
   )

   binned = bin_average(ds, bin_variable="prdM", bin_size=1.0)

   labeled = split(binned, "prdM", cast_type=CastType.BOTH,
   exclude_bad_scans=True, drop=False)

   down = labeled.where(labeled["cast_type"] == CastType.DOWNCAST.value,
   drop=True)
