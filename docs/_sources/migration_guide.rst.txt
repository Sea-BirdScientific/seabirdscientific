.. _migration_guide:

Migration Guide
###############

v2.7.3 → v3.0.0
****************

Breaking Changes Summary
========================

*   **Pandas → Xarray**:

    *   `read_hex_file(...)` now returns :class:`xarray.Dataset` (previously :class:`pandas.DataFrame`).
    *   `processing.bin_average(...)` now *accepts* and *returns* :class:`xarray.Dataset`.
    *   New `instrument_data.read_cnv_file(...) -> xarray.Dataset` provides a direct CNV→xarray path, replacing the old custom container approach.

*   **Snake-case reader names + deprecation shims**:

    *   `read_SBE19plus_format_0` → `read_sbe19plus_format_0`
    *   `read_SBE37SM_format_0` → `read_sbe37sm_format_0`
    *   `read_SBE39plus_format_0` → `read_sbe39plus_format_0`
    *   `read_SBE911plus_format_0` → `read_sbe911plus_format_0`  
        Old names remain as *wrappers* that emit `DeprecationWarning`.

*   **Enums at call sites deprecated → use string literals**:

    *   For processing functions, pass strings (e.g., `"fixed"`, `"percent"`, `"boxcar"`, `"median"`) instead of `MinVelocityType` or `WindowFilterType`. Using the Enums still works but warns.

*   **HEX metadata refactor**:

    *   Prefer new `HEX_TYPE_*` and `HEX_LEN_*` constants for type names and field lengths.
    *   `HexDataTypes` Enum is retained but emits warnings on member access.

*   **Parameter cleanup / deprecations removed**:

    *   SeaFET pH converters removed `ph_counts`; use `raw_ph` with `ph_units="counts"|"volts"`.
    *   SeaFET coefficient constructors dropped deprecated aliases (`k0`, `k2`, `int_k0`, `int_k2`, `ext_k0`, `ext_k2`).
    *   Legacy dummy `hex=...` parameters removed from low-level readers.

*   **New EOS-80 module and deprecation routing**:

    *   New `eos80_conversion` module contains EOS-80 implementations (`density`, `potential_temperature`, `adiabatic_temperature_gradient`, and EOS-80 `bouyancy_frequency`).
    *   `eos80_processing` now delegates to `eos80_conversion` and marks old functions as deprecated.

*   **Loop-edit API clarified**:

    *   New canonical `processing.loop_edit(...)` operates on *depth* or *pressure* (with `units`/`latitude`), returns updated flag array.
    *   `loop_edit_depth(...)` and `loop_edit_pressure(...)` retained as **deprecated** compatibility wrappers.

*   **Cast splitting now labels data**:

    *   `processing.split(...)` *adds* a `cast_type` coordinate to a dataset (`"downcast"`, `"upcast"`, or empty string) and can optionally `drop=True` to subset. It no longer returns a list of frames by default.


Module-level API Changes
========================

`cal_coefficients` module
-------------------------

Removed deprecated constructor parameters:

*   `PHSeaFETInternalCoefficients`: no longer accepts `k0`, `k2`, `int_k0`, `int_k2`.
*   `PHSeaFETExternalCoefficients`: no longer accepts `ext_k0`, `ext_k2`.

**Migration Example**

.. code-block:: python

   # Before
   ci = PHSeaFETInternalCoefficients(k0=1.0, k2=2.0)

   # After
   ci = PHSeaFETInternalCoefficients(kdf0=1.0, kdf2=2.0)


.. code-block:: python

   # Before
   ce = PHSeaFETExternalCoefficients(ext_k0=1.0, ext_k2=2.0)

   # After
   ce = PHSeaFETExternalCoefficients(k0=1.0, k2=2.0)


`conversion` module
-------------------

Parameter removals and doc clarifications:

*   `convert_internal_seafet_ph(...)`: removed `ph_counts`; use `raw_ph` with `ph_units`.
*   `convert_external_seafet_ph(...)`: removed `ph_counts`.

**Migration Example**

.. code-block:: python

   # Before
   ph = convert_internal_seafet_ph(ph_counts=counts)

   # After
   ph = convert_internal_seafet_ph(raw_ph=counts, ph_units="counts")


.. code-block:: python

   # Before
   ph = convert_external_seafet_ph(ph_counts=counts)

   # After
   ph = convert_external_seafet_ph(raw_ph=counts, ph_units="counts")


New buoyancy utilities added here (moved from `processing`):

*   `buoyancy_frequency(temperature, salinity, pressure, gravity)`
*   `buoyancy(temperature, salinity, pressure, latitude, longitude, window_size, ...)`

`eos80_conversion` module
-------------------------

EOS-80 implementations moved from `eos80_procesing`:

*   `density(salinity, temperature, pressure)`
*   `potential_temperature(salinity, temperature, pressure, mean_pressure)`
*   `adiabatic_temperature_gradient(salinity, temperature, pressure)`
*   `bouyancy_frequency(temperature, salinity, pressure, gravity)`


`eos80_processing` module
-------------------------

Functions now delegate to `eos80_conversion` equivalents and emit deprecation warnings:

*   `bouyancy_frequency(...)` → calls `eos80_conversion.bouyancy_frequency(...)`
*   `density(...)` → `eos80_conversion.density(...)`
*   `potential_temperature(...)` → `eos80_conversion.potential_temperature(...)`
*   `adiabatic_temperature_gradient(...)` → `eos80_conversion.adiabatic_temperature_gradient(...)`

**Migration Example**

.. code-block:: python

   # Before
   from seabirdscientific import eos80_processing as ep
   n2 = ep.bouyancy_frequency(
       temp_ITS_subset=t,
       salinity_prac_subset=s,
       pressure_dbar_subset=p,
       gravity=g
   )

   # After
   from seabirdscientific import eos80_conversion as ec
   n2 = ec.bouyancy_frequency(temperature=t, salinity=s, pressure=p, gravity=g)


`instrument_data` module
------------------------


CNV Reading Now Returns Xarray

*   `read_cnv_file(filepath) -> xarray.Dataset`  
    Replaces `InstrumentData` / `MeasurementSeries` / `cnv_to_instrument_data`.

**Migration Example**

.. code-block:: python

   # Before
   inst = cnv_to_instrument_data("file.cnv")
   dataframe = inst._to_dataframe()

   # After
   dataset = read_cnv_file("file.cnv")
   dataframe = ds.to_dataframe()

HEX Reading and Decoding

*   `read_hex_file` now returns `xarray.Dataset` instead of `pandas.DataFrame`.
*   Added `_preallocate_dataset` helper.
*   Introduced `HEX_TYPE_*` and `HEX_LEN_*` constants.
*   `HexDataTypes` still present but every member access triggers a `DeprecationWarning`.

**Migration Example**

.. code-block:: python

   # Before
   dataframe = read_hex_file("cast.hex", itype, enabled_sensors=[...])

   # After
   dataset = read_hex_file("cast.hex", itype, enabled_sensors=[...])
   dataframe = dataset.to_dataframe()

Snake-case Format Readers

Old → new names:

*   `read_SBE19plus_format_0` → `read_sbe19plus_format_0`
*   `read_SBE37SM_format_0` → `read_sbe37sm_format_0`
*   `read_SBE39plus_format_0` → `read_sbe39plus_format_0`
*   `read_SBE911plus_format_0` → `read_sbe911plus_format_0`

`processing` module
-------------------

Loop Edit

*   `loop_edit` replaces `loop_edit_depth` and `loop_edit_pressure`, which remain only as compatibility shims and emit deprecation warnings. Use `units="depth"` (default) or `units="pressure"`
*   The `min_velocity_type` argument now expects strings (e.g., `"fixed"`, `"percent"`) instead of the `MinVelocityType` enum. Passing the enum is accepted but deprecated.

**Migration Example (prefer string literal over enum and use the core API)**

.. code-block:: python

      # Before
      m = loop_edit_depth(depth, flag, dt, MinVelocityType.FIXED, 0.1)

      # After
      m = loop_edit_depth(depth, flag, dt, "fixed", 0.1)


Window Filter

*   `window_filter(...)` now expects a string literal for the window type: `"boxcar"`, `"cosine"`, `"gaussian"`, `"median"`, or `"triangle"`.
*   Passing `WindowFilterType` still works but warns.

**Migration Example**

.. code-block:: python

  # Before
  y = window_filter(x, flags, WindowFilterType.MEDIAN, 5, dt)

  # After
  y = window_filter(x, flags, "median", 5, dt)


Xarray bin-average

- ``bin_average(dataset, ...)`` now takes and return an :class:`xarray.Dataset`.
- Output uses a ``bin_number`` dimension; dataset attributes are preserved.
- Supports options like ``include_scan_count``, ``interpolate``, cast selection (via ``cast_type`` alongside ``split(...)``), and flag handling.

**Migration Example**

.. code-block:: python

   dataset = read_hex_file(...)
   binned_dataset = bin_average(dataset, bin_variable="prdM", bin_size=1.0)

Cast splitting via coordinate

*   `split(...)` now labels the data by adding a `cast_type` coordinate (values: `"downcast"`, `"upcast"`, or `""`).
*   Use `drop=True` to immediately subset.

**Migration Example**

.. code-block:: python

   labeled = split(
       binned_dataset,
       "prdM",
       cast_type=CastType.BOTH,
       exclude_bad_scans=True,
       drop=False
   )
   down = labeled.where(labeled["cast_type"] == CastType.DOWNCAST.value, drop=True)
   up = labeled.where(labeled["cast_type"] == CastType.UPCAST.value, drop=True)

TEOS-10 buoyancy relocation

- ``processing.bouyancy_frequency(...)`` and ``processing.buoyancy(...)`` are deprecated wrappers that now call the new implementations in ``conversion``. Prefer calling ``conversion`` directly.

**Migration Example**

.. code-block:: python

   from seabirdscientific.conversion import buoyancy_frequency
   n2 = buoyancy_frequency(temperature=ct, salinity=sa, pressure=p, gravity=g)

Enum value adjustments

- ``CastType`` values are now strings: ``BOTH="both"``, ``DOWNCAST="downcast"``, ``UPCAST="upcast"``, ``NONE=""``.
- Logic based on equality remains the same; avoid numeric comparisons.

``utils`` module
----------------

- ``get_tolerance(data, flag_value=-9.99e-29)`` now accepts an explicit ``flag_value`` argument for tests.
- New: ``profile`` decorator (uses ``line_profiler``) to print per-line timings during development.
- New: ``WarnAllMembersMeta`` metaclass used to warn when deprecated Enum members are accessed (e.g., ``HexDataTypes.temperature``).


End-to-End Examples
===================

HEX → bin → cast labeling (xarray)
----------------------------------

.. code-block:: python

   from seabirdscientific.instrument_data import read_hex_file, Sensors, InstrumentType
   from seabirdscientific.processing import bin_average, split, CastType

   dataset = read_hex_file(
       filepath="cast.hex",
       instrument_type=InstrumentType.SBE37SM,
       enabled_sensors=[Sensors.Temperature, Sensors.Conductivity, Sensors.Pressure],
       moored_mode=False,
   )

   binned_dataset = bin_average(dataset, bin_variable="prdM", bin_size=1.0)

   labeled = split(
       binned_dataset,
       "prdM",
       cast_type=CastType.BOTH,
       exclude_bad_scans=True,
       drop=False
   )

   down = labeled.where(
       labeled["cast_type"] == CastType.DOWNCAST.value,
       drop=True
   )
