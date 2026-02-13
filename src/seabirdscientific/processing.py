"""Functions for processing converted SBS instrument data."""

# Native imports
import warnings
from enum import Enum
from logging import getLogger
from typing import Literal

# Third-party imports
import numpy as np
import pandas as pd
import xarray as xr
from scipy import signal

# Sea-Bird imports
from seabirdscientific import conversion as c


logger = getLogger(__name__)


FLAG_VALUE = -9.99e-29


class MinVelocityType(Enum):
    """The minimum velocity type used with loop edit
    DEPRECATED. Use Literal defined in loop_edit"""

    FIXED = 0
    PERCENT = 1


class WindowFilterType(Enum):
    """The window filter type. See CDT Processing in the docs for
    details.
    DEPRECATED. Use literals defined in windowd_filter
    """

    BOXCAR = "boxcar"
    COSINE = "cosine"
    TRIANGLE = "triangle"
    GAUSSIAN = "gaussian"
    MEDIAN = "median"


class CastType(Enum):
    """The subsection of data to use when splitting by upcast and/or
    downcast
    """

    BOTH = "both"
    DOWNCAST = "downcast"
    UPCAST = "upcast"
    NONE = ""


def butterworth_filter(x: np.ndarray, time_constant: float, sample_interval: float) -> np.ndarray:
    """Applies a butterworth filter to a dataset.

    :param x: A numpy array of floats
    :param time_constant: 1 / (2 * pi * cutoff_frequency)
    :param sample_interval: 1 / sampling_frequency

    :return: the filtered data
    """

    cutoff_frequency = 1 / (2 * np.pi * time_constant)
    sampling_frequency = 1 / sample_interval
    b, a = signal.butter(1, cutoff_frequency, fs=sampling_frequency)
    y = signal.filtfilt(b, a, x)
    return y


def low_pass_filter(x: np.ndarray, time_constant: float, sample_interval: float) -> np.ndarray:
    """Applies a low pass filter as defined in the SeaSoft manual
    v7.26.8 page 101.

    :param x: A numpy array of floats
    :param time_constant: 1 / (2 * pi * cutoff_frequency)
    :param sample_interval: 1 / sampling_frequency

    :return: the filtered data
    """

    a = 1 / (1 + 2 * time_constant / sample_interval)
    b = a * (1 - 2 * time_constant / sample_interval)

    x2 = x.copy()
    for n in range(1, len(x)):
        x2[n] = a * x[n] + a * x[n - 1] - b * x2[n - 1]

    y = x2.copy()
    for n in range(len(x) - 2, -1, -1):
        y[n] = a * x2[n] + a * x2[n + 1] - b * y[n + 1]

    return y


def align_ctd(
    x: np.ndarray, offset: float, sample_interval: float, flag_value=FLAG_VALUE
) -> np.ndarray:
    """Takes an ndarray object of data for a single variable and applies
    a time offset to the series.

    Performs interpolation when offset is not a factor of sampleinterval

    :param x: array of values to apply shift to
    :param offset: offset value to shift by (s)
    :param sample_interval: time between samples (s)

    :return: the aligned data
    """

    sample_times = np.arange(0, len(x) * sample_interval, sample_interval)
    interp_res = np.interp(sample_times + offset, sample_times, x, flag_value, flag_value)
    return interp_res


def cell_thermal_mass(
    temperature: np.ndarray = [],
    conductivity: np.ndarray = [],
    amplitude: float = 1,
    time_constant: float = 1,
    sample_interval: float = 1,
    temperature_C: np.ndarray = None,  # Deprecated
    conductivity_Sm: np.ndarray = None,  # Deprecated
) -> np.ndarray:
    """Removes conductivity cell thermal mass effects from measured
    conductivity.

    Cell Thermal Mass uses a recursive filter to remove conductivity
    cell thermal mass effects from the measured conductivity [From the
    SeaSoft manual, page 92]

    :param temperature_C: temperature in degrees C
    :param conductivity_Sm: conductivity in S/m
    :param amplitude: thermal anomaly amplitude (alpha)
    :param time_constant: thermal anomoly time constant (1/beta)
    :param sample_interval: time between samples

    :return: the corrected conductivity in S/m
    """

    if temperature_C is not None:
        warnings.warn("Deprecated, use temperature", DeprecationWarning)
        temperature = temperature_C

    if conductivity_Sm is not None:
        warnings.warn("Deprecated, use salinity", DeprecationWarning)
        conductivity = conductivity_Sm

    a = 2 * amplitude / (sample_interval / time_constant + 2)
    b = 1 - (2 * a / amplitude)
    ctm = np.zeros(len(temperature))  # cell thermal mass
    corrected_conductivity = conductivity.copy()

    for n in range(1, len(ctm)):
        dc_dt = 0.1 * (1.0 + 0.006 * (temperature[n] - 20.0))
        dt = temperature[n] - temperature[n - 1]
        ctm[n] = -1.0 * b * ctm[n - 1] + a * dc_dt * dt
        corrected_conductivity[n] += ctm[n]

    return corrected_conductivity


def loop_edit_pressure(
    pressure: np.ndarray,
    latitude: float,
    flag: np.ndarray,
    sample_interval: float,
    min_velocity_type: MinVelocityType,
    min_velocity: float,
    window_size: float,
    mean_speed_percent: float,
    remove_surface_soak: bool,
    min_soak_depth: float,
    max_soak_depth: float,
    use_deck_pressure_offset: bool,
    exclude_flags: bool,
    flag_value=FLAG_VALUE,
) -> np.ndarray:
    """Deprecated. Use loop_edit"""

    warnings.warn("Deprecated. Use loop_edit", DeprecationWarning)

    depth = c.depth_from_pressure(pressure, latitude)
    return loop_edit_depth(
        depth,
        flag,
        sample_interval,
        min_velocity_type,
        min_velocity,
        window_size,
        mean_speed_percent,
        remove_surface_soak,
        min_soak_depth,
        max_soak_depth,
        use_deck_pressure_offset,
        exclude_flags,
        flag_value,
    )


def loop_edit_depth(
    depth: np.ndarray,
    flag: np.ndarray,
    sample_interval: float,
    min_velocity_type: Literal["fixed", "percent"],
    min_velocity: float,
    window_size: float,
    mean_speed_percent: float,
    remove_surface_soak: bool,
    min_soak_depth: float,
    max_soak_depth: float,
    use_deck_pressure_offset: bool,
    exclude_flags: bool,
    flag_value=FLAG_VALUE,
) -> np.ndarray:
    """Deprecated. Use loop_edit"""

    warnings.warn("Deprecated. Use loop_edit", DeprecationWarning)

    return loop_edit(
        depth,
        flag,
        sample_interval,
        min_velocity_type,
        min_velocity,
        window_size,
        mean_speed_percent,
        remove_surface_soak,
        min_soak_depth,
        max_soak_depth,
        use_deck_pressure_offset,
        exclude_flags,
        flag_value,
        units="depth"
    )


def loop_edit(
    measurand: np.ndarray,
    flag: np.ndarray,
    sample_interval: float,
    min_velocity_type: Literal["fixed", "percent"] = "fixed",
    min_velocity: float = 0.25,
    window_size: float = 300,
    mean_speed_percent: float = 20,
    remove_surface_soak: bool = True,
    min_soak_depth: float = 5,
    max_soak_depth: float = 20,
    use_deck_pressure_offset: bool = True,
    exclude_flags: bool = True,
    flag_value=FLAG_VALUE,
    latitude: float = 0,
    units: Literal["depth", "pressure"] = "depth"
) -> np.ndarray:
    """Marks scans determined to be part of a pressure loop as bad.

    Loop Edit marks scans bad by setting the flag value associated with
    the scan to badflag in data that has pressure slowdowns or reversals
    (typically caused by ship heave).

    :param measurand: Salt water depth or pressure as an array of
        positive numbers
    :param flag: Array of flag values. The flag for a good value is 0.0.
        The flag for a detected loop value defaults to -9.99e-29
    :param sample_interval: Time interval between samples in seconds
    :param min_velocity_type: Sets whether flags are based on min
        velocity or a percentage of mean speed
    :param min_velocity: The minimum velocity for data to be considered
        good
    :param window_size: Time interval to include in mean speed
        calculation
    :param mean_speed_percent: Percentage of mean speed for data to be
        considered good
    :param remove_surface_soak: If true, data that occur before the
        minimum soak depth is reached are marked as bad
    :param min_soak_depth: The depth that must be reached before data
        can be considered good
    :param max_soak_depth: The maximum depth that can be considered the
        start of a downcast
    :param use_deck_pressure_offset: If true, min and max soak depths
        are offset by the first value in the depth array
    :param exclude_flags: If true, existing bad flags are preserved
    :param flag_value: Passing is 0.0, failing defaults to -9.99e-29.

    Returns:
        np.ndarray: the input data with updated flags
    """

    if isinstance(min_velocity_type, MinVelocityType):
        warnings.warn("MinVelocityType Enum is deprecated, use Literals", DeprecationWarning)

    if units == "pressure":
        depth = c.depth_from_pressure(measurand, latitude)
    else:
        depth = measurand.copy()
    
    _flag = flag.copy()

    if not exclude_flags:
        _flag[:] = 0.0

    if use_deck_pressure_offset:
        min_soak_depth -= depth[0]
        max_soak_depth -= depth[0]

    (min_depth_n, max_depth_n) = _find_depth_peaks(
        depth, _flag, remove_surface_soak, flag_value, min_soak_depth, max_soak_depth
    )

    if min_velocity_type in ["fixed", MinVelocityType.FIXED]:
        downcast_mask = _min_velocity_mask(
            depth, sample_interval, min_velocity, min_depth_n, max_depth_n + 1, False
        )
        upcast_mask = _min_velocity_mask(
            depth, sample_interval, min_velocity, max_depth_n, len(depth), True
        )

    elif min_velocity_type in ["percent", MinVelocityType.PERCENT]:
        diff_length = int(window_size / sample_interval)
        downcast_mask = _mean_speed_percent_mask(
            depth,
            sample_interval,
            min_velocity,
            mean_speed_percent,
            min_depth_n,
            max_depth_n,
            False,
            diff_length,
        )
        upcast_mask = _mean_speed_percent_mask(
            depth,
            sample_interval,
            min_velocity,
            mean_speed_percent,
            max_depth_n,
            len(depth),
            True,
            diff_length,
        )
    else:
        raise ValueError

    _flag[~downcast_mask & ~upcast_mask] = flag_value

    _flag_by_minima_maxima(depth, _flag, min_depth_n, max_depth_n, flag_value)

    return _flag


def _find_depth_peaks(
    depth: np.ndarray,
    flag: np.ndarray,
    remove_surface_soak: bool,
    flag_value: float,
    min_soak_depth: float,
    max_soak_depth: float,
) -> tuple[int, int]:
    """Finds the global depth minima and maxima.

    This determinines the earliest points where the downcast and upcast
    can begin

    :param depth: depth data
    :param flag: flag data
    :param remove_surface_soak: If true, scans before a minimum depth
        are mrked as bad
    :param flag_value: the flag value (typically -9.99e-29)
    :param min_soak_depth: the minimum depth that must be reached before
        a series can be considered a downcast
    :param max_soak_depth: maximumm depth that can be considered the
        start of a downcast

    :return: minimum and maximum index corresponding to minimum and
        maximum depth
    """

    # first index where min_soak_depth < depth < max_soak_depth
    if remove_surface_soak:
        min_soak_depth_n = min(
            n
            for n, d in enumerate(depth)
            if flag[n] != flag_value and min_soak_depth < d < max_soak_depth and d != flag_value
        )
    else:
        min_soak_depth_n = 0

    # beginning of possible upcast domain
    max_depth = max([d for n, d in enumerate(depth) if flag_value not in [d, flag[n]]])
    max_depth_n = np.where(depth == max_depth)[0][0]

    # beginning of possible downcast domain
    min_depth = min(
        [
            d
            for n, d in enumerate(depth)
            if flag_value not in [d, flag[n]] and min_soak_depth_n < n < max_depth_n
        ]
    )
    min_depth_n = np.where(depth == min_depth)[0][0]

    return (min_depth_n, max_depth_n)


def _min_velocity_mask(
    depth: np.ndarray,
    interval: float,
    min_velocity: float,
    domain_start: int,
    domain_end: int,
    is_upcast: bool,
) -> np.ndarray:
    """Creates a mask to assign bad flags where velocity is less that
    the provided minimum.

    :param depth: depth data
    :param interval: sampling interval
    :param min_velocity: minimum velocity (only used for samples within
        first window)
    :param domain_start: earliest possible beginning of downcast or
        upcast
    :param domain_end: earliest possible end of downcast or upcast
    :param is_upcast: inverts velocity sign for upcast

    :return: true/false mask where true means velocity is at least the
        minimum
    """

    sign = -1 if is_upcast else 1
    mask0 = sign * np.diff(depth, prepend=depth[0]) / interval >= min_velocity
    mask1 = sign * np.diff(depth, append=depth[-1]) / interval >= min_velocity

    mask = mask0 & mask1

    mask[:domain_start] = False
    mask[domain_end:] = False

    return mask


def _mean_speed_percent_mask(
    depth: np.ndarray,
    interval: float,
    min_velocity: float,
    mean_speed_percent: float,
    domain_start: int,
    domain_end: int,
    is_upcast: bool,
    diff_length: int,
) -> np.ndarray:
    """Assigns bad flags to scans below a specified velocity.

    This creates a mask to assign bad flags where the velocity is less
    than a given percentage of the mean accross a number of samples
    according to diff_length (determined by the window and sample_rate)

    :param depth: depth data
    :param interval: sampling interval in seconds
    :param min_velocity: minimum velocity
    :param mean_speed_percent: the minimum percentage of the mean speed
        that qualifies a sample as good
    :param domain_start: earliest possible beginning of downcast
    :param domain_end: earliest possible beginning of upcast
    :param is_upcast: inverts velocity sign for upcast
    :param diff_length: averaging window divided by sample interval

    :return: true/false mask where true means the sample velocity is at
        least the given percentage of the mean
    """

    sign = -1 if is_upcast else 1
    mask0 = sign * np.diff(depth[0:diff_length]) / interval > min_velocity
    mask1 = sign * np.diff(depth[1 : diff_length + 1]) / interval > min_velocity
    first_window_mask = np.concatenate([[False], (mask0 & mask1)])

    mean_speed = sign * (depth[diff_length:] - depth[:-diff_length]) / diff_length / interval
    speed = sign * np.diff(depth[diff_length:], prepend=depth[diff_length - 1]) / interval

    mask = np.concatenate((first_window_mask, (speed > (mean_speed * mean_speed_percent / 100.0))))
    mask[:domain_start] = False
    mask[domain_end:] = False

    return mask


def _flag_by_minima_maxima(
    depth: np.ndarray,
    flag: np.ndarray,
    min_depth_n: int,
    max_depth_n: int,
    flag_value: float,
):
    """Flags data that is less than the most recent local minima or
    maxima.

    This condition occurs following a ship heave, where velocity
    temporarily inverts

    :param depth: depth data
    :param flag: flag data.  Will be modified by this function.
    :param min_depth_n: index of global minimum depth
    :param max_depth_n: index of global maximum depth
    :param flag_value: value assigned to bad flags
    """

    local_max = -10000.0
    local_min = 10000.0

    # flag values that don't exceed most recent valid local minima/maxima
    for n, d in enumerate(depth):
        if n >= max_depth_n and d < local_min and flag[n] != flag_value:
            local_min = d
        elif n >= min_depth_n and d > local_max and flag[n] != flag_value:
            local_max = d
        else:
            flag[n] = flag_value


def bin_average(
    dataset: xr.Dataset,
    bin_variable: str,
    bin_size: float,
    include_scan_count: bool = True,
    min_scans: int = 1,
    max_scans: int = 999999,
    exclude_bad_scans: bool = True,
    interpolate: bool = False,
    cast_type: CastType = CastType.BOTH,
    trim_start: int = 0,
    trim_end: int = 0,
    include_surface_bin: bool = False,
    surface_bin_min: float = 0,
    surface_bin_max: float = 5,
    surface_bin_value: float = 2.5,
    flag_value=FLAG_VALUE,
) -> xr.Dataset:
    """Averages data into bins, using intervals based on bin_variable.
    Returns a new dataframe with the binned data

    :param dataset: Dataset containing all data to group into bins
    :param bin_variable: The variable that will control binning,
        typically pressure, depth, time, or scan number. For scan number,
        use 'nScan'.
    :param bin_size: The bin width or range of data for each bin
    :param include_scan_count: If True includes a column (nbin) in the
        returned dataframe for the number of scans in each bin. Defaults
        to True
    :param min_scans: The minimum number of scans required in a bin for
        it to be included in the dataset. Defaults to 1
    :param max_scans: the maximum number of scans allowed to for a bin
        to be included in the dataset. Defaults to 999999
    :param exclude_bad_scans: If True, removes scans marked bad by loop
        edit (scans marked bad by wild edit are always excluded).
        Defaults to True
    :param interpolate: If True interpolates bins after averaging.
        Defaults to False
    :param cast_type: Sets which data to include. When binning by depth
        or pressure use, UPCAST, DOWNCAST, or BOTH. When binning by time
        or other variables use NA. Defaults to CastType.BOTH
    :param trim_start: Remove this number of scans from the beginning
        of the initial dataset. Defaults to 0
    :param trim_end: Remove this number of scans fro mteh end of the
        initial dataset. Defaults to 0
    :param include_surface_bin: Includes a surface bin at the beginning
        of the downcast and/or at the end of the upcast. Defaults to
        False
    :param surface_bin_min: The minimum value of the bin_variable to
        include in the surface bin. Defaults to 0 (value less than 0
        will be set to 0)
    :param surface_bin_max: The maximum value of the bin_variable to
        include in the surface bin. Defaults to 5 (this will be
        constrained to surface_bin_min and bin_size / 2)
    :param surface_bin_value: The target value for interpolating the
        surface bin at the beginning of the downcast. Defaults to 2.5
        (the target value for the upcast surface bin is always 0)
    :param flag_value: The magical number indicating bad data.
        Defaults to -9.99e-29
    :return: A new Dataset with binned data
    """

    df = dataset.to_dataframe().iloc[trim_start : len(dataset["sample"]) - trim_end]

    # remove scans marked as bad during loop edit
    if exclude_bad_scans and "flag" in df.columns:
        df = df.drop(df[df["flag"] == flag_value].index)

    # always remove scans marked bad during wild edit
    for column in df.columns.difference(["flag"]):
        df = df.drop(df[df[column] == flag_value].index)

    if bin_variable == "nScan" and "nScan" not in df.columns:
        # We want to bin by scans, ensure it's a column
        df.insert(0, "nScan", range(0, len(df)))

    # pd series containing the variable we want to bin for, converted to ndarray
    control = df[bin_variable].to_numpy()

    bin_min = bin_size / 2.0  # min value of first bin

    if bin_variable == "nScan":
        # when binning by scan number, start at 0
        bin_min = 0

    control_max = np.amax(control)
    bin_max = control_max - ((bin_min + control_max) % bin_size) + bin_size

    # split into descending and ascending, including peak in both
    peak_index = np.argmax(control)
    control_desc = control[: peak_index + 1]
    control_asc = control[peak_index:]

    # create the bins to sort into
    desc_bin_edges = np.arange(start=bin_min, stop=bin_max + bin_size, step=bin_size)
    asc_bin_adges = np.arange(start=bin_max, stop=bin_min - bin_size, step=-bin_size)

    # setup bins to indicate where each index should be sorted into
    desc_bins = np.digitize(x=control_desc, bins=desc_bin_edges)
    asc_bins = np.digitize(x=control_asc, bins=asc_bin_adges, right=True)
    asc_bins += np.amax(desc_bins) - 1
    df["bin_number"] = np.concat((desc_bins[:-1], asc_bins))

    if interpolate:
        desc_midpoints = (desc_bin_edges[:-1] + desc_bin_edges[1:]) / 2
        asc_midpoints = (asc_bin_adges[1:-1] + asc_bin_adges[2:]) / 2
        midpoints = np.concat((desc_midpoints, asc_midpoints))
        df["midpoint"] = df["bin_number"].map(
            lambda n: midpoints[n - 1] if 0 < n <= len(midpoints) else np.nan
        )

    if include_surface_bin:
        if surface_bin_min < 0:
            logger.warning("Surface bin min set to 0")
        _surface_bin_min = max(surface_bin_min, 0)

        if not surface_bin_min <= surface_bin_max <= bin_min:
            logger.warning(f"Surface bin max set to {surface_bin_max}")
        _surface_bin_max = max(surface_bin_min, min(surface_bin_max, bin_min))

        bins = (_surface_bin_min, _surface_bin_max)
        surface_desc_bin = np.digitize(x=control_desc, bins=bins) == 1
        surface_asc_bin = np.digitize(x=control_asc, bins=bins, right=True) == 1

        # these will get added back in depending on the cast type
        surface_desc = df[: len(surface_desc_bin)][surface_desc_bin]
        surface_asc = df[len(surface_desc_bin) - 1 :][surface_asc_bin]

        if interpolate:
            surface_desc["midpoint"] = surface_desc["midpoint"].fillna(surface_bin_value)
            surface_asc["midpoint"] = surface_asc["midpoint"].fillna(0)

    if cast_type == CastType.DOWNCAST:
        # keeping one past the peak index to match SBE data processing
        df = df[: peak_index + 2]
        min_bin_number = 1
        below_min = df["bin_number"] < min_bin_number
        df = df.drop(df[below_min].index)

        if include_surface_bin:
            df = pd.concat((surface_desc, df))

    elif cast_type == CastType.UPCAST:
        df = df[peak_index:]
        # discarding the first bin to match SBE data processing
        min_bin_number = np.amin(asc_bins) + 1
        max_bin_number = np.amax(asc_bins) - 1
        below_min = df["bin_number"] < min_bin_number
        over_max = df["bin_number"] > max_bin_number
        df = df.drop(df[below_min | over_max].index)

        if include_surface_bin:
            df = pd.concat((df, surface_asc))

    elif cast_type == CastType.BOTH:
        # drop first and last these since they're not necessarily the same as surface bin
        min_bin_number = 1
        max_bin_number = np.amax(asc_bins) - 1
        below_min = df["bin_number"] < min_bin_number
        over_max = df["bin_number"] > max_bin_number
        df = df.drop(df[below_min | over_max].index)

        if include_surface_bin:
            df = pd.concat((surface_desc, df, surface_asc))

    # else cast_type == CastType.NA:
    # do nothing

    # get the number of scans in each bin
    scans_per_bin = np.bincount(df["bin_number"])
    df["nbin"] = df["bin_number"].map(lambda x: scans_per_bin[x])

    if exclude_bad_scans:
        df = df.groupby("bin_number", as_index=False).mean()
    else:
        # need to handle the flag column differently
        not_flag = df[df.columns.difference(["flag"])].groupby("bin_number").mean()
        # if all the values in a group are the flag value the assign the
        # flag value to the group, otherwise 0
        flag = df[["bin_number", "flag"]].groupby("bin_number").mean()
        flag.loc[flag["flag"] != flag_value] = 0
        df = pd.concat([not_flag, flag], axis=1).reset_index()

    df = df.drop(df[df["nbin"] < min_scans].index)
    df = df.drop(df[df["nbin"] > max_scans].index)

    if interpolate:

        def interp(p_p, x_p, p_c, x_c, p_i):
            """Interpolate according to SBE Data Processing manual
            version 7.26.8, page 89
            """
            x_i = ((x_c - x_p) * (p_i - p_p) / (p_c - p_p)) + x_p
            return x_i

        excluded_columns = ["nbin", "flag", "bin_number", bin_variable, "midpoint"]
        for column in (df.columns).difference(excluded_columns):
            interp_result = []
            for n in range(len(df[column])):
                n_p = 1 if n == 0 else n - 1
                p_p = df[bin_variable].iloc[n_p]
                x_p = df[column].iloc[n_p]
                p_c = df[bin_variable].iloc[n]
                x_c = df[column].iloc[n]
                p_i = df["midpoint"].iloc[n]
                x_i = interp(p_p, x_p, p_c, x_c, p_i)
                interp_result.append(x_i)

            df[column] = pd.Series(interp_result, index=df.index)

        df[bin_variable] = df["midpoint"]
        df = df.drop("midpoint", axis=1)

    df = df.drop("bin_number", axis=1)

    if not include_scan_count:
        df = df.drop("nbin", axis=1)

    _dataset = xr.Dataset(df, attrs=dataset.attrs).rename({"dim_0": "bin_number"})

    return _dataset


def wild_edit(
    data: np.ndarray,
    flags: np.ndarray,
    std_pass_1: float,
    std_pass_2: float,
    scans_per_block: int,
    distance_to_mean: float,
    exclude_bad_flags: bool,
    flag_value=FLAG_VALUE,
) -> np.ndarray:
    """Flags outliers in a dataset.

    Outliers are flagged by iterating over the data in blocks, taking
    the mean and standard deviation, and flagging data outside the
    combined variance (standard deviation of the block multiplied by the
    standard deviation argument). Each block is processed three times:
    first to remove loop edit flags if applicable, second to temporarily
    remove outliers outside std_pass_1, third to remove outliers outside
    std_pass_2 from the returned data.

    If the final block contains fewer samples than scans_per_block, data
    is backfilled from the previous block before computing the mean and
    standard deviation.

    This algorithm may introduce nonlinearities in the data because it
    runs one block at a time instead of a rolling window. This is done
    to maintain parity with SBE Data Processing.

    :param data: The data to be flagged, such as temperature or pressure
    :param flags: Flag data from loop edit
    :param std_pass_1: Standard deviation for the first pass
    :param std_pass_2: Standard deviation for the second pass
    :param scans_per_block: Number of samples to process in each block
    :param distance_to_mean: Minimum threshod to flag data. Values
        within this range to mean are kept even if their standard
        deviation is outside the limit
    :param exclude_bad_flags: Excludes all loop edit flags
    :param flag_value: The flag value written in place
        of data. Defaults to -9.99e-29.

    :return: The data with flag values in place of outliers
    """

    lower_index = 0
    upper_index = scans_per_block
    flagged_data = data.copy()

    while upper_index <= len(data):
        flagged_data[lower_index:upper_index] = _flag_data(
            data[lower_index:upper_index],
            flags,
            std_pass_1,
            std_pass_2,
            distance_to_mean,
            exclude_bad_flags,
            flag_value,
        )
        lower_index += scans_per_block
        upper_index += scans_per_block

    if lower_index < len(data):
        lower_index = len(data) - scans_per_block
        upper_index = len(data)
        flagged_data[len(data) - len(data) % scans_per_block :] = _flag_data(
            data[lower_index:upper_index],
            flags,
            std_pass_1,
            std_pass_2,
            distance_to_mean,
            exclude_bad_flags,
            flag_value,
        )[scans_per_block - len(data) % scans_per_block :]

    return flagged_data


def _flag_data(
    data: np.ndarray,
    flags: np.ndarray,
    std_pass_1: float,
    std_pass_2: float,
    distance_to_mean: float,
    exclude_bad_flags: bool,
    flag_value=FLAG_VALUE,
) -> np.ndarray:
    """Helper function for wild_edit() that handles the three main loops

    :param data: The data to be flagged, such as temperature or pressure
    :param flags: Flag data from loop edit
    :param std_pass_1: Standard deviation for the first pass
    :param std_pass_2: Standard deviation for the second pass
    :param distance_to_mean: Minimum threshod to flag data. Values
        within this range to mean are kept even if their standard
        deviation is outside the limit
    :param exclude_bad_flags: Excludes all loop edit flags
    :param flag_value: The flag value written in place of data. Defaults
        to -9.99e-29.

    :return: The data with flag values in place of outliers
    """

    data_copy = pd.Series(data.copy())
    flagged_data = data.copy()

    for n, value in enumerate(data_copy):
        if exclude_bad_flags and flags[n] == flag_value:
            data_copy[n] = np.nan

    mean = data_copy.mean()
    std = data_copy.std()

    for n, value in enumerate(data_copy):
        if abs(value - mean) >= std * std_pass_1:
            data_copy[n] = np.nan

    mean = data_copy.mean()
    std = data_copy.std()

    for n, value in enumerate(flagged_data):
        if abs(value - mean) > std * std_pass_2 and abs(value - mean) > distance_to_mean:
            flagged_data[n] = flag_value

    return flagged_data


def window_filter(
    data_in: np.ndarray,
    flags: np.ndarray,
    window_type: Literal["boxcar", "cosine", "gaussian", "median", "triangle"],
    window_width: int,
    sample_interval: float,
    half_width=1,
    offset=0.0,
    exclude_flags=False,
    flag_value=FLAG_VALUE,
) -> np.ndarray:
    """Filters a dataset by convolving it with an array of weights.

    The available window filter types are boxcar, cosine, triangle,
    gaussian and median. Refer to the SeaSoft data processing manual
    version 7.26.8, page 108

    :param data_in: data to be filtered
    :param flags: flagged data defined by loop edit
    :param window_type: the filter type (boxcar, cosine, triangle, or
        gaussian)
    :param window_width: width of the window filter (must be odd)
    :param sample_interval: sample interval of the dataset. Defaults to
        1.0.
    :param half_width: width of the guassian curve. Defaults to 1.
    :param offset: shifts the center point of the gaussian. Defaults to
        0.0.
    :param exclude_flags: exclude values from the dataset that are
        flagged in flags. Also excludes corresponding weights from
        normalization. Defaults to False.
    :param flag_value: the flag value in flags. Defaults to -9.99e-29.

    :return: the convolution of data_in and the window filter
    """

    if isinstance(window_type, WindowFilterType):
        warnings.warn("WindowFilterType Enum is deprecated, use Literals", DeprecationWarning)

    # convert flags to nan for processing
    data = [d if d != flag_value else np.nan for d in data_in]

    if exclude_flags:
        data = [data[n] if f != flag_value else np.nan for n, f in enumerate(flags)]

    # for simplicity later, window indices start negative, end positive, and are centered on 0
    window_start = int(-(window_width - 1) / 2)
    window_end = int((window_width - 1) / 2 + 1)
    window = np.array([])

    # define the window filter
    if window_type in ["boxcar", WindowFilterType.BOXCAR]:
        window = signal.windows.boxcar(window_width)

    elif window_type in ["cosine", WindowFilterType.COSINE]:
        for n in range(window_start, window_end):
            window = np.append(window, np.cos((n * np.pi) / (window_width + 1)))

    elif window_type in ["triangle", WindowFilterType.TRIANGLE]:
        window = signal.windows.triang(window_width)

    elif window_type in ["gaussian", WindowFilterType.GAUSSIAN]:
        phase = offset / sample_interval
        # the manual defines scale with sample rate, but seasoft uses sample interval
        scale = np.log(2) * (2 * sample_interval / half_width) ** 2
        for n in range(window_start, window_end):
            window = np.append(window, np.exp(-((n - phase) ** 2) * scale))

    data_valid = np.nan_to_num(data)
    data_out = data_valid.copy()
    data_start = int((window_width - 1) / 2)
    data_end = len(data) - data_start

    # run a convolution, renormalizing weights as necessary
    for n in range(len(data)):
        if data_start <= n < data_end:
            value = 0

            if window_type in ["median", WindowFilterType.MEDIAN]:
                window = np.array(data[n + window_start : n + window_end])
                value = np.nanmedian(window)
            else:
                window_valid = window.copy()

                # exclude weights from normalization if they correspond to nans in the data
                for n2 in range(len(window)):
                    if np.isnan(data[n + data_start - n2]):
                        window_valid[n2] = 0

                window_valid /= window_valid.sum()

                for n2 in range(len(window)):
                    value += window_valid[n2] * data_valid[n + data_start - n2]

            data_out[n] = value

    return data_out


def bouyancy_frequency(
    temp_conservative_subset: np.ndarray = None,
    salinity_abs_subset: np.ndarray = None,
    pressure_dbar_subset: np.ndarray = None,
    gravity: float = 0,
):
    """Deprecated. Use conversion.buoyancy_frequency"""

    warnings.warn("Deprecated. Use conversion.buoyancy_frequency", DeprecationWarning)

    return c.buoyancy_frequency(
        temp_conservative_subset, salinity_abs_subset, pressure_dbar_subset, gravity
    )


def buoyancy(
    temperature_c: np.ndarray,
    salinity_prac: np.ndarray,
    pressure_dbar: np.ndarray,
    latitude: np.ndarray = 0,
    longitude: np.ndarray = 0,
    window_size: float = 11,
    use_modern_formula=True,
    flag_value=FLAG_VALUE,
):
    """Deprecated. Use conversion.buoyancy"""

    warnings.warn("Deprecated. Use conversion.buoyancy", DeprecationWarning)

    return c.buoyancy(
        temperature_c,
        salinity_prac,
        pressure_dbar,
        latitude,
        longitude,
        window_size,
        use_modern_formula,
        flag_value,
    )


def _get_downcast_mask(
    dataset: xr.Dataset,
    control_variable: str,
    min_value: float = -np.inf,
    exclude_bad_scans=False,
) -> np.ndarray:
    """Creates a mask of the downcast of a dataset. The min index starts
    at the first value greater than min_value, and the max index is at
    the greatest value in the control_variable array. If flags are
    excluded those samples will be excluded but they won't be removed
    from the dataframe. In that case, the last valid sample on the
    downcast side of the peak will be determine the index to end on

    :param dataset: The dataset to create a downcast mask from
    :param control_variable: The control variable, typically pressure or
        depth
    :param min_value: Exclude values less than this from the result,
        defaults to -np.inf
    :param exclude_bad_scans: If True, these samples wil be skipped when
        determining the index of the max value, defaults to False
    :return: and array of booleans where true is a downcast
    """

    n_min = 0
    n_max = dataset[control_variable].idxmax(dim="sample").item()

    if exclude_bad_scans and "flag" in list(dataset.data_vars):
        # this diverges slightly from seasoft, where seasoft will take
        # the index of the max value that isn't a flag, this version
        # first finds the peak (regardless of flags) then finds the
        # first non-flagged value on the downcast side of the peak
        mask = (dataset["flag"].values != FLAG_VALUE) & (dataset["sample"].values < n_max)
        n_max = dataset[control_variable].where(mask).idxmax(dim="sample").item()

    if min_value > -np.inf:
        mask = (min_value < dataset[control_variable].values) & (dataset["sample"].values < n_max)
        n_min = dataset[control_variable].where(mask).idxmin(dim="sample").item()

    downcast_mask = (n_min <= dataset["sample"].values) & (dataset["sample"].values <= n_max)

    return downcast_mask


def _get_upcast_mask(
    dataset: xr.Dataset,
    control_variable: str,
    min_value: float = -np.inf,
    exclude_bad_scans=False,
) -> np.ndarray:
    """Creates a mask of the upcast of a dataset. The max index is at
    the greatest value of the control_variable values. If flags are
    excluded those samples will be excluded when determining the max
    index, but they won't be removed them from the dataframe. In that
    case, the first valid sample on the upcast side of the peak will be
    determine the index to start on. The min index is at the last value
    greater than min_value

    :param dataset: The dataset to create an upcast mask from
    :param control_variable: The control variable, typically pressure or
        depth
    :param min_value: Exclude values less than this from the result,
        defaults to -np.inf
    :param exclude_bad_scans: If True, these samples wil be skipped when
        determining the index of the max value, defaults to False
    :return: The upcast portion of the input dataframe
    """

    # n_min and n_max refer to the index of min/max value of the control
    # variable, meaning n_max will be less than n_min for an upcast
    n_min = len(dataset["sample"]) - 1
    n_max = dataset[control_variable].idxmax(dim="sample").item() + 1

    if exclude_bad_scans and "flag" in list(dataset.data_vars):
        # this diverges slightly from seasoft, where seasoft will take
        # the index of the max value that isn't a flag, this version
        # first finds the peak (regardless of flags) then finds the
        # first non-flagged value on the upcast side of the peak
        mask = (n_max < dataset["sample"].values) & (dataset["flag"].values != FLAG_VALUE)
        n_max = dataset[control_variable].where(mask).idxmax(dim="sample").item() + 1

    if min_value > -np.inf:
        # seasoft data processing doesn't include the max value
        mask = (n_max < dataset["sample"].values) & (min_value < dataset[control_variable].values)
        n_min = dataset[control_variable].where(mask).idxmin(dim="sample").item() - 1

    upcast_mask = (n_max <= dataset["sample"].values) & (dataset["sample"].values <= n_min)

    return upcast_mask


def split(
    dataset: xr.Dataset,
    control_variable: str,
    cast_type=CastType.BOTH,
    exclude_bad_scans=False,
    min_value=-np.inf,
    drop=False,
) -> xr.Dataset:
    """Adds a cast_type coordinate to the dataset that labels samples
    as "upcast", "downcast", or ""

    :param dataset: The dataset to add the cast_type to
    :param control_variable: The variable to determine where to split,
        typically pressure or depth
    :param split_mode: Determines which casts to include in the result,
        defaults to CastType.BOTH
    :param exclude_bad_scans: If True, rows with bad flags will be
        excluded when determining the max depth or pressure. Defaults
        to False
    :param min_value: Values below this will be excluded from the
        beginning of the downcast and end of the upcast
    :return: A new dataset with the cast_type coordinate
    """

    ds = dataset.copy()

    cast_type_coord = np.full(dataset.sizes["sample"], None, dtype=object)
    if cast_type in [CastType.BOTH, CastType.DOWNCAST]:
        downcast_mask = _get_downcast_mask(ds, control_variable, min_value, exclude_bad_scans)
        cast_type_coord[downcast_mask] = CastType.DOWNCAST.value

    if cast_type in [CastType.BOTH, CastType.UPCAST]:
        upcast_mask = _get_upcast_mask(ds, control_variable, min_value, exclude_bad_scans)
        cast_type_coord[upcast_mask] = CastType.UPCAST.value

    ds = ds.assign_coords(cast_type=("sample", cast_type_coord))

    if drop:
        if cast_type is CastType.BOTH:
            ds = ds.where(
                (ds["cast_type"] == CastType.UPCAST.value)
                | (ds["cast_type"] == CastType.DOWNCAST.value),
                drop=True,
            )
        if cast_type is CastType.DOWNCAST:
            ds = ds.where(ds["cast_type"] == CastType.DOWNCAST.value, drop=True)
        if cast_type is CastType.UPCAST:
            ds = ds.where(ds["cast_type"] == CastType.UPCAST.value, drop=True)

    return ds
