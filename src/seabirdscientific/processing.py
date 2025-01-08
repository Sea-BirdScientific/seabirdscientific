"""A collection of functions for processing converted SBS instrument
data.
"""

# Classes:

#     MinVelocityType (Enum)
#     WindowFilterType (Enum)

# Functions:

#     butterworth_filter (np.ndarray, float, float) -> np.ndarray
#     low_pass_filter (np.ndarray, float, float) -> np.ndarray
#     align_ctd (np.ndarray, float, float, float) -> np.ndarray
#     cell_thermal_mass (np.ndarray, np.ndarray, float, float, float) -> np.ndarray
#     loop_edit_pressure (
#         np.ndarray, float, np.ndarray, float, MinVelocityType, float,
#         float, float, bool, float, float, bool, bool,float
#     ) -> np.ndarray
#     loop_edit_depth (
#         np.ndarray, np.ndarray, float, MinVelocityType, float, float,
#         float, bool, float, float, bool, bool,float
#     ) -> np.ndarray
#     find_depth_peaks (np.ndarray, np.ndarray, bool, float, float, float) -> tuple[int, int]
#     min_velocity_mask (np.ndarray, float, float, int, int, bool) -> np.ndarray
#     mean_speed_percent_mask (np.ndarray, float, float, float, int, int, bool, int) -> np.ndarray
#     flag_by_minima_maxima (np.ndarray, np.ndarray, int, int, float)
#     bin_average (pd.DataFrame, str, float, bool, int, int, bool, bool)
#     wild_edit (np.ndarray, np.ndarray, float, float, int, float, bool,float) -> np.ndarray
#     flag_data (np.ndarray, np.ndarray, float, float, float, bool,float) -> np.ndarray
#     window_filter (
#         np.ndarray, np.ndarray, WindowFilterType, int, float, int, float, bool, float
#     ) -> np.ndarray
#     bouyancy_frequency (np.ndarray, np.ndarray, np.ndarray, float)
#     buoyancy (np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, bool, float)

# Native imports
from enum import Enum
import math

# Third-party imports
import gsw
import numpy as np
import pandas as pd
from scipy import signal, stats

# Sea-Bird imports

# Internal imports
from .conversion import depth_from_pressure
from . import eos80_processing as eos80


class MinVelocityType(Enum):
    """The minimum velocity type used with loop edit"""
    FIXED = 0
    PERCENT = 1


class WindowFilterType(Enum):
    """The window filter type. See CDT Processing in the docs for
    details.
    """
    BOXCAR = "boxcar"
    COSINE = "cosine"
    TRIANGLE = "triangle"
    GAUSSIAN = "gaussian"
    MEDIAN = "median"


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
    x: np.ndarray, offset: float, sample_interval: float, flag_value=-9.99e-29
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
    temperature_C: np.ndarray,
    conductivity_Sm: np.ndarray,
    amplitude: float,
    time_constant: float,
    sample_interval: float,
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

    a = 2 * amplitude / (sample_interval / time_constant + 2)
    b = 1 - (2 * a / amplitude)
    ctm = np.zeros(len(temperature_C))  # cell thermal mass
    corrected_conductivity = conductivity_Sm.copy()

    for n in range(1, len(ctm)):
        dc_dt = 0.1 * (1.0 + 0.006 * (temperature_C[n] - 20.0))
        dt = temperature_C[n] - temperature_C[n - 1]
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
    flag_value=-9.99e-29,
) -> np.ndarray:
    """Variation of loop_edit_depth that derives depth from pressure
    and latitude.

    :param pressure: A pressure array in dbar
    :param latitude: A single latitude value where the cast occurred
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

    :return: the input data with updated flags
    """

    depth = depth_from_pressure(pressure, latitude)
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
    min_velocity_type: MinVelocityType,
    min_velocity: float,
    window_size: float,
    mean_speed_percent: float,
    remove_surface_soak: bool,
    min_soak_depth: float,
    max_soak_depth: float,
    use_deck_pressure_offset: bool,
    exclude_flags: bool,
    flag_value=-9.99e-29,
) -> np.ndarray:
    """Marks scans determined to be part of a pressure loop as bad.

    Loop Edit marks scans bad by setting the flag value associated with
    the scan to badflag in data that has pressure slowdowns or reversals
    (typically caused by ship heave).

    :param depth: Salt water depth as an array of positive numbers.
        Colloquially used interchangeably with pressure in dbar
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

    if not exclude_flags:
        flag[:] = 0.0

    if use_deck_pressure_offset:
        min_soak_depth -= depth[0]
        max_soak_depth -= depth[0]

    (min_depth_n, max_depth_n) = find_depth_peaks(
        depth, flag, remove_surface_soak, flag_value, min_soak_depth, max_soak_depth
    )

    if min_velocity_type == MinVelocityType.FIXED:
        downcast_mask = min_velocity_mask(
            depth, sample_interval, min_velocity, min_depth_n, max_depth_n + 1, False
        )
        upcast_mask = min_velocity_mask(
            depth, sample_interval, min_velocity, max_depth_n, len(depth), True
        )

    elif min_velocity_type == MinVelocityType.PERCENT:
        diff_length = int(window_size / sample_interval)
        downcast_mask = mean_speed_percent_mask(
            depth,
            sample_interval,
            min_velocity,
            mean_speed_percent,
            min_depth_n,
            max_depth_n,
            False,
            diff_length,
        )
        upcast_mask = mean_speed_percent_mask(
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

    flag[~downcast_mask & ~upcast_mask] = flag_value

    flag_by_minima_maxima(depth, flag, min_depth_n, max_depth_n, flag_value)

    cast = np.array([0 for _ in range(len(flag))])
    cast[downcast_mask] = -1
    cast[upcast_mask] = 1
    return cast  # TODO: refactor to handle cast and flag in seperate functions


def find_depth_peaks(
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
            if flag[n] != flag_value and min_soak_depth < d < max_soak_depth
        )
    else:
        min_soak_depth_n = 0

    max_soak_depth_n = min(
        n for n, d in enumerate(depth) if flag[n] != flag_value and d > max_soak_depth
    )

    # beginning of possible downcast domain
    min_depth_n = min(
        (d, n)
        for n, d in enumerate(depth)
        if flag[n] != flag_value and min_soak_depth_n <= n < max_soak_depth_n
    )[1]

    # beginning of possible upcast domain
    if min_depth_n == len(depth) - 1:
        max_depth_n = -1
    else:
        max_depth_n = [n for n, d in enumerate(depth) if d == max(depth) and n > min_depth_n][0]

    return (min_depth_n, max_depth_n)


def min_velocity_mask(
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


def mean_speed_percent_mask(
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


def flag_by_minima_maxima(
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
    dataset: pd.DataFrame,
    bin_variable: str,
    bin_size: float,
    include_scan_count: bool,
    min_scans: int,
    max_scans: int,
    exclude_bad_scans: bool,
    interpolate: bool,
):
    """Averages data into a number of set bins, averaging values within
    each bin.

    Bins with fewer than min_scans or larger than max_scans get removed.
    Updates the dataset argument to be in bins.

    :param dataset: Data frame containing all data to group into bins.
        Will be modified by this function.
    :param bin_variable: name of the variable to use bins for. Expected
        to be one of "pressure", "depth", "scan_number", or "time"
    :param bin_size: size for each bin of data
    :param include_scan_count: indicates whether to include a column in
        the updated dataset for the number of scans in each bin
    :param min_scans: indicates the minimum number of scans required in
        a bin for it to be included in the dataset. Note that bins with
        0 scans are always dropped.
    :param max_scans: indicates the maximum number of scans allowed to
        be in a bin for it to be included in the dataset. Bins with more
        scans than max_Scans will be dropped.
    :param exclude_bad_scans: indicates whether to include bad scans
        within the bins.
    :param interpolate: indicates whether to interpolate the balues in
        the bins after averaging. Only possible for pressure or depth
        bins

    """
    # TODO: Currently unimplemented: surface bin, begin/end scans to omit,
    #     processing upcast vs downcast vs both, time bins, scan number bins

    if exclude_bad_scans and "flag" in dataset.columns:
        # remove all scans marked as bad (i.e. flag greater than 0)
        dataset.drop(dataset[dataset["flag"] > 0].index, inplace=True)
    bin_row = dataset[
        bin_variable
    ].to_numpy()  # pd series containing the variable we want to bin for, converted to ndarray
    max_val = np.amax(bin_row)

    bin_min = bin_size - (bin_size / 2.0)  # min value of first bin

    # create the bins to sort into
    # note that we create an inital "dummy" bin from -bin_min to bin_min to catch extra low values
    # TODO: Above issue can be addressed by surface bin likely
    bin_targets = np.arange(start=bin_min - bin_size, stop=max_val + bin_size, step=bin_size)

    # setup bins to indicate where each index should be sorted into
    bins = np.digitize(x=bin_row, bins=bin_targets, right=True)

    dataset["bin_number"] = bins

    dataset = dataset.groupby("bin_number").mean()

    # get the number of scans in each bin
    nbin_unfiltered = np.bincount(bins)
    nbin = np.delete(nbin_unfiltered, np.where(nbin_unfiltered == 0))
    dataset["nbin"] = nbin

    dataset.drop(dataset[dataset["nbin"] < min_scans].index, inplace=True)
    dataset.drop(dataset[dataset["nbin"] > max_scans].index, inplace=True)

    print(dataset.to_string())

    # TODO: validate that this is running correctly
    if interpolate:
        new_dataset = dataset.copy()
        prev_row = None
        first_row = None
        first_row_index = None
        second_row = None
        for index, row in dataset.iterrows():
            if prev_row is None:
                # we'll come back to the first row at the end
                first_row = row
                first_row_index = index
            else:
                prev_presure = prev_row[
                    bin_variable
                ]  # use bin_variable since this could be pressure or depth
                curr_pressure = row[
                    bin_variable
                ]  # use bin_variable since this could be pressure or depth
                center_pressure = index * bin_size

                for col_name, col_data in dataset.items():
                    if col_name == "nbin" or col_name == "flag" or col_name == "bin_number":
                        continue
                    prev_val = prev_row[col_name]
                    curr_val = row[col_name]

                    # formula from the seasoft data processing manual, page 89, version 7.26.8-3
                    interpolated_val = (
                        (curr_val - prev_val)
                        * (center_pressure - prev_presure)
                        / (curr_pressure - prev_presure)
                    ) + prev_val
                    new_dataset.loc[index, col_name] = interpolated_val

            prev_row = row
            if index == 2:
                # save this so we can reference it at the end for interpolating the first row
                second_row = row

        # back to the first row now
        prev_presure = second_row[bin_variable]  # reference second row's value
        curr_pressure = first_row[
            bin_variable
        ]  # use bin_variable since this could be pressure or depth
        center_pressure = first_row_index * bin_size

        for col_name, col_data in dataset.items():
            if col_name == "nbin" or col_name == "flag" or col_name == "bin_number":
                continue
            prev_val = second_row[col_name]  # reference second row's value
            curr_val = first_row[col_name]

            # formula from the seasoft manual, page 89
            interpolated_val = (
                (curr_val - prev_val)
                * (center_pressure - prev_presure)
                / (curr_pressure - prev_presure)
            ) + prev_val
            new_dataset.loc[first_row_index, col_name] = interpolated_val

        dataset = new_dataset
        print(dataset.to_string())
    if not include_scan_count:
        dataset.drop("nbin", axis=1, inplace=True)


def wild_edit(
    data: np.ndarray,
    flags: np.ndarray,
    std_pass_1: float,
    std_pass_2: float,
    scans_per_block: int,
    distance_to_mean: float,
    exclude_bad_flags: bool,
    flag_value=-9.99e-29,
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
        flagged_data[lower_index:upper_index] = flag_data(
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
        flagged_data[len(data) - len(data) % scans_per_block :] = flag_data(
            data[lower_index:upper_index],
            flags,
            std_pass_1,
            std_pass_2,
            distance_to_mean,
            exclude_bad_flags,
            flag_value,
        )[scans_per_block - len(data) % scans_per_block :]

    return flagged_data


def flag_data(
    data: np.ndarray,
    flags: np.ndarray,
    std_pass_1: float,
    std_pass_2: float,
    distance_to_mean: float,
    exclude_bad_flags: bool,
    flag_value=-9.99e-29,
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
    window_type: WindowFilterType,
    window_width: int,
    sample_interval: float,
    half_width=1,
    offset=0.0,
    exclude_flags=False,
    flag_value=-9.99e-29,
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

    # convert flags to nan for processing
    data = [d if d != flag_value else np.nan for d in data_in]

    if exclude_flags:
        data = [data[n] if f != flag_value else np.nan for n, f in enumerate(flags)]

    # for simplicity later, window indices start negative, end positive, and are centered on 0
    window_start = int(-(window_width - 1) / 2)
    window_end = int((window_width - 1) / 2 + 1)
    window = np.array([])

    # define the window filter
    if window_type == WindowFilterType.BOXCAR:
        window = signal.windows.boxcar(window_width)

    elif window_type == WindowFilterType.COSINE:
        for n in range(window_start, window_end):
            window = np.append(window, np.cos((n * np.pi) / (window_width + 1)))

    elif window_type == WindowFilterType.TRIANGLE:
        window = signal.windows.triang(window_width)

    elif window_type == WindowFilterType.GAUSSIAN:
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
        if data_start <= n and n < data_end:
            value = 0

            if window_type == WindowFilterType.MEDIAN:
                window = data[n + window_start : n + window_end]
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
    temp_conservative_subset: np.ndarray,
    salinity_abs_subset: np.ndarray,
    pressure_dbar_subset: np.ndarray,
    gravity: float,
):
    """Calculates an N^2 value (buoyancy frequency) for the given window
    of temperature, salinity, and pressure, at the given latitude.

    Expect temperature as conservative temperature, salinity as abslute
    salinity, and pressure as dbar, all of the same length. Performs the
    calculation using TEOS-10 and specific volume.

    :param temp_conservative_subset: temperature values for the given
        window
    :param salinity_abs_subset: salinity values for the given window
    :param pressure_dbar_subset: pressure values for the given window
    :param gravity: gravity value

    :return: A single N^2 [Brunt-Väisälä (buoyancy) frequency]
    """

    db_to_pa = 1e4
    # Wrap these as a length-1 array so that GSW accepts them
    pressure_bar = [np.mean(pressure_dbar_subset)]
    temperature_bar = [np.mean(temp_conservative_subset)]
    salinity_bar = [np.mean(salinity_abs_subset)]

    # Compute average specific volume, temp expansion ceoff,
    # and saline contraction coeff over window
    (v_bar, alpha_bar, beta_bar) = gsw.specvol_alpha_beta(
        salinity_bar, temperature_bar, pressure_bar
    )

    # Estimate vertical gradient of conservative temp
    dct_dp = stats.linregress(pressure_dbar_subset, temp_conservative_subset)
    # TODO: error handling with r, p, std_error

    # Estimate vertical gradient of absolute salinity
    dsa_dp = stats.linregress(pressure_dbar_subset, salinity_abs_subset)
    # TODO: error handling with r, p, std_error

    # Compute N2 combining computed ceofficients and vertical gradients.
    # we index into v_bar, alpha_bar, and beta_bar as they are all arrays of len 1
    n2 = gravity**2 / (v_bar[0] * db_to_pa)
    n2 *= beta_bar[0] * dsa_dp.slope - alpha_bar[0] * dct_dp.slope
    return n2


def buoyancy(
    temperature_c: np.ndarray,
    salinity_prac: np.ndarray,
    pressure_dbar: np.ndarray,
    latitude: np.ndarray,
    longitude: np.ndarray,
    window_size: float,
    use_modern_formula=True,
    flag_value=-9.99e-29,
):
    """Calculates the 4 buoyancy values based off the incoming data.

    Data is expected to have already been binned via Bin_Average using
    decibar pressure bins. All arrays are expected to be the same
    length, except for latitude and longitude, which can be length 1.
    Optionally can use the former calculation for N^2 from the SBE Data
    Processing Manual, but defaults to a newer formula using TEOS-10.

    :param temperature_c: Temperature in ITS-90 degrees C
    :param salinity_prac: Practical salinity in PSU
    :param pressure_dbar: Pressure in dbar
    :param latitude: latitude values. If length 1, gets applied to all
        values.
    :param longitude: longitude values. If length 1, gets applied to all
        values.
    :param window_size: window size to use. If this number is smaller
        than the binned window size, round up to a minium of 3 scans.
        I.E. uses the center scan and one scan on each side of it at the
        very least
    :param use_modern_formula: Whether to use a modern formula for
        calculating N^2. Defaults to true.
    :param flag_value: Bad Flag value to use for marking bad scans.
        Defaults to -9.99e-29

    :return: dataframe with 4 columns, one with each calculated
        variable. The columns are as follows:
    """

    salinity_prac, temperature_c, pressure_dbar, latitude, longitude = np.broadcast_arrays(
        salinity_prac, temperature_c, pressure_dbar, latitude, longitude
    )
    pressure_dbar = pressure_dbar.astype(np.double)

    # Get the original bin size that we're working with, using the
    # second and third bin so we don't have to worry about the surface
    # bin
    original_bin_size = abs(pressure_dbar[2] - pressure_dbar[1])

    # Calculates how many scans to have on either side of our median
    # point, but need at least 1 (for a total of 3 scans)
    scans_per_side = max(math.floor(window_size / original_bin_size / 2), 1)

    salinity_abs = gsw.SA_from_SP(salinity_prac, pressure_dbar, longitude, latitude)
    temperature_conservative = gsw.CT_from_t(salinity_abs, temperature_c, pressure_dbar)

    result = pd.DataFrame()

    # create our result np.ndarrays with the flag value as default
    result["N2"] = np.full(len(temperature_c), flag_value)
    result["N"] = np.full(len(temperature_c), flag_value)
    result["E"] = np.full(len(temperature_c), flag_value)
    result["E10^-8"] = np.full(len(temperature_c), flag_value)

    # start loop at scans_per_side
    for i in range(scans_per_side, len(temperature_conservative) - scans_per_side):
        min_index = i - scans_per_side
        max_index = (
            i + scans_per_side + 1
        )  # add + 1 because slicing does not include the max_index

        pressure_subset = pressure_dbar[min_index:max_index]
        temperature_cons_subset = temperature_conservative[min_index:max_index]
        temperature_its_subset = temperature_c[min_index:max_index]
        salinity_subset = salinity_abs[min_index:max_index]

        pressure_bar = [np.mean(pressure_subset)]
        gravity = gsw.grav([latitude[i]], pressure_bar)[0]

        if use_modern_formula:
            salinity_subset = salinity_abs[min_index:max_index]
            n2 = bouyancy_frequency(
                temperature_cons_subset, salinity_subset, pressure_subset, gravity
            )
        else:
            salinity_subset = salinity_prac[min_index:max_index]
            n2 = eos80.bouyancy_frequency(
                temperature_its_subset, salinity_subset, pressure_subset, gravity
            )

        result.at[i, "N2"] = n2
        if n2 >= 0:
            result.at[i, "N"] = math.sqrt(n2) * 3600 / (2 * np.pi)
        else:
            result.at[i, "N"] = np.nan
        result.at[i, "E"] = n2 / gravity
        result.at[i, "E10^-8"] = result.at[i, "E"] * 1e8

    return result
