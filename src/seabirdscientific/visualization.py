"""Functions to support data visualization."""

# Native imports
import json
import warnings
from dataclasses import dataclass
from datetime import timedelta
from logging import getLogger
from pathlib import Path
from typing import Dict, List, Literal, Optional, Union

# Third-party imports
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import xarray as xr
from plotly import subplots

# Internal imports
from seabirdscientific.interpret_sbs_variable import interpret_sbs_variable

# TODO: check python version
# TODO: translate warning and error messages

logger = getLogger(__name__)


@dataclass
class Transform:
    slope: float
    offset: float


@dataclass
class AxisSettings:
    name: str  # e.g. "tv290C", "c0S/m"
    title: str  # e.g. "Temperature", "Conductivity"
    units: str  # e.g. "ITS-90 deg C", "S/m"
    range: List[float]  # e.g. [0, 50], [20, 100] one pair per variable
    marker: str  # e.g. "circle", "square"
    color: str  #  e.g. "#00576d", "#0083a4", "#8ae7ff"
    color_scale: List[
        List[Union[float, str]]
    ]  # e.g. [[0, "#00576d"],[0.4, "#0083a4"],[1, "#8ae7ff"]]
    transform: Transform


@dataclass
class Axes:
    x: List[AxisSettings]
    y: List[AxisSettings]
    z: List[AxisSettings]


@dataclass
class PlotSettings:
    axes: Axes
    title: str
    layout_type: Literal["multiple-x", "multiple-y", "contour"]
    show_loop_edit_flags: bool
    connect_gaps: bool
    flag_value: float
    invert_y_axis: bool
    title_font_size: float
    label_font_size: float
    marker_size: float
    series_names: List[str]
    annotate: bool
    background_color: str
    height: int
    width: int
    disabled: bool
    plot_previous_casts: str
    plot_all_scans: bool
    previous_scans_to_plot: int
    previous_interval_to_plot: timedelta


class ChartConfig:
    """Dataclass to contain chart information and plotly settings"""

    def __init__(
        self,
        title: str,
        x_names: List[str],
        y_names: List[str],
        z_names: List[str],
        chart_type: Literal["overlay", "subplots"],
        bounds: Optional[Dict[Literal["x", "y", "z"], Dict[int, List[int]]]] = None,
        x_titles: Optional[List[str]] = None,
        y_titles: Optional[List[str]] = None,
        z_titles: Optional[List[str]] = None,
        plot_loop_edit_flags=False,
        lift_pen_over_bad_data=False,
        flag_value=-9.99e-29,
    ):
        """Initializes a chart config object to store the names, units,
        and data format.

        Config parameters must be in the same order as the data.

        :param title: Title of the chart
        :param x_names: X axis names
        :param y_names: Y axis names
        :param z_names: Z axis names
        :param chart_type: string to select type of chart "overlay" for
            multiple datasets sharing one axis, or "subplots" for
            multiple datasets on separate subplots.
        :param bounds: Chart axis bounds, for example:
            {'x': {0: [2, 10], 1: [1, 100]}
        :param x_titles: List of titles corresponding to x axis dataset.
            Defaults to units in line charts and empty string in ts plot
        :param y_titles: List of titles corresponding to y axis dataset.
            Defaults to units in line charts and empty string in ts plot
        :param z_titles: List of titles corresponding to z axis dataset.
            Defaults to units in line charts and empty string in ts plot
        :param plot_loop_edit_flags: If true, data are plotted
            regardless of loop edit flag values, otherwise data are not
            plotted where loop edit flags are equal to the flag_value
        :param lift_pen_over_bad_data: If true, flagged data are not
            drawn and instead leave gaps in plotted lines. Otherwise,
            lines are drawn between points surrounding flagged data (see
            SeaSoft manual, page 122)
        :param flag_value: A user configurable value that identifies
            flagged data. Defaults to -9.99e-29
        """

        self.title = title
        x_info = [interpret_sbs_variable(name) for name in x_names]
        y_info = [interpret_sbs_variable(name) for name in y_names]
        z_info = [interpret_sbs_variable(name) for name in z_names]
        self.x_names = [info["name"] for info in x_info]
        self.y_names = [info["name"] for info in y_info]
        self.z_names = [info["name"] for info in z_info]
        self.x_units = [info["units"] for info in x_info]
        self.y_units = [info["units"] for info in y_info]
        self.z_units = [info["units"] for info in z_info]
        self.chart_type = chart_type
        self.bounds = bounds if bounds is not None else {}
        self.x_titles = x_titles if x_titles is not None else x_names
        self.y_titles = y_titles if y_titles is not None else y_names
        self.z_titles = z_titles if z_titles is not None else z_names
        self.plot_loop_edit_flags = plot_loop_edit_flags
        self.lift_pen_over_bad_data = lift_pen_over_bad_data
        self.flag_value = flag_value

        axes: List[Literal["x", "y", "z"]] = ["x", "y", "z"]
        for axis in axes:
            if axis not in self.bounds.keys():
                self.bounds[axis] = {}


def _mask_dataset(dataset: xr.Dataset, config: ChartConfig) -> xr.Dataset:
    """Apply plotting masks to the dataset based on chart configuration."""

    masked = dataset.copy(deep=True)

    for name, data_array in masked.data_vars.items():
        if np.issubdtype(data_array.dtype, np.number):
            masked[name] = data_array.where(data_array != config.flag_value)

    if not config.plot_loop_edit_flags and "flag" in masked.data_vars:
        masked = masked.where(masked["flag"].notnull())

    return masked


def parse_instrument_data(source: Union[str, Path, pd.DataFrame, xr.Dataset]) -> xr.Dataset:
    """Top level function for converting instrument data to numpy array.

    Currently supports pandas dataframes, json strings, or a Path to the
    following file types: .csv, .asc (comma separated only), .json.

    :param source: A JSON string, file path (.csv, .asc, .json), or
        pandas DataFrame

    :return: xarray dataset containing field names and data

    """

    try:
        if isinstance(source, xr.Dataset):
            data = source.to_dataframe().reset_index(drop=True)

        elif isinstance(source, pd.DataFrame):
            data = source.copy()

        elif isinstance(source, Path):
            suffix = source.suffix.lower()
            if suffix in (".csv", ".asc"):
                data = pd.read_csv(source)
            elif suffix == ".json":
                with open(source, mode="r") as js_data:
                    data = pd.DataFrame.from_dict(json.load(js_data), orient="columns")

        elif isinstance(source, str):
            data = pd.DataFrame.from_dict(json.loads(source), orient="columns")

        else:
            raise TypeError

        if "data" not in locals():
            raise NameError

        if not isinstance(data, pd.DataFrame):
            raise TypeError

        columns = [interpret_sbs_variable(column)["name"] for column in data.columns]
        for old_column, new_column in zip(data.columns, columns):
            data.rename(columns={old_column: new_column}, inplace=True)
        scan_coords = np.arange(len(data))
        dataset = xr.Dataset(coords={"scan": scan_coords})
        for column in data.columns:
            dataset[column] = xr.DataArray(data[column].to_numpy(), dims=["scan"])

        return dataset

    except (NameError, TypeError) as e:
        logger.error(e)
        return None


def select_subset(axis_names: list[str], data: xr.Dataset) -> xr.Dataset:
    """Takes a list of axis names and returns a data set for each name
    in the list.

    If axis_names is empty the function will return a Dataset of
    integers representing the scan count of the data. This could be
    used in a single series chart for example.

    Otherwise, the function will return a Dataset for each name in the
    list. This would be for a single xy chart or an overlay/subplot
    chart.

    Example:

    data = read_data("./example.csv")

    subset = select_subset(["T090C", "C0Sm"], data)

    :param axis_names: List of axis names corresponding to the data
    :param data: xarray dataset containing one-dimensional instrument
        variables

    :return: A tuple with the axis name and data

    """

    if len(axis_names) == 0:
        sample_count = data.sizes.get("scan", 0)
        return xr.Dataset(
            {"Scan Count": xr.DataArray(np.arange(sample_count), dims=["scan"])},
            coords={"scan": np.arange(sample_count)},
        )

    missing_names = [name for name in axis_names if name not in data.data_vars]
    if missing_names:
        raise KeyError(missing_names[0])

    return data[axis_names]


def plot_xy_chart(data: xr.Dataset, config: ChartConfig) -> go.Figure:
    """Takes instrument data and a config and plots an XY chart with one
    or more data sets.

    :param data: xarray dataset, such as that returned by read_cnv_file
        or read_hex_file
    :param config: Config object with various plotly
        settings

    :return: A plotly.graph_objects.Figure
    """

    figure = go.Figure()

    masked_data = _mask_dataset(data, config)
    x_data = select_subset(config.x_names, masked_data)
    y_data = select_subset(config.y_names, masked_data)

    # single data set
    if len(x_data.data_vars) == 1 and len(y_data.data_vars) == 1:
        figure = create_single_plot(x_data, y_data, config)

    # too many data sets
    elif len(x_data.data_vars) > 1 and len(y_data.data_vars) > 1:
        message = "Only one axis can support multiple data sets"
        logger.warning(message)
        warnings.warn(message)

    # multiple data sets
    elif len(x_data.data_vars) > 1 or len(y_data.data_vars) > 1:
        if config.chart_type == "overlay":
            figure = create_overlay(x_data, y_data, config)
        elif config.chart_type == "subplots":
            figure = create_subplots(x_data, y_data, config)

    else:
        # getting here should not be possible unless data and config are
        # altered outside their init functions
        raise ValueError

    if not config.lift_pen_over_bad_data:
        figure.update_traces(connectgaps=True)

    return figure


def create_single_plot(x: xr.Dataset, y: xr.Dataset, config: ChartConfig) -> go.Figure:
    """Creates a single XY plot, with one or more data sets.

    If there are multiple datasets for the x or y axis, an overlay plot
    will be generated

    :param x: xarray dataset of data for the x axis
    :param y: xarray dataset of data for the y axis
    :param config: Dataclass with settings for the plotly chart

    :return: A plotly.graph_objects.Figure displaying the provided x y
        data
    """

    x_names = list(x.data_vars)
    y_names = list(y.data_vars)

    if any(name in x_names for name in y_names):
        logger.warning("Duplicate data names will be omitted")

    x_name = x_names[0]
    y_name = y_names[0]
    figure = go.Figure()
    figure.add_trace(
        go.Scatter(
            x=x[x_name].data,
            y=y[y_name].data,
            mode="lines",
            name=y_name,
        )
    )
    figure.update_layout(title=config.title)

    apply_single_config(figure, config)

    return figure


def create_subplots(x: xr.Dataset, y: xr.Dataset, config: ChartConfig) -> go.Figure:
    """Creates a chart with multiple subplots.

    :param x: xarray dataset of data for the x axis
    :param y: xarray dataset of data for the y axis
    :param config: Dataclass with settings for the plotly chart

    :return: A plotly.graph_objects.Figure with multiple subplots
    """

    x_names = list(x.data_vars)
    y_names = list(y.data_vars)

    if any(name in x_names for name in y_names):
        logger.warning("Duplicate data names will be omitted")

    figure = go.Figure()
    figure.layout.title = config.title

    if len(x_names) > 1 and len(y_names) == 1:
        figure = subplots.make_subplots(
            rows=1,
            cols=len(x_names),
            column_titles=x_names,
            y_title=y_names[0],
            figure=figure,
        )
        column = 1
        for x_column in x_names:
            figure.add_trace(
                go.Scatter(
                    x=x[x_column].data,
                    y=y[y_names[0]].data,
                    name=x_column,
                ),
                row=1,
                col=column,
            )
            column += 1
        apply_subplots_x_config(figure, config)

    elif len(x_names) == 1 and len(y_names) > 1:
        figure = subplots.make_subplots(
            rows=len(y_names),
            cols=1,
            row_titles=y_names,
            x_title=x_names[0],
            figure=figure,
        )
        row = 1
        for y_column in y_names:
            figure.add_trace(
                go.Scatter(
                    x=x[x_names[0]].data,
                    y=y[y_column].data,
                    name=y_column,
                ),
                row=row,
                col=1,
            )
            row += 1
        apply_subplots_y_config(figure, config)

    return figure


def create_overlay(x: xr.Dataset, y: xr.Dataset, config: ChartConfig) -> go.Figure:
    """Creates a chart with multiple datasets overlayed on one axis.

    :param x: xarray dataset of data for the x axis
    :param y: xarray dataset of data for the y axis
    :param config: Dataclass with settings for the plotly chart

    :return: A plotly.graph_objects.Figure
    """

    x_names = list(x.data_vars)
    y_names = list(y.data_vars)

    if any(name in x_names for name in y_names):
        logger.warning("Duplicate data names will be omitted")

    figure = go.Figure()
    figure.layout.title = config.title

    if len(x_names) > 1 and len(y_names) == 1:
        x_axis = 1
        for x_column in x_names:
            figure.add_trace(
                go.Scatter(
                    x=x[x_column].data,
                    y=y[y_names[0]].data,
                    name=x_column,
                    xaxis=f"x{x_axis}",
                )
            )
            x_axis += 1

        apply_overlay_config(figure, config)

    elif len(x_names) == 1 and len(y_names) > 1:
        y_axis = 1
        for y_column in y_names:
            figure.add_trace(
                go.Scatter(
                    x=x[x_names[0]].data,
                    y=y[y_column].data,
                    name=y_column,
                    yaxis=f"y{y_axis}",
                )
            )
            y_axis += 1

        apply_overlay_config(figure, config)

    return figure


def apply_single_config(figure: go.Figure, config: ChartConfig):
    """Updates various chart settings for single plots.

    :param figure: The figure being updated
    :param config: The user defined config being applied to the figure
    """

    figure.update_layout(
        xaxis={
            "title": "" if len(config.x_titles) < 1 else config.x_titles[0],
            "domain": [max(0, 0.1 * (len(config.y_titles) - 1)), 1],
            "range": None if len(config.bounds["x"]) < 1 else config.bounds["x"][0],
        },
        yaxis={
            "title": "" if len(config.y_titles) < 1 else config.y_titles[0],
            "domain": [max(0, 0.1 * (len(config.x_titles) - 1)), 1],
            "range": None if len(config.bounds["y"]) < 1 else config.bounds["y"][0],
        },
    )


def apply_subplots_x_config(figure: go.Figure, config: ChartConfig):
    """Updates various chart settings for charts with multiple x axes.

    Config parameters may contain upto 4 arguments per axis, and must be
    in the same order as the data. Hence all of the magic number
    indexing below

    :param figure: The figure being updated
    :param config: The user defined config being applied to the figure
    """

    y_range = None if len(config.bounds["y"]) < 1 else config.bounds["y"][0]

    figure.update_layout(
        xaxis={
            "title": "" if len(config.x_titles) < 1 else config.x_titles[0],
            "range": None if len(config.bounds["x"]) < 1 else config.bounds["x"][0],
        },
        xaxis2={
            "title": "" if len(config.x_titles) < 2 else config.x_titles[1],
            "range": None if len(config.bounds["x"]) < 2 else config.bounds["x"][1],
        },
        xaxis3={
            "title": "" if len(config.x_titles) < 3 else config.x_titles[2],
            "range": None if len(config.bounds["x"]) < 3 else config.bounds["x"][2],
        },
        xaxis4={
            "title": "" if len(config.x_titles) < 4 else config.x_titles[3],
            "range": None if len(config.bounds["x"]) < 4 else config.bounds["x"][3],
        },
        yaxis={
            "title": "" if len(config.y_titles) < 1 else config.y_titles[0],
            "range": y_range,
        },
        yaxis2={"range": y_range},
        yaxis3={"range": y_range},
        yaxis4={"range": y_range},
    )


def apply_subplots_y_config(figure: go.Figure, config: ChartConfig):
    """Updates various chart settings for charts with multiple y axes.

    Config parameters may contain upto 4 arguments per axis, and must be
    in the same order as the data. Hence all of the magic number
    indexing below

    :param figure: The figure being updated
    :param config: The user defined config being applied to the figure
    """

    x_range = None if len(config.bounds["x"]) < 1 else config.bounds["x"][0]

    figure.update_layout(
        xaxis={"range": x_range},
        xaxis2={"range": x_range},
        xaxis3={"range": x_range},
        xaxis4={
            "title": "" if len(config.x_titles) < 1 else config.x_titles[0],
            "range": x_range,
        },
        yaxis={
            "title": "" if len(config.y_titles) < 1 else config.y_titles[0],
            "range": None if len(config.bounds["y"]) < 1 else config.bounds["y"][0],
        },
        yaxis2={
            "title": "" if len(config.y_titles) < 2 else config.y_titles[1],
            "range": None if len(config.bounds["y"]) < 2 else config.bounds["y"][1],
        },
        yaxis3={
            "title": "" if len(config.y_titles) < 3 else config.y_titles[2],
            "range": None if len(config.bounds["y"]) < 3 else config.bounds["y"][2],
        },
        yaxis4={
            "title": "" if len(config.y_titles) < 4 else config.y_titles[3],
            "range": None if len(config.bounds["y"]) < 3 else config.bounds["y"][3],
        },
    )


def apply_overlay_config(figure: go.Figure, config: ChartConfig):
    """Updates various chart settings for charts with multiple y axes.

    Config parameters may contain upto 4 arguments per axis, and must be
    in the same order as the data.
    Hence all of the magic number indexing below

    :param figure: The figure being updated
    :param config: The user defined config being applied to the figure
    """

    figure.update_layout(
        xaxis={
            "title": "" if len(config.x_titles) < 1 else config.x_titles[0],
            "domain": [max(0, 0.1 * (len(config.y_titles) - 1)), 1],
            "range": None if len(config.bounds["x"]) < 1 else config.bounds["x"][0],
            "title_standoff": 5,
        },
        xaxis2={
            "title": "" if len(config.x_titles) < 2 else config.x_titles[1],
            "overlaying": "x",
            "position": 0.1,
            "range": None if len(config.bounds["x"]) < 2 else config.bounds["x"][1],
            "title_standoff": 5,
        },
        xaxis3={
            "title": "" if len(config.x_titles) < 3 else config.x_titles[2],
            "overlaying": "x",
            "position": 0.2,
            "range": None if len(config.bounds["x"]) < 3 else config.bounds["x"][2],
            "title_standoff": 5,
        },
        xaxis4={
            "title": "" if len(config.x_titles) < 4 else config.x_titles[3],
            "overlaying": "x",
            "position": 0.3,
            "range": None if len(config.bounds["x"]) < 4 else config.bounds["x"][3],
            "title_standoff": 5,
        },
        yaxis={
            "title": "" if len(config.y_titles) < 1 else config.y_titles[0],
            "domain": [0.1 * (len(config.x_titles)), 1],
            "position": 0,
            "range": None if len(config.bounds["y"]) < 1 else config.bounds["y"][0],
            "title_standoff": 5,
        },
        yaxis2={
            "title": "" if len(config.y_titles) < 2 else config.y_titles[1],
            "overlaying": "y",
            "position": 0.1,
            "range": None if len(config.bounds["y"]) < 2 else config.bounds["y"][1],
            "title_standoff": 5,
        },
        yaxis3={
            "title": "" if len(config.y_titles) < 3 else config.y_titles[2],
            "overlaying": "y",
            "position": 0.2,
            "range": None if len(config.bounds["y"]) < 3 else config.bounds["y"][2],
            "title_standoff": 5,
        },
        yaxis4={
            "title": "" if len(config.y_titles) < 4 else config.y_titles[3],
            "overlaying": "y",
            "position": 0.3,
            "range": None if len(config.bounds["y"]) < 3 else config.bounds["y"][3],
            "title_standoff": 5,
        },
    )


def plot_ts_chart(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    x_vec: np.ndarray,
    y_vec: np.ndarray,
    z_mat: np.ndarray,
    config: ChartConfig,
) -> go.Figure:
    """Overlays a scatter plot onto a contour plot to create a TS plot.
    Takes as args the xyz properties on a Contour object. In a future
    version these will be replaced with a single contour object

    :param x: absolute salinity
    :param y: conservative temperature
    :param z: potential density
    :param x_vec: absolute salinity vector
    :param y_vec: conservative temperature vector
    :param z_mat: potential density matrix
    :param config: Config object with key/values required by conversion
        function
    :return: A plotly.graph_objects.Figure
    """

    if len(config.x_names) > 1 or len(config.y_names) > 1 or len(config.z_names) > 1:
        logger.warning(
            "plot_ts_chart expects one data set for each axis. Extra data sets are ignored"
        )

    # Create 2 plots using plotly (not plotly express)
    # T-S diagram with x-y values colored by z
    trace1 = go.Scatter(
        x=x,
        y=y,
        mode="markers",
        marker={
            "color": z,
            "showscale": True,
            "size": 2,
            "colorbar": {
                # sigma_theta kg m-3
                "title": {
                    "text": config.z_titles[0],
                    "side": "top",
                }
            },
        },
    )

    # Contours in gray
    colorscale = [[0, "gray"], [1, "gray"]]
    trace2 = go.Contour(
        x=x_vec,
        y=y_vec,
        z=z_mat,
        showscale=False,
        colorscale=colorscale,
        contours={"coloring": "lines", "showlabels": True},
    )

    # Overlay the plots
    fig = subplots.make_subplots()
    fig.add_trace(trace1)
    fig.add_trace(trace2)

    # Centered title, axis labels
    fig.update_layout(
        height=600,
        width=600,
        title_text=config.title,
        title_x=0.5,
        xaxis_title=config.x_titles[0],
        yaxis_title=config.y_titles[0],
    )

    return fig
