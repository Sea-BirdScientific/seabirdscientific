#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""A collection of classes and functions to support visualization of data

Classes:

    ChartConfig
    ChartData

Functions:

    parse_instrument_data (Union[str, Path, pd.DataFrame]) -> pd.DataFrame
    select_subset (list[str], pd.DataFrame) -> pd.DataFrame
    plot_xy_chart (ChartData, ChartConfig) -> go.Figure:
    create_single_plot (pd.DataFrame, pd.DataFrame, ChartConfig) -> go.Figure
    create_subplots (pd.DataFrame, pd.DataFrame, ChartConfig) -> go.Figure:
    create_overlay (pd.DataFrame, pd.DataFrame, ChartConfig) -> go.Figure
    apply_single_config (go.Figure, ChartConfig)
    apply_subplots_x_config (go.Figure, ChartConfig)
    apply_subplots_y_config (go.Figure, ChartConfig)
    apply_overlay_config (go.Figure, ChartConfig)
    plot_ts_chart (np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, ChartConfig) -> go.Figure

"""

# Native imports
import json
from logging import getLogger
from typing import Union
from pathlib import Path

# Third-party imports
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly import subplots

# Internal imports
from .interpret_sbs_variable import interpret_sbs_variable

# TODO: check python version
# TODO: translate warning and error messages

logger = getLogger(__name__)

PLOTLY_PATH = "/node_modules/plotly.js-dist/plotly.js"


class ChartConfig:
    """Dataclass to contain chart information and plotly settings"""

    # TODO: replace chart_type with enum
    def __init__(
        self,
        title,
        x_names,
        y_names,
        z_names,
        chart_type,
        bounds={},
        x_titles=[""],
        y_titles=[""],
        z_titles=[""],
        plot_loop_edit_flags=False,
        lift_pen_over_bad_data=False,
        flag_value=-9.99e-29,
    ):
        """Initializes a chart config object to store the names, units, and data format.

        Config parameters must be in the same order as the data.

        Args:
            title (str): Title of the chart
            x_names (list[str]): X axis names
            y_names (list[str]): Y axis names
            z_names (list[str]): Z axis names
            chart_type (str): Optional string to select type of chart "overlay" for multiple datasets sharing one axis,
                or "subplots" for multiple datasets on separate subplots.
            bounds (dict): Chart axis bounds, structure as:
                bounds = {
                    'x': {0: [n, m], 1: [n, m], 2: [n, m], 3: [n, m]}
                    'y': {0: [n, m], 1: [n, m], 2: [n, m], 3: [n, m]}
                    'z': {0: [n, m], 1: [n, m], 2: [n, m], 3: [n, m]}
                    }
            x/y/z_titles (list[str]): List of titles corresponding to each dataset.
                Defaults to units in line charts and empty string in ts plots.
            plot_loop_edit_flags (bool): If true, data are plotted regardless of
                loop edit flag values, otherwise data are not plotted where loop
                edit flags are equal to the flag_value (-9.99e-29 by default)
            lift_pen_over_bad_data (bool): If true, flagged data are not drawn and
                instead leave gaps in plotted lines. Otherwise, lines are drawn
                between points surrounding flagged data (see SeaSoft manual, page 122)
            flag_value (float): A user configurable value that identifies flagged data
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
        self.bounds = bounds
        self.chart_type = chart_type
        self.x_titles = x_titles
        self.y_titles = y_titles
        self.z_titles = z_titles
        self.plot_loop_edit_flags = plot_loop_edit_flags
        self.lift_pen_over_bad_data = lift_pen_over_bad_data
        self.flag_value = flag_value

        for axis in ["x", "y", "z"]:
            if axis not in self.bounds.keys():
                self.bounds[axis] = {}


class ChartData:
    """Class to contain chart data and helper functions"""

    # TODO: move helper functions into ChartData

    def __init__(self, data_source, config):
        """Initializes an object to store chart data.

        Args:
            data_source (str | pd.DataFrame): A file path (.csv, .asc, .json),
            JSON string, or pandas DataFrame
            config (ChartConfig): A ChartConfig object to configure plotly
        """

        data = parse_instrument_data(data_source)
        if data is not None:
            data.mask(data == config.flag_value, inplace=True)
            if not config.plot_loop_edit_flags:
                mask = data["flag"].isnull()
                data.loc[mask, :] = np.nan
            self.x = select_subset(config.x_names, data)
            self.y = select_subset(config.y_names, data)
            self.z = select_subset(config.z_names, data)


def parse_instrument_data(source: Union[str, Path, pd.DataFrame]) -> pd.DataFrame:
    """Top level function for converting instrument data to numpy array.

    Currently supports pandas dataframes, json strings, and a string path to
    the following file types: .csv, .asc (comma separated), and .json.

    Args:
        source (str | | pathlib.Path | pd.DataFrame): A JSON string, file path (.csv, .asc, .json),
        or pandas DataFrame

    Returns:
        pd.DataFrame: Numpy ndarray containing field names and data

    Example:
        array = parse_instrument_data("./example.csv")
    """

    try:
        if isinstance(source, pd.DataFrame):
            data = source.copy()

        elif isinstance(source, Path):
            suffix = source.suffix.lower()
            if suffix == ".csv" or suffix == ".asc":
                data = pd.read_csv(source)
            elif suffix == ".json":
                with open(source, encoding="utf-8") as js_data:
                    data = pd.DataFrame.from_dict(json.load(js_data), orient="columns")

        elif isinstance(source, str):
            data = pd.DataFrame.from_dict(json.loads(source), orient="columns")

        else:
            raise Exception

        if "data" not in locals() or not isinstance(data, pd.DataFrame):
            raise Exception

        columns = [interpret_sbs_variable(column)["name"] for column in data.columns]
        for old_column, new_column in zip(data.columns, columns):
            data.rename(columns={old_column: new_column}, inplace=True)
        return data

    except Exception as e:
        logger.error(e)


def select_subset(axis_names: list[str], data: pd.DataFrame) -> pd.DataFrame:
    """Takes a list of axis names and returns a data set for each name in the list.

    If axis_names is empty the function will return a DataFrame
    of integers representing the sample count of the data.
    This could be used in a single series chart for example.

    Otherwise, the function will return a DataFrame for each name in the list.
    This would be for a single xy chart or an overlay/subplot chart.

    Args:
        axis_names (list of strings): List of axis names corresponding to the data
        data (DataFrame): The numpy DataFrame returned from read_data()

    Returns:
        tuple(string, DataFrame): A tuple with the axis name and data

    Example:
        data = read_data("./example.csv")

        subset = select_subset(["T090C", "C0Sm"], data)
    """

    if len(axis_names) == 0:
        return pd.DataFrame({"Sample Count": list(range(0, len(data)))})
    else:
        return data[axis_names]


def plot_xy_chart(data: ChartData, config: ChartConfig) -> go.Figure:
    """Takes instrument data and a config and plots an XY chart with one or more data sets.

    Args:
        data (ChartData): Data object with x, y, z, data selected according to the config
        config (ChartConfig): Config object with various plotly settings

    Returns:
        Figure: A plotly.graph_objects.Figure
    """

    figure = go.Figure()

    # single data set
    if len(data.x.columns) == 1 and len(data.y.columns) == 1:
        figure = create_single_plot(data.x, data.y, config)

    # too many data sets
    elif len(data.x.columns) > 1 and len(data.y.columns) > 1:
        logger.warning("Only one axis can support multiple data sets")
        # return go.Figure()

    # multiple data sets
    elif len(data.x.columns) > 1 or len(data.y.columns) > 1:
        if config.chart_type == "overlay":
            figure = create_overlay(data.x, data.y, config)
        elif config.chart_type == "subplots":
            figure = create_subplots(data.x, data.y, config)

    else:
        # getting here should not be possible unless data and config are altered outside their init functions
        raise Exception

    if not config.lift_pen_over_bad_data:
        figure.update_traces(connectgaps=True)

    return figure


def create_single_plot(x: pd.DataFrame, y: pd.DataFrame, config: ChartConfig) -> go.Figure:
    """Creates a single XY plot, with one or more data sets.

    If there are multiple datasets for the x or y axis, an overlay plot will be generated

    Args:
        x (ndarray): Numpy array of data for the x axis
        y (ndarray): Numpy array of data for the y axis
        config (ChartConfig): Dataclass with settings for the plotly chart

    Returns:
        go.Figure: A plotly.graph_objects.Figure displaying the provided x y data
    """

    if any(name in x.columns for name in y.columns):
        logger.warning("Duplicate data names will be omitted")

    figure = px.line(
        data_frame=pd.concat([x, y], axis=1),
        x=x.columns[0] if len(x.columns) == 1 else x.columns,
        y=y.columns[0] if len(y.columns) == 1 else y.columns,
        title=config.title,
        # markers=True,
    )

    apply_single_config(figure, config)

    return figure


def create_subplots(x: pd.DataFrame, y: pd.DataFrame, config: ChartConfig) -> go.Figure:
    """Creates a chart with multiple subplots.

    Args:
        x (DataFrame): Pandas DataFrame of data for the x axis
        y (DataFrame): Pandas DataFrame of data for the y axis
        config (ChartConfig): Dataclass with settings for the plotly chart

    Returns:
        go.Figure: A plotly.graph_objects.Figure with multiple subplots
    """

    if any(name in x.columns for name in y.columns):
        logger.warning("Duplicate data names will be omitted")

    figure = go.Figure()
    figure.layout.title = config.title

    if len(x.columns) > 1 and len(y.columns) == 1:
        figure = subplots.make_subplots(
            rows=len(y.columns),
            cols=len(x.columns),
            column_titles=list(x.columns),
            y_title=y.columns[0],
            figure=figure,
        )
        column = 1
        for x_column in x.columns:
            figure.add_trace(
                go.Scatter(
                    x=x[x_column],
                    y=y[y.columns[0]],
                    name=x_column,
                ),
                row=1,
                col=column,
            )
            column += 1
        apply_subplots_x_config(figure, config)

    elif len(x.columns) == 1 and len(y.columns) > 1:
        figure = subplots.make_subplots(
            rows=len(y.columns),
            cols=len(x.columns),
            row_titles=list(y.columns),
            x_title=x.columns[0],
            figure=figure,
        )
        row = 1
        for y_column in y.columns:
            figure.add_trace(
                go.Scatter(
                    x=x[x.columns[0]],
                    y=y[y_column],
                    name=y_column,
                ),
                row=row,
                col=1,
            )
            row += 1
        apply_subplots_y_config(figure, config)

    return figure


def create_overlay(x: pd.DataFrame, y: pd.DataFrame, config: ChartConfig) -> go.Figure:
    """Creates a chart with multiple datasets overlayed on one axis.

    Args:
        x (DataFrame): Pandas DataFrame of data for the x axis
        y (DataFrame): Pandas DataFrame of data for the y axis
        config (ChartConfig): Dataclass with settings for the plotly chart

    Returns:
        Figure: A plotly.graph_objects.Figure
    """

    if any(name in x.columns for name in y.columns):
        logger.warning("Duplicate data names will be omitted")

    figure = go.Figure()
    figure.layout.title = config.title

    if len(x.columns) > 1 and len(y.columns) == 1:
        x_axis = 1
        for x_column in x.columns:
            figure.add_trace(go.Scatter(x=x[x_column], y=y[y.columns[0]], name=x_column, xaxis=f"x{x_axis}"))
            x_axis += 1

        apply_overlay_config(figure, config)

    elif len(x.columns) == 1 and len(y.columns) > 1:
        y_axis = 1
        for y_column in y.columns:
            figure.add_trace(go.Scatter(x=x[x.columns[0]], y=y[y_column], name=y_column, yaxis=f"y{y_axis}"))
            y_axis += 1

        apply_overlay_config(figure, config)

    return figure


def apply_single_config(figure: go.Figure, config: ChartConfig):
    """Updates various chart settings for single plots.

    Args:
        figure (go.Figure): The figure being updated
        config (ChartConfig): The user defined config being applied to the figure
    """

    figure.update_layout(
        xaxis=dict(
            title="" if len(config.x_units) < 1 else config.x_units[0],
            domain=[max(0, 0.1 * (len(config.y_units) - 1)), 1],
            range=None if len(config.bounds["x"]) < 1 else config.bounds["x"][0],
        ),
        yaxis=dict(
            title="" if len(config.y_units) < 1 else config.y_units[0],
            domain=[max(0, 0.1 * (len(config.x_units) - 1)), 1],
            range=None if len(config.bounds["y"]) < 1 else config.bounds["y"][0],
        ),
    )


def apply_subplots_x_config(figure: go.Figure, config: ChartConfig):
    """Updates various chart settings for charts with multiple x axes.

    Config parameters may contain upto 4 arguments per axis, and must be in the same order as the data. Hence all of the magic number indexing below

    Args:
        figure (go.Figure): The figure being updated
        config (ChartConfig): The user defined config being applied to the figure
    """

    y_range = None if len(config.bounds["y"]) < 1 else config.bounds["y"][0]

    figure.update_layout(
        xaxis=dict(
            title="" if len(config.x_units) < 1 else config.x_units[0],
            range=None if len(config.bounds["x"]) < 1 else config.bounds["x"][0],
        ),
        xaxis2=dict(
            title="" if len(config.x_units) < 2 else config.x_units[1],
            range=None if len(config.bounds["x"]) < 2 else config.bounds["x"][1],
        ),
        xaxis3=dict(
            title="" if len(config.x_units) < 3 else config.x_units[2],
            range=None if len(config.bounds["x"]) < 3 else config.bounds["x"][2],
        ),
        xaxis4=dict(
            title="" if len(config.x_units) < 4 else config.x_units[3],
            range=None if len(config.bounds["x"]) < 4 else config.bounds["x"][3],
        ),
        yaxis=dict(title="" if len(config.y_units) < 1 else config.y_units[0], range=y_range),
        yaxis2=dict(range=y_range),
        yaxis3=dict(range=y_range),
        yaxis4=dict(range=y_range),
    )


def apply_subplots_y_config(figure: go.Figure, config: ChartConfig):
    """Updates various chart settings for charts with multiple y axes.

    Config parameters may contain upto 4 arguments per axis, and must be in the same order as the data. Hence all of the magic number indexing below

    Args:
        figure (go.Figure): The figure being updated
        config (ChartConfig): The user defined config being applied to the figure
    """

    x_range = None if len(config.bounds["x"]) < 1 else config.bounds["x"][0]

    figure.update_layout(
        xaxis=dict(range=x_range),
        xaxis2=dict(range=x_range),
        xaxis3=dict(range=x_range),
        xaxis4=dict(title="" if len(config.x_units) < 1 else config.x_units[0], range=x_range),
        yaxis=dict(
            title="" if len(config.y_units) < 1 else config.y_units[0],
            range=None if len(config.bounds["y"]) < 1 else config.bounds["y"][0],
        ),
        yaxis2=dict(
            title="" if len(config.y_units) < 2 else config.y_units[1],
            range=None if len(config.bounds["y"]) < 2 else config.bounds["y"][1],
        ),
        yaxis3=dict(
            title="" if len(config.y_units) < 3 else config.y_units[2],
            range=None if len(config.bounds["y"]) < 3 else config.bounds["y"][2],
        ),
        yaxis4=dict(
            title="" if len(config.y_units) < 4 else config.y_units[3],
            range=None if len(config.bounds["y"]) < 3 else config.bounds["y"][3],
        ),
    )


def apply_overlay_config(figure: go.Figure, config: ChartConfig):
    """Updates various chart settings for charts with multiple y axes.

    Config parameters may contain upto 4 arguments per axis, and must be in the same order as the data.
    Hence all of the magic number indexing below

    Args:
        figure (go.Figure): The figure being updated
        config (ChartConfig): The user defined config being applied to the figure
    """

    figure.update_layout(
        xaxis=dict(
            title="" if len(config.x_units) < 1 else config.x_units[0],
            domain=[max(0, 0.1 * (len(config.y_units) - 1)), 1],
            range=None if len(config.bounds["x"]) < 1 else config.bounds["x"][0],
        ),
        xaxis2=dict(
            title="" if len(config.x_units) < 2 else config.x_units[1],
            overlaying="x",
            position=0.1,
            range=None if len(config.bounds["x"]) < 2 else config.bounds["x"][1],
        ),
        xaxis3=dict(
            title="" if len(config.x_units) < 3 else config.x_units[2],
            overlaying="x",
            position=0.2,
            range=None if len(config.bounds["x"]) < 3 else config.bounds["x"][2],
        ),
        xaxis4=dict(
            title="" if len(config.x_units) < 4 else config.x_units[3],
            overlaying="x",
            position=0.3,
            range=None if len(config.bounds["x"]) < 4 else config.bounds["x"][3],
        ),
        yaxis=dict(
            title="" if len(config.y_units) < 1 else config.y_units[0],
            domain=[0.1 * (len(config.x_units)), 1],
            position=0,
            range=None if len(config.bounds["y"]) < 1 else config.bounds["y"][0],
        ),
        yaxis2=dict(
            title="" if len(config.y_units) < 2 else config.y_units[1],
            overlaying="y",
            position=0.1,
            range=None if len(config.bounds["y"]) < 2 else config.bounds["y"][1],
        ),
        yaxis3=dict(
            title="" if len(config.y_units) < 3 else config.y_units[2],
            overlaying="y",
            position=0.2,
            range=None if len(config.bounds["y"]) < 3 else config.bounds["y"][2],
        ),
        yaxis4=dict(
            title="" if len(config.y_units) < 4 else config.y_units[3],
            overlaying="y",
            position=0.3,
            range=None if len(config.bounds["y"]) < 3 else config.bounds["y"][3],
        ),
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
    """Calls various conversion functions (such as ct_sa_pd_from_t_s_p) and creates a TS plot.

    Args:
        data (ChartData): Data object with x, y, z, data selected according to the config
        config (ChartConfig): Config object with key/values required by conversion function
        min_salinity (optional, float): Minimum salinity to display
        lat, lon (optional, float) : Used to determine absolute salinity (SA)

    Returns:
        TSPlotValues: dataclass with various data for creating a TS plot
    """

    if len(config.x_names) > 1 or len(config.y_names) > 1 or len(config.z_names) > 1:
        logger.warning("plot_ts_chart expects one data set for each x, y, and z parameter. Extra data sets are ignored")

    # Create 2 plots using plotly (not plotly express)
    # T-S diagram with x-y values colored by z
    trace1 = go.Scatter(
        x=x,
        y=y,
        mode="markers",
        marker=dict(
            color=z,
            showscale=True,
            size=2,
            colorbar=dict(
                # sigma_theta kg m-3
                title=config.z_titles[0],
                titleside="top",
            ),
        ),
    )

    # Contours in gray
    colorscale = [[0, "gray"], [1, "gray"]]
    trace2 = go.Contour(
        x=x_vec,
        y=y_vec,
        z=z_mat,
        showscale=False,
        colorscale=colorscale,
        contours=dict(coloring="lines", showlabels=True),
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
