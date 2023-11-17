"""TODO: test_visualization docstring"""

# Native imports
import importlib.resources

# Third-party imports
import numpy as np
import pandas as pd
import pytest
import plotly.graph_objects as go

# Sea-Bird imports

# Internal imports
from sbs.visualize import visualization as dv

test_resources = importlib.resources.files('resources')


class TestParseChartData:
    @pytest.mark.parametrize("path", [
        test_resources / "example_pass.csv",
        test_resources / "example_pass.asc"
    ])
    def test_parse_instrument_data_csv_pass(self, path):
        assert isinstance(dv.parse_instrument_data(path), pd.DataFrame)

    def test_parse_instrument_data_json_pass(self):
        assert isinstance(dv.parse_instrument_data(test_resources / "example_pass.json"), pd.DataFrame)

    def test_parse_instrument_data_txt_error(self, caplog):
        dv.parse_instrument_data(test_resources / "example_fail_filetype.txt")
        assert('ERROR' in caplog.text)


class TestSelectSubset:
    def test_select_subset_empty_pass(self):
        data = dv.parse_instrument_data(test_resources / "example_pass.asc")
        assert np.all(dv.select_subset([], data) == pd.DataFrame({'Sample Count': [0, 1, 2]}))

    def test_select_subset_single_pass(self):
        data = dv.parse_instrument_data(test_resources / "example_pass.asc")
        subset = dv.select_subset(['Col1'], data)
        assert subset.columns == ['Col1']
        assert np.all(subset == pd.DataFrame({'Col1': [1, 4, 7]}))


class TestPlotXYChart:
    data_path = test_resources / 'example_pass.asc'
    config = dv.ChartConfig(
        title=data_path.name,
        x_names=["Col1"],
        y_names=[],
        z_names=[],
        chart_type=""
    )
    data = dv.ChartData(data_path, config)

    def test_plot_xy_chart_pass(self):
        assert isinstance(dv.plot_xy_chart(self.data, self.config), go.Figure)

    @pytest.mark.parametrize("chart_type", ["overlay", "subplots"])
    def test_plot_xy_chart_multiple_y_pass(self, chart_type):
        self.config.x_names = []
        self.config.y_names = ["Col1", "Col2", "Col3"]
        self.config.chart_type = chart_type
        data = dv.ChartData(self.data_path, self.config)
        assert isinstance(dv.plot_xy_chart(data, self.config), go.Figure)

    @pytest.mark.parametrize("chart_type", ["overlay", "subplots"])
    def test_plot_xy_chart_multiple_x_pass(self, chart_type):
        self.config.x_names = ["Col1", "Col2", "Col3"]
        self.config.y_names = []
        self.config.chart_type = chart_type
        assert isinstance(dv.plot_xy_chart(self.data, self.config), go.Figure)

    @pytest.mark.parametrize("chart_type", ["overlay", "subplots"])
    def test_plot_xy_chart_warn(self, chart_type):
        self.config.x_names = ["Col1", "Col2"]
        self.config.y_names = ["Col1", "Col2"]
        self.config.chart_type = chart_type
        data = dv.ChartData(self.data_path, self.config)
        with pytest.warns(UserWarning, match="Only one axis can support multiple data sets"):
            dv.plot_xy_chart(data, self.config)


class TestChartData:
    data_path = test_resources / 'example_pass.asc'
    config = dv.ChartConfig(
        title=data_path.name,
        x_names=["Col1"],
        y_names=[],
        z_names=[],
        chart_type="overlay"
    )
    data = dv.ChartData(data_path, config)

    def test_chart_data_x_error(self, caplog):
        with pytest.raises(KeyError):
            config = self.config
            config.x_names = ["Col999"]
            data = dv.ChartData(self.data_path, config)

    def test_chart_data_y_error(self, caplog):
        with pytest.raises(KeyError):
            config = self.config
            config.y_names = ["Col999"]
            data = dv.ChartData(self.data_path, config)

    def test_chart_data_z_error(self, caplog):
        with pytest.raises(KeyError):
            config = self.config
            config.z_names = ["Col999"]
            data = dv.ChartData(self.data_path, config)

    def test_chart_data_slash_pass(self):
        data_path = test_resources / 'example_fail_slash.csv'
        config = dv.ChartConfig(
            title=data_path,
            x_names = ["Col/1"],
            y_names = ["Col\\2"],
            z_names = ["Col-3"],
            chart_type=""
        )
        data = dv.ChartData(data_path, config)
        print(data)
        assert isinstance(data, dv.ChartData)


class TestPlotXYChart:
    def test_plot_xy_chart(self):
        config = dv.ChartConfig(
            title='title',
            x_names=["Col1"],
            y_names=["Col2"],
            z_names=[],
            chart_type='subplots'
        )
        source = '{"Col1": [1.0, 4.0, 7.0], "Col2": [2.0, 5.0, 8.0], "Col3": [3.0, 6.0, 9.0], "flag": [0, 0, 0]}'
        data = dv.ChartData(source, config)
        figure = dv.plot_xy_chart(data, config)
        assert isinstance(figure, go.Figure)
    # TODO: add failing tests


class TestPlotTSChart:
    def test_plot_ts_chart(self):
        config = dv.ChartConfig(
            title='title',
            x_names=["Col1"],
            y_names=["Col2"],
            z_names=["Col3"],
            chart_type=''
        )

        figure = dv.plot_ts_chart(
            np.array([1.91865655, 4.61772767, 6.98555581]),
            np.array([1.07317135, 4.21589043, 7.32884086]),
            np.array([1.45629369, 3.65123201, 5.38958412]),
            np.array([1.99946998, 2.0975469 , 2.19562383, 2.29370075, 2.39177767, 2.4898546 , 2.58793152, 2.68600844, 2.78408537, 2.88216229, 2.98023921, 3.07831613, 3.17639306, 3.27446998, 3.3725469 , 3.47062383, 3.56870075, 3.66677767, 3.7648546 , 3.86293152, 3.96100844, 4.05908537, 4.15716229, 4.25523921, 4.35331613, 4.45139306, 4.54946998, 4.6475469 , 4.74562383, 4.84370075, 4.94177767, 5.0398546 , 5.13793152, 5.23600844, 5.33408537, 5.43216229, 5.53023921, 5.62831613, 5.72639306, 5.82446998, 5.9225469 , 6.02062383, 6.11870075, 6.21677767, 6.3148546 , 6.41293152, 6.51100844, 6.60908537, 6.70716229, 6.80523921, 6.90331613, 7.00139306, 7.09946998]),
            np.array([1.96585421, 2.82299707, 3.68013993, 4.53728278, 5.39442564, 6.2515685 , 7.10871136, 7.96585421]),
            np.array([
                [1.55093735, 1.62975141, 1.70854623, 1.78732225, 1.86607991, 1.94481962, 2.02354179, 2.10224685, 2.18093517, 2.25960715, 2.33826317, 2.41690359, 2.49552878, 2.57413909, 2.65273487, 2.73131647, 2.8098842 , 2.88843839, 2.96697937, 3.04550744, 3.1240229 , 3.20252605, 3.28101719, 3.35949659, 3.43796454, 3.51642129, 3.59486713, 3.67330231, 3.75172708, 3.83014169, 3.90854638, 3.9869414 , 4.06532696, 4.14370331, 4.22207065, 4.30042922, 4.37877921, 4.45712084, 4.5354543 , 4.61377981, 4.69209754, 4.7704077 , 4.84871045, 4.927006  , 5.0052945 , 5.08357614, 5.16185108, 5.24011949, 5.31838153, 5.39663735, 5.47488712, 5.55313098, 5.63136907],
                [1.567872  , 1.64639465, 1.72489841, 1.80338371, 1.88185098, 1.96030064, 2.0387331 , 2.11714877, 2.19554802, 2.27393125, 2.35229883, 2.43065113, 2.50898851, 2.58731131, 2.66561989, 2.74391457, 2.82219568, 2.90046355, 2.97871849, 3.0569608 , 3.13519079, 3.21340874, 3.29161495, 3.3698097 , 3.44799326, 3.5261659 , 3.60432788, 3.68247946, 3.76062088, 3.8387524 , 3.91687425, 3.99498667, 4.07308989, 4.15118413, 4.22926961, 4.30734654, 4.38541513, 4.46347559, 4.54152812, 4.61957291, 4.69761016, 4.77564005, 4.85366276, 4.93167847, 5.00968735, 5.08768958, 5.16568533, 5.24367475, 5.321658  , 5.39963524, 5.47760662, 5.55557229, 5.6335324 ],
                [1.57398767, 1.65222672, 1.73044719, 1.80864951, 1.88683412, 1.96500142, 2.04315181, 2.12128571, 2.19940349, 2.27750554, 2.35559223, 2.43366392, 2.51172096, 2.5897637 , 2.66779249, 2.74580766, 2.82380952, 2.9017984 , 2.97977461, 3.05773846, 3.13569023, 3.21363022, 3.29155872, 3.36947599, 3.44738232, 3.52527797, 3.6031632 , 3.68103826, 3.75890339, 3.83675885, 3.91460487, 3.99244168, 4.07026951, 4.14808857, 4.22589909, 4.30370128, 4.38149534, 4.45928148, 4.53705989, 4.61483076, 4.69259429, 4.77035066, 4.84810005, 4.92584263, 5.00357859, 5.08130808, 5.15903126, 5.23674831, 5.31445938, 5.39216461, 5.46986417, 5.54755819, 5.62524683],
                [1.56954418, 1.64750705, 1.72545163, 1.80337835, 1.88128764, 1.95917989, 2.03705551, 2.1149149 , 2.19275845, 2.27058652, 2.34839949, 2.42619772, 2.50398156, 2.58175135, 2.65950743, 2.73725013, 2.81497977, 2.89269666, 2.97040112, 3.04809345, 3.12577393, 3.20344287, 3.28110052, 3.35874719, 3.43638312, 3.51400859, 3.59162385, 3.66922915, 3.74682474, 3.82441086, 3.90198774, 3.97955561, 4.0571147 , 4.13466522, 4.2122074 , 4.28974143, 4.36726752, 4.44478588, 4.5222967 , 4.59980016, 4.67729646, 4.75478578, 4.83226829, 4.90974418, 4.9872136 , 5.06467673, 5.14213373, 5.21958475, 5.29702996, 5.3744695 , 5.45190352, 5.52933217, 5.60675559],
                [1.55479071, 1.63248447, 1.7101602 , 1.78781833, 1.86545926, 1.94308342, 2.02069121, 2.098283  , 2.17585919, 2.25342014, 2.33096623, 2.40849781, 2.48601523, 2.56351882, 2.64100893, 2.71848588, 2.79595   , 2.87340158, 2.95084094, 3.02826837, 3.10568417, 3.18308863, 3.26048201, 3.3378646 , 3.41523665, 3.49259844, 3.56995021, 3.64729221, 3.72462469, 3.80194788, 3.87926202, 3.95656733, 4.03386404, 4.11115235, 4.1884325 , 4.26570467, 4.34296907, 4.42022591, 4.49747537, 4.57471765, 4.65195292, 4.72918137, 4.80640318, 4.88361851, 4.96082754, 5.03803042, 5.11522733, 5.19241841, 5.26960382, 5.34678371, 5.42395823, 5.50112751, 5.5782917 ],
                [1.5299663 , 1.60739764, 1.6848112 , 1.76220738, 1.83958661, 1.91694929, 1.99429581, 2.07162657, 2.14894194, 2.2262423 , 2.303528  , 2.3807994 , 2.45805686, 2.53530069, 2.61253125, 2.68974884, 2.76695379, 2.84414641, 2.921327  , 2.99849585, 3.07565326, 3.1527995 , 3.22993486, 3.30705961, 3.384174  , 3.46127829, 3.53837275, 3.6154576 , 3.69253311, 3.76959949, 3.84665699, 3.92370582, 4.00074621, 4.07777837, 4.15480251, 4.23181884, 4.30882756, 4.38582885, 4.46282293, 4.53980996, 4.61679013, 4.69376363, 4.77073063, 4.84769129, 4.92464578, 5.00159428, 5.07853692, 5.15547388, 5.2324053 , 5.30933133, 5.38625211, 5.46316779, 5.5400785 ],
                [1.49530023, 1.57247552, 1.64963323, 1.72677379, 1.8038976 , 1.88100507, 1.95809659, 2.03517254, 2.11223331, 2.18927926, 2.26631075, 2.34332813, 2.42033175, 2.49732194, 2.57429903, 2.65126334, 2.72821519, 2.80515489, 2.88208273, 2.958999  , 3.035904  , 3.11279801, 3.18968129, 3.26655412, 3.34341676, 3.42026947, 3.49711249, 3.57394607, 3.65077045, 3.72758586, 3.80439253, 3.88119069, 3.95798055, 4.03476232, 4.11153621, 4.18830244, 4.26506118, 4.34181265, 4.41855702, 4.49529449, 4.57202523, 4.64874942, 4.72546724, 4.80217886, 4.87888443, 4.95558412, 5.03227808, 5.10896648, 5.18564946, 5.26232716, 5.33899973, 5.41566732, 5.49233005],
                [1.45101254, 1.5279378 , 1.60484567, 1.68173658, 1.75861094, 1.83546915, 1.91231159, 1.98913866, 2.06595071, 2.14274813, 2.21953127, 2.29630047, 2.37305608, 2.44979843, 2.52652785, 2.60324465, 2.67994916, 2.75664167, 2.83332248, 2.90999188, 2.98665017, 3.06329761, 3.13993448, 3.21656105, 3.29317757, 3.3697843 , 3.44638148, 3.52296937, 3.5995482 , 3.67611819, 3.75267958, 3.82923258, 3.90577742, 3.98231431, 4.05884344, 4.13536503, 4.21187926, 4.28838634, 4.36488645, 4.44137977, 4.51786648, 4.59434677, 4.67082079, 4.74728872, 4.82375072, 4.90020695, 4.97665757, 5.05310273, 5.12954257, 5.20597724, 5.28240689, 5.35883165, 5.43525166]
            ]),
            config
        )
        assert isinstance(figure, go.Figure)
    # TODO: add failing tests


class TestChartBounds:
    data_path = test_resources / 'example_pass.asc'
    
    def test_chart_bounds_multiple_x_pass(self):
        config = dv.ChartConfig(
            title=self.data_path.name,
            x_names=["Col1", "Col2", "Col3", "Col4"],
            y_names=[],
            z_names=[],
            chart_type = "overlay",
            bounds = {
                'x': { 0: [1, 5], 1: [6, 8], 2: [5, 9], 3: [9, 11] }
            }
        )
        data = dv.ChartData(self.data_path, config)
        assert isinstance(dv.plot_xy_chart(data, config), go.Figure)
    
    def test_chart_bounds_multiple_y_pass(self):
        config = dv.ChartConfig(
            title=self.data_path.name,
            x_names=[],
            y_names=["Col1", "Col2", "Col3", "Col4"],
            z_names=[],
            chart_type = "overlay",
            bounds = {
                'y': { 0: [1, 5], 1: [6, 8], 2: [5, 9], 3: [9, 11] }
            }
        )
        data = dv.ChartData(self.data_path, config)
        assert isinstance(dv.plot_xy_chart(data, config), go.Figure)
