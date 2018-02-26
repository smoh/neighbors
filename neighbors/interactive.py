from bokeh.plotting import figure, show
from bokeh.models import TapTool, CustomJS, ColumnDataSource
from bokeh.models.widgets import DataTable, TableColumn, Button
from bokeh.layouts import row, column


class SelectPoints(object):
    """Bokeh layout for selecting data points interactively"""
    def __init__(self, x, y, radius=0.5, x_axis_label='x', y_axis_label='y'):
        """
        Make dashboard for selecting data points on x-y scatter plot

        x, y : array-like
            data to plot as x and y
        radius : float
            circle radius in data units
        x_axis_label, y_axis_label : str
            axis labels
        """
        p = figure(x_axis_label=x_axis_label, y_axis_label=y_axis_label)

        data = dict(idx=[], x=[], y=[])
        source = ColumnDataSource(data)
        columns = [
                TableColumn(field='idx', title='idx'),
                TableColumn(field="x", title=x_axis_label),
                TableColumn(field="y", title=y_axis_label),
            ]
        data_table = DataTable(source=source, columns=columns, width=400, height=400)
        # add reset button to clear DataTable
        buttoncallback = CustomJS(args={'source':source}, code="""
        source.data = {'idx':[], 'x':[], 'y':[]};
        """)
        buttonTableReset = Button(label="Reset table", button_type="success", callback=buttoncallback)

        cr = p.circle(x, y, radius=radius)
        # callback to add selected points to DataTable
        code = """
        var idx = circle.selected['1d'].indices;
        var table_data = {'idx':[], 'x': [], 'y': []};

        idx.forEach(function(i) {
            table_data['idx'].push(i);
            table_data['x'].push(circle.data['x'][i]);
            table_data['y'].push(circle.data['y'][i]);
        });

        source.data = table_data;
        """
        callback = CustomJS(args={'circle':cr.data_source, 'source':source}, code=code)
        p.add_tools(TapTool(callback=callback, renderers=[cr]))

        # show the results
        show(row([p, column([data_table, buttonTableReset])]))
        self.data_table = data_table
