#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions to ...

Code documentation
------------------
"""


import plotly
import plotly.graph_objs as go


def scatter_tracer(x_data, y_data, tracer_name, tracer_mode,
                   ref_yaxis, marker_symbol, marker_size,
                   marker_color, line_dash):
    """ Creates a tracer object for a scatter plot.

        Parameters
        ----------
        x_data : list
            xaxis values corresponding to the maximum
            number of missing genes threshold values.
        y_data : list
            yaxis values corresponding to either the
            number of genomes included in the analysis
            at each threshold or to the number of genes
            present in 95%, 99%, 99.5%, or 100% of the
            genomes at each threshold.
        tracer_name : str
            name to show in the plot legend.
        tracer_mode : str
            type of symbols used used to represent data
            (lines, markers, lines+markers...).
        ref_yaxis : str
            the yaxis that will be used as reference.
        marker_symbol : str
            the type of symbol used for the markers.
        marker_size : int
            the size of the marker symbol.
        marker_color : str
            color of the markers and of the line (if any).
        line_dash : str
            type of line (solid, dash...).

        Returns
        -------
        tracer : plotly.graph_objs.Scatter
            a Plotly tracer with the data for a scatter
            plot (a group of data points).
    """

    tracer = go.Scatter(x=x_data,
                        y=y_data,
                        name=tracer_name,
                        mode=tracer_mode,
                        yaxis=ref_yaxis,
                        marker=dict(symbol=marker_symbol,
                                    size=marker_size,
                                    color=marker_color),
                        line=dict(dash=line_dash)
                        )

    return tracer


def create_figure(figure_data, figure_layout):
	"""
	"""

	fig = go.Figure(data=figure_data, layout=figure_layout)

	return fig


def plot_figure(figure, output_file, auto_open=False):
	"""
	"""

	plotly.offline.plot(figure, filename=output_file, auto_open=auto_open)
