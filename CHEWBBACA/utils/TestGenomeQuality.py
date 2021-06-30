#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module ...

Code documentation
------------------
"""


import os
import argparse
from copy import deepcopy
from itertools import chain

import plotly
import numpy as np
import pandas as pd
import plotly.graph_objs as go

try:
    from utils import (file_operations as fo,
                       matrix_manipulation as mm)
except:
    from CHEWBBACA.utils import (file_operations as fo,
                                 matrix_manipulation as mm)


def presence_absence_iterator(matrix, threshold):
    """ Determines the number of genes present in 0%, 95%, 99%,
        99.5% and 100% of genomes, excluding the genomes that
        have a number of missing genes greater than a defined
        threshold.

        Parameters
        ----------
        matrix : numpy.ndarray
            Array with the profiles determined by the AlleleCall
            process of chewBBACCA in a binary format
            (Presence/Absence matrix).
        threshold : int
            Maximum acceptable number of missing genes. Genomes
            with a number of missing genes above this threshold
            will be excluded from the calculations to determine
            the number of present genes.

        Returns
        -------
        exclude_genomes : list
            Genomes that had a number of missing loci above the
            threshold.
        vector : list
            Updated vector with the results for the current
            threshold.
    """

    iteration_matrix = deepcopy(matrix)

    # stores genes that were not in any of the genomes
    # stores genes that were in all genomes

    # determine genomes with number of missing loci above threshold
    missing = mm.missing_data_table(iteration_matrix)
    # select identifiers for genomes with missing loci above threshold
    exclude_genomes = missing.index[missing['missing'] > threshold].tolist()
    # remove genomes above threshold
    iteration_matrix = mm.remove_genomes(iteration_matrix, exclude_genomes)

    genomes = list(iteration_matrix.index)

    # this list includes 'FILE' header from matrix
    genes = list(iteration_matrix.columns)

    # determine cgMLST values
    genes_values = {}
    for g in genes:
        # get values counts
        presence_absence_counts = iteration_matrix[g].value_counts()
        total_presences = presence_absence_counts.get(1, 0)
        total_absences = presence_absence_counts.get(0, 0)

        # calculate percentage of genomes that have the locus
        gene_presence = float(len(genomes) - total_absences) / len(genomes)
        genes_values[g] = gene_presence

    # determine number of loci above each cgMLST threshold
    cgMLST_95 = [k for k, v in genes_values.items() if v >= 0.95]
    cgMLST_99 = [k for k, v in genes_values.items() if v >= 0.99]
    cgMLST_995 = [k for k, v in genes_values.items() if v >= 0.995]
    cgMLST_100 = [k for k, v in genes_values.items() if v >= 1.0]

    vector = [len(genomes),
              (cgMLST_95, len(cgMLST_95)),
              (cgMLST_99, len(cgMLST_99)),
              (cgMLST_995, len(cgMLST_995)),
              (cgMLST_100, len(cgMLST_100)),
              exclude_genomes]

    return [exclude_genomes, vector]


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


input_file = '/home/rfm/Desktop/rfm/bugfixing/chewBBACA_tutorial/chewBBACA_tutorial_copy/cgMLST_all.tsv'
missing_loci_threshold = 290
step = 50
output_directory = '/home/rfm/Desktop/rfm/bugfixing/chewBBACA_tutorial/chewBBACA_tutorial_copy/novel_TestGenomeQuality'
def main(input_file, missing_loci_threshold, step, output_directory):

    # create output directory, if it does not exist
    fo.create_directory(output_directory)

    # import matrix with allelic profiles
    matrix = pd.read_csv(input_file, header=0, index_col=0,
                         sep='\t', low_memory=False)

    # mask missing data
    print('Masking missing data...', end='')
    masked_matrix = matrix.apply(mm.replace_chars)
    print('done.')

    # build presence/absence matrix
    print('Building presence and absence matrix...', end='')
    presence_absence = mm.presAbs(masked_matrix, output_directory)
    print('done.')

    # generate list with all threshold values
    threshold_list = list(chain(range(0, missing_loci_threshold, step),
                                [missing_loci_threshold]))

    results = []
    for threshold in threshold_list:
        print('\r', 'Determining core-genome for threshold '
              '{0}/{1}'.format(threshold, threshold_list[-1]),
              end='')

        removed_genomes, vector = presence_absence_iterator(presence_absence,
                                                            threshold)

        results.append(vector)

    # plot legend labels
    labels = ['Genomes used',
              'Loci present in 95% genomes',
              'Loci present in 99% genomes',
              'Loci present in 99.5% genomes',
              'Loci present in 100% genomes']

    colors = ['#3690c0', '#ec7014', '#66bd63', '#807dba', '#fdbb84']
    # xaxis data is the list of thresholds
    x_data = threshold_list

    # number of genomes used per threshold
    y_genomes = [res[0] for res in results]
    # number of genes per threshold and per
    # genome presence percentage
    y_95 = [res[1][-1] for res in results]
    y_99 = [res[2][-1] for res in results]
    y_995 = [res[3][-1] for res in results]
    y_100 = [res[4][-1] for res in results]

    # group all yaxis datasets into list
    y_datasets = [y_genomes, y_95, y_99, y_995, y_100]

    # create all tracers
    tracers = []
    for d in range(len(y_datasets)):
        # tracer with used genomes data
        if d == 0:
            tracer = scatter_tracer(x_data,
                                    y_datasets[d],
                                    labels[d],
                                    'lines+markers',
                                    'y2',
                                    'diamond-dot',
                                    10,
                                    colors[d],
                                    'solid')
        # tracers for number of genes per threshold
        else:
            tracer = scatter_tracer(x_data,
                                    y_datasets[d],
                                    labels[d],
                                    'lines+markers',
                                    'y1',
                                    'circle',
                                    10,
                                    colors[d],
                                    'dash')

        tracers.append(tracer)

    # define layout attributes
    fig_layout = go.Layout(title='Test genomes quality',
                           xaxis=dict(title='Missing loci threshold',
                                      showline=True,
                                      mirror=True,
                                      linecolor='#EBEBEB',
                                      ticks='outside',
                                      tickcolor='#EBEBEB',
                                      showgrid=True,
                                      gridcolor='rgb(255,255,255)',
                                      range=[-5, missing_loci_threshold+5]),
                           yaxis=dict(title='Number of loci<br>'
                                            'in the core-genome',
                                      showline=True,
                                      linecolor='#EBEBEB',
                                      ticks='outside',
                                      tickcolor='#EBEBEB',
                                      showgrid=True,
                                      gridcolor='rgb(255,255,255)'),
                           yaxis2=dict(title='Number of genomes used<br>'
                                             'to determine the core-genome',
                                       showline=True,
                                       linecolor='#EBEBEB',
                                       ticks='outside',
                                       tickcolor='#EBEBEB',
                                       showgrid=True,
                                       gridcolor='rgb(255,255,255)',
                                       overlaying='y',
                                       side='right'),
                           plot_bgcolor='#EBEBEB'
                           )

    fig = go.Figure(data=tracers, layout=fig_layout)
    plot_file = os.path.join(output_directory, 'GenomeQualityPlot.html')
    plotly.offline.plot(fig, filename=plot_file, auto_open=False)

    print('\nResults available at: {0}'.format(output_directory))


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input-file', type=str,
                        required=True, dest='input_file',
                        help='Raw allele call matrix file (results_alleles.tsv).')

    parser.add_argument('-t', '--missing-loci-threshold', type=int,
                        required=True, dest='missing_loci_threshold',
                        help='Threshold value for the maximum number of '
                             'missing genes per genome.')

    parser.add_argument('-s', '--step', type=int,
                        required=True, dest='step',
                        help='Step between each threshold analysis. '
                             'The process will start with a threshold '
                             'value of 0 and will increment the threshold '
                             'value based on this step value.')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='Path to output directory.')

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = parse_arguments()
    main(**vars(args))
