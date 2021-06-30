#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions related to ...

Code documentation
------------------
"""


import os
import numpy as np
import pandas as pd


def replace_chars(column):
    """ Replaces all non-numeric characters
        in a column with allele identifiers.

        Parameters
        ----------
        column : pandas.core.series.Series
            Pandas dataframe column.

        Returns
        -------
        replace_missing : pandas.core.series.Series
            Input column with cells that only contain
            numeric characters.
    """

    # remove 'INF-' from inferred alleles
    replace_inf = column.replace(to_replace='INF-',
                                 value='', regex=True)
    # replace '*' in novel alleles from schemas in Chewie-NS
    # before replacing missing data cases to avoid replacing '*' with '0'
    replace_inf = replace_inf.replace(to_replace='\*',
                                      value='', regex=True)
    # replace missing data with '0'
    replace_missing = replace_inf.replace(to_replace='\D+.*',
                                          value='0', regex=True)

    return replace_missing


def binarize_matrix(column):
    """ Converts a Pandas dataframe column values
        into numeric values.

        Parameters
        ----------
        column : pandas.core.series.Series
            Pandas dataframe column.

        Returns
        -------
        Numpy array corresponding to the input column
        with numeric values equal to 1 for cells that
        had valid allele identifiers and equal to 0
        for cells that had missing data.
    """

    coln = pd.to_numeric(column)

    return np.int64(coln > 0)


def remove_genomes(matrix, genomesToRemove):
    """ Removes rows from a Pandas dataframe if the
        index identifier matches the identifier of
        a genome to remove.

        Parameters
        ----------
        matrix : pandas.core.frame.DataFrame
            Pandas dataframe with allelic profiles.
            Each row has the allelic profile of a genome
            and each column has the allele identifiers
            determined for a single gene.
        genomesToRemove : list
            List with the set of genomes to remove.

        Returns
        -------
        pruned_matrix : pandas.core.frame.DataFrame
            Input dataframe without the rows whose
            index matched an identifier of a genome
            to remove.
    """

    to_remove_bool = matrix.index.isin(genomesToRemove)
    pruned_matrix = matrix.loc[~ to_remove_bool]

    return pruned_matrix


def presAbs(matrix, output_directory):
    """ Creates a presence absence matrix.

        Parameters
        ----------
        matrix : pandas.core.frame.DataFrame
            Pandas dataframe with allelic profiles.
            Each row has the allelic profile of a genome
            and each column has the allele identifiers
            determined for a single gene.
        output_directory : str
            Path to the directory where the TSV file with
            the presence absence matrix will be stored.

        Returns
        -------
        presence_absence : pandas.core.frame.DataFrame
            Pandas dataframe with numeric values equal to
            1 for the cells that had valid allele identifiers
            and equal to 0 for missing data.
    """

    presence_absence = matrix.apply(binarize_matrix)

    pa_path = os.path.join(output_directory, 'Presence_Absence.tsv')
    presence_absence.to_csv(pa_path, sep='\t')

    return presence_absence


def missing_data_table(presence_absence):
    """ Determines missing data per genome.

        Parameters
        ----------
        presence_absence : pandas.core.frame.DataFrame
            Pandas dataframe with numeric values equal to
            1 for the cells that have valid allele identifiers
            and equal to 0 for missing data.

        Returns
        -------
        missing_data_df : pandas.core.frame.DataFrame
            Dataframe with number of missing genes and
            percentage of missing genes per genome.
    """

    _, n_genes = presence_absence.shape
    genes_present = presence_absence.apply(np.count_nonzero, axis=1)

    missing_data = {'FILE': presence_absence.index,
                    'missing': n_genes - genes_present,
                    'percentage': 1 - (genes_present / n_genes)}

    missing_data_df = pd.DataFrame(missing_data,
                                   columns=['FILE',
                                            'missing',
                                            'percentage'])

    return missing_data_df

