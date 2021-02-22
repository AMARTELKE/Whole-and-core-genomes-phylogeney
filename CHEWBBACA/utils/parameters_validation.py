#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR

    Pedro Cerqueira
    github: @pedrorvc

    Rafael Mamede
    github: @rfm-targa

DESCRIPTION

"""


import os
import re
import sys
import argparse
import platform
import subprocess
import multiprocessing

try:
    from utils import (constants as ct,
                       file_operations as fo,
                       chewiens_requests as cr,
                       fasta_operations as fao)
except:
    from CHEWBBACA.utils import (constants as ct,
                                 file_operations as fo,
                                 chewiens_requests as cr,
                                 fasta_operations as fao)


class ModifiedHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):

    # prog is the name of the program 'ex: chewBBACA.py'
    def __init__(self, prog, indent_increment=2, max_help_position=56, width=100):
        super().__init__(prog, indent_increment, max_help_position, width)

    # override split lines method
    def _split_lines(self, text, width):
        lines = super()._split_lines(text, width) + ['']
        return lines

    def _format_action_invocation(self, action):
        if not action.option_strings:
            default = self._get_default_metavar_for_positional(action)
            metavar, = self._metavar_formatter(action, default)(1)
            return metavar

        else:
            parts = []

            # if the Optional doesn't take a value, format is:
            #    -s, --long
            if action.nargs == 0:
                parts.extend(action.option_strings)

            # if the Optional takes a value, format is:
            #    -s ARGS, --long ARGS
            else:
                default = self._get_default_metavar_for_optional(action)
                args_string = self._format_args(action, default)
                for option_string in action.option_strings:
                    parts.append(option_string)

                return '%s %s' % (', '.join(parts), args_string)

            return ', '.join(parts)

    def _get_default_metavar_for_optional(self, action):
        return action.dest.upper()


def arg_list(arg, arg_name):
    """
    """

    if isinstance(arg, list) is True:
        if len(arg) > 1:
            sys.exit('\nMultiple {0} values.'.format(arg_name))
        else:
            arg = arg[0]

    return arg


# custom functions to validate arguments type and value
def bsr_type(arg, min_value=ct.BSR_MIN, max_value=ct.BSR_MAX):
    """
    """

    arg = arg_list(arg, 'BLAST Score Ratio')

    try:
        schema_bsr = float(arg)
        if schema_bsr >= min_value and schema_bsr <= max_value:
            valid = schema_bsr
        elif schema_bsr < min_value or schema_bsr > max_value:
            sys.exit('\nBSR value is not contained in the '
                     '[0.0, 1.0] interval.')
    except Exception:
        sys.exit('\nInvalid BSR value of {0}. BSR value must be contained'
                 ' in the [0.0, 1.0] interval.'.format(arg))

    return valid


def minimum_sequence_length_type(arg, min_value=ct.MSL_MIN, max_value=ct.MSL_MAX):
    """
    """

    arg = arg_list(arg, 'minimum sequence length')

    try:
        schema_ml = int(arg)
        if schema_ml >= min_value and schema_ml <= max_value:
            valid = schema_ml
        elif schema_ml < min_value or schema_ml > max_value:
            sys.exit('\nInvalid minimum sequence length value. '
                     'Must be equal or greater than 0.')
    except Exception:
        sys.exit('\nInvalid minimum sequence length value used to '
                 'create schema. Value must be a positive integer.')

    return valid


def size_threshold_type(arg, min_value=ct.ST_MIN, max_value=ct.ST_MAX):
    """
    """

    arg = arg_list(arg, 'size threshold')

    try:
        schema_st = float(arg)
        if schema_st >= min_value and schema_st <= max_value:
            valid = schema_st
        elif schema_st < min_value or schema_st > max_value:
            sys.exit('\nInvalid size threshold value. '
                     'Must be contained in the [0.0, 1.0] interval.')
    except Exception:
        if arg in [None, 'None']:
            valid = None
        else:
            sys.exit('\nInvalid size threshold value used to '
                     'create schema. Value must be None or a '
                     'positive float in the [0.0, 1.0] interval.')

    return valid


def translation_table_type(arg, genetic_codes=ct.GENETIC_CODES):
    """
    """

    arg = arg_list(arg, 'translation table')

    try:
        schema_gen_code = int(arg)
        if schema_gen_code in genetic_codes:
            valid = schema_gen_code
        else:
            valid = False
    except Exception:
        valid = False

    if valid is False:
        # format available genetic codes into list
        lines = ['\t{0}: {1}'.format(k, v) for k, v in genetic_codes.items()]
        gc_table = '\n{0}\n'.format('\n'.join(lines))

        sys.exit('\nInvalid genetic code value.\nValue must correspond to '
                 'one of the accepted genetic codes\n\nAccepted genetic '
                 'codes:\n{0}'.format(gc_table))

    return valid


def validate_ws(arg, min_value=ct.WORD_SIZE_MIN, max_value=ct.WORD_SIZE_MAX):
    """
    """

    arg = arg_list(arg, 'word size')

    try:
        if arg == None:
            valid = 'None'
        else:
            word_size = int(arg)
            if word_size >= min_value and word_size <= max_value:
                valid = word_size
            else:
                sys.exit('\nWord size for the clustering step '
                         'must be equal or greater than {0} and '
                         'equal or smaller than {1}.'.format(min_value, max_value))
    except Exception:
        sys.exit('\nSchema created with invalid clustering word '
                 'size value.')

    return valid


def validate_cs(arg, min_value=ct.CLUSTERING_SIMILARITY_MIN, max_value=ct.CLUSTERING_SIMILARITY_MAX):
    """
    """

    arg = arg_list(arg, 'clustering similarity')

    try:
        if arg == None:
            valid = 'None'
        else:
            cluster_sim = float(arg)
            if cluster_sim >= min_value and cluster_sim <= max_value:
                valid = cluster_sim
            else:
                sys.exit('\nClustering similarity threshold value '
                         'must be contained in the [0.0, 1.0] '
                         'interval.')
    except Exception:
        sys.exit('\nSchema created with invalid clustering '
                 'threshold value.')

    return valid


def validate_rf(arg, min_value=ct.REPRESENTATIVE_FILTER_MIN, max_value=ct.REPRESENTATIVE_FILTER_MAX):
    """
    """

    arg = arg_list(arg, 'representative filter')

    try:
        if arg == None:
            valid = 'None'
        else:
            representative_filter = float(arg)
            if representative_filter >= min_value and representative_filter <= max_value:
                valid = representative_filter
            else:
                sys.exit('\nRepresentative filter threshold value '
                         'must be contained in the [0.0, 1.0] '
                         'interval.')
    except Exception:
        sys.exit('\nSchema created with invalid representative filter value.')

    return valid


def validate_if(arg, min_value=ct.INTRA_CLUSTER_MIN, max_value=ct.INTRA_CLUSTER_MAX):
    """
    """

    arg = arg_list(arg, 'intra-cluster filter')

    try:
        if arg == None:
            valid = 'None'
        else:
            intraCluster_filter = float(arg)
            if intraCluster_filter >= min_value and intraCluster_filter <= max_value:
                valid = intraCluster_filter
            else:
                sys.exit('\nIntra-cluster filter value '
                         'must be contained in the [0.0, 1.0] '
                         'interval.')
    except Exception:
        sys.exit('\nSchema created with invalid intra-cluster filter '
                 'value.')

    return valid


def validate_ptf(arg, input_path):
    """
    """

    arg = arg_list(arg, 'Prodigal training file')

    # perform this check in LoadSchema!!!
    #if arg == '':
    #    sys.exit('Cannot upload a schema that was created '
    #             'without a Prodigal training file.')

    schema_ptfs = [os.path.join(input_path, file)
                   for file in os.listdir(input_path) if '.trn' in file]
    if len(schema_ptfs) == 1:
        schema_ptf = schema_ptfs[0]
        ptf_hash = fo.hash_file(schema_ptf, 'rb')
        if ptf_hash == arg:
            valid = [schema_ptf, arg]
        else:
            sys.exit('Training file in schema directory is not the original.')
    elif len(schema_ptfs) > 1:
        sys.exit('More than one training file in schema directory.')
    elif len(schema_ptfs) == 0:
        sys.exit('Could not find a valid training file in schema directory.')

    return valid


def validate_ns_url(arg):
    """
    """

    if arg in ct.HOST_NS:
        ns_url = ct.HOST_NS[arg]
    else:
        ns_url = arg

    # sync schema has None by default to get ns_url in schema URI
    if ns_url is not None:
        # check if server is up
        conn = cr.check_connection(ns_url)
        if conn is False:
            sys.exit('Failed to establish a connection to the Chewie-NS '
                     'at {0}.'.format(ns_url))

    return ns_url


def validate_python_version():
    """
    """

    python_version = platform.python_version()

    try:
        assert tuple(map(int, python_version.split('.'))) >= (3, 4, 0)
    except AssertionError:
        print('Python version found: {} '.format(python_version))
        print('Please use version Python >= 3.4')
        sys.exit(0)

    return python_version


def verify_cpu_usage(cpu_to_use):
    """ Verify if the value provided for the number of CPU
        cores/threads does not exceed system limit or affect
        system performance.

        Parameters
        ----------
        cpu_to_use : int
            Value provided for the number of CPU cores/threads.

        Returns
        -------
        cpu_to_use : int
            Value of CPU cores/threads that will be used after
            determining if the provided value was safe.

    """
    total_cpu = multiprocessing.cpu_count()

    cpu_to_use = int(cpu_to_use)

    # do not allow a value greater than the number of cores
    if cpu_to_use >= total_cpu:
        print('Warning! You have provided a CPU core count value '
              'that is equal to or exceeds the number of CPU '
              'cores in your system!')
        # define a value that is safe according to the number of
        # available cores/threads
        if total_cpu > 2:
            cpu_to_use = total_cpu - 2
        elif total_cpu == 2:
            cpu_to_use = 1
        print('Resetting to: {0}'.format(cpu_to_use))
    elif cpu_to_use == (total_cpu - 1):
        print('Warning! You have provided a CPU core count value '
              'that is close to the maximum core count of your '
              'machine ({0}/{1}). This may affect your system '
              'responsiveness.'.format(cpu_to_use, total_cpu))

    return cpu_to_use


def is_exe(fpath):
    """ Determines if path points to file and
        if the file is executable.
    """

    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)


def which(program):
    """
    """

    proc = subprocess.Popen(['which', program],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)

    stdout, stderr = proc.communicate()

    program_path = stdout.decode('utf8').strip()

    return program_path


def check_blast(blast_path, major=ct.BLAST_MAJOR, minor=ct.BLAST_MINOR):
    """
    """

    # search for BLAST in PATH
    if not blast_path:
        blastp_path = which(ct.BLASTP_ALIAS)
        if not blastp_path:
            sys.exit('Could not find BLAST executables in PATH.')
        else:
            blast_path = os.path.dirname(blastp_path)
    # validate user-provided path
    else:
        blastp_path = os.path.join(blast_path, ct.BLASTP_ALIAS)
        executable = is_exe(blastp_path)
        if executable is False:
            sys.exit('Provided path does not contain BLAST executables.')

    # check BLAST version
    try:
        proc = subprocess.Popen([blastp_path, '-version'],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    except:
        sys.exit('Could not determine BLAST version.\n'
                 'Please verify that BLAST is installed and added to PATH.')

    stdout, stderr = proc.communicate()

    version_string = stdout.decode('utf8')
    version_pattern = r'^blastp:\s(?P<MAJOR>\d+).(?P<MINOR>\d+).(?P<REV>\d+).*'
    blast_version_pat = re.compile(version_pattern)

    match = blast_version_pat.search(version_string)
    if match is None:
        sys.exit('Could not determine BLAST version.')

    version = {k: int(v) for k, v in match.groupdict().items()}
    if version['MAJOR'] < major or (version['MAJOR'] >= major and version['MINOR'] < minor):
        sys.exit('Please update BLAST to version >= {0}.{1}.0'.format(major, minor))

    return blast_path


def check_prodigal(prodigal_path):
    """
    """

    # check Prodigal version
    try:
        proc = subprocess.Popen([prodigal_path, '-v'],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    except:
        sys.exit('Could not determine Prodigal version '
                 '({0} -v).\nPlease verify that Prodigal '
                 'is installed and added to PATH.'.format(prodigal_path))

    # prodigal version is captured in stderr
    stdout, stderr = proc.communicate()

    return True



def prompt_arguments(ptf_path, blast_score_ratio, translation_table,
                     minimum_length, size_threshold, version):
    """
    """

    prompt = ('It seems that your schema was created with chewBBACA '
              '2.1.0 or lower.\nIt is highly recommended that you run '
              'the PrepExternalSchema process to guarantee full '
              'compatibility with the new chewBBACA version.\nIf you '
              'wish to continue, the AlleleCall process will convert '
              'the schema to v{0}, but will not determine if schema '
              'structure respects configuration values.\nDo you wish '
              'to proceed?\n'.format(version))
    proceed = fo.input_timeout(prompt, ct.prompt_timeout)
    if proceed.lower() not in ['y', 'yes']:
        sys.exit('Exited.')

    print('\n')
    # determine parameters default values to include in config file
    if ptf_path is None:
        prompt = ('Full path to the Prodigal training file:\n')
        ptf_path = fo.input_timeout(prompt, ct.prompt_timeout)
        if ptf_path in ['', 'None']:
            ptf_path = False
        else:
            if os.path.isfile(ptf_path) is False:
                sys.exit('Provided path is not a valid file.')
    else:
        if os.path.isfile(ptf_path) is False:
            sys.exit('Provided path is not a valid file.')

    if blast_score_ratio is None:
        prompt = ('BLAST score ratio value:\n')
        blast_score_ratio = fo.input_timeout(prompt, ct.prompt_timeout)
    
    blast_score_ratio = bsr_type(blast_score_ratio)

    if translation_table is None:
        prompt = ('Translation table value:\n')
        translation_table = fo.input_timeout(prompt, ct.prompt_timeout)

    translation_table = translation_table_type(translation_table)

    if minimum_length is None:
        prompt = ('Minimum length value:\n')
        minimum_length = fo.input_timeout(prompt, ct.prompt_timeout)

    minimum_length = minimum_sequence_length_type(minimum_length)

    if size_threshold is None:
        prompt = ('Size threshold value:\n')
        size_threshold = fo.input_timeout(prompt, ct.prompt_timeout)

    size_threshold = size_threshold_type(size_threshold)

    return [ptf_path, blast_score_ratio, translation_table,
            minimum_length, size_threshold]


def auto_arguments(ptf_path, blast_score_ratio, translation_table,
                   minimum_length, size_threshold):
    """
    """

    if ptf_path is not None:
        if os.path.isfile(ptf_path) is False:
            sys.exit('Provided path to Prodigal training file '
                     'is not valid.')
    else:
        ptf_path = False

    blast_score_ratio = (blast_score_ratio
                         if blast_score_ratio is not None
                         else ct.DEFAULT_BSR)
    blast_score_ratio = bsr_type(blast_score_ratio)

    translation_table = (translation_table
                         if translation_table is not None
                         else ct.GENETIC_CODES_DEFAULT)
    translation_table = translation_table_type(translation_table)

    minimum_length = (minimum_length
                      if minimum_length is not None
                      else ct.MINIMUM_LENGTH_DEFAULT)
    minimum_length = minimum_sequence_length_type(minimum_length)

    size_threshold = (size_threshold
                      if size_threshold is not None
                      else ct.SIZE_THRESHOLD_DEFAULT)
    size_threshold = size_threshold_type(size_threshold)

    return [ptf_path, blast_score_ratio, translation_table,
            minimum_length, size_threshold]


def upgrade_legacy_schema(ptf_path, schema_directory, blast_score_ratio,
                          translation_table, minimum_length, version,
                          size_threshold, force_continue):
    """
    """

    if force_continue is False:
        values = prompt_arguments(ptf_path, blast_score_ratio,
                                  translation_table, minimum_length,
                                  size_threshold, version)
    # forced method, create with default args defined in constants
    else:
        values = auto_arguments(ptf_path, blast_score_ratio,
                                translation_table, minimum_length,
                                size_threshold)

    ptf_path, blast_score_ratio,\
    translation_table, minimum_length, size_threshold = values

    # copy training file to schema directory
    if ptf_path is not False:
        print('\nAdding Prodigal training file to schema...')
        shutil.copy(ptf_path, schema_directory)
        print('Created {0}'.format(os.path.join(schema_directory,
                                                os.path.basename(ptf_path))))

    # determine PTF hash
    ptf_hash = (fo.hash_file(ptf_path, 'rb')
                if ptf_path is not False
                else False)

    print('\nCreating file with schema configs...')
    # write schema config file
    schema_config = write_schema_config(blast_score_ratio, ptf_hash,
                                        translation_table, minimum_length,
                                        version, size_threshold,
                                        None, None, None, None,
                                        schema_directory)
    print('Created {0}'.format(os.path.join(schema_directory,
                                            '.schema_config')))

    print('\nCreating file with list of genes...')
    # create hidden file with genes/loci list
    genes_list_file = write_gene_list(schema_directory)
    print('Created {0}\n'.format(os.path.join(schema_directory,
                                              '.genes_list')))

    return [ptf_path, blast_score_ratio, translation_table,
            minimum_length, size_threshold]


def validate_ptf(ptf_path, schema_directory, schema_params,
                 unmatch_params):
    """
    """

    if ptf_path is None:
        # deal with multiple training files
        schema_ptfs = [file
                       for file in os.listdir(schema_directory)
                       if file.endswith('.trn')]
        if len(schema_ptfs) > 1:
            sys.exit('Found more than one Prodigal training '
                     'file in schema directory.\nPlease maintain '
                     'only the training file used in the schema '
                     'creation process.')
        elif len(schema_ptfs) == 1:
            ptf_path = os.path.join(schema_directory, schema_ptfs[0])
        elif len(schema_ptfs) == 0:
            print('There is no Prodigal training file in schema\'s '
                  'directory.')
            ptf_path = False
    # if user provides a training file
    elif ptf_path == 'False':
        ptf_path = False
    else:
        if os.path.isfile(ptf_path) is False:
            sys.exit('Invalid path for Prodigal training file.')

    # determine PTF checksum
    if ptf_path is not False:
        ptf_hash = fo.hash_file(ptf_path, 'rb')
    else:
        ptf_hash = False

    if ptf_hash not in schema_params['prodigal_training_file']:
        ptf_num = len(schema_params['prodigal_training_file'])
        if force_continue is False:
            if ptf_num == 1:
                print('Prodigal training file is not the one '
                      'used to create the schema.')
                prompt = ('Using this training file might lead to '
                          'results not consistent with previous runs '
                          'and invalidate the schema for usage with '
                          'the NS.\nContinue process?\n')
                ptf_answer = fo.input_timeout(prompt, ct.prompt_timeout)
            if ptf_num > 1:
                print('Prodigal training file is not any of the {0} '
                      'used in previous runs.'.format(ptf_num))
                prompt = ('Continue?\n')
                ptf_answer = fo.input_timeout(prompt, ct.prompt_timeout)
        else:
            ptf_answer = 'yes'

        if ptf_answer.lower() not in ['y', 'yes']:
            sys.exit('Exited.')
        else:
            schema_params['prodigal_training_file'].append(ptf_hash)
            unmatch_params['prodigal_training_file'] = ptf_hash

    return [ptf_path, schema_params, unmatch_params]


def solve_conflicting_arguments(schema_params, ptf_path, blast_score_ratio,
                                translation_table, minimum_length,
                                size_threshold, force_continue, config_file,
                                schema_directory):

    # run parameters values
    run_params = {'bsr': blast_score_ratio,
                  'translation_table': translation_table,
                  'minimum_locus_length':minimum_length,
                  'size_threshold': size_threshold}

    # determine user provided arguments values that differ from default
    unmatch_params = {k: v
                      for k, v in run_params.items()
                      if v not in schema_params[k] and v is not None}
    # determine arguments values not provided by user
    default_params = {k: schema_params[k][0]
                      for k, v in run_params.items()
                      if v is None}

    # update arguments for current run
    for k in run_params:
        if k in default_params:
            run_params[k] = default_params[k]

    if len(unmatch_params) > 0:
        print('Provided arguments values differ from arguments '
              'values used for schema creation:\n')
        params_diffs = [[p, ':'.join(map(str, schema_params[p])),
                         str(unmatch_params[p])]
                        for p in unmatch_params]
        params_diffs_text = ['{:^20} {:^20} {:^10}'.format('Argument', 'Schema', 'Provided')]
        params_diffs_text += ['{:^20} {:^20} {:^10}'.format(p[0], p[1], p[2]) for p in params_diffs]
        print('\n'.join(params_diffs_text))
        if force_continue is False:
            prompt = ('\nContinuing might lead to results not '
                      'consistent with previous runs.\nProviding '
                      'parameters values that differ from the ones '
                      'used for schema creation will also invalidate '
                      'the schema for uploading and synchronization '
                      'with Chewie-NS.\nContinue?\n')
            params_answer = fo.input_timeout(prompt, ct.prompt_timeout)
        else:
            params_answer = 'yes'

        if params_answer.lower() not in ['y', 'yes']:
            sys.exit('Exited.')
        else:
            # append new arguments values to configs values
            for p in unmatch_params:
                schema_params[p].append(unmatch_params[p])

    # default is to get the training file in schema directory
    ptf_path, schema_params, unmatch_params = validate_ptf(ptf_path, schema_directory,
                                                           schema_params, unmatch_params)
    run_params['ptf_path'] = ptf_path

    # save updated schema config file
    if len(unmatch_params) > 0:
        fo.pickle_dumper(schema_params, config_file)

    return [schema_params, unmatch_params, run_params]


def write_gene_list(schema_dir):
    """ Creates list with gene files in a schema and
        uses the pickle module to save the list to a file.

        Parameters
        ----------
        schema_dir : str
            Path to the directory with schema files.

        Returns
        -------
        A list with two elements. A boolean value that
        is True if the file with the list of genes was
        created and False otherwise. The second element
        is the path to the created file.
    """

    schema_files = [file for file in os.listdir(schema_dir) if '.fasta' in file]
    schema_list_file = fo.join_paths(schema_dir, ['.genes_list'])
    fo.pickle_dumper(schema_files, schema_list_file)

    return [os.path.isfile(schema_list_file), schema_list_file]


def write_schema_config(blast_score_ratio, ptf_hash, translation_table,
                        minimum_sequence_length, chewie_version, size_threshold,
                        word_size, window_size, clustering_sim, representative_filter,
                        intra_filter, output_directory):
    """ Writes chewBBACA's parameters values used to create
        a schema to a file.

        Parameters
        ----------
        blast_score_ratio : float
            BLAST Score Ratio value used to create the
            schema.
        ptf_hash : str
            BLAKE2 hash of the Prodigal training file
            content.
        translation_table : int
            Genetic code used to predict and translate
            coding sequences.
        minimum_sequence_length : int
            Minimum sequence length, sequences with a
            length value lower than this value are not
            included in the schema.
        chewie_version : str
            Version of the chewBBACA suite used to create
            the schema.
        size_threshold : float
            Sequence size variation percentage threshold,
            new alleles cannot have a length value that
            deviates +/- than this value in relation to the
            locus's representative sequence.
        word_size : int
            Word/k value used to cluster protein sequences
            during schema creation and allele calling.
        clustering_sim : float
            Proportion of k-mers/minimizers that two proteins
            need to have in common to be clustered together.
        representative_filter : float
            Proportion of k-mers/minimizers that a clustered
            protein has to have in common with the representative
            protein of the cluster to be considered the same gene.
        intra_filter : float
            Proportion of k-mers/minimizers that clustered
            proteins have to have in common to be considered
            of the same gene.
        output_directory : str
            Path to the output directory where the file with
            schema parameters values will be created.

        Returns
        -------
        A list with two elements. A boolean value that
        is True if the file with the parameters values was
        created and False otherwise. The second element
        is the path to the created file.
    """

    size_threshold = None if size_threshold in [None, 'None'] else float(size_threshold)

    params = {}
    params['bsr'] = [float(blast_score_ratio)]
    params['prodigal_training_file'] = [ptf_hash]
    params['translation_table'] = [int(translation_table)]
    params['minimum_locus_length'] = [int(minimum_sequence_length)]
    params['chewBBACA_version'] = [chewie_version]
    params['size_threshold'] = [size_threshold]
    params['word_size'] = [word_size]
    params['window_size'] = [window_size]
    params['cluster_sim'] = [clustering_sim]
    params['representative_filter'] = [representative_filter]
    params['intraCluster_filter'] = [intra_filter]

    config_file = os.path.join(output_directory, '.schema_config')
    fo.pickle_dumper(params, config_file)

    return [os.path.isfile(config_file), config_file]


def read_configs(schema_path, filename):
    """ Reads file with schema config values.

        Parameters
        ----------
        schema_path : str
            Path to the schema's directory.
        filename : str
            Name of the file that contains the config values.

        Returns
        -------
        configs : dict
            Dictionary with config names as keys and config
            values as values.
    """

    config_file = os.path.join(schema_path, filename)
    if os.path.isfile(config_file):
        # Load configs dictionary
        configs = fo.pickle_loader(config_file)
    else:
        sys.exit('Could not find a valid config file.')

    return configs


def check_input_type(input_path, output_file, parent_dir=None):
    """ Checks if the input path is for a file or for a
        directory. If the path is for a directory, the
        function creates a file with the list of paths
        to FASTA files in the directory.

        Parameters
        ----------
        input_path : str
            Path to file or directory.
        output_file : str
            Path to the output file with the list of FASTA
            files.

        Returns
        -------
        list_files : str
            Path to a file with a list of paths for FASTA
            files.

        Raises
        ------
        SystemExit
            - If there were no FASTA files in the directory.
            - If input path is not a valid path for a file or
              for a directory.
    """

    # check if input argument is a file or a directory
    if os.path.isfile(input_path):
        if parent_dir is not None:

            with open(input_path, 'r') as infile:
                lines = list(csv.reader(infile))

            lines = [os.path.join(parent_dir, f)
                     if parent_dir not in f
                     else f
                     for f in lines]

            with open(output_file, 'w') as outfile:
                outfile.write('\n'.join(lines))

            list_files = output_file

    elif os.path.isdir(input_path):
        # we need to get only files with FASTA extension
        files = os.listdir(input_path)
        files = fo.filter_files(files, ct.FASTA_SUFFIXES)
        # get absolute paths
        files = [os.path.join(input_path, file) for file in files]
        # filter any directories that might end with FASTA extension
        files = [file for file in files if os.path.isdir(file) is False]

        # only keep files whose content is typical of a FASTA file
        fasta_files = fao.filter_non_fasta(files)

        # if there are FASTA files
        if len(fasta_files) > 0:
            # store full paths to FASTA files
            with open(output_file, 'w') as f:
                for file in fasta_files:
                    f.write(file + '\n')
        else:
            sys.exit('\nCould not get input files. Please '
                     'provide a directory with FASTA files '
                     'or a file with the list of full paths '
                     'to the FASTA files and ensure that '
                     'filenames end with one of the '
                     'following suffixes: {0}.'
                     ''.format(ct.FASTA_SUFFIXES))

        list_files = output_file
    else:
        sys.exit('\nInput argument is not a valid directory or '
                 'file with a list of paths. Please provide a '
                 'valid input, either a folder with FASTA files '
                 'or a file with the list of full paths to FASTA '
                 'files (one per line).')

    return list_files


def process_header(process):
    """
    """

    header = 'chewBBACA - {0}'.format(process)
    hf = '='*(len(header)+4)
    print('{0}\n  {1}\n{0}'.format(hf, header, hf))


def check_ptf(ptf_path):
    """ Determines if path to Prodigal training file exists.

        Parameters
        ----------
        ptf_path : str
            Path to the Prodigal training file.

        Returns
        -------
        A list with a bool value that is True if the Prodigal
        training file exists or False otherwise and the path
        to the file if it exists or a message if it does not
        exist.
    """

    if os.path.isfile(ptf_path) is False:
        message = ('Cannot find specified Prodigal training file.'
                   '\nPlease provide a valid training file.\n\nYou '
                   'can create a training file for a species of '
                   'interest with the following command:\n  prodigal '
                   '-i <reference_genome> -t <training_file.trn> -p '
                   'single\n\nIt is strongly advised to provide a '
                   'high-quality and closed genome for the training '
                   'process.')
        return [False, message]
    else:
        return [True, ptf_path]
