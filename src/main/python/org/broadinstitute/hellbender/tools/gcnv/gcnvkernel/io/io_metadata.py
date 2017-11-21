import os
import csv
import numpy as np
from typing import List
import pandas as pd
import logging

from ..structs.metadata import SampleReadDepthMetadata, SamplePloidyMetadata, SampleCoverageMetadata,\
    SampleMetadataCollection
from .. import types
from . import io_commons
from . import io_consts

_logger = logging.getLogger(__name__)


def write_sample_coverage_metadata(sample_metadata_collection: SampleMetadataCollection,
                                   sample_names: List[str],
                                   output_file: str):
    assert len(sample_names) > 0
    assert sample_metadata_collection.all_samples_have_coverage_metadata(sample_names)
    contig_list = sample_metadata_collection.sample_coverage_metadata_dict[sample_names[0]].contig_list
    for sample_name in sample_names:
        assert sample_metadata_collection.sample_coverage_metadata_dict[sample_name].contig_list == contig_list
    parent_path = os.path.dirname(output_file)
    io_commons.assert_output_path_writable(parent_path)
    with open(output_file, 'w') as tsv_file:
        writer = csv.writer(tsv_file, delimiter='\t')
        header = [io_consts.sample_name_column_name] + [contig for contig in contig_list]
        writer.writerow(header)
        for sample_name in sample_names:
            sample_coverage_metadata = sample_metadata_collection.get_sample_coverage_metadata(sample_name)
            row = ([sample_name] + [repr(sample_coverage_metadata.n_j[j]) for j in range(len(contig_list))])
            writer.writerow(row)


def read_sample_coverage_metadata(sample_metadata_collection: SampleMetadataCollection,
                                  input_file: str) -> List[str]:
    with open(input_file, 'r') as tsv_file:
        reader = csv.reader(tsv_file, delimiter='\t')
        row_num = 0
        contig_list = None
        sample_names = []
        for row in reader:
            row_num += 1
            if row_num == 1:  # header
                num_header_elems = len(row)
                assert num_header_elems > 1, "malformed sample coverage metadata file"
                assert row[0] == io_consts.sample_name_column_name, "malformed sample ploidy metadata file"
                num_contigs = num_header_elems - 1
                contig_list = row[1:]
                continue

            assert len(row) == num_header_elems
            sample_name = row[0]
            n_j = np.asarray([int(row[k + 1]) for k in range(num_contigs)], dtype=types.big_uint)
            sample_metadata_collection.add_sample_coverage_metadata(SampleCoverageMetadata(
                sample_name, n_j, contig_list))
            sample_names.append(sample_name)

    return sample_names


def update_sample_metadata_collection_from_ploidy_determination_calls(
        sample_metadata_collection: SampleMetadataCollection,
        input_calls_path: str,
        comment='#',
        delimiter='\t'):

    def get_sample_name(input_path: str) -> str:
        sample_name_file = os.path.join(input_path, io_consts.default_sample_name_txt_filename)
        assert os.path.exists(sample_name_file), \
            "Sample name could not be found in the ploidy results located at \"{0}\"".format(input_path)
        with open(sample_name_file, 'r') as f:
            for line in f:
                return line.strip()

    def get_sample_read_depth_metadata(input_path: str) -> SampleReadDepthMetadata:
        sample_read_depth_file = os.path.join(input_path, io_consts.default_sample_read_depth_tsv_filename)
        assert os.path.exists(sample_read_depth_file), \
            "Sample read depth could not be found in the ploidy results located at \"{0}\"".format(input_path)
        _sample_name = io_commons.extract_sample_name_from_header(sample_read_depth_file)
        sample_read_depth_pd = pd.read_csv(sample_read_depth_file, delimiter=delimiter, comment=comment)
        assert io_consts.global_read_depth_column_name in sample_read_depth_pd.columns.values,\
            "Read depth file \"{0}\" does not contain the mandatory column \"{1}\"".format(
                sample_read_depth_file, io_consts.global_read_depth_column_name)
        read_depth = sample_read_depth_pd[io_consts.global_read_depth_column_name].values[0]
        return SampleReadDepthMetadata(_sample_name, read_depth)

    def get_sample_ploidy_metadata(input_path: str) -> SamplePloidyMetadata:
        sample_ploidy_file = os.path.join(input_path, io_consts.default_sample_contig_ploidy_tsv_filename)
        assert os.path.exists(sample_ploidy_file), \
            "Sample ploidy results could not be found in the ploidy results located at \"{0}\"".format(input_path)
        _sample_name = io_commons.extract_sample_name_from_header(sample_ploidy_file)
        sample_ploidy_pd = pd.read_csv(sample_ploidy_file, delimiter=delimiter, comment=comment)
        for mandatory_column in [io_consts.contig_column_name,
                                 io_consts.ploidy_column_name,
                                 io_consts.ploidy_gq_column_name]:
            assert mandatory_column in sample_ploidy_pd.columns.values,\
                "Sample contig ploidy file \"{0}\" does not contain the mandatory column \"{1}\"".format(
                    sample_ploidy_file, mandatory_column)
        contig_list = [str(x) for x in sample_ploidy_pd[io_consts.contig_column_name].values]
        ploidy_list = [int(x) for x in sample_ploidy_pd[io_consts.ploidy_column_name].values]
        ploidy_gq_list = [float(x) for x in sample_ploidy_pd[io_consts.ploidy_gq_column_name].values]
        return SamplePloidyMetadata(_sample_name,
                                    np.asarray(ploidy_list, dtype=types.small_uint),
                                    np.asarray(ploidy_gq_list, dtype=types.floatX),
                                    contig_list)

    _logger.info("Loading germline contig ploidy and global read depth metadata...")
    assert os.path.exists(input_calls_path) and os.path.isdir(input_calls_path), \
        "The provided path to ploidy determination results \"{0}\" is not a directory".format(input_calls_path)
    subdirs = os.listdir(input_calls_path)
    for subdir in subdirs:
        if subdir.find(io_consts.sample_folder_prefix) >= 0:
            sample_ploidy_results_dir = os.path.join(input_calls_path, subdir)
            sample_name = get_sample_name(sample_ploidy_results_dir)
            sample_read_depth_metadata = get_sample_read_depth_metadata(sample_ploidy_results_dir)
            sample_ploidy_metadata = get_sample_ploidy_metadata(sample_ploidy_results_dir)
            assert sample_read_depth_metadata.sample_name == sample_name, \
                "Inconsistency detected in the ploidy determination results; cannot continue"
            assert sample_ploidy_metadata.sample_name == sample_name, \
                "Inconsistency detected in the ploidy determination results; cannot continue"
            sample_metadata_collection.add_sample_read_depth_metadata(sample_read_depth_metadata)
            sample_metadata_collection.add_sample_ploidy_metadata(sample_ploidy_metadata)
