import numpy as np
from typing import List, Set, Dict
from .interval import Interval
from .. import types
import logging
import csv
import os

_logger = logging.getLogger(__name__)


class TargetsIntervalListMetadata:
    def __init__(self, targets_interval_list: List[Interval]):
        _logger.info("Generating targets metadata...")
        self.targets_interval_list = targets_interval_list
        self.num_targets = len(targets_interval_list)
        self.contig_set = self._get_contig_set_from_interval_list(targets_interval_list)
        self.contig_list = sorted(list(self.contig_set))
        self.num_contigs = len(self.contig_list)

        # map from contig to indices in the target list
        self.contig_target_indices: Dict[str, List[int]] = \
            {contig: [ti for ti in range(len(targets_interval_list))
                      if targets_interval_list[ti].contig == contig]
             for contig in self.contig_set}

        # number of targets per contig
        self.t_j = np.asarray([len(self.contig_target_indices[self.contig_list[j]])
                               for j in range(self.num_contigs)], dtype=types.big_uint)

    @staticmethod
    def _get_contig_set_from_interval_list(targets_interval_list: List[Interval]) -> Set[str]:
        return {target.contig for target in targets_interval_list}


class SampleCoverageMetadata:
    """ Represents essential metadata collected from a sample's coverage profile """
    def __init__(self,
                 sample_name: str,
                 n_j: np.ndarray,
                 contig_list: List[str]):
        """
        :param sample_name: a string identifier
        :param n_j: total read count per contig array
        :param contig_list: list of contigs
        """
        assert n_j.ndim == 1
        assert n_j.size == len(contig_list)

        self.sample_name = sample_name
        self.contig_list = contig_list

        # total count per contig
        self.n_j = n_j.astype(types.med_uint)

        # total count
        self.n_total = np.sum(self.n_j)
        self._contig_map: Dict[str, int] = {contig: j for j, contig in enumerate(contig_list)}

    def get_contig_total_count(self, contig: str):
        return self.n_j[self._contig_map[contig]]

    def get_total_count(self):
        return self.n_total

    @staticmethod
    def generate_sample_coverage_metadata(sample_name,
                                          n_t: np.ndarray,
                                          targets_metadata: TargetsIntervalListMetadata):
        n_j = np.zeros((len(targets_metadata.contig_list),), dtype=types.big_uint)
        for j, contig in enumerate(targets_metadata.contig_list):
            n_j[j] = np.sum(n_t[targets_metadata.contig_target_indices[contig]])
        return SampleCoverageMetadata(sample_name, n_j, targets_metadata.contig_list)


class SamplePloidyMetadata:
    """ Estimated contig ploidy for a sample """
    def __init__(self,
                 sample_name: str,
                 ploidy_j: np.ndarray,
                 ploidy_genotyping_quality_j: np.ndarray,
                 contig_list: List[str]):
        """
        :param sample_name: a string identifier
        :param ploidy_j: a vector of ploidy per contig
        :param contig_list: list of contigs
        """
        assert ploidy_j.ndim == 1
        assert ploidy_j.size == len(contig_list)
        assert ploidy_genotyping_quality_j.ndim == 1
        assert ploidy_genotyping_quality_j.size == len(contig_list)

        self.sample_name = sample_name
        self.contig_list = contig_list
        self.ploidy_j = ploidy_j.astype(types.small_uint)
        self.ploidy_genotyping_quality_j = ploidy_genotyping_quality_j.astype(types.floatX)
        self._contig_map: Dict[str, int] = {contig: j for j, contig in enumerate(contig_list)}

    def get_contig_ploidy(self, contig: str):
        return self.ploidy_j[self._contig_map[contig]]

    def get_contig_ploidy_genotyping_quality(self, contig: str):
        return self.ploidy_genotyping_quality_j[self._contig_map[contig]]


class SampleReadDepthMetadata:
    def __init__(self,
                 sample_name: str,
                 read_depth: float):
        self.sample_name = sample_name
        self.read_depth = read_depth

    def get_read_depth(self):
        return self.read_depth

    @staticmethod
    def generate_sample_read_depth_metadata(sample_coverage_metadata: SampleCoverageMetadata,
                                            sample_ploidy_metadata: SamplePloidyMetadata,
                                            targets_metadata: TargetsIntervalListMetadata) -> 'SampleReadDepthMetadata':
        assert sample_coverage_metadata.sample_name == sample_ploidy_metadata.sample_name
        assert targets_metadata.contig_list == sample_ploidy_metadata.contig_list
        sample_name = sample_ploidy_metadata.sample_name
        n_total = sample_coverage_metadata.n_total
        t_j = targets_metadata.t_j
        ploidy_j = sample_ploidy_metadata.ploidy_j
        effective_total_copies = float(np.sum(t_j * ploidy_j))
        read_depth = float(n_total) / effective_total_copies
        return SampleReadDepthMetadata(sample_name, read_depth)


class SampleMetadataCollection:
    # for input/output tables
    _sample_name_column_name = 'SAMPLE_NAME'
    _ploidy_contig_header_prefix = 'PLOIDY_CONTIG_'
    _ploidy_gq_contig_header_prefix = 'PLOIDY_GQ_CONTIG_'
    _read_depth_column_name = 'READ_DEPTH'
    _total_count_column_name = 'TOTAL_COUNT'
    _contig_header_prefix = ''

    def __init__(self):
        self.sample_coverage_metadata_dict: Dict[str, SampleCoverageMetadata] = dict()
        self.sample_ploidy_metadata_dict: Dict[str, SamplePloidyMetadata] = dict()
        self.sample_read_depth_metadata_dict: Dict[str, SampleReadDepthMetadata] = dict()

    def add_sample_coverage_metadata(self, sample_coverage_metadata: SampleCoverageMetadata):
        sample_name = sample_coverage_metadata.sample_name
        if sample_name in self.sample_coverage_metadata_dict.keys():
            raise SampleAlreadyInCollectionException(
                'sample "{0}" already has coverage metadata annotations'.format(sample_name))
        else:
            self.sample_coverage_metadata_dict[sample_name] = sample_coverage_metadata

    def add_sample_ploidy_metadata(self, sample_ploidy_metadata: SamplePloidyMetadata):
        sample_name = sample_ploidy_metadata.sample_name
        if sample_name in self.sample_ploidy_metadata_dict.keys():
            raise SampleAlreadyInCollectionException(
                'sample "{0}" already has ploidy metadata annotations'.format(sample_name))
        else:
            self.sample_ploidy_metadata_dict[sample_name] = sample_ploidy_metadata

    def add_sample_read_depth_metadata(self, sample_read_depth_metadata: SampleReadDepthMetadata):
        sample_name = sample_read_depth_metadata.sample_name
        if sample_name in self.sample_read_depth_metadata_dict.keys():
            raise SampleAlreadyInCollectionException(
                'sample "{0}" already has read depth metadata annotations'.format(sample_name))
        else:
            self.sample_read_depth_metadata_dict[sample_name] = sample_read_depth_metadata

    def all_samples_have_coverage_metadata(self, sample_names: List[str]):
        return all([sample_name in self.sample_coverage_metadata_dict.keys()
                    for sample_name in sample_names])

    def all_samples_have_ploidy_metadata(self, sample_names: List[str]):
        return all([sample_name in self.sample_ploidy_metadata_dict.keys()
                    for sample_name in sample_names])

    def all_samples_have_read_depth_metadata(self, sample_names: List[str]):
        return all([sample_name in self.sample_read_depth_metadata_dict.keys()
                    for sample_name in sample_names])

    def get_sample_coverage_metadata(self, sample_name: str) -> SampleCoverageMetadata:
        return self.sample_coverage_metadata_dict[sample_name]

    def get_sample_ploidy_metadata(self, sample_name: str) -> SamplePloidyMetadata:
        return self.sample_ploidy_metadata_dict[sample_name]

    def get_sample_read_depth_metadata(self, sample_name: str) -> SampleReadDepthMetadata:
        return self.sample_read_depth_metadata_dict[sample_name]

    def get_sample_read_depth_array(self, sample_names: List[str]) -> np.ndarray:
        return np.asarray([self.sample_read_depth_metadata_dict[sample_name].get_read_depth()
                           for sample_name in sample_names], dtype=types.floatX)

    def get_sample_contig_ploidy_array(self, contig: str, sample_names: List[str]) -> np.ndarray:
        return np.asarray([self.get_sample_ploidy_metadata(sample_name).get_contig_ploidy(contig)
                           for sample_name in sample_names], dtype=types.small_uint)

    def write_sample_coverage_metadata(self,
                                       sample_names: List[str],
                                       output_file: str):
        from ..utils.io import assert_output_path_writable
        assert len(sample_names) > 0
        assert self.all_samples_have_coverage_metadata(sample_names)
        contig_list = self.sample_coverage_metadata_dict[sample_names[0]].contig_list
        for sample_name in sample_names:
            assert self.sample_coverage_metadata_dict[sample_name].contig_list == contig_list
        parent_path = os.path.dirname(output_file)
        assert_output_path_writable(parent_path)
        with open(output_file, 'w') as tsv_file:
            writer = csv.writer(tsv_file, delimiter='\t')
            header = ([self._sample_name_column_name]
                      + [self._contig_header_prefix + contig for contig in contig_list])
            writer.writerow(header)
            for sample_name in sample_names:
                sample_coverage_metadata = self.get_sample_coverage_metadata(sample_name)
                row = ([sample_name] + [repr(sample_coverage_metadata.n_j[j]) for j in range(len(contig_list))])
                writer.writerow(row)

    def write_sample_read_depth_metadata(self,
                                         sample_names: List[str],
                                         output_file: str):
        from ..utils.io import assert_output_path_writable
        assert len(sample_names) > 0
        assert self.all_samples_have_read_depth_metadata(sample_names)
        parent_path = os.path.dirname(output_file)
        assert_output_path_writable(parent_path)
        with open(output_file, 'w') as tsv_file:
            writer = csv.writer(tsv_file, delimiter='\t')
            header = [self._sample_name_column_name, self._read_depth_column_name]
            writer.writerow(header)
            for sample_name in sample_names:
                sample_read_depth_metadata = self.get_sample_read_depth_metadata(sample_name)
                row = [sample_name, repr(sample_read_depth_metadata.get_read_depth())]
                writer.writerow(row)

    def write_sample_contig_ploidy_metadata(self,
                                            sample_names: List[str],
                                            output_file: str):
        from ..utils.io import assert_output_path_writable
        assert len(sample_names) > 0
        assert self.all_samples_have_ploidy_metadata(sample_names)
        contig_list = self.sample_ploidy_metadata_dict[sample_names[0]].contig_list
        for sample_name in sample_names:
            assert self.sample_ploidy_metadata_dict[sample_name].contig_list == contig_list

        parent_path = os.path.dirname(output_file)
        assert_output_path_writable(parent_path)
        with open(output_file, 'w') as tsv_file:
            writer = csv.writer(tsv_file, delimiter='\t')
            header = ([self._sample_name_column_name]
                      + [self._ploidy_contig_header_prefix + contig for contig in contig_list]
                      + [self._ploidy_gq_contig_header_prefix + contig for contig in contig_list])
            writer.writerow(header)
            for sample_name in sample_names:
                sample_ploidy_metadata = self.get_sample_ploidy_metadata(sample_name)
                row = ([sample_name]
                       + [repr(sample_ploidy_metadata.ploidy_j[j]) for j in range(len(contig_list))]
                       + [repr(sample_ploidy_metadata.ploidy_genotyping_quality_j[j]) for j in range(len(contig_list))])
                writer.writerow(row)

    def read_sample_coverage_metadata(self, input_file: str) -> List[str]:
        sample_names = []
        with open(input_file, 'r') as tsv_file:
            reader = csv.reader(tsv_file, delimiter='\t')
            row_num = 0
            contig_list = []
            for row in reader:
                row_num += 1
                num_header_elems = None
                if row_num == 1:  # header
                    num_header_elems = len(row)
                    assert num_header_elems > 1, "malformed sample coverage metadata file"
                    assert row[0] == self._sample_name_column_name, "malformed sample ploidy metadata file"
                    num_contigs = num_header_elems - 1
                    assert all(
                        len(row[k + 1]) > len(self._contig_header_prefix)
                        and row[k + 1][:len(self._contig_header_prefix)] == self._contig_header_prefix
                        for k in range(num_contigs)), "malformed sample ploidy metadata file"
                    for k in range(num_contigs):
                        contig_list.append(row[k + 1][len(self._contig_header_prefix):])
                else:
                    assert len(row) == num_header_elems
                    sample_name = row[0]
                    n_j = np.asarray([int(row[k + 1]) for k in range(num_contigs)], dtype=types.big_uint)
                    self.add_sample_coverage_metadata(SampleCoverageMetadata(sample_name, n_j, contig_list))
                    sample_names.append(sample_name)
        return sample_names

    def read_sample_read_depth_metadata(self, input_file: str):
        with open(input_file, 'r') as tsv_file:
            reader = csv.reader(tsv_file, delimiter='\t')
            expected_header = [self._sample_name_column_name, self._read_depth_column_name]
            row_num = 0
            for row in reader:
                row_num += 1
                if row_num == 1:
                    assert row == expected_header, "malformed sample read depth metadata file"
                    continue
                sample_name = row[0]
                read_depth = float(row[1])
                self.add_sample_read_depth_metadata(SampleReadDepthMetadata(sample_name, read_depth))

    def read_sample_ploidy_metadata(self, input_file: str):
        with open(input_file, 'r') as tsv_file:
            reader = csv.reader(tsv_file, delimiter='\t')
            row_num = 0
            contig_list = []
            for row in reader:
                row_num += 1
                if row_num == 1:
                    assert row[0] == self._sample_name_column_name, "malformed sample ploidy metadata file"
                    num_header_elems = len(row)
                    assert num_header_elems % 2 == 1, "malformed sample ploidy metadata file"
                    num_contigs = (num_header_elems - 1) // 2
                    assert all(
                        len(row[k + 1]) > len(self._ploidy_contig_header_prefix)
                        and row[k + 1][:len(self._ploidy_contig_header_prefix)] == self._ploidy_contig_header_prefix
                        for k in range(num_contigs)), "malformed sample ploidy metadata file"
                    assert all(
                        len(row[k + num_contigs + 1]) > len(self._ploidy_gq_contig_header_prefix)
                        and row[k + num_contigs + 1][:len(self._ploidy_gq_contig_header_prefix)] ==
                        self._ploidy_gq_contig_header_prefix
                        for k in range(num_contigs)), "malformed sample ploidy metadata file"
                    for k in range(num_contigs):
                        contig_list.append(row[k + 1][len(self._ploidy_contig_header_prefix):])
                    continue

                sample_name = row[0]
                ploidy_j = np.asarray([int(row[k + 1]) for k in range(num_contigs)], dtype=types.small_uint)
                ploidy_gq_j = np.asarray([float(row[k + num_contigs + 1]) for k in range(num_contigs)],
                                         dtype=types.floatX)
                self.add_sample_ploidy_metadata(SamplePloidyMetadata(sample_name, ploidy_j, ploidy_gq_j, contig_list))

    @staticmethod
    def read_sample_names(input_file: str) -> List[str]:
        sample_names = []
        with open(input_file, 'r') as tsv_file:
            reader = csv.reader(tsv_file, delimiter='\t')
            row_num = 0
            for row in reader:
                row_num += 1
                if row_num == 1:  # header
                    num_header_elems = len(row)
                    assert num_header_elems == 1, "malformed sample names file"
                    assert row[0] == SampleMetadataCollection._sample_name_column_name, "malformed sample names file"
                else:
                    assert len(row) == 1, "malformed sample names file"
                    sample_names.append(row[0])
        return sample_names


class SampleAlreadyInCollectionException(Exception):
    def __init__(self, msg):
        super().__init__(msg)
