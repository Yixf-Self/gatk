import numpy as np
from typing import List, Set, Dict
from .interval import Interval
from .. import types
import logging
import csv
import os

_logger = logging.getLogger(__name__)


class IntervalListMetadata:
    def __init__(self, interval_list: List[Interval]):
        _logger.info("Generating intervals metadata...")
        self.interval_list = interval_list
        self.num_intervals = len(interval_list)
        self.contig_set = self._get_contig_set_from_interval_list(interval_list)
        self.contig_list = sorted(list(self.contig_set))
        self.num_contigs = len(self.contig_list)

        # map from contig to indices in the interval list
        self.contig_interval_indices: Dict[str, List[int]] = \
            {contig: [ti for ti in range(len(interval_list))
                      if interval_list[ti].contig == contig]
             for contig in self.contig_set}

        # number of intervals per contig
        self.t_j = np.asarray([len(self.contig_interval_indices[self.contig_list[j]])
                               for j in range(self.num_contigs)], dtype=types.big_uint)

    @staticmethod
    def _get_contig_set_from_interval_list(interval_list: List[Interval]) -> Set[str]:
        return {interval.contig for interval in interval_list}


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
        self._contig_map = {contig: j for j, contig in enumerate(contig_list)}

    def _assert_contig_exists(self, contig: str):
        assert contig in self._contig_map, \
            "Sample ({0}) does not have coverage metadata for contig ({1})".format(self.sample_name, contig)

    def get_contig_total_count(self, contig: str):
        self._assert_contig_exists(contig)
        return self.n_j[self._contig_map[contig]]

    def get_total_count(self):
        return self.n_total

    @staticmethod
    def generate_sample_coverage_metadata(sample_name,
                                          n_t: np.ndarray,
                                          interval_list_metadata: IntervalListMetadata):
        n_j = np.zeros((len(interval_list_metadata.contig_list),), dtype=types.big_uint)
        for j, contig in enumerate(interval_list_metadata.contig_list):
            n_j[j] = np.sum(n_t[interval_list_metadata.contig_interval_indices[contig]])
        return SampleCoverageMetadata(sample_name, n_j, interval_list_metadata.contig_list)


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
        self._contig_map = {contig: j for j, contig in enumerate(contig_list)}

    def _assert_contig_exists(self, contig: str):
        assert contig in self._contig_map, \
            "Sample ({0}) does not have ploidy metadata for contig ({1})".format(self.sample_name, contig)

    def get_contig_ploidy(self, contig: str):
        self._assert_contig_exists(contig)
        return self.ploidy_j[self._contig_map[contig]]

    def get_contig_ploidy_genotyping_quality(self, contig: str):
        self._assert_contig_exists(contig)
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
                                            interval_list_metadata: IntervalListMetadata) -> 'SampleReadDepthMetadata':
        assert sample_coverage_metadata.sample_name == sample_ploidy_metadata.sample_name
        assert interval_list_metadata.contig_list == sample_ploidy_metadata.contig_list
        sample_name = sample_ploidy_metadata.sample_name
        n_total = sample_coverage_metadata.n_total
        t_j = interval_list_metadata.t_j
        ploidy_j = sample_ploidy_metadata.ploidy_j
        effective_total_copies = float(np.sum(t_j * ploidy_j))
        read_depth = float(n_total) / effective_total_copies
        return SampleReadDepthMetadata(sample_name, read_depth)


class SampleMetadataCollection:
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


class SampleAlreadyInCollectionException(Exception):
    def __init__(self, msg):
        super().__init__(msg)
