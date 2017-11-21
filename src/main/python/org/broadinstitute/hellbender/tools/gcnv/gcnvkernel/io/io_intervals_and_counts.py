import numpy as np
import pandas as pd
import logging
from typing import Optional, List, Tuple, Set

from ..structs.interval import Interval, IntervalAnnotation, interval_annotations_dtypes, interval_annotations_dict
from .. import types
from . import io_commons
from . import io_consts

interval_dtypes_dict = {
    io_consts.contig_column_header: np.str,
    io_consts.start_column_header: types.med_uint,
    io_consts.end_column_header: types.med_uint
}

read_count_dtypes_dict = {
    **interval_dtypes_dict,
    io_consts.count_column_header: types.med_uint
}

_logger = logging.getLogger(__name__)


def load_read_counts_tsv_file(read_counts_tsv_file: str,
                              max_rows: Optional[int] = None,
                              output_targets: bool = False) \
        -> Tuple[str, np.ndarray, Optional[List[Interval]]]:
    sample_name = io_commons.extract_sample_name_from_header(read_counts_tsv_file)
    counts_pd = pd.read_csv(read_counts_tsv_file, delimiter='\t', comment='#', nrows=max_rows,
                            dtype={**read_count_dtypes_dict})
    if output_targets:
        targets_pd = counts_pd[list(interval_dtypes_dict.keys())]
        targets_interval_list = _convert_targets_pd_to_interval_list(targets_pd)
        return sample_name, counts_pd[io_consts.count_column_header].as_matrix(), targets_interval_list
    else:
        return sample_name, counts_pd[io_consts.count_column_header].as_matrix(), None


def load_interval_list_tsv_file(targets_tsv_file: str) -> List[Interval]:
    targets_pd = pd.read_csv(targets_tsv_file, delimiter='\t',
                             dtype={**interval_dtypes_dict, **interval_annotations_dtypes})
    return _convert_targets_pd_to_interval_list(targets_pd)


def load_counts_in_the_modeling_zone(read_count_file_list: List[str],
                                     modeling_interval_list: List[Interval]):
    """ Note: it is assumed that all read counts have the same intervals; this is not asserted for speed """
    num_targets = len(modeling_interval_list)
    num_samples = len(read_count_file_list)
    assert num_samples > 0
    assert num_targets > 0

    sample_names: List[str] = []
    n_st = np.zeros((num_samples, num_targets), dtype=types.med_uint)
    master_interval_list = None
    interval_to_index_map = None
    for si, read_count_file in enumerate(read_count_file_list):
        if master_interval_list is None:
            sample_name, n_t, master_interval_list = load_read_counts_tsv_file(read_count_file, output_targets=True)
            interval_to_index_map = {interval: ti for ti, interval in enumerate(master_interval_list)}
            assert all([interval in interval_to_index_map for interval in modeling_interval_list]), \
                "Some of the modeling intervals are absent in the provided read counts .tsv file"
        else:  # we do not load targets again for speed, assume it is the same as the first sample
            sample_name, n_t, _ = load_read_counts_tsv_file(read_count_file, output_targets=False)
        n_st[si, :] = np.asarray([n_t[interval_to_index_map[interval]]
                                  for interval in modeling_interval_list], dtype=types.med_uint)
        sample_names.append(sample_name)
    return sample_names, n_st


def _convert_targets_pd_to_interval_list(targets_pd: pd.DataFrame) -> List[Interval]:
    """
    Converts a pandas dataframe targets intervals to list(Interval). Annotations will be parsed
    and added to the intervals as well.
    """
    interval_list: List[Interval] = []
    columns = [str(x) for x in targets_pd.columns.values]
    assert all([required_column in columns
                for required_column in interval_dtypes_dict.keys()]), "Some columns missing"
    for contig, start, end in zip(targets_pd[io_consts.contig_column_header],
                                  targets_pd[io_consts.start_column_header],
                                  targets_pd[io_consts.end_column_header]):
        interval = Interval(contig, start, end)
        interval_list.append(interval)

    for annotation_key in set(columns).intersection(interval_annotations_dict.keys()):
        bad_annotations_found = False
        for ti, raw_value in enumerate(targets_pd[annotation_key]):
            try:
                annotation: IntervalAnnotation = interval_annotations_dict[annotation_key](raw_value)
                interval_list[ti].add_annotation(annotation_key, annotation)
            except ValueError:
                bad_annotations_found = True
        if bad_annotations_found:
            _logger.warning("Some of the annotations for {0} contained bad values and were ignored".format(
                annotation_key))

    return interval_list


def write_interval_list_to_tsv_file(output_file: str, interval_list: List[Interval]):
    assert len(interval_list) > 0, "can not write an empty interval list to disk"
    annotation_found_keys: Set[str] = set()
    for interval in interval_list:
        for key in interval.annotations.keys():
            annotation_found_keys.add(key)
    mutual_annotation_key_list: List[str] = []
    for key in annotation_found_keys:
        if all(key in interval.annotations.keys() for interval in interval_list):
            mutual_annotation_key_list.append(key)
        else:
            _logger.warning("Only some targets have annotation \"{0}\" and others do not; "
                            "cannot write this annotation to disk; proceeding...")
    with open(output_file, 'w') as out:
        header = '\t'.join([io_consts.contig_column_header,
                            io_consts.start_column_header,
                            io_consts.end_column_header]
                           + mutual_annotation_key_list)
        out.write(header + '\n')
        for interval in interval_list:
            row = '\t'.join([interval.contig, repr(interval.start), repr(interval.end)] +
                            [str(interval.annotations[key]) for key in mutual_annotation_key_list])
            out.write(row + '\n')
