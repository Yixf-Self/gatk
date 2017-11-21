contig_column_header = 'CONTIG'
start_column_header = 'START'
end_column_header = 'END'
count_column_header = 'COUNT'
sample_folder_prefix = "SAMPLE_"
copy_number_column_prefix = "COPY_NUMBER_"
sample_name_column_name = 'SAMPLE_NAME'
global_read_depth_column_name = 'GLOBAL_READ_DEPTH'
contig_column_name = 'CONTIG'
ploidy_column_name = 'PLOIDY'
ploidy_gq_column_name = 'PLOIDY_GQ'
sample_name_header_regexp = "^[# ]+SAMPLE_NAME[\\s\t]*=[\\s\t]*(.*\\S)[\\s\t]*$"
sample_name_header_prefix = "SAMPLE_NAME="

default_sample_read_depth_tsv_filename = 'global_read_depth.tsv'
default_sample_name_txt_filename = 'sample_name.txt'
default_sample_contig_ploidy_tsv_filename = 'contig_ploidy.tsv'
default_copy_number_log_posterior_tsv_filename = "log_q_c_tc.tsv"
default_class_log_posterior_tsv_filename = "log_q_tau_tk.tsv"

default_denoising_config_json_filename = "denoising_config.json"
default_calling_config_json_filename = "calling_config.json"
default_ploidy_config_json_filename = "ploidy_config.json"
default_gcnvkernel_version_json_filename = "gcnvkernel_version.json"
default_contig_ploidy_prior_tsv_filename = 'contig_ploidy_prior.tsv'
default_interval_list_filename = "interval_list.tsv"
