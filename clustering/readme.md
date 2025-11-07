
File summary:
`file_summary.py`, `gnps_new_summary.py`

----------
`get_ms2_from_raw_data.py`
- on minglab's server: Generate a mgf file for each dataset (only singly charged spectra) and transfer them to Dlab server


`gnps_download_merge.py`, `get_ms2_from_raw_data.py`
- on Dlab's server: Download new GNPS datasets and generate mgf files


> Why to generate mgf files for each dataset? 
1. read polarity correctly for each spectrum 
2. easy to split by mz values for clustering 
3. some mzML files may corrupt

----------
`ms2_summary.py`
- Generate MS/MS summary for all datasets 
- Re-split all mgf files according to mz values (easy for clustering, `unclustered_mzgrouped`)

----------
`falcon_unclustered_mgf.py`
- Run falcon on all unclustered mgf files
- Generate clustered mgf files & pkl files (dict of representative cluster_id: [all cluster_ids])