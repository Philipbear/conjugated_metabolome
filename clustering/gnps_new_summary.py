import pandas as pd
import os

from tqdm import tqdm


def gnps_new_datasets():
    df = pd.read_csv('data/all_gnps_datasets.tsv', sep='\t', low_memory=False)
    df = df[df['dataset_id'] > 95251].reset_index(drop=True)

    df['metabolomics_data'] = df.apply(lambda x: _recognize_metabolomics_dataset(x), axis=1)
    df = df[df['metabolomics_data']].reset_index(drop=True)

    unique_instruments = df['instrument_resolved'].value_counts().reset_index()

    unique_instruments.to_csv('data/gnps_new_model_summary.tsv', sep='\t', index=False)


def assign_ms_type_and_filter():
    ms_df = pd.read_csv('data/gnps_new_ms_type.tsv', sep='\t')  # manually curated

    ms_df['instrument_resolved'] = ms_df['instrument_resolved'].apply(lambda x: x.strip())
    ms_df['MS_type'] = ms_df['MS_type'].apply(lambda x: x.strip())

    # dict from model to ms type
    model_to_ms = dict(zip(ms_df['instrument_resolved'], ms_df['MS_type']))

    df = pd.read_csv('data/gnps_new_datasets.tsv', sep='\t')
    df = df[df['dataset_id'] > 95251].reset_index(drop=True)

    df['metabolomics_data'] = df.apply(lambda x: _recognize_metabolomics_dataset(x), axis=1)
    df = df[df['metabolomics_data']].reset_index(drop=True)
    df['MS_type'] = df['instrument_resolved'].apply(lambda x: model_to_ms.get(x, ''))

    # filter
    df = df[(df['MS_type'].isin(['orbitrap', 'tof']))].reset_index(drop=True)

    df.to_csv('data/gnps_new_datasets_filtered.tsv', sep='\t', index=False)


def count_ms2():
    # run on server
    mgfs = os.listdir('new_msv_unclustered')
    mgfs = [x for x in mgfs if x.endswith('.mgf')]
    mgfs = [os.path.join('new_msv_unclustered', x) for x in mgfs]

    ms2_count = 0
    for mgf in tqdm(mgfs):
        with open(mgf) as f:
            lines = f.readlines()
            ms2_count += len([x for x in lines if x.strip() == 'BEGIN IONS'])

    print('Total number of MS2s:', ms2_count)



def _recognize_metabolomics_dataset(x):

    modif_str = x['modification_resolved'].lower() if pd.notnull(x['modification_resolved']) else ''

    if 'unimod' in modif_str or 'pride' in modif_str:
        return False

    metabolo_keywords = ['metabolomics', 'metabolome', 'metabolite']
    proteo_keywords = ['proteomics', 'proteome', 'protein', 'translational']

    title = x['title'].lower() if pd.notnull(x['title']) else ''
    description = x['description'].lower() if pd.notnull(x['description']) else ''
    keywords = x['keywords'].lower() if pd.notnull(x['keywords']) else ''
    all_description = title + ' ' + description + ' ' + keywords

    if any([keyword in all_description for keyword in metabolo_keywords]):
        return True
    elif any([keyword in all_description for keyword in proteo_keywords]):
        return False
    else:
        return True


if __name__ == '__main__':
    gnps_new_datasets()

    assign_ms_type_and_filter()

    # count_ms2()
