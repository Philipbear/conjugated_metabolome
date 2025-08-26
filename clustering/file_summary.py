"""
file summaries obtained from Yasin. fill in the polarities by Shipei (get_file_polarity_centroid_server.py)
select files with MS2 data, fill in ms instrument type.
"""
import pickle

import pandas as pd


def print_stats():
    all_models = []
    df = pd.read_csv('data/mwb_results_with_polarity_centroided.tsv', sep='\t')
    all_models.extend(df['Model'].tolist())
    print(df['Polarity'].value_counts())
    print(df['Centroided'].value_counts())
    print(df['Model'].value_counts())

    df = pd.read_csv('data/mtbls_results_with_polarity_centroided.tsv', sep='\t')
    all_models.extend(df['Model'].tolist())
    print(df['Polarity'].value_counts())
    print(df['Centroided'].value_counts())
    print(df['Model'].value_counts())

    df = pd.read_csv('data/msv_results_with_polarity_centroided.tsv', sep='\t')
    all_models.extend(df['Model'].tolist())
    print(df['Polarity'].value_counts())
    print(df['Centroided'].value_counts())
    print(df['Model'].value_counts())

    print(pd.Series(all_models).value_counts()[:20])


def merge_and_filter_file_summary():
    df1 = pd.read_csv('data/mwb_results_with_polarity_centroided.tsv', sep='\t')
    df2 = pd.read_csv('data/mtbls_results_with_polarity_centroided.tsv', sep='\t')
    df3 = pd.read_csv('data/msv_results_with_polarity_centroided.tsv', sep='\t')
    df = pd.concat([df1, df2, df3])

    print('Total number of files:', df.shape[0])

    # generate summary of the models as a csv file
    model_summary = df['Model'].value_counts().reset_index()

    model_summary.to_csv('data/model_summary.tsv', sep='\t', index=False)
    df.to_csv('data/all_with_polarity.tsv', sep='\t', index=False)


def assign_ms_type():
    ms_df = pd.read_csv('data/ms_type.csv')  # manually curated

    ms_df['Model'] = ms_df['Model'].apply(lambda x: x.strip())
    ms_df['MS_type'] = ms_df['MS_type'].apply(lambda x: x.strip())

    # dict from model to ms type
    model_to_ms = dict(zip(ms_df['Model'], ms_df['MS_type']))

    df = pd.read_csv('data/all_with_polarity.tsv', sep='\t')
    df['MS_type'] = df['Model'].apply(lambda x: model_to_ms.get(x, ''))

    df.to_csv('data/all_with_polarity_and_ms_type.tsv', sep='\t', index=False)


def filter_data_for_clustering():

    df = pd.read_csv('data/all_with_polarity_and_ms_type.tsv', sep='\t')

    df['Repo'] = df['Filename'].apply(lambda x: x.split('/')[4])
    df['Dataset'] = df['Filename'].apply(lambda x: x.split('/')[5])

    print('Total number of files:', df.shape[0])
    print('Total number of datasets:', df['Dataset'].nunique())
    print('Total number of MS2s:', df['MS2s'].sum())

    # filter
    df = df[(df['MS2s'] > 0) &
            (df['MS_type'].isin(['orbitrap', 'tof', 'fticr'])) &
            (df['Centroided'] == 'centroided') &
            (~df['MS_type'].str.contains('gc'))].reset_index(drop=True)

    # all gnps datasets
    gnps_df = pd.read_csv('data/all_gnps_datasets.tsv', sep='\t', low_memory=False)
    gnps_df['metabolomics_data'] = gnps_df.apply(lambda x: _recognize_metabolomics_dataset(x), axis=1)
    proteomics_datasets = gnps_df[~gnps_df['metabolomics_data']]['dataset'].tolist()
    df = df[~df['Dataset'].isin(proteomics_datasets)].reset_index(drop=True)

    df.to_csv('data/all_filtered.tsv', sep='\t', index=False)

    print('After filtering:')
    print('Total number of files:', df.shape[0])
    print('Total number of datasets:', df['Dataset'].nunique())
    print('Total number of MS2s:', df['MS2s'].sum())
    print('MS types:', df['MS_type'].value_counts())
    # repo: repos of unique datasets
    df = df.drop_duplicates(subset='Dataset')
    print('Repos:', df['Repo'].value_counts())


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


def _recognize_dia_dataset(x):

    dia_keywords = ['(dia)', '(dia ', ' dia)', ' dia ', 'data-independent', 'data independent']

    title = x['title'].lower() if pd.notnull(x['title']) else ''
    description = x['description'].lower() if pd.notnull(x['description']) else ''
    keywords = x['keywords'].lower() if pd.notnull(x['keywords']) else ''
    all_description = title + ' ' + description + ' ' + keywords

    if any([keyword in all_description for keyword in dia_keywords]):
        return True

    return False


def recognize_dia_datasets():

    gnps_df = pd.read_csv('data/all_gnps_datasets.tsv', sep='\t', low_memory=False)
    gnps_df['dia_data'] = gnps_df.apply(lambda x: _recognize_dia_dataset(x), axis=1)

    dia_datasets = gnps_df[gnps_df['dia_data']]['dataset'].tolist()
    print('GNPS DIA datasets:', dia_datasets)

    # save as pickle
    with open('data/dia_datasets.pkl', 'wb') as f:
        pickle.dump(dia_datasets, f)


if __name__ == '__main__':

    # print_stats()
    #
    # merge_and_filter_file_summary()
    #
    # # perform MS type determination on Claude 3.5
    #
    # assign_ms_type()

    # filter_data_for_clustering()
    '''
    Total number of files: 632919
    Total number of datasets: 3059
    Total number of MS2s: 1213366329
    
    After filtering:
    Total number of files: 200820
    Total number of datasets: 1248
    Total number of MS2s: 661981944
    MS types: MS_type
    orbitrap    189120
    tof          11472
    fticr          228
    Name: count, dtype: int64
    Repos: Repo
    MassIVE    1107
    MTBLS        80
    ST           61
    
    
    GNPS new datasets data:
    281 datasets
    
    After filtering:
    253 datasets
    Total number of MS2s: 110519249
    '''

    recognize_dia_datasets()
