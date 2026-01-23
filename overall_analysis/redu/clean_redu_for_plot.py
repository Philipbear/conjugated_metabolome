import pickle
import pandas as pd
from collections import Counter



def clean_redu():
    redu = pd.read_csv('data/all_sampleinformation.tsv', sep='\t', low_memory=False)
    redu = redu[['ATTRIBUTE_DatasetAccession', 'filename', 'SampleType', 'NCBITaxonomy', 'NCBIDivision', 'UBERONBodyPartName', 'HealthStatus']]

    redu['filename'] = redu['filename'].apply(lambda x: x.split('/')[-1].split('.mz')[0])

    redu['mri'] = redu['ATTRIBUTE_DatasetAccession'] + ':' + redu['filename']

    # drop cols
    redu = redu.drop(columns=['ATTRIBUTE_DatasetAccession', 'filename'])

    # save
    redu.to_csv('data/redu.tsv', sep='\t', index=False)


def load_data():
    with open('data/pos_all_mri.pkl', 'rb') as f:
        pos_all = pickle.load(f)
    with open('data/neg_all_mri.pkl', 'rb') as f:
        neg_all = pickle.load(f)
    with open('data/pos_annotated.pkl', 'rb') as f:
        pos_annotated = pickle.load(f)
    with open('data/neg_annotated.pkl', 'rb') as f:
        neg_annotated = pickle.load(f)

    print('All data loaded')
    return pos_all, neg_all, pos_annotated, neg_annotated


def merge_data():
    # Load the data
    redu = pd.read_csv('data/redu.tsv', sep='\t', low_memory=False)

    pos_all, neg_all, pos_annotated, neg_annotated = load_data()

    # Pre-compute counts using Counter objects (much faster than manual counting)
    print("Creating count dictionaries...")
    pos_all_counts = Counter(pos_all)
    neg_all_counts = Counter(neg_all)
    pos_annotated_counts = Counter(pos_annotated)
    neg_annotated_counts = Counter(neg_annotated)

    # Create combined count dictionaries
    all_counts = pos_all_counts + neg_all_counts
    annotated_counts = pos_annotated_counts + neg_annotated_counts

    print("Applying counts to dataframe...")
    redu['all'] = redu['mri'].map(all_counts).fillna(0).astype(int)
    redu['annotated'] = redu['mri'].map(annotated_counts).fillna(0).astype(int)

    # Save
    print("Saving results...")
    redu.to_csv('data/redu_with_data.tsv', sep='\t', index=False)


if __name__ == '__main__':
    import os
    os.chdir(os.path.dirname(__file__))

    clean_redu()
    merge_data()
