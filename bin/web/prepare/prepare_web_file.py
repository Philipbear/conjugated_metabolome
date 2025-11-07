import pandas as pd

'''
qry_id	qry_mz	ref_1_id	ref_1_score	ref_1_peak	ref_1_usage	ref_1_prec_frag_int	ref_1_prec_frag_water_loss_int	ref_1_prec_nl_int	ref_1_prec_nl_water_loss_int	
ref_2_id	ref_2_score	ref_2_peak	ref_2_usage	total_usage	
ref_1_inchi	ref_1_name	ref_1_formula	ref_1_mass	ref_1_cls_superclass	ref_1_cls_class	ref_1_cls_subclass	ref_1_np_superclass	ref_1_np_class	ref_1_np_pathway	
ref_2_inchi	ref_2_name	ref_2_formula	ref_2_mass	ref_2_cls_superclass	ref_2_cls_class	ref_2_cls_subclass	ref_2_np_superclass	ref_2_np_class	ref_2_np_pathway	
delta_mass
'''

def prepare_search_results(result_path='/home/shipei/projects/revcos/search/results/pos_stitched/pos_best_raw.tsv',
                           out_path='/home/shipei/projects/revcos/search/results/pos_stitched/pos_web.parquet'):

    print(f"Preparing search results from {result_path} to {out_path}")
    df = pd.read_csv(result_path, sep='\t', low_memory=False)
    print('All columns:')
    print(df.columns.tolist())

    df['ref_1_prec_int'] = df.apply(lambda x: max(x['ref_1_prec_frag_int'], x['ref_1_prec_frag_water_loss_int']), axis=1)
    
    df = df[df['ref_1_prec_int'] >= 0.05].reset_index(drop=True)

    # drop cols ref_1_prec_frag_int, ref_1_prec_frag_water_loss_int, ref_1_prec_nl_int, ref_1_prec_nl_water_loss_int, total_usage
    df = df.drop(columns=['ref_1_prec_frag_int', 'ref_1_prec_frag_water_loss_int', 'ref_1_prec_nl_int', 'ref_1_prec_nl_water_loss_int', 'total_usage']).reset_index(drop=True)

    df.loc[df['ref_2_score'] > 0.5, 'spec_spec'] = True
    df.loc[df['ref_2_score'] < 0.5, 'spec_spec'] = False
    
    # load metadata
    print("Loading metadata...")
    # metadata = pd.read_pickle('/Users/shipei/Documents/projects/conjugated_metabolome/db/ms2db/all/all_ms2db_metadata.pkl')
    metadata = pd.read_pickle('/home/shipei/projects/revcos/search/all_ms2db_metadata.pkl')
    db_id_to_db = metadata.set_index('db_id')['db'].to_dict()
    db_id_to_mass = metadata.set_index('db_id')['monoisotopic_mass'].to_dict()
    db_id_to_inchi14 = metadata.set_index('db_id')['inchikey_14'].to_dict()
    
    # add metadata to df
    df['ref_1_mass'] = df['ref_1_id'].map(db_id_to_mass)
    df['delta_mass'] = df['qry_mz'] - df['ref_1_mass'] - 1.007276
    # delta_mass, round to 2 decimal places
    df['delta_mass'] = df['delta_mass'].round(2)
    
    df['ref_1_db'] = df['ref_1_id'].map(db_id_to_db)
    df['ref_1_db'] = df['ref_1_db'].apply(lambda x: 1 if x == 'gnps' or x == 'mona' else 0)
    
    def _fill_ref_2_inchi(row):
        if pd.isnull(row['ref_2_id']):
            return str(row['delta_mass'])
        else:
            id = str(row['ref_2_id'])
            if id.endswith('.0'):
                id = id[:-2]
            return db_id_to_inchi14[id]
    def _fill_ref_2_db(row):
        if pd.isnull(row['ref_2_id']):
            return 0
        else:
            id = str(row['ref_2_id'])
            if id.endswith('.0'):
                id = id[:-2]
            db = db_id_to_db[id]
            if db == 'gnps' or db == 'mona':
                return 1
            else:
                return 0
    df['ref_2_inchi'] = df.apply(_fill_ref_2_inchi, axis=1)
    df['ref_2_db'] = df.apply(_fill_ref_2_db, axis=1)
        
    # sort by ref_1_inchi, ref_2_inchi, ref_1_db, ref_2_db, ref_1_prec_int, ref_1_usage, ref_1_score, ref_2_score
    df = df.sort_values(by=['ref_1_inchi', 'ref_2_inchi', 'ref_1_db', 'ref_2_db', 'ref_1_prec_int', 'ref_1_score',
                            'ref_1_usage', 'ref_2_score'], ascending=[True, True, False, False, False, False, False, False]).reset_index(drop=True)
    
    # dereplicate by dataset, ref_1_inchi and ref_2_inchi, keep the first
    df['dataset'] = df['qry_id'].apply(lambda x: x.split(':')[0])
    df = df.drop_duplicates(subset=['dataset', 'ref_1_inchi', 'ref_2_inchi'], keep='first').reset_index(drop=True)
    
    # First get the group counts and collect all qry_ids
    print("Grouping by ref_1_inchi and ref_2_inchi...")
    grouped = df.groupby(['ref_1_inchi', 'ref_2_inchi']).agg({
        'qry_id': ['first', 'count', list]  # first for primary, count for total, list for all qry_ids
    }).reset_index()

    # Flatten column names
    grouped.columns = ['ref_1_inchi', 'ref_2_inchi', 'primary_qry_id', 'count', 'all_qry_ids']

    # Use drop_duplicates for other columns
    agg_df = df.drop_duplicates(subset=['ref_1_inchi', 'ref_2_inchi'], keep='first')[
        ['ref_1_inchi', 'ref_2_inchi', 'qry_id', 'qry_mz', 'ref_1_id', 'ref_2_id', 'delta_mass', 'spec_spec']
    ].reset_index(drop=True)

    # Merge the grouped data
    df = agg_df.merge(grouped, on=['ref_1_inchi', 'ref_2_inchi'])

    # Create other_qry_ids by removing the primary qry_id from the list
    df['other_qry_ids'] = df.apply(lambda row: [qid for qid in row['all_qry_ids'] if qid != row['qry_id']], axis=1)

    # Drop unnecessary columns
    df = df.drop(columns=['all_qry_ids', 'primary_qry_id'])
    
    # for spec_spec, if False, then set ref_2_id to None
    df['ref_2_id'] = df.apply(lambda x: x['ref_2_id'] if x['spec_spec'] else None, axis=1)
    
    def remove_prefix(db_id):
        if pd.isnull(db_id):
            return db_id
        return db_id.replace('CCMSLIB000', '').replace('MSBNK-', '')
    
    df['ref_1_id'] = df['ref_1_id'].apply(remove_prefix)
    df['ref_2_id'] = df['ref_2_id'].apply(remove_prefix)

    # drop cols
    df = df.drop(columns=['spec_spec', 'ref_1_inchi', 'ref_2_inchi']).reset_index(drop=True)
    
    # save as parquet file
    df.to_parquet(out_path, index=False)  # qry_id, qry_mz, ref_1_id, ref_2_id, delta_mass, other_qry_ids, count
    

def prepare_ms2db_metadata(path='/Users/shipei/Documents/projects/conjugated_metabolome/db/ms2db/all/all_ms2db_metadata.pkl'):
    
    metadata = pd.read_pickle(path)
    
    metadata = metadata[['name', 'db', 'db_id', 'inchikey_14', 'monoisotopic_mass']]
    # change the column 'db' to categorical
    metadata['db'] = metadata['db'].astype('category')
    # change the column 'monoisotopic_mass' to 4 decimal places
    metadata['monoisotopic_mass'] = metadata['monoisotopic_mass'].round(4)
    
    # show the first 5 rows of gnps db
    print(metadata['db'].cat.categories)
    # for each category, show the first 5 rows
    for category in metadata['db'].cat.categories:
        print(f"Category: {category}")
        print(metadata[metadata['db'] == category].head(5))

    # verify if all gnps spectra, db_id starts with 'CCMSLIB000'
    print(metadata[metadata['db'] == 'gnps']['db_id'].str.startswith('CCMSLIB000').all())
    # verify if all massbank spectra, db_id starts with 'MSBNK-'
    print(metadata[metadata['db'] == 'massbank']['db_id'].str.startswith('MSBNK-').all())
    
    # some rules
    # 1. if db is 'gnps', then db_id removes 'CCMSLIB000'
    # 2. if db is 'massbank', then db_id removes 'MSBNK-'
    # 3. if db is mona or nist20, then db_id remains the same
    def remove_prefix(db_id):
        if pd.isnull(db_id):
            return db_id
        return db_id.replace('CCMSLIB000', '').replace('MSBNK-', '')
    metadata['db_id'] = metadata['db_id'].apply(remove_prefix)
    
    # check if all db_id are unique
    if metadata['db_id'].is_unique:
        print("All db_id are unique")
    
    # save as parquet file
    metadata.to_parquet('/Users/shipei/Documents/projects/conjugated_metabolome/bin/web/prepare/ms2db.parquet', index=False)
    

def compress_ms2db_metadata(path='/Users/shipei/Documents/projects/conjugated_metabolome/bin/web/prepare/ms2db.parquet'):
    """
    Compress the ms2db metadata file, by db ids that are in the result files, to reduce its size.
    """
    db = pd.read_parquet(path)

    # load the search results
    print("Loading search results to filter db_ids...")
    neg = pd.read_parquet('/Users/shipei/Documents/projects/conjugated_metabolome/bin/web/prepare/neg_web.parquet')
    pos = pd.read_parquet('/Users/shipei/Documents/projects/conjugated_metabolome/bin/web/prepare/pos_web.parquet')
    all_ids = set(neg['ref_1_id']).union(set(neg['ref_2_id'])).union(set(pos['ref_1_id'])).union(set(pos['ref_2_id']))
    del neg, pos
    
    print(f"Total unique db_ids in search results: {len(all_ids)}")
    # filter the db by db_id
    db = db[db['db_id'].isin(all_ids)].reset_index(drop=True)
    print(f"Total unique db_ids in ms2db metadata after filtering: {len(db)}")
    # save as parquet file
    db.to_parquet('/Users/shipei/Documents/projects/conjugated_metabolome/bin/web/prepare/ms2db.parquet', index=False)



if __name__ == "__main__":
     
    # prepare_ms2db_metadata() # local
    
    ######
    # server
    # prepare_search_results(result_path='/home/shipei/projects/revcos/search/results/neg_stitched/neg_best_raw.tsv',
    #                        out_path='/home/shipei/projects/revcos/search/results/neg_stitched/neg_web.parquet')
    
    # prepare_search_results(result_path='/home/shipei/projects/revcos/search/results/pos_stitched/pos_best_raw.tsv',
    #                        out_path='/home/shipei/projects/revcos/search/results/pos_stitched/pos_web.parquet')
    
    ##############    
    compress_ms2db_metadata(path='/Users/shipei/Documents/projects/conjugated_metabolome/bin/web/prepare/ms2db.parquet')