import os
import pickle
import pandas as pd

from flash_cos import FlashCos


def prepare_ms2_lib(mode='p', mz_tol=0.025,
                    peak_intensity_power=0.5,
                    gnps_only=False,
                    new_path=None):
    """
    prepare ms2 db
    """

    db = []

    all_df = pd.read_pickle('//db/ms2db/all/all_ms2db_metadata.pkl')
    all_df = all_df[all_df['ion_mode'] == mode].reset_index(drop=True)

    if gnps_only:
        all_df = all_df[all_df['db'] == 'gnps'].reset_index(drop=True)

    for i, row in all_df.iterrows():

        db.append({
            'id': row['db_id'],
            'precursor_mz': float(row['recalc_prec_mz']),
            'peaks': row['peaks']
        })

    print('Number of spectra in the database:', len(db))

    print('initializing search engine')
    search_engine = FlashCos(max_ms2_tolerance_in_da=mz_tol * 1.005,
                             mz_index_step=0.0001,
                             peak_scale_k=0.0,
                             peak_intensity_power=peak_intensity_power)
    print('building index')
    search_engine.build_index(db,
                              max_indexed_mz=2000.,
                              precursor_ions_removal_da=-0.5,
                              noise_threshold=0.0,
                              min_ms2_difference_in_da=mz_tol * 2.02,
                              clean_spectra=True)

    # save as pickle
    with open(new_path, 'wb') as file:
        pickle.dump(search_engine, file)

    print(f"Pickle file saved to: {new_path}")
    return search_engine


if __name__ == "__main__":
    # prepare_ms2_lib(mode='p', mz_tol=0.025,
    #                 peak_intensity_power=0.5,
    #                 gnps_only=False,
    #                 new_path='data/indexed_lib_pos.pkl')
    # prepare_ms2_lib(mode='n', mz_tol=0.025,
    #                 peak_intensity_power=0.5,
    #                 gnps_only=False,
    #                 new_path='data/indexed_lib_neg.pkl')

    prepare_ms2_lib(mode='p', mz_tol=0.025,
                    peak_intensity_power=0.5,
                    gnps_only=True,
                    new_path='data/indexed_lib_gnps_pos.pkl')
