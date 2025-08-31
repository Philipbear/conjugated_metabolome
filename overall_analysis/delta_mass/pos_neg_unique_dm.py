import pickle
import numpy as np
import matplotlib.pyplot as plt


def main():
    with open('data/pos_delta_mass.pkl', 'rb') as f:
        pos_all = pickle.load(f)
    with open('data/neg_delta_mass.pkl', 'rb') as f:
        neg_all = pickle.load(f)

    print('All data loaded')

    # round to 2 decimals
    pos_all = [round(x, 2) for x in pos_all]
    neg_all = [round(x, 2) for x in neg_all]
    
    # unique rounded value
    unique_pos = set(pos_all)
    unique_neg = set(neg_all)
    print(f'all unique pos: {len(unique_pos)}')
    print(f'all unique neg: {len(unique_neg)}')
    
    shared_count = len(unique_pos.intersection(unique_neg))
    print(f'shared: {shared_count}')
    print(f'exclusively unique pos: {len(unique_pos) - shared_count}')
    print(f'exclusively unique neg: {len(unique_neg) - shared_count}')

if __name__ == '__main__':
    import os
    os.chdir(os.path.dirname(__file__))
    
    main()
    '''
    All data loaded
    all unique pos: 101526
    all unique neg: 41172
    shared: 38901
    exclusively unique pos: 62625
    exclusively unique neg: 2271
    '''