from flask.cli import F
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import normalize
from sklearn.decomposition import PCA
from umap import UMAP
import os


def load_ms2db_metadata():
    df = pd.read_pickle('/Users/shipei/Documents/projects/conjugated_metabolome/db/ms2db/all/all_ms2db_metadata.pkl')
        
    # inchikey_to_classifier = df.set_index('inchikey_14')['npclassifier_class_results'].to_dict()
    
    # inchikey_to_classifier = df.set_index('inchikey_14')['npclassifier_superclass_results'].to_dict()
    
    inchikey_to_classifier = df.set_index('inchikey_14')['npclassifier_pathway_results'].to_dict()
    
    # inchikey_to_classifier = df.set_index('inchikey_14')['classyfire_superclass'].to_dict()
    
    # inchikey_to_classifier = df.set_index('inchikey_14')['classyfire_class'].to_dict()
    
    # inchikey_to_classifier = df.set_index('inchikey_14')['classyfire_subclass'].to_dict()
    
    return inchikey_to_classifier


def get_compound_classes(compound_ids, inchikey_to_classifier):
    """Get chemical classes for compounds"""
    classes = []
    unknown_count = 0
    
    for inchikey_14 in compound_ids:
        
        if inchikey_14 and inchikey_14 in inchikey_to_classifier:
            class_info = inchikey_to_classifier[inchikey_14]
            # Handle NaN, None, empty string, and 'Unclassified' cases
            if pd.isna(class_info) or class_info is None or class_info == '' or class_info == 'Unclassified':
                classes.append('Unclassified')
                unknown_count += 1
            else:
                # Convert to string and strip whitespace in case of any formatting issues
                class_str = str(class_info).strip()
                if class_str == '' or class_str.lower() == 'nan':
                    classes.append('Unclassified')
                    unknown_count += 1
                else:
                    # Split by ';' and take the first part
                    if ';' in class_str:
                        # class_str = 'Unclassified'
                        class_str = class_str.split(';')[0].strip()
                    
                    # Check again if after splitting it's empty
                    if class_str == '':
                        classes.append('Unclassified')
                        unknown_count += 1
                    else:
                        classes.append(class_str)
        else:
            classes.append('Unclassified')
            unknown_count += 1
    
    print(f"Found classes for {len(compound_ids) - unknown_count}/{len(compound_ids)} compounds")
    print(f"Unclassified compounds: {unknown_count}")

    return classes


def check_processed_data_exists():
    """Check if processed data files exist"""
    data_folder = 'data'
    required_files = [
        'umap_feature_matrix_clean.npy',
        'compound_ids_clean.txt',
        'umap_embedding.npy'
    ]
    
    for file in required_files:
        if not os.path.exists(os.path.join(data_folder, file)):
            return False
    return True

def load_raw_data():
    """Load the prepared data from prepare_data.py output"""
    data_folder = 'data'
    
    # Load feature matrix
    feature_matrix = np.load(os.path.join(data_folder, 'umap_feature_matrix.npy'))
    
    # Load compound IDs
    with open(os.path.join(data_folder, 'compound_ids.txt'), 'r') as f:
        compound_ids = [line.strip() for line in f]
    
    # Load delta mass array
    delta_mass_array = np.load(os.path.join(data_folder, 'delta_mass_array.npy'))
    
    print(f"Loaded raw data: {feature_matrix.shape[0]} compounds, {feature_matrix.shape[1]} features")
    
    return feature_matrix, compound_ids, delta_mass_array

def load_processed_data():
    """Load already processed data"""
    data_folder = 'data'
    
    # Load cleaned feature matrix
    feature_matrix_clean = np.load(os.path.join(data_folder, 'umap_feature_matrix_clean.npy'))
    
    # Load cleaned compound IDs
    with open(os.path.join(data_folder, 'compound_ids_clean.txt'), 'r') as f:
        compound_ids_clean = [line.strip() for line in f]
    
    # Load UMAP embedding
    embedding = np.load(os.path.join(data_folder, 'umap_embedding.npy'))
    
    print(f"Loaded processed data: {len(compound_ids_clean)} compounds, {feature_matrix_clean.shape[1]} features")
    print(f"UMAP embedding shape: {embedding.shape}")
    
    return feature_matrix_clean, compound_ids_clean, embedding

def save_processed_data(feature_matrix_clean, compound_ids_clean, embedding):
    """Save processed data for future use"""
    data_folder = 'data'
    os.makedirs(data_folder, exist_ok=True)
    
    # Save cleaned feature matrix
    np.save(os.path.join(data_folder, 'umap_feature_matrix_clean.npy'), feature_matrix_clean)
    
    # Save cleaned compound IDs
    with open(os.path.join(data_folder, 'compound_ids_clean.txt'), 'w') as f:
        for compound_id in compound_ids_clean:
            f.write(f"{compound_id}\n")
    
    # Save UMAP embedding
    np.save(os.path.join(data_folder, 'umap_embedding.npy'), embedding)
    
    print("Processed data saved for future use")

def clean_and_normalize_data(feature_matrix, min_count=3, min_compounds=20, min_features=10):
    """
    Clean and normalize the feature matrix:
    1. Remove counts that are less than min_count
    2. Remove delta masses that appear in fewer than min_compounds compounds
    3. Remove compounds that have fewer than min_features features remaining
    4. Normalize each compound's vector to sum to 1
    """
    print(f"Original matrix shape: {feature_matrix.shape}")
    print(f"Original sparsity: {1 - np.count_nonzero(feature_matrix) / feature_matrix.size:.2%}")

    # Step 1: Remove counts that are less than min_count
    feature_matrix_count_filtered = np.where(feature_matrix >= min_count, feature_matrix, 0)
    
    # Count how many compounds have at least one feature above min_count
    valid_compounds_after_count = np.any(feature_matrix_count_filtered > 0, axis=1)
    print(f"After count filtering: keeping {np.sum(valid_compounds_after_count)} compounds with at least one feature >= {min_count}")
    
    # Apply the count filter
    feature_matrix_count_filtered = feature_matrix_count_filtered[valid_compounds_after_count, :]
    
    # Step 2: Remove features (delta masses) that appear in fewer than min_compounds compounds
    feature_counts = np.count_nonzero(feature_matrix_count_filtered, axis=0)
    valid_features = feature_counts >= min_compounds
    
    print(f"Removing {np.sum(~valid_features)} features that appear in < {min_compounds} compounds")
    feature_matrix_cleaned = feature_matrix_count_filtered[:, valid_features]
    
    # Step 3: Remove compounds that have fewer than min_features features remaining
    compound_feature_counts = np.count_nonzero(feature_matrix_cleaned, axis=1)
    valid_compounds_after_features = compound_feature_counts >= min_features

    print(f"Removing {np.sum(~valid_compounds_after_features)} compounds with fewer than {min_features} features")
    feature_matrix_cleaned = feature_matrix_cleaned[valid_compounds_after_features, :]
    
    # Create the final valid_compounds mask for the original matrix
    valid_compounds = np.zeros(len(feature_matrix), dtype=bool)
    
    # First apply count filter
    temp_indices = np.where(valid_compounds_after_count)[0]
    # Then apply feature filter to the remaining compounds
    final_indices = temp_indices[valid_compounds_after_features]
    valid_compounds[final_indices] = True
    
    print(f"Final valid compounds: {np.sum(valid_compounds)}")
    
    # Step 4: Normalization
    feature_matrix_cleaned = np.log1p(feature_matrix_cleaned)  # log(1+x) to handle zeros
    feature_matrix_normalized = normalize(feature_matrix_cleaned, norm='l1', axis=1)
    
    print(f"Final matrix shape: {feature_matrix_normalized.shape}")
    print(f"Final sparsity: {1 - np.count_nonzero(feature_matrix_normalized) / feature_matrix_normalized.size:.2%}")
    
    return feature_matrix_normalized, valid_compounds, valid_features

def reduce_dimensionality_pca(feature_matrix, n_components=50):
    """Apply PCA to reduce dimensionality before UMAP (optional preprocessing)"""
    print(f"Applying PCA to reduce from {feature_matrix.shape[1]} to {n_components} dimensions")
    
    pca = PCA(n_components=n_components, random_state=42)
    pca_features = pca.fit_transform(feature_matrix)
    
    explained_variance_ratio = np.sum(pca.explained_variance_ratio_)
    print(f"PCA explained variance ratio: {explained_variance_ratio:.3f}")
    
    return pca_features, pca

def create_umap_embedding(feature_matrix, n_neighbors=15, min_dist=0.1, n_components=2, metric='euclidean'):
    """Create UMAP embedding"""
    print(f"Creating UMAP embedding with n_neighbors={n_neighbors}, min_dist={min_dist}, metric={metric}")
    
    reducer = UMAP(  # Use UMAP directly instead of umap.UMAP
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        n_components=n_components,
        metric=metric,
        random_state=42,
        verbose=True
    )
    
    embedding = reducer.fit_transform(feature_matrix)
    
    return embedding, reducer

def plot_umap_basic_by_class(embedding, compound_ids, classes, highlight_classes=None, output_folder='plots'):
    """Create UMAP plots colored by chemical classes"""
    os.makedirs(output_folder, exist_ok=True)
    
    # Get unique classes and their counts
    unique_classes = list(set(classes))
    # Sort alphabetically but keep Unclassified last
    classified_classes = [cls for cls in unique_classes if cls != 'Unclassified']
    classified_classes.sort()  # Sort alphabetically
    unique_classes = classified_classes + ['Unclassified']
        
    # Custom color palette based on sector_colors
    sector_colors = {
        "Fatty acids": "#4e639e", 
        "Shikimates and Phenylpropanoids": "#e54616", 
        "Terpenoids": "#dba053",
        "Alkaloids": "#ff997c", 
        "Amino acids and Peptides": "#7fbfdd", 
        "Carbohydrates": "#96a46b",
        "Polyketides": "#760f00",
        "Unclassified": "#808080"  # grey for unclassified
    }
    
    # Fallback colors for classes not in sector_colors
    fallback_colors = [
        '#17becf',  # cyan
        '#bcbd22',  # olive
        '#9467bd',  # purple
        '#8c564b',  # brown
        '#e377c2',  # pink
        '#2ca02c',  # green
        '#d62728',  # red
    ]
    
    # Create color mapping
    class_to_color = {}
    fallback_index = 0
    
    for cls in unique_classes:
        if cls in sector_colors:
            class_to_color[cls] = sector_colors[cls]
        else:
            # Use fallback colors for classes not in the predefined list
            if fallback_index < len(fallback_colors):
                class_to_color[cls] = fallback_colors[fallback_index]
                fallback_index += 1
            else:
                class_to_color[cls] = '#808080'  # grey for additional classes
    
    # font for all text elements
    plt.rcParams['font.family'] = 'Helvetica'
    # Plot colored by class
    fig, ax = plt.subplots(figsize=(2.9, 1.6))

    for i, cls in enumerate(unique_classes):
        mask = np.array(classes) == cls
        if np.any(mask):
            ax.scatter(embedding[mask, 0], embedding[mask, 1], 
                       c=[class_to_color[cls]], label=f'{cls}',
                       alpha=1, s=0.5, edgecolor='none')

    ax.set_xlabel('UMAP 1', fontsize=4.5, labelpad=-6, color='0.2')
    ax.set_ylabel('UMAP 2', fontsize=4.5, labelpad=-7, color='0.2')

    # no ticks
    ax.set_xticks([])
    ax.set_yticks([])

    ax.legend(bbox_to_anchor=(0.95, 0.45), loc='center left', fontsize=5, markerscale=4.5, frameon=False)
    # plt.grid(True, alpha=0.3)    
    
    # # remove top and right spines
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    # for spine in ax.spines.values():
    #     spine.set_color('0.2')
    #     spine.set_linewidth(0.5)
    
    # Get data range for half-length spines
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    
    # Calculate half ranges
    x_range = x_max - x_min
    y_range = y_max - y_min
    
    # Set custom spine positions (half length)
    x_spine_start = x_min
    x_spine_end = x_min + x_range * 0.35
    
    # Y-axis spine (left) - half length, centered  
    y_spine_start = y_min
    y_spine_end = y_min + y_range * 0.35
    
    # Remove default spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    # Add custom half-length spines
    ax.plot([x_spine_start, x_spine_end], [y_min, y_min], 
            color='0.2', linewidth=0.5, clip_on=False)  # bottom spine
    ax.plot([x_min, x_min], [y_spine_start, y_spine_end], 
            color='0.2', linewidth=0.5, clip_on=False)  # left spine
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, 'umap_by_class.png'), dpi=600, bbox_inches='tight', transparent=True)
    # plt.show()


def plot_umap_more_by_class(embedding, compound_ids, classes, highlight_classes=None, output_folder='plots'):
    os.makedirs(output_folder, exist_ok=True)
    
    # Get unique classes and their counts
    unique_classes = list(set(classes))
    class_counts = {cls: classes.count(cls) for cls in unique_classes}
    # Create a simplified plot with top classes only
    top_classes = sorted(class_counts.items(), key=lambda x: x[1], reverse=True)[:10]
    top_class_names = [cls for cls, _ in top_classes]
    
    plt.figure(figsize=(15, 12))
    
    # Plot "Other" category first (background)
    other_mask = ~np.isin(classes, top_class_names)
    if np.any(other_mask):
        plt.scatter(embedding[other_mask, 0], embedding[other_mask, 1], 
                   c='lightgray', label=f'Other (n={np.sum(other_mask)})',
                   alpha=0.3, s=6, edgecolor='none')
    
    # Plot top classes
    colors_top = plt.cm.tab10(np.linspace(0, 1, len(top_class_names)))
    for i, cls in enumerate(top_class_names):
        mask = np.array(classes) == cls
        if np.any(mask):
            plt.scatter(embedding[mask, 0], embedding[mask, 1], 
                       c=[colors_top[i]], label=f'{cls} (n={class_counts[cls]})',
                       alpha=0.8, s=6, edgecolor='none')
    
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    plt.title(f'UMAP of Chemical Compounds by Top 10 Classes\n({len(compound_ids)} compounds)')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10, markerscale=3)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, 'umap_by_top_classes.png'), dpi=600, bbox_inches='tight')
    plt.show()
    
    # Create highlighted classes plot if highlight_classes is provided
    if highlight_classes:
        plt.figure(figsize=(15, 12))
        
        # Find which classes are actually present in the data
        present_highlight_classes = [cls for cls in highlight_classes if cls in unique_classes]
        missing_highlight_classes = [cls for cls in highlight_classes if cls not in unique_classes]
        
        if missing_highlight_classes:
            print(f"\nWarning: The following highlight classes were not found in the data: {missing_highlight_classes}")
        
        # Plot non-highlighted classes first (background)
        non_highlight_mask = ~np.isin(classes, present_highlight_classes)
        if np.any(non_highlight_mask):
            plt.scatter(embedding[non_highlight_mask, 0], embedding[non_highlight_mask, 1], 
                       c='lightgray', label=f'Other classes (n={np.sum(non_highlight_mask)})',
                       alpha=0.3, s=3, edgecolor='none')
        
        # Plot highlighted classes with distinct colors
        colors_highlight = plt.cm.tab10(np.linspace(0, 1, len(present_highlight_classes)))
        
        for i, cls in enumerate(present_highlight_classes):
            mask = np.array(classes) == cls
            if np.any(mask):
                plt.scatter(embedding[mask, 0], embedding[mask, 1], 
                           c=[colors_highlight[i]], label=f'{cls} (n={class_counts[cls]})',
                           alpha=0.8, s=12, edgecolor='black', linewidth=0.3)
        
        plt.xlabel('UMAP 1')
        plt.ylabel('UMAP 2')
        plt.title(f'UMAP of Chemical Compounds - Highlighted Classes\n({len(compound_ids)} compounds)')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10, markerscale=2)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(output_folder, 'umap_highlighted_classes.png'), dpi=600, bbox_inches='tight')
        plt.show()
        
        # Print summary of highlighted classes
        print(f"\nHighlighted classes summary:")
        for cls in present_highlight_classes:
            count = class_counts[cls]
            print(f"  {cls}: {count} compounds ({count/len(classes)*100:.1f}%)")


def plot_umap(embedding, compound_ids, output_folder='plots', highlight_inchikeys=None):
    """Create basic UMAP plots (without classes)"""
    os.makedirs(output_folder, exist_ok=True)
    
    # Basic scatter plot
    plt.figure(figsize=(12, 10))
    
    if highlight_inchikeys is not None and len(highlight_inchikeys) > 0:
        # Create masks for highlighted and non-highlighted compounds
        highlight_mask = np.isin(compound_ids, highlight_inchikeys)
        
        # Plot non-highlighted compounds first (background)
        non_highlight_mask = ~highlight_mask
        if np.any(non_highlight_mask):
            plt.scatter(embedding[non_highlight_mask, 0], embedding[non_highlight_mask, 1], 
                       alpha=0.3, s=6, edgecolor='none', c='lightgray', 
                       label=f'Other compounds (n={np.sum(non_highlight_mask)})')
        
        # Plot highlighted compounds on top
        if np.any(highlight_mask):
            plt.scatter(embedding[highlight_mask, 0], embedding[highlight_mask, 1], 
                       alpha=0.8, s=15, edgecolor='black', linewidth=0.5, c='red',
                       label=f'Highlighted compounds (n={np.sum(highlight_mask)})')
        
        plt.legend(loc='upper right', markerscale=2)
        print(f"Highlighted {np.sum(highlight_mask)} out of {len(highlight_inchikeys)} requested compounds")
        
        # Print which compounds were found/not found
        found_compounds = set(compound_ids) & set(highlight_inchikeys)
        missing_compounds = set(highlight_inchikeys) - set(compound_ids)
        
        if found_compounds:
            print(f"Found compounds: {len(found_compounds)}")
        if missing_compounds:
            print(f"Missing compounds: {len(missing_compounds)}")
            for compound in list(missing_compounds)[:5]:  # Show first 5 missing
                print(f"  - {compound}")
            if len(missing_compounds) > 5:
                print(f"  ... and {len(missing_compounds) - 5} more")
    else:
        # Regular plot without highlighting
        plt.scatter(embedding[:, 0], embedding[:, 1], alpha=0.6, s=3, edgecolor='none')
    
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    title = f'UMAP of Chemical Compounds\n({len(compound_ids)} compounds)'
    if highlight_inchikeys is not None and len(highlight_inchikeys) > 0:
        title += f' - {len(highlight_inchikeys)} highlighted'
    plt.title(title)
    
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Add suffix to filename if highlighting
    filename = 'umap_basic.png'
    if highlight_inchikeys is not None and len(highlight_inchikeys) > 0:
        filename = 'umap_basic_highlighted.png'
    
    plt.savefig(os.path.join(output_folder, filename), dpi=600, bbox_inches='tight')
    plt.show()


def process_data_pipeline(min_count=3, min_compounds=20, min_features=10, pca_n_components=50):
    """Complete data processing pipeline"""
    print("Processing raw data...")
    
    # Load raw data
    feature_matrix, compound_ids, delta_mass_array = load_raw_data()
    
    # Clean and normalize data
    feature_matrix_clean, valid_compounds, valid_features = clean_and_normalize_data(
        feature_matrix, min_count, min_compounds, min_features
    )
    
    # Update compound IDs to match cleaned data
    compound_ids_clean = [compound_ids[i] for i in range(len(compound_ids)) if valid_compounds[i]]
    delta_mass_clean = delta_mass_array[valid_features]
    
    print(f"\nFinal data: {len(compound_ids_clean)} compounds, {len(delta_mass_clean)} features")
    
    # Optional: Apply PCA for dimensionality reduction (useful for very high-dimensional data)
    pca_features, pca = reduce_dimensionality_pca(feature_matrix_clean, pca_n_components)
    input_features = pca_features
    print("Using PCA-reduced features for UMAP")
    
    # Create UMAP embedding
    embedding, reducer = create_umap_embedding(
        input_features, 
        n_neighbors=15, 
        min_dist=0.1, 
        metric='euclidean'
    )
    
    # Save processed data
    save_processed_data(feature_matrix_clean, compound_ids_clean, embedding)
    
    return feature_matrix_clean, compound_ids_clean, embedding

def main(reprocess_data=False, highlight_inchikeys=None, highlight_classes=None, min_count=3,
         min_compounds=20, min_features=10, pca_n_components=50):
    """Main analysis pipeline"""
    
    # Load metadata for chemical classes
    print("Loading chemical class metadata...")
    inchikey_to_classifier = load_ms2db_metadata()
    
    # Check if processed data already exists
    if reprocess_data:
        feature_matrix_clean, compound_ids_clean, embedding = process_data_pipeline(min_count, min_compounds, min_features, pca_n_components)
    elif check_processed_data_exists():
        print("Processed data found. Loading existing data...")
        feature_matrix_clean, compound_ids_clean, embedding = load_processed_data()
    else:
        print("No processed data found. Processing raw data...")
        feature_matrix_clean, compound_ids_clean, embedding = process_data_pipeline(min_count, min_compounds, min_features, pca_n_components)

    # Get chemical classes for compounds
    print("\nMapping compounds to chemical classes...")
    classes = get_compound_classes(compound_ids_clean, inchikey_to_classifier)
    
    # Create plots
    print("\nGenerating plots...")
    # plot_umap(embedding, compound_ids_clean, highlight_inchikeys=highlight_inchikeys)  # Basic plots with highlighting
    plot_umap_basic_by_class(embedding, compound_ids_clean, classes, highlight_classes=highlight_classes)  # Class-colored plots
    # plot_umap_more_by_class(embedding, compound_ids_clean, classes, highlight_classes=highlight_classes)  # Class-colored plots

    # Save final results with classes
    results_df = pd.DataFrame({
        'compound_id': compound_ids_clean,
        'umap_1': embedding[:, 0],
        'umap_2': embedding[:, 1],
        'chemical_class': classes
    })
    
    results_df.to_csv('data/umap_results.csv', index=False)
    print("\nResults saved to 'data/umap_results.csv'")


if __name__ == "__main__":
    import os
    os.chdir(os.path.dirname(__file__))

    main(reprocess_data=False, min_count=5, min_compounds=200, min_features=20, pca_n_components=50)
    