"""
Create scatter plots for annotation frequency analysis based on unique datasets
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
from collections import defaultdict

def load_annotation_data():
    """Load the pickle files containing annotation data"""
    
    # Load merged data
    spec_spec = pd.read_pickle('overall_analysis/count/data/all_unique_spec_annotation_merged.pkl')
    spec_delta = pd.read_pickle('overall_analysis/count/data/all_unique_delta_annotation_merged.pkl')
    
    return spec_spec, spec_delta

def count_datasets_per_annotation(df):
    """
    Count unique datasets for each annotation
    
    Args:
        df: DataFrame with 'dataset' and 'annotations' columns
    
    Returns:
        dict: annotation -> number of unique datasets
    """
    annotation_datasets = defaultdict(set)
    
    for _, row in df.iterrows():
        annotation_datasets[row['annotations']].add(row['dataset'])
    
    # Convert sets to counts
    return {annotation: len(datasets) for annotation, datasets in annotation_datasets.items()}

def create_scatter_plot(annotation_counts, x_label, output_path, color='#c7522a'):
    """
    Create a scatter plot for annotation dataset frequencies
    
    Args:
        annotation_counts: Dictionary of annotations and their dataset counts
        x_label: Label for x-axis
        output_path: Path to save the plot
        color: Color for the scatter points
    """
    # Convert to sorted list (high to low frequency)
    sorted_annotations = sorted(annotation_counts.items(), key=lambda x: x[1], reverse=True)
    
    # Extract dataset counts
    dataset_counts = [count for _, count in sorted_annotations]
    
    # Create x-axis positions (annotation rank)
    x_positions = range(1, len(dataset_counts) + 1)
    
    # Set Arial font
    plt.rcParams['font.family'] = 'Arial'
    
    # Create figure
    fig, ax = plt.subplots(figsize=(1.65, 1.15))
    
    # Create scatter plot
    ax.scatter(x_positions, dataset_counts, 
              color=color, 
              alpha=1, 
              s=1,  # Point size
              edgecolors='none')
    
    # Set log scale for x-axis, linear for y-axis (since dataset counts are typically smaller)
    ax.set_xscale('log')
    # Use log scale for y-axis only if there's a wide range of values
    if max(dataset_counts) > 50:
        ax.set_yscale('log')
    
    # Set labels
    ax.set_xlabel(x_label, fontsize=5, color='0.2', labelpad=1.25)
    ax.set_ylabel('Number of unique datasets', fontsize=5, color='0.2', labelpad=1.2)
    
    # Style the plot
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.25)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('0.4')
    ax.spines['bottom'].set_color('0.4')
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)
    
    # Style ticks
    ax.tick_params(axis='both', which='major', labelsize=3.8, length=1.25, width=0.5, pad=1, color='0.4')
    ax.tick_params(axis='both', which='minor', length=0.5, width=0.25, color='0.4')
        
    # Set axis limits with some padding
    ax.set_xlim(0, len(dataset_counts) * 1.02)
    ax.set_ylim(0, max(dataset_counts) * 1.2)
    
    plt.tight_layout()
    
    # Save the figure
    plt.savefig(output_path + '.png', dpi=1000, bbox_inches='tight', transparent=True)
    plt.close()
    
    print(f"Saved plot: {output_path}")
    
    # Print statistics about dataset coverage
    dataset_counts_values = list(dataset_counts)
    high_coverage_3 = sum(1 for count in dataset_counts_values if count >= 3)
    high_coverage_5 = sum(1 for count in dataset_counts_values if count >= 5)
    high_coverage_10 = sum(1 for count in dataset_counts_values if count >= 10)
    
    print(f"Number of annotations appearing in >= 3 datasets: {high_coverage_3}")
    print(f"Number of annotations appearing in >= 5 datasets: {high_coverage_5}")
    print(f"Number of annotations appearing in >= 10 datasets: {high_coverage_10}")
    print(f"Max datasets for single annotation: {max(dataset_counts_values)}")
    print(f"Total unique annotations: {len(dataset_counts_values)}")


def plot_all_annotation_frequencies():
    """Create all scatter plots for annotation dataset frequencies"""
    
    # Load data
    spec_spec_df, spec_delta_df = load_annotation_data()
    
    # Count datasets per annotation
    spec_spec_counts = count_datasets_per_annotation(spec_spec_df)
    spec_delta_counts = count_datasets_per_annotation(spec_delta_df)
    
    # Create output directory
    output_dir = 'overall_analysis/count/plots'
    os.makedirs(output_dir, exist_ok=True)
    
    # Define colors for different plot types
    spec_spec_color = '0.7'
    spec_delta_color = '0.7'
    
    print("Creating annotation dataset frequency scatter plots...")
    
    # spec-spec annotations
    print("\nStructure-structure pairs:")
    create_scatter_plot(
        spec_spec_counts,
        'Structure-structure pairs',
        f'{output_dir}/spec_spec_dataset_frequencies',
        color=spec_spec_color
    )
    
    # spec-delta annotations
    print("\nStructure-delta mass pairs:")
    create_scatter_plot(
        spec_delta_counts,
        'Structure-delta mass pairs',
        f'{output_dir}/spec_delta_dataset_frequencies',
        color=spec_delta_color
    )
    
    print(f"\nAll plots saved to: {output_dir}")


if __name__ == '__main__':
        
    # Create individual plots
    plot_all_annotation_frequencies()
    