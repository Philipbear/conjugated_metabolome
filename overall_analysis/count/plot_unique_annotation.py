"""
Create scatter plots for annotation frequency analysis
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
import os

def load_annotation_data():
    """Load the pickle files containing annotation data"""
    
    # Load data
    with open('overall_analysis/count/data/all_unique_annotation.pkl', 'rb') as f:
        data = pickle.load(f)
    
    return data

def create_scatter_plot(annotations_dict, x_label, output_path, color='#c7522a'):
    """
    Create a scatter plot for annotation frequencies
    
    Args:
        annotations_dict: Dictionary of annotations and their counts
        title: Title for the plot
        output_path: Path to save the plot
        color: Color for the scatter points
    """
    # Convert to sorted list (high to low frequency)
    sorted_annotations = sorted(annotations_dict.items(), key=lambda x: x[1], reverse=True)
    
    # Extract frequencies
    frequencies = [count for _, count in sorted_annotations]
    
    # Create x-axis positions (annotation rank)
    x_positions = range(1, len(frequencies) + 1)
    
    # Set Arial font
    plt.rcParams['font.family'] = 'Arial'
    
    # Create figure
    fig, ax = plt.subplots(figsize=(1.65, 1.15))
    
    # Create scatter plot
    ax.scatter(x_positions, frequencies, 
              color=color, 
              alpha=1, 
              s=1,  # Point size
              edgecolors='none')
    
    # Set log scale for y-axis to handle wide range of frequencies
    ax.set_yscale('log')
    ax.set_xscale('log')
    
    # Set labels and title
    ax.set_xlabel(x_label, fontsize=5, color='0.2', labelpad=1.25)
    ax.set_ylabel('Frequency', fontsize=5, color='0.2', labelpad=1.2)
    # ax.set_title(title, fontsize=5, color='0.2', pad=5)
    
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
    ax.set_xlim(0, len(frequencies) * 1.02)
    ax.set_ylim(0, max(frequencies) * 1.2)
    
    plt.tight_layout()
    
    # Save the figure
    plt.savefig(output_path + '.png', dpi=1000, bbox_inches='tight', transparent=True)
    # plt.savefig(output_path + '.svg', format='svg', bbox_inches='tight', transparent=True)
    plt.close()
    
    print(f"Saved plot: {output_path}")
    
    # print out how many annotations have frequency >= 3
    high_freq_count = sum(1 for freq in frequencies if freq >= 3)
    print(f"Number of annotations with frequency >= 3: {high_freq_count}")
    # print out how many annotations have frequency >= 5
    high_freq_count_5 = sum(1 for freq in frequencies if freq >= 5)
    print(f"Number of annotations with frequency >= 5: {high_freq_count_5}")
    # print out how many annotations have frequency >= 10
    high_freq_count_10 = sum(1 for freq in frequencies if freq >= 10)
    print(f"Number of annotations with frequency >= 10: {high_freq_count_10}")


def plot_all_annotation_frequencies():
    """Create all scatter plots for annotation frequencies"""
    
    # Load data
    data = load_annotation_data()
    
    # Create output directory
    output_dir = 'overall_analysis/count/plots'
    os.makedirs(output_dir, exist_ok=True)
    
    # Define colors for different plot types
    spec_spec_color = '0.7'
    spec_delta_color = '0.7'
    
    print("Creating annotation frequency scatter plots...")
    
    # spec-spec annotations
    create_scatter_plot(
        data['spec_spec'],
        'Structure-structure pairs',
        f'{output_dir}/spec_spec_frequencies',
        color=spec_spec_color
    )
    
    # spec-delta annotations
    create_scatter_plot(
        data['spec_delta'],
        'Structure-delta mass pairs',
        f'{output_dir}/spec_delta_frequencies',
        color=spec_delta_color
    )
    
    print(f"\nAll plots saved to: {output_dir}")


if __name__ == '__main__':
        
    # Create individual plots
    plot_all_annotation_frequencies()        
    print("\nAll plots completed!")
    