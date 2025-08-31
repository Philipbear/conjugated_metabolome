import pickle
import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_indicators
import numpy as np


def load_data():
    """
    Load unique annotation data for mouse body parts
    """
    print('Loading unique annotation data...')
    # load data
    with open('overall_analysis/redu/data/all_unique_annotation_mouse_bodypart.pkl', 'rb') as f:
        all = pickle.load(f)

    return all  # DataFrame with columns ['annotation', 'UBERONBodyPartName', 'count']


def prepare_upset_data(df):
    """
    Convert annotation-bodypart data to format suitable for UpSet plot
    
    Args:
        df: DataFrame with columns ['annotation', 'UBERONBodyPartName', 'count']
    
    Returns:
        pandas Series suitable for UpSet plotting
    """
    print('Preparing data for UpSet plot...')
    
    # Create a pivot table where annotations are rows and body parts are columns
    # Fill with 1 if annotation exists in that body part, 0 otherwise
    pivot_df = df.pivot_table(
        index='annotation', 
        columns='UBERONBodyPartName', 
        values='count', 
        fill_value=0
    )
    
    # Convert counts to binary (1 if annotation exists, 0 if not)
    binary_df = (pivot_df > 0).astype(bool)
    
    # Convert to format expected by upsetplot
    upset_data = from_indicators(binary_df)
    
    print(f"Found {len(binary_df)} unique annotations across {len(binary_df.columns)} body parts")
    print(f"Body parts: {list(binary_df.columns)}")
    
    return upset_data, binary_df


def create_upset_plot(upset_data, binary_df, figsize=(8.8, 3.5), min_subset_size=5):
    """
    Create and save UpSet plot
    """
    print('Creating UpSet plot...')
    
    # Set Arial font and figure size
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['figure.figsize'] = figsize
    plt.rcParams['font.size'] = 4.25
    
    # Create UpSet plot with visible elements and custom styling
    upset = UpSet(
        upset_data,
        subset_size='count',
        intersection_plot_elements=6,  # Show intersection bars
        totals_plot_elements=2,        # Show totals bars
        min_subset_size=min_subset_size,
        show_counts=True,
        sort_by='cardinality',         # Sort by number of sets in intersection
        sort_categories_by='-cardinality',  # Sort categories by total count
        element_size=6.5,               # Make dots smaller (default is usually around 32)
        facecolor='0.5',          # Optional: customize dot color
        connecting_line_width=1
    )
    
    # Generate the plot - this returns a dict of axes, not a figure
    axes_dict = upset.plot()
    
    # Get the figure from one of the axes
    fig = list(axes_dict.values())[0].figure
    
    # Now set the figure size
    fig.set_size_inches(figsize)
    
    # Adjust tick parameters and spine line width for bar plots
    for ax_name, ax in axes_dict.items():
        if ax_name == 'intersections':            
            # # set log scale for y-axis
            # ax.set_yscale('log')
            
            # # Set y-axis limits for intersection bars
            # ax.set_ylim(1, 13000)  # Adjust as needed
            
            # Intersection bar plot (top)
            ax.tick_params(
                axis='both',
                which='major',
                length=1,      # Tick length (shorter)
                width=0.5,     # Tick width (thinner)
                labelsize=3.5,   # Label size
                color='0.3',   # Tick color
                pad=1          # Distance between tick and label
            )
            # Remove minor ticks if present
            ax.tick_params(axis='both', which='minor', length=0, width=0)
            
            # grid
            ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.25, axis='y')
            
            # Adjust spine line width for intersection bars
            for spine in ax.spines.values():
                spine.set_linewidth(0.5)  # Make spine lines thinner
                spine.set_color('1')    # Optional: set spine color
            
            for patch in ax.patches:
                patch.set_facecolor('0.15')      # Bar fill color
                patch.set_edgecolor('0.15')  # Bar edge color
                
        elif ax_name == 'totals':
            # Totals bar plot (right side)
            ax.tick_params(
                axis='both',
                which='major',
                length=1,
                width=0.5,
                labelsize=3.5,
                color='0.3',
                pad=1
            )
            ax.tick_params(axis='both', which='minor', length=0, width=0)
            
            # grid
            ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.25, axis='x')
            
            # Adjust spine line width for totals bars
            for spine in ax.spines.values():
                spine.set_linewidth(0.5)  # Make spine lines thinner
                spine.set_color('0.4')    # Optional: set spine color
            for patch in ax.patches:
                patch.set_facecolor('0.15')      # Bar fill color
                patch.set_edgecolor('0.15')  # Bar edge color
    
    # Adjust position and format for intersection bar labels
    if 'intersections' in axes_dict:
        ax = axes_dict['intersections']
        for text in ax.texts:
            # Get current position
            x, y = text.get_position()
            
            # Adjust vertical position (move up/down)
            y_offset = 5500  # Adjust this value to move text up (+) or down (-)
            text.set_position((x, y + y_offset))
            
            # Adjust horizontal alignment
            text.set_ha('center')  # 'left', 'center', 'right'
            text.set_va('center')  # 'bottom', 'center', 'top'
            
            # Set font size
            text.set_fontsize(1.8)
            
            # Add comma formatting for thousands
            current_text = text.get_text()
            if current_text.isdigit():
                formatted_text = f"{int(current_text):,}"
                text.set_text(formatted_text)
    
    # Adjust position and format for totals bar labels
    if 'totals' in axes_dict:
        ax = axes_dict['totals']
        for text in ax.texts:
            # Get current position
            x, y = text.get_position()
            
            # Adjust horizontal position (move left/right)
            x_offset = 4800  # Adjust this value to move text left (-) or right (+)
            text.set_position((x + x_offset, y))
            
            # Adjust alignment
            text.set_ha('right')    # 'left', 'center', 'right'
            text.set_va('center')  # 'bottom', 'center', 'top'
            
            # Set font size
            text.set_fontsize(3.5)
            
            # Add comma formatting for thousands
            current_text = text.get_text()
            if current_text.isdigit():
                formatted_text = f"{int(current_text):,}"
                text.set_text(formatted_text)
        
    # Save the plot
    plt.tight_layout()
    plt.savefig('overall_analysis/redu/plots/upset_mouse_bodypart_annotations.svg', 
                format='svg', bbox_inches='tight', transparent=True)
    
    print('UpSet plot saved to overall_analysis/redu/plots/')
    
    return fig



def analyze_intersections(binary_df, top_n=20):
    """
    Analyze and print the top intersections
    
    Args:
        binary_df: Binary DataFrame
        top_n: Number of top intersections to analyze
    """
    print(f'\nAnalyzing top {top_n} intersections...')
    
    # Count annotations per body part
    body_part_counts = binary_df.sum(axis=0).sort_values(ascending=False)
    print(f'\nAnnotations per body part:')
    for bp, count in body_part_counts.head(10).items():
        print(f'  {bp}: {count:,} annotations')
    
    # Find intersections
    from upsetplot import from_indicators
    upset_data = from_indicators(binary_df)
    
    # Get top intersections
    top_intersections = upset_data.sort_values(ascending=False).head(top_n)
    
    print(f'\nTop {top_n} intersections:')
    
    # Get body part names from the binary_df columns
    body_part_names = list(binary_df.columns)
    
    for intersection, count in top_intersections.items():
        # intersection is a tuple of boolean values corresponding to each body part
        body_parts = [body_part_names[i] for i, present in enumerate(intersection) if present]
        
        if len(body_parts) == 1:
            print(f'  {body_parts[0]} only: {count:,} unique annotations')
        else:
            body_parts_str = ' ∩ '.join(body_parts)
            print(f'  {body_parts_str}: {count:,} unique annotations')
    
    return top_intersections


def create_summary_statistics(df, binary_df):
    """
    Create summary statistics for the dataset
    """
    print('\n=== Summary Statistics ===')
    
    total_annotations = len(binary_df)
    total_body_parts = len(binary_df.columns)
    total_spectra = df['count'].sum()
    
    print(f'Total unique annotations: {total_annotations:,}')
    print(f'Total body parts: {total_body_parts}')
    print(f'Total annotated spectra: {total_spectra:,}')
    
    # Annotations found in multiple body parts
    annotations_per_bodypart = binary_df.sum(axis=1)
    multi_bodypart_annotations = (annotations_per_bodypart > 1).sum()
    single_bodypart_annotations = (annotations_per_bodypart == 1).sum()
    
    print(f'Annotations found in multiple body parts: {multi_bodypart_annotations:,} ({multi_bodypart_annotations/total_annotations*100:.1f}%)')
    print(f'Annotations found in single body part: {single_bodypart_annotations:,} ({single_bodypart_annotations/total_annotations*100:.1f}%)')
    
    # Most ubiquitous annotations
    max_bodyparts = annotations_per_bodypart.max()
    most_ubiquitous = annotations_per_bodypart[annotations_per_bodypart == max_bodyparts].index
    print(f'Most ubiquitous annotation(s) found in {max_bodyparts} body parts: {len(most_ubiquitous)} annotation(s)')


def plot_body_part_distribution(binary_df, figsize=(10, 6)):
    """
    Create a bar plot showing annotation distribution per body part
    """
    print('Creating body part distribution plot...')
    
    # Count annotations per body part
    body_part_counts = binary_df.sum(axis=0).sort_values(ascending=False)
    
    plt.figure(figsize=figsize)
    
    # Create bar plot
    bars = plt.bar(range(len(body_part_counts)), body_part_counts.values, 
                   color='#c7522a', alpha=0.8)
    
    plt.xlabel('Body Parts', fontsize=12)
    plt.ylabel('Number of Unique Annotations', fontsize=12)
    plt.title('Distribution of Unique Annotations Across Mouse Body Parts', fontsize=14)

    # Set x-axis labels
    plt.xticks(range(len(body_part_counts)), body_part_counts.index, 
               rotation=45, ha='right', fontsize=10)
    
    # Add value labels on bars
    for i, bar in enumerate(bars):
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height):,}',
                ha='center', va='bottom', fontsize=8)
    
    plt.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    
    # Save plot
    plt.savefig('overall_analysis/redu/plots/bodypart_annotation_distribution.svg', 
                format='svg', bbox_inches='tight', transparent=True)
    
    print('Body part distribution plot saved to overall_analysis/redu/plots/')
    plt.show()


def main():
    """
    Main function to run the complete analysis
    """
    import os
    
    # Create output directory
    os.makedirs('overall_analysis/redu/plots', exist_ok=True)
    
    # Load data
    df = load_data()
    
    # Prepare data for UpSet plot
    upset_data, binary_df = prepare_upset_data(df)
    
    # Create summary statistics
    create_summary_statistics(df, binary_df)
    
    # Analyze intersections
    top_intersections = analyze_intersections(binary_df)
    
    # Create UpSet plot
    fig = create_upset_plot(upset_data, binary_df, figsize=(8.6, 3.5), min_subset_size=101)
    
    # # Create body part distribution plot
    # plot_body_part_distribution(binary_df)
    
    print('\nAnalysis complete! All plots saved to overall_analysis/redu/plots/')
    
    return df, binary_df, upset_data, top_intersections


if __name__ == '__main__':
    # Run the complete analysis
    df, binary_df, upset_data, top_intersections = main()
    