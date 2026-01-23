import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.patches import Rectangle


def create_human_data():
    # Load the data
    df = pd.read_csv('data/redu_with_data.tsv', sep='\t', low_memory=False)
    df = df[df['all'] > 0].reset_index(drop=True)

    # Calculate annotation ratio
    df['annotation_ratio'] = df['annotated'] / df['all'] * 100  # Convert to percentage

    # Filter for human samples only
    df = df[df['NCBITaxonomy'].str.contains('9606|Homo sapiens', na=False)]

    # Filter out missing values from UBERONBodyPartName
    df = df[df['UBERONBodyPartName'] != 'missing value']
    df = df.dropna(subset=['UBERONBodyPartName'])
    
    # save
    df.to_csv('data/human_body_part_annotation_ratios.tsv', sep='\t', index=False)


def create_human_body_part_boxplot():

    # Load the data
    df = pd.read_csv('data/human_body_part_annotation_ratios.tsv', sep='\t', low_memory=False)

    # Get median values for each body part to rank them
    body_part_medians = df.groupby('UBERONBodyPartName')['annotation_ratio'].median().sort_values(ascending=False)

    # Keep only body parts with at least 3 samples
    counts = df['UBERONBodyPartName'].value_counts()
    valid_body_parts = counts[counts >= 3].index.tolist()

    # Filter dataframe to include only valid body parts
    df = df[df['UBERONBodyPartName'].isin(valid_body_parts)]

    # Re-calculate medians for the filtered data
    body_part_medians = df.groupby('UBERONBodyPartName')['annotation_ratio'].median().sort_values(ascending=False)

    # Create ordered categorical data type for proper sorting
    ordered_body_parts = body_part_medians.index.tolist()
    df['UBERONBodyPartName'] = pd.Categorical(df['UBERONBodyPartName'],
                                              categories=ordered_body_parts,
                                              ordered=True)

    # Set Arial font for all text elements
    plt.rcParams['font.family'] = 'Arial'

    # Create figure with appropriate size for vertical orientation
    fig, ax = plt.subplots(figsize=(2.9, 3.))  # Taller figure for vertical layout

    # Create a dictionary to store boxplot statistics for each body part
    boxplot_stats = {}

    # Calculate boxplot statistics manually for each body part
    for body_part in ordered_body_parts:
        data = df[df['UBERONBodyPartName'] == body_part]['annotation_ratio']
        q1 = data.quantile(0.25)
        median = data.median()
        q3 = data.quantile(0.75)
        iqr = q3 - q1
        
        whisker_bottom = data.quantile(0.05)
        whisker_top = data.quantile(0.95)

        boxplot_stats[body_part] = {
            'q1': q1,
            'median': median,
            'q3': q3,
            'whisker_bottom': whisker_bottom,
            'whisker_top': whisker_top
        }

    # First plot individual data points with jitter (BEHIND boxes)
    for i, body_part in enumerate(ordered_body_parts):
        # Get data for this body part
        values = df[df['UBERONBodyPartName'] == body_part]['annotation_ratio']

        # Add small random jitter to y position (vertical jitter now)
        jitter = np.random.uniform(-0.1, 0.1, size=len(values))

        # Plot points with lower zorder to ensure they're behind boxes
        # Note: x and y are swapped for vertical orientation
        ax.scatter(values, [i + j for j in jitter],
                   s=0.0012, color='0.8', alpha=1, zorder=1)

    # Then plot boxplots with higher zorder to ensure they're in front
    for i, body_part in enumerate(ordered_body_parts):
        stats = boxplot_stats[body_part]

        # Draw box (horizontal orientation now)
        box_height = 0.3
        box_bottom = i - box_height / 2
        box = Rectangle((stats['q1'], box_bottom),
                        stats['q3'] - stats['q1'], box_height,
                        facecolor='#F1DED7', edgecolor='0.4',
                        linewidth=0.5, zorder=2)
        ax.add_patch(box)

        # Draw median line (vertical line now)
        ax.plot([stats['median'], stats['median']],
                [box_bottom, box_bottom + box_height],
                color='#c7522a', linewidth=0.5, zorder=3)

        # Draw whiskers (horizontal now)
        ax.plot([stats['q1'], stats['whisker_bottom']], [i, i],
                color='0.4', linewidth=0.5, zorder=2)
        ax.plot([stats['q3'], stats['whisker_top']], [i, i],
                color='0.4', linewidth=0.5, zorder=2)

        # Draw caps on whiskers (vertical caps now)
        cap_height = 0.05
        ax.plot([stats['whisker_bottom'], stats['whisker_bottom']],
                [i - cap_height, i + cap_height],
                color='0.4', linewidth=0.5, zorder=2)
        ax.plot([stats['whisker_top'], stats['whisker_top']],
                [i - cap_height, i + cap_height],
                color='0.4', linewidth=0.5, zorder=2)

    # Style the plot
    ax.set_xlabel('Annotated metabolite conjugates (%)', fontsize=5, labelpad=3, color='0.2')
    ax.set_ylabel('')

    # Set y-ticks and labels (body parts on y-axis now)
    ax.set_yticks(range(len(ordered_body_parts)))
    ax.set_yticklabels(ordered_body_parts)

    # Invert y-axis so highest median is at top
    ax.invert_yaxis()

    ax.tick_params(axis='y', which='major', labelsize=4.5, colors='0.2',
                   length=1, pad=0.5, width=0.5)

    # Set x-axis limits
    ax.set_xlim(-2, 67)

    # Style the x-axis and frame
    _color = 0.4
    ax.tick_params(axis='x', which='major', length=1, width=0.5, pad=1,
                   colors=str(_color), labelsize=4)

    for spine in ax.spines.values():
        spine.set_color((_color, _color, _color))
        spine.set_linewidth(0.5)

    # Add grid lines (vertical now)
    ax.grid(axis='x', linestyle='--', alpha=0.7, color='0.9', linewidth=0.5)

    # Save the plot
    plt.tight_layout()
    plt.savefig('data/human_body_part_annotation_boxplot_vertical.svg', format='svg',
                bbox_inches='tight', transparent=True)

    print(f"Number of body parts plotted: {len(ordered_body_parts)}")


if __name__ == "__main__":
    import os
    os.chdir(os.path.dirname(__file__))
    
    create_human_data()
    # create_human_body_part_boxplot()