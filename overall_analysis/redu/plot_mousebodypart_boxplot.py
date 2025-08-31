def create_body_part_boxplot():
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np
    from matplotlib.patches import Rectangle

    # Load the data
    df = pd.read_csv('data/redu_with_data.tsv', sep='\t', low_memory=False)
    df = df[df['all'] > 0].reset_index(drop=True)

    # Calculate annotation ratio
    df['annotation_ratio'] = df['annotated'] / df['all'] * 100  # Convert to percentage

    # Filter for human samples only
    df = df[df['NCBITaxonomy'].isin(['10088|Mus', '10090|Mus musculus', '10105|Mus minutoides', '10114|Rattus', '10116|Rattus norvegicus'])]

    # Filter out missing values from UBERONBodyPartName
    df = df[df['UBERONBodyPartName'] != 'missing value']
    df = df.dropna(subset=['UBERONBodyPartName'])

    # Get median values for each body part to rank them
    body_part_medians = df.groupby('UBERONBodyPartName')['annotation_ratio'].median().sort_values(ascending=False)

    # Keep only body parts with at least 1 samples
    counts = df['UBERONBodyPartName'].value_counts()
    valid_body_parts = counts[counts >= 1].index.tolist()

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
    fig, ax = plt.subplots(figsize=(7, 1.5))

    # Create a dictionary to store boxplot statistics for each body part
    boxplot_stats = {}

    # Calculate boxplot statistics manually for each body part
    for body_part in ordered_body_parts:
        data = df[df['UBERONBodyPartName'] == body_part]['annotation_ratio']
        q1 = data.quantile(0.25)
        median = data.median()
        q3 = data.quantile(0.75)
        iqr = q3 - q1
        
        # whisker_bottom = max(data.min(), q1 - 1.5 * iqr)
        # whisker_top = min(data.max(), q3 + 1.5 * iqr)
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

        # Add small random jitter to x position
        jitter = np.random.uniform(-0.1, 0.1, size=len(values))

        # Plot points with lower zorder to ensure they're behind boxes
        ax.scatter([i + j for j in jitter], values,
                   s=0.0012, color='0.8', alpha=1, zorder=1)

    # Then plot boxplots with higher zorder to ensure they're in front
    for i, body_part in enumerate(ordered_body_parts):
        stats = boxplot_stats[body_part]

        # Draw box
        box_width = 0.3
        box_left = i - box_width / 2
        box = Rectangle((box_left, stats['q1']),
                        box_width, stats['q3'] - stats['q1'],
                        facecolor='#F1DED7', edgecolor='0.4',
                        linewidth=0.5, zorder=2)
        ax.add_patch(box)

        # Draw median line
        ax.plot([box_left, box_left + box_width],
                [stats['median'], stats['median']],
                color='#c7522a', linewidth=0.5, zorder=3)

        # Draw whiskers
        ax.plot([i, i], [stats['q1'], stats['whisker_bottom']],
                color='0.4', linewidth=0.5, zorder=2)
        ax.plot([i, i], [stats['q3'], stats['whisker_top']],
                color='0.4', linewidth=0.5, zorder=2)

        # Draw caps on whiskers
        cap_width = 0.05
        ax.plot([i - cap_width, i + cap_width],
                [stats['whisker_bottom'], stats['whisker_bottom']],
                color='0.4', linewidth=0.5, zorder=2)
        ax.plot([i - cap_width, i + cap_width],
                [stats['whisker_top'], stats['whisker_top']],
                color='0.4', linewidth=0.5, zorder=2)

    # Style the plot
    ax.set_xlabel('')
    ax.set_ylabel('Annotated metabolite\nconjugates (%)', fontsize=5, labelpad=2, color='0.2')

    # Set x-ticks and labels
    ax.set_xticks(range(len(ordered_body_parts)))
    ax.set_xticklabels(ordered_body_parts)

    ax.tick_params(axis='x', which='major', labelsize=4.5, colors='0.2',
                   length=1, rotation=30, pad=0.5, width=0.5)

    # Update x-tick labels
    labels = ax.get_xticklabels()
    for label in labels:
        label.set_ha('right')
    ax.set_xticklabels(labels)

    # Set y-axis limits
    ax.set_ylim(0, 70)

    # Style the y-axis and frame
    _color = 0.4
    ax.tick_params(axis='y', which='major', length=1, width=0.5, pad=1,
                   colors=str(_color), labelsize=4)

    for spine in ax.spines.values():
        spine.set_color((_color, _color, _color))
        spine.set_linewidth(0.5)

    # Add grid lines
    ax.grid(axis='y', linestyle='--', alpha=0.7, color='0.9', linewidth=0.5)

    # Save the plot
    plt.tight_layout()
    # plt.savefig('data/mouse_body_part_annotation_boxplot.png', format='png', dpi=600,
    #             bbox_inches='tight')
    plt.savefig('data/mouse_body_part_annotation_boxplot.svg', format='svg',
                bbox_inches='tight', transparent=True)
    plt.show()

    print(f"Number of body parts plotted: {len(ordered_body_parts)}")


if __name__ == "__main__":
    import os
    os.chdir(os.path.dirname(__file__))
    create_body_part_boxplot()