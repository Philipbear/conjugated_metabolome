import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def create_ranked_annotation_plot(fig_size=(3.45, 1.15)):

    # Load the data
    df = pd.read_csv('data/redu_with_data.tsv', sep='\t', low_memory=False)
    df = df[df['all'] > 0].reset_index(drop=True)

    # Calculate annotation ratio
    df['annotation_ratio'] = df['annotated'] / df['all'] * 100  # Convert to percentage

    # Filter out the specified categories
    excluded_categories = ['missing value', 'Synthetic and Chimeric']
    df = df[~df['NCBIDivision'].isin(excluded_categories)]

    # Merge vertebrates and invertebrates into animals
    df['NCBIDivision'] = df['NCBIDivision'].replace(['Vertebrates', 'Invertebrates'], 'Animals')
    df['NCBIDivision'] = df['NCBIDivision'].replace(['Primates', 'Rodents', 'Mammals'], 'Animals')

    # Get unique categories (NCBIDivision)
    categories = df['NCBIDivision'].unique()

    # Set Arial font for all text elements
    plt.rcParams['font.family'] = 'Arial'

    # Define colors and variables used in both plots
    # annotated_color = '#F1DED7'
    # non_annotated_color = '0.94'
    annotated_color = '#c7522a'
    non_annotated_color = '#809bce'

    # Process data once to get positions and statistics
    current_x = 0
    tick_positions = []
    tick_labels = []
    category_widths = []
    annotation_stats = []

    for category in categories:
        category_df = df[df['NCBIDivision'] == category].copy()
        category_df = category_df.sort_values('annotation_ratio', ascending=False)

        total_annotated = category_df['annotated'].sum()
        total_all = category_df['all'].sum()
        ratio = total_annotated / total_all * 100
        annotation_stats.append(f"Total: {ratio:.1f}%")

        middle_x = current_x + len(category_df) / 2
        tick_positions.append(middle_x)
        tick_labels.append(category)
        category_widths.append(len(category_df))

        current_x += len(category_df) + 2  # Add gap between categories

    # Get the final x limit
    x_limit = current_x - 2

    # STEP 1: Generate the data visualization only (PNG)
    # ==================================================
    fig_data, ax_data = plt.subplots(figsize=fig_size)

    current_x = 0
    for category in categories:
        category_df = df[df['NCBIDivision'] == category].copy()
        category_df = category_df.sort_values('annotation_ratio', ascending=False)

        for _, row in category_df.iterrows():
            annotated_height = row['annotation_ratio']
            non_annotated_height = 100 - annotated_height

            ax_data.bar(current_x, annotated_height, width=1, color=annotated_color, edgecolor='none', alpha=0.7)
            ax_data.bar(current_x, non_annotated_height, width=1, bottom=annotated_height,
                        color=non_annotated_color, edgecolor='none', alpha=0.2)

            current_x += 1

        current_x += 2  # Gap between categories

    # Remove all elements except the bars
    ax_data.set_xticks([])
    ax_data.set_yticks([])
    ax_data.set_xlim(-1, x_limit)
    ax_data.set_ylim(0, 100)

    for spine in ax_data.spines.values():
        spine.set_visible(False)

    # Save the data visualization as PNG
    plt.tight_layout()
    plt.savefig('data/annotation_ratio_data.png', format='png', dpi=600,
                bbox_inches='tight', transparent=True)
    plt.close(fig_data)

    # STEP 2: Create text and frame elements only (SVG)
    # =================================================
    fig_text, ax_text = plt.subplots(figsize=fig_size)

    # Set up the empty plot with correct dimensions
    ax_text.set_xlim(-1, x_limit)
    ax_text.set_ylim(0, 100)

    # Add y-axis label
    ax_text.set_ylabel('Annotated metabolite\nconjugates (%)', fontsize=5, labelpad=3, color='0.2')
    ax_text.set_xlabel('Public LC-MS/MS files (ordered)', fontsize=5, labelpad=2.5, color='0.2')

    # Remove default x-axis ticks and labels
    ax_text.set_xticks([])

    # Add category labels and stats
    for pos, label, stat, width in zip(tick_positions, tick_labels, annotation_stats, category_widths):
        ax_text.text(pos, 90, f'{label}', ha='center', va='top', fontfamily='Arial',
                    color='0.1', fontsize=5)
        ax_text.text(pos, 75, f'(n = {width:,})', ha='center', va='top', fontfamily='Arial',
                    color='0.1', fontsize=4.5)
        ax_text.text(pos, 60, stat, ha='center', va='top', fontfamily='Arial',
                    color='0.1', fontsize=4.5)

    # Style the y-axis and frame
    _color = 0.4
    ax_text.tick_params(axis='y', which='major', length=1, width=0.5, pad=1,
                        colors=str(_color), labelsize=4)

    for spine in ax_text.spines.values():
        spine.set_color((_color, _color, _color))
        spine.set_linewidth(0.5)

    # Add the legend
    import matplotlib.patches as mpatches
    annotated_patch = mpatches.Patch(color=annotated_color, label='Annotated metabolite conjugates', alpha=0.7)
    non_annotated_patch = mpatches.Patch(color=non_annotated_color, label='Other metabolites', alpha=0.2)

    legend = ax_text.legend(
        handles=[annotated_patch, non_annotated_patch],
        loc='upper center',
        bbox_to_anchor=(0.5, 1.27),
        ncol=2,
        frameon=False,
        prop={'family': 'Arial', 'size': 5},
    )

    for text in legend.get_texts():
        text.set_color('0.2')

    # Make the plot area transparent (remove any background fill)
    ax_text.patch.set_alpha(0)
    fig_text.patch.set_alpha(0)

    # Save the text and frame as SVG
    plt.tight_layout()
    plt.savefig('data/annotation_ratio_text.svg', format='svg',
                bbox_inches='tight', transparent=True)
    plt.close(fig_text)

    print(f"Total MRIs plotted: {len(df)}")
    print(f"Number of categories plotted: {len(categories)}")
    print("Files created:")
    print("  - data/annotation_ratio_data.png (raster data visualization)")
    print("  - data/annotation_ratio_text.svg (vector text and frame)")
    print("  - You can now combine these in PowerPoint slides")


if __name__ == "__main__":
    import os
    os.chdir(os.path.dirname(__file__))
    create_ranked_annotation_plot((4.05, 1.15))