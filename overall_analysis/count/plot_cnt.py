"""
Count annotated spectra and unique annotations using parallel processing
"""

import os
import matplotlib.pyplot as plt


def plot_annotation_counts():
    """Create a single stacked bar plot showing annotation types"""
    
    # Data
    ref_2_annotated = 3546343
    ref_2_delta_mass = 20681096

    # Set Arial font
    plt.rcParams['font.family'] = 'Arial'
    
    # Create figure
    fig, ax = plt.subplots(figsize=(0.5, 0.98))
        
    # Create single stacked bar
    ax.bar(0, ref_2_annotated, bottom=ref_2_delta_mass, color='#c7522a', width=0.6)
    ax.bar(0, ref_2_delta_mass, color='#c7522a', width=0.6, alpha=0.4)
    
    # Add value labels on the bars
    total = ref_2_annotated + ref_2_delta_mass

    # Remove all axes and ticks
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    # Set limits to center the bar
    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(0, total * 1.1)
    
    plt.tight_layout()
    
    # Save figure as SVG    
    plt.savefig('data/annotation_stacked_bar.svg', format='svg', bbox_inches='tight', transparent=True)
    
    return fig

    
if __name__ == '__main__':
    
    import os
    os.chdir(os.path.dirname(__file__))
    # Create the stacked bar plot
    plot_annotation_counts()