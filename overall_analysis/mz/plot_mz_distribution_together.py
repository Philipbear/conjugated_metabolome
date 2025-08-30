import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import random  # Added for downsampling


def load_data():
    with open('data/pos_all_prec_mz.pkl', 'rb') as f:
        pos_all = pickle.load(f)
    with open('data/neg_all_prec_mz.pkl', 'rb') as f:
        neg_all = pickle.load(f)
    with open('data/pos_annotated_mz.pkl', 'rb') as f:
        pos_annotated = pickle.load(f)
    with open('data/neg_annotated_mz.pkl', 'rb') as f:
        neg_annotated = pickle.load(f)

    print('All data loaded')

    return pos_all, neg_all, pos_annotated, neg_annotated


def plot_mz_density_combined(pos_all, neg_all, pos_annotated, neg_annotated, save_path=None, sample_size=10000):
    """
    Create combined density plot for positive and negative mode m/z values.

    Parameters:
    -----------
    pos_all : list
        List of all m/z values from positive ionization mode
    neg_all : list
        List of all m/z values from negative ionization mode
    pos_annotated : list
        List of annotated m/z values from positive ionization mode
    neg_annotated : list
        List of annotated m/z values from negative ionization mode
    save_path : str, optional
        Path to save the figure, if None the figure is displayed
    sample_size : int, optional
        Number of points to sample for KDE computation

    Returns:
    --------
    fig : matplotlib.figure.Figure
        The generated figure object
    """
    # Downsample data for faster KDE computation
    pos_all_sample = random.sample(pos_all, min(sample_size, len(pos_all)))
    pos_annotated_sample = random.sample(pos_annotated, min(sample_size, len(pos_annotated)))
    neg_all_sample = random.sample(neg_all, min(sample_size, len(neg_all)))
    neg_annotated_sample = random.sample(neg_annotated, min(sample_size, len(neg_annotated)))

    print(f"Downsampled from {len(pos_all)} to {len(pos_all_sample)} points for all positive MS/MS")
    print(f"Downsampled from {len(pos_annotated)} to {len(pos_annotated_sample)} points for annotated positive MS/MS")
    print(f"Downsampled from {len(neg_all)} to {len(neg_all_sample)} points for all negative MS/MS")
    print(f"Downsampled from {len(neg_annotated)} to {len(neg_annotated_sample)} points for annotated negative MS/MS")

    # Set Arial font for all text elements
    plt.rcParams['font.family'] = 'Arial'

    # Create figure with two subplots sharing the x-axis
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(2, 1.1), sharex=True, gridspec_kw={'hspace': 0})

    # Compute x grid for KDE
    x_grid = np.linspace(0, 1500, 1000)

    # Compute KDEs with downsampled data for positive mode
    kde_pos_all = gaussian_kde(pos_all_sample)
    kde_pos_annotated = gaussian_kde(pos_annotated_sample)

    # Compute KDEs with downsampled data for negative mode
    kde_neg_all = gaussian_kde(neg_all_sample)
    kde_neg_annotated = gaussian_kde(neg_annotated_sample)

    # Plot densities for positive mode on the first subplot
    ax1.plot(x_grid, kde_pos_all(x_grid), color='#809bce', alpha=1, label='All public MS/MS', linewidth=1)
    ax1.fill_between(x_grid, kde_pos_all(x_grid), alpha=0.3, color='#809bce')

    ax1.plot(x_grid, kde_pos_annotated(x_grid), color='#c7522a', alpha=1, label='Annotated MS/MS', linewidth=1)
    ax1.fill_between(x_grid, kde_pos_annotated(x_grid), alpha=0.1, color='#c7522a')

    # Add "Positive Mode" text inside the plot
    ax1.text(0.6, 0.75, 'Positive ion mode', transform=ax1.transAxes, fontsize=5, color='0.2')

    # Plot densities for negative mode on the second subplot
    ax2.plot(x_grid, kde_neg_all(x_grid), color='#809bce', alpha=1, label='All public MS/MS', linewidth=1)
    ax2.fill_between(x_grid, kde_neg_all(x_grid), alpha=0.3, color='#809bce')

    ax2.plot(x_grid, kde_neg_annotated(x_grid), color='#c7522a', alpha=1, label='Annotated MS/MS', linewidth=1)
    ax2.fill_between(x_grid, kde_neg_annotated(x_grid), alpha=0.1, color='#c7522a')

    # Add "Negative Mode" text inside the plot
    ax2.text(0.6, 0.75, 'Negative ion mode', transform=ax2.transAxes, fontsize=5, color='0.2')

    # Set common x-label for both subplots
    ax2.set_xlabel('MS/MS precursor $\mathit{m/z}$', fontsize=5, labelpad=2, color='0.2')

    # Set y-labels for each subplot
    ax1.set_ylabel('Density', fontsize=5, labelpad=2, color='0.2')
    ax2.set_ylabel('Density', fontsize=5, labelpad=2, color='0.2')

    # Set common x-limit
    ax1.set_xlim(0, 1400)  # ax2 will share the same limit because of sharex=True

    # Adjust y-axis to start from 0 for both subplots
    y_max_pos = max(np.max(kde_pos_all(x_grid)), np.max(kde_pos_annotated(x_grid)))
    y_max_neg = max(np.max(kde_neg_all(x_grid)), np.max(kde_neg_annotated(x_grid)))

    ax1.set_ylim(0, y_max_pos * 1.1)
    ax2.set_ylim(0, y_max_neg * 1.1)

    # Add legend to both subplots
    legend1 = ax1.legend(frameon=False, fontsize=4.5,
                        loc='lower right',  # Position at lower right
                        bbox_to_anchor=(0.98, 0.01),  # Move legend lower
                        handlelength=1.0,  # Shorter line length
                        handletextpad=0.5,  # Space between line and text
                        labelspacing=0.2)

    legend2 = ax2.legend(frameon=False, fontsize=4.5,
                        loc='lower right',  # Position at lower right
                        bbox_to_anchor=(0.98, 0.01),  # Move legend lower
                        handlelength=1.0,  # Shorter line length
                        handletextpad=0.5,  # Space between line and text
                        labelspacing=0.2)

    # Change the color of the legend text for both legends
    for legend in [legend1, legend2]:
        for text in legend.get_texts():
            text.set_color('0.2')  # Set text color to 20% black

    # Format ticks only for the bottom subplot
    _color = 0.4
    ax2.tick_params(axis='x', which='major', length=1, width=0.5, pad=1,
                    colors=str(_color), labelsize=4)

    # Set custom x-ticks at 0, 200, 400, etc.
    x_ticks = np.arange(0, 1501, 200)  # Generate ticks from 0 to 1500 in steps of 200
    ax2.set_xticks(x_ticks)

    # Remove x-axis completely from top plot
    ax1.spines['bottom'].set_visible(False)
    ax1.xaxis.set_visible(False)

    # Hide y-ticks for both subplots
    ax1.set_yticks([])
    ax2.set_yticks([])

    # Set spine colors and width for both subplots
    for ax in [ax1, ax2]:
        for spine in ax.spines.values():
            spine.set_color((_color, _color, _color))
            spine.set_linewidth(0.5)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, format='svg', bbox_inches='tight', transparent=True)

    plt.show()

    return fig


if __name__ == '__main__':
    import os
    os.chdir(os.path.dirname(__file__))
    
    pos_all, neg_all, pos_annotated, neg_annotated = load_data()
    plot_mz_density_combined(pos_all, neg_all, pos_annotated, neg_annotated,
                             save_path='data/combined_mz_density.svg', sample_size=500000)