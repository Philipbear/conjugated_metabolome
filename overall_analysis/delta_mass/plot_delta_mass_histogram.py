import pickle
import numpy as np
import matplotlib.pyplot as plt


def load_data():
    with open('data/pos_delta_mass.pkl', 'rb') as f:
        pos_all = pickle.load(f)
    with open('data/neg_delta_mass.pkl', 'rb') as f:
        neg_all = pickle.load(f)

    print('All data loaded')

    return pos_all, neg_all


def plot_delta_mass_histogram(delta_mass_list, mode='pos', save_name=None):
    """
    Create bar histogram for delta mass values.
    """
    # Round values to 2 decimal places
    rounded_values = [round(x, 2) for x in delta_mass_list]

    # Count frequency of each rounded value
    unique_values, counts = np.unique(rounded_values, return_counts=True)

    print(f"Processing {len(rounded_values)} delta mass values for {mode} mode")
    print(f"Found {len(unique_values)} unique rounded values")

    # Set Arial font for all text elements
    plt.rcParams['font.family'] = 'Arial'

    # Create figure
    fig, ax = plt.subplots(figsize=(2.85, 1.1))

    # Plot histogram
    color = '#809bce' if mode == 'pos' else '#c7522a'
    ax.bar(unique_values, counts, width=0.25, color=color, alpha=1,
           label=f"{'Positive' if mode == 'pos' else 'Negative'} mode")

    # Set labels and style
    ax.set_xlabel('Delta mass (Da)', fontsize=5, labelpad=2, color='0.2')
    ax.set_ylabel('Frequency', fontsize=5, labelpad=3, color='0.2')

    # Set x-axis limits
    x_min, x_max = 0, 180
    ax.set_xlim(x_min, x_max)

    # Set y-axis limit based on max frequency within the x-axis range
    mask = (unique_values >= x_min) & (unique_values <= x_max)
    if np.any(mask):
        max_freq_in_range = np.max(counts[mask])
        ax.set_ylim(0, max_freq_in_range * 1.1)  # Set y-limit to 1.1 times max frequency

        # print top freq and their corresponding delta mass
        top_indices = np.argsort(counts[mask])[-10:]
        top_values = unique_values[mask][top_indices]
        top_counts = counts[mask][top_indices]
        print(f"Top delta mass values in range {x_min} to {x_max}:")
        for value, count in zip(top_values, top_counts):
            print(f"Delta mass: {value}, Frequency: {count}")
            
            

    # Add legend
    legend = ax.legend(frameon=False, fontsize=5, loc='upper right', handlelength=1.0, handletextpad=0.5)
    for text in legend.get_texts():
        text.set_color('0.2')

    # Style adjustments
    _color = 0.4
    ax.tick_params(axis='both', which='major', length=1, width=0.5, pad=1,
                   colors=str(_color), labelsize=4)

    for spine in ax.spines.values():
        spine.set_color((_color, _color, _color))
        spine.set_linewidth(0.5)

    plt.tight_layout()

    # Save figure if path is provided
    if save_name:
        # plt.savefig(f'{save_name}.svg', format='svg', bbox_inches='tight', transparent=True)
        plt.savefig(f'{save_name}.png', format='png', dpi=1000, bbox_inches='tight', transparent=True)
        print(f"Figure saved to {save_name}")

    # plt.show()
    return fig


def plot_delta_mass_histogram_both_modes(pos_delta_mass_list, neg_delta_mass_list, save_name=None):
    """
    Create merged bar histogram for delta mass values with both ion modes.
    Uses single color and combines all frequencies.
    """
    # Round values to 2 decimal places for both modes
    pos_rounded = [round(x, 2) for x in pos_delta_mass_list]
    neg_rounded = [round(x, 2) for x in neg_delta_mass_list]
    
    # Combine both datasets
    all_rounded = pos_rounded + neg_rounded
    
    # Count frequency of each rounded value
    unique_values, counts = np.unique(all_rounded, return_counts=True)
    
    print(f"Processing {len(pos_rounded)} positive and {len(neg_rounded)} negative delta mass values")
    print(f"Combined total: {len(all_rounded)} values")
    print(f"Found {len(unique_values)} unique delta mass values")
    
    # Set Arial font for all text elements
    plt.rcParams['font.family'] = 'Arial'
    
    # Create figure
    fig, ax = plt.subplots(figsize=(3.05, 1.1))
    
    # Plot bars with single color
    color = '#c7522a'
    width = 0.35
    
    ax.bar(unique_values, counts, width=width, color=color, alpha=1, label='Both modes combined')
    
    # Set labels and style
    ax.set_xlabel('Delta mass (Da)', fontsize=5, labelpad=2, color='0.2')
    ax.set_ylabel('Frequency', fontsize=5, labelpad=3, color='0.2')
    
    # Set x-axis limits
    x_min, x_max = 0, 480
    ax.set_xlim(x_min, x_max)
    
    # Set y-axis limit based on max frequency within the x-axis range
    if len(unique_values) > 0:
        mask = (unique_values >= x_min) & (unique_values <= x_max)
        if np.any(mask):
            max_freq_in_range = np.max(counts[mask])
            ax.set_ylim(0, max_freq_in_range * 1.1)
            
            # Print top 2 values for each interval
            print(f"Top delta mass values in 20 Da intervals from {x_min} to {x_max}:")
            
            for interval_start in range(x_min, x_max, 20):
                interval_end = min(interval_start + 20, x_max)
                interval_mask = (unique_values >= interval_start) & (unique_values < interval_end)
                
                if np.any(interval_mask):
                    interval_values = unique_values[interval_mask]
                    interval_counts = counts[interval_mask]
                    
                    # Get top in this interval
                    top_indices = np.argsort(interval_counts)[-3:]
                    top_values = interval_values[top_indices]
                    top_counts = interval_counts[top_indices]
                    
                    print(f"  Interval [{interval_start}-{interval_end}) Da:")
                    for value, count in zip(reversed(top_values), reversed(top_counts)):
                        print(f"    Delta mass: {value}, Frequency: {count}")
    
    # # Add legend
    # legend = ax.legend(frameon=False, fontsize=4, loc='upper right', 
    #                   handlelength=1.0, handletextpad=0.5)
    # for text in legend.get_texts():
    #     text.set_color('0.2')
    
    # Style adjustments
    _color = 0.4
    ax.tick_params(axis='both', which='major', length=1, width=0.5, pad=1,
                   colors=str(_color), labelsize=4)
    
    for spine in ax.spines.values():
        spine.set_color((_color, _color, _color))
        spine.set_linewidth(0.5)
    
    plt.tight_layout()
    
    # Save figure if path is provided
    if save_name:
        plt.savefig(f'{save_name}.png', format='png', dpi=1000, bbox_inches='tight', transparent=True)
        print(f"Figure saved to {save_name}")
    
    return fig


if __name__ == '__main__':
    import os
    os.chdir(os.path.dirname(__file__))
    
    pos_delta_mass, neg_delta_mass = load_data()
    # plot_delta_mass_histogram(pos_delta_mass, mode='pos', save_name='data/pos_delta_mass_histogram')
    # plot_delta_mass_histogram(neg_delta_mass, mode='neg', save_name='data/neg_delta_mass_histogram')
    
    # Plot combined histogram for both modes
    plot_delta_mass_histogram_both_modes(pos_delta_mass, neg_delta_mass, save_name='data/both_modes_delta_mass_histogram')