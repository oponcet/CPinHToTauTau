import uproot
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep

# Use CMS style for mplhep
plt.style.use(hep.style.CMS)

# Define colors and labels
COLORS = {
    'tt': "#9e9ac8",  # Violet
    'dy': "#feb24c",  # Orange
    'diboson': "#a96b59",  # Brown for Diboson
    'wj': "#d73027",  # Red for W+jets
    'higgs': "#253494",  # Dark blue for Higgs
    'fake': "#a1d99b"   # Green for Fake
}

LABELS = {
    'tt': "t/#bar{t}",
    'dy': "Drell-Yan",
    'diboson': "Diboson/Triboson",
    'wj': "W+jets",
    'higgs': "Higgs",
    'fake': "Fake"
}

def assign_color_and_label(dataset_name):
    """
    Assign color and label to histograms based on dataset name.
    """
    for key in COLORS.keys():
        if key in dataset_name.lower():
            return COLORS[key], LABELS[key]
    return "#000000", "Unknown"  # Default color and label if no match

def create_stack_plot_hep(mc_histograms, data_histogram, var, cat_dir):
    """
    Create a Matplotlib stack plot for the MC histograms and overlay the data histogram.

    Parameters:
    - mc_histograms (dict): Dictionary of MC histograms as {name: (counts, bin_edges)}.
    - data_histogram (tuple): Data histogram as (counts, bin_edges).
    - var (str): Variable name to label the plot.
    - cat_dir (str): Directory where the plot will be saved.
    """
    # Prepare data for stacking
    stack_counts = []
    stack_labels = []
    stack_colors = []
    all_bin_edges = None

    for hist_name, (counts, bin_edges) in mc_histograms.items():
        if all_bin_edges is None:
            all_bin_edges = bin_edges  # Use the first histogram's bin edges
        color, label = assign_color_and_label(hist_name)
        stack_counts.append(counts)
        stack_labels.append(label)
        stack_colors.append(color)

    # Create the plot
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot stacked MC histograms
    hep.histplot(
        np.array(stack_counts),
        all_bin_edges,
        label=stack_labels,
        stack=True,
        color=stack_colors,
        histtype="fill",
        edgecolor="black"
    )

    # Overlay data histogram, if available
    if data_histogram:
        data_counts, data_bin_edges = data_histogram
        bin_centers = 0.5 * (data_bin_edges[:-1] + data_bin_edges[1:])
        ax.errorbar(
            bin_centers,
            data_counts,
            yerr=np.sqrt(data_counts),  # Poisson uncertainty
            fmt='o',
            color="black",
            label="Data"
        )

    # Set axis labels and title
    ax.set_xlabel(var)
    ax.set_ylabel("Events")

    # Add a legend
    ax.legend(loc="upper right", fontsize=10)

    # Save the plot
    output_path = f"{cat_dir}/{var}_mc_stack_with_data.png"
    plt.savefig(output_path, dpi=300)
    print(f"Plot saved to: {output_path}")

    # Show the plot (optional)
    plt.close(fig)


# Example of usage
# Assuming `mc_histograms` is a dictionary with MC histograms, and `data_histogram` is provided
# Structure: {histogram_name: (counts, bin_edges)}, where counts = np.array, bin_edges = np.array

# Replace the following paths with actual file paths and names
root_file_path = "path/to/your_file.root"
thstack_name = "path/to/your_THStack"
data_histogram_name = "path/to/your_Data_Histogram"
output_dir = "output/plots"

# Extract histograms from the ROOT file
with uproot.open(root_file_path) as root_file:
    # Extract MC histograms
    mc_stack = root_file[thstack_name]
    mc_histograms = {
        hist.name: (hist.values(), hist.axis().edges()) for hist in mc_stack.values()
    }

    # Extract Data histogram (if available)
    data_histogram = None
    if data_histogram_name in root_file:
        data_hist = root_file[data_histogram_name]
        data_histogram = (data_hist.values(), data_hist.axis().edges())

# Variable name for labeling
variable_name = "MET"  # Example variable

# Create the plot
create_stack_plot_and_summary(mc_histograms, data_histogram, variable_name, output_dir)
