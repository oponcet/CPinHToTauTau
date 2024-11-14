import pickle
import hist
import numpy as np
import matplotlib.pyplot as plt
import mplhep
import scinum as sn



def load_histogram_data(pickle_file_path):
    """
    Load histogram data from a specified pickle file.

    Parameters:
    pickle_file_path (str): The path to the pickle file.

    Returns:
    dict or object: The loaded histogram data from the pickle file.

    Raises:
    FileNotFoundError: If the specified file is not found.
    Exception: For other errors during the loading process.
    """
    try:
        with open(pickle_file_path, 'rb') as f:
            hist_data = pickle.load(f)
            print("Histogram data loaded successfully!")
            return hist_data
    except FileNotFoundError:
        print(f"Error: The file '{pickle_file_path}' was not found.")
    except Exception as e:
        print(f"An unexpected error occurred while loading the file: {e}")
        raise



def hist_to_num(h: hist.hist, unc_name=str(sn.DEFAULT)) -> sn.Number:
    # Return an `sn.Number` object, where:
    # - The nominal values (i.e., bin contents) are `h.values()`.
    # - The uncertainties are stored in a dictionary where the key is `unc_name` and the value 
    #   is the square root of the variances (h.variances()**0.5).
    return sn.Number(h.values(), {unc_name: h.variances()**0.5})
    return sn.Number(h.values(), {unc_name: h.variances()**0.5})

# helper to integrate values stored in an array based number object
# Define a helper function `integrate_num` that integrates the values of an `sn.Number` object `num`.
# 

mplhep.style.use("CMS")

# Load histogram data from the pickle file
cf_store_path = "/eos/user/o/oponcet2/analysis/FakeFactor/analysis_httcp"

variable = "tau_1_pt"
category_num = "FFDRIso_tautau" # Numerator category
category_den = "FFDRantiIso_tautau" # Denominator category

print(f"Variable Name: {variable}")
print(f"Category Num: {category_num}")
print(f"Category Den: {category_den}")


pickle_file_path_num = cf_store_path + '/cf.DataDrivenEstimation/run3_2022_preEE_nano_cp_tau_v12/nominal/calib__main/sel__main_FF/prod__main_FF/datasets_19_130ad08fc8/htcondor2/qcd_histogram__'+category_num+'_'+variable+'.pickle'
print(f"Loading histogram data from: {pickle_file_path_num}")

pickle_file_path_den = cf_store_path + '/cf.DataDrivenEstimation/run3_2022_preEE_nano_cp_tau_v12/nominal/calib__main/sel__main_FF/prod__main_FF/datasets_19_130ad08fc8/htcondor2/qcd_histogram__'+category_den+'_'+variable+'.pickle'
print(f"Loading histogram data from: {pickle_file_path_den}")


ss_iso_qcd_hist = load_histogram_data(pickle_file_path_num)
ss_noniso_qcd_hist = load_histogram_data(pickle_file_path_den)

# print(f"ss_iso_qcd_hist : {ss_iso_qcd_hist}")
# print(f"ss_noniso_qcd_hist : {ss_noniso_qcd_hist}")


## Calculate ratio of histograms = Fake factors
ss_noniso_qcd = hist_to_num(ss_noniso_qcd_hist) # ss_noniso
ss_iso_qcd = hist_to_num(ss_iso_qcd_hist) # ss_iso

# print(f"ss_noniso_qcd : {ss_noniso_qcd}")
# print(f"ss_iso_qcd : {ss_iso_qcd}")

# calculate the pt-dependent fake factor
fake_factor = (ss_iso_qcd / ss_noniso_qcd)[:, None]
print(f"Fake Factor : {fake_factor}")
fake_factor_values = np.squeeze(np.nan_to_num(fake_factor()), axis=0)
fake_factor_variances = fake_factor(sn.UP, sn.ALL, unc=True)**2

# from IPython import embed; embed()

# get pt edges : 
pt_bin_edges = ss_iso_qcd_hist.axes['tau_1_pt'].edges
pt_centers = (pt_bin_edges[:-1] + pt_bin_edges[1:]) / 2  # Get the bin centers for plotting

# # Plot using Matplotlib
plt.figure(figsize=(8, 6))
plt.errorbar(pt_centers, fake_factor_values[0], yerr=fake_factor_variances[0], fmt='o', label='Fake Factor', color='black')
plt.xlabel('pT [GeV]')  # Adjust x-axis label as necessary
plt.ylabel('Fake Factor')
plt.title(f'Fake Factor vs pT for group "tautau"')
# plt.grid(True)
plt.legend()
plt.savefig(f"script_FF/plots/FakeFactors_script_tautau.png")  # save


### PLOT 
## Define custom axes and create a figure with subplots
fig = plt.figure()
grid = fig.add_gridspec(2, 1, hspace=0, height_ratios=[3, 1])
main_ax = fig.add_subplot(grid[0])
subplot_ax = fig.add_subplot(grid[1], sharex=main_ax)
plt.setp(main_ax.get_xticklabels(), visible=False)

# Define axis dictionary
ax_dict = {
    'main_ax': main_ax,
    'ratio_ax': subplot_ax  # This should match the argument in the plot_ratio function
}

# Plot ratio with custom axes
main_ax_artists, subplot_ax_artists = ss_iso_qcd_hist.project(variable).plot_ratio(
    ss_noniso_qcd_hist.project(variable),
    rp_ylabel=r"Iso/AntiIso",
    rp_num_label="Iso",
    rp_denom_label="AntiIso",
    ax_dict=ax_dict
)

#Set axis titles
main_ax.set_xlabel('Tau $p_{T}$ / GeV')  # X-axis title for the main axis
main_ax.set_ylabel("Events")  # Y-axis title for the main axis

# Set axis ranges
main_ax.set_xlim(0, 200)  # Set x-axis range for the main axis

subplot_ax.set_xlabel('Tau $p_{T}$ / GeV')  # X-axis title for the subplot axis
subplot_ax.set_ylabel(r"Iso/AntiIso")  # Y-axis title for the subplot axis

# Adjust subplot axis ranges if needed
subplot_ax.set_xlim(0, 200)  # Set x-axis range for the subplot axis
subplot_ax.set_ylim(-0.1, 1)  # Set y-axis range for the subplot axis (example range for pull)


plt.savefig(f"script_FF/plots/plot_{category_num}_{category_den}_{variable}_v2.png", format='png')