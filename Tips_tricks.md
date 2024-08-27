# CF_inspect Cheat Sheet

## In Selection

### Basic Usage

```bash
cf_inspect /eos/user/o/oponcet/code/analysis/CP/cf_store/analysis_httcp/cf.SelectEvents/run3_2022_preEE_nano_cp_tau_v12_limited/data_tau_C/nominal/calib__main/sel__main_FF/DRiso_DRantiIso-4/results_0.parquet

# Access the loaded data
data = objects[0]

# Check the structure
print(ak.fields(data))

# Show the first few entries
print(ak.to_list(data[:5]))

# Check the number of events
len(data)

# Display specific fields/columns
print(ak.to_list(data["objects"]["Tau"][:10]))

# Display nested fields
print(ak.to_list(data["objects"]["Tau"]["Tau"][:10]))

# Summarize the data
total_sum = ak.sum(data["objects"]["Tau"]["Tau"][:])
print(total_sum)  # Replace "field_name" with the actual name of the field

import matplotlib.pyplot as plt
import numpy as np
import awkward as ak

# Flatten and convert data for visualization
Tau = ak.flatten(data["objects"]["Tau"]["Tau"])
Tau = ak.to_numpy(Tau_flat)

# Create a histogram
plt.figure(figsize=(10, 6))
plt.hist(Tau_numpy, bins=50, alpha=0.7, color='blue', edgecolor='black')
plt.title("Distribution of Tau")
plt.xlabel("Value")
plt.ylabel("Number of Events")
plt.grid(True)
plt.show()
