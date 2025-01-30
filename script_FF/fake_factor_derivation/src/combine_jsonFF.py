'''
Author : @oponcet
Date : 29-01-2025
Script to combine the json files into one.
'''
import json
import glob

# Define the files and the structure
input_files = glob.glob("script_FF/fake_factor_derivation/outputs/FakeFactors/json/*.json")  # Replace with the path to your files
merged_structure = {
    "schema_version": 2,
    "description": "Fake factors for the httcp analysis",
    "corrections": [
        {
            "name": "fake_factors_fit",
            "description": "Fit fake factors for all decay modes and jet categories",
            "version": 1,
            "inputs": [
                {
                    "name": "pt",
                    "type": "real",
                    "description": "Transverse momentum of the tau"
                },
                {
                    "name": "dm",
                    "type": "int",
                    "description": "Reconstructed tau decay mode of leading tau: pi_1, rho_1, a1dm11_1, a1dm10_1, a1dm2_1"
                },
                {
                    "name": "njets",
                    "type": "int",
                    "description": "Number of jets in the event (has_0j, has_1j, has_2j)"
                },
                { 
                    "name": "syst",
                    "type": "string",
                    "description": "Systematic variations: 'nom', 'up', 'down'"
                }
            ],
            "output": {
                "name": "fake_factor",
                "type": "real",
                "description": "Fake factor to apply to data-MC"
            },
            "data": {
                "nodetype": "category",
                "input": "dm",
                "content": []
            }
        }
    ]
}

# Process each file and merge into the main structure
dm_content = {}

for file in input_files:
    with open(file, "r") as f:
        data = json.load(f)

    correction_data = data["corrections"][0]["data"]
    dm = correction_data["content"][0]["key"]  # Get the decay mode
    njets_content = correction_data["content"][0]["value"]["content"]  # Get jet-specific data

    # Add the decay mode and its associated jet categories to the main content
    if dm not in dm_content:
        dm_content[dm] = {"nodetype": "category", "input": "njets", "content": njets_content}
    else:
        # If the decay mode exists, extend its `njets` content
        dm_content[dm]["content"].extend(njets_content)

# Add all decay modes to the merged structure
merged_structure["corrections"][0]["data"]["content"] = [
    {"key": dm, "value": value} for dm, value in dm_content.items()
]

year = "2022_preEE"
# Save the merged JSON
output_file = f"script_FF/fake_factor_derivation/outputs/FakeFactors/json/fake_factor_{year}.json"
with open(output_file, "w") as f:
    json.dump(merged_structure, f, indent=4)

print(f"Merged JSON saved to {output_file}")
