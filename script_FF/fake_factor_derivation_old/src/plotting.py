import ROOT
from ROOT import TCanvas, THStack, TLegend, TColor

# Define colors for specific datasets
COLORS = {
    'tt': "#9e9ac8",  # Violet
    'dy': "#feb24c",  # Orange
    'diboson': "#a96b59",  # Brown for Diboson (WW)
    'wj': "#d73027",  # Red for W+jets
    'higgs': "#253494",  # Dark blue for Higgs
    'fake': "#a1d99b"  # Green for Fake
}

LABELS = {
    'tt': "t/#bar{t}",
    'dy': "Drell-Yan",
    'diboson': "Diboson/Triboson",
    'wj': "W+jets",
    'higgs': "Higgs",
    'fake': "Fake"
}

def assign_color_and_label(dataset):
    """
    Assigns the appropriate color and label for the given dataset.
    """
    if dataset.startswith("tt") or dataset.startswith("st"):
        # print(f"Found ttbar dataset: {dataset}")
        return COLORS['tt'], LABELS['tt']
    elif dataset in ["ww", "wz", "zz", "www", "wwz", "wzz", "zzz"]:
        # print(f"Found diboson dataset: {dataset}")
        return COLORS['diboson'], LABELS['diboson']
    elif dataset.startswith("dy"):
        return COLORS['dy'], LABELS['dy']
    elif dataset.startswith("wj"):
        return COLORS['wj'], LABELS['wj']
    elif dataset.startswith("higgs"):
        return COLORS['higgs'], LABELS['higgs']
    elif dataset.startswith("fake"):
        return COLORS['fake'], LABELS['fake']
    else:
        return "#000000", "Unknown"