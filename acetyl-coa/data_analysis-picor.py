from pathlib import Path

import pandas as pd
import picor as ic

molecule = "acetyl-coa"
spect = pd.read_csv("acetyl-coa/acetylcoa_data.csv", index_col=0)
spect.drop(columns=["compound", "replicate"], inplace=True)
#spect.columns = ["No Label"] + [f"{i}C13" for i in range(1,8)]

# Get abundances relative to heighest peak
m = spect.max().max()
spect = spect / m * 100
res = ic.calc_isotopologue_correction(
    spect,
    molecule_name=molecule,
    resolution_correction=True,
    resolution=60000,
    molecules_file="metabolites.csv",
)

res.to_csv("data_iso_corr.csv")
