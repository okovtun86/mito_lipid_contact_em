# -*- coding: utf-8 -*-
"""
Created on Sun Jun  1 13:40:13 2025

@author: Oleg Kovtun, PhD
"""

#%%

import os
import numpy as np
import tifffile
from scipy.ndimage import binary_erosion, binary_dilation, generate_binary_structure
import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm

pixel_size_nm = 8
dilate_radii_nm = [24, 56]

mito_dir = "WT_Fasted_9428_mitochondria_instance_segmentation"
ld_dir = "lipid droplet segmentation slices/lipid droplet segmentation slices"

# output_csv = "mito_ld_contact_analysis.csv"
n_jobs = -1

footprint = generate_binary_structure(2, 1)
dilate_radii_px = [round(r / pixel_size_nm) for r in dilate_radii_nm]

all_mito_slices = sorted([f for f in os.listdir(mito_dir) if f.endswith(".tiff")])
all_ld_slices = sorted([f for f in os.listdir(ld_dir) if f.endswith(".tiff")])

if len(all_mito_slices) != len(all_ld_slices):
    raise ValueError("Mismatch in number of mitochondrial and lipid droplet slices.")

def process_slice(z, mito_file, ld_file):
    mito_path = os.path.join(mito_dir, mito_file)
    ld_path = os.path.join(ld_dir, ld_file)

    mito_slice = tifffile.imread(mito_path)  # 32-bit labeled mask
    ld_slice = tifffile.imread(ld_path) > 0   # 8-bit binary

    mito_labels = np.unique(mito_slice)
    mito_labels = mito_labels[mito_labels != 0]  

    results = []

    for label in mito_labels:
        mito_mask = mito_slice == label
        eroded = binary_erosion(mito_mask, structure=footprint)
        surface = mito_mask & ~eroded
        surface_area = surface.sum()

        result_row = {
            "z": z + 1,
            "mito_id": int(label),
            "surface_area": surface_area,
        }

        for r_nm, r_px in zip(dilate_radii_nm, dilate_radii_px):
            ld_expanded = ld_slice.copy()
            for _ in range(r_px):
                ld_expanded = binary_dilation(ld_expanded, structure=footprint)

            contact = surface & ld_expanded
            contact_area = contact.sum()
            pct_contact = (contact_area / surface_area) * 100 if surface_area > 0 else 0

            result_row[f"contact_area_nm{r_nm}"] = contact_area
            result_row[f"percent_contact_nm{r_nm}"] = pct_contact

        results.append(result_row)

    return results

print("\nStarting parallel per-slice analysis...\n")

parallel_results = Parallel(n_jobs=n_jobs)(
    delayed(process_slice)(z, mito_file, ld_file)
    for z, (mito_file, ld_file) in enumerate(zip(all_mito_slices, all_ld_slices))
)

results = [row for group in parallel_results for row in group]

df = pd.DataFrame(results)
# df.to_csv(output_csv, index=False)

# print(f"\nSaved results to: {output_csv}")
