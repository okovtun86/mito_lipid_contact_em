# -*- coding: utf-8 -*-
"""
Created on Sun Jun  1 13:14:22 2025

@author: Oleg kovtun, PhD
"""

#%%

import os
import numpy as np
import tifffile
from scipy.ndimage import binary_erosion, binary_dilation, generate_binary_structure
import pandas as pd


pixel_size_nm = 8
dilate_radii_nm = [24, 56]
dilate_radii_px = [round(r / pixel_size_nm) for r in dilate_radii_nm]

mito_dir = "WT_Fasted_9428_mitochondria_instance_segmentation"
ld_dir = "lipid droplet segmentation slices/lipid droplet segmentation slices"

# structuring element (disk for 2D)
footprint = generate_binary_structure(2, 1)

results = []

for z in range(0, 9500):
    
    print('Processing z slice ', z)
    mito_path = os.path.join(mito_dir, f"{z:04d}.tiff")
    ld_path = os.path.join(ld_dir, f"{z:04d}.tiff")

    mito_slice = tifffile.imread(mito_path)  # 32-bit labeled instances
    ld_slice = tifffile.imread(ld_path) > 0   # 8-bit binary

    mito_labels = np.unique(mito_slice)
    mito_labels = mito_labels[mito_labels != 0]  

    for label in mito_labels:
        mito_mask = mito_slice == label
        eroded = binary_erosion(mito_mask, structure=footprint)
        surface = mito_mask & ~eroded

        surface_area = surface.sum()

        for r_nm, r_px in zip(dilate_radii_nm, dilate_radii_px):
            ld_expanded = ld_slice.copy()
            for _ in range(r_px):
                ld_expanded = binary_dilation(ld_expanded, structure=footprint)

            contact = surface & ld_expanded
            contact_area = contact.sum()
            pct_contact = (contact_area / surface_area) * 100 if surface_area > 0 else 0

            results.append({
                "z": z,
                "mito_id": int(label),
                "surface_area": surface_area,
                "contact_area_nm{}".format(r_nm): contact_area,
                "percent_contact_nm{}".format(r_nm): pct_contact
            })

df = pd.DataFrame(results)

# df.to_csv("mito_ld_contact_analysis.csv", index=False)
