# mito_lipid_contact_em

This script quantifies lipid droplet (LD) interactions with individual mitochondria (Mito) from per-slice Mito instance and semantic LD segmentations

For each 2D plane:

1. A 1-pixel-thick outer surface is computed for each Mito instance via binary erosion and subtraction.
2. The binary lipid droplet mask is dilated to simulate proximity zones at specified distances (e.g., 24 nm, 56 nm)
3. Overlap between each Mito surface and expanded LD region is computed to define the “interaction surface.”
4. The percentage of each Mito surface involved in contact is calculated as interaction_area / surface_area × 100

Test data: https://www.ebi.ac.uk/empiar/EMPIAR-12017/
