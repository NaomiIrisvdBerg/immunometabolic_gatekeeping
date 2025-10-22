##This Python script binarises exported PNGs of whole slides loaded in the Xenium Explorer 3 software.
#We binarise cell segmentation maps and transcription maps of whole slides as follows:
#Cell segmentation maps were set to a single colour (yellow) in Xenium Explorer prior to PNG export via: Cells --> K-Means Clustering (K=2) --> Cluster 1 == yellow & Cluster 2 == yellow
#Transcript density maps (per panel as described in the manuscript) were set to 10uM bins at a Density map scale threshold of 0 to 0.01, and visualised using Inferno coloring,
#so that every 10uM bin containing at least one transcript positive for the respective panel (and meeting the transcript QC criteria as defaulted in Xenium Explorer 3 software) is shown as pale yellow as well
##hence the below code is able to binarise both the cell segmentation maps exported to PNGs as well as the transcript density maps exported to PNGs:

#Note that the base_dir needs to be changed in accordance with the tissue slide analysed and binarised.
#All PNGs exported from the respective slide were stored in the same folder that is generated when downloading from the Xenium database (see links in Manuscript, Supporting Materials IIIB)

import sys
import os
from PIL import Image
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

print("Using interpreter:", sys.executable)  ###interpreter used by us: miniconda3/bin/python  (v 3.10)
print("script started")

##path to the folder with the all PNGs (incl. PNG of cell segmentation map (made yellow) and the transcription denisty maps)
base_dir = "/Users/vandenn/Downloads/Xenium/Xenium_Preview_Human_Non_diseased_Lung_With_Add_on_FFPE_outs (1)"

##list all PNG files in the folder
png_files = [f for f in os.listdir(base_dir) if f.lower().endswith(".png")]

##Loop through each PNG
for fname in png_files:
    full_path = os.path.join(base_dir, fname)
    print(f"Processing: {fname}")
    
    try:
        ##load RGB image
        img = Image.open(full_path).convert("RGB")
        rgb = np.array(img)
        
        ##extract RGB channels
        R, G, B = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]
        R_norm = R / 255.0
        G_norm = G / 255.0
        B_norm = B / 255.0
        
        ##brightness filter
        eps = 1e-6
        brightness = (R + G + B) / (3 * 255.0)
        b_thresh = np.quantile(brightness, 0.995)
        is_bright = brightness >= b_thresh
        
        ##yellowness score
        yellowness = ((R_norm + G_norm) - 2 * B_norm) / ((R_norm + G_norm) + 2 * B_norm + eps)
        
        ##threshold for binarisation (adjustable)
        threshold = 0
        binary_mask = (yellowness > threshold).astype(int)

        print("Signal pixels (value = 1):", binary_mask.sum())
        
        ##convert to DataFrame
        rows, cols = binary_mask.shape
        df = pd.DataFrame([
            {"x": x, "y": y, "value": binary_mask[y, x]}
            for y in range(rows)
            for x in range(cols)
        ])
        
        ##output filename
        gene_name = os.path.splitext(fname)[0]
        output_path = os.path.join(base_dir, f"{gene_name}_filtered.csv")
        df.to_csv(output_path, index=False)
        
        print(f"Saved to: {output_path}")
        
    except Exception as e:
        print(f"Error processing {fname}: {e}")

print("script finished")
