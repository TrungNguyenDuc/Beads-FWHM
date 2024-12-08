# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 21:26:38 2024

@author: S233755
"""

import os
import numpy as np
import tifffile as tiff
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from skimage.measure import label, regionprops
from skimage.filters import threshold_otsu
from tkinter import Tk, filedialog

# Gaussian function for curve fitting
def gaussian(x, a, b, c):
    return a * np.exp(-((x - b) / c) ** 2)

# Input parameters
lateral_pixel_size = 147
axial_pixel_size = 147
pad = 1
roi_radius = 15
threshold_value = 1000
data_number = [3]
fit_quality = 0.85


def process_tiff_files(lateral_pixel_size=147, axial_pixel_size=147, pad=1, roi_radius=15, threshold_value=1000, fit_quality=0.85):
   
    # Select multiple files
    root = Tk()
    root.withdraw()
    file_paths = filedialog.askopenfilenames(filetypes=[("TIFF files", "*.tif")])
    if not file_paths:
        print("File selection canceled.")
        exit()
    
    print(f"Selected files: {file_paths}")
    
    # Initialize final results container
    final_results_all = []
    
    xRes_mean_all = []
    yRes_mean_all = []
    zRes_mean_all = []
    
    xRes_std_all = []
    yRes_std_all = []
    zRes_std_all = []
    
    for file_path in file_paths:
        print(f"\nProcessing file: {file_path}")
        # Load TIFF data
        im_data = tiff.imread(file_path).astype(np.float64)
        im_data = np.pad(im_data, ((pad, pad), (pad, pad), (pad, pad)), mode="constant")  # Added padding uniformly
        print("Data loaded.")
        
        # Determine the dimensions
        dims = np.array(im_data.shape)
        sorted_dims = np.argsort(dims)  # Sort dimensions to identify shortest, second, and longest
        z_idx, x_idx, y_idx = sorted_dims  # Assign indices for z, x, and y
        
        # Rearrange data for consistent axis ordering: (z, x, y)
        im_data = np.transpose(im_data, (z_idx, x_idx, y_idx))
        
        # Thresholding
        im_bw = im_data > threshold_value
        plt.figure()
        plt.imshow(np.max(im_bw, axis=0), cmap="gray")  # Project along z-axis for visualization
        plt.title("Thresholded Data")
        plt.show()
        
        # Particle identification
        stats = regionprops(label(im_bw))
        print("Particles identified.")
        
        final_results = {
            "file_name": os.path.basename(file_path),
            "fwhmX": [],
            "fwhmY": [],
            "fwhmZ": [],
            "xR2": [],
            "yR2": [],
            "zR2": [],
            "PeakIn": [],
            "xAxis": [],
            "yAxis": [],
            "xRes": [],
            "yRes": [],
            "zRes": [],
            "PEAKY": [],
        }
        
        # Process each particle
        for idx, particle in enumerate(stats):
            print(f"Processing particle {idx + 1}/{len(stats)}")
            
            centroid = np.round(particle.centroid).astype(int)
            start_idx = centroid - roi_radius
            end_idx = centroid + roi_radius + 1
        
            if np.any(start_idx < 0) or np.any(end_idx >= im_data.shape):
                continue
        
            # Extract region of interest (ROI)
            im_roi = im_data[start_idx[0]:end_idx[0], start_idx[1]:end_idx[1], start_idx[2]:end_idx[2]]
        
            # Find the brightest voxel
            max_idx = np.unravel_index(np.argmax(im_roi), im_roi.shape)
            line_z = im_roi[:, max_idx[1], max_idx[2]]
            line_x = im_roi[max_idx[0], :, max_idx[2]]
            line_y = im_roi[max_idx[0], max_idx[1], :]
        
            # Normalize lines
            line_x = (line_x - np.min(line_x)) / (np.max(line_x) - np.min(line_x))
            line_y = (line_y - np.min(line_y)) / (np.max(line_y) - np.min(line_y))
            line_z = (line_z - np.min(line_z)) / (np.max(line_z) - np.min(line_z))
        
            # Fit Gaussian and calculate FWHM
            try:
                x = np.arange(line_x.size)
                y = np.arange(line_y.size)
                z = np.arange(line_z.size)
                
                # Fit Gaussian to x data
                popt_x, _ = curve_fit(gaussian, x, line_x, bounds=([0, 0.25 * x.size, 0], [1, 0.75 * x.size, 20]))
                fitted_x = gaussian(x, *popt_x)
        
                # Fit Gaussian to y data
                popt_y, _ = curve_fit(gaussian, y, line_y, bounds=([0, 0.25 * y.size, 0], [1, 0.75 * y.size, 20]))
                fitted_y = gaussian(y, *popt_y)
        
                # Fit Gaussian to z data
                popt_z, _ = curve_fit(gaussian, z, line_z, bounds=([0, 0.25 * z.size, 0], [1, 0.75 * z.size, 20]))
                fitted_z = gaussian(z, *popt_z)
        
                # Calculate R^2 for x
                ss_res_x = np.sum((line_x - fitted_x) ** 2)  # Residual sum of squares
                ss_tot_x = np.sum((line_x - np.mean(line_x)) ** 2)  # Total sum of squares
                r2_x = 1 - (ss_res_x / ss_tot_x)
        
                # Calculate R^2 for y
                ss_res_y = np.sum((line_y - fitted_y) ** 2)  # Residual sum of squares
                ss_tot_y = np.sum((line_y - np.mean(line_y)) ** 2)  # Total sum of squares
                r2_y = 1 - (ss_res_y / ss_tot_y)
        
                # Calculate R^2 for z
                ss_res_z = np.sum((line_z - fitted_z) ** 2)  # Residual sum of squares
                ss_tot_z = np.sum((line_z - np.mean(line_z)) ** 2)  # Total sum of squares
                r2_z = 1 - (ss_res_z / ss_tot_z)
        
                scale_factor = 2 * np.sqrt(2 * np.log(2))
                fwhm_x = popt_x[2] * lateral_pixel_size * scale_factor / np.sqrt(2)
                fwhm_y = popt_y[2] * lateral_pixel_size * scale_factor / np.sqrt(2)
                fwhm_z = popt_z[2] * axial_pixel_size * scale_factor / np.sqrt(2)
            except RuntimeError:
                continue
        
            # Store results
            final_results["fwhmX"].append(fwhm_x)
            final_results["fwhmY"].append(fwhm_y)
            final_results["fwhmZ"].append(fwhm_z)
            final_results["xR2"].append(r2_x)
            final_results["yR2"].append(r2_y)
            final_results["zR2"].append(r2_z)
            final_results["PeakIn"].append(np.max(im_roi))
        
        # Plot resolution scatter plots
        x_axis = [stat.centroid[1] * lateral_pixel_size/1000 for stat in stats]
        y_axis = [stat.centroid[0] * lateral_pixel_size/1000 for stat in stats]
        
        plt.figure()
        plt.subplot(3, 1, 1)
        plt.scatter(x_axis[0:len(final_results["fwhmX"])], final_results["fwhmX"], color="k")
        plt.title("X-Resolution")
        plt.xlabel("Lateral Position (µm)")
        plt.ylabel("Resolution (nm)")
        
        plt.subplot(3, 1, 2)
        plt.scatter(x_axis[0:len(final_results["fwhmY"])], final_results["fwhmY"], color="k")
        plt.title("Y-Resolution")
        plt.xlabel("Lateral Position (µm)")
        plt.ylabel("Resolution (nm)")
        
        plt.subplot(3, 1, 3)
        plt.scatter(x_axis[0:len(final_results["fwhmZ"])], final_results["fwhmZ"], color="k")
        plt.title("Z-Resolution")
        plt.xlabel("Lateral Position (µm)")
        plt.ylabel("Resolution (nm)")
        
        plt.tight_layout()
        plt.show()
        
        
        # Initialize counters and arrays
        counter = 0
        xAxis = []
        yAxis = []
        xRes = []
        yRes = []
        zRes = []
        PEAKY = []
        
        # Iterate through each bead
        for beadIdx in range(len(final_results["zR2"])):
            # Quality of fit for each dimension
            qualityOfFit = [
                final_results["xR2"][beadIdx],
                final_results["yR2"][beadIdx],
                final_results["zR2"][beadIdx]
            ]
            Peak = final_results["PeakIn"][beadIdx]
            
            # If all fits are good (> 0.95) and other criteria are met, proceed
            if all(np.array(qualityOfFit) > fit_quality):  # You can add additional conditions here
                counter += 1
                xAxis.append(stats[beadIdx]["Centroid"][1]* lateral_pixel_size/1000)
                yAxis.append(stats[beadIdx]["Centroid"][0]* lateral_pixel_size/1000)
                xRes.append(final_results["fwhmX"][beadIdx])
                yRes.append(final_results["fwhmY"][beadIdx])
                zRes.append(final_results["fwhmZ"][beadIdx])
                PEAKY.append(Peak)
        
        # Convert results to NumPy arrays for further processing
        xAxis = np.array(xAxis)
        yAxis = np.array(yAxis)
        xRes = np.array(xRes)
        yRes = np.array(yRes)
        zRes = np.array(zRes)
        PEAKY = np.array(PEAKY)
        
        xRes_mean = np.mean(xRes)
        yRes_mean = np.mean(yRes)
        zRes_mean = np.mean(zRes)
        
        xRes_std = np.std(xRes)
        yRes_std = np.std(yRes)
        zRes_std = np.std(zRes)
        
        
        final_results["xAxis"] = xAxis
        final_results["yAxis"] = yAxis
        final_results["xRes"] = xRes
        final_results["yRes"] = yRes
        final_results["zRes"] = zRes
        final_results["PEAKY"] = PEAKY
        
        
        plt.figure()
        plt.subplot(3, 1, 1)
        plt.scatter(xAxis, xRes, color="k")
        plt.title("X-Resolution")
        plt.xlabel("Lateral Position (µm)")
        plt.ylabel("Resolution (nm)")
        
        plt.subplot(3, 1, 2)
        plt.scatter(xAxis, yRes, color="k")
        plt.title("Y-Resolution")
        plt.xlabel("Lateral Position (µm)")
        plt.ylabel("Resolution (nm)")
        
        plt.subplot(3, 1, 3)
        plt.scatter(xAxis, zRes, color="k")
        plt.title("Z-Resolution")
        plt.xlabel("Lateral Position (µm)")
        plt.ylabel("Resolution (nm)")
        
        plt.tight_layout()
        plt.show()
        
        
        # Debugging: print filtered results
        print(f"Filtered xAxis: {xAxis}")
        print(f"Filtered yAxis: {yAxis}")
        print(f"Filtered xRes: {xRes}")
        print(f"Filtered yRes: {yRes}")
        print(f"Filtered zRes: {zRes}")
        print(f"Filtered PEAKY: {PEAKY}")
        print("Resolution analysis complete.")
        # Append current file results
        final_results_all.append(final_results)
        
        xRes_mean_all.append(xRes_mean)
        yRes_mean_all.append(yRes_mean)
        zRes_mean_all.append(zRes_mean)
    
        xRes_std_all.append(xRes_std)
        yRes_std_all.append(yRes_std)
        zRes_std_all.append(zRes_std)
        
        print("\nProcessing completed for all files.")    
        
    ind_ = np.arange(len(xRes_mean_all))    
    xRes_mean_all = np.array(xRes_mean_all)
    yRes_mean_all = np.array(yRes_mean_all)
    zRes_mean_all = np.array(zRes_mean_all)
    
       
    xRes_std_all = np.array(xRes_std_all)
    yRes_std_all = np.array(yRes_std_all)
    zRes_std_all = np.array(zRes_std_all)
    
    plt.figure()
    plt.subplot(3, 1, 1)
    plt.plot(ind_, xRes_mean_all, color="k")
    plt.title("X-FWHM")
    plt.xlabel("Time")
    plt.ylabel("FWHM (nm)")
    
    plt.subplot(3, 1, 2)
    plt.plot(ind_, yRes_mean_all, color="k")
    plt.title("Y-FWHM")
    plt.xlabel("Time")
    plt.ylabel("FWHM (nm)")
    
    plt.subplot(3, 1, 3)
    plt.plot(ind_, zRes_mean_all, color="k")
    plt.title("FWHM")
    plt.xlabel("Time")
    plt.ylabel("Resolution (nm)")
    
    plt.tight_layout()
    plt.show()
    
    return final_results_all, xRes_mean_all, yRes_mean_all, zRes_mean_all, xRes_std_all, yRes_std_all, zRes_std_all 

print('Input auto image')
    
auto_final_results_all, xRes_auto, yRes_auto, zRes_auto, xRes_std_auto, yRes_std_auto, zRes_std_auto = process_tiff_files(
    lateral_pixel_size=147,
    axial_pixel_size=147,
    pad=1,
    roi_radius=15,
    threshold_value=1000,
    fit_quality=0.85
)

print('Input without auto image')

without_auto_final_results_all, xRes_without_auto, yRes_without_auto, zRes_without_auto, xRes_std_without_auto, yRes_std_without_auto, zRes_std_without_auto = process_tiff_files(
    lateral_pixel_size=147,
    axial_pixel_size=147,
    pad=1,
    roi_radius=15,
    threshold_value=1000,
    fit_quality=0.85
)



ind_ = np.arange(len(xRes_auto))*2  


plt.figure()
plt.subplot(3, 1, 1)
plt.errorbar(ind_, xRes_auto,yerr= xRes_std_auto, label='With autofocus', color="r")
plt.errorbar(ind_, xRes_without_auto, yerr= xRes_std_without_auto, label='Without autofocus', color="k")
plt.title("X-FWHM")
plt.xlabel("Time (minute)")
plt.ylabel("FWHM (nm)")

plt.subplot(3, 1, 2)
plt.errorbar(ind_, yRes_auto, yerr= yRes_std_auto, label='With autofocus', color="r")
plt.errorbar(ind_, yRes_without_auto, yerr= yRes_std_without_auto, label='Without autofocus', color="k")
plt.title("Y-FWHM")
plt.xlabel("Time (minute)")
plt.ylabel("FWHM (nm)")

plt.subplot(3, 1, 3)
plt.errorbar(ind_, zRes_auto, yerr= zRes_std_auto, label='With autofocus', color="r")
plt.errorbar(ind_, zRes_without_auto,yerr= zRes_std_without_auto,label='Without autofocus', color="k")
plt.title("Z-FWHM")
plt.xlabel("Time (minute)")
plt.ylabel("Resolution (nm)")

plt.tight_layout()
plt.show()
