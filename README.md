# Bead Image Resolution Analysis

This Python script evaluates the Full Width at Half Maximum (FWHM) of 100 nm fluorescence bead images in the x, y, and z dimensions to estimate the resolution of an imaging system. 
It performs Gaussian fitting on the intensity profiles of individual beads to calculate the lateral and axial resolution and assess the system's imaging quality. See our bioRxiv paper [1] (https://doi.org/10.1101/2024.11.29.626121) for full optical setup instructions.

![OPM-Autofocus-Fig 4](https://github.com/user-attachments/assets/27e0122d-49aa-464d-a5e4-3b5ff4cd2925)

Fig. 1. Time lapse volumetric imaging using OPM with and without remote focus stabilization. A-B 100nm fluorescent nanospheres, as imaged with OPM without remote focus stabilization at the beginning 
and end of an one hour time lapse series. Cross sectional maximum intensity projections are shown. C-D OPM imaging of the same volume as shown in A-B, but with remote focus stabilization active. 
E Three timepoints of a timelapse of A375 cancer cells labeled with tractin-mRuby, using OPM with remote focus stabilization active. Cross sectional maximum intensity projections are shown.
## Features

- **TIFF file processing**: Handles multiple 3D TIFF files with bead images.
- **Thresholding and Particle Detection**: Identifies individual bead particles using a specified intensity threshold.
- **Gaussian Fitting**: Fits a Gaussian function to the intensity profiles in x, y, and z dimensions for each bead.
- **Resolution Calculation**: Calculates FWHM and estimates resolution in nanometers.
- **Quality Filtering**: Filters results based on the fit quality (R² values).
- **Visualization**: Generates scatter plots of resolution values versus lateral positions for x, y, and z dimensions.
- **Summary Statistics**: Computes mean and standard deviation of resolutions across all analyzed beads.

## Dependencies

- Python 3.x
- Libraries: `numpy`, `tifffile`, `scipy`, `matplotlib`, `skimage`, `tkinter`

Install the required libraries using:

```bash
pip install numpy tifffile scipy matplotlib scikit-image
```

## Usage

1. Clone the repository and navigate to the script directory.
2. Run the script:

   ```bash
   python resolution_analysis.py
   ```

3. A file selection dialog will appear. Select one or more 3D TIFF files containing bead images.
4. The script will process the selected files and display resolution scatter plots for x, y, and z dimensions.
5. The results, including mean and standard deviations of resolution, are printed to the console.

## Parameters

The script uses the following adjustable parameters:

- `lateral_pixel_size`: Pixel size in the lateral (x, y) dimensions (default: 147 nm).
- `axial_pixel_size`: Pixel size in the axial (z) dimension (default: 147 nm).
- `pad`: Padding added to image borders (default: 1 pixel).
- `roi_radius`: Radius of the region of interest (ROI) around detected bead centroids (default: 15 pixels).
- `threshold_value`: Intensity threshold for particle detection (default: 1000).
- `fit_quality`: Minimum R² value for Gaussian fitting to accept results (default: 0.85).

Modify these parameters directly in the script to adapt it to your imaging data.

## Output

- **Scatter Plots**: Resolution versus lateral position for x, y, and z dimensions.
- **Console Log**: Summary statistics (mean and standard deviation) of x, y, and z resolutions.
- **Debugging Info**: Coordinates, FWHM values, and quality metrics of individual beads.

## Example

After running the script and selecting files, the scatter plots for x, y, and z resolutions will look similar to the example below:

```
X-Resolution Mean: 148.3 nm, STD: 5.2 nm
Y-Resolution Mean: 150.1 nm, STD: 6.1 nm
Z-Resolution Mean: 400.8 nm, STD: 12.3 nm
```

We have used the code to assess the imaging performance of the autofocusing system, we measured the full width at half maximum (FWHM) of 100 nm diameter fluorescent nanospheres 
and performed image decorrelation analysis on subcellular features. These metrics were used to evaluate resolution and the preservation of structural details. 
Additionally, we compared image quality under conditions with and without the autofocus mechanism to demonstrate its effectiveness in maintaining optimal focal plane alignment during acquisition.

<img width="424" alt="OPM-Autofocus-Fig S3" src="https://github.com/user-attachments/assets/e76a26b8-b536-4374-858a-05700b5ff571">

Fig. 2. Comparison of the XYZ FWHM of 100 nm fluorescence bead with and without autofocusing system.

## Notes

- Ensure the TIFF files contain isotropic bead images, ideally from a standardized calibration sample.
- Adjust the `threshold_value` if beads are not detected or too many artifacts are included.
