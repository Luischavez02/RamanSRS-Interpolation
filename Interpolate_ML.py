import numpy as np
from scipy.interpolate import interp1d 
import os
from scipy.ndimage import gaussian_filter1d

def Gaussian(subtype_file_name, input_folder, fwhm=7):
    # Load data
    filepath = os.path.join(input_folder, subtype_file_name)

    # Load Data
    data = np.loadtxt(filepath)
    wavenumbers = data[:, 0]
    intensity = data[:, 1]

    # Calculate Gaussian
    dx = np.mean(np.diff(wavenumbers))
    fwhm = 7
    sigma = fwhm / 2.35 / dx
    gaus = gaussian_filter1d(intensity, sigma)
    
    # Make the wavenumbers into a column stack
    stack = np.column_stack([wavenumbers, gaus])
    return stack

def interp(num_range, stack, file_name, output_folder):
    # Loads Wavenumbers and Intensities
    wavenumbers = stack[:, 0]
    intensity = stack[:, 1]
    # Interpolate function
    data_interpolate = interp1d(wavenumbers, intensity, bounds_error=False, fill_value="extrapolate")
    # New intensity after interpolating
    new_intensity = data_interpolate(num_range)
    # Directory
    os.makedirs(output_folder, exist_ok=True)
    output_path = os.path.join(output_folder, f"{file_name}.txt")
    # Save file
    np.savetxt(output_path, np.column_stack([num_range,new_intensity]), delimiter=' ', fmt='%.4f')
    

# Output & Input Folder
output_folder = "./Interpolated_References/Interpolated"

input_folder = "Interpolated_References/Subtype_Reference"

# Steps and number range
num_range = np.arange(784, 800.5, 0.5)
new_num_range = []
for i in num_range:
    wavenumber = ((1/i) - (1/1031)) * (1 * 10**7)
    new_num_range.append(wavenumber)
new_num_range = new_num_range[::-1]

# Run functions
gaussian_stack = Gaussian("MilkProtein.txt",input_folder)
interp(new_num_range, gaussian_stack, "MilkProtein_Interp", output_folder)