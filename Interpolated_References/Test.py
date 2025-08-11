import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d 

# Load Data
data = np.loadtxt("MilkProtein.txt")
data2 = np.loadtxt("MilkProtein_Interp.txt")

# Grab wavenumbers and intensities
wavenumbers = data[:, 0]
intensity = data[:, 1]

wn = data2[:, 0]
inten = data2[:, 1]

# Calculate Gaussian
mean = np.mean(wavenumbers)
print(mean)
fwhm = 7
sigma = fwhm/2.35
top = (wavenumbers-mean)**2
gaussian = np.exp(-(top)/(2*(sigma**2)))

# Gaussian v2
dx = np.mean(np.diff(wavenumbers))
fwhm = 7
sigma = fwhm / 2.35 / dx
gaus = gaussian_filter1d(intensity, sigma)


# Interpolate
data_interpolate = interp1d(wavenumbers, gaus, bounds_error=False)

# Get the range of numbers we want
num_range = np.arange(784, 800.5, 0.5)
new_num_range = []
for i in num_range:
    wavenumber = ((1/i) - (1/1031)) * (1 * 10**7)
    new_num_range.append(wavenumber)
new_num_range = new_num_range[::-1]

# Interpolated values
smoothed_intensity = data_interpolate(new_num_range)
# New stack with the smoothed intensity over the range we wanted
gaussian_stack = np.column_stack([new_num_range, smoothed_intensity])

# Plot data
plt.plot(wavenumbers, gaus, label='Original')
plt.plot(wn, inten, label='Interpolated')
plt.xlabel('Wavenumbers')
plt.ylabel('Intensity')
plt.title('Gum')
plt.legend()
plt.show()