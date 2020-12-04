import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from itertools import product

# fname = 'ap_cutoff3475_rap12_ran3'
fname = 'final_scan_cutoff3490_v3'
fitsname = 'A1_pic.fits'
bins_number = 25

data = np.loadtxt(fname, skiprows=1)  # Data structure 0)ypos 1)xpos 2)radius 
Magerror = 2.000E-02
Magzpt = 2.530E+01

image_data = fits.getdata(fitsname, ext=0)
shapes = image_data.shape

def fit_galaxy(ypos, xpos, r_in, r_out = 4):
    """ fit the galaxy to a circle of radius r_in and use the outer ring to calculate the local background """
    count_out = []
    count_in = []
    for j, i in product(np.arange(ypos - (r_out + r_in), ypos + r_out + r_in + 1),np.arange(xpos - (r_out + r_in), xpos + 1 + r_out + r_in)):  # Create square
        if (j - ypos) ** 2 + (i - xpos) ** 2 <= r_in ** 2 and 0<= j <= shapes[0] - 1 and 0<= i <= shapes[1] - 1: # make sure points are in a circle
            j,i = [int(j),int(i)]
            if 3419<image_data[j,i]:
                count_in.append(image_data[j,i])
        if r_in ** 2 < (j - ypos) ** 2 + (i - xpos) ** 2 <= (r_in + r_out)**2 and 0<= j <= (shapes[0] - 1) and 0<= i <= shapes[1] - 1: # in the outer ring
            j,i = [int(j),int(i)]
            if 3419 - 13.5 <image_data[j,i]< 3419 + 13.5:       # Make sure it is within 1 sd away from the mean 
                count_out.append(image_data[j][i]) 
    n = len(count_in)
    count_in = np.array(count_in).sum()
    if len(count_out)>4:
        count_out = np.array(count_out).sum() / len(count_out)    
    else:
        count_out = 3419
    return count_in, count_out, n

count = np.zeros(len(data))
for c in range(len(data)):
    co = fit_galaxy(data[c,0], data[c,1], data[c,2])
    count[c] = co[0] - co[1] * co[2]
       
counts = np.array(count)
magnitude = Magzpt - 2.5 * np.log10(counts)
magnitudes = []
for i in range(len(magnitude)):
    if magnitude[i] < 18:
        magnitudes.append(magnitude[i])
    else:
        print(magnitude[i])
magnitudes = np.array(magnitudes)

numbers, xbin = np.histogram(magnitudes, bins = bins_number)
number = []

# Get cumulative sum
for i in range(len(numbers)):
    number.append(numbers[: (i+1)].sum())
number = np.array(number)

# print("number", number, "xbin", xbin)
x = np.linspace(xbin.min(), xbin.max(), bins_number)
log_n = np.array(np.log10(number))

# Plot fitting linear
plt.figure(1)
plt.errorbar(x[1:], log_n[1:], xerr= Magerror, yerr = np.array(1/np.sqrt(number))[1:], fmt = 'x')

coeff, cov = np.polyfit(x, log_n, 1, cov=True)
formula = np.poly1d(coeff)
print("coeff =", coeff[0], "coeff error =", cov[0])
plt.plot(x, formula(x))

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("magnitude", fontsize=16)
plt.ylabel("log10 (N)", fontsize=16)
plt.title("log(N) against magnitude cut off at 3490 counts", fontsize=20)

# Plot histogram
plt.figure(2)
plt.hist(magnitudes, bins = bins_number)
plt.xlabel("magnitude", fontsize=16)
plt.ylabel("Number of galaxies", fontsize=16)
plt.title("Histogram of the galaxies", fontsize=20)
plt.show()
