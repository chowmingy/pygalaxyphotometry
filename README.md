# pygalaxyphotometry
This code is written for the A1 image processing for physics year 3 project but the code can be adapted for various photometry by changing the parameters.



Uploaded file: 23/11/20

# Remember to install astropy before running the code.


Example:
________________________________________________________________


from scan_and_fit_v3.1 import analysis, galaxy

ana = analysis("fits_file_name")

Choose one:
ana.scan_ap(cut_off= 3475, r_ap = 12) 	# Apecture method
ana.scan(3475)				# Varting apecture method

fname = "file_name"

info = []
for g in range(len(ana.galaxies)):
    info.append(ana.galaxies[g].get_info())	# Getting the information of galaxies
info = np.array(info)
np.savetxt(fname, info, header = "ypos xpos radius count_sum bg_galaxy no_count")
