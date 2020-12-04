import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from itertools import product
import copy
import time

class galaxy:
    def __init__(self, ypos, xpos, radius, count_sum, bg_galaxy=3419, no_count = 0):  ### xdim, ydim??? and Angle
        self.ypos = ypos
        self.xpos = xpos
        self.radius = radius
        self.count_sum = count_sum
        self.bg_galaxy = bg_galaxy
        self.no_count = no_count
        # Maybe some loss function to calculate counting

    def get_info(self):
        """ Return x, y position, radius, sum of pixels and local background"""
        return np.array([self.ypos, self.xpos, self.radius, self.count_sum, self.bg_galaxy, self.no_count])
    
class analysis:
    def __init__(self, fname):
        self.fname = fname
        self.raw_image_data = fits.getdata(self.fname, ext=0)###[3820:3975,422:570]
        self.shapes = [self.raw_image_data.shape[0], self.raw_image_data.shape[1]]
        self.dimension = self.shapes[0] * self.shapes[1]
        # creating masked array
        self.masked = np.zeros(self.shapes)
        for j, i in product(np.arange(0, self.shapes[0]), np.arange(0, self.shapes[1])):
            self.masked[j, i] = 1
        # ranking the image_data from the highest counts to the lowest
        self.image_data = copy.deepcopy(self.raw_image_data)  # Creating a deep copy of the original data, not advisable for a large data file
        self.image_data = self.image_data.reshape([self.shapes[0] * self.shapes[1]])  # ranked image
        self.rank = np.argsort(self.image_data)[::-1]         # the corresponding position of the data can be found from the rank
        self.image_data.sort(axis=0)
        self.image_data = self.image_data[::-1]
        self.galaxies = []
        print(self.raw_image_data[116,79])


    def rank_yx(self, rankyx, rank_to_yx=1):
        """ Transform rank (single value) to x,y position if rank_to_yx = 1 (True)"""
        if rank_to_yx == 1:
            x = int(rankyx) % int(self.shapes[1])
            y = (rankyx - x) / int(self.shapes[1])
            return [y, x]  # More convenient to return y, x
        
        if rank_to_yx == 0:  # that means transform yx to rank, expecting rankyx to be a list, may not be necessary
            rankyx = rankyx[0] * int(self.shapes[1]) + rankyx[1]
            return rankyx  # returns back a float

    def mask_region(self, ypos, xpos, r):
        """ masking the circular region at radius r"""
        for j, i in product(np.arange(ypos - r, ypos + r + 1), np.arange(xpos - r, xpos + 1 + r)):  # Create square
            if (j - ypos) ** 2 + (i - xpos) ** 2 <= r ** 2 and 0 <= j<= self.shapes[0] - 1 and 0<= i <=self.shapes[1] - 1:
                j = int(j)
                i = int(i)
                self.masked[j, i] = 0

    def pick_largest(self, cut_off):
        """ Look for the largest unmasked value"""
        for i in range(self.dimension):
            m = self.masked[int(self.rank_yx(self.rank[i])[0])      # locating the corresponding mark array
                            ,int(self.rank_yx(self.rank[i])[1])]
            if m * self.image_data[i] == self.image_data[i]:
                if self.image_data[i] <= cut_off:
                    print("Surveying completed")
                    return -1,-1  # returns -1,-1 if scan is completed
                else:
                    return self.image_data[i], np.array(self.rank[i])

    def fit_galaxy(self, ypos, xpos, r_in, r_out = 0):
        """ fit the galaxy to a circle of radius r """
        count_out = []
        count_in = []
        for j, i in product(np.arange(ypos - (r_out + r_in), ypos + r_out + r_in + 1),np.arange(xpos - (r_out + r_in), xpos + 1 + r_out + r_in)):  # Create square
            if (j - ypos) ** 2 + (i - xpos) ** 2 <= r_in ** 2 and 0<= j <= self.shapes[0] - 1 and 0<= i <= self.shapes[1] - 1: # make sure points are in a circle
                j = int(j)
                i = int(i)
                if self.raw_image_data[j,i] * self.masked[j,i] == self.raw_image_data[j,i]:
                    count_in.append(self.raw_image_data[j,i])
                self.masked[j,i] = 0        # self.mask_region runs the for loop again
            if r_in ** 2 < (j - ypos) ** 2 + (i - xpos) ** 2 <= (r_in + r_out)**2 and 0<= j <= (self.shapes[0] - 1) and 0<= i  <= self.shapes[1] - 1: # in the outer ring
                j = int(j)
                i = int(i)
                if self.raw_image_data[j,i] * self.masked[j,i] == self.raw_image_data[j,i]:        
                    count_out.append(self.raw_image_data[j][i]) 
                self.masked[j,i]
        return count_in, count_out

    def scan_ap(self, cut_off = 3480, r_ap = 12, r_an = 3):
        """ Do the aperture method scan by setting radius for the aperture and the annulus"""
        for trial in range(self.dimension):
            max_count, max_rank = self.pick_largest(cut_off = cut_off)
            if max_count >= 0:
                y,x = self.rank_yx(max_rank)
                print("Scan  pos", y,x," scanning",trial,"counts", max_count)
                count_in, count_out = self.fit_galaxy(y,x,r_ap, r_an)
                count_sum = []
                local_bg = []
                for c in range(len(count_in)):
                    if count_in[c] >= cut_off:
                        count_sum.append(count_in[c])
                for c in range(len(count_out)):
                    if count_out[c] <= cut_off:
                        local_bg.append(count_out[c])
                no_c = len(count_in)
                if len(count_sum) >= int(np.pi * (r_ap **2) / 2):  # Make sure it is not noise
                    count_sum = np.array(count_sum).sum()
                    if len(local_bg) != 0:
                        total = 0
                        for c in range(len(local_bg)):
                            if 3*13.8 <= abs(local_bg[c] - 3419):
                                total += local_bg[c]
                        local_bg = total / len(local_bg)
                    else:
                        local_bg = 3419
                    print("galaxy founded at ", y, x)
                    self.galaxies.append(galaxy(y, x, r_ap, count_sum, bg_galaxy=local_bg, no_count = no_c))
           
            elif max_count == -1:
                print("aperture scan completed, number of galaxies found is", len(self.galaxies))
                break



###########################################################################################################

    def scan(self, cut_off, r0 = 4, dr = 2, rmax = 80, bg_global = 3419, sigma = 3.5 * 13.8):
        """ A variation on the aperture method.The idea is to
        start at radius 2 and increase the radius to get a better fit on the galaxy"""
        
        start_time = time.time()
        for trial in range(self.dimension):
            
            max_count, max_rank = self.pick_largest(cut_off = cut_off)
            fitted = 0                   # galaxy not fitted to a circle
            point = []                   # 
            circle_number = 0
            
            if max_count == -1:
                print("Scan completed, number of galaxies found is ", len(self.galaxies), "run time is", time.time() - start_time)
                break
            
            if max_count >= 0:      # That means a value that is larger than cut_off exists
                ypos,xpos = self.rank_yx(max_rank)
                # print("max_count, y, x", max_count, ypos, xpos)    
            
            
            for r in range(r0, rmax, dr):     # r = radius, we know the largest radius can't be >80 by inspecting the pic.
                print("locating the galaxy position at", ypos, xpos, "at a radius r =", r)               
                if fitted == 1 or r == 80:
                    # print("max_count, yx",max_count, ypos, xpos, "cut =", no_cut, len(new_point)/2)
                    
                    self.mask_region(ypos, xpos, r - dr)
                    # print(bg_local, no_bg)
                    if r - dr > r0:         # it must be a galaxy
                        
                        # getting the local bg
                        if no_bg > 3:   # check if enough data to deduce bg
                            bg = bg_local / no_bg   
                        else:
                            bg = bg_global
                        
                        if circle_number >= np.pi* 2* (r - dr)**2:
                            self.galaxies.append(galaxy(ypos, xpos, r - dr, np.array(point).sum(), bg, circle_number))

                    else:
                        if no_cut < circle_number/2 and circle_number*2 >= np.pi * (r0**2):
                            self.galaxies.append(galaxy(ypos, xpos, r - dr, np.array(point).sum(), bg, circle_number))

                    print("\ngalaxy scan completed,radius =", r - dr, "position =", ypos, xpos, "max count =", max_count, "time = ", time.time() - start_time)

                    break
                
                ##### Resetting parameters ####
                no_bg = 0              # number 0f background pixels
                bg_local = 0           # sum of local background noises
                no_cut = 0             # number of counts below cut off
                new_point = []         # pending pixels to be added
                ###############################
                
                if fitted == 0:
                    
                    if r == r0:
                        
                        for j, i in product(np.arange(ypos - r, ypos + r + 1), np.arange(xpos - r, xpos + 1 + r)):  # Create square
            
                            # Check if it is within the circle radius = 2
                            if int((i - xpos) ** 2 + (j - ypos) ** 2) <= r ** 2 and 0<= j <= (self.shapes[0] - 1) and 0<= i  <= self.shapes[1] - 1:
                                
                                i,j =[int(i), int(j)]
                                if self.raw_image_data[j,i] == self.raw_image_data[j,i] *self.masked[j,i]:    # Append the ppoint if not masked (masked = 1)

                                    if self.raw_image_data[j,i] <= cut_off:
                                        no_cut += 1

                                    if self.raw_image_data[j,i] > cut_off:
                                        point.append(self.raw_image_data[j, i])
                                        circle_number += 1

                                    if abs(self.raw_image_data[j,i] - bg_global) <= sigma:
                                        bg_local += self.raw_image_data[j,i]
                                        no_bg += 1
                                        
                        if no_cut >= len(point)/2 or 2*circle_number < np.pi * r0**2:
                            fitted = 1
                        else:
                            pass
                       
#######################################################################################################
                    if r > r0:
                        for j, i in product(np.arange(ypos - r, ypos + r + 1), np.arange(xpos - r, xpos + 1 + r)):              
                            
                            # Check if data are in between the previous and the new circle
                            if (r - dr)**2 < int((i - xpos) ** 2 + (j - ypos) ** 2) <= r ** 2 and 0<= j <= (self.shapes[0] - 1) and 0<= i  <= self.shapes[1] - 1:
                                i,j =[int(i), int(j)]               # just incase    
                                
                                if self.raw_image_data[j,i] * self.masked[j,i] == self.raw_image_data[j,i]:

                                    if self.raw_image_data[j,i] <= cut_off:
                                        no_cut += 1
                                    if self.raw_image_data[j,i] > cut_off:
                                        new_point.append(self.raw_image_data[j, i])  # points are pending to be added in
                                        circle_number += 1

                                    if abs(self.raw_image_data[j,i] - bg_global) <= sigma:
                                        bg_local += self.raw_image_data[j,i]
                                        no_bg += 1
                        
                        # Check if half of the new data points are inside cut off region
                        if no_cut <= int(len(new_point)/2) or circle_number*2 < np.pi * r**2:
                            for rannk in range(len(new_point)):
                                point.append(new_point[rannk])

                        else:
                            fitted = 1

#######################################################################################################

start_time = time.time()
ana = analysis("A1_pic_bleed_filtered.fits")

# ana.scan_ap(cut_off= 3490, r_ap = 12)
ana.scan(3490)

#%%
fname_ap = "test_ap_cutoff3490_rap12_ran3"
fname_sc = "final_scan_cutoff3490_v3"
print("Run time =",time.time() - start_time)

info = []
for g in range(len(ana.galaxies)):
    info.append(ana.galaxies[g].get_info())
info = np.array(info)

# np.savetxt(fname_ap, info, header = "ypos xpos radius count_sum bg_galaxy no_count")
np.savetxt(fname_sc, info, header = "ypos xpos radius count_sum bg_galaxy no_count")