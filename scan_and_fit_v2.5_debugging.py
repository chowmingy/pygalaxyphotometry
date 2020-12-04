import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from itertools import product
import time
        
############ Debuugging  #############
def rank_yx(rankyx, rank_to_yx = 1):
    """ Transform rank (single value) to x,y position if rank_to_yx == 1 (True)"""
    if rank_to_yx == 1:
        x = int(rankyx) % int(shapes[1])
        y = (rankyx - x) / int(shapes[1])
        return [y, x]                   # More convenient to return y, x
    if rank_to_yx == 0:                 # that means transfer yx to rank
        rankyx = rankyx[0] * int(shapes[1]) + rankyx[1]
        return rankyx                   # returns back to array

def pick_largest():
    """ Look for the largest unmasked value"""
    for i in range(len(image_data)):
        if masked[i] * image_data[i] == image_data[i]:
            if image_data[i]<= cut_off:
                print("Surveying completed")
                survey_completed = True
                break
            else:
                return image_data[i], np.array(rank[i])
            
def mask_region(ypos, xpos, r):
    """ Given an array of positions, the function masks the corresponding region (a circle centered at position [xpos, ypos] with radius r"""

    circle_rank = []
    # Get
    for j, i in product(np.arange(ypos - r, ypos + r + 1), np.arange(xpos - r, xpos + 1 + r)):  # Create square
        if (j - ypos)**2 + (i - xpos)**2 <=r**2:
            circle_rank.append(rank_yx([j, i], rank_to_yx = 0))

    for temp in range(len(circle_rank)):  # temporary
            for dim in range(len(rank)):
                if circle_rank[temp] == rank[dim] * masked[dim]:  # Locating the specific rank
                    masked[dim] = 0
                    break
        
data = []        
##########################################
def fit_galaxy(ypos, xpos):

    """ Put in the location and attempt to fit the galaxy with a circle, stops when more than half of the newly added region are bg,
    Count it as counting error or weird bg if most surrounding is in the bg region for radius <= 2 """

    fitted = 0
    circle_rank = []
    point = []                   # list of pixels within the region "galaxy"
    start_time = time.time()

    dimension = shapes[0]*shapes[1]

    
    for r in range(2, 500,2):    # r = radius, we know the largest radius can't be >500 as it will be 1/4 of the pic.
        print("locating the galaxy position at", ypos, xpos, "at a radius r =", r)
        
        if fitted == 1 or r == 80:
            mask_region(ypos, xpos, r - 2)
            if r - 2 > 2: # Too small so assume counting error
                data.append([ypos, xpos, r -2, np.array(point).sum(), bg_local / no_bg])
            print("run time = ", time.time() - start_time)
            print("\nscan of galaxy completed at radius", r - 2, "position =", ypos, xpos)
            break

        # Resetting parameters
        no_bg = 0              # number 0f background pixels
        bg_local = 0           # sum of local background noises
        new_rank = []
        new_point = []         # pending pixels to be added
    
        if fitted == 0:
            
            if r == 2:
                for j, i in product(np.arange(ypos - r, ypos + r + 1), np.arange(xpos - r, xpos + 1 + r)):  # Create square
    
                    # Check if it is within the circle radius = 2
                    if int((i - xpos) ** 2 + (j - ypos) ** 2) <= r ** 2:
                        circle_rank.append(rank_yx([j, i], rank_to_yx = 0))
                        
                # Append the region to the circle_rank
                for temp in range(len(circle_rank)):                                 # temporary data pending to be added in
                    for dim in range(dimension):
                        if circle_rank[temp] == rank[dim] * masked[dim]:             # i.e. not masked
                            point.append(image_data[dim])

                            # print("dim = ", dim, "rank = ", rank[dim], "pixel1 =", image_data[dim])
                            if image_data[dim] <= cut_off:
                                    no_bg += 1
                                    bg_local += image_data[dim]

                if no_bg >= int(len(circle_rank))/2:
                    fitted = 1
                else:
                    pass
#########################################################
            if r > 2:
                # print("largest pixel =", image_data[0:5], "Corresponding rank", rank[0:5])
                for j, i in product(np.arange(ypos - r, ypos + r + 1), np.arange(xpos - r, xpos + 1 + r)):
                    # print(j,i)
                    
                    # Check if data are in between the previous and the new circle
                    if (r - 2)**2 < int((i - xpos) ** 2 + (j - ypos) ** 2) <= r ** 2:
                        new_rank.append(rank_yx([j, i], rank_to_yx = 0))
                # print("largest pixel =", image_data[0:5], "Corresponding rank", rank[0:5])        
                # Getting the temporary new_point
                for temp in range(len(new_rank)):  # temporary
                    for dim in range(dimension):
                        if new_rank[temp] == rank[dim] * masked[dim]:  # i.e. not masked

                            new_point.append(image_data[dim])
                            if image_data[dim] <= cut_off:
                                
                                # print("dim = ", dim, "rank = ", new_rank[temp], "pixel =", image_data[dim])
                                no_bg += 1
                                bg_local += image_data[dim]
                                
                # Check if half of the new data points are inside cut off region
                if no_bg < int(len(new_rank))/2:
                    for rannk in range(len(new_point)):
                        circle_rank.append(new_rank[rannk])
                        point.append(new_point[rannk])
                else:
                    fitted = 1
                    
                
            
        print("run time = ", time.time() - start_time)

##########################################
        
fname = "A1_pic_bleed_filtered.fits"
image_file = fits.open(fname)
whole_image_data = fits.getdata(fname, ext=0)

rank = np.argsort(image_data)[::-1]
print(rank[0])
############ Change Boundary #########
# ylb = 0
# yub = 150
### Left side 0<x<1435, rhs 1435<x<2610

cut_off = 3440
data_with_actual_position = []

for y,x in product(np.arange(0,2), np.arange(0,1)):
    print("At region", y, x)
    image_data =whole_image_data[y*150:(y+1)*150,x*150:(x+1)*150]
    shapes = [image_data.shape[0],image_data.shape[1]]
    # for i,j in product(np.arange(0,shapes[0]-1), np.arange(0,shapes[1]-1)):
    #     if image_data[i][j]<3420:
    #         image_data[i][j] = 1
    
    masked = np.zeros(shapes[0] * shapes[1])
    for i in range(len(masked)):
        masked[i] = 1
    
    image_data = image_data.reshape([shapes[0] * shapes[1]])
    rank = np.argsort(image_data)[::-1] 
    rank=np.array(rank)
    image_data.sort(axis=0)
    image_data = image_data[::-1] 
    # print("largest pixel =", image_data[0:5], "Corresponding rank", rank[0:5])
    for i in range(100000):         # doesmn't matter as long as it is large
        a = pick_largest()
        if a == None:
            break
        else:
            fit_galaxy(rank_yx(a[1])[0],rank_yx(a[1])[1])
            data_with_actual_position.append([y,x])
            
    
#%%
np.savetxt("galaxies_cut_off_at_3440",data)
np.savetxt("data_with_actual_region_3440", data_with_actual_position)
