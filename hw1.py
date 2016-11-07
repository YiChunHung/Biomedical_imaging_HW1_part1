import csv
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from scipy.signal import hilbert

#parameter
time_offset = 6.48 * (10**-6) # in usec
fs = 50 * (10 ** 6)    # sampling rate, in MHz
aper_size = 6  # aperture size, in mm
#focal_pt = ;   # focal point, in mm
dx = 0.050    # distance between two successive scanning positions, in mm
soundv = 1.54 * (10 ** 6)   # speed of sound, in mm/us

#(a)
DR = 40  #dynamic range in dB

#####read file and switch to list object
points_rf_data = csv.reader(open('points_rf_data','r',newline=''),delimiter=' ',skipinitialspace=True)
rf_list_str = list(points_rf_data)

#####create the array of data
rf_list_num = [[float(col) for col in row] for row in rf_list_str]
rf_array = np.array(rf_list_num)

#####set the x and z axis of image
size = np.shape(rf_array)   #size = (zdm,xdm)
z_axis = np.array(range(0,size[0]))
z_axis = (z_axis/fs+time_offset) * soundv/2
x_axis = np.array(range(0,size[1])) * dx

#####Detect the envolope of the signal
rf_array_trans = np.transpose(rf_array)
envolope = abs(hilbert(rf_array_trans))
envolope_dB = 20 * np.log10(np.transpose(envolope)/np.amax(envolope))

#####Show the image
cmap = cm.get_cmap('Greys_r',40);
plt.imshow(envolope_dB+DR, cmap = cmap, vmax = 40, vmin = 0, 
			extent = [x_axis.min(),x_axis.max(),z_axis.max(),z_axis.min()])
plt.xlabel('Lateral position (mm)')
plt.ylabel('Depth (mm)')
plt.title('Point targets')
plt.colorbar()
plt.show()
	



