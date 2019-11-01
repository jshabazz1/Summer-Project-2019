# For querying for data
import requests
from astroquery.mast import Tesscut, Catalogs

# For manipulating data
import numpy as np

from astropy.table import Table
from astropy.coordinates import SkyCoord

import re

# For matplotlib plotting
import matplotlib

import matplotlib.pyplot as plt
import matplotlib.animation as animation

# For animation display
from matplotlib import rc
from IPython.display import HTML
rc('animation', html='jshtml')

# For bokeh plotting
from bokeh import plotting


ticid = 167602025
istart = 1616.0
istop = 1620.0



def my_animation (ticid, istart, istop, flare=None):
	Writer = animation.writers['imagemagick']
	Writer = Writer(fps=20, metadata=dict(artist='Me'), bitrate=1800)

	starName = "TIC " + str(ticid)
	radSearch = 4/60 #radius in degrees

	#Querying RA and DEC from Catalogs given a TIC ID
	catalogData = Catalogs.query_object(starName, radius = radSearch, catalog = "TIC")
	ra = catalogData[0]['ra']
	dec = catalogData[0]['dec']

	#Gets Skycoord given RA and DEC then calls to TESScut
	obj_coord = SkyCoord(ra,dec,unit="deg")
	print(obj_coord)
	Tesscut.get_sectors(obj_coord)
	cutout_hdu = Tesscut.get_cutouts(obj_coord, size=20)[0]
	cutout_hdu.info()
	cutout_table = cutout_hdu[1].data
	cutout_table.columns


	def find_index(btjd):
		return (np.abs(cutout_table['TIME']-btjd)).argmin()

	start = find_index(istart)
	end = find_index(istop)

	print(f"Frames {istart}-{istop} ({istop-istart} frames)")

	def make_animation(data_array, start_frame=start, end_frame=end, vmin=None, vmax=None, delay=50):
	    """
	    Function that takes an array where each frame is a 2D image array and make an animated plot
	    that runs through the frames.
	    
	    Note: This can take a long time to run if you have a lot of frames.    
	    Parameters
	    ----------
	    data_array : array
	        Array of 2D images.
	    start_frame : int
	        The index of the initial frame to show. Default is the first frame.
	    end_frame : int
	        The index of the final frame to show. Default is the last frame.
	    vmin : float
	        Data range min for the colormap. Defaults to data minimum value.
	    vmax : float
	        Data range max for the colormap. Defaults to data maximum value.
	    delay: 
	        Delay before the next frame is shown in milliseconds.

	    Returns
	    -------
	    response : `animation.FuncAnimation`
	    """
	    
	    if not vmin:
	        vmin = np.min(data_array)
	    if not vmax:
	        vmax = np.max(data_array)
	        
	    if not end_frame:
	        end_frame = len(data_array) - 1 # set to the end of the array
	        
	    num_frames = end_frame - start_frame + 1 # include the end frame
	        
	    def animate(i, fig, ax, binarytab, start=start):
	        """Function used to update the animation"""
	        ax.set_title("Epoch #" + str(i+start))
	        fig.suptitle(f'TIC {ticid}')
	        im = ax.imshow(binarytab[i+start], cmap=plt.cm.YlGnBu_r, vmin=vmin, vmax=vmax)
	        return im
	    
	    # Create initial plot.
	    fig, ax = plt.subplots(figsize=(10,10))
	    ax.imshow(data_array[start_frame], cmap=plt.cm.YlGnBu_r, vmin=vmin, vmax=vmax)

	    ani = animation.FuncAnimation(fig, animate, fargs=(fig, ax, data_array, start_frame), frames=num_frames, 
	                                  interval=delay, repeat_delay=1000)
	    
	    ani.save('/Users/jshabazz/Work/TESScut_anims/' + str(ticid) +'_flareevent'+ str(flare) +'.gif', writer=Writer);print('File created')
	    plt.close()

	    return ani
	make_animation(cutout_table['FLUX'], vmax=500)

my_animation(ticid, istart, istop)




