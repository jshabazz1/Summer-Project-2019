# %%
'''
# Finding Flares and Variable Stars: A MAST Case Study

In this tutorial we will learn about MAST's programmatic tools for accessing TESS time series data while exploring a flaring star from the literature.  We will view the flaring star and expore its immediate neighborhood using the MAST API in Python to access and view TESS time series, FFI data, and catalog entries.

Topics to be covered include:
- Using the MAST API to get mission pipeline and TASOC light curves
- Plotting TESS light curves in Python
- Using the MAST API to make an FFI cutout
- Creating a movie of TPF frames in Python
- Using the MAST API to get a list of TESS Input Catalog (TIC) sources
- Over plotting TIC sources on TESS images

See the __[MAST TESS site](http://archive.stsci.edu/tess/)__ for more information and examples of how to access and use TESS data.
'''

# %%
'''
## Terminology

- **TESS:** The Transiting Exoplanet Survey Satellite
- **TASOC:** The TESS Asteroseismic Science Operations Center
- **Sector:** TESS observed the sky in regions of 24x96 degrees along the southern, then northern, ecliptic hemispheres. Each of these regions is referred to as a "sector", starting with Sector 1.
- **TIC:** The TESS input catalog.
- **FFI:** TESS periodically reads out the entire frame of all four cameras, nominally every 30 minutes, and stores them as full frame images (FFIs). 
- **HDU:** Header Data Unit. A FITS file is made up of HDUs that contain data and metadata relating to the file. The first HDU is called the primary HDU, and anything that follows is considered an "extension", e.g., "the first FITS extension", "the second FITS extension", etc.
- **HDUList:** A list of HDUs that comprise a fits file.
- **BJD:** Barycentric Julian Date, the Julian Date that has been corrected for differences in the Earth's position with respect to the Solar System center of mass.
- **BTJD:** Barycentric TESS Julian Date, the timestamp measured in BJD, but offset by 2457000.0. I.e., BTJD = BJD - 2457000.0
- **WCS:** World Coordinate System, the coordinates that locate an astronomical object on the sky. 
'''

# %%
'''
## Imports

In this tutorial we will use the MAST module of Astroquery to query and download data.

We will use both the matplotlib and bokeh packages to visualize our data as they have different strengths and weaknesses.
'''

# %%
# For querying for data
from astroquery.mast import Tesscut, Observations, Catalogs

# For manipulating data
import numpy as np

from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.timeseries import LombScargle
from astropy.time import Time
import astropy.units as u

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
from bokeh.models import Span


# %%
'''
## Exploring a stellar flare

### Selecting the flare

We will start with a known flare from the literature, in this case from [G&uuml;nther, M. N., Zhan, Z., Seager, S.,
et al. 2019, arXiv e-prints, arXiv:1901.00443](https://arxiv.org/abs/1901.00443). We picked a particularly long flare to give use the best chance of finding it in the half hour cadence data as well as the 2 minute cadence data.

We've made note of the TIC ID and sector for our flare of interest, as well as its peak time in BJD format:
'''

# %%
tic_id = 141914082
sector = 1

tpeak = 2458341.89227 # Julian Day

# %%
'''
### Querying MAST

#### Mission light curves

Here we choose the TESS mission (`obs_collection`) and query on our TIC ID and sector.
'''

# %%
mission_res = Observations.query_criteria(obs_collection="TESS", 
                                          target_name=tic_id, 
                                          sequence_number=sector)
mission_res

# %%
'''
#### TASOC light curves

MAST also hosts a variety of community contributed High Level Science Products (HLSPs), all of which are given the mission "HLSP". In this case we will specifically search for HLSPs in the TESS project, which will return the light curves provided by the TASOC (note the `provenance_name` of "TASOC").
'''

# %%
tasoc_res = Observations.query_criteria(target_name=tic_id, 
                                        obs_collection="HLSP", 
                                        project="TESS",
                                        sequence_number=sector)
tasoc_res['dataproduct_type',"obs_collection","target_name","t_exptime","filters",
          "provenance_name","project","sequence_number","instrument_name"]

# %%
'''
In this case there are two light curves, to understand the difference between the two light curves we look to the `t_exptime` column, and note the different values. These exposure times correspond to 2 minutes (short cadence) and 30 minutes (long cadence). We will explore both light curves.
'''

# %%
'''
### Downloading the data products

From here on we will work with the TASOC light curves only, although we could do the same with the mission pipeline light curves as well.

#### Querying for the list of associated data products

Each observation may have one or more associated data products. In the case of the TASOC light curves, there is simply a single light curve file for each observation. 
'''

# %%
tasoc_prod = Observations.get_product_list(tasoc_res)
tasoc_prod["dataproduct_type", "description", "dataURI", "size"]

# %%
'''
#### Downloading files

We can choose do download some or all of the associated data files, in this case since we just have the two light curves, we will download all of the products.
'''

# %%
tasoc_manifest = Observations.download_products(tasoc_prod)
tasoc_manifest

# %%
'''
### Plotting the light curves

We will use bokeh for plotting so that we can have interactivity, and will plot both the 2 minute and 30 mintue cadence data. 

We can tell which is which by examining the filenames and noting that one contains `c0120` (2 min) and the other `c1800` (30 min).
'''

# %%
# Loading the short cadence light curve
hdu = fits.open(tasoc_manifest["Local Path"][0])
short_cad_lc = Table(hdu[1].data)
hdu.close()

# Loading the long cadence light curve
hdu = fits.open(tasoc_manifest["Local Path"][1])
long_cad_lc = Table(hdu[1].data)
hdu.close()

# %%
bfig = plotting.figure(plot_width=850, plot_height=300, title=f"Detrended Lightcurve (TIC{tic_id})")

# Short cadence
bfig.circle(short_cad_lc["TIME"],short_cad_lc["FLUX_RAW"], fill_color="black",size=2, line_color=None)
bfig.line(short_cad_lc["TIME"],short_cad_lc["FLUX_RAW"], line_color='black')

# Long cadence
bfig.circle(long_cad_lc["TIME"],long_cad_lc["FLUX_RAW"], fill_color="#553be7",size=6, line_color=None)
bfig.line(long_cad_lc["TIME"],long_cad_lc["FLUX_RAW"], line_color='#553be7')

# Marking the flare (tpeak is in BJD, while the time column in the light curve is BTJD, so we must convert)
vline = Span(location=(tpeak - 2457000), dimension='height', line_color='#8c0051', line_width=1)
bfig.renderers.extend([vline])

# Labeling the axes
bfig.xaxis.axis_label = "Time (BTJD)"
bfig.yaxis.axis_label = "Flux"

plotting.show(bfig)

# %%
'''
### Making a video

Looking at the above plot we can see the flare event in both the long and short cadence light curves. Since we can see it even in the half hour cadence data, we should be able to make an animation of the area around the flaring star and see the flare happen.

We will use TESScut, the MAST cutout tool for full-frame images to cutout the area around the flaring star across the entire sector, and then make a movie that shows how it changes over time.

We will use the `astroquery.mast` __[Tesscut](https://astroquery.readthedocs.io/en/latest/mast/mast.html#tesscut)__ class to make this cutout.  
We will use two functions:
- Find the sky coordinate of our flare star: `Observations._resolve_object`\*
- Query for cutouts and get the result as a list of HDUList objects: `Tesscut.get_cutouts` \*\*

\* `Observations._resolve_object` is a private (not documented) function which is being removed in favor of the public function `Observations.resolve_object` in the next Astroquery release.

\*\* We must start by finding the sky coordinate of our star, however starting with the next Astroquery release, `Tesscut` functions will be able to take an object name such as a TIC ID as well.
'''

# %%
coord = Observations._resolve_object(f"TIC {tic_id}")

# %%
'''
**Requesting a cutout target pixel file. **

This query will return a list of `HDUList` objects, each of which is the cutout target pixel file for a single sector. In this case, because we specified a single sector we know that the resulting list will only have one element and can pull it out directly.
'''

# %%
cutout_hdu = Tesscut.get_cutouts(coordinates=coord, size=40, sector=1)[0]

# %%
cutout_hdu.info()

# %%
cutout_table = Table(cutout_hdu[1].data)
cutout_table.columns

# %%
'''
#### Exploring the cutout time series

We want to explore what is happening with in our cutout area over the time that the flare occurs, so we will make an animated plot of the cutout frames.

We can't make a movie of the whole sector (it would take too long), so we will choose only the time range around the flare.
'''

# %%
def find_index(btjd):
    """
    Given a time as a Barycentric TESS Julian Date (BTJD) timestamp, return the closest index in a table
    that is assumed to have a TIME column that is also in BTJD"""
    
    return (np.abs(cutout_table['TIME'] - btjd)).argmin()

# %%
start = find_index(1341.5)
end = find_index(1342.5)

print(f"Frames {start}-{end} ({end-start} frames)")

# %%
'''
#### Looking at the animated cutout
'''

# %%
def make_animation(data_array, start_frame=0, end_frame=None, vmin=None, vmax=None, delay=50):
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
        
    def animate(i, fig, ax, binarytab, start=0):
        """Function used to update the animation"""
        ax.set_title("Epoch #" + str(i+start))
        im = ax.imshow(binarytab[i+start], cmap=plt.cm.YlGnBu_r, vmin=vmin, vmax=vmax)
        return im,
    
    # Create initial plot.
    fig, ax = plt.subplots(figsize=(10,10))
    ax.imshow(data_array[start_frame], cmap=plt.cm.YlGnBu_r, vmin=vmin, vmax=vmax)

    ani = animation.FuncAnimation(fig, animate, fargs=(fig, ax, data_array, start_frame), frames=num_frames, 
                                  interval=delay, repeat_delay=1000)
    
    plt.close()
    
    return ani

# %%
make_animation(cutout_table['FLUX'], start, end, vmax=700, delay=150)

# %%
'''
We can see three things in this plot:
- The flare that occures in frames 740-743
- An abberition that appears in frame 754
- A variable star pulsing in the lower right corner
'''

# %%
'''
## Exploring the variable star

Now we will look more closely at the variable star we can see in the animation. 

### Querying the TESS Input Catalog

To start with we will overlay the nearby TIC sources onto the image so we can identify the star in question. To do this we will use the `astroquery.mast` Catalog clas to search the TIC.
'''

# %%
sources = Catalogs.query_object(catalog="TIC", objectname=f"TIC {tic_id}", radius=10*u.arcmin)
sources = sources[sources["Tmag"] < 12]
print(f"Number of sources: {len(sources)}")
print(sources)

# %%
'''
### Overlaying the sources on a single cutout image

We will get the WCS infomation associated with our cutout so that we can make a WCS-aware plot, and identify a single cutout image to show. Then we display the image and sources together, and label the sources with their row number in the catalog table.
'''

# %%
cutout_wcs = WCS(cutout_hdu[2].header)
cutout_img = cutout_table["FLUX"][start]

# %%
fig, ax = plt.subplots(subplot_kw={'projection':cutout_wcs})
fig.set_size_inches(10,10)
plt.grid(color='white', ls='solid')
    
# Setup WCS axes.
xcoords = ax.coords[0]
ycoords = ax.coords[1]
xcoords.set_major_formatter('d.ddd')
ycoords.set_major_formatter('d.ddd')
xcoords.set_axislabel("RA (deg)")
ycoords.set_axislabel("Dec (deg)")
ax.imshow(cutout_img, cmap=plt.cm.YlGnBu_r,vmin=0,vmax=700)
ax.plot(sources['ra'],sources['dec'],'x',transform=ax.get_transform('icrs'),color="red")

# Annotating the sources with their row nnumber in the sources table
for i,star in enumerate(sources):
    ax.text(star['ra']+0.01,star['dec'],i,transform=ax.get_transform('icrs'))

ax.set_xlim(0,cutout_img.shape[1]-1)
ax.set_ylim(cutout_img.shape[0]-1,0)

plt.show()

# %%
'''
The variable star is row 4 in the catalog sources table.
'''

# %%
sources["ID","ra","dec"][4]

# %%
'''
### Getting the variable star light curve

Again, we will look specifically for the TASOC light curve(s) associated with this star, rather than the mission pipeline ones. Below we go through the same process to serch for the observation, then find the associated data products, and download them.
'''

# %%
variable_tic_id = sources["ID"][4]

variable_res = Observations.query_criteria(target_name=variable_tic_id, 
                                        obs_collection="HLSP", 
                                        filters="TESS")
print(f"Number of tasoc light curves for {variable_tic_id}: {len(variable_res)}\n")

        
variable_prod = Observations.get_product_list(variable_res[0])
variable_manifest = Observations.download_products(variable_prod)

# %%
'''
Note that this time there is only one TASOC light curve, and it is at the 30 minute cadence.  This was not a star that TESS observed at the short cadence.
'''

# %%
hdu = fits.open(variable_manifest["Local Path"][0])
variable_lc = Table(hdu[1].data)
hdu.close()

# %%
'''
### Plotting the variable star light curve

We wil again plot the light curve using bokeh, for the interactive tools.
'''

# %%
bfig = plotting.figure(plot_width=850, plot_height=300, title=f"Detrended Lightcurve (TIC{variable_tic_id})")

bfig.circle(variable_lc["TIME"],variable_lc["FLUX_RAW"], fill_color="black",size=4, line_color=None)
bfig.line(variable_lc["TIME"],variable_lc["FLUX_RAW"], line_color='black')

# Labeling the axes
bfig.xaxis.axis_label = "Time (BTJD)"
bfig.yaxis.axis_label = "Flux"

plotting.show(bfig)

# %%
'''
That looks variable all right!

### Finding the period

We'll run a quick Lomb Scargle priodogram on this light curve to see if we can quantify the periodic behavior. To do this we will use the `astropy.timeseries` class LombScargle (LINK THIS).
'''

# %%
lomb = LombScargle(variable_lc["TIME"], variable_lc["FLUX_RAW"])
frequency, power = lomb.autopower()

# %%
'''
#### Plotting the periodogram
'''

# %%
bfig = plotting.figure(plot_width=850, plot_height=300, x_axis_type="log", x_range=(0.008, 1),
                       title=f"Periodogram (TIC{variable_tic_id})")

bfig.line(1/frequency, power, line_color='black')

# Labeling the axes
bfig.xaxis.axis_label = "Period"
bfig.yaxis.axis_label = "Power"

plotting.show(bfig)

# %%
'''
#### Phasing on the highest power period

We will pick out the highest powered period in the abover periodogram and phase the stellar light curve on that period.
'''

# %%
period = 1/frequency[np.argmax(power)].value
period

# %%
bfig = plotting.figure(plot_width=850, plot_height=300, title=f"Phased Lightcurve (TIC{variable_tic_id})")

# Plotting the phased light curve
bfig.circle(variable_lc["TIME"]%period,variable_lc["FLUX_RAW"], fill_color="black",size=4, line_color=None)

# Plotting the periodic fit
t_fit = np.linspace(0,period,100)
bfig.line(t_fit, lomb.model(t_fit, 1/period), color='#1b9f00', line_width=2)

# Labeling the axes
bfig.xaxis.axis_label = "Phase (days)"
bfig.yaxis.axis_label = "Flux"

plotting.show(bfig)

# %%
