import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from Time_and_flux import gettimeflux_1800

time = gettimeflux_1800('167602025', '1')[0]
flux = gettimeflux_1800('167602025', '1')[1]



title = 'SAP_FLUX'
x = np.array(time).byteswap().newbyteorder()
y = np.array(flux).byteswap().newbyteorder()

total = len(x)
flare = pd.DataFrame(y, x)

Writer = animation.writers['imagemagick']
Writer = Writer(fps=60, metadata=dict(artist='Me'), bitrate=1800)

for i in range(50):
	fig, ax = plt.subplots(figsize=(10,6))
	plt.xlim(1618.0, 1619.0)
	plt.ylim(np.min(flare)[0], np.max(flare)[0])
	plt.xlabel('TBJD',fontsize=20)
	plt.ylabel('FLUX',fontsize=20)
	plt.title('TIC 167602025 flare',fontsize=20)



def animate(i):
	ax.scatter(x[:int(i+250)], y[:int(i+250)])

    #p.tick_params(labelsize=17)
    #plt.setp(p.lines,linewidth=7)

ani = matplotlib.animation.FuncAnimation(fig, animate, frames=1000, interval=500, repeat=True)
ani.save('TIC167602025_flare.mp4', writer=Writer)



