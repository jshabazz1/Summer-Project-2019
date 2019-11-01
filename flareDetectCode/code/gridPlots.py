import matplotlib.pyplot as plt



def flare_focus (timestart,fluxstart, timestop, fluxstop):

    try:
        fig, axs = plt.subplots(2, 2)
        axs[0, 0].scatter(timestart[0], fluxstart[0])
        axs[0, 0].scatter(timestop[0], fluxstop[0])
        axs[0, 0].set_title('Flare candidate 1')
        axs[0, 1].plot(timestart[1], fluxstart[1], 'tab:orange')
        axs[0, 1].plot(timestop[1], fluxstop[1], 'tab:orange')
        axs[0, 1].set_title('Flare candidate 2')
        axs[1, 0].plot(timestart[2], fluxstart[2], 'tab:green')
        axs[1, 0].plot(timestop[2], fluxstop[2], 'tab:green')
        axs[1, 0].set_title('Flare candidate 3')
        axs[1, 1].plot(timestart[3], fluxstart[3], 'tab:red')
        axs[1, 1].plot(timestop[3], fluxstop[3], 'tab:red')
        axs[1, 1].set_title('Flare candidate 4')
    except:
        pass