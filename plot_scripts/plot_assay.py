import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

# Specify pCa data file, hill fit file, color and marker type for plotting.
experiments = [ ('log/serca_alone.log', 'fit/serca_alone.fit.txt', 'black', 'o'),
                ('log/serca_afa.log', 'fit/serca_afa.fit.txt', 'blue', 's'),
                ('log/serca_afa_ps16.log', 'fit/serca_afa_ps16.fit.txt', 'red', '^')]

# Plot Normalized activity assays
fig,ax = plt.subplots(figsize=(3,2))

for e in experiments:
    x = np.genfromtxt(e[0],usecols=(0),dtype=float)
    y = np.genfromtxt(e[0],usecols=(8),dtype=float)
    y_err = np.genfromtxt(e[0],usecols=(9),dtype=float)
    x_fit = np.genfromtxt(e[1],usecols=(0),dtype=float)
    y_fit = np.genfromtxt(e[1],usecols=(7),dtype=float)

    # Plot data
    ax.plot(x_fit,y_fit,linewidth=0.5,color=e[2])
    ax.errorbar(x,y,yerr=y_err,linewidth=0.0,ecolor=e[2], elinewidth=0.5, capsize=1)
    ax.scatter(x,y,color=e[2],marker=e[3],s=10)

    # Set axis titles
    ax.set_xlabel('pCa',size=8)
    ax.set_ylabel('Normalized Activity',size=8)
    
    # Set tick params
    ax.set_xlim(8.5,4.5)
    ax.set_ylim(0,1.1)
    ax.xaxis.set_tick_params(labelsize=8)
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.yaxis.set_major_locator(MultipleLocator(0.25))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.yaxis.set_minor_locator(MultipleLocator(0.25))
    
    
plt.tight_layout()
plt.savefig('plot_assay.ps',dpi=300)
#plt.show()
