import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

# Name of columns. Manually enter pkcas from .txt files.
names = ['SERCA', '+AFA', '+pSer16-AFA']
pkcas = [6.8413, 6.4668, 6.6277]
erros = [0.0111, 0.0121, 0.0471]

# Plot bar
fig,ax = plt.subplots(figsize=(1.5,2))

# Bar params
index = np.arange(len(names))
bar_width = 0.5
opacity = 0.5

# Plot
ax.bar(index, pkcas, bar_width, alpha=opacity,tick_label=names,color='blue')
ax.errorbar(index,pkcas,yerr=erros,linewidth=0.0,ecolor='black', elinewidth=0.5, capsize=1)

# Set tick params
ax.set_ylim(6.2,7.0)
for tick in ax.get_xticklabels():
    tick.set_rotation(45)
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=8)
    
# Set axis titles
ax.set_ylabel('pKCa',size=8)

plt.tight_layout()
plt.savefig('plot_assay_bar.ps',dpi=300)
plt.show()
