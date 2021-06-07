import pandas as pd
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

df=pd.read_csv('reffs_to_plot.csv')

fig=plt.figure(figsize=(18,9))
plt.axhline(1, color='grey', lw=2, ls= '-.')
plt.axvline(0.514, color = 'grey', lw=2, ls ='-.')
plt.plot(df[['weighted_npi']], df[['R_perc']], label = "Model with percolation ", color = 'black', lw = 4)
plt.plot(df[['weighted_npi']], df[['R_noperc']], label ="Standard model", color = 'black', ls = 'dashed', lw = 4)
plt.xlabel('$W_{NPI}$', fontsize = 44)
plt.ylabel('$R_0$', fontsize = 44)
plt.xticks([0,0.25, 0.514, 0.75, 1], ('$0$', '$0.25$', '$T_{perc}$', '0.75', '1'), fontsize = 32)
plt.yticks(fontsize = 32)
plt.legend(loc = 'best', fontsize=44)
#plt.show()
plt.tight_layout()
plt.savefig('r0_wnpi.pdf')

