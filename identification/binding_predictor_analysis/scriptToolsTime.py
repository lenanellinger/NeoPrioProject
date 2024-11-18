import matplotlib.pyplot as plt
import os
import numpy as np

dragen = {'MHCflurry': 130.42013883590698, 'MHCnuggetsI': 35.42291522026062, 'NetMHC': 23.87152600288391, 'NetMHCpan': 42.37384915351868, 'NetMHCpanEL': 38.5609495639801}
strelka = {'MHCflurry': 169.65769600868225, 'MHCnuggetsI': 50.71703600883484, 'NetMHC': 35.499252796173096, 'NetMHCpan': 64.75278925895691, 'NetMHCpanEL': 58.31923747062683}

X = np.arange(len(dragen))
ax = plt.subplot(111)
ax.bar(X, dragen.values(), width=0.2, color='b', align='center')
ax.bar(X-0.2, strelka.values(), width=0.2, color='g', align='center')
ax.legend(('DRAGEN','STRELKA'))
plt.xticks(X, dragen.keys())
plt.title("Time in seconds", fontsize=17)
plt.show()

plt.bar(range(len(dragen)), list(dragen.values()), align='center')
plt.xticks(range(len(dragen)), list(dragen.keys()))

plt.savefig(os.path.join("./", 'dragen.png'), bbox_inches='tight', dpi=100)

plt.bar(range(len(strelka)), list(strelka.values()), align='center')
plt.xticks(range(len(strelka)), list(strelka.keys()))

plt.savefig(os.path.join("./images", "time", 'together.png'), bbox_inches='tight', dpi=100)
plt.clf()
mean = {}
for tool in dragen:
    mean[tool] = (dragen[tool] + strelka[tool]) / 2

plt.bar(range(len(mean)), list(mean.values()), align='center')
plt.xticks(range(len(mean)), list(mean.keys()))

plt.savefig(os.path.join("./images", "time", 'mean.png'), bbox_inches='tight', dpi=100)
