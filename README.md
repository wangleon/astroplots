Astroplots
==========

A collection of astronomical plots.

Skymap of *Kepler* field of view
---------------------------------

```python
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import keplermap

fig = plt.figure(figsize=(8,8), dpi=150)
kmap = keplermap.KeplerMap(fig, [0.1,0.1,0.8,0.8], season=1)
fig.add_axes(kmap)
kmap.plot_kepler_field(lw=0.5, fc='C0', ec='k', alpha=0.3)
kmap.plot_field_stars()
kmap.plot_clusters()
kmap.grid(color='k', ls='--', alpha=0.2)
kmap.set_xlabel('Right Ascension')
kmap.set_ylabel('Declination')
fig.savefig('kepler_map.png')
plt.show()
```

<img src="https://github.com/wangleon/astroplots/blob/master/keplermap/kepler_map.png" width=650>
