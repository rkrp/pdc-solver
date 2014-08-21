from matplotlib import pyplot as plt
from scipy.signal import argrelextrema
import numpy as np

splot = 111

y = np.array([56.8,41.3,53.5,62.7,52.9,44.8,58.7,60,42,50.3,60.9,56.8,42.2,56.2,62.3,44.5,47.6,59.3,61,41.9,55.6,64.4,46.8,47.3,59.9,61,41.7,55.2,64.7,47.5,46.5,58.5,60.8,41.4,53.6,63.8,48.7,44.3,56.4,63.1,41.8,50.3,63.3,52.6,41.6,53.5,63.9,44.6,47.5,63,58.3,41.4,51.3,64.5,47.9,46.2,62,58.8,41.9,50.8,64.5,49.2,44.8,61.1,59.8,42.1,50.1,64.3,49.7,44.8,60.8,60.1,41.8,48.9,64,50.4,44.8,60.4,60.8,42.1,48.6,63.6,51.9,43.3,58.7,62.6,43.8,47.7,62.9,53.1,42.6,58.5,63,43.8,47.3,62.8,53.3,42.7,58.2,63.1,44,47.1,62.6,53.4,42.5,58,63.2,44.2,46.2,62.1,54,42.6,58,63.3,44.3,45.7,61.7,54.9,42,57,63.8,45.3,45.5,61.1,55.3,42.4,57.1,64.5,45.7,44.8,60.6,56.9,41,55.5,65,47.2,44.2,59.6,58.1,41.5,54.5,65.2,48.3,43.1,58.9,58.9,41.8,53.7,65.1,49.5,42.5,57.7,60.5,42.1,52.2,65.1,51.2,42.3,57.1,61.5,42.9,51.9,65.3,51.9,41.4,55.7,62.3,43.6,50.7,64.7,53.1,41.3,55.2,62.9,44.1,50,64.2,54.1,41.6,54.8,63.3,44.2,49.4,63.6,54.9,41.9,54,63.8,45.1,48.5,63.5,55.8,42,53.2])
x = np.array(range(len(y)))

pos = argrelextrema(np.array(y), np.greater)
pos = pos[0]
posdiff = [pos[n] - pos[n-1] for n in range(2, len(pos))]

print np.mean(posdiff)


#~ plt.subplot(splot)
#~ plt.plot(x, y, 'b-')
#~ plt.xlabel('Time')
#~ plt.ylabel('PV')
#~ plt.show()
