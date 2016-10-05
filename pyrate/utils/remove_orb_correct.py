'''
* removes the orbital error from orbital corrected ifgs to get the original ifgs
* USAGE: remove_orb_correct.py ifglist error_mat
'''

# firstly just try show original - orbital corrected to see if it is working correcly

from PIL import Image
import sys
import numpy as np

import matplotlib.pyplot as pp
from matplotlib import figure

tif1 = Image.open(fp=sys.argv[1])
tif2 = Image.open(fp=sys.argv[2])

tif1_dat = np.array(tif1.getdata()).reshape(tif1.size[::-1])
tif2_dat = np.array(tif2.getdata()).reshape(tif2.size[::-1])

orb_error = np.add(tif1_dat, -tif2_dat)
#orb_error = tif1_dat

plt = pp.imshow(orb_error)
#pp.title('my plot')
#pp.title('my plot')
fig = pp.gcf()
fig.canvas.set_window_title(sys.argv[1])
pp.axis('off')
pp.colorbar()
pp.show()