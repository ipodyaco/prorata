#!/usr/bin/python
import sys, getopt, warnings, os, re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from matplotlib._png import read_png

X_values = [1, 2, 3, 4, 5]
Y_values = [0.9, 0.23, 1.5, 1.27, 0.4]
plt.bar(X_values, Y_values, 0.1, color='red', edgecolor='none')
plt.xlabel('m/z')
plt.ylabel('intensity')

arr_hand = read_png('/auto/proteomics/yingfeng/projects/C100/scripts/visulization/example3.png')
imagebox = OffsetImage(arr_hand, zoom=1.0)
xy = [3, 1.5]               # coordinates to position this image

ab = AnnotationBbox(imagebox, xy,
    xybox=(0, 30),
    xycoords='data',
    #bbox = dict (boxstyle = "pad =0"),
    boxcoords="offset points")  

ab. patch.set_boxstyle("round")
ax = plt.gca()
                                
ax.add_artist(ab)

plt.axis([0, 5, 0, 1.8])

plt.draw()
plt.show()

plt.savefig("bar_sample", format="png")
