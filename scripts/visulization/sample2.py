
#http://stackoverflow.com/questions/4860417/placing-custom-images-in-a-plot-window-as-custom-data-markers-or-to-annotate-t

import matplotlib.pyplot as PLT
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from matplotlib._png import read_png

fig = PLT.gcf()
fig.clf()
ax = PLT.subplot(111)

# add a first image
arr_hand = read_png('/auto/proteomics/yingfeng/projects/C100/scripts/visulization/example2.png')
imagebox = OffsetImage(arr_hand, zoom=.5)
xy = [0.25, 0.45]               # coordinates to position this image

ab = AnnotationBbox(imagebox, xy,
    xybox=(30., -30.),
    xycoords='data',
    boxcoords="offset points")                                  
ax.add_artist(ab)

# add second image
arr_vic = read_png('/auto/proteomics/yingfeng/projects/C100/scripts/visulization/b.png')
imagebox = OffsetImage(arr_vic, zoom=.1)
xy = [.6, .3]                  # coordinates to position 2nd image

ab = AnnotationBbox(imagebox, xy,
    xybox=(30, -30),
    xycoords='data',
    boxcoords="offset points")
ax.add_artist(ab)

# rest is just standard matplotlib boilerplate
ax.grid(True)
PLT.draw()
PLT.show()
