#!/usr/bin/env python

### Draft script to display thumbnails of Pyrate tifs and the option to display a 'zoomed in' image when a thumbnail is clicked.
### Creates a set of jpgs (thumbnail and zoomed in versions) that are stored with the original tifs. If a jpg exists, it isn't re-processed.
### Very slow at generating jpgs for some reason. 

## Need to manually load these modules first
    # module load gdal/1.11.1-python
    # module load gmt/5.1.0

## change the img_dir path at the bottom of the script to point to the required directory

import os
import sys
import math
from Tkinter import *
from PIL import Image, ImageFont, ImageDraw
from PIL.ImageTk import PhotoImage
import gdal
import re
import tkFont

# create thumbnail and jpeg directories
def makeThumbs(imgdir, subdir1="thumbnails", subdir2="jpegs"):
    """
    get thumbnail images for all PyRate geotif images in a directory; 
    for each image, create and save a new thumb, or load and return an existing thumb; 
    makes thumb  dir if needed;
    returns a list of (image filename, thumb image object)
    """
    thumb_dir = os.path.join(imgdir, subdir1)
    if not os.path.exists(thumb_dir):
        os.mkdir(thumb_dir)

    jpeg_dir = os.path.join(imgdir, subdir2)
    if not os.path.exists(jpeg_dir):
        os.mkdir(jpeg_dir)

    # get list of tif files and output jpg files
    in_tifs = []
    out_jpgs = []
    out_thumbs = []
    for file in os.listdir(imgdir):
        if file.endswith(".tif"):
            infile = os.path.join(imgdir, file)
            filename = os.path.splitext(file)[0]
            jpg = os.path.join(filename + ".jpg")
            outjpg = os.path.join(jpeg_dir, jpg)
            thumb_jpg = os.path.join(filename + "_thumb.jpg")
            outthumb = os.path.join(thumb_dir, thumb_jpg)
            in_tifs.append(infile)
            out_jpgs.append(outjpg)
            out_thumbs.append(outthumb)

    in_tifs.sort()
    out_jpgs.sort()
    out_thumbs.sort()

    # create thumbnail
    for i, file in enumerate(out_thumbs):
        if os.path.isfile(file):
            print "A jpeg file already exists for %s" % file
        else:
            try:
                # change no data value to -9999
                params = "gdalwarp -srcnodata 0 -dstnodata -9999 %s out1.tif" %(in_tifs[i])
                os.system(params)
                # get min-max data values
                params = "grdinfo out1.tif > temp"
                os.system(params)
                values = open('temp')
                for line in values:
                    line = line.rstrip()
                    if re.search('z_min:', line):
                        line1 = line.strip().split()
                        z_min = float(line1[2])
                        z_max = float(line1[4])
                        max_inc = z_max / 3
                        max_1 = max_inc
                        max_2 = max_inc * 2
                        max_3 = z_max
                        min_inc = z_min / 3
                        min_1 = min_inc
                        min_2 = min_inc * 2
                        min_3 = z_min
                # create colour table with min-max data values
                f = open('col.txt', 'w')
                line1 = "%f 0 0 255" %(max_3)
                f.write(line1 + "\n")
                f.close()
                f = open('col.txt', 'a')
                line2 = "%f 0 170 255" %(max_2)
                f.write(line2 + "\n")    
                line3 = "%f 0 234 255" %(max_1)
                f.write(line3 + "\n") 
                line4 = "0 0 255 0"
                f.write(line4 + "\n") 
                line5 = "%f 255 255 0" %(min_1)
                f.write(line5 + "\n") 
                line6 = "%f 255 170 0" %(min_2)
                f.write(line6 + "\n") 
                line7 = "%f 255 0 0" %(min_3)
                f.write(line7 + "\n") 
                line8 = "-9999 255 255 255"
                f.write(line8 + "\n")
                f.close()
                # apply colour table to tif
                params = "gdaldem color-relief out1.tif col.txt out2.tif"
                os.system(params)
                # use ImageMagick to create smaller jpeg for viewing when thumbnail is clicked
                params = "convert -resize 10% out2.tif out3.tif"
                os.system(params)
                jpg_file = out_jpgs[i]
                params = "convert out3.tif %s" %(jpg_file)
                os.system(params)
                # use ImageMagick to create jpeg thumbnail (200px width)
                params = "convert -thumbnail 200 out2.tif out4.tif"
                os.system(params)
                thumbs_file = out_thumbs[i]
                params = "convert out4.tif %s" %(thumbs_file)
                os.system(params)
                # file cleanup
                os.remove("out1.tif")
                os.remove("temp")
                os.remove("col.txt")
                os.remove("out2.tif")
                os.remove("out3.tif")
                os.remove("out4.tif")
            except Exception, e:
                print e  
    # make list of image filename, thumb image object and jpeg image object
    thumbs = []
    for i, file in enumerate(out_thumbs):
        jpg = out_jpgs[i]
        thumb_jpg = out_thumbs[i]
        thumbobj = Image.open(thumb_jpg)
        jpegobj = Image.open(jpg)
        thumbs.append((thumb_jpg, thumbobj, jpg, jpegobj))
    return thumbs
    return jpeg_dir

# Single pop up window to view individual thumbnail (scrolling around image enabled)
class ScrolledCanvas(Toplevel):
     def __init__(self, jpg, parent=None):
          Toplevel.__init__(self, parent)
          self.master.title(jpg)
          canv = Canvas(self, relief=SUNKEN)
          canv.config(width=1000, height=900)
          canv.config(highlightthickness=0)
          sbarV = Scrollbar(self, orient=VERTICAL)
          sbarH = Scrollbar(self, orient=HORIZONTAL)

          sbarV.config(command=canv.yview)
          sbarH.config(command=canv.xview)

          canv.config(yscrollcommand=sbarV.set)
          canv.config(xscrollcommand=sbarH.set)
          canv.bind('<Button-4>', lambda event: event.widget.yview_scroll(-1, UNITS))
          canv.bind('<Button-5>', lambda event: event.widget.yview_scroll(1, UNITS))

          sbarV.pack(side=RIGHT, fill=Y)
          sbarH.pack(side=BOTTOM, fill=X)

          canv.pack(side=LEFT, expand=YES, fill=BOTH)
          self.im=Image.open(jpg)
          width,height=self.im.size
          canv.config(scrollregion=(0,0,width,height))
          self.im2=PhotoImage(self.im)
          self.imgtag=canv.create_image(0,0,anchor="nw",image=self.im2)
          font = tkFont.Font(family='Helvetica', size=14, weight='bold')
          name = os.path.basename(jpg)
          canv.create_text(150, 30, font=font, text=name)

# View thumbnails
def viewer(imgdir, kind=Toplevel, numcols=None, height=900, width=1050):
    """
    use fixed-size buttons, scrollable canvas;
    sets scrollable (full) size, and places thumbs at absolute x,y 
    coordinates in canvas;  caveat: assumes all thumbs are same size
    """
    win = kind()
    win.title('Viewing: ' + imgdir)
    #quit = Button(win, text='Quit', command=win.quit, bg='white')
    #quit.pack(side=BOTTOM)

    canvas = Canvas(win, borderwidth=0)

    ## Vertical scroll bar
    vbar = Scrollbar(win)
    vbar.pack(side=RIGHT, fill=Y)
    vbar.config(command=canvas.yview) 
    canvas.config(yscrollcommand=vbar.set) 

    ## Horizonal scroll bar
    #hbar = Scrollbar(win, orient='horizontal')
    #hbar.pack(side=BOTTOM, fill=X)          
    #hbar.config(command=canvas.xview)
    #canvas.config(xscrollcommand=hbar.set)        

    canvas.pack(side=LEFT)
    canvas.config(height=height, width=width)       # init viewable area size
                                                    # changes if user resizes
    thumbs = makeThumbs(imgdir)                     # [(thumb_jpg, thumbobj, jpg, jpegobj)]
    numthumbs = len(thumbs)

    numcols = 5
    #if not numcols:
    #    numcols = int(math.ceil(math.sqrt(numthumbs)))  # fixed or N x N
    numrows = int(math.ceil(numthumbs / numcols))       # 3.x true div

    linksize = max(thumbs[0][1].size)                   # (width, height)
    fullsize = (0, 0,                                   # upper left  X,Y
        (linksize * numcols), (linksize * numrows) )    # lower right X,Y
    canvas.config(scrollregion=fullsize)                # scrollable area size

    rowpos = 0
    savephotos = []
    while thumbs:
        thumbsrow, thumbs = thumbs[:numcols], thumbs[numcols:]
        colpos = 0
        for (thumb_jpg, thumbobj, jpg, jpegobj) in thumbsrow:
            photo = PhotoImage(thumbobj)
            jpg1 = os.path.basename(jpg)
            dates = jpg1.split('_')[0]
            font = tkFont.Font(family='Helvetica', size=10, weight='normal')
            link = Button(canvas, text=dates, image=photo, compound="center", font=font, foreground='black')
            handler = lambda savefile=jpg: ScrolledCanvas(savefile)
            link.config(command=handler, width=linksize, height=linksize)
            link.pack(side="top")
            #link.pack(side=LEFT, expand=YES)
            canvas.create_window(colpos, rowpos, anchor=NW,
                    window=link, width=linksize, height=linksize)
            colpos += linksize
            canvas.bind('<Button-4>', lambda event: event.widget.yview_scroll(-1, UNITS))
            canvas.bind('<Button-5>', lambda event: event.widget.yview_scroll(1, UNITS))
            savephotos.append(photo)
        rowpos += linksize
    return win, savephotos

if __name__ == '__main__':
    #imgdir = os.getcwd()
    #imgdir="/g/data1/dg9/INSAR_ANALYSIS/LOWER_DARLING/RSAT2/PYRATE/T210D_out_prepifg/"
    imgdir= "/g/data1/dg9/INSAR_ANALYSIS/ORD/RSAT2/PYRATE/T125D_out_prepifg/"
    main, save = viewer(imgdir, kind=Tk)
    main.mainloop()
