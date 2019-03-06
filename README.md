# ThreshSeg
This repository contains several MATLAB scripts used to do image segmentation via an MBO-type iterative thresholding method.

Also the manuscript is stored in `manuscript/`.

The code is to accompany the paper:  
Dong Wang, Haohan Li, Xiaoyu Wei, Xiaoping Wang.  
An efficient iterative thresholding method for image segmentation. (2016).  
[arXiv:1608.01431 [cs.CV]](https://arxiv.org/abs/1608.01431).

## Code Usage
 - Run ThreshSeg.m to use the GUI interface (somehow it runs much slower in GUI).
 - For better performance, you may call the library directly. To do that, copy main_template.m to a new file (e.g. main.m) and edit the parameters. In the same directory as main.m, put all input files under ./input. Then execute main.m.

## Examples
Run examples/demo_XXXX.m and read comments therein.

## Using GUI
The initials can be set through GUI with a mouse.

 - For rectangular regions, left click on the image twice to pick a rectangle.

 - For polygonal regions, left click to add vertices, and right click to add the last vertex as well as connect the last
   vertex with the first one.

## Initial File Format

### Rectangles
 - The file contains (n_phases-1) lines. Each line consists of four real numbers
   (xmin, xmax, ymin, ymax), which specifies a rectangular region.
 - The entries should between 0 and 1, otherwise they will be normalized to be so.
 - When the rectangular regions intersect, the intersection parts are assigned to none
   of them, and is left to be the lase phase.

### Polygons
 - The file contains 2*(n_phases-1) lines. Every two lines hold information for a
   polygon, the first of which consists of x-coordinates, and the second of which 
   consists of y-coordinates.
 - Since the polygons may have different numbers of vertices, the coordinates may be
   padded with zeros to make the length of lines even. Note that this is optional. You
   may write the file without zero padding and it also works.
 - The rules for out-of-image points and intersecting polygons are the same as
   rectanglular cases.
