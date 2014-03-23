lkopticalflow
=============

Implementation is based on the following paper (section 2.4):
http://robots.stanford.edu/cs223b04/algo_tracking.pdf

Usage:
V = lk(I,J,U,lkwin,levels);

Inputs:
I: Previous Image. It should be a grayscale image with pixel values scaled to [0,1] in double format.
J: Next Image. It should be a grayscale image with pixel values scaled to [0,1] in double format.
U: 2xN array containing the [col,row]' indices or [X,Y]' coordinates of the points of interest in image I. Origin at top left. N is the number of points.
lkwin: window size for lk algorithm. recommended value : 11.
levels: number of levels to use for lk algorithm. recommended value : 2 to 5.

Outputs:
V: 2xN array containing the output points in image J. Same format as V.
