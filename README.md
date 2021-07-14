# tangential_fill_stitch
Python script which takes a polygon as input and calculates a tangential fill pattern for stitching.

Entry point is Start.py. There you can type in a filename for svg polygon import or -if you leave it empty- a rectangle with hole will be generated for demonstration purposes. The script uses svgtoshapes.py which was taken from https://gist.github.com/un1tz3r0/9f473e4de65787d336ca60681bc6fcbd for importing svg files. You should not use splines but straight lines for your input shape.

As stitching strategy you can use "Stitcher.StitchingStrategy.CLOSEST_POINT" (minimizes the connection length between adjacent offsetted curves) or "Stitcher.StitchingStrategy.INNER_TO_OUTER" (stitches from the starting point as fast as possible to the innermost curve to stitch afterwards curve by curve from inner to outer).

