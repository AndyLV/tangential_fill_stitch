import svgwrite
import svgtoshapes as svgtoshapes
from stitches import StitchPattern
from shapely.geometry.polygon import LinearRing
from shapely.geometry import Polygon,JOIN_STYLE, Point
import matplotlib.pyplot as plt
import numpy as np

def drawresult(resultcoords, resultcoords_Origin, colorString):
    fig, axs = plt.subplots(1, 1)
    axs.axis('equal')
    plt.gca().invert_yaxis()
    plt.plot(*zip(*resultcoords), colorString)

    colormap = np.array(['r', 'g', 'b', 'c', 'm', 'y', 'k', 'gray', 'm','r'])
    labelmap = np.array(['MUST_USE', 'REGULAR_SPACING', 'INITIAL_RASTERING', 'EDGE_NEEDED', 'NOT_NEEDED',
                        'ALREADY_TRANSFERRED', 'ADDITIONAL_TRACKING_POINT_NOT_NEEDED', 'EDGE_RASTERING_ALLOWED', 'EDGE_PREVIOUSLY_SHIFTED', 'ENTER_LEAVING_POINT'])

    for i in range(0, 9+1):
        # if i != Sampler.PointSource.EDGE_NEEDED and i != Sampler.PointSource.INITIAL_RASTERING:
        #    continue
        selection = []
        for j in range(len(resultcoords)):
            if i == resultcoords_Origin[j]:
                selection.append(resultcoords[j])
        if len(selection) > 0:
            plt.scatter(*zip(*selection), c=colormap[i], label=labelmap[i])

  #  plt.scatter(*zip(*resultcoords),
  #              c=colormap[resultcoords_Origin])
    axs.legend()
    plt.show(block=True)



if __name__ == "__main__":
    mmtopx = 3.77376            #factor from mm to pixel
    offset = 0.125  # mm         #used offset
    stitchlength = 3.0  # mm    #used maximum stitch distance
    filename_input = './Debug2.svg' #Stitch a svg file (note that the svg file must not contain splines - only straight line segments)
    #If filename_input is empty a rectangle with a hole will be stitched as demo
    filename_output = './output.svg'
    joint_style = JOIN_STYLE.mitre #https://shapely.readthedocs.io/en/stable/_images/parallel_offset.png
    interlaced = True
    strategy = StitchPattern.StitchingStrategy.INNER_TO_OUTER

    mypoly = Polygon()
    if filename_input=='':
        #Stitch a rect with a hole:
        outerrectsize = 40  # mm
        innerrectsize = 20  # mm
        innerdistance = (outerrectsize-innerrectsize)/4.0+0.1
        outerring = LinearRing([(0, 0), (0, outerrectsize*mmtopx),
                            (outerrectsize*mmtopx, outerrectsize*mmtopx), (outerrectsize*mmtopx, 0)])
        innerring = LinearRing([(innerdistance*mmtopx, innerdistance*mmtopx), ((outerrectsize-innerdistance)*mmtopx, innerdistance*mmtopx), ((
            outerrectsize-innerdistance)*mmtopx, (outerrectsize-innerdistance)*mmtopx), (innerdistance*mmtopx, (outerrectsize-innerdistance)*mmtopx)])

        mypoly = Polygon(outerring, holes=[innerring])
    else:
        result = list(svgtoshapes.svgtoshapes(filename_input))
        if len(result) == 1:
            mypoly = Polygon(result[0])
        else:
            mypoly = Polygon(result[0], [result[1]])
    
 #   fig, axs = plt.subplots(1, 1)
 #   axs.axis('equal')
 #   plt.gca().invert_yaxis()
 #   x2, y2 = mypoly.exterior.coords.xy
 #   plt.plot(x2, y2, 'r')
 #   plt.show(block=True)

    
    connectedLine, connectedLineOrigin = StitchPattern.offsetPoly(
        mypoly, -offset*mmtopx, joint_style, stitchlength*mmtopx, interlaced,strategy, Point(0,0))

    dwg = svgwrite.Drawing('output.svg', profile='full')
    dwg.add(dwg.polyline(connectedLine,
                         fill="none",
                         stroke='#000000',
                         stroke_width=0.1,
                         stroke_linejoin="round",
                         stroke_linecap="round"))
    dwg.save()

    drawresult(connectedLine, connectedLineOrigin, 'r-')


