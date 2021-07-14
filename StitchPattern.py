from shapely.geometry.polygon import LinearRing, LineString
from shapely.geometry import Polygon, MultiLineString, Point
from shapely.ops import  polygonize
from shapely.geometry import MultiPolygon
from anytree import AnyNode, PreOrderIter
from shapely.geometry.polygon import orient
import ConnectAndSamplePattern as ConnectAndSample
from depq import DEPQ
import constants
from enum import IntEnum


# Problem: When shapely offsets a LinearRing the start/end point might be handled wrongly since they are only treated as LineString.
# (See e.g. https://i.stack.imgur.com/vVh56.png as a problematic example)
# This method checks first whether the start/end point form a problematic edge with respect to the offset side. If it is not a problematic
# edge we can use the normal offset_routine. Otherwise we need to perform two offsets:
# -offset the ring
# -offset the start/end point + its two neighbors left and right
# Finally both offsets are merged together to get the correct offset of a LinearRing
def offsetLinearring(ring, offset, side, resolution, join_style, mitre_limit):
    coords = ring.coords[:]
    # check whether edge at index 0 is concave or convex. Only for concave edges we need to spend additional effort
    dx_seg1 = dy_seg1 = 0
    if coords[0] != coords[-1]:
        dx_seg1 = coords[0][0]-coords[-1][0]
        dy_seg1 = coords[0][1]-coords[-1][1]
    else:
        dx_seg1 = coords[0][0]-coords[-2][0]
        dy_seg1 = coords[0][1]-coords[-2][1]
    dx_seg2 = coords[1][0]-coords[0][0]
    dy_seg2 = coords[1][1]-coords[0][1]
    # use cross product:
    crossvalue = dx_seg1*dy_seg2-dy_seg1*dx_seg2
    sidesign = 1
    if side == 'left':
        sidesign = -1

    # We do not need to take care of the joint n-0 since we offset along a concave edge:
    if sidesign*offset*crossvalue <= 0:
        return ring.parallel_offset(offset, side, resolution, join_style, mitre_limit)

    # We offset along a convex edge so we offset the joint n-0 separately:
    if coords[0] != coords[-1]:
        coords.append(coords[0])
    offsetring1 = ring.parallel_offset(
        offset, side, resolution, join_style, mitre_limit)
    offsetring2 = LineString((coords[-2], coords[0], coords[1])).parallel_offset(
        offset, side, resolution, join_style, mitre_limit)

    # Next we need to merge the results:
    if offsetring1.geom_type == 'LineString':
        return LinearRing(offsetring2.coords[:]+offsetring1.coords[1:-1])
    else:
        # We have more than one resulting LineString for offset of the geometry (ring) = offsetring1.
        # Hence we need to find the LineString which belongs to the offset of element 0 in coords =offsetring2
        # in order to add offsetring2 geometry to it:
        resultList = []
        thresh = constants.offset_factor_for_adjacent_geometry*abs(offset)
        for offsets in offsetring1:
            if(abs(offsets.coords[0][0]-coords[0][0]) < thresh and abs(offsets.coords[0][1]-coords[0][1]) < thresh):
                resultList.append(LinearRing(
                    offsetring2.coords[:]+offsets.coords[1:-1]))
            else:
                resultList.append(LinearRing(offsets))
        return MultiLineString(resultList)


# Removes all geometries which do not form a "valid" LinearRing (meaning a ring which does not form a straight line)
def takeonlyvalidLinearRings(rings):
    if(rings.geom_type == 'MultiLineString'):
        newList = []
        for ring in rings:
            if len(ring.coords) > 3 or (len(ring.coords) == 3 and ring.coords[0] != ring.coords[-1]):
                newList.append(ring)
        if len(newList) == 1:
            return LinearRing(newList[0])
        else:
            return MultiLineString(newList)
    else:
        if len(rings.coords) <= 2:
            return LinearRing()
        elif len(rings.coords) == 3 and rings.coords[0] == rings.coords[-1]:
            return LinearRing()
        else:
            return rings


#Since naturally holes have the opposite point ordering than non-holes we make 
#all lines within the tree "root" uniform (having all the same ordering direction)
def Make_tree_uniform_cw_ccw(root):
    for node in PreOrderIter(root):
        if(node.id == 'hole'):
            node.val.coords = list(node.val.coords)[::-1]


#Used to define which stitching strategy shall be used
class StitchingStrategy(IntEnum):
    CLOSEST_POINT = 0
    INNER_TO_OUTER = 1

# Takes a polygon (which can have holes) as input and creates offsetted versions until the polygon is filled with these smaller offsets.
# These created geometries are afterwards connected to each other and resampled with a maximum stitchdistance.
# The return value is a LineString which should cover the full polygon.
#Input:
#-poly: The shapely polygon which can have holes
#-offset: The used offset for the curves
#-joinstyle: Join style for the offset - can be round, mitered or bevel (https://shapely.readthedocs.io/en/stable/manual.html#shapely.geometry.JOIN_STYLE)
#For examples look at https://shapely.readthedocs.io/en/stable/_images/parallel_offset.png
#-stitchdistance maximum allowed stitch distance between two points
#-offset_by_half: True if the points shall be interlaced
#-strategy: According to StitchingStrategy you can select between different strategies for the connection between parent and childs
#Output:
#-List of point coordinate tuples
#-Tag (origin) of each point to analyze why a point was placed at this position
def offsetPoly(poly, offset, joinstyle, stitchdistance, offset_by_half, strategy):
    orderedPoly = orient(poly, -1)
    orderedPoly = orderedPoly.simplify(
        constants.simplification_threshold, False)
    root = AnyNode(id="node", val=orderedPoly.exterior, already_rastered=False, transferred_point_priority_deque=DEPQ(
        iterable=None, maxlen=None))
    activePolys = [root]
    activeHoles = [[]]

    for holes in orderedPoly.interiors:
        #print("hole: - is ccw: ", LinearRing(holes).is_ccw)
        activeHoles[0].append(
            AnyNode(id="hole", val=holes, already_rastered=False, transferred_point_priority_deque=DEPQ(
                iterable=None, maxlen=None)))

    # counter = 0
    while len(activePolys) > 0:  # and counter < 20:
        # counter += 1
        # print("New iter")
        curPoly = activePolys.pop()
        curHoles = activeHoles.pop()
        polyinners = []

        # outer = curPoly.val.parallel_offset(offset,'left', 5, joinstyle, 10)
        outer = offsetLinearring(curPoly.val, offset, 'left', 5, joinstyle, 10)
        outer = outer.simplify(0.01, False)
        outer = takeonlyvalidLinearRings(outer)

        for j in range(len(curHoles)):
            # inner = closeLinearRing(curHoles[j].val,offset/2.0).parallel_offset(offset,'left', 5, joinstyle, 10)
            inner = offsetLinearring(
                curHoles[j].val, offset, 'left', 5, joinstyle, 10)
            inner = inner.simplify(0.01, False)
            inner = takeonlyvalidLinearRings(inner)
            if not inner.is_empty:
                polyinners.append(Polygon(inner))
        if not outer.is_empty:
            if len(polyinners) == 0:
                if outer.geom_type == 'LineString':
                    result = Polygon(outer)
                else:
                    result = MultiPolygon(polygonize(outer))
            else:
                if outer.geom_type == 'LineString':
                    result = Polygon(outer).difference(
                        MultiPolygon(polyinners))
                else:
                    result = MultiPolygon(outer).difference(
                        MultiPolygon(polyinners))

            if not result.is_empty and result.area > offset*offset/10:
                resultlist = []
                if result.geom_type == 'Polygon':
                    resultlist = [result]
                else:
                    resultlist = list(result)
                # print("New resultlist: ", len(resultlist))
                for pol in resultlist:
                    pol = orient(pol, -1)
                    if pol.area < offset*offset/10:
                        continue
                    # if pol.exterior.is_ccw:
                    #   hole.coords = list(hole.coords)[::-1]
                    polcoords = pol.exterior.simplify(0.01, False)
                    polcoords = takeonlyvalidLinearRings(polcoords)
                    if polcoords.is_empty:
                        continue
                    #print("node: - is ccw: ", LinearRing(polcoords).is_ccw)
                    # if(LinearRing(polcoords).is_ccw):
                    #    print("Fehler!")
                    node = AnyNode(id="node", parent=curPoly,
                                   val=polcoords, already_rastered=False, transferred_point_priority_deque=DEPQ(
                                       iterable=None, maxlen=None))
                    activePolys.append(node)
                    holenodelist = []
                    for hole in pol.interiors:
                        holenode = AnyNode(
                            id="hole", val=hole, already_rastered=False, transferred_point_priority_deque=DEPQ(
                                iterable=None, maxlen=None))
                        for previoushole in curHoles:
                            if Polygon(hole).contains(Polygon(previoushole.val)):
                                previoushole.parent = holenode
                        holenodelist.append(holenode)
                    activeHoles.append(holenodelist)
        for previoushole in curHoles:  # if the previous holes are not contained in the new holes they have been merged with the outer polygon
            if previoushole.parent == None:
                previoushole.parent = curPoly


   # DBG.drawPoly(root, 'r-')

    Make_tree_uniform_cw_ccw(root)
    # print(RenderTree(root))
    if strategy == StitchingStrategy.CLOSEST_POINT:
        connectedLine, connectedLineOrigin = ConnectAndSample.connect_raster_tree_nearest_neighbor(
                root, offset, stitchdistance, Point(0, 0), offset_by_half)
    elif strategy == StitchingStrategy.INNER_TO_OUTER:
        connectedLine, connectedLineOrigin = ConnectAndSample.connect_raster_tree_from_inner_to_outer(
            root, offset, stitchdistance, Point(0, 0), offset_by_half)
    else:
        print("Invalid strategy!")
        assert(0)

    return connectedLine, connectedLineOrigin
<<<<<<< HEAD
=======


if __name__ == "__main__":

    # polygon = Polygon([(0, 0), (1, 1), (1, 0)])
    # test2 = list(polygon.interiors)

    # mmtopx = 3.77376
    # offset = 0.25 #mm
    # stitchlength = 1.5 #mm

    # outerrectsize = 40 #mm
    # innerrectsize = 20 #mm

    # outerring = LinearRing([(0,0), (0,outerrectsize*mmtopx),(outerrectsize*mmtopx,outerrectsize*mmtopx),(outerrectsize*mmtopx,0)])
    # innerring = LinearRing([(innerrectsize*mmtopx,innerrectsize*mmtopx), ((outerrectsize-innerrectsize)*mmtopx,innerrectsize*mmtopx),((outerrectsize-innerrectsize)*mmtopx,(outerrectsize-innerrectsize)*mmtopx),(innerrectsize*mmtopx,74.5)])
    # #innerring = LinearRing([(7,7), (93,7),(93,93),(7,93)])

    # #innerring = LinearRing([(5,5), (95,5),(95,95),(5,95)])

    # #innerring = LinearRing([(69.5,30.5),(80,40),(69.5,69.5),(40,80),(30.5,69.5),(20,40),(30.5,30.5),(40,20)])

    # mypoly = Polygon(outerring,holes=[innerring])

    # #fig, axs = plt.subplots(1, 1)

    # #axs.axis('equal')
    # #plt.gca().invert_yaxis()

    # #plot_MultiPolygon(mypoly,plt,'go-')

    # #myroot = offsetPoly(mypoly,-1,1)
    # #drawPoly(myroot, 'r-')
    # #test = subtractResult(mypoly,myroot,0.6)

    # connectedLine = offsetPoly(mypoly,-offset*mmtopx,1,stitchlength*mmtopx)

    # dwg = svgwrite.Drawing('test.svg', profile='full')
    # dwg.add(dwg.polyline(connectedLine,
    #                             fill="none",
    #                             stroke='rgb(0,0,0)',
    #                             stroke_width=0.5,
    #                             stroke_linejoin="round",
    #                             stroke_linecap="round"))
    # dwg.save()

    mmtopx = 3.77376
    offset = 0.25  # mm
    stitchlength = 3.0  # mm

    outerrectsize = 40  # mm
    innerrectsize = 20  # mm

    innerdistance = (outerrectsize-innerrectsize)/4.0+0.1

    outerring = LinearRing([(0, 0), (0, outerrectsize*mmtopx),
                           (outerrectsize*mmtopx, outerrectsize*mmtopx), (outerrectsize*mmtopx, 0)])
    innerring = LinearRing([(innerdistance*mmtopx, innerdistance*mmtopx), ((outerrectsize-innerdistance)*mmtopx, innerdistance*mmtopx), ((
        outerrectsize-innerdistance)*mmtopx, (outerrectsize-innerdistance)*mmtopx), (innerdistance*mmtopx, (outerrectsize-innerdistance)*mmtopx)])
    #
    # outerring = LinearRing([(0,0), (0,100),(100,100),(100,0)])
    # innerring = LinearRing([(7,7), (93,7),(93,93),(7,93)])

    #mypoly = Polygon(outerring, holes=[innerring])
    filename = "kaninchenscaled.svg"
    result = list(svgtoshapes.svgtoshapes(filename))
    mypoly = Polygon(result[0])
    connectedLine, connectedLineOrigin = offsetPoly(
        mypoly, -offset*mmtopx, 2, stitchlength*mmtopx, True)

    DBG.drawresult(connectedLine, connectedLineOrigin, 'r-')

    dwg = svgwrite.Drawing('test.svg', profile='full')
    dwg.add(dwg.polyline(connectedLine,
                         fill="none",
                         stroke='#000000',
                         stroke_width=0.1,
                         stroke_linejoin="round",
                         stroke_linecap="round"))
    dwg.save()

    # fig, axs = plt.subplots(1, 1)

    # axs.axis('equal')
    # plt.gca().invert_yaxis()

    # plot_MultiPolygon(mypoly,plt,'go-')

    # myroot = offsetPoly(mypoly,-1,1,1)

    # fig, ax = plt.subplots()
    # patch = PolygonPatch(test)
    # ax.add_patch(patch)
    # ax.axis('equal')
    # plt.gca().invert_yaxis()
    # plt.show()

    # for i in range(19):
    #    myresult = offset_polygons(mypoly,-i,1)
    #    plot_MultiPolygon(myresult,print(RenderTree(root, style=AsciiStyle())) plt, 'r-')

    # plt.show()

    # test3 = list(mypoly.interiors)

    # mypolyinner = Polygon(innerring)
    # mypolyouter = Polygon(outerring)

    # test = mypolyouter.difference(mypolyinner)
    # #polygons = [mypolyinner, mypolyouter]
    # #test = unary_union(polygons)
    # test4 = list(test.interiors)

    # filename = "kaninchentest2.svg"
    # result = svgtoshapes.svgtoshapes(filename)

    # fig, axs = plt.subplots(1, 1)
    # axs.axis('equal')
    # plt.gca().invert_yaxis()
    # counter = 0
    # for geom in result:
    #     x,y = geom.coords.xy
    #     plt.plot(x,y,'go-')
    #     geom2 = geom
    #     while not geom2.is_empty and counter < 40:
    #         if geom2.geom_type == 'LineString':
    #             x2,y2 = geom2.coords.xy
    #             plt.plot(x2,y2,'ro-')
    #         else:
    #             for lines in geom2:
    #                 x2,y2 = lines.coords.xy
    #                 plt.plot(x2,y2,'ro-')
    #         geom2 = geom2.parallel_offset(1,'left', 5, 1, 1)
    #         geom2 = geom2.simplify(0.01)
    #         counter = counter+1

    # plt.show()

    # print("Goodbye, World!")
    # #np.absolute()
>>>>>>> Fixed bug for line string equidistant rastering
