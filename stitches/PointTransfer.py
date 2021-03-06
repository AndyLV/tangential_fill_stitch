from stitches import constants
from stitches import LineStringSampling
from shapely.geometry import  Point, MultiPoint
from shapely.geometry.polygon import LineString, LinearRing
from collections import namedtuple
from shapely.ops import nearest_points
import math

projected_point_tuple = namedtuple('projected_point_tuple', ['point', 'point_source'])

#Calculated the nearest interserction point of "bisectorline" with the coordinates of child (child.val). 
#It returns the intersection point and its distance along the coordinates of the child or "None, None" if no
#intersection was found.
def calc_transferred_point(bisectorline, child):
    result = bisectorline.intersection(child.val)
    if result.is_empty:
        return None, None
    desired_point = Point()
    if result.geom_type == 'Point':
        desired_point = result
    elif result.geom_type == 'LineString':
        desired_point = Point(result.coords[0])
    else:
        resultlist = list(result)
        desired_point = resultlist[0]
        if len(resultlist) > 1:
             desired_point = nearest_points(result, Point(bisectorline.coords[0]))[0]

    priority = child.val.project(desired_point)
    point = desired_point
    return point, priority


#Takes the current tree item and its rastered points (to_transfer_points) and transfers these points to its parent, siblings and childs
# To do so it calculates the current normal and determines its intersection with the neighbors which gives the transferred points.
#Input:
#-treenode: Tree node whose points stored in "to_transfer_points" shall be transferred to its neighbors.
#-used_offset: The used offset when the curves where offsetted
#-offset_by_half: True if the transferred points shall be interlaced with respect to the points in "to_transfer_points"
#-max_stitching_distance: The maximum allowed stitch distance between two points
#-to_transfer_points: List of points belonging to treenode which shall be transferred - it is assumed that to_transfer_points can be handled as closed ring
#-to_transfer_points_origin: The origin tag of each point in to_transfer_points
#-overnext_neighbor: Transfer the points to the overnext neighbor (gives a more stable interlacing)
#-transfer_forbidden_points: Only allowed for interlacing (offset_by_half): Might be used to transfer points unshifted as forbidden points to the neighbor to avoid a point placing there
#-transfer_to_parent: If True, points will be transferred to the parent
#-transfer_to_sibling: If True, points will be transferred to the siblings
#-transfer_to_child: If True, points will be transferred to the childs
#Output:
#-Fills the attribute "transferred_point_priority_deque" of the siblings and parent in the tree datastructure. An item of the deque
#is setup as follows: ((projected point on line, LineStringSampling.PointSource), priority=distance along line)
#index of point_origin is the index of the point in the neighboring line
def transfer_points_to_surrounding(treenode, used_offset, offset_by_half, max_stitching_distance, to_transfer_points,  to_transfer_points_origin=[], 
                                   overnext_neighbor = False, transfer_forbidden_points = False, transfer_to_parent=True, transfer_to_sibling=True, transfer_to_child=True):

    assert(len(to_transfer_points)==len(to_transfer_points_origin) or len(to_transfer_points_origin) == 0)
    assert((overnext_neighbor and not offset_by_half) or not overnext_neighbor)
    assert(not transfer_forbidden_points or transfer_forbidden_points and (offset_by_half or not offset_by_half and overnext_neighbor))

    if len(to_transfer_points) == 0:
        return

    # Get a list of all possible adjacent nodes which will be considered for transferring the points of treenode:
    childs_tuple = treenode.children
    siblings_tuple = treenode.siblings
    # Take only neighbors which have not rastered before
    # We need to distinguish between childs (project towards inner) and parent/siblings (project towards outer)
    child_list = []
    child_list_forbidden = []
    neighbor_list = []
    neighbor_list_forbidden = []

    if transfer_to_child:
        for child in childs_tuple:
            if child.already_rastered == False:
                if not overnext_neighbor:
                    child_list.append(child)
                if transfer_forbidden_points:
                    child_list_forbidden.append(child)
            if overnext_neighbor:
                for subchild in child.children:
                    if subchild.already_rastered == False:
                        child_list.append(subchild)

    if transfer_to_sibling:
        for sibling in siblings_tuple:
            if sibling.already_rastered == False:
                if not overnext_neighbor:
                    neighbor_list.append(sibling)
                if transfer_forbidden_points:
                    neighbor_list_forbidden.append(sibling)
            if overnext_neighbor:
                for subchild in sibling.children:
                    if subchild.already_rastered == False:
                        neighbor_list.append(subchild)

    if transfer_to_parent and treenode.parent != None:
        if treenode.parent.already_rastered == False:
            if not overnext_neighbor:
                    neighbor_list.append(treenode.parent)
            if transfer_forbidden_points:
                    neighbor_list_forbidden.append(treenode.parent)
        if overnext_neighbor:
            if treenode.parent.parent != None:
                if treenode.parent.parent.already_rastered == False:
                    neighbor_list.append(treenode.parent.parent)

    if not neighbor_list and not child_list:
        return

    # Go through all rastered points of treenode and check where they should be transferred to its neighbar
    # The source of each point is stored in treenode.pointsourcelist
    point_list = list(MultiPoint(to_transfer_points))
    point_source_list = to_transfer_points_origin.copy()

    # For a linear ring the last point is the same as the starting point which we delete
    # since we do not want to transfer the starting and end point twice
    closedLine = LineString(to_transfer_points)
    if point_list[0].distance(point_list[-1]) < constants.point_spacing_to_be_considered_equal:
        point_list.pop()
        if(point_source_list):
            point_source_list.pop()
        if len(point_list) == 0:
            return
    else:
        # closed line is needed if we offset by half since we need to determine the line
        # length including the closing segment
        closedLine = LinearRing(to_transfer_points)

    bisectorline_length = abs(used_offset) * \
        constants.transfer_point_distance_factor*(2.0 if overnext_neighbor else 1.0)

    bisectorline_length_forbidden_points = abs(used_offset) * \
        constants.transfer_point_distance_factor

    linesign_child = math.copysign(1, used_offset)


    i = 0
    currentDistance = 0
    while i < len(point_list):
        assert(point_source_list[i] != LineStringSampling.PointSource.ENTER_LEAVING_POINT)
        #if abs(point_list[i].coords[0][0]-47) < 0.3 and abs(point_list[i].coords[0][1]-4.5) < 0.3:
        #    print("HIIIIIIIIIIIERRR")

        # We create a bisecting line through the current point
        normalized_vector_prev_x = (
            point_list[i].coords[0][0]-point_list[i-1].coords[0][0])  # makes use of closed shape
        normalized_vector_prev_y = (
            point_list[i].coords[0][1]-point_list[i-1].coords[0][1])
        prev_spacing = math.sqrt(normalized_vector_prev_x*normalized_vector_prev_x +
                                 normalized_vector_prev_y*normalized_vector_prev_y)

        normalized_vector_prev_x /= prev_spacing
        normalized_vector_prev_y /= prev_spacing


        normalized_vector_next_x = normalized_vector_next_y = 0
        next_spacing = 0
        while True:
            normalized_vector_next_x = (
                point_list[i].coords[0][0]-point_list[(i+1) % len(point_list)].coords[0][0])
            normalized_vector_next_y = (
                point_list[i].coords[0][1]-point_list[(i+1) % len(point_list)].coords[0][1])
            next_spacing = math.sqrt(normalized_vector_next_x*normalized_vector_next_x +
                                     normalized_vector_next_y*normalized_vector_next_y)
            if next_spacing < constants.line_lengh_seen_as_one_point:
                point_list.pop(i)
                if(point_source_list):
                    point_source_list.pop(i)
                currentDistance += next_spacing
                continue

            normalized_vector_next_x /= next_spacing
            normalized_vector_next_y /= next_spacing
            break

        vecx = (normalized_vector_next_x+normalized_vector_prev_x)
        vecy = (normalized_vector_next_y+normalized_vector_prev_y)
        vec_length = math.sqrt(vecx*vecx+vecy*vecy)

        vecx_forbidden_point = vecx
        vecy_forbidden_point = vecy

        # The two sides are (anti)parallel - construct normal vector (bisector) manually:
        # If we offset by half we are offseting normal to the next segment
        if(vec_length < constants.line_lengh_seen_as_one_point or offset_by_half):
            vecx = linesign_child*bisectorline_length*normalized_vector_next_y
            vecy = -linesign_child*bisectorline_length*normalized_vector_next_x

            if transfer_forbidden_points:
                vecx_forbidden_point = linesign_child*bisectorline_length_forbidden_points*normalized_vector_next_y
                vecy_forbidden_point = -linesign_child*bisectorline_length_forbidden_points*normalized_vector_next_x

        else:
            vecx *= bisectorline_length/vec_length
            vecy *= bisectorline_length/vec_length
            
            if (vecx*normalized_vector_next_y-vecy * normalized_vector_next_x)*linesign_child < 0:
                vecx = -vecx
                vecy = -vecy
            vecx_forbidden_point = vecx
            vecy_forbidden_point = vecy

        assert((vecx*normalized_vector_next_y-vecy *
               normalized_vector_next_x)*linesign_child >= 0)

        originPoint = point_list[i]
        originPoint_forbidden_point = point_list[i]
        if(offset_by_half):
            off = currentDistance+next_spacing/2
            if off > closedLine.length:
                off -= closedLine.length
            originPoint = closedLine.interpolate(off)

        bisectorline_child = LineString([(originPoint.coords[0][0],
                                  originPoint.coords[0][1]),
                                   (originPoint.coords[0][0]+vecx,
                                  originPoint.coords[0][1]+vecy)])

        bisectorline_neighbor = LineString([(originPoint.coords[0][0],
                                  originPoint.coords[0][1]),
                                   (originPoint.coords[0][0]-vecx,
                                  originPoint.coords[0][1]-vecy)])

        bisectorline_forbidden_point_child = LineString([(originPoint_forbidden_point.coords[0][0],
                                        originPoint_forbidden_point.coords[0][1]),
                                        (originPoint_forbidden_point.coords[0][0]+vecx_forbidden_point,
                                        originPoint_forbidden_point.coords[0][1]+vecy_forbidden_point)])

        bisectorline_forbidden_point_neighbor = LineString([(originPoint_forbidden_point.coords[0][0],
                                        originPoint_forbidden_point.coords[0][1]),
                                        (originPoint_forbidden_point.coords[0][0]-vecx_forbidden_point,
                                        originPoint_forbidden_point.coords[0][1]-vecy_forbidden_point)])

        for child in child_list:
            point, priority = calc_transferred_point(bisectorline_child,child)
            if point==None:
                continue
            child.transferred_point_priority_deque.insert(projected_point_tuple(point = point, point_source=LineStringSampling.PointSource.OVERNEXT if overnext_neighbor else LineStringSampling.PointSource.DIRECT), priority)
        for child in child_list_forbidden:
            point, priority = calc_transferred_point(bisectorline_forbidden_point_child,child)
            if point == None:
                continue
            child.transferred_point_priority_deque.insert(projected_point_tuple(point=point, point_source=LineStringSampling.PointSource.FORBIDDEN_POINT), priority)
        
        for neighbor in neighbor_list:
            point, priority = calc_transferred_point(bisectorline_neighbor,neighbor)
            if point==None:
                continue
            neighbor.transferred_point_priority_deque.insert(projected_point_tuple(point = point, point_source=LineStringSampling.PointSource.OVERNEXT if overnext_neighbor else LineStringSampling.PointSource.DIRECT), priority)
        for neighbor in neighbor_list_forbidden:
            point, priority = calc_transferred_point(bisectorline_forbidden_point_neighbor,neighbor)
            if point == None:
                continue
            neighbor.transferred_point_priority_deque.insert(projected_point_tuple(point=point, point_source=LineStringSampling.PointSource.FORBIDDEN_POINT), priority)

        i += 1
        currentDistance += next_spacing

    assert(len(point_list) == len(point_source_list))
