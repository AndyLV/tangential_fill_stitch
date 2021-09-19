from shapely.geometry.polygon import LineString, LinearRing
from shapely.geometry import  Point, MultiPoint, linestring
from shapely.ops import nearest_points, polygonize
from collections import namedtuple
import math
from stitches import LineStringSampling
from stitches import PointTransfer
from depq import DEPQ
from stitches import constants

nearest_neighbor_tuple = namedtuple('nearest_neighbor_tuple', ['nearest_point_parent', 'nearest_point_child', 'projected_distance_parent', 'child_node'])


# Cuts a closed line so that the new closed line starts at the point with "distance" to the beginning of the old line.
def cut(line, distance):
    if distance <= 0.0 or distance >= line.length:
        return [LineString(line)]
    coords = list(line.coords)
    for i, p in enumerate(coords):
        if i > 0 and p == coords[0]:
            pd = line.length
        else:
            pd = line.project(Point(p))
        if pd == distance:
            if coords[0] == coords[-1]:
                return LineString(coords[i:]+coords[1:i+1])
            else:
                return LineString(coords[i:]+coords[:i])
        if pd > distance:
            cp = line.interpolate(distance)
            if coords[0] == coords[-1]:
                return LineString([(cp.x, cp.y)] + coords[i:]+coords[1:i]+[(cp.x, cp.y)])
            else:
                return LineString([(cp.x, cp.y)] + coords[i:]+coords[:i])


def cyclic_distance(val1, val2, line_length):
    absdiff = abs(val1-val2)
    # if val1 < val2:
    #print("Needed cyclic_distance")
    return min(absdiff, line_length-absdiff)



    


#Takes the offsetted curves organized as tree, connects and samples them.
#Strategy: A connection from parent to child is made where both curves come closest together.
#Input:
#-tree: contains the offsetted curves in a hierachical organized data structure.
#-usedoffset: used offset when the offsetted curves were generated
#-stitchdistance: maximum allowed distance between two points after sampling
#-closePoint: defines the beginning point for stitching (stitching starts always from the undisplaced curve)
#-offset_by_half: If true the resulting points are interlaced otherwise not.
#Returnvalues:
#-All offsetted curves connected to one line and sampled with points obeying stitchdistance and offset_by_half
#-Tag (origin) of each point to analyze why a point was placed at this position
def connect_raster_tree_nearest_neighbor(tree, usedoffset, stitchdistance, closePoint, offset_by_half):
    # We cut the current item so that its index 0 is closest to closePoint
    currentCoords = tree.val
    absoffset = abs(usedoffset)
    resultCoords = []
    resultCoords_Origin = []

    startDistance = tree.val.project(closePoint)
    if startDistance > 0:
        currentCoords = cut(currentCoords, startDistance)
        tree.val = currentCoords

        if not tree.transferred_point_priority_deque.is_empty():
            newDEPQ = DEPQ(iterable=None, maxlen=None)
            for item,priority in tree.transferred_point_priority_deque:
                newDEPQ.insert(item, math.fmod(
                    priority-startDistance+currentCoords.length, currentCoords.length))
            tree.transferred_point_priority_deque = newDEPQ
        #print("Gecutted")

    ownCoords = []
    ownCoordsOrigin = []
    stitching_direction = 1

    # This list should contain a tuple of nearest points between the current geometry
    # and the subgeometry, the projected distance along the current geometry,
    # and the belonging subtree node
    nearestPointsList = []
    
    for subnode in tree.children:
        point_parent, point_child = nearest_points(currentCoords, subnode.val)
        #if point_pair[0].distance(point_pair[1]) > 1.5*absoffset:
            #print("WARNING: sub geometry closest point is too far apart!")
        projDistance = currentCoords.project(point_parent)
        nearestPointsList.append(nearest_neighbor_tuple(nearest_point_parent = point_parent, 
                                                        nearest_point_child = point_child,
                                                        projected_distance_parent = projDistance,
                                                        child_node=subnode))
    nearestPointsList.sort(reverse=False, key=lambda tup: tup.projected_distance_parent)

    if nearestPointsList:
        start_distance = min(absoffset*constants.factor_offset_starting_points, nearestPointsList[0].projected_distance_parent)
    else:
        start_distance = absoffset*constants.factor_offset_starting_points

    temp, temp_origin = LineStringSampling.rasterLineString2_Priority2(currentCoords, start_distance,  # We add/subtract absoffset/2 to not sample the same point again (avoid double points for start and end)
                                                            currentCoords.length, stitchdistance, stitching_direction, tree.transferred_point_priority_deque, absoffset)
    assert(len(temp) == len(temp_origin))
    temp_origin[0] = LineStringSampling.PointSource.ENTER_LEAVING_POINT

    # We might add the last point wich might not rastered by the previous method
    lastpoint_distance = currentCoords.project(Point(temp[-1]))
    # We use cyclic distance since the point can be assigned distance 0 if the last point is also the first point
    lastpoint_distance = cyclic_distance(
        lastpoint_distance, currentCoords.length, currentCoords.length)
    assert(len(ownCoords) == len(ownCoordsOrigin))
    delta = min(lastpoint_distance, absoffset)
    if delta > constants.line_lengh_seen_as_one_point:
        temp.append(currentCoords.interpolate(
            currentCoords.length-delta/2).coords[0])
        temp_origin.append(LineStringSampling.PointSource.ENTER_LEAVING_POINT)
        # if (abs(temp[-1][0]-9.8) < 0.2 and abs(temp[-1][1]-141.1) < 0.2):
        #    print("HIIER FOUNDED3")
    else:
        temp_origin[-1] = LineStringSampling.PointSource.ENTER_LEAVING_POINT


    if len(temp) > 4:
        ownCoords.extend(temp[1:-1]) # Take not the first an the last one since they are specific to the entering and leaving the geometry
        ownCoordsOrigin.extend(temp_origin[1:-1])
    else:
        ownCoords.extend(temp)
        ownCoordsOrigin.extend(temp_origin)

    assert(len(ownCoords) == len(ownCoordsOrigin))
    
    tree.val = LineString(ownCoords)
    tree.pointsourcelist = ownCoordsOrigin
    tree.stitching_direction = stitching_direction
    tree.already_rastered = True
    assert(len(ownCoords) == len(ownCoordsOrigin))

    #Next we need to transfer our rastered points to siblings and childs
    #First we start with the siblings:
    PointTransfer.transfer_inner_points_to_surrounding2(
        tree, usedoffset, offset_by_half, stitchdistance, False, False)  # Only used to transfer to possible siblings
    assert(len(ownCoords) == len(ownCoordsOrigin))
    to_transfer_point_list = []
    to_transfer_point_list_origin = []
    for k in range(1, len(temp)-1):
        # if abs(temp[k][0]-5.25) < 0.5 and abs(temp[k][1]-42.9) < 0.5:
        #    print("HIER gefunden!")

        #if(
        #   temp_origin[k] == LineStringSampling.PointSource.NOT_NEEDED or
        #   temp_origin[k] == LineStringSampling.PointSource.ADDITIONAL_TRACKING_POINT_NOT_NEEDED):
        #    continue
        #if (offset_by_half and (temp_origin[(k+1) % len(temp_origin)] == LineStringSampling.PointSource.NOT_NEEDED or
        #                        temp_origin[(k+1) % len(temp_origin)] == LineStringSampling.PointSource.ADDITIONAL_TRACKING_POINT_NOT_NEEDED)):
        #    continue
        if (not offset_by_half and temp_origin[k] == LineStringSampling.PointSource.EDGE_NEEDED):
            continue
        if temp_origin[k] == LineStringSampling.PointSource.ENTER_LEAVING_POINT:
            continue
        #if (offset_by_half and ((temp_origin[k] == LineStringSampling.PointSource.EDGE_NEEDED and temp_origin[(k+1) % len(temp_origin)] != LineStringSampling.PointSource.EDGE_NEEDED) or
        #                        (temp_origin[k] != LineStringSampling.PointSource.EDGE_NEEDED and temp_origin[(k+1) % len(temp_origin)] == LineStringSampling.PointSource.EDGE_NEEDED))):
        #    continue
        to_transfer_point_list.append(Point(temp[k]))
        point_origin = temp_origin[k]
        #if (offset_by_half and point_origin != LineStringSampling.PointSource.EDGE_PREVIOUSLY_SHIFTED):
        #    if(temp_origin[(k+1)%len(temp_origin)]==LineStringSampling.PointSource.EDGE_NEEDED):
        #        point_origin = LineStringSampling.PointSource.EDGE_NEEDED
        #    elif temp_origin[(k+1)%len(temp_origin)]==LineStringSampling.PointSource.EDGE_PREVIOUSLY_SHIFTED:
        #        point_origin = LineStringSampling.PointSource.EDGE_PREVIOUSLY_SHIFTED

            
        to_transfer_point_list_origin.append(point_origin)

    assert(len(ownCoords) == len(ownCoordsOrigin))
    #Next we transfer to the childs
    PointTransfer.transfer_transferred_points_also_to_childs(
        tree, usedoffset, offset_by_half, stitchdistance, to_transfer_point_list, False,to_transfer_point_list_origin)
    assert(len(ownCoords) == len(ownCoordsOrigin))
    

    if not nearestPointsList:
        resultCoords = temp
        resultCoords_Origin = temp_origin
    else:
        # if len(nearestPointsList) > 1:
        #    print("HIERRR!")

        #To create a closed ring
        temp.append(temp[0])
        temp_origin.append(temp_origin[0])


        # temp does not start with currentCoords but has an offset (see call of rasterLineString2_Priority2)
        total_distance = absoffset*constants.factor_offset_starting_points
        current_item_index = 0
        resultCoords = [temp[0]]
        resultCoords_Origin = [LineStringSampling.PointSource.ENTER_LEAVING_POINT]
        for i in range(1, len(temp)):
            next_distance = math.sqrt((temp[i][0]-temp[i-1][0])**2 +
                                      (temp[i][1]-temp[i-1][1])**2)
            while (current_item_index < len(nearestPointsList) and
                    total_distance+next_distance+constants.eps > nearestPointsList[current_item_index].projected_distance_parent):
               # was_inside = True
                item = nearestPointsList[current_item_index]
                temp3, temp3_origin = connect_raster_tree_nearest_neighbor(
                    item.child_node, usedoffset, stitchdistance, item.nearest_point_child, offset_by_half)

                delta = item.nearest_point_parent.distance(Point(temp[i-1]))
                if delta > absoffset*constants.factor_offset_starting_points:
                    resultCoords.append(item.nearest_point_parent.coords[0])
                    resultCoords_Origin.append(LineStringSampling.PointSource.ENTER_LEAVING_POINT)
                # reversing avoids crossing when entering and leaving the child segment
                resultCoords.extend(temp3[::-1])
                resultCoords_Origin.extend(temp3_origin[::-1])

                delta = item.nearest_point_parent.distance(Point(temp[i]))
                if current_item_index < len(nearestPointsList)-1:
                    delta = min(delta, abs(
                        nearestPointsList[current_item_index+1].projected_distance_parent-item.projected_distance_parent))

                if delta > absoffset*constants.factor_offset_starting_points:
                    resultCoords.append(currentCoords.interpolate(
                        item.projected_distance_parent+absoffset*constants.factor_offset_starting_points).coords[0])
                    resultCoords_Origin.append(LineStringSampling.PointSource.ENTER_LEAVING_POINT)

                current_item_index += 1
            if i < len(temp)-1:
                if(Point(resultCoords[-1]).distance(Point(temp[i])) > absoffset*constants.factor_offset_remove_points):
                    resultCoords.append(temp[i])
                    resultCoords_Origin.append(temp_origin[i])

            # Since currentCoords and temp are rastered differently there accumulate errors regarding the current distance.
            # Since a projection of each point in temp would be very time consuming we project only every n-th point which resets the accumulated error every n-th point.
            if i % 20 == 0:
                total_distance = currentCoords.project(Point(temp[i]))
            else:
                total_distance += next_distance

    assert(len(ownCoords) == len(ownCoordsOrigin))
    assert(len(resultCoords) == len(resultCoords_Origin))

    return resultCoords, resultCoords_Origin
    #return remove_dense_points(resultCoords, resultCoords_Origin, stitchdistance,absoffset)

#Takes a line and calculates the nearest distance along this line to enter the next_line
#Input:
#-travel_line: The "parent" line for which the distance should be minimized to enter next_line
#-next_line: contains the next_line which need to be entered
#-thresh: The distance between travel_line and next_line needs to below thresh to be a valid point for entering
#Output:
#-tuple - the tuple structure is: (nearest point in travel_line, nearest point in next_line)
def get_nearest_points_closer_than_thresh(travel_line, next_line,thresh):
    point_list = list(MultiPoint(travel_line.coords)) 

    if point_list[0].distance(next_line) < thresh:
        return nearest_points(point_list[0], next_line)

    for i in range(len(point_list)-1):
        line_segment = LineString([point_list[i], point_list[i+1]])
        result = nearest_points(line_segment,next_line)

        if result[0].distance(result[1])< thresh:
            return result
    line_segment = LineString([point_list[-1], point_list[0]])
    result = nearest_points(line_segment,next_line)

    if result[0].distance(result[1])< thresh:
        return result
    else:
        return None


#Takes a line and calculates the nearest distance along this line to enter the childs in children_list
#The method calculates the distances along the line and along the reversed line to find the best direction
#which minimizes the overall distance for all childs.
#Input:
#-travel_line: The "parent" line for which the distance should be minimized to enter the childs
#-children_list: contains the childs of travel_line which need to be entered
#-threshold: The distance between travel_line and a child needs to below threshold to be a valid point for entering
#-preferred_direction: Put a bias on the desired travel direction along travel_line. If equals zero no bias is applied.
# preferred_direction=1 means we prefer the direction of travel_line; preferred_direction=-1 means we prefer the opposite direction.
#Output:
#-stitching direction for travel_line
#-list of tuples (one tuple per child). The tuple structure is: ((nearest point in travel_line, nearest point in child), distance along travel_line, belonging child)
def create_nearest_points_list(travel_line, children_list, threshold, threshold_hard,preferred_direction=0):
    result_list_in_order = []
    result_list_reversed_order = []

    travel_line_reversed = LinearRing(travel_line.coords[::-1])

    weight_in_order = 0
    weight_reversed_order = 0
    for child in children_list:
        result = get_nearest_points_closer_than_thresh(travel_line, child.val, threshold)
        if result == None: #where holes meet outer borders a distance up to 2*used offset can arise
            result = get_nearest_points_closer_than_thresh(travel_line, child.val, threshold_hard)
            assert(result != None)
        proj = travel_line.project(result[0])
        weight_in_order += proj
        result_list_in_order.append(nearest_neighbor_tuple(nearest_point_parent = result[0],
                                                           nearest_point_child = result[1],
                                                           projected_distance_parent = proj,
                                                           child_node = child))

        result = get_nearest_points_closer_than_thresh(travel_line_reversed, child.val, threshold)
        if result == None: #where holes meet outer borders a distance up to 2*used offset can arise
            result = get_nearest_points_closer_than_thresh(travel_line_reversed, child.val, threshold_hard)
            assert(result != None)
        proj = travel_line_reversed.project(result[0])
        weight_reversed_order += proj
        result_list_reversed_order.append(nearest_neighbor_tuple(nearest_point_parent = result[0],
                                                                nearest_point_child = result[1],
                                                                projected_distance_parent = proj,
                                                                child_node = child))

    if preferred_direction == 1:
        weight_in_order=min(weight_in_order/2, max(0, weight_in_order-10*threshold))
        if weight_in_order == weight_reversed_order:
            return (1, result_list_in_order)
    elif preferred_direction == -1:
        weight_reversed_order=min(weight_reversed_order/2, max(0, weight_reversed_order-10*threshold))
        if weight_in_order == weight_reversed_order:
            return (-1, result_list_reversed_order)


    if weight_in_order < weight_reversed_order:
        return (1, result_list_in_order)
    else:
        return (-1, result_list_reversed_order)


#Takes the offsetted curves organized as tree, connects and samples them.
#Strategy: A connection from parent to child is made as fast as possible to reach the innermost child as fast as possible in order
# to stich afterwards from inner to outer. 
#Input:
#-tree: contains the offsetted curves in a hierachical organized data structure.
#-usedoffset: used offset when the offsetted curves were generated
#-stitchdistance: maximum allowed distance between two points after sampling
#-closePoint: defines the beginning point for stitching (stitching starts always from the undisplaced curve)
#-offset_by_half: If true the resulting points are interlaced otherwise not.
#Returnvalues:
#-All offsetted curves connected to one line and sampled with points obeying stitchdistance and offset_by_half
#-Tag (origin) of each point to analyze why a point was placed at this position
def connect_raster_tree_from_inner_to_outer(tree, usedoffset, stitchdistance, closePoint, offset_by_half):
    # We cut the current item so that its index 0 is closest to closePoint
    currentCoords = tree.val
    absoffset = abs(usedoffset)
    resultCoords = []
    resultCoords_Origin = []

    startDistance = tree.val.project(closePoint)
    if startDistance > 0:
        currentCoords = cut(currentCoords, startDistance)
        tree.val = currentCoords
        # print("After cut")
        # checkTree(tree)
        if not tree.transferred_point_priority_deque.is_empty():
            newDEPQ = DEPQ(iterable=None, maxlen=None)
            for item, priority in tree.transferred_point_priority_deque:
                newDEPQ.insert(item, math.fmod(
                    priority-startDistance+currentCoords.length, currentCoords.length))
            tree.transferred_point_priority_deque = newDEPQ
        #print("Gecutted")

    ownCoords = []
    ownCoordsOrigin = []

    #We try to use always the opposite stitching direction with respect to the parent to avoid crossings when entering and leaving the child
    parent_stitching_direction = -1
    if tree.parent != None:
        parent_stitching_direction = tree.parent.stitching_direction
    stitching_direction, nearestPointsList = create_nearest_points_list(currentCoords, tree.children, 1.5*absoffset,2.05*absoffset,-parent_stitching_direction)
    nearestPointsList.sort(reverse=False, key=lambda tup: tup.projected_distance_parent)


    temp = []
    temp_origin = []
    if nearestPointsList:
        startoffset = min(absoffset*constants.factor_offset_starting_points, nearestPointsList[0].projected_distance_parent)
        endoffset = max(currentCoords.length-absoffset*constants.factor_offset_starting_points, nearestPointsList[-1].projected_distance_parent)
    else:
        startoffset = absoffset*constants.factor_offset_starting_points
        endoffset = currentCoords.length-absoffset*constants.factor_offset_starting_points
    

    if stitching_direction == 1:
        temp, temp_origin = LineStringSampling.rasterLineString2_Priority2(currentCoords, startoffset,  # We add startoffset to not sample the same point again (avoid double points for start and end)
                                                            endoffset, stitchdistance, stitching_direction, tree.transferred_point_priority_deque, absoffset)
    else:
        temp, temp_origin = LineStringSampling.rasterLineString2_Priority2(currentCoords, currentCoords.length-startoffset,  # We subtract startoffset to not sample the same point again (avoid double points for start and end)
                                                            currentCoords.length-endoffset, stitchdistance, stitching_direction, tree.transferred_point_priority_deque, absoffset)
        currentCoords.coords = currentCoords.coords[::-1] 

    # We add the last point wich might not rastered by the previous method
    lastpoint_distance = currentCoords.project(Point(temp[-1]))
    # We use cyclic distance since the point can be assigned distance 0 if the last point is also the first point
    lastpoint_distance = cyclic_distance(
        lastpoint_distance, currentCoords.length, currentCoords.length)
    assert(len(ownCoords) == len(ownCoordsOrigin))
    delta = min(lastpoint_distance, absoffset)
    if delta > constants.line_lengh_seen_as_one_point:
        temp.append(currentCoords.interpolate(
            currentCoords.length-delta/2).coords[0]) 
        temp_origin.append(LineStringSampling.PointSource.ENTER_LEAVING_POINT)
        #if (abs(temp[-1][0]-61.7) < 0.2 and abs(temp[-1][1]-105.1) < 0.2):
        #    print("HIIER FOUNDED3")
    else:
        temp_origin[-1] = LineStringSampling.PointSource.ENTER_LEAVING_POINT

    temp_origin[0] = LineStringSampling.PointSource.ENTER_LEAVING_POINT
    assert(len(temp) == len(temp_origin))

    if len(temp) > 4:
        ownCoords.extend(temp[1:-1]) #Do not take the Enter and Leaving point into account
        ownCoordsOrigin.extend(temp_origin[1:-1])
    else:
        ownCoords.extend(temp)
        ownCoordsOrigin.extend(temp_origin)

    assert(len(ownCoords) == len(ownCoordsOrigin))
    
    #Next we need to transfer our rastered points to siblings and childs
    #First we start with the siblings:
    tree.val = LineString(ownCoords)
    tree.pointsourcelist = ownCoordsOrigin
    tree.stitching_direction = stitching_direction
    tree.already_rastered = True
    assert(len(ownCoords) == len(ownCoordsOrigin))
    PointTransfer.transfer_inner_points_to_surrounding2(
        tree, usedoffset, offset_by_half, stitchdistance, False, False)  # Only used to transfer to possible siblings
    assert(len(ownCoords) == len(ownCoordsOrigin))
    to_transfer_point_list = []
    to_transfer_point_list_origin = []
    for k in range(1, len(temp)-1):
        # if abs(temp[k][0]-5.25) < 0.5 and abs(temp[k][1]-42.9) < 0.5:
        #    print("HIER gefunden!")

        #if(
        #   temp_origin[k] == LineStringSampling.PointSource.NOT_NEEDED or
        #   temp_origin[k] == LineStringSampling.PointSource.ADDITIONAL_TRACKING_POINT_NOT_NEEDED):
        #    continue
        #if (offset_by_half and (temp_origin[(k+1) % len(temp_origin)] == LineStringSampling.PointSource.NOT_NEEDED or
        #                        temp_origin[(k+1) % len(temp_origin)] == LineStringSampling.PointSource.ADDITIONAL_TRACKING_POINT_NOT_NEEDED)):
        #    continue
        if (not offset_by_half and temp_origin[k] == LineStringSampling.PointSource.EDGE_NEEDED):
            continue
        if temp_origin[k] == LineStringSampling.PointSource.ENTER_LEAVING_POINT:
            continue
        #if (offset_by_half and ((temp_origin[k] == LineStringSampling.PointSource.EDGE_NEEDED and temp_origin[(k+1) % len(temp_origin)] != LineStringSampling.PointSource.EDGE_NEEDED) or
        #                        (temp_origin[k] != LineStringSampling.PointSource.EDGE_NEEDED and temp_origin[(k+1) % len(temp_origin)] == LineStringSampling.PointSource.EDGE_NEEDED))):
        #    continue
        to_transfer_point_list.append(Point(temp[k]))
        point_origin = temp_origin[k]
        #if (offset_by_half and point_origin != LineStringSampling.PointSource.EDGE_PREVIOUSLY_SHIFTED):
        #    if(temp_origin[(k+1)%len(temp_origin)]==LineStringSampling.PointSource.EDGE_NEEDED):
        #        point_origin = LineStringSampling.PointSource.EDGE_NEEDED
        #    elif temp_origin[(k+1)%len(temp_origin)]==LineStringSampling.PointSource.EDGE_PREVIOUSLY_SHIFTED:
        #        point_origin = LineStringSampling.PointSource.EDGE_PREVIOUSLY_SHIFTED

            
        to_transfer_point_list_origin.append(point_origin)

    assert(len(ownCoords) == len(ownCoordsOrigin))
    assert(len(to_transfer_point_list) == len(to_transfer_point_list_origin))
    #since the projection is only in ccw direction towards inner we need to use "-usedoffset" for stitching_direction==-1
    #print("1: ", len(to_transfer_point_list), ' - ', len(to_transfer_point_list_origin))
    PointTransfer.transfer_transferred_points_also_to_childs(
            tree, stitching_direction*usedoffset, False, stitchdistance, to_transfer_point_list, 
            overnext_child=offset_by_half,transfer_forbidden_points= offset_by_half, to_transfer_points_origin=to_transfer_point_list_origin)
    #print("1: ", len(to_transfer_point_list), ' - ', len(to_transfer_point_list_origin))
    assert(len(to_transfer_point_list) == len(to_transfer_point_list_origin))
    if offset_by_half:# and tree.parent == None:
        PointTransfer.transfer_transferred_points_also_to_childs(
            tree, stitching_direction*usedoffset, True, stitchdistance, to_transfer_point_list,
            overnext_child=False,transfer_forbidden_points=False, to_transfer_points_origin=to_transfer_point_list_origin)

    assert(len(ownCoords) == len(ownCoordsOrigin))
   
    
    if not nearestPointsList:
        resultCoords = temp
        resultCoords_Origin = temp_origin
    else:
        # if len(nearestPointsList) > 1:
        #    print("HIERRR!")

        #Create a closed ring for the following code
        temp.append(temp[0])
        temp_origin.append(temp_origin[0])

        #print("Nearest Distance Proj: ", nearestPointsList[0][1])
        # temp does not start with currentCoords but has an offset (see call of rasterLineString2_Priority2)
        total_distance = startoffset

        current_item_index = 0
        resultCoords = [temp[0]]
        resultCoords_Origin = [temp_origin[0]]

        for i in range(1, len(temp)):
            next_distance = math.sqrt((temp[i][0]-temp[i-1][0])**2 +
                                      (temp[i][1]-temp[i-1][1])**2)
            while (current_item_index < len(nearestPointsList) and
                    total_distance+next_distance+constants.eps > nearestPointsList[current_item_index].projected_distance_parent):
               # was_inside = True
                item = nearestPointsList[current_item_index]
                temp3, temp3_origin = connect_raster_tree_from_inner_to_outer(
                    item.child_node, usedoffset, stitchdistance, item.nearest_point_child, offset_by_half)

                if(Point(resultCoords[-1]).distance(item.nearest_point_parent) > constants.factor_offset_starting_points*absoffset):    
                    resultCoords.append(item.nearest_point_parent.coords[0])
                    resultCoords_Origin.append(LineStringSampling.PointSource.ENTER_LEAVING_POINT)
                #if (abs(resultCoords[-1][0]-61.7) < 0.2 and abs(resultCoords[-1][1]-105.1) < 0.2):
                #    print("HIIER FOUNDED3")
                
                # reversing avoids crossing when entering and leaving the child segment
                resultCoords.extend(temp3)
                resultCoords_Origin.extend(temp3_origin)

                delta = item.nearest_point_parent.distance(Point(temp[i]))
                if current_item_index < len(nearestPointsList)-1:
                    delta = min(delta, abs(
                        nearestPointsList[current_item_index+1].projected_distance_parent-item.projected_distance_parent))

                if delta > constants.factor_offset_starting_points*absoffset:
                    resultCoords.append(currentCoords.interpolate(
                        item.projected_distance_parent+2*constants.factor_offset_starting_points*absoffset).coords[0]) 
                    resultCoords_Origin.append(LineStringSampling.PointSource.ENTER_LEAVING_POINT)
                #if (abs(resultCoords[-1][0]-61.7) < 0.2 and abs(resultCoords[-1][1]-105.1) < 0.2):
                #    print("HIIER FOUNDED3")

                current_item_index += 1
            if i < len(temp)-1:
                if(Point(resultCoords[-1]).distance(Point(temp[i])) > absoffset*constants.factor_offset_remove_points):
                    resultCoords.append(temp[i])
                    resultCoords_Origin.append(temp_origin[i])

            # Since currentCoords and temp are rastered differently there accumulate errors regarding the current distance.
            # Since a projection of each point in temp would be very time consuming we project only every n-th point which resets the accumulated error every n-th point.
            if i % 20 == 0:
                total_distance = currentCoords.project(Point(temp[i]))
            else:
                total_distance += next_distance

    assert(len(ownCoords) == len(ownCoordsOrigin))
    assert(len(resultCoords) == len(resultCoords_Origin))

    #return remove_dense_points(resultCoords, resultCoords_Origin, stitchdistance,absoffset)
    return resultCoords, resultCoords_Origin
