from shapely.geometry.polygon import LineString, LinearRing
from shapely.geometry import  Point, MultiPoint
from shapely.ops import nearest_points, polygonize
import math
import LineStringSampling as Sampler
import PointTransfer as PointTransferrer
from depq import DEPQ
import constants

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


<<<<<<< HEAD
=======
# Receives an AnyTree where each node contains a closed path. The goal is to create one new path which contains all paths within AnyTree as well as their connection
# closePoint is used to calculate the closest point on the outest closed path which is used as starting point.
def connect_raster_tree(tree, usedoffset, stitchdistance, closePoint):
    # We cut the current item so that its index 0 is closest to closePoint
    currentCoords = tree.val
    absoffset = abs(usedoffset)
    resultCoords = []
    # print("Bin HIER!")
    startDistance = tree.val.project(closePoint)
    if startDistance > 0:
        currentCoords = cut(currentCoords, startDistance)

    # This list should contain a tuple of nearest points between the current geometry
    # and the subgeometry, the projected distance along the current geometry,
    # and the belonging subtree node
    nearestPointsList = []
    for subnode in tree.children:
        point_pair = nearest_points(currentCoords, subnode.val)
        if point_pair[0].distance(point_pair[1]) > 1.5*absoffset:
            print("WARNING: sub geometry closest point is too far apart!")
        projDistance = currentCoords.project(point_pair[0])
        nearestPointsList.append((point_pair, projDistance, subnode))

    if len(nearestPointsList) > 0:
        if len(nearestPointsList) > 1:
            print("HIERRR!")
        # Use reverse since pop() removes the last element which should be the nearest
        nearestPointsList.sort(reverse=True, key=lambda tup: tup[1])
        lastDistance = 0
        while nearestPointsList:
            item = nearestPointsList.pop()
            resultCoords.extend(Sampler.rasterLineString2(
                substring(currentCoords, lastDistance, item[1]), stitchdistance))
            resultCoords.extend(connect_raster_tree(
                item[2], usedoffset, stitchdistance, item[0][1]))
            lastDistance = item[1]
        resultCoords.extend(Sampler.rasterLineString2(
            substring(currentCoords, lastDistance, currentCoords.length), stitchdistance))
    else:
        resultCoords = Sampler.rasterLineString2(currentCoords, stitchdistance)

    return resultCoords

# Uses nearest point of the adjacent neighbors for transfer


def transfer_inner_points_to_surrounding(treenode, abs_used_offset, offset_by_half):
    # Get a list of all possible adjacent nodes which will be considered for transferring the points of treenode:
    siblings_tuple = treenode.siblings
    # Take only siblings which have not rastered before
    neighbor_list = []
    for sibling in siblings_tuple:
        if sibling.already_rastered == False:
            neighbor_list.append(sibling)

    # Take also the parend into account:
    if treenode.parent != None:
        neighbor_list.append(treenode.parent)

    # Go through all rastered points of treenode and check whether they should be transferred to its neighbar
    point_list = list(MultiPoint(treenode.val.coords))
    # For a linear ring the last point is the same as the starting point
    if point_list[0] == point_list[-1]:
        point_list.pop()

    for i in range(len(point_list)):
        prev_spacing = point_list[i].distance(
            point_list[i-1])  # makes use of closed shape
        # ahead_spacing = point_list[i].distance(
        #    point_list[i+1])  # makes use of closed shape
        for neighbor in neighbor_list:
            point_pair = nearest_points(neighbor.val, point_list[i])
            distance = point_pair[0].distance(point_pair[1])
            if distance < constants.offset_factor_for_adjacent_geometry*abs_used_offset:
                priority = neighbor.val.project(point_pair[0])
                point = point_pair[0]
                if(offset_by_half and priority > prev_spacing/2):
                    priority -= prev_spacing/2
                    point = neighbor.val.interpolate(priority)
                neighbor.transferred_point_priority_deque.insert(
                    (point, i, id(treenode)), priority)  # We add also i as well as the source id as information which is helpful later to detect interrupted transferred points to the parent (e.g. because another sibling was closer)


# def checkTree(tree):
#    for node in PreOrderIter(tree):
#        if(node.id == 'hole' and not LinearRing(node.val).is_ccw):
#            print("Error in Check!")
#        if(node.id == 'node' and LinearRing(node.val).is_ccw):
#            print("ERROR2 in checktree")


def transfer_transferred_points_also_to_childs(treenode, used_offset, offset_by_half, max_stitching_distance, to_transfer_points, recursive=True):
   # checkTree(treenode)

    if len(to_transfer_points) == 0:
        return

    # Get a list of all possible adjacent nodes which will be considered for transferring the points of treenode:
    childs_tuple = treenode.children
    # Take only siblings which have not rastered before
    child_list = []
    for child in childs_tuple:
        if child.already_rastered == False:
            # child.to_transfer = []
            child_list.append(child)

    to_transfer_points_per_child = [[] for i in range(len(child_list))]

    # Go through all rastered points of treenode and check where they should be transferred to its neighbar
    # The source of each point is stored in treenode.pointsourcelist
    point_list = list(MultiPoint(to_transfer_points))

    # For a linear ring the last point is the same as the starting point which we delete
    # since we do not want to transfer the starting and end point twice
    closedLine = treenode.val
    if point_list[0].distance(point_list[-1]) < constants.point_spacing_to_be_considered_equal:
        point_list.pop()
        if len(point_list) == 0:
            return
    else:
        # closed line is needed if we offset by half since we need to determine the line
        # length including the closing segment
        closedLine = LinearRing(treenode.val)

    bisectorline_length = abs(used_offset) * \
        constants.transfer_point_distance_factor
   # if hasattr(treenode, 'stitching_direction'):
   #     linesign = -treenode.stitching_direction
  #  else:
    linesign = math.copysign(1, used_offset)
   # print(RenderTree(treenode))
   # if treenode.id == 'hole':
   #     linesign *= -1

    i = 0
    while i < len(point_list):
        currentDistance = closedLine.project(point_list[i])

        if abs(point_list[i].coords[0][0]-9.2) < 0.2 and abs(point_list[i].coords[0][1]-150.5) < 0.2:
            print("HIIIIIIIIIIIERRR")

        # We create a bisecting line through the current point
        normalized_vector_prev_x = (
            point_list[i].coords[0][0]-point_list[i-1].coords[0][0])  # makes use of closed shape
        normalized_vector_prev_y = (
            point_list[i].coords[0][1]-point_list[i-1].coords[0][1])
        prev_spacing = math.sqrt(normalized_vector_prev_x*normalized_vector_prev_x +
                                 normalized_vector_prev_y*normalized_vector_prev_y)

        # If spacing is too large we interpolate for a better bisecting line
        if(prev_spacing-constants.line_lengh_seen_as_one_point > max_stitching_distance):
            prev_spacing = min(max_stitching_distance, closedLine.length/10.0)
            prev_point_distance = currentDistance-prev_spacing
            if currentDistance < prev_spacing:
                prev_point_distance += closedLine.length  # makes use of closed shape
            prev_point = closedLine.interpolate(prev_point_distance)

            normalized_vector_prev_x = (
                point_list[i].coords[0][0]-prev_point.coords[0][0])
            normalized_vector_prev_y = (
                point_list[i].coords[0][1]-prev_point.coords[0][1])

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
                continue

            # If spacing is too large we interpolate for a better bisecting line
            if(next_spacing-constants.line_lengh_seen_as_one_point > max_stitching_distance):
                next_spacing = min(max_stitching_distance,
                                   closedLine.length/10.0)
                next_point_distance = currentDistance+next_spacing
                if next_point_distance > closedLine.length:
                    next_point_distance -= closedLine.length  # makes use of closed shape
                next_point = closedLine.interpolate(next_point_distance)

                normalized_vector_next_x = (
                    point_list[i].coords[0][0]-next_point.coords[0][0])
                normalized_vector_next_y = (
                    point_list[i].coords[0][1]-next_point.coords[0][1])

            normalized_vector_next_x /= next_spacing
            normalized_vector_next_y /= next_spacing
            break

        vecx = (normalized_vector_next_x+normalized_vector_prev_x)
        vecy = (normalized_vector_next_y+normalized_vector_prev_y)
        vec_length = math.sqrt(vecx*vecx+vecy*vecy)

        # The two sides are (anti)parallel - construct normal vector (bisector) manually:
        # If we offset by half we are offseting normal to the next segment
        if(vec_length < constants.line_lengh_seen_as_one_point or offset_by_half):
            vecx = linesign*bisectorline_length*normalized_vector_next_y
            vecy = -linesign*bisectorline_length*normalized_vector_next_x
        else:
            vecx *= bisectorline_length/vec_length
            vecy *= bisectorline_length/vec_length
            if (vecx*normalized_vector_next_y-vecy * normalized_vector_next_x)*linesign < 0:
                vecx = -vecx
                vecy = -vecy

        assert((vecx*normalized_vector_next_y-vecy *
               normalized_vector_next_x)*linesign >= 0)

        originPoint = point_list[i]
        if(offset_by_half):
            off = currentDistance+next_spacing/2
            # negative priority values start from end
            # if off < 0:
            #    off += treenode.val.length
            if off > closedLine.length:
                off -= closedLine.length
            originPoint = closedLine.interpolate(off)

        bisectorline = LineString([(originPoint.coords[0][0],
                                  originPoint.coords[0][1]),
                                   (originPoint.coords[0][0]+vecx,
                                  originPoint.coords[0][1]+vecy)])

        for child, transfer_point_list_child in zip(child_list, to_transfer_points_per_child):
            result = bisectorline.intersection(child.val)
            if result.is_empty:
                continue
            desired_point = Point()
            if result.geom_type == 'Point':
                desired_point = result
            else:
                resultlist = list(result)
                desired_point = resultlist[0]
                if len(resultlist) > 1:
                    desired_point = nearest_points(result, point_list[i])[0]

            priority = child.val.project(desired_point)
            point = desired_point

            if abs(point.coords[0][0]-34) < 0.5 and abs(point.coords[0][1]-165.57) < 0.5:
                print("HIER gefunden!")
            if(recursive):
                transfer_point_list_child.append(point)
            child.transferred_point_priority_deque.insert(
                (point, i, id(treenode), Sampler.PointSource.ALREADY_TRANSFERRED), priority)  # We add also i as well as the source id as information which is helpful later to detect interrupted transferred points to the parent (e.g. because another sibling was closer)
        i += 1

    if(recursive):
        for child, transfer_point_list_child in zip(child_list, to_transfer_points_per_child):
            transfer_transferred_points_also_to_childs(
                child, used_offset, offset_by_half, max_stitching_distance, transfer_point_list_child)

# Calculates the current normal and determines its intersection with the neighbors
# to determine the transferred points


def transfer_inner_points_to_surrounding2(treenode, used_offset, offset_by_half, max_stitching_distance, transfer_to_siblings_childs=True):
    # if(treenode.id == 'hole' and not LinearRing(treenode.val).is_ccw):
    #    print("Error in transfer!")
    # print("Hier in Transfer!")
    # checkTree(treenode)

    # Get a list of all possible adjacent nodes which will be considered for transferring the points of treenode:
    siblings_tuple = treenode.siblings
    # Take only siblings which have not rastered before
    neighbor_list = []
    siblings_list = []
    for sibling in siblings_tuple:
        if sibling.already_rastered == False:
            # sibling.to_transfer = []
            neighbor_list.append(sibling)
            siblings_list.append(sibling)

    # to_transfer_points_per_sibling = [[] for i in range(len(siblings_list))]
    to_transfer_points_per_sibling = {
        id(sibling): [] for sibling in siblings_list}

    # Take also the parend into account:
    if treenode.parent != None:
        neighbor_list.append(treenode.parent)

    # if treenode.id == 'hole' and treenode.parent.id == 'node':
    #    print("HIER IB")

    # Go through all rastered points of treenode and check where they should be transferred to its neighbar
    # The source of each point is stored in treenode.pointsourcelist
    point_list = list(MultiPoint(treenode.val.coords))
    point_source_list = treenode.pointsourcelist.copy()
    assert(len(point_list) == len(point_source_list))

    # For a linear ring the last point is the same as the starting point which we delete
    # since we do not want to transfer the starting and end point twice
    closedLine = treenode.val
    if point_list[0].distance(point_list[-1]) < constants.point_spacing_to_be_considered_equal:
        point_list.pop()
        point_source_list.pop()
    else:
        # closed line is needed if we offset by half since we need to determine the line
        # length including the closing segment
        closedLine = LinearRing(treenode.val)

    bisectorline_length = abs(used_offset) * \
        constants.transfer_point_distance_factor
    linesign = treenode.stitching_direction
   # print(RenderTree(treenode))
   # if treenode.id == 'hole':
   #     linesign *= -1

    i = 0
    currentDistance = 0
    while i < len(point_list):
        if abs(point_list[i].coords[0][0]-6.1) < 0.5 and abs(point_list[i].coords[0][1]-42) < 0.5:
            print("HIER gefunden!")
        need_additional_nearest_neighbor = False

        point_source = Sampler.PointSource.MUST_USE

        # We do not transfer points which were only characteristic ("fix") to the geometry of the source
        if(offset_by_half and point_source_list[i] == Sampler.PointSource.EDGE_NEEDED and
                point_source_list[(i+1) % len(point_list)] == Sampler.PointSource.EDGE_NEEDED):
            point_source = Sampler.PointSource.EDGE_RASTERING_ALLOWED
        elif(point_source_list[i] == Sampler.PointSource.EDGE_NEEDED):
            point_source = Sampler.PointSource.NOT_NEEDED
        elif(point_source_list[i] == Sampler.PointSource.ALREADY_TRANSFERRED or
             point_source_list[i] == Sampler.PointSource.NOT_NEEDED or
             point_source_list[i] == Sampler.PointSource.ADDITIONAL_TRACKING_POINT_NOT_NEEDED):
            currentDistance += point_list[i].distance(
                point_list[(i+1) % len(point_list)])
            i += 1
            continue

        # if(abs(point_list[i].coords[0][0]-134) < 0.5 and abs(point_list[i].coords[0][1]-75.5) < 0.5):
        #    print("HIERRR point")
        # if(abs(point_list[i].coords[0][0]-133) < 0.5 and abs(point_list[i].coords[0][1]-75.5) < 0.5):
        #    print("HIERRR point2")
        # prev_spacing = point_list[i].distance(
        #    point_list[i-1])
        # if(abs(point_list[i].coords[0][0]-54) < 0.5 and abs(point_list[i].coords[0][1]-16.5) < 0.5):
        #    print("HIERRR point2")

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
                point_source_list.pop(i)
                currentDistance += next_spacing
                continue
            normalized_vector_next_x /= next_spacing
            normalized_vector_next_y /= next_spacing
            break

        # We do not transfer points which were only characteristic ("fix") to the geometry of the source
        # Since offset_by_half depends also on the following point, we need to check this point, too
        if((offset_by_half and point_source_list[(i+1) % len(point_list)] == Sampler.PointSource.EDGE_NEEDED) and
                point_source_list[i] != Sampler.PointSource.EDGE_NEEDED):
            point_source = Sampler.PointSource.NOT_NEEDED
        elif(offset_by_half and point_source_list[(i+1) % len(point_list)] == Sampler.PointSource.ALREADY_TRANSFERRED or
             offset_by_half and point_source_list[(i+1) % len(point_list)] == Sampler.PointSource.NOT_NEEDED or
             offset_by_half and point_source_list[(i+1) % len(point_list)] == Sampler.PointSource.ADDITIONAL_TRACKING_POINT_NOT_NEEDED):
            i += 1
            currentDistance += next_spacing
            continue

        vecx = (normalized_vector_next_x+normalized_vector_prev_x)
        vecy = (normalized_vector_next_y+normalized_vector_prev_y)
        vec_length = math.sqrt(vecx*vecx+vecy*vecy)

        # The two sides are (anti)parallel - construct normal vector (bisector) manually:
        # If we offset by half we are offseting normal to the next segment
        if(vec_length < 1E-2 or offset_by_half):
            vecx = linesign*bisectorline_length*normalized_vector_next_y
            vecy = -linesign*bisectorline_length*normalized_vector_next_x
            need_additional_nearest_neighbor = offset_by_half
        else:
            vecx *= bisectorline_length/vec_length
            vecy *= bisectorline_length/vec_length
            if (vecx*normalized_vector_next_y-vecy * normalized_vector_next_x)*linesign < 0:
                vecx = -vecx
                vecy = -vecy
            need_additional_nearest_neighbor = True

        assert((vecx*normalized_vector_next_y-vecy *
               normalized_vector_next_x)*linesign >= 0)

        originPoint = point_list[i]
        if(offset_by_half):
            off = currentDistance+next_spacing/2
            # negative priority values start from end
            # if off < 0:
            #    off += treenode.val.length
            if off > closedLine.length:
                off -= closedLine.length
            originPoint = closedLine.interpolate(off)

        bisectorline = LineString([(originPoint.coords[0][0],
                                  originPoint.coords[0][1]),
                                   (originPoint.coords[0][0]+vecx,
                                  originPoint.coords[0][1]+vecy)])

        for neighbor in neighbor_list:
            result = bisectorline.intersection(neighbor.val)
            if not result.is_empty:
                desired_point = Point()
                if result.geom_type == 'Point':
                    desired_point = result
                else:
                    resultlist = list(result)
                    desired_point = resultlist[0]
                    if len(resultlist) > 1:
                        desired_point = nearest_points(
                            result, point_list[i])[0]

                priority = neighbor.val.project(desired_point)
                point = desired_point

                # if abs(point.coords[0][0]-9.4) < 0.2 and abs(point.coords[0][1]-13.3) < 0.2:
                #    print("HIIIIIIIIIIIERRR")

                # if abs(point.coords[0][0]-58) < 0.5 and abs(point.coords[0][1]-132) < 0.5:
                #    print("HIER gefunden!")
                # its a sibling - the transferred point to this sibling needs to be transferred to the siblings childs afterwards
                if transfer_to_siblings_childs and neighbor.parent == treenode.parent and point_source != Sampler.PointSource.NOT_NEEDED:
                    to_transfer_points_per_sibling[id(neighbor)].append(point)
                    # neighbor.to_transfer.append(point)
                neighbor.transferred_point_priority_deque.insert(
                    (point, i, id(treenode), point_source), priority)  # We add also i as well as the source id as information which is helpful later to detect interrupted transferred points to the parent (e.g. because another sibling was closer)

            if need_additional_nearest_neighbor:
                # if abs(point_list[i].coords[0][0]-10.4) < 0.2 and abs(point_list[i].coords[0][1]-14) < 0.2:
                #    print("HIIIIIIIIIIIERRR")
                add_point = nearest_points(neighbor.val, point_list[i])[0]

                # if abs(add_point.coords[0][0]-9.4) < 0.2 and abs(add_point.coords[0][1]-13.3) < 0.2:
                #    print("HIIIIIIIIIIIERRR")

                if add_point.distance(point_list[i]) < constants.offset_factor_for_adjacent_geometry*abs(used_offset):
                    priority = neighbor.val.project(add_point)
                    neighbor.transferred_point_priority_deque.insert(
                        (add_point, i, id(treenode), Sampler.PointSource.ADDITIONAL_TRACKING_POINT_NOT_NEEDED), priority)

        i += 1
        currentDistance += next_spacing

    assert(len(point_list) == len(point_source_list))

    if transfer_to_siblings_childs:
        # Transfer the points transferred to the siblings also to their childs
        for sibling in siblings_list:
            transfer_transferred_points_also_to_childs(
                sibling, used_offset, offset_by_half, max_stitching_distance, to_transfer_points_per_sibling[id(sibling)][::linesign])


>>>>>>> Fixed bug for line string equidistant rastering
def cyclic_distance(val1, val2, line_length):
    absdiff = abs(val1-val2)
    # if val1 < val2:
    #print("Needed cyclic_distance")
    return min(absdiff, line_length-absdiff)


#Takes the offsetted curves organized as tree, connects and samples them.
#Strategy: A connection from parent to child is made where both curves come closest together.
#Input:
#-tree: contains the offsetted curves in a hierachical organized data structure.
#-usedoffset: used offset when there offsetted curves were generated
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
            for item in tree.transferred_point_priority_deque:
                newDEPQ.insert(item[0], math.fmod(
                    item[1]-startDistance+currentCoords.length, currentCoords.length))
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
        point_pair = nearest_points(currentCoords, subnode.val)
        if point_pair[0].distance(point_pair[1]) > 1.5*absoffset:
            print("WARNING: sub geometry closest point is too far apart!")
        projDistance = currentCoords.project(point_pair[0])
        nearestPointsList.append((point_pair, projDistance, subnode))

    temp, temp_origin = Sampler.rasterLineString2_Priority2(currentCoords, absoffset*constants.factor_offset_starting_points,  # We add/subtract absoffset/2 to not sample the same point again (avoid double points for start and end)
                                                            currentCoords.length, stitchdistance, stitching_direction, tree.transferred_point_priority_deque, absoffset)
    assert(len(temp) == len(temp_origin))
    if len(temp) > 3:
        ownCoords.extend(temp[1:])
        ownCoordsOrigin.extend(temp_origin[1:])
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
    PointTransferrer.transfer_inner_points_to_surrounding2(
        tree, usedoffset, offset_by_half, stitchdistance, False, False)  # Only used to transfer to possible siblings
    assert(len(ownCoords) == len(ownCoordsOrigin))
    to_transfer_point_list = []
    to_transfer_point_list_origin = []
    for k in range(1, len(temp)):
        # if abs(temp[k][0]-5.25) < 0.5 and abs(temp[k][1]-42.9) < 0.5:
        #    print("HIER gefunden!")

        #if(
        #   temp_origin[k] == Sampler.PointSource.NOT_NEEDED or
        #   temp_origin[k] == Sampler.PointSource.ADDITIONAL_TRACKING_POINT_NOT_NEEDED):
        #    continue
        #if (offset_by_half and (temp_origin[(k+1) % len(temp_origin)] == Sampler.PointSource.NOT_NEEDED or
        #                        temp_origin[(k+1) % len(temp_origin)] == Sampler.PointSource.ADDITIONAL_TRACKING_POINT_NOT_NEEDED)):
        #    continue
        if (not offset_by_half and temp_origin[k] == Sampler.PointSource.EDGE_NEEDED):
            continue
        #if (offset_by_half and ((temp_origin[k] == Sampler.PointSource.EDGE_NEEDED and temp_origin[(k+1) % len(temp_origin)] != Sampler.PointSource.EDGE_NEEDED) or
        #                        (temp_origin[k] != Sampler.PointSource.EDGE_NEEDED and temp_origin[(k+1) % len(temp_origin)] == Sampler.PointSource.EDGE_NEEDED))):
        #    continue
        to_transfer_point_list.append(Point(temp[k]))
        point_origin = temp_origin[k]
        #if (offset_by_half and point_origin != Sampler.PointSource.EDGE_PREVIOUSLY_SHIFTED):
        #    if(temp_origin[(k+1)%len(temp_origin)]==Sampler.PointSource.EDGE_NEEDED):
        #        point_origin = Sampler.PointSource.EDGE_NEEDED
        #    elif temp_origin[(k+1)%len(temp_origin)]==Sampler.PointSource.EDGE_PREVIOUSLY_SHIFTED:
        #        point_origin = Sampler.PointSource.EDGE_PREVIOUSLY_SHIFTED

            
        to_transfer_point_list_origin.append(point_origin)

    assert(len(ownCoords) == len(ownCoordsOrigin))
    #Next we transfer to the childs
    PointTransferrer.transfer_transferred_points_also_to_childs(
        tree, usedoffset, offset_by_half, stitchdistance, to_transfer_point_list, False,to_transfer_point_list_origin)
    assert(len(ownCoords) == len(ownCoordsOrigin))
    
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
        temp_origin.append(Sampler.PointSource.EDGE_NEEDED)
        # if (abs(temp[-1][0]-9.8) < 0.2 and abs(temp[-1][1]-141.1) < 0.2):
        #    print("HIIER FOUNDED3")

    if len(nearestPointsList) == 0:
        resultCoords = temp
        resultCoords_Origin = temp_origin
    else:
        # if len(nearestPointsList) > 1:
        #    print("HIERRR!")
        temp.append(temp[0])
        temp_origin.append(temp_origin[0])

        nearestPointsList.sort(reverse=False, key=lambda tup: tup[1])

        # temp does not start with currentCoords but has an offset (see call of rasterLineString2_Priority2)
        total_distance = absoffset*constants.factor_offset_starting_points
        current_item_index = 0
        resultCoords = [temp[0]]
        resultCoords_Origin = [temp_origin[0]]
        for i in range(1, len(temp)):
            next_distance = math.sqrt((temp[i][0]-temp[i-1][0])**2 +
                                      (temp[i][1]-temp[i-1][1])**2)
            while (current_item_index < len(nearestPointsList) and
                    total_distance+next_distance+constants.eps > nearestPointsList[current_item_index][1]):
               # was_inside = True
                item = nearestPointsList[current_item_index]
                temp3, temp3_origin = connect_raster_tree_nearest_neighbor(
                    item[2], usedoffset, stitchdistance, item[0][1], offset_by_half)
                resultCoords.append(item[0][0].coords[0])
                resultCoords_Origin.append(Sampler.PointSource.EDGE_NEEDED)
                # reversing avoids crossing when entering and leaving the child segment
                resultCoords.extend(temp3[::-1])
                resultCoords_Origin.extend(temp3_origin[::-1])

                delta = min(absoffset, item[0][0].distance(Point(temp[i])))
                if current_item_index < len(nearestPointsList)-1:
                    delta = min(delta, abs(
                        nearestPointsList[current_item_index+1][1]-item[1]))

                if delta > constants.line_lengh_seen_as_one_point:
                    resultCoords.append(currentCoords.interpolate(
                        item[1]+delta/2).coords[0])
                    resultCoords_Origin.append(Sampler.PointSource.EDGE_NEEDED)

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

#Takes a line and calculates the nearest distance along this line to enter the next_line
#Input:
#-travel_line: The "parent" line for which the distance should be minimized to enter next_line
#-next_line: contains the next_line which need to be entered
#-thresh: The distance between travel_line and next_line needs to below thresh to be a valid point for entering
#Output:
#-tuple - the tuple structure is: (nearest point in travel_line, nearest point in next_line)
def get_nearest_points_closer_than_thresh(travel_line, next_line,thresh):
    point_list = list(MultiPoint(travel_line.coords)) 

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
def create_nearest_points_list(travel_line, children_list, threshold, preferred_direction=0):
    result_list_in_order = []
    result_list_reversed_order = []

    travel_line_reversed = LinearRing(travel_line.coords[::-1])

    weight_in_order = 0
    weight_reversed_order = 0
    for child in children_list:
        result = get_nearest_points_closer_than_thresh(travel_line, child.val, threshold)
        assert(result != None)
        proj = travel_line.project(result[0])
        weight_in_order += proj
        result_list_in_order.append((result,proj, child))

        result = get_nearest_points_closer_than_thresh(travel_line_reversed, child.val, threshold)
        assert(result != None)
        proj = travel_line_reversed.project(result[0])
        weight_reversed_order += proj
        result_list_reversed_order.append((result,proj ,child))

    if preferred_direction == 1:
        weight_in_order=min(weight_in_order/2, max(0, weight_in_order-10*threshold))
    elif preferred_direction == -1:
        weight_reversed_order=min(weight_reversed_order/2, max(0, weight_reversed_order-10*threshold))

    if weight_in_order <= weight_reversed_order:
        return (1, result_list_in_order)
    else:
        return (-1, result_list_reversed_order)


#Takes the offsetted curves organized as tree, connects and samples them.
#Strategy: A connection from parent to child is made as fast as possible to reach the innermost child as fast as possible in order
# to stich afterwards from inner to outer. 
#Input:
#-tree: contains the offsetted curves in a hierachical organized data structure.
#-usedoffset: used offset when there offsetted curves were generated
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
            for item in tree.transferred_point_priority_deque:
                newDEPQ.insert(item[0], math.fmod(
                    item[1]-startDistance+currentCoords.length, currentCoords.length))
            tree.transferred_point_priority_deque = newDEPQ
        #print("Gecutted")

    ownCoords = []
    ownCoordsOrigin = []

    #We try to use always the opposite stitching direction with respect to the parent to avoid crossings when entering and leaving the child
    parent_stitching_direction = -1
    if tree.parent != None:
        parent_stitching_direction = tree.parent.stitching_direction
    stitching_direction, nearestPointsList = create_nearest_points_list(currentCoords, tree.children, 1.5*absoffset,-parent_stitching_direction)

    temp = []
    temp_origin = []
    startoffset = absoffset*constants.factor_offset_starting_points
    if stitching_direction == 1:
        temp, temp_origin = Sampler.rasterLineString2_Priority2(currentCoords, startoffset,  # We add startoffset to not sample the same point again (avoid double points for start and end)
                                                            currentCoords.length, stitchdistance, stitching_direction, tree.transferred_point_priority_deque, absoffset)
    else:
        temp, temp_origin = Sampler.rasterLineString2_Priority2(currentCoords, currentCoords.length-startoffset,  # We subtract startoffset to not sample the same point again (avoid double points for start and end)
                                                            0, stitchdistance, stitching_direction, tree.transferred_point_priority_deque, absoffset)
        currentCoords.coords = currentCoords.coords[::-1] 

<<<<<<< HEAD
=======
    temp, temp_origin = Sampler.rasterLineString2_Priority2(currentCoords, absoffset*constants.factor_offset_starting_points,  # We add/subtract absoffset/2 to not sample the same point again (avoid double points)
                                                            currentCoords.length, stitchdistance, stitching_direction, tree.transferred_point_priority_deque)
>>>>>>> Fixed bug for line string equidistant rastering
    assert(len(temp) == len(temp_origin))
    if len(temp) > 3:
        ownCoords.extend(temp[1:])
        ownCoordsOrigin.extend(temp_origin[1:])
    else:
        ownCoords.extend(temp)
        ownCoordsOrigin.extend(temp_origin)

    assert(len(ownCoords) == len(ownCoordsOrigin))
    
<<<<<<< HEAD
    #Next we need to transfer our rastered points to siblings and childs
    #First we start with the siblings:
=======
>>>>>>> Fixed bug for line string equidistant rastering
    tree.val = LineString(ownCoords)
    tree.pointsourcelist = ownCoordsOrigin
    tree.stitching_direction = stitching_direction
    tree.already_rastered = True
    assert(len(ownCoords) == len(ownCoordsOrigin))
<<<<<<< HEAD
    PointTransferrer.transfer_inner_points_to_surrounding2(
        tree, usedoffset, offset_by_half, stitchdistance, False, False)  # Only used to transfer to possible siblings
=======
    transfer_inner_points_to_surrounding2(
        tree, usedoffset, offset_by_half, stitchdistance, False)  # Only used to transfer to possible siblings
>>>>>>> Fixed bug for line string equidistant rastering
    assert(len(ownCoords) == len(ownCoordsOrigin))
    to_transfer_point_list = []
    to_transfer_point_list_origin = []
    for k in range(1, len(temp)):
        # if abs(temp[k][0]-5.25) < 0.5 and abs(temp[k][1]-42.9) < 0.5:
        #    print("HIER gefunden!")

        #if(
        #   temp_origin[k] == Sampler.PointSource.NOT_NEEDED or
        #   temp_origin[k] == Sampler.PointSource.ADDITIONAL_TRACKING_POINT_NOT_NEEDED):
        #    continue
        #if (offset_by_half and (temp_origin[(k+1) % len(temp_origin)] == Sampler.PointSource.NOT_NEEDED or
        #                        temp_origin[(k+1) % len(temp_origin)] == Sampler.PointSource.ADDITIONAL_TRACKING_POINT_NOT_NEEDED)):
        #    continue
        if (not offset_by_half and temp_origin[k] == Sampler.PointSource.EDGE_NEEDED):
            continue
        #if (offset_by_half and ((temp_origin[k] == Sampler.PointSource.EDGE_NEEDED and temp_origin[(k+1) % len(temp_origin)] != Sampler.PointSource.EDGE_NEEDED) or
        #                        (temp_origin[k] != Sampler.PointSource.EDGE_NEEDED and temp_origin[(k+1) % len(temp_origin)] == Sampler.PointSource.EDGE_NEEDED))):
        #    continue
        to_transfer_point_list.append(Point(temp[k]))
        point_origin = temp_origin[k]
        #if (offset_by_half and point_origin != Sampler.PointSource.EDGE_PREVIOUSLY_SHIFTED):
        #    if(temp_origin[(k+1)%len(temp_origin)]==Sampler.PointSource.EDGE_NEEDED):
        #        point_origin = Sampler.PointSource.EDGE_NEEDED
        #    elif temp_origin[(k+1)%len(temp_origin)]==Sampler.PointSource.EDGE_PREVIOUSLY_SHIFTED:
        #        point_origin = Sampler.PointSource.EDGE_PREVIOUSLY_SHIFTED

<<<<<<< HEAD
            
        to_transfer_point_list_origin.append(point_origin)

    assert(len(ownCoords) == len(ownCoordsOrigin))
    if stitching_direction == -1: #since the projection is only in ccw direction towards inner we need to use "-usedoffset"
        PointTransferrer.transfer_transferred_points_also_to_childs(
            tree, -usedoffset, offset_by_half, stitchdistance, to_transfer_point_list, False,to_transfer_point_list_origin)
    else:
        PointTransferrer.transfer_transferred_points_also_to_childs(
            tree, usedoffset, offset_by_half, stitchdistance, to_transfer_point_list, False,to_transfer_point_list_origin)
    assert(len(ownCoords) == len(ownCoordsOrigin))
   
    # We add the last point wich might not rastered by the previous method
=======
    assert(len(ownCoords) == len(ownCoordsOrigin))
    transfer_transferred_points_also_to_childs(
        tree, usedoffset, offset_by_half, stitchdistance, to_transfer_point_list, False)
    assert(len(ownCoords) == len(ownCoordsOrigin))
    # We add the last point wich might not rastered by the previous method
    # resultCoords.append(currentCoords.interpolate(
    #   currentCoords.length-absoffset/2).coords[0])
    # temp.append((currentCoords.coords[0]))
>>>>>>> Fixed bug for line string equidistant rastering
    lastpoint_distance = currentCoords.project(Point(temp[-1]))
    # We use cyclic distance since the point can be assigned distance 0 if the last point is also the first point
    lastpoint_distance = cyclic_distance(
        lastpoint_distance, currentCoords.length, currentCoords.length)
    assert(len(ownCoords) == len(ownCoordsOrigin))
    delta = min(lastpoint_distance, absoffset)
    if delta > constants.line_lengh_seen_as_one_point:
        temp.append(currentCoords.interpolate(
            currentCoords.length-delta/2).coords[0]) 
        temp_origin.append(Sampler.PointSource.EDGE_NEEDED)
<<<<<<< HEAD
        #if (abs(temp[-1][0]-61.7) < 0.2 and abs(temp[-1][1]-105.1) < 0.2):
=======
        # if (abs(temp[-1][0]-9.8) < 0.2 and abs(temp[-1][1]-141.1) < 0.2):
>>>>>>> Fixed bug for line string equidistant rastering
        #    print("HIIER FOUNDED3")

    if len(nearestPointsList) == 0:
        resultCoords = temp
        resultCoords_Origin = temp_origin
    else:
        # if len(nearestPointsList) > 1:
        #    print("HIERRR!")
        temp.append(temp[0])
        temp_origin.append(temp_origin[0])

        nearestPointsList.sort(reverse=False, key=lambda tup: tup[1])
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
                    total_distance+next_distance+constants.eps > nearestPointsList[current_item_index][1]):
               # was_inside = True
                item = nearestPointsList[current_item_index]
                temp3, temp3_origin = connect_raster_tree_from_inner_to_outer(
                    item[2], usedoffset, stitchdistance, item[0][1], offset_by_half)

                if(Point(resultCoords[-1]).distance(item[0][0]) > absoffset*constants.factor_offset_remove_points):    
                    resultCoords.append(item[0][0].coords[0])
                    resultCoords_Origin.append(Sampler.PointSource.EDGE_NEEDED)
                #if (abs(resultCoords[-1][0]-61.7) < 0.2 and abs(resultCoords[-1][1]-105.1) < 0.2):
                #    print("HIIER FOUNDED3")
                
                # reversing avoids crossing when entering and leaving the child segment
                resultCoords.extend(temp3)
                resultCoords_Origin.extend(temp3_origin)

                delta = min(absoffset, item[0][0].distance(Point(temp[i])))
                if current_item_index < len(nearestPointsList)-1:
                    delta = min(delta, abs(
                        nearestPointsList[current_item_index+1][1]-item[1]))

                if delta > constants.line_lengh_seen_as_one_point:
                    resultCoords.append(currentCoords.interpolate(
                        item[1]+delta/2).coords[0]) 
                    resultCoords_Origin.append(Sampler.PointSource.EDGE_NEEDED)
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

    return resultCoords, resultCoords_Origin
