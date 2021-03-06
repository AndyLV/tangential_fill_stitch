from shapely.geometry.polygon import LineString
from shapely.geometry import Point
from shapely.ops import substring
import math
import numpy as np
from enum import IntEnum
from stitches import constants
from stitches import PointTransfer

#Used to tag the origin of a rastered point
class PointSource(IntEnum):
    #MUST_USE = 0  # Legacy
    REGULAR_SPACING = 1  # introduced to not exceed maximal stichting distance
    #INITIAL_RASTERING = 2  #Legacy
    EDGE_NEEDED = 3 # point which must be stitched to avoid to large deviations to the desired path
    #NOT_NEEDED = 4 #Legacy
    #ALREADY_TRANSFERRED = 5 #Legacy
    #ADDITIONAL_TRACKING_POINT_NOT_NEEDED = 6 #Legacy
    #EDGE_RASTERING_ALLOWED = 7 #Legacy
    #EDGE_PREVIOUSLY_SHIFTED = 8  #Legacy
    ENTER_LEAVING_POINT = 9 #Whether this point is used to enter or leave a child
    SOFT_EDGE_INTERNAL = 10 #If the angle at a point is <= constants.limiting_angle this point is marked as SOFT_EDGE
    HARD_EDGE_INTERNAL = 11 #If the angle at a point is > constants.limiting_angle this point is marked as HARD_EDGE (HARD_EDGES will always be stitched)
    PROJECTED_POINT = 12 #If the point was created by a projection (transferred point) of a neighbor it is marked as PROJECTED_POINT
    REGULAR_SPACING_INTERNAL = 13 # introduced to not exceed maximal stichting distance
    #FORBIDDEN_POINT_INTERNAL=14  #Legacy
    SOFT_EDGE = 15 #If the angle at a point is <= constants.limiting_angle this point is marked as SOFT_EDGE
    HARD_EDGE = 16 #If the angle at a point is > constants.limiting_angle this point is marked as HARD_EDGE (HARD_EDGES will always be stitched)
    FORBIDDEN_POINT=17  #Only relevant for desired interlacing - non-shifted point positions at the next neighbor are marked as forbidden
    REPLACED_FORBIDDEN_POINT=18 #If one decides to avoid forbidden points new points to the left and to the right as replacement are created
    DIRECT = 19     #Calculated by next neighbor projection
    OVERNEXT = 20   #Calculated by overnext neighbor projection


# Calculates the angles between adjacent edges at each interior point
#Note that the first and last values in the return array are zero since for the boundary points no angle calculations were possible
def calculate_line_angles(line):
    Angles = np.zeros(len(line.coords))
    for i in range(1, len(line.coords)-1):
        vec1 = np.array(line.coords[i])-np.array(line.coords[i-1])
        vec2 = np.array(line.coords[i+1])-np.array(line.coords[i])
        vec1length = np.linalg.norm(vec1)
        vec2length = np.linalg.norm(vec2)
        scalar_prod=np.dot(vec1, vec2)/(vec1length*vec2length)
        scalar_prod = min(max(scalar_prod,-1),1)
        #if scalar_prod > 1.0:
        #    scalar_prod = 1.0
        #elif scalar_prod < -1.0:
        #    scalar_prod = -1.0
        Angles[i] = math.acos(scalar_prod)
    return Angles

#Rasters a line between start_distance and end_distance.
#Input:
#-line: The line to be rastered
#-start_distance: The distance along the line from which the rastering should start
#-end_distance: The distance along the line until which the rastering should be done
#-maxstitch_distance: The maximum allowed stitch distance
#-stitching_direction: =1 is stitched along line direction, =-1 if stitched in reversed order. Note that
# start_distance > end_distance for stitching_direction = -1
#-must_use_points_deque: deque with projected points on line from its neighbors. An item of the deque
#is setup as follows: ((projected point on line, LineStringSampling.PointSource), priority=distance along line)
#index of point_origin is the index of the point in the neighboring line
#-abs_offset: used offset between to offsetted curves
#Output:
#-List of tuples with the rastered point coordinates
#-List which defines the point origin for each point according to the PointSource enum.
def raster_line_string_with_priority_points(line, start_distance, end_distance, maxstitch_distance, stitching_direction, must_use_points_deque, abs_offset):
    if (abs(end_distance-start_distance) < constants.line_lengh_seen_as_one_point):
        return [line.interpolate(start_distance).coords[0]], [PointSource.HARD_EDGE]

    assert (stitching_direction == -1 and start_distance >= end_distance) or (
        stitching_direction == 1 and start_distance <= end_distance)
   
    deque_points = list(must_use_points_deque)

    linecoords = line.coords

    if start_distance > end_distance:
        start_distance, end_distance = line.length - \
            start_distance, line.length-end_distance
        linecoords = linecoords[::-1]
        for i in range(len(deque_points)):
            deque_points[i] = (deque_points[i][0],
                               line.length-deque_points[i][1])
    else:
        deque_points = deque_points[::-1]

    # Remove all points from the deque which do not fall in the segment [start_distance; end_distance]
    while (len(deque_points) > 0 and deque_points[0][1] <= start_distance+min(maxstitch_distance/20, constants.point_spacing_to_be_considered_equal)):
        deque_points.pop(0)
    while (len(deque_points) > 0 and deque_points[-1][1] >= end_distance-min(maxstitch_distance/20, constants.point_spacing_to_be_considered_equal)):
        deque_points.pop()


# Ordering in priority queue:
#   (point, LineStringSampling.PointSource), priority)
    aligned_line = LineString(linecoords)
    path_coords = substring(aligned_line,
                            start_distance, end_distance)
    angles = calculate_line_angles(path_coords)

    current_distance = start_distance

    #Next we merge the line points and the projected (deque) points into one list
    merged_point_list = []
    dq_iter = 0
    for point,angle in zip(path_coords.coords,angles):
        #if abs(point[0]-40.4) < 0.2 and abs(point[1]-2.3)< 0.2:
        #    print("GEFUNDEN")
        current_distance = start_distance+path_coords.project(Point(point))
        while dq_iter < len(deque_points) and deque_points[dq_iter][1] < current_distance:
            #We want to avoid setting points at soft edges close to forbidden points
            if deque_points[dq_iter][0].point_source == PointSource.FORBIDDEN_POINT:
                #Check whether a previous added point is a soft edge close to the forbidden point
                if (merged_point_list[-1][0].point_source == PointSource.SOFT_EDGE_INTERNAL and 
                   abs(merged_point_list[-1][1]-deque_points[dq_iter][1] < abs_offset*constants.factor_offset_forbidden_point)):
                    item = merged_point_list.pop()
                    merged_point_list.append((PointTransfer.projected_point_tuple(point=item[0].point, point_source=\
                        PointSource.FORBIDDEN_POINT),item[1]))
            else:
                merged_point_list.append(deque_points[dq_iter])
            dq_iter+=1
        #Check whether the current point is close to a forbidden point
        if (dq_iter < len(deque_points) and 
            deque_points[dq_iter-1][0].point_source == PointSource.FORBIDDEN_POINT and
            angle < constants.limiting_angle and
            abs(deque_points[dq_iter-1][1]-current_distance) < abs_offset*constants.factor_offset_forbidden_point):
            point_source = PointSource.FORBIDDEN_POINT
        else:
            if angle < constants.limiting_angle:
                point_source = PointSource.SOFT_EDGE_INTERNAL
            else:
                point_source = PointSource.HARD_EDGE_INTERNAL
        merged_point_list.append((PointTransfer.projected_point_tuple(point=Point(point), point_source=point_source),current_distance))

    result_list = [merged_point_list[0]]
   
    #General idea: Take one point of merged_point_list after another into the current segment until this segment is not simplified to a straight line by shapelys simplify method.
    #Then, look at the points within this segment and choose the best fitting one (HARD_EDGE > OVERNEXT projected point > DIRECT projected point) as termination of this segment 
    # and start point for the next segment (so we do not always take the maximum possible length for a segment)
    segment_start_index = 0
    segment_end_index = 1
    forbidden_point_list = []
    while segment_end_index < len(merged_point_list): 
        #if abs(merged_point_list[segment_end_index-1][0].point.coords[0][0]-67.9) < 0.2 and abs(merged_point_list[segment_end_index-1][0].point.coords[0][1]-161.0)< 0.2:
        #    print("GEFUNDEN")

        #Collection of points for the current segment
        current_point_list = [merged_point_list[segment_start_index][0].point]
       
        while segment_end_index < len(merged_point_list):
            segment_length = merged_point_list[segment_end_index][1]-merged_point_list[segment_start_index][1]
            if segment_length > maxstitch_distance+constants.point_spacing_to_be_considered_equal:
                new_distance = merged_point_list[segment_start_index][1]+maxstitch_distance
                merged_point_list.insert(segment_end_index,(PointTransfer.projected_point_tuple(point=aligned_line.interpolate(new_distance), point_source=\
                PointSource.REGULAR_SPACING_INTERNAL),new_distance))
                if abs(merged_point_list[segment_end_index][0].point.coords[0][0]-12.2) < 0.2 and abs(merged_point_list[segment_end_index][0].point.coords[0][1]-0.9)< 0.2:
                    print("GEFUNDEN")
                segment_end_index+=1
                break
            #if abs(merged_point_list[segment_end_index][0].point.coords[0][0]-93.6) < 0.2 and abs(merged_point_list[segment_end_index][0].point.coords[0][1]-122.7)< 0.2:
            #    print("GEFUNDEN")
        
            current_point_list.append(merged_point_list[segment_end_index][0].point)
            simplified_len = len(LineString(current_point_list).simplify(constants.factor_offset_remove_dense_points*abs_offset,preserve_topology=False).coords)
            if simplified_len > 2: #not all points have been simplified - so we need to add it
                break

            if merged_point_list[segment_end_index][0].point_source ==PointSource.HARD_EDGE_INTERNAL:
                segment_end_index+=1
                break
            segment_end_index+=1

        segment_end_index-=1

        #Now we choose the best fitting point within this segment
        index_overnext = -1
        index_direct = -1
        index_hard_edge = -1

        iter = segment_start_index+1 
        while (iter <= segment_end_index):
            if merged_point_list[iter][0].point_source == PointSource.OVERNEXT:
                index_overnext = iter
            elif merged_point_list[iter][0].point_source == PointSource.DIRECT:
                index_direct = iter
            elif merged_point_list[iter][0].point_source == PointSource.HARD_EDGE_INTERNAL:
                index_hard_edge = iter
            iter += 1
        if index_hard_edge != -1:
            segment_end_index = index_hard_edge
        else:
            if index_overnext != -1:
                if (index_direct != -1 and index_direct > index_overnext and 
                        (merged_point_list[index_direct][1]-merged_point_list[index_overnext][1]) >= 
                        constants.factor_segment_length_direct_preferred_over_overnext*
                        (merged_point_list[index_overnext][1]-merged_point_list[segment_start_index][1])):
                    #We allow to take the direct projected point instead of the overnext projected point if it would result in a
                    #significant longer segment length
                    segment_end_index = index_direct
                else:
                    segment_end_index = index_overnext
            elif index_direct != -1:
                segment_end_index = index_direct

        #Usually OVERNEXT and DIRECT points are close to each other and in some cases both were selected as segment edges
        #If they are too close (<abs_offset) we remove one of it
        if (((merged_point_list[segment_start_index][0].point_source == PointSource.OVERNEXT and
            merged_point_list[segment_end_index][0].point_source == PointSource.DIRECT) or
            (merged_point_list[segment_start_index][0].point_source == PointSource.DIRECT and
            merged_point_list[segment_end_index][0].point_source == PointSource.OVERNEXT)) and
            abs(merged_point_list[segment_end_index][1] - merged_point_list[segment_start_index][1]) < abs_offset):
            result_list.pop()

        result_list.append(merged_point_list[segment_end_index])
        #To have a chance to replace all forbidden points afterwards
        if merged_point_list[segment_end_index][0].point_source == PointSource.FORBIDDEN_POINT:
            forbidden_point_list.append(len(result_list)-1)

        segment_start_index = segment_end_index
        segment_end_index+=1

    return_point_list = [] #[result_list[0][0].point.coords[0]]
    return_point_source_list = []#[result_list[0][0].point_source]

    #Currently replacement of forbidden points not satisfying
    #result_list = replace_forbidden_points(aligned_line, result_list, forbidden_point_list,abs_offset)

    #Finally we create the final return_point_list and return_point_source_list
    for i in range(len(result_list)):
        return_point_list.append(result_list[i][0].point.coords[0])
        #if abs(result_list[i][0].point.coords[0][0]-91.7) < 0.2 and abs(result_list[i][0].point.coords[0][1]-106.15)< 0.2:
        #    print("GEFUNDEN")
        if result_list[i][0].point_source == PointSource.HARD_EDGE_INTERNAL:
            point_source = PointSource.HARD_EDGE        
        elif result_list[i][0].point_source == PointSource.SOFT_EDGE_INTERNAL:
            point_source = PointSource.SOFT_EDGE
        elif result_list[i][0].point_source == PointSource.REGULAR_SPACING_INTERNAL:
            point_source = PointSource.REGULAR_SPACING
        elif result_list[i][0].point_source == PointSource.FORBIDDEN_POINT:
            point_source = PointSource.FORBIDDEN_POINT
        else:
            point_source = PointSource.PROJECTED_POINT

        return_point_source_list.append(point_source)


    assert(len(return_point_list) == len(return_point_source_list))

    #return remove_dense_points(returnpointlist, returnpointsourcelist, maxstitch_distance,abs_offset)
    return return_point_list, return_point_source_list


if __name__ == "__main__":
    line = LineString([(0,0), (1,0), (2,1),(3,0),(4,0)])

    print(calculate_line_angles(line)*180.0/math.pi)
