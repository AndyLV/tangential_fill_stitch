#from sys import path
#from typing import NoReturn
#import matplotlib.pyplot as plt
from shapely.geometry.polygon import LineString
from shapely.geometry import Point
from shapely.ops import substring
import math
import numpy as np
from depq import DEPQ
from enum import IntEnum
from stitches import constants
from stitches import PointTransfer

#Used to tag the origin of a rastered point
class PointSource(IntEnum):
    #MUST_USE = 0  # Legacy
    REGULAR_SPACING = 1  # introduced to not exceed maximal stichting distance
    INITIAL_RASTERING = 2  # No transferred points for this segment were available
    EDGE_NEEDED = 3 # point which must be stitched to avoid to large deviations to the desired path (relevant for interlaced mode)
    #NOT_NEEDED = 4 #Legacy
    ALREADY_TRANSFERRED = 5 #The point transfer can be done recursively down to all sub-childs. 
    #These points are tagged as ALREADY_TRANSFERRED to avoid that they are transferred again to all sub-childs if the next child in the order is rastered.
    
    #ADDITIONAL_TRACKING_POINT_NOT_NEEDED = 6 #Legacy
    #EDGE_RASTERING_ALLOWED = 7 #Legacy
    EDGE_PREVIOUSLY_SHIFTED = 8 #If an edge_needed was previously shifted the next edge (childs edge) should not be shifted in interlaced mode
    ENTER_LEAVING_POINT = 9 #Whether this point is used to enter or leave a child
    SOFT_EDGE = 10
    HARD_EDGE = 11
    PROJECTED_POINT = 12
    MAX_DISTANCE = 13
    FORBIDDEN_POINT=14
    SSOFT_EDGE = 15
    HHARD_EDGE = 16
    FFORBIDDEN_POINT=17
    REPLACED_FORBIDDEN_POINT=18
    DIRECT = 19
    OVERNEXT = 20


# Calculates the angles between adjacent edges at each interior point
#Note that the first and last values in the return array are zero since for the boundary points no angle calculations were possible
def calculate_line_angles(line):
    Angles = np.zeros(len(line.coords))
    for i in range(1, len(line.coords)-1):
        vec1 = np.array(line.coords[i])-np.array(line.coords[i-1])
        vec2 = np.array(line.coords[i+1])-np.array(line.coords[i])
        vec1length = np.linalg.norm(vec1)
        vec2length = np.linalg.norm(vec2)
        Angles[i] = math.acos(np.dot(vec1, vec2)/(vec1length*vec2length))
    return Angles

# Calculates the angles between adjacent edges at each interior point with a sign indicating whether the segment turns to the left (positive) or to the right (negative)
#Note that the first and last values in the return array are zero since for the boundary points no angle calculations were possible
def calculate_signed_line_angles(line):
    Angles = np.zeros(len(line.coords))
    for i in range(1, len(line.coords)-1):
        vec1 = np.array(line.coords[i])-np.array(line.coords[i-1])
        vec2 = np.array(line.coords[i+1])-np.array(line.coords[i])
        vec1length = np.linalg.norm(vec1)
        vec2length = np.linalg.norm(vec2)
        z = np.cross(vec1,vec2)
        Angles[i] = math.copysign(math.acos(np.dot(vec1, vec2)/(vec1length*vec2length)),z)
    return Angles


#Takes all points within "line" and adds additional points if the spacing between two points in line exceeds maxstitchdistance
#Output:
#-List of tuples with the rastered point coordinates
#-List which defines the point origin for each point according to the PointSource enum.
def rasterLineString2(line, maxstitchdistance):
    if line.length < constants.line_lengh_seen_as_one_point:
        return [line.coords[0]], [PointSource.EDGE_NEEDED]
   
    returnpointlist = []
    returnpointsourcelist = []
    overall_distance = 0
    skipped_segment_length = 0
    skipped_segment_start_index = -1

    for i in range(len(line.coords)-1):
        startindex = i
        segmentlength = math.sqrt(pow(
            line.coords[startindex][0]-line.coords[startindex+1][0], 2)+pow(line.coords[startindex][1]-line.coords[startindex+1][1], 2))

        if segmentlength < constants.line_lengh_seen_as_one_point:
            overall_distance += segmentlength
            if skipped_segment_length == 0:
                skipped_segment_start_index = i
                skipped_segment_length = segmentlength
                continue
            else:
                skipped_segment_length += segmentlength
                if skipped_segment_length < constants.line_lengh_seen_as_one_point:
                    continue
                else:
                    startindex = skipped_segment_start_index
                    segmentlength = skipped_segment_length
                    skipped_segment_length = 0
        else:
            skipped_segment_length = 0
            # we subtract "eps=constants.line_lengh_seen_as_one_point" to account for numerical inaccuracies -
            # e.g. segment_length = 2*maxstitchdistance would otherwise sometimes cause 1 and sometimes cause 2 subpoints
        numberofsubpoints = math.ceil(
            (segmentlength-constants.line_lengh_seen_as_one_point)/maxstitchdistance)-1
        subsegmentlength = segmentlength/(numberofsubpoints+1)
        returnpointlist.append((line.coords[startindex]))
        returnpointsourcelist.append(PointSource.INITIAL_RASTERING)
        #if abs(returnpointlist[-1][0]-29)< 0.2 and abs(returnpointlist[-1][1]-11)<0.2:
        #    print("Initial Rastering gefunden!")
        for j in range(1, numberofsubpoints+1):
            returnpointlist.append(
                (line.interpolate(overall_distance+j*subsegmentlength).coords[0]))
            #if abs(returnpointlist[-1][0]-29)< 0.2 and abs(returnpointlist[-1][1]-11)<0.2:
            #    print("Initial Rastering gefunden!")
            returnpointsourcelist.append(PointSource.INITIAL_RASTERING)
        overall_distance = overall_distance+segmentlength
    return returnpointlist, returnpointsourcelist


#Checks whether an edge (unshifted point) can be replaced by shifted_point without have too strong deviations of unshifted_point to the resulting line segment.
#line contains both points (unshifted_point and shifted_point)
#The threshold whether the shift is allowed is a comparison of the resulting distance and a factor times absoffset
#To create a line segment to which unshifted_point is compared we need one point left (e.g. shifted_point) and one point right to unshifted point on the line "line"-
#This second point to the right is calculated by traveling the same spacing between unshifted_point and shifted_point in the opposite direction. This traveling distance is cropped by "maxdelta_next_point"
#Output:
#-True or false whether the shift of the edge is allowed (not a too large deviation)
def check_edge_needed_shift_allowed(line, unshifted_point, shifted_point, absoffset, maxdelta_next_point):
    #assert(abs(line.coords[0][0]-line.coords[-1][0])<constants.eps and abs(line.coords[0][1]-line.coords[-1][1])<constants.eps)
    proj1 = line.project(unshifted_point)
    proj2 = line.project(shifted_point)
    delta = proj2-proj1
    if maxdelta_next_point != -1:
        delta = math.copysign(min(abs(delta), maxdelta_next_point),delta)
    newproj = proj1-delta
    if newproj < 0:
        #newproj+=line.length
        return False
    if newproj > line.length:
        #newproj-=line.length
        return False
    thirdPoint = line.interpolate(newproj)
    bisection = LineString([shifted_point, thirdPoint])
    distance = bisection.distance(unshifted_point)

    return (distance<=constants.fac_offset_edge_shift*absoffset)


def remove_dense_points(connectedLine, connectedLineOrigin,stitchdistance,absoffset):
    i = 1
    while i < len(connectedLine)-1:
        line_segment = LineString([connectedLine[i-1], connectedLine[i+1]])
        if line_segment.length > stitchdistance:
            i+=1
            continue
        distance = Point(connectedLine[i]).distance(line_segment)
        if distance < absoffset*constants.factor_offset_remove_dense_points:
            connectedLine.pop(i)
            connectedLineOrigin.pop(i)
        else:
            i+=1
    return connectedLine, connectedLineOrigin

#Rasters a line between start_distance and end_distance.
#Input:
#-line: The line to be rastered
#-start_distance: The distance along the line from which the rastering should start
#-end_distance: The distance along the line until which the rastering should be done
#-maxstitchdistance: The maximum allowed stitch distance
#-stitching_direction: =1 is stitched along line direction, =-1 if stitched in reversed order. Note that
# start_distance > end_distance for stitching_direction = -1
#-must_use_points_deque: deque with projected points on line from its neighbors. An item of the deque
#is setup as follows: ((projected point on line, LineStringSampling.PointSource), priority=distance along line)
#index of point_origin is the index of the point in the neighboring line
#-absoffset: used offset between to offsetted curves
#Output:
#-List of tuples with the rastered point coordinates
#-List which defines the point origin for each point according to the PointSource enum.
def rasterLineString2_Priority3(line, start_distance, end_distance, maxstitchdistance, stitching_direction, must_use_points_deque, absoffset):
    if (abs(end_distance-start_distance) < constants.line_lengh_seen_as_one_point):
        return [line.interpolate(start_distance).coords[0]], [PointSource.EDGE_NEEDED]

    assert (stitching_direction == -1 and start_distance >= end_distance) or (
        stitching_direction == 1 and start_distance <= end_distance)
   
    deque_points = list(must_use_points_deque)

    if not deque_points:
        return rasterLineString2(substring(line, start_distance, end_distance), maxstitchdistance)

    linecoords = line.coords
    # for item in linecoords:
    #    if abs(item[0]-135) < 0.5 and abs(item[1]-16) < 0.5:
    #        print("GEFUNDEN!")

    if start_distance > end_distance:
        start_distance, end_distance = line.length - \
            start_distance, line.length-end_distance
        linecoords = linecoords[::-1]
        for i in range(len(deque_points)):
            deque_points[i] = (deque_points[i][0],
                               line.length-deque_points[i][1])
    else:
        deque_points = deque_points[::-1]

    # Remove all points from the "must use point list" which do not fall in the segment [start_distance; end_distance]
    while (len(deque_points) > 0 and deque_points[0][1] <= start_distance+min(maxstitchdistance/20, constants.point_spacing_to_be_considered_equal)):
        deque_points.pop(0)
    while (len(deque_points) > 0 and deque_points[-1][1] >= end_distance-min(maxstitchdistance/20, constants.point_spacing_to_be_considered_equal)):
        deque_points.pop()

    if not deque_points:
        return rasterLineString2(substring(line, start_distance, end_distance), maxstitchdistance)

# Ordering in priority queue:
#   (point, LineStringSampling.PointSource), priority)

    returnpointlist = []
    returnpointsourcelist = []

    path_coords = substring(LineString(linecoords),
                            start_distance, end_distance)
    signed_Angles = calculate_signed_line_angles(path_coords)

    

    segment_start_index = 0
    segment_end_index = 1
    dq_start_index = 0
    dq_end_index_right = 0 #index within the deque whose point is straight "behind" the segment end
    dq_end_index_left = 0 #index within the deque whose point is straight "before" the segment end
    startpoint_proj = 0

    startpoint = path_coords.coords[segment_start_index]
    startpoint_source = PointSource.EDGE_NEEDED

    returnpointsourcelist.append(startpoint_source)
    returnpointlist.append(startpoint)
    #if (abs(returnpointlist[-1][0]-13.2)< 0.2 and abs(returnpointlist[-1][1]-141.4)<0.2):
    #            print("HIIER FOUNDED1")  

    #last_shifted_edge_points = [None, None]

    #We divide the line into segments separated by angles which are not negligble straight
    while segment_end_index < len(path_coords.coords):
        current_angle_sum = 0

        #Get the starting index within deque_points belonging to the current segment
        while dq_start_index < len(deque_points) and deque_points[dq_start_index][1]-start_distance < startpoint_proj:
            dq_start_index += 1

        current_point_list = [path_coords.coords[segment_start_index], path_coords.coords[segment_end_index]]
        current_proj_endpoint = path_coords.project(Point(current_point_list[-1]))
        dq_iter = dq_start_index
        found_edge_shifted = False
        while dq_iter < len(deque_points) and deque_points[dq_iter][1]-start_distance < current_proj_endpoint:
            if deque_points[dq_iter][0].point_source == PointSource.EDGE_PREVIOUSLY_SHIFTED:
                    found_edge_shifted=True
            dq_iter+=1
        min_deviation = min(abs(current_proj_endpoint+start_distance-deque_points[min(dq_iter,len(deque_points)-1)][1]),abs(current_proj_endpoint+start_distance-deque_points[max(0,dq_iter-1)][1]))
        min_index = segment_end_index

 
        while segment_end_index < len(path_coords.coords)-1 and len(LineString(current_point_list).simplify(constants.factor_offset_remove_dense_points*absoffset,preserve_topology=False).coords)==2:
            segment_end_index+=1
            current_point_list.append(path_coords.coords[segment_end_index])
            current_proj_endpoint = path_coords.project(Point(current_point_list[-1]))
            dq_iter-=1
            while dq_iter < len(deque_points) and deque_points[dq_iter][1]-start_distance < current_proj_endpoint:
                if deque_points[dq_iter][0].point_source == PointSource.EDGE_PREVIOUSLY_SHIFTED:
                    found_edge_shifted=True
                dq_iter+=1
            deviation = min(abs(current_proj_endpoint+start_distance-deque_points[max(0,dq_iter-1)][1]),abs(current_proj_endpoint+start_distance-deque_points[min(dq_iter,len(deque_points)-1)][1]))
            if  deviation < min_deviation:
                min_index = segment_end_index
                min_deviation = deviation
        if found_edge_shifted: #Shifted points are by definition projected points - so we want to take only for unshifted points the one which is closest to a projected one
            segment_end_index = min_index

        #while segment_end_index < len(signed_Angles)-1 and abs(current_angle_sum+signed_Angles[segment_end_index]) <= constants.limiting_angle:
        #    current_angle_sum+=signed_Angles[segment_end_index]
        #    segment_end_index += 1



        #Get end index within deque_points for the current segment:
        dq_end_index_right = dq_start_index
        end_point_proj = path_coords.project(Point(path_coords.coords[segment_end_index]))
        while dq_end_index_right < len(deque_points) and deque_points[dq_end_index_right][1]-start_distance < end_point_proj:
            dq_end_index_right += 1    
        dq_end_index_left = dq_end_index_right-1

        #if (abs(path_coords.coords[segment_end_index][0]-82.259859301875)< 0.2 and abs(path_coords.coords[segment_end_index][1]-98.03874253852013)<0.2):
       #     print("HIIER FOUNDED1")  

        if (abs(path_coords.coords[segment_end_index][0]-77.5)< 0.2 and abs(path_coords.coords[segment_end_index][1]-105.4)<0.2):
            print("HIIER FOUNDED1")  

        
        #We check the deque_points close to the end point (edge) of the segment. If there is no point source coming from an shifted edge we are allowed to shift the edge of this segment.
        #By this we create an interlaced structure - every second row will be shifted.
        #TODO TODO check whether starting with dq_end_index_left or dq_end_index_right
        dq_iter = dq_end_index_left
        found_EDGE_Needed_and_not_edge_shifted = False
        found_index = -1
        
        if not(dq_end_index_right < len(deque_points) and deque_points[dq_end_index_right][1]-start_distance <= end_point_proj+absoffset and deque_points[dq_end_index_right][0].point_source == PointSource.EDGE_PREVIOUSLY_SHIFTED):
            while dq_iter > 0 and deque_points[dq_iter][1]-start_distance >= max(startpoint_proj,end_point_proj- maxstitchdistance/2.0)  and dq_end_index_right-dq_iter < 4: #look max 3 points to the left of the end segment
                if deque_points[dq_iter][0].point_source == PointSource.EDGE_PREVIOUSLY_SHIFTED:
                    found_EDGE_Needed_and_not_edge_shifted = False
                    break
                elif not found_EDGE_Needed_and_not_edge_shifted and (deque_points[dq_iter][0].point_source == PointSource.EDGE_NEEDED or deque_points[dq_iter][0].point_source == PointSource.INITIAL_RASTERING):
                    found_index = dq_iter
                    found_EDGE_Needed_and_not_edge_shifted = True
                dq_iter -= 1

        endpoint = path_coords.coords[segment_end_index]
        endpoint_source = PointSource.EDGE_NEEDED

        if  found_EDGE_Needed_and_not_edge_shifted:
            if segment_end_index < len(path_coords.coords)-1:
                found_EDGE_Needed_and_not_edge_shifted = check_edge_needed_shift_allowed(path_coords,Point(endpoint), Point(deque_points[found_index][0].point.coords[0]),absoffset, Point(endpoint).distance(Point(path_coords.coords[segment_end_index+1])))
            else:
                found_EDGE_Needed_and_not_edge_shifted = check_edge_needed_shift_allowed(path_coords,Point(endpoint), Point(deque_points[found_index][0].point.coords[0]),absoffset, -1)
            if found_EDGE_Needed_and_not_edge_shifted:
                #if (abs( deque_points[found_index].point[0].coords[0][0]-90.0)< 0.2 and abs( deque_points[found_index].point[0].coords[0][1]-174.5)<0.2):
                #    print("HIIER FOUNDED1")
                #if(deque_points[found_index].point[0].distance(Point(startpoint)) > constants.fact_offset_edge_shift_remaining_line_length*absoffset):
                #last_shifted_edge_points.pop(0)
                #last_shifted_edge_points.append(Point(endpoint))
                
                #If the shift is allowed we replace the segments endpoint by a projected point which lies within the segment (shifting towards inner)       
                endpoint =  deque_points[found_index][0].point.coords[0]
                end_point_proj = path_coords.project(Point(endpoint))
                endpoint_source = PointSource.EDGE_PREVIOUSLY_SHIFTED

        #Collect all deque points within the segment in an extra list 
        proj_list = []
        while dq_start_index < len(deque_points) and deque_points[dq_start_index][1]-start_distance < end_point_proj:
            # TODO maybe take not all points (e.g. EDGE NEEDED)
            if deque_points[dq_start_index][0].point_source != PointSource.EDGE_PREVIOUSLY_SHIFTED:
                proj_list.append(
                    deque_points[dq_start_index][1]-startpoint_proj-start_distance)
            dq_start_index += 1


        #Distribute the points within the current segment using the projected points from the deque
        distributed_proj_list = distribute_points_proj2(
            Point(startpoint), Point(endpoint), proj_list, maxstitchdistance)

        # if the end point of the previous segment (= start point of this segment) was shifted (shift is always towards the inner of the segment),
        # we need to check here whether this shift was ok from the perspective of this segment
        if returnpointsourcelist[-1] == PointSource.EDGE_PREVIOUSLY_SHIFTED and distributed_proj_list:
            last_shifted_edge_point = Point(path_coords.coords[segment_start_index])
            next_point = path_coords.interpolate(distributed_proj_list[0]+startpoint_proj) #the first point after the shifted edge
            bisectorline = LineString([next_point, Point(returnpointlist[-1])])
            distance = bisectorline.distance(last_shifted_edge_point)
            if distance > constants.fac_offset_edge_shift*absoffset: #The shift induced deviation to the path is too large - we either unshift the shifted point or add the unshifted point additionally
                
                if len(returnpointlist) > 1 and last_shifted_edge_point.distance(Point(returnpointlist[-2])) > maxstitchdistance :
                    #We cannot simply unshift the point since this would induce a too large distance between the next-to-last point in the previous segement and the unshifted point (which moved away from this point compared to the shifted one)
                    #Hence we add the unshifted point in addition to fulfill the maxstitchdistance constraint.
                    returnpointlist.append(last_shifted_edge_point.coords[0])
                    returnpointsourcelist.append(PointSource.EDGE_NEEDED)
                else: 
                    #Unshift the point
                    returnpointlist.pop()
                    returnpointsourcelist.pop()
                    returnpointlist.append(last_shifted_edge_point.coords[0])
                    returnpointsourcelist.append(PointSource.EDGE_NEEDED)
                
        #Add the points to the return list
        for item in distributed_proj_list:
            returnpointlist.append(path_coords.interpolate(
                item+startpoint_proj).coords[0])
            #if (abs(returnpointlist[-1][0]-34)< 0.2 and abs(returnpointlist[-1][1]-165.6)<0.2):
            #    print("HIIER FOUNDED1")    
            returnpointsourcelist.append(PointSource.REGULAR_SPACING)
        
        #Add the end point only if it is not too close
        if (Point(endpoint).distance(Point(returnpointlist[-1])) > absoffset*constants.factor_offset_remove_points 
                    or endpoint_source == PointSource.EDGE_PREVIOUSLY_SHIFTED): #The "or endpoint_source == PointSource.EDGE_PREVIOUSLY_SHIFTED" is necessary since shifting can bring points very close together 
                                                                                #and if we do not add it here we do not know that the previous point was shifted in the next round
            returnpointlist.append(endpoint)
            returnpointsourcelist.append(endpoint_source)
        startpoint_proj = end_point_proj
        startpoint = endpoint
        segment_start_index = segment_end_index
        segment_end_index += 1
        #if (abs(returnpointlist[-1][0]-10)< 0.2 and abs(returnpointlist[-1][1]-143.8)<0.2):
        #    print("HIIER FOUNDED1")  


    # TODO: Check whether the angle splitting was too lax so that we might need to add further EDGE_NEEDED afterwards to minimize the deviation to the original path

    assert(len(returnpointlist) == len(returnpointsourcelist))

    #return remove_dense_points(returnpointlist, returnpointsourcelist, maxstitchdistance,absoffset)
    return returnpointlist, returnpointsourcelist


#Returns the index in arr whose value is closest to target
def get_closest_value_index(arr, target):
    n = len(arr)
    left = 0
    right = n - 1
    mid = 0

    # edge case - last or above all
    if target >= arr[n - 1]:
        return arr[n - 1], n-1
    # edge case - first or below all
    if target <= arr[0]:
        return arr[0], 0
    # BSearch solution: Time & Space: Log(N)

    while left < right:
        mid = (left + right) // 2  # find the mid
        if target < arr[mid]:
            right = mid
        elif target > arr[mid]:
            left = mid + 1
        else:
            return arr[mid], mid

    if target < arr[mid]:
        if target - arr[mid-1] >= arr[mid] - target:
            return arr[mid], mid
        else:
            return arr[mid-1], mid-1
        # return find_closest(arr[mid - 1], arr[mid], target)
    else:
        if target - arr[mid] >= arr[mid+1] - target:
            return arr[mid+1], mid+1
        else:
            return arr[mid], mid
        # return find_closest(arr[mid], arr[mid + 1], target)



# findClosest
# We find the closest by taking the difference
# between the target and both values. It assumes
# that val2 is greater than val1 and target lies
# between these two.
def find_closest(val1, val2, target):
    return val2 if target - val1 >= val2 - target else val1


#Tries to optimally distribute the points between start_point and end_point 
# - thereby considering projected points on this linesegment (transferPointList_proj) and
# - the maximum allowed stitch distance (maxStitchdistance)
#Output:
#-List of the distance of the calculated, distributed points with respect to the start_point (projection on the linesegment)
def distribute_points_proj2(start_point, end_point, transferPointList_proj, maxStitchdistance):

   # if abs(end_point.coords[0][0]-58.4) < 0.2 and abs(end_point.coords[0][1]-45.5)< 0.2:
   #     print("GEFUNDEN")

    L = start_point.distance(end_point)

    min_number_of_points_in_between = math.floor(L/maxStitchdistance)

    #We remove transfer points which are very close at the boundary since they sometimes are projection artifacts and 
    #should not influence the decisions below with len(transferPointList_proj) == 0, 1 or 2
    limit3_fac = 0.1
    while len(transferPointList_proj) > 0 and transferPointList_proj[0] <= limit3_fac*min(L,maxStitchdistance):
        transferPointList_proj.pop(0)

    while len(transferPointList_proj) > 0 and L-transferPointList_proj[-1] <= limit3_fac*min(L,maxStitchdistance):
        transferPointList_proj.pop()

    if len(transferPointList_proj) == 0:
        resultlist = []
        delta_step = L/(min_number_of_points_in_between+1)
        for j in range(1, min_number_of_points_in_between+1):
            resultlist.append(delta_step*j)
        return resultlist

    limit1 = 0.5*maxStitchdistance
    limit2_fac = 0.25

    if L <= limit1:
        return []

    limit_assymetric_fac = 0.25

    if len(transferPointList_proj) == 1 and min_number_of_points_in_between <= 1:
        if L <= maxStitchdistance and abs((L/2)-transferPointList_proj[0]) > limit_assymetric_fac*L:
            return []
        else:
            return transferPointList_proj

    if len(transferPointList_proj) == 2 and L <= maxStitchdistance:
        return []

    # if min_number_of_points_in_between == 0 and len(transferPointList_proj) != 1:
    #    return []
    # elif min_number_of_points_in_between == 0:
    #    return transferPointList_proj

    val, index = get_closest_value_index(transferPointList_proj, L/2)
    #point = LineString([start_point, end_point]).interpolate(val)

    delta = 0
    if index == len(transferPointList_proj)-1:
        delta = val-transferPointList_proj[index-1]
    elif index == 0:
        delta = transferPointList_proj[index+1]-val
    else:
        delta = max(
            val-transferPointList_proj[index-1],  transferPointList_proj[index+1]-val)

    if delta < constants.point_spacing_to_be_considered_equal:
        delta = maxStitchdistance

    number_of_points_with_delta = 1+math.floor((val-limit2_fac*delta)/delta)+math.floor((L-val-limit2_fac*delta)/delta)

    if number_of_points_with_delta %2 != 0:#len(transferPointList_proj) % 2 != 0:
        if delta < limit1:
            i = 2
            while(delta*i < limit1):
                i += 1
            delta *= i

        if delta > maxStitchdistance:
            i = 2
            while(delta/i > maxStitchdistance):
                i += 1
            delta /= i

    #residuum_left = val % delta
    #residuum_right = (L-val) % delta
    #NLeft = math.floor(val/delta)
    #NRight = math.floor((L-val)/delta)

    delta_left = delta_right = delta
    #delta_new_left = delta_new_right = delta
    #if residuum_left <= delta/2.0 and NLeft != 0:
    #        delta_new_left = min(maxStitchdistance, val/NLeft)
    #elif residuum_left > delta/2.0:
    #        delta_new_left = min(maxStitchdistance, val/(NLeft+1))

    #if residuum_right <= delta/2.0 and NRight != 0:
    #        delta_new_right = min(maxStitchdistance,(L-val)/NRight)
    #elif residuum_right > delta/2.0:
    #        delta_new_right = min(maxStitchdistance,(L-val)/(NRight+1))

    #mixing_factor_left = 1.0-math.exp(-30*(residuum_left-delta/2.0)**2/delta**2)
    #mixing_factor_right = 1.0-math.exp(-30*(residuum_right-delta/2.0)**2/delta**2)
            
    #delta_left = mixing_factor_left*delta_new_left+(1.0-mixing_factor_left)*delta
    #delta_right = mixing_factor_right*delta_new_right+(1.0-mixing_factor_right)*delta   

    returnlist = [val]

    current_val = val-delta_left
    while current_val > 0:
        if current_val < limit2_fac*delta_left and (delta_left+current_val) <= maxStitchdistance:
            break
        returnlist.insert(0, current_val)
        current_val -= delta_left

    current_val = val+delta_right
    while current_val < L:
        if (L-current_val) < limit2_fac*delta_right and (delta_right+L-current_val) <= maxStitchdistance:
            break
        returnlist.append(current_val)
        current_val += delta_right

    return returnlist



#Rasters a line between start_distance and end_distance.
#Input:
#-line: The line to be rastered
#-start_distance: The distance along the line from which the rastering should start
#-end_distance: The distance along the line until which the rastering should be done
#-maxstitchdistance: The maximum allowed stitch distance
#-stitching_direction: =1 is stitched along line direction, =-1 if stitched in reversed order. Note that
# start_distance > end_distance for stitching_direction = -1
#-must_use_points_deque: deque with projected points on line from its neighbors. An item of the deque
#is setup as follows: ((projected point on line, LineStringSampling.PointSource), priority=distance along line)
#index of point_origin is the index of the point in the neighboring line
#-absoffset: used offset between to offsetted curves
#Output:
#-List of tuples with the rastered point coordinates
#-List which defines the point origin for each point according to the PointSource enum.
def rasterLineString2_Priority2(line, start_distance, end_distance, maxstitchdistance, stitching_direction, must_use_points_deque, absoffset):
    if (abs(end_distance-start_distance) < constants.line_lengh_seen_as_one_point):
        return [line.interpolate(start_distance).coords[0]], [PointSource.EDGE_NEEDED]

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

    # Remove all points from the "must use point list" which do not fall in the segment [start_distance; end_distance]
    while (len(deque_points) > 0 and deque_points[0][1] <= start_distance+min(maxstitchdistance/20, constants.point_spacing_to_be_considered_equal)):
        deque_points.pop(0)
    while (len(deque_points) > 0 and deque_points[-1][1] >= end_distance-min(maxstitchdistance/20, constants.point_spacing_to_be_considered_equal)):
        deque_points.pop()


# Ordering in priority queue:
#   (point, LineStringSampling.PointSource), priority)
    aligned_line = LineString(linecoords)
    path_coords = substring(aligned_line,
                            start_distance, end_distance)
    angles = calculate_line_angles(path_coords)

    current_distance = start_distance

    #Next we merge the line points and the projected points into one list
    merged_point_list = []
    dq_iter = 0
    for point,angle in zip(path_coords.coords,angles):
        if abs(point[0]-40.4) < 0.2 and abs(point[1]-2.3)< 0.2:
            print("GEFUNDEN")
        current_distance = start_distance+path_coords.project(Point(point))
        while dq_iter < len(deque_points) and deque_points[dq_iter][1] < current_distance:
            if deque_points[dq_iter][0].point_source == PointSource.FORBIDDEN_POINT:
                if (merged_point_list[-1][0].point_source == PointSource.SOFT_EDGE and 
                   abs(merged_point_list[-1][1]-deque_points[dq_iter][1] < absoffset*constants.factor_offset_forbidden_point)):
                    item = merged_point_list.pop()
                    merged_point_list.append((PointTransfer.projected_point_tuple(point=item[0].point, point_source=\
                        PointSource.FORBIDDEN_POINT),item[1]))
            else:
                merged_point_list.append(deque_points[dq_iter])
            dq_iter+=1
        if (dq_iter < len(deque_points) and 
            deque_points[dq_iter-1][0].point_source == PointSource.FORBIDDEN_POINT and
            angle < constants.limiting_angle and
            abs(deque_points[dq_iter-1][1]-current_distance) < absoffset*constants.factor_offset_forbidden_point):
            point_source = PointSource.FORBIDDEN_POINT
        else:
            if angle < constants.limiting_angle:
                point_source = PointSource.SOFT_EDGE
            else:
                point_source = PointSource.HARD_EDGE
        merged_point_list.append((PointTransfer.projected_point_tuple(point=Point(point), point_source=point_source),current_distance))

    result_list = [merged_point_list[0]]
   
    segment_start_index = 0
    segment_end_index = 1
    forbidden_point_list = []
    while segment_end_index < len(merged_point_list): 
        print(segment_end_index)
        #if segment_end_index == 101:
        #    print("HIER")

        if abs(merged_point_list[segment_end_index-1][0].point.coords[0][0]-67.9) < 0.2 and abs(merged_point_list[segment_end_index-1][0].point.coords[0][1]-161.0)< 0.2:
            print("GEFUNDEN")



        current_point_list = [merged_point_list[segment_start_index][0].point]
       
        while segment_end_index < len(merged_point_list):
            segment_length = merged_point_list[segment_end_index][1]-merged_point_list[segment_start_index][1]
            if segment_length > maxstitchdistance+constants.point_spacing_to_be_considered_equal:
                new_distance = merged_point_list[segment_start_index][1]+maxstitchdistance
                merged_point_list.insert(segment_end_index,(PointTransfer.projected_point_tuple(point=aligned_line.interpolate(new_distance), point_source=\
                PointSource.MAX_DISTANCE),new_distance))
                #new_distance = (merged_point_list[segment_start_index][1]+merged_point_list[segment_end_index][1])/2.0
                #merged_point_list.insert(segment_end_index,(PointTransfer.projected_point_tuple(point=path_coords.interpolate(new_distance), point_source=\
                #PointSource.MAX_DISTANCE),new_distance))
                if abs(merged_point_list[segment_end_index][0].point.coords[0][0]-64.6) < 0.2 and abs(merged_point_list[segment_end_index][0].point.coords[0][1]-161)< 0.2:
                    print("GEFUNDEN")
                segment_end_index+=1
                break
            if abs(merged_point_list[segment_end_index][0].point.coords[0][0]-93.6) < 0.2 and abs(merged_point_list[segment_end_index][0].point.coords[0][1]-122.7)< 0.2:
                print("GEFUNDEN")
            #if abs(merged_point_list[segment_end_index][0].point.coords[0][0]-89.15) < 0.2 and abs(merged_point_list[segment_end_index][0].point.coords[0][1]-107.42)< 0.2:
            #    print("GEFUNDEN")
            

            current_point_list.append(merged_point_list[segment_end_index][0].point)
            simplified_len = len(LineString(current_point_list).simplify(constants.factor_offset_remove_dense_points*absoffset,preserve_topology=False).coords)
            if simplified_len > 2: #not all points have been simplified - so we need to add it
                break

            if merged_point_list[segment_end_index][0].point_source ==PointSource.HARD_EDGE:
                segment_end_index+=1
                break
            segment_end_index+=1


        if True:
            segment_end_index-=1
            index_overnext = -1
            index_direct = -1
            index_hard_edge = -1

            iter = segment_start_index+1 
            while (iter <= segment_end_index):# and 
                #(merged_point_list[iter][0].point_source == PointSource.SOFT_EDGE or 
                # merged_point_list[iter][0].point_source == PointSource.MAX_DISTANCE) or
                # merged_point_list[iter][0].point_source == PointSource.FORBIDDEN_POINT):
                if merged_point_list[iter][0].point_source == PointSource.OVERNEXT:
                    index_overnext = iter
                elif merged_point_list[iter][0].point_source == PointSource.DIRECT:
                    index_direct = iter
                elif merged_point_list[iter][0].point_source == PointSource.HARD_EDGE:
                    index_hard_edge = iter
                iter += 1
                print ("iter: ",iter)
            print ("Segment start: ",segment_start_index, " Segment end: ", segment_end_index)
            #if iter > segment_start_index:
            #    segment_end_index = iter
            if index_hard_edge != -1:
                #if index_overnext != -1:
                #    segment_end_index = min(index_overnext,index_hard_edge)
                #elif index_direct != -1:
                #    segment_end_index = min(index_direct,index_hard_edge)
                #else:
                segment_end_index = index_hard_edge
            else:
                if index_overnext != -1:
                    if (index_direct != -1 and index_direct > index_overnext and 
                            (merged_point_list[index_direct][1]-merged_point_list[index_overnext][1]) >= 
                            constants.factor_segment_length_direct_preferred_over_overnext*
                            (merged_point_list[index_overnext][1]-merged_point_list[segment_start_index][1])):
                        segment_end_index = index_direct
                    else:
                        segment_end_index = index_overnext
                elif index_direct != -1:
                    segment_end_index = index_direct

        else:
            segment_end_index-=1
            iter = segment_end_index
            while (iter >= segment_start_index and 
                (merged_point_list[iter][0].point_source == PointSource.SOFT_EDGE or 
                merged_point_list[iter][0].point_source == PointSource.MAX_DISTANCE) or
                merged_point_list[iter][0].point_source == PointSource.FORBIDDEN_POINT):
                iter -= 1

            if iter > segment_start_index:
                segment_end_index = iter

        if (merged_point_list[segment_start_index][0].point_source == PointSource.OVERNEXT and
            merged_point_list[segment_end_index][0].point_source == PointSource.DIRECT and
                (merged_point_list[segment_end_index][1] - merged_point_list[segment_start_index][1]) < absoffset):
            result_list.pop()

        result_list.append(merged_point_list[segment_end_index])
        if merged_point_list[segment_end_index][0].point_source == PointSource.FORBIDDEN_POINT:
            forbidden_point_list.append(len(result_list)-1)

        segment_start_index = segment_end_index
        segment_end_index+=1

    return_point_list = [result_list[0][0].point.coords[0]]
    return_point_source_list = [PointSource.EDGE_NEEDED]

    #result_list = replace_forbidden_points(aligned_line, result_list, forbidden_point_list,absoffset)
    #Finally we need to traverse the result_list and add points where the maxstitchdistance constraint is violated:
    for i in range(1,len(result_list)):
        #segment_length = result_list[i][1]-result_list[i-1][1]
        #if segment_length > maxstitchdistance+constants.eps:
        #    numberofsubpoints = math.ceil(
        #    (segment_length-constants.line_lengh_seen_as_one_point)/maxstitchdistance)-1
        #    subsegmentlength = segment_length/(numberofsubpoints+1)
        #    for j in range(1, numberofsubpoints+1):
        #        return_point_list.append(
        #            (path_coords.interpolate(result_list[i-1][1]+j*subsegmentlength).coords[0]))
        #        return_point_source_list.append(PointSource.REGULAR_SPACING)

        return_point_list.append(result_list[i][0].point.coords[0])
        #if abs(result_list[i][0].point.coords[0][0]-91.7) < 0.2 and abs(result_list[i][0].point.coords[0][1]-106.15)< 0.2:
        #    print("GEFUNDEN")
        if result_list[i][0].point_source == PointSource.HARD_EDGE:
            point_source = PointSource.HHARD_EDGE        
        elif result_list[i][0].point_source == PointSource.SOFT_EDGE:
            point_source = PointSource.SSOFT_EDGE
        elif result_list[i][0].point_source == PointSource.MAX_DISTANCE:
            point_source = PointSource.REGULAR_SPACING
        elif result_list[i][0].point_source == PointSource.FORBIDDEN_POINT:
            point_source = PointSource.FORBIDDEN_POINT
        else:
            point_source = PointSource.PROJECTED_POINT

        return_point_source_list.append(point_source)


    assert(len(return_point_list) == len(return_point_source_list))

    #return remove_dense_points(returnpointlist, returnpointsourcelist, maxstitchdistance,absoffset)
    return return_point_list, return_point_source_list

def replace_forbidden_points(line, result_list, forbidden_point_list_indices, absoffset):
    current_index_shift = 0 #since we add and remove points in the result_list, we need to adjust the indices stored in forbidden_point_list_indices
    for index in forbidden_point_list_indices:
        if abs(result_list[index][0].point.coords[0][0]-40.7) < 0.2 and abs(result_list[index][0].point.coords[0][1]-1.3)< 0.2:
            print("GEFUNDEN")
        index+=current_index_shift
        distance_left = result_list[index][0].point.distance(result_list[index-1][0].point)/2.0
        distance_right = result_list[index][0].point.distance(result_list[(index+1)%len(result_list)][0].point)/2.0
        while distance_left > constants.point_spacing_to_be_considered_equal and distance_right > constants.point_spacing_to_be_considered_equal:
            new_point_left_proj = result_list[index][1]-distance_left
            if new_point_left_proj < 0:
                new_point_left_proj += line.length
            new_point_right_proj = result_list[index][1]+distance_right
            if new_point_right_proj > line.length:
                new_point_right_proj-=line.length
            point_left = line.interpolate(new_point_left_proj)
            point_right = line.interpolate(new_point_right_proj)
            forbidden_point_distance =  result_list[index][0].point.distance(LineString([point_left, point_right]))
            if forbidden_point_distance < constants.factor_offset_remove_dense_points*absoffset:
                del result_list[index]
                result_list.insert(index, (PointTransfer.projected_point_tuple(point=point_right, point_source=\
                PointSource.REPLACED_FORBIDDEN_POINT),new_point_right_proj))
                result_list.insert(index, (PointTransfer.projected_point_tuple(point=point_left, point_source=\
                PointSource.REPLACED_FORBIDDEN_POINT),new_point_left_proj))
                current_index_shift+=1
                break
            else:
                distance_left/=2.0
                distance_right/=2.0
    return result_list


if __name__ == "__main__":
    line = LineString([(0,0), (1,0), (2,1),(3,0),(4,0)])

    print(calculate_line_angles(line)*180.0/math.pi)

    print(calculate_signed_line_angles(line)*180.0/math.pi)