import math

# Used in the simplify routine of shapely
simplification_threshold = 0.01

# If a transferred point is closer than this value to one of its neighbors, it will be checked whether it can be removed
distance_thresh_remove_transferred_point = 0.15

# If a line segment is shorter than this threshold it is handled as a single point
line_lengh_seen_as_one_point = 0.05

# E.g. to check whether a point is already present in a point list, the point is allowed to be this value in distance apart
point_spacing_to_be_considered_equal = 0.05

# Adjacent geometry should have points closer than offset*offset_factor_for_adjacent_geometry to be considered adjacent
offset_factor_for_adjacent_geometry = 1.5

# Transfer point distance is used for projecting points from already rastered geometry to adjacent geometry
# (max spacing transfer_point_distance_factor*offset) to get a more regular pattern
transfer_point_distance_factor = 1.5

# Used to handle numerical inaccuracies during comparisons
eps = 1E-3

factor_offset_starting_points=0.5 #When entering and leaving a child from a parent we introduce an offset of absoffset*factor_offset_starting_points so 
                                  #that entering and leaving points are not lying above each other.

factor_offset_remove_points=0.75 #if points are closer than absoffset*factor_offset_remove_points one of it is removed

fac_offset_edge_shift = 0.25 #if an unshifted relevant edge is closer than absoffset*fac_offset_edge_shift to the line segment created by the shifted edge,
                            #the shift is allowed - otherwise the edge must not be shifted.

limiting_angle = math.pi*3.5/180.0 #angles straighter (smaller) than this are considered as more or less straight (no concrete edges required for path segments having only angles <= this value)

factor_offset_remove_dense_points=0.05 #if a point distance to the connected line of its two neighbors is smaller than absoffset times this factor, this point will be removed if the stitching distance will not be exceeded