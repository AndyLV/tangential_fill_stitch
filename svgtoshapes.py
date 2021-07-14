''' Read and parse SVG files returning a shapely GeometryCollection containing the paths and primitives found.
Transforms are applied and all shapes are returned in the SVG root elements coordinate space. Primitives (rect,
circle and ellipse) are converted to equivalent paths, and only shape outlines are rendered, stroke properties
have no affect and are discarded. 

Curved paths (arc and bezier commands and rounded rect, ellipse and circle elements) are linearized using a 
naive approach that attempts to create enough segments to reduce maximum distance to the curve to below a 
certain threshold, although this is probably not very efficient and most likely implemented in a fairly
broken way.'''

import math
import re


def matmul(a, b):
    return [[sum([a[row][i]*b[i][col] for i in range(0, 3)]) for col in range(0, 3)] for row in range(0, 3)]


def matident():
    return [[1, 0, 0], [0, 1, 0], [0, 0, 1]]


def matmat(a, b, c, d, e, f):
    return [[a, c, e], [b, d, f], [0, 0, 1]]


def matoffs(x, y):
    return [[1, 0, x], [0, 1, y], [0, 0, 1]]


def matrot(theta, xo=0, yo=0):
    cosp = math.cos(theta * (math.pi/180))
    sinp = math.sin(theta * (math.pi/180))
    xoff = xo - xo * cosp + yo * sinp
    yoff = yo - xo * cosp - yo * cosp
    return [[cosp, -sinp, xoff], [sinp, cosp, yoff], [0, 0, 1]]


def matscale(xfact=1, yfact=1, xo=0, yo=0):
    xoff = xo - xo * xfact
    yoff = yo - yo * yfact
    return [[xfact, 0, xoff], [0, yfact, yoff], [0, 0, 1]]


def matskew(xs=0, ys=0, xo=0, yo=0):
    xoff = -yo * math.tan(xs * (math.pi/180))
    yoff = -xo * math.tan(ys * (math.pi/180))
    return [[1, math.tan(xs*(math.pi/180)), xoff], [math.tan(ys*(math.pi/180)), 1, yoff], [0, 0, 1]]


def matpoint(p, m):
    ''' apply the affine 3x3 transform matrix represented by m[][] to
    the point p[] (a tuple of (x, y)), returning the transformed point
    as a new tuple (newx,newy) '''
    px, py = p
    nx = m[0][0] * px + m[0][1] * py + m[0][2]
    ny = m[1][0] * px + m[1][1] * py + m[1][1]
    return (nx, ny)


def parsesvgtransform(s, initmatrix):
    done = False
    offset = 0
    mats = []
    while not done:
        m = re.match("^\\s*(\\w+?)\\s*\\(([^\\)]*)\\)\\s*", s[offset:])
        if m is None:
            if len(s[offset:].strip()) > 0:
                raise("Error parsing transform string at: '{}'".format(
                    s[offset:].strip()))
            done = True
            continue
        offset = offset + len(m.group(0))
        op = m.group(1)
        #params = [float(p.strip()) for p in re.split("\\s*,\\s*", m.group(2))]
        params = [float(p.strip()) for p in re.split("\\s+", m.group(2))]
        #print("{}({})".format(op, repr(params)))

        def nthparam(n, defval=None):
            nonlocal params
            if len(params) > n:
                return params[n]
            return defval
        if op == "translate":
            mats.append(matoffs(nthparam(0, 0), nthparam(1, 0)))
        elif op == "rotate":
            mats.append(matrot(nthparam(0, 0), nthparam(1, 0), nthparam(2, 0)))
        elif op == "scale":
            mats.append(matscale(nthparam(0, 1), nthparam(
                1, 1), nthparam(2, 0), nthparam(3, 0)))
        elif op == "skewX":
            mats.append(matskew(nthparam(0, 0), 0))
        elif op == "skewY":
            mats.append(matskew(0, nthparam(0, 0)))
        elif op == "matrix":
            mats.append(
                matmat(*[nthparam(n, [1, 0, 0, 1, 0, 0][n]) for n in range(0, 6)]))
        else:
            raise RuntimeError(
                "unknown transform function '{}'".format(m.group(0)))
    # print(repr(m))
    # print(repr(mats))

    result = matident()
    for cur in reversed(mats):
        result = matmul(result, cur)
    if initmatrix != None:
        result = matmul(initmatrix, result)
    return result


def svgtoshapes(svg_file):
    from xml.dom import minidom
    from svg.path import parse_path
    import svg.path
    from shapely import geometry as g
    from shapely import ops

    with open(svg_file, "r") as f:
        svg_string = f.read()

    svg_dom = minidom.parseString(svg_string)

    def composetransform(elem):
        print(elem)
        mat = matident()
        while not (elem.parentNode is None):
            if elem.hasAttribute('transform'):
                mat = matmul(mat, parsesvgtransform(
                    elem.getAttribute('transform'), None))
            elem = elem.parentNode
        return mat

    def attrvalue(elem, attr, default=None, required=False):
        if elem.hasAttributes(attr):
            return elem.getAttribute(attr)
        if default is None and required:
            # print("Error")
            raise ValueError(
                f"cannot get attribute '{attr}' of element {elem}")
        return default

    def rectpath(elem):
        x0, y0 = float(attrvalue(elem, "x", 0)), float(attrvalue(elem, "y", 0))
        w, h = float(attrvalue(elem, "width", required=True)), float(
            attrvalue(elem, "height", required=True))
        x1, y1 = x0 + w, y0
        x2, y2 = x0 + w, y0 + h
        x3, y3 = x0, y0 + h

        d = ("M{} {} L {} {} L {} {} L {} {} z"
             "".format(x0, y0, x1, y1, x2, y2, x3, y3))
        return d

    def elipsepath(elem):
        """converts the parameters from an ellipse or a circle to a string for a 
            Path object d-attribute"""

        cx = attrvalue(elem, 'cx', 0, required=False)
        cy = attrvalue(elem, 'cy', 0, required=False)
        rx = attrvalue(elem, 'rx', None, required=False)
        ry = attrvalue(elem, 'ry', None, required=False)
        r = attrvalue(elem, 'r', None, required=False)

        if r is not None:
            rx = ry = float(r)
        else:
            rx = float(rx)
            ry = float(ry)

        cx = float(cx)
        cy = float(cy)

        d = ''
        d += 'M' + str(cx - rx) + ',' + str(cy)
        d += 'a' + str(rx) + ',' + str(ry) + ' 0 1,0 ' + str(2 * rx) + ',0'
        d += 'a' + str(rx) + ',' + str(ry) + ' 0 1,0 ' + str(-2 * rx) + ',0'

        return d

   # for path in svg_dom.getElementsByTagName('path'):
   #     print(path)
   #     composetransform(path)
    path_strings = [(path.getAttribute('d'), composetransform(path))
                    for path in svg_dom.getElementsByTagName('path')]
    path_strings.extend([(rectpath(rect), composetransform(rect))
                        for rect in svg_dom.getElementsByTagName('rect')])
    path_strings.extend([(ellipsepath(circle), composetransform(circle))
                        for circle in svg_dom.getElementsByTagName('circle')])

    shapes = []
    for path_string, matrix in path_strings:
        path_data = parse_path(path_string)
        # print(repr(path_data))
        #path_shape = path_to_points(path_data)
        # paths.append(path_data)
        points = []

        def commitpoints():
            nonlocal shapes
            nonlocal points
            if len(points) > 1:
                shapes.append(g.asLineString(
                    [(p.real, p.imag) for p in points]))
            elif len(points) > 0:
                shapes.append(g.asPoint((points[0].real, points[0].imag)))
            points = []

        def addpoint(p):
            nonlocal points
            if len(points) < 1 or points[-1] != p:
                points.append(p)

        for seg in path_data:
            isbreak = False
            if len(points) > 0 and points[-1] != seg.start:
                isbreak = True
            if isinstance(seg, svg.path.Move):
                isbreak = True
            if isbreak:
                commitpoints()
            addpoint(seg.start)
            if isinstance(seg, svg.path.Arc):
                numdivs = max(1, int((seg.theta / 360) * 16))
                for i in range(0, numdivs):
                    addpoint(seg.point((i+1)/(numdivs+1)))
            elif isinstance(seg, (svg.path.QuadraticBezier, svg.path.CubicBezier)):
                numdivs = min(10, max(1, int(seg.length() * 10)))
                for i in range(0, numdivs):
                    addpoint(seg.point((i+1)/(numdivs+1)))
            addpoint(seg.end)
        commitpoints()

    return g.GeometryCollection(shapes)
