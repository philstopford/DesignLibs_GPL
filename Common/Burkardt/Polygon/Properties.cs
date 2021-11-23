using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Polygon;

public static class Properties
{
    public static double[] angle_half(double[] p1, double[] p2, double[] p3, int p1Index = 0, int p2Index = 0, int p3Index = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ANGLE_HALF finds half an angle.
        //
        //  Discussion:
        //
        //    The original angle is defined by the sequence of points P1, P2 and P3.
        //
        //    The point P4 is calculated so that:
        //
        //      (P1,P2,P4) = (P1,P2,P3) / 2
        //
        //        P1
        //        /
        //       /   P4
        //      /  .
        //     / .
        //    P2--------->P3
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 May 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[2], P2[2], P3[2], points defining the angle.
        //
        //    Input, double ANGLE_HALF[2], a point P4 defining the half angle.
        //    The vector P4 - P2 will have unit norm.
        //
    {
        int i;

        double[] p4 = new double[2];

        double norm = Math.Sqrt((p1[p1Index + 0] - p2[p2Index + 0]) * (p1[p1Index + 0] - p2[p2Index + 0])
                                + (p1[p1Index + 1] - p2[p2Index + 1]) * (p1[p1Index + 1] - p2[p2Index + 1]));

        for (i = 0; i < 2; i++)
        {
            p4[i] = (p1[p1Index + i] - p2[p2Index + i]) / norm;
        }

        norm = Math.Sqrt((p3[p3Index + 0] - p2[p2Index + 0]) * (p3[p3Index + 0] - p2[p2Index + 0])
                         + (p3[p3Index + 1] - p2[p2Index + 1]) * (p3[p3Index + 1] - p2[p2Index + 1]));

        for (i = 0; i < 2; i++)
        {
            p4[i] += (p3[p3Index + i] - p2[p2Index + i]) / norm;
        }

        for (i = 0; i < 2; i++)
        {
            p4[i] = 0.5 * p4[i];
        }

        norm = typeMethods.r8vec_norm(2, p4);

        for (i = 0; i < 2; i++)
        {
            p4[i] = p2[p2Index + i] + p4[i] / norm;
        }

        return p4;
    }

    public static double angle_rad(double[] p1, double[] p2, double[] p3, int p1Index = 0, int p2Index = 0, int p3Index = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ANGLE_RAD returns the angle in radians swept out between two rays.
        //
        //  Discussion:
        //
        //      ANGLE_RAD ( P1, P2, P3 ) + ANGLE_RAD ( P3, P2, P1 ) = 2 * PI
        //
        //        P1
        //        /
        //       /
        //      /
        //     /
        //    P2--------->P3
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[2], P2[2], P3[2], define the rays
        //    P1 - P2 and P3 - P2 which define the angle.
        //
        //    Output, double ANGLE_RAD, the angle between the two rays,
        //    in radians.  This value will always be between 0 and 2*PI.  If either 
        //    ray has zero length, then the angle is returned as zero.
        //
    {
        double[] p = new double[2];
            
        double value;

        p[0] = (p3[p3Index + 0] - p2[p2Index + 0]) * (p1[p1Index + 0] - p2[p2Index + 0])
               + (p3[p3Index + 1] - p2[p2Index + 1]) * (p1[p1Index + 1] - p2[p2Index + 1]);


        p[1] = (p3[p3Index + 0] - p2[p2Index + 0]) * (p1[p1Index + 1] - p2[p2Index + 1])
               - (p3[p3Index + 1] - p2[p2Index + 1]) * (p1[p1Index + 0] - p2[p2Index + 0]);

        switch (p[0])
        {
            case 0.0 when p[1] == 0.0:
                value = 0.0;
                return value;
        }

        value = Math.Atan2(p[1], p[0]);

        switch (value)
        {
            case < 0.0:
                value += 2.0 * Math.PI;
                break;
        }

        return value;
    }

    public static bool between(double xa, double ya, double xb, double yb, double xc,
            double yc)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BETWEEN is TRUE if vertex C is between vertices A and B.
        //
        //  Discussion:
        //
        //    The points must be (numerically) collinear.
        //
        //    Given that condition, we take the greater of XA - XB and YA - YB
        //    as a "scale" and check where C's value lies.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 May 2014
        //
        //  Author:
        //
        //    Original C version by Joseph ORourke.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Joseph ORourke,
        //    Computational Geometry in C,
        //    Cambridge, 1998,
        //    ISBN: 0521649765,
        //    LC: QA448.D38.
        //
        //  Parameters:
        //
        //    Input, double XA, YA, XB, YB, XC, YC, the coordinates of 
        //    the vertices.
        //
        //    Output, bool BETWEEN, is TRUE if C is between A and B.
        //
    {
        bool value;

        if (!collinear(xa, ya, xb, yb, xc, yc))
        {
            value = false;
        }
        else if (Math.Abs(ya - yb) < Math.Abs(xa - xb))
        {
            double xmax = Math.Max(xa, xb);
            double xmin = Math.Min(xa, xb);
            value = xmin <= xc && xc <= xmax;
        }
        else
        {
            double ymax = Math.Max(ya, yb);
            double ymin = Math.Min(ya, yb);
            value = ymin <= yc && yc <= ymax;
        }

        return value;
    }

    public static bool collinear(double xa, double ya, double xb, double yb, double xc,
            double yc)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COLLINEAR returns a measure of collinearity for three points.
        //
        //  Discussion:
        //
        //    In order to deal with collinear points whose coordinates are not
        //    numerically exact, we compare the area of the largest square
        //    that can be created by the line segment between two of the points
        //    to (twice) the area of the triangle formed by the points.
        //
        //    If the points are collinear, their triangle has zero area.
        //    If the points are close to collinear, then the area of this triangle
        //    will be small relative to the square of the longest segment.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 May 2014
        //
        //  Author:
        //
        //    Original C version by Joseph ORourke.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Joseph ORourke,
        //    Computational Geometry in C,
        //    Cambridge, 1998,
        //    ISBN: 0521649765,
        //    LC: QA448.D38.
        //
        //  Parameters:
        //
        //    Input, double XA, YA, XB, YB, XC, YC, the coordinates of 
        //    the vertices.
        //
        //    Output, bool COLLINEAR, is TRUE if the points are judged 
        //    to be collinear.
        //
    {
        const double r8_eps = 2.220446049250313E-016;
        bool value;

        double area = 0.5 * (
            (xb - xa) * (yc - ya)
            - (xc - xa) * (yb - ya));

        double side_ab_sq = Math.Pow(xa - xb, 2) + Math.Pow(ya - yb, 2);
        double side_bc_sq = Math.Pow(xb - xc, 2) + Math.Pow(yb - yc, 2);
        double side_ca_sq = Math.Pow(xc - xa, 2) + Math.Pow(yc - ya, 2);

        double side_max_sq = Math.Max(side_ab_sq, Math.Max(side_bc_sq, side_ca_sq));

        switch (side_max_sq)
        {
            case <= r8_eps:
                value = true;
                break;
            default:
            {
                value = 2.0 * Math.Abs(area) <= r8_eps * side_max_sq;

                break;
            }
        }

        return value;
    }

    public static bool diagonal(int im1, int ip1, int n, int[] prev, int[] next, double[] x,
            double[] y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIAGONAL: VERTEX(IM1) to VERTEX(IP1) is a proper internal diagonal.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 May 2014
        //
        //  Author:
        //
        //    Original C version by Joseph ORourke.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Joseph ORourke,
        //    Computational Geometry in C,
        //    Cambridge, 1998,
        //    ISBN: 0521649765,
        //    LC: QA448.D38.
        //
        //  Parameters:
        //
        //    Input, int IM1, IP1, the indices of two vertices.
        //
        //    Input, int N, the number of vertices.
        //
        //    Input, int PREV[N], the previous neighbor of each vertex.
        //
        //    Input, int NEXT[N], the next neighbor of each vertex.
        //
        //    Input, double X[N], Y[N], the coordinates of each vertex.
        //
        //    Output, bool DIAGONAL, the value of the test.
        //
    {
        bool value1 = in_cone(im1, ip1, n, prev, next, x, y);
        bool value2 = in_cone(ip1, im1, n, prev, next, x, y);
        bool value3 = diagonalie(im1, ip1, n, next, x, y);

        bool value = value1 && value2 && value3;

        return value;
    }

    public static bool diagonalie(int im1, int ip1, int n, int[] next, double[] x, double[] y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIAGONALIE is true if VERTEX(IM1):VERTEX(IP1) is a proper diagonal.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 May 2014
        //
        //  Author:
        //
        //    Original C version by Joseph ORourke.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Joseph ORourke,
        //    Computational Geometry in C,
        //    Cambridge, 1998,
        //    ISBN: 0521649765,
        //    LC: QA448.D38.
        //
        //  Parameters:
        //
        //    Input, int IM1, IP1, the indices of two vertices.
        //
        //    Input, int N, the number of vertices.
        //
        //    Input, int NEXT[N], the next neighbor of each vertex.
        //
        //    Input, double X[N], Y[N], the coordinates of each vertex.
        //
        //    Output, bool DIAGONALIE, the value of the test.
        //
    {
        int j = im1;
        int jp1 = next[im1];

        bool value = true;
        //
        //  For each edge VERTEX(J):VERTEX(JP1) of the polygon:
        //
        while (true)
        {
            //
            //  Skip any edge that includes vertex IM1 or IP1.
            //
            if (j == im1 || j == ip1 || jp1 == im1 || jp1 == ip1)
            {
            }
            else
            {
                bool value2 = intersect(x[im1], y[im1], x[ip1], y[ip1], x[j], y[j],
                    x[jp1], y[jp1]);

                if (value2)
                {
                    value = false;
                    break;
                }
            }

            j = jp1;
            jp1 = next[j];

            if (j == im1)
            {
                break;
            }
        }

        return value;
    }

    public static bool in_cone(int im1, int ip1, int n, int[] prev, int[] next, double[] x,
            double[] y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    IN_CONE is TRUE if the diagonal VERTEX(IM1):VERTEX(IP1) is strictly internal.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 May 2014
        //
        //  Author:
        //
        //    Original C version by Joseph ORourke.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Joseph ORourke,
        //    Computational Geometry in C,
        //    Cambridge, 1998,
        //    ISBN: 0521649765,
        //    LC: QA448.D38.
        //
        //  Parameters:
        //
        //    Input, int IM1, IP1, the indices of two vertices.
        //
        //    Input, int N, the number of vertices.
        //
        //    Input, int PREV[N], the previous neighbor of each vertex.
        //
        //    Input, int NEXT[N], the next neighbor of each vertex.
        //
        //    Input, double X[N], Y[N], the coordinates of each vertex.
        //
        //    Output, bool IN_CONE, the value of the test.
        //
    {
        bool value;

        int im2 = prev[im1];
        int i = next[im1];

        double t1 = triangle_area(x[im1], y[im1], x[i], y[i], x[im2], y[im2]);

        switch (t1)
        {
            case >= 0.0:
                double t2 = triangle_area(x[im1], y[im1], x[ip1], y[ip1], x[im2], y[im2]);
                double t3 = triangle_area(x[ip1], y[ip1], x[im1], y[im1], x[i], y[i]);
                value = 0.0 < t2 && 0.0 < t3;
                break;
            default:
                double t4 = triangle_area(x[im1], y[im1], x[ip1], y[ip1], x[i], y[i]);
                double t5 = triangle_area(x[ip1], y[ip1], x[im1], y[im1], x[im2], y[im2]);
                value = !(0.0 <= t4 && 0.0 <= t5);
                break;
        }

        return value;
    }

    public static bool intersect(double xa, double ya, double xb, double yb, double xc,
            double yc, double xd, double yd)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INTERSECT is true if lines VA:VB and VC:VD intersect.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 May 2014
        //
        //  Author:
        //
        //    Original C version by Joseph ORourke.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Joseph ORourke,
        //    Computational Geometry in C,
        //    Cambridge, 1998,
        //    ISBN: 0521649765,
        //    LC: QA448.D38.
        //
        //  Parameters:
        //
        //    Input, double XA, YA, XB, YB, XC, YC, XD, YD, the X and Y 
        //    coordinates of the four vertices.
        //
        //    Output, bool INTERSECT, the value of the test.
        //
    {
        bool value;

        if (intersect_prop(xa, ya, xb, yb, xc, yc, xc, yd))
        {
            value = true;
        }
        else if (between(xa, ya, xb, yb, xc, yc))
        {
            value = true;
        }
        else if (between(xa, ya, xb, yb, xd, yd))
        {
            value = true;
        }
        else if (between(xc, yc, xd, yd, xa, ya))
        {
            value = true;
        }
        else if (between(xc, yc, xd, yd, xb, yb))
        {
            value = true;
        }
        else
        {
            value = false;
        }

        return value;
    }

    public static bool intersect_prop(double xa, double ya, double xb, double yb, double xc,
            double yc, double xd, double yd)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INTERSECT_PROP is TRUE if lines VA:VB and VC:VD have a proper intersection.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 May 2014
        //
        //  Author:
        //
        //    Original C version by Joseph ORourke.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Joseph ORourke,
        //    Computational Geometry in C,
        //    Cambridge, 1998,
        //    ISBN: 0521649765,
        //    LC: QA448.D38.
        //
        //  Parameters:
        //
        //    Input, double XA, YA, XB, YB, XC, YC, XD, YD, the X and Y 
        //    coordinates of the four vertices.
        //
        //    Output, bool INTERSECT_PROP, the result of the test.
        //
    {
        bool value;

        if (collinear(xa, ya, xb, yb, xc, yc))
        {
            value = false;
        }
        else if (collinear(xa, ya, xb, yb, xd, yd))
        {
            value = false;
        }
        else if (collinear(xc, yc, xd, yd, xa, ya))
        {
            value = false;
        }
        else if (collinear(xc, yc, xd, yd, xb, yb))
        {
            value = false;
        }
        else
        {
            double t1 = triangle_area(xa, ya, xb, yb, xc, yc);
            double t2 = triangle_area(xa, ya, xb, yb, xd, yd);
            double t3 = triangle_area(xc, yc, xd, yd, xa, ya);
            double t4 = triangle_area(xc, yc, xd, yd, xb, yb);

            bool value1 = 0.0 < t1;
            bool value2 = 0.0 < t2;
            bool value3 = 0.0 < t3;
            bool value4 = 0.0 < t4;

            value = l4_xor(value1, value2) && l4_xor(value3, value4);
        }

        return value;
    }

    public static bool l4_xor(bool l1, bool l2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    L4_XOR returns the exclusive OR of two L4's.
        //
        //  Discussion:
        //
        //    An L4 is a logical value.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 May 2014
        //
        //  Author:
        //
        //   John Burkardt
        //
        //  Parameters:
        //
        //    Input, bool L1, L2, two values whose exclusive OR is needed.
        //
        //    Output, bool L4_XOR, the exclusive OR of L1 and L2.
        //
    {
        bool value1 = l1 && !l2;
        bool value2 = !l1 && l2;

        bool value = value1 || value2;

        return value;
    }

    public static double[] polygon_angles(int n, double[] v)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_ANGLES computes the interior angles of a polygon.
        //
        //  Discussion:
        //
        //    The vertices should be listed in counter clockwise order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 March 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of vertices of the polygon.
        //
        //    Input, double V[2*N], the vertices.
        //
        //    Output, double POLYGON_ANGLES[N], the angles of the polygon,
        //    in radians.
        //
    {
        int i;

        switch (n)
        {
            case < 1:
                return null;
        }

        double[] angle = new double[n];

        switch (n)
        {
            case <= 2:
            {
                for (i = 0; i < n; i++)
                {
                    angle[i] = 0.0;
                }

                return angle;
            }
        }

        for (i = 0; i < n; i++)
        {
            int im1 = typeMethods.i4_wrap(i - 1, 0, n - 1);
            int ip1 = typeMethods.i4_wrap(i + 1, 0, n - 1);

            angle[i] = angle_rad(v, v, v, p1Index: + im1 * 2, p2Index: + i * 2, p3Index: + ip1 * 2);
        }

        return angle;
    }

    public static double polygon_area(int n, double[] v)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_AREA computes the area of a polygon.
        //
        //  Discussion:
        //
        //    AREA = ABS ( 0.5 * SUM ( 1 <= I <= N ) X(I) * ( Y(I+1)-Y(I-1) ) )
        //    where Y[N] should be replaced by Y[0], and Y[N+1] by Y[1].
        //
        //    If the vertices are given in counter clockwise order, the area
        //    will be positive.  If the vertices are given in clockwise order,
        //    the area will be negative.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of vertices of the polygon.
        //
        //    Input, double V[2*N], the coordinates of the vertices.
        //
        //    Output, double POLYGON_AREA, the area of the polygon.
        //
    {
        int i;

        double area = 0.0;

        for (i = 0; i < n; i++)
        {
            int im1 = i - 1;
            int ip1 = i + 1;

            im1 = im1 switch
            {
                < 0 => n - 1,
                _ => im1
            };

            if (n <= ip1)
            {
                ip1 = 0;
            }

            area += v[0 + i * 2] * (v[1 + ip1 * 2] - v[1 + im1 * 2]);
        }

        area = 0.5 * area;

        return area;
    }

    public static double polygon_area_2(int n, double[] v)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_AREA_2 computes the area of a polygon.
        //
        //  Discussion:
        //
        //    The area is the sum of the areas of the triangles formed by
        //    node N with consecutive pairs of nodes.
        //
        //    If the vertices are given in counter clockwise order, the area
        //    will be positive.  If the vertices are given in clockwise order,
        //    the area will be negative.
        //
        //    Thanks to Martin Pineault for noticing that an earlier version
        //    of this routine would not correctly compute the area of a nonconvex
        //    polygon.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Adrian Bowyer, John Woodwark,
        //    A Programmer's Geometry,
        //    Butterworths, 1983.
        //
        //  Parameters:
        //
        //    Input, int N, the number of vertices of the polygon.
        //
        //    Input, double V[2*N], the coordinates of the vertices.
        //
        //    Output, double POLYGON_AREA_2, the area of the polygon.
        //
    {
        int i;

        double area = 0.0;

        for (i = 0; i < n - 2; i++)
        {
            double area_triangle = triangle_area(
                v[0 + i * 2], v[1 + i * 2],
                v[0 + (i + 1) * 2], v[1 + (i + 1) * 2],
                v[0 + (n - 1) * 2], v[1 + (n - 1) * 2]);

            area += area_triangle;
        }

        return area;
    }

    public static double[] polygon_centroid(int n, double[] v)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_CENTROID computes the centroid of a polygon.
        //
        //  Discussion:
        //
        //    Denoting the centroid coordinates by (CX,CY), then
        //
        //      CX = Integral ( Polygon interior ) x dx dy / Area ( Polygon )
        //      CY = Integral ( Polygon interior ) y dx dy / Area ( Polygon ).
        //
        //    Green's theorem states that
        //
        //      Integral ( Polygon boundary ) ( M dx + N dy ) =
        //      Integral ( Polygon interior ) ( dN/dx - dM/dy ) dx dy.
        //
        //    Using M = 0 and N = x^2/2, we get:
        //
        //      CX = 0.5 * Integral ( Polygon boundary ) x^2 dy,
        //
        //    which becomes
        //
        //      CX = 1/6 SUM ( 1 <= I <= N ) 
        //        ( X[I+1] + X[I] ) * ( X[I] * Y[I+1] - X[I+1] * Y[I] )
        //
        //    where, when I = N, the index "I+1" is replaced by 1.
        //
        //    A similar calculation gives us a formula for CY.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Gerard Bashein, Paul Detmer,
        //    Centroid of a Polygon,
        //    Graphics Gems IV, edited by Paul Heckbert,
        //    AP Professional, 1994, T385.G6974.
        //
        //  Parameters:
        //
        //    Input, int N, the number of sides of the polygon.
        //
        //    Input, double V[2*N], the coordinates of the vertices.
        //
        //    Output, double[] POLYGON_CENTROID[2], the coordinates of the centroid.
        //
    {
        int i;
        int j;

        double area = 0.0;
        double[] centroid = new double[2];

        for (j = 0; j < 2; j++)
        {
            centroid[j] = 0.0;
        }

        for (i = 0; i < n; i++)
        {
            int ip1;
            if (i < n - 1)
            {
                ip1 = i + 1;
            }
            else
            {
                ip1 = 0;
            }

            double temp = v[0 + i * 2] * v[1 + ip1 * 2] - v[0 + ip1 * 2] * v[1 + i * 2];

            area += temp;

            centroid[0] += (v[0 + ip1 * 2] + v[0 + i * 2]) * temp;
            centroid[1] += (v[1 + ip1 * 2] + v[1 + i * 2]) * temp;
        }

        area /= 2.0;

        for (j = 0; j < 2; j++)
        {
            centroid[j] /= 6.0 * area;
        }

        return centroid;
    }

    public static double[] polygon_centroid_2(int n, double[] v)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_CENTROID_2 computes the centroid of a polygon.
        //
        //  Method:
        //
        //    The centroid is the area-weighted sum of the centroids of
        //    disjoint triangles that make up the polygon.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Adrian Bowyer, John Woodwark,
        //    A Programmer's Geometry,
        //    Butterworths, 1983.
        //
        //  Parameters:
        //
        //    Input, int N, the number of vertices of the polygon.
        //
        //    Input, double V[2*N], the coordinates of the vertices. 
        //
        //    Output, double POLYGON_CENTROID_2[2], the coordinates of the centroid.
        //
    {
        int i;
        int j;

        double area = 0.0;
        double[] centroid = new double[2];

        for (j = 0; j < 2; j++)
        {
            centroid[j] = 0.0;
        }

        for (i = 0; i < n - 2; i++)
        {
            double area_triangle = triangle_area(
                v[0 + i * 2], v[1 + i * 2],
                v[0 + (i + 1) * 2], v[1 + (i + 1) * 2],
                v[0 + (n - 1) * 2], v[1 + (n - 1) * 2]);

            area += area_triangle;
            centroid[0] += area_triangle * (v[0 + i * 2] + v[0 + (i + 1) * 2] + v[0 + (n - 1) * 2]) / 3.0;
            centroid[1] += area_triangle * (v[1 + i * 2] + v[1 + (i + 1) * 2] + v[1 + (n - 1) * 2]) / 3.0;
        }

        for (j = 0; j < 2; j++)
        {
            centroid[j] /= area;
        }

        return centroid;
    }

    public static bool polygon_contains_point(int n, double[] v, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_CONTAINS_POINT finds if a point is inside a simple polygon .
        //
        //  Discussion:
        //
        //    A simple polygon is one whose boundary never crosses itself.
        //    The polygon does not need to be convex.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    M Shimrat,
        //    Position of Point Relative to Polygon,
        //    ACM Algorithm 112,
        //    Communications of the ACM,
        //    Volume 5, Number 8, page 434, August 1962.
        //
        //  Parameters:
        //
        //    Input, int N, the number of nodes or vertices in the polygon.
        //    N must be at least 3.
        //
        //    Input, double V[2*N], the coordinates of the vertices.
        //
        //    Input, double P[2], the coordinates of the point to be tested.
        //
        //    Output, bool POLYGON_CONTAINS_POINT, is TRUE if the point
        //    is inside the polygon or on its boundary, and FALSE otherwise.
        //
    {
        int i;

        bool value = false;

        for (i = 0; i < n; i++)
        {
            double x1 = v[0 + i * 2];
            double y1 = v[1 + i * 2];

            double x2;
            double y2;
            if (i < n - 1)
            {
                x2 = v[0 + (i + 1) * 2];
                y2 = v[1 + (i + 1) * 2];
            }
            else
            {
                x2 = v[0 + 0 * 2];
                y2 = v[1 + 0 * 2];
            }

            if (y1 < p[1] && p[1] <= y2 ||
                p[1] <= y1 && y2 < p[1])
            {
                value = (p[0] - x1 - (p[1] - y1) * (x2 - x1) / (y2 - y1)) switch
                {
                    < 0 => !value,
                    _ => value
                };
            }
        }

        return value;
    }

    public static bool polygon_contains_point_2(int n, double[] v, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_CONTAINS_POINT_2 finds if a point is inside a convex polygon.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of nodes or vertices in the polygon.
        //    N must be at least 3.
        //
        //    Input, double V[2*N], the coordinates of the vertices.
        //
        //    Input, double P[2], the coordinates of the point to be tested.
        //
        //    Output, bool POLYGON_CONTAINS_POINT_2, is TRUE if the point
        //    is inside the polygon or on its boundary, and FALSE otherwise.
        //
    {
        int i;
        double[] t = new double[2 * 3];
        //
        //  A point is inside a convex polygon if and only if it is inside
        //  one of the triangles formed by the first vertex and any two consecutive
        //  vertices.
        //
        t[0 + 0 * 2] = v[0 + 0 * 2];
        t[1 + 0 * 2] = v[1 + 0 * 2];

        for (i = 1; i < n - 1; i++)
        {
            t[0 + 1 * 2] = v[0 + i * 2];
            t[1 + 1 * 2] = v[1 + i * 2];
            t[0 + 2 * 2] = v[0 + (i + 1) * 2];
            t[1 + 2 * 2] = v[1 + (i + 1) * 2];

            if (triangle_contains_point_1(t, p))
            {
                return true;
            }
        }

        return false;
    }

    public static double polygon_diameter(int n, double[] v)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_DIAMETER computes the diameter of a polygon.
        //
        //  Discussion:
        //
        //    The diameter of a polygon is the maximum distance between any
        //    two points on the polygon.  It is guaranteed that this maximum
        //    distance occurs between two vertices of the polygon.  It is
        //    sufficient to check the distance between all pairs of vertices.
        //    This is an N^2 algorithm.  There is an algorithm by Shamos which
        //    can compute this quantity in order N time instead.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of vertices of the polygon.
        //
        //    Input, double V[2*N], the coordinates of the vertices.
        //
        //    Output, double POLYGON_DIAMETER, the diameter of the polygon.
        //
    {
        int i;

        double diameter = 0.0;

        for (i = 0; i < n; i++)
        {
            int j;
            for (j = i + 1; j < n; j++)
            {
                double t = Math.Sqrt((v[0 + i * 2] - v[0 + j * 2]) * (v[0 + i * 2] - v[0 + j * 2])
                                     + (v[1 + i * 2] - v[1 + j * 2]) * (v[1 + i * 2] - v[1 + j * 2]));
                if (diameter < t)
                {
                    diameter = t;
                }
            }
        }

        return diameter;
    }

    public static double[] polygon_expand(int n, double[] v, double h)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_EXPAND expands a polygon.
        //
        //  Discussion:
        //
        //    This routine simple moves each vertex of the polygon outwards
        //    in such a way that the sides of the polygon advance by H.
        //
        //    This approach should always work if the polygon is convex, or
        //    star-shaped.  But for general polygons, it is possible
        //    that this procedure, for large enough H, will create a polygon
        //    whose sides intersect.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //   21 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of sides of the polygon.
        //
        //    Input, double V[2*N], the coordinates of the vertices.
        //
        //    Input, double H, the expansion amount.
        //
        //    Output, double POLYGON_EXPAND[2*N], the "expanded" coordinates.
        //
    {
        int j;

        double[] w = new double[2 * n];
        //
        //  Consider each angle, formed by the nodes P(I-1), P(I), P(I+1).
        //
        for (j = 0; j < n; j++)
        {
            int jm1 = typeMethods.i4_wrap(j - 1, 0, n - 1);
            int jp1 = typeMethods.i4_wrap(j + 1, 0, n - 1);
            //
            //        P1
            //        /
            //       /   P4
            //      /  .
            //     / .
            //    P2--------->P3
            //
            double[] p4 = angle_half(v, v, v, p1Index: + jm1 * 2, p2Index: + j * 2, p3Index: + jp1 * 2);
            //
            //  Compute the value of the half angle.
            //
            double angle = angle_rad(v, v, p4, p1Index: + jm1 * 2, p2Index: + j * 2);
            //
            //  The stepsize along the ray must be adjusted so that the sides
            //  move out by H.
            //
            double h2 = h / Math.Sin(angle);

            int i;
            for (i = 0; i < 2; i++)
            {
                w[i + j * 2] = v[i + j * 2] - h2 * (p4[i] - v[i + j * 2]);
            }

        }

        return w;
    }

    public static void polygon_inrad_data(int n, double radin, ref double area, ref double radout,
            ref double side)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_INRAD_DATA determines polygonal data from its inner radius.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of sides of the polygon.
        //    N must be at least 3.
        //
        //    Input, double RADIN, the inner radius of the polygon, that is,
        //    the radius of the largest circle that can be inscribed within
        //    the polygon.
        //
        //    Output, double &AREA, the area of the regular polygon.
        //
        //    Output, double &RADOUT, the outer radius of the polygon, that is,
        //    the radius of the smallest circle that can be described about
        //    the polygon.
        //
        //    Output, double &SIDE, the length of one side of the polygon.
        //
    {
        switch (n)
        {
            case < 3:
                Console.WriteLine("");
                Console.WriteLine("POLYGON_INRAD_DATA - Fatal error!");
                Console.WriteLine("  Input value of N must be at least 3,");
                Console.WriteLine("  but your input value was N = " + n + "");
                return;
        }

        double angle = Math.PI / n;
        area = n * radin * radin * Math.Tan(angle);
        side = 2.0 * radin * Math.Tan(angle);
        radout = 0.5 * side / Math.Sin(angle);

    }

    public static double polygon_integral_1(int n, double[] v)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_INTEGRAL_1 integrates the function 1 over a polygon.
        //
        //  Discussion:
        //
        //    INTEGRAL = 0.5 * SUM ( 1 <= I <= N ) (X(I)+X(I-1)) * (Y(I)-Y(I-1))
        //
        //    where X[N] and Y[N] should be replaced by X[0] and Y[0].
        //
        //    The integral of 1 over a polygon is the area of the polygon.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    SF Bockman,
        //    Generalizing the Formula for Areas of Polygons to Moments,
        //    American Mathematical Society Monthly,
        //    1989, pages 131-132.
        //
        //  Parameters:
        //
        //    Input, int N, the number of vertices of the polygon.
        //    N should be at least 3 for a nonzero result.
        //
        //    Input, double V[2*N], the coordinates of the vertices
        //    of the polygon.  These vertices should be given in
        //    counter clockwise order.
        //
        //    Output, double POLYGON_INTEGRAL_1, the value of the integral.
        //
    {
        int i;

        double result = 0.0;

        switch (n)
        {
            case < 3:
                Console.WriteLine("");
                Console.WriteLine("POLYGON_INTEGRAL_1 - Fatal error!");
                Console.WriteLine("  The number of vertices must be at least 3.");
                Console.WriteLine("  The input value of N = " + n + "");
                return 1;
        }

        for (i = 0; i < n; i++)
        {
            int im1 = i switch
            {
                0 => n - 1,
                _ => i - 1
            };

            result += 0.5 * (v[0 + i * 2] + v[0 + im1 * 2])
                          * (v[1 + i * 2] - v[1 + im1 * 2]);
        }

        return result;
    }

    public static double polygon_integral_x(int n, double[] v)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_INTEGRAL_X integrates the function X over a polygon.
        //
        //  Discussion:
        //
        //    INTEGRAL = (1/6) * SUM ( I = 1 to N )
        //      ( X[I]^2 + X[I] * X[I-1] + X[I-1]^2 ) * ( Y[I] - Y[I-1] )
        //
        //    where X[N] and Y[N] should be replaced by X[0] and Y[0].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    SF Bockman,
        //    Generalizing the Formula for Areas of Polygons to Moments,
        //    American Mathematical Society Monthly,
        //    1989, pages 131-132.
        //
        //  Parameters:
        //
        //    Input, int N, the number of vertices of the polygon.
        //    N should be at least 3 for a nonzero result.
        //
        //    Input, double V[2*N], the coordinates of the vertices.
        //
        //    Output, double POLYGON_INTEGRAL_X, the value of the integral.
        //
    {
        int i;

        double result = 0.0;

        switch (n)
        {
            case < 3:
                Console.WriteLine("");
                Console.WriteLine("POLYGON_INTEGRAL_X - Fatal error!");
                Console.WriteLine("  The number of vertices must be at least 3.");
                Console.WriteLine("  The input value of N = " + n + "");
                return 1;
        }

        for (i = 0; i < n; i++)
        {
            int im1 = i switch
            {
                0 => n - 1,
                _ => i - 1
            };

            result += (v[0 + i * 2] * v[0 + i * 2]
                       + v[0 + i * 2] * v[0 + im1 * 2]
                       + v[0 + im1 * 2] * v[0 + im1 * 2])
                      * (v[1 + i * 2] - v[1 + im1 * 2]);
        }

        result /= 6.0;

        return result;
    }

    public static double polygon_integral_xx(int n, double[] v)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_INTEGRAL_XX integrates the function X*X over a polygon.
        //
        //  Discussion:
        //
        //    INTEGRAL = (1/12) * SUM ( I = 1 to N )
        //      ( X[I]^3 + X[I]^2 * X[I-1] + X[I] * X[I-1]^2 + X[I-1]^3 )
        //      * ( Y[I] - Y[I-1] )
        //
        //    where X[N] and Y[N] should be replaced by X[0] and Y[0].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    SF Bockman,
        //    Generalizing the Formula for Areas of Polygons to Moments,
        //    American Mathematical Society Monthly,
        //    1989, pages 131-132.
        //
        //  Parameters:
        //
        //    Input, int N, the number of vertices of the polygon.
        //    N should be at least 3 for a nonzero result.
        //
        //    Input, double V[2*N], the coordinates of the vertices.
        //
        //    Output, double POLYGON_INTEGRAL_XX, the value of the integral.
        //
    {
        int i;

        double result = 0.0;

        switch (n)
        {
            case < 3:
                Console.WriteLine("");
                Console.WriteLine("POLYGON_INTEGRAL_XX - Fatal error!");
                Console.WriteLine("  The number of vertices must be at least 3.");
                Console.WriteLine("  The input value of N = " + n + "");
                return 1;
        }

        for (i = 0; i < n; i++)
        {
            int im1 = i switch
            {
                0 => n - 1,
                _ => i - 1
            };

            result += (
                v[0 + i * 2] * v[0 + i * 2] * v[0 + i * 2]
                + v[0 + i * 2] * v[0 + i * 2] * v[0 + im1 * 2]
                + v[0 + i * 2] * v[0 + im1 * 2] * v[0 + im1 * 2]
                + v[0 + im1 * 2] * v[0 + im1 * 2] * v[0 + im1 * 2]) * (v[1 + i * 2] - v[1 + im1 * 2]);
        }

        result /= 12.0;

        return result;
    }

    public static double polygon_integral_xy(int n, double[] v)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_INTEGRAL_XY integrates the function X*Y over a polygon.
        //
        //  Discussion:
        //
        //    INTEGRAL = (1/24) * SUM (I=1 to N)
        //      ( Y[I] * 
        //        ( 3 * X[I]^2 + 2 * X[I] * X[I-1] + X[I-1]^2 )
        //      + Y[I-1] *
        //        ( X[I]^2 + 2 * X[I] * X[I-1] + 3 * X[I-1]^2 ) 
        //      ) * ( Y[I] - Y[I-1] )
        //
        //    where X[N] and Y[N] should be replaced by X[0] and Y[0].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    SF Bockman,
        //    Generalizing the Formula for Areas of Polygons to Moments,
        //    American Mathematical Society Monthly,
        //    1989, pages 131-132.
        //
        //  Parameters:
        //
        //    Input, int N, the number of vertices of the polygon.
        //    N should be at least 3 for a nonzero result.
        //
        //    Input, double V[2*N], the coordinates of the vertices.
        //
        //    Output, double POLYGON_INTEGRAL_XY, the value of the integral.
        //
    {
        int i;

        double result = 0.0;

        switch (n)
        {
            case < 3:
                Console.WriteLine("");
                Console.WriteLine("POLYGON_INTEGRAL_XY - Fatal error!");
                Console.WriteLine("  The number of vertices must be at least 3.");
                Console.WriteLine("  The input value of N = " + n + "");
                return 1;
        }

        for (i = 0; i < n; i++)
        {
            int im1 = i switch
            {
                0 => n - 1,
                _ => i - 1
            };

            result += (
                v[1 + i * 2] * (3.0 * v[0 + i * 2] * v[0 + i * 2]
                                + 2.0 * v[0 + i * 2] * v[0 + im1 * 2]
                                + v[0 + im1 * 2] * v[0 + im1 * 2])
                + v[1 + im1 * 2] * (v[0 + i * 2] * v[0 + i * 2]
                                    + 2.0 * v[0 + i * 2] * v[0 + im1 * 2]
                                    + 3.0 * v[0 + im1 * 2] * v[0 + im1 * 2])
            ) * (v[1 + i * 2] - v[1 + im1 * 2]);
        }

        result /= 24.0;

        return result;
    }

    public static double polygon_integral_y(int n, double[] v)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_Y integrates the function Y over a polygon.
        //
        //  Discussion:
        //
        //    INTEGRAL = (1/6) * SUM ( I = 1 to N )
        //      - ( Y[I]^2 + Y[I] * Y[I-1] + Y[I-1]^2 ) * ( X[I] - X[I-1] )
        //
        //    where X[N] and Y[N] should be replaced by X[0] and Y[0].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    SF Bockman,
        //    Generalizing the Formula for Areas of Polygons to Moments,
        //    American Mathematical Society Monthly,
        //    1989, pages 131-132.
        //
        //  Parameters:
        //
        //    Input, int N, the number of vertices of the polygon.
        //    N should be at least 3 for a nonzero result.
        //
        //    Input, double V[2*N], the coordinates of the vertices.
        //
        //    Output, double POLYGON_INTEGRAL_Y, the value of the integral.
        //
    {
        int i;

        double result = 0.0;

        switch (n)
        {
            case < 3:
                Console.WriteLine("");
                Console.WriteLine("POLYGON_INTEGRAL_Y - Fatal error!");
                Console.WriteLine("  The number of vertices must be at least 3.");
                Console.WriteLine("  The input value of N = " + n + "");
                return 1;
        }

        for (i = 0; i < n; i++)
        {
            int im1 = i switch
            {
                0 => n - 1,
                _ => i - 1
            };

            result -= (v[1 + i * 2] * v[1 + i * 2]
                       + v[1 + i * 2] * v[1 + im1 * 2]
                       + v[1 + im1 * 2] * v[1 + im1 * 2])
                      * (v[0 + i * 2] - v[0 + im1 * 2]);
        }

        result /= 6.0;

        return result;
    }

    public static double polygon_integral_yy(int n, double[] v)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_INTEGRAL_YY integrates the function Y*Y over a polygon.
        //
        //  Discussion:
        //
        //    INTEGRAL = (1/12) * SUM ( I = 1 to N )
        //      - ( Y[I]^3 + Y[I]^2 * Y[I-1] + Y[I] * Y[I-1]^2 + Y[I-1]^3 )
        //      * ( X[I] - X[I-1] )
        //
        //    where X[N] and Y[N] should be replaced by X[0] and Y[0].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    SF Bockman,
        //    Generalizing the Formula for Areas of Polygons to Moments,
        //    American Mathematical Society Monthly,
        //    1989, pages 131-132.
        //
        //  Parameters:
        //
        //    Input, int N, the number of vertices of the polygon.
        //    N should be at least 3 for a nonzero result.
        //
        //    Input, double V[2*N], the coordinates of the vertices.
        //
        //    Output, double POLYGON_INTEGRAL_YY, the value of the integral.
        //
        //
    {
        int i;

        double result = 0.0;

        switch (n)
        {
            case < 3:
                Console.WriteLine("");
                Console.WriteLine("POLYGON_INTEGRAL_YY - Fatal error!");
                Console.WriteLine("  The number of vertices must be at least 3.");
                Console.WriteLine("  The input value of N = " + n + "");
                return 1;
        }

        for (i = 0; i < n; i++)
        {
            int im1 = i switch
            {
                0 => n - 1,
                _ => i - 1
            };

            result -= (v[1 + i * 2] * v[1 + i * 2] * v[1 + i * 2]
                       + v[1 + i * 2] * v[1 + i * 2] * v[1 + im1 * 2]
                       + v[1 + i * 2] * v[1 + im1 * 2] * v[1 + im1 * 2]
                       + v[1 + im1 * 2] * v[1 + im1 * 2] * v[1 + im1 * 2])
                      * (v[0 + i * 2] - v[0 + im1 * 2]);
        }

        result /= 12.0;

        return result;
    }

    public static int polygon_is_convex(int n, double[] v)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_IS_CONVEX determines whether a polygon is convex in 2D.
        //
        //  Discussion:
        //
        //    If the polygon has less than 3 distinct vertices, it is
        //    classified as convex degenerate.
        //
        //    If the polygon "goes around" more than once, it is classified
        //    as NOT convex.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Peter Schorn, Frederick Fisher,
        //    Testing the Convexity of a Polygon,
        //    Graphics Gems IV,
        //    edited by Paul Heckbert,
        //    AP Professsional, 1994, T385.G6974.
        //
        //  Parameters
        //
        //    Input, int N, the number of vertices.
        //
        //    Input/output, double V[2*N], the coordinates of the vertices of the
        //    polygon.  On output, duplicate consecutive points have been deleted,
        //    and the vertices have been reordered so that the lexicographically
        //    least point comes first.
        //
        //    Output, int POLYGON_IS_CONVEX:
        //    -1, the polygon is not convex;
        //     0, the polygon has less than 3 vertices; it is "degenerately" convex;
        //     1, the polygon is convex and counter clockwise;
        //     2, the polygon is convex and clockwise.
        //
    {
        const int NOT_CONVEX = -1;
        const int DEGENERATE_CONVEX = 0;
        const int CONVEX_CCW = 1;
        const int CONVEX_CW = 2;

        int i;
        const double TOL = 1.0;
        int value = 0;

        double exterior_total = 0.0;
        switch (n)
        {
            //
            //  If there are not at least 3 distinct vertices, we are done.
            //
            case < 3:
                return DEGENERATE_CONVEX;
        }

        double sense = 0.0;
        //
        //  Consider each polygonal vertex I.
        //
        for (i = 0; i < n; i++)
        {
            int ip1 = typeMethods.i4_wrap(i + 1, 0, n - 1);
            int ip2 = typeMethods.i4_wrap(i + 2, 0, n - 1);

            double dot = (v[0 + ip2 * 2] - v[0 + ip1 * 2]) * (v[0 + i * 2] - v[0 + ip1 * 2])
                         + (v[1 + ip2 * 2] - v[1 + ip1 * 2]) * (v[1 + i * 2] - v[1 + ip1 * 2]);

            double cross = (v[0 + ip2 * 2] - v[0 + ip1 * 2]) * (v[1 + i * 2] - v[1 + ip1 * 2])
                           - (v[0 + i * 2] - v[0 + ip1 * 2]) * (v[1 + ip2 * 2] - v[1 + ip1 * 2]);

            double angle = Math.Atan2(cross, dot);
            switch (sense)
            {
                //
                //  See if the turn defined by this vertex is our first indication of
                //  the "sense" of the polygon, or if it disagrees with the previously
                //  defined sense.
                //
                case 0.0 when angle < 0.0:
                    sense = -1.0;
                    break;
                case 0.0:
                {
                    sense = angle switch
                    {
                        > 0.0 => +1.0,
                        _ => sense
                    };

                    break;
                }
                case 1.0 when angle < 0.0:
                    return NOT_CONVEX;
                case -1.0 when 0.0 < angle:
                    return NOT_CONVEX;
            }

            //
            //  If the exterior total is greater than 360, then the polygon is
            //  going around again.
            //
            angle = Math.Atan2(-cross, -dot);

            exterior_total += angle;

            if (360.0 + TOL < typeMethods.r8_degrees(Math.Abs(exterior_total)))
            {
                return NOT_CONVEX;
            }

        }

        value = sense switch
        {
            +1.0 => CONVEX_CCW,
            -1.0 => CONVEX_CW,
            _ => value
        };

        return value;

    }

    public static double polygon_lattice_area(int i, int b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_LATTICE_AREA computes the area of a lattice polygon.
        //
        //  Discussion:
        //
        //    We define a lattice to be the plane, in which the points
        //    whose coordinates are both integers are given a special
        //    status as "lattice points".
        //
        //    A lattice polygon is a polygon whose vertices are lattice points.
        //
        //    The area of a lattice polygon can be computed by Pick's Theorem:
        //
        //      Area = I + B / 2 - 1
        //
        //    where
        //
        //      I = the number of lattice points contained strictly inside the polygon;
        //
        //      B = the number of lattice points that lie exactly on the boundary.
        //  
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Branko Gruenbaum, G C Shephard,
        //    Pick's Theorem,
        //    The American Mathematical Monthly,
        //    Volume 100, 1993, pages 150-161.
        //
        //  Parameters:
        //
        //    Input, int I, the number of interior lattice points.
        //
        //    Input, int B, the number of boundary lattice points.
        //
        //    Output, double POLYGON_LATTICE_AREA, the area of the lattice polygon.
        //
    {
        double value = 0;

        value = i + b / 2.0 - 1.0;

        return value;
    }

    public static void polygon_outrad_data(int n, double radout, ref double area, ref double radin,
            ref double side)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_OUTRAD_DATA determines polygonal data from its outer radius.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of sides of the polygon.
        //    N must be at least 3.
        //
        //    Input, double RADOUT, the outer radius of the polygon, that is,
        //    the radius of the smallest circle that can be described
        //    around the polygon.
        //
        //    Output, double &AREA, the area of the regular polygon.
        //
        //    Output, double &RADIN, the inner radius of the polygon, that is,
        //    the radius of the largest circle that can be inscribed
        //    within the polygon.
        //
        //    Output, double &SIDE, the length of one side of the polygon.
        //
    {
        switch (n)
        {
            case < 3:
                Console.WriteLine("");
                Console.WriteLine("POLYGON_OUTRAD_DATA - Fatal error!");
                Console.WriteLine("  Input value of N must be at least 3,");
                Console.WriteLine("  but your input value was N = " + n + "");
                return;
        }

        double angle = Math.PI / n;
        area = 0.5 * n * radout * radout * Math.Sin(2.0 * angle);
        side = 2.0 * radout * Math.Sin(angle);
        radin = 0.5 * side / Math.Tan(angle);

    }

    public static double polygon_perimeter(int n, double[] v)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_PERIMETER computes the perimeter of a polygon.
        //
        //  Discussion:
        //
        //    The perimeter is simply the sum of the side lengths.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of vertices of the polygon.
        //
        //    Input, double V(2,N), the vertices.
        //
        //    Output, double POLYGON_PERIMETER, the perimeter.
        //
    {
        int i;

        double value = 0.0;

        int im1 = n - 1;

        for (i = 0; i < n; i++)
        {
            double l = Math.Sqrt(Math.Pow(v[0 + im1 * 2] - v[0 + i * 2], 2)
                                 + Math.Pow(v[1 + im1 * 2] - v[1 + i * 2], 2));
            value += l;
            im1 = i;
        }

        return value;
    }

    public static double polygon_perimeter_quad(int n, double[] v, double hmax,
            Func<double, double, double> f)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_PERIMETER_QUAD estimates an integral over the perimeter of a polygon.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of vertices of the polygon.
        //
        //    Input, double V(2,N), the vertices.
        //
        //    Input, double HMAX, the maximum length of a quadrature interval.
        //
        //    Input, double F ( double X, double Y ), a function whose integral 
        //    over the perimeter is desired.
        //
        //    Output, double POLYGON_PERIMETER_QUAD, the estimated integral.
        //
    {
        int i;

        double value = 0.0;

        for (i = 0; i < n; i++)
        {
            int ip1 = typeMethods.i4_wrap(i + 1, 0, n - 1);
            double l = Math.Sqrt(Math.Pow(v[0 + ip1 * 2] - v[0 + i * 2], 2) + Math.Pow(v[1 + ip1 * 2] - v[1 + i * 2], 2));
            int m = (int)Math.Ceiling(l / hmax);
            double dxy = l / m;

            int j;
            for (j = 1; j <= 2 * m - 1; j += 2)
            {
                double x = ((2 * m - j) * v[0 + i * 2]
                            + j * v[0 + ip1 * 2])
                           / (2 * m);
                double y = ((2 * m - j) * v[1 + i * 2]
                            + j * v[1 + ip1 * 2])
                           / (2 * m);
                double fxy = f(x, y);
                value += fxy * dxy;
            }
        }

        return value;
    }

    public static double polygon_point_dist(int n, double[] v, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_POINT_DIST: distance ( polygon, point ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of vertices.
        //
        //    Input, double V[2*N], the triangle vertices.
        //
        //    Input, double P[2], the point to be checked.
        //
        //    Output, double POLYGON_POINT_DIST, the distance from the point to the
        //    polygon.
        //
    {
        int j;
        //
        //  Find the distance to each of the line segments.
        //
        double dist = typeMethods.r8_huge();

        for (j = 0; j < n; j++)
        {
            int jp1 = typeMethods.i4_wrap(j + 1, 0, n - 1);

            double dist2 = segment_point_dist(v, v, p, p1Index: + j * 2, p2Index: + jp1 * 2);

            if (dist2 < dist)
            {
                dist = dist2;
            }
        }

        return dist;
    }

    public static double[] polygon_point_near(int n, double[] v, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_POINT_NEAR computes the nearest point on a polygon.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 February 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of vertices.
        //
        //    Input, double V[2*N], the polygon vertices.
        //
        //    Input, double P[2], the point whose nearest polygon point
        //    is to be determined.
        //
        //    Output, double POLYGON_POINT_NEAR[2], the nearest point to P.
        //
    {
        int j;
        //
        //  Find the distance to each of the line segments that make up the edges
        //  of the polygon.
        //
        double dist = typeMethods.r8_huge();

        double[] pn = new double[2];
        pn[0] = 0.0;
        pn[1] = 0.0;

        for (j = 0; j < n; j++)
        {
            int jp1 = typeMethods.i4_wrap(j + 1, 0, n - 1);

            double[] pn2 = segment_point_near(v, v, p, p1Index: + j * 2, p2Index: + jp1 * 2);

            double dist2 = Math.Sqrt(Math.Pow(pn2[0] - p[0], 2) + Math.Pow(pn2[1] - p[1], 2));

            if (dist2 < dist)
            {
                dist = dist2;
                pn[0] = pn2[0];
                pn[1] = pn2[1];
            }
        }

        return pn;
    }

    public static double[] polygon_sample(int nv, double[] v, int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_SAMPLE uniformly samples a polygon.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NV, the number of vertices.
        //
        //    Input, double V[2*NV], the vertices of the polygon, listed in
        //    counterclockwise order.
        //
        //    Input, int N, the number of points to create.
        //
        //    Input/output, int *SEED, a seed for the random
        //    number generator.
        //
        //    Output, double POLYGON_SAMPLE[2*N], the points.
        //
    {
        int i;
        int j;
        //
        //  Triangulate the polygon.
        //
        double[] x = new double[nv];
        double[] y = new double[nv];
        for (i = 0; i < nv; i++)
        {
            x[i] = v[0 + i * 2];
            y[i] = v[1 + i * 2];
        }

        int[] triangles = Triangulate.polygon_triangulate(nv, x, y);
        //
        //  Determine the areas of each triangle.
        //
        double[] area_triangle = new double[nv - 2];

        for (i = 0; i < nv - 2; i++)
        {
            area_triangle[i] = triangle_area(
                v[0 + triangles[0 + i * 3] * 2], v[1 + triangles[0 + i * 3] * 2],
                v[0 + triangles[1 + i * 3] * 2], v[1 + triangles[1 + i * 3] * 2],
                v[0 + triangles[2 + i * 3] * 2], v[1 + triangles[2 + i * 3] * 2]);
        }

        //
        //  Normalize the areas.
        //
        double area_polygon = typeMethods.r8vec_sum(nv - 2, area_triangle);

        double[] area_relative = new double[nv - 2];

        for (i = 0; i < nv - 2; i++)
        {
            area_relative[i] = area_triangle[i] / area_polygon;
        }

        //
        //  Replace each area by the sum of itself and all previous ones.
        //
        double[] area_cumulative = new double[nv - 2];

        area_cumulative[0] = area_relative[0];
        for (i = 1; i < nv - 2; i++)
        {
            area_cumulative[i] = area_relative[i] + area_cumulative[i - 1];
        }

        double[] s = new double[2 * n];

        for (j = 0; j < n; j++)
        {
            //
            //  Choose triangle I at random, based on areas.
            //
            double area_percent = UniformRNG.r8_uniform_01(ref seed);

            int k;
            for (k = 0; k < nv - 2; k++)
            {
                i = k;

                if (area_percent <= area_cumulative[k])
                {
                    break;
                }
            }

            //
            //  Now choose a point at random in triangle I.
            //
            double[] r = UniformRNG.r8vec_uniform_01_new(2, ref seed);

            switch (r[0] + r[1])
            {
                case > 1.0:
                    r[0] = 1.0 - r[0];
                    r[1] = 1.0 - r[1];
                    break;
            }

            s[0 + j * 2] = (1.0 - r[0] - r[1]) * v[0 + triangles[0 + i * 3] * 2]
                           + r[0] * v[0 + triangles[1 + i * 3] * 2]
                           + r[1] * v[0 + triangles[2 + i * 3] * 2];

            s[1 + j * 2] = (1.0 - r[0] - r[1]) * v[1 + triangles[0 + i * 3] * 2]
                           + r[0] * v[1 + triangles[1 + i * 3] * 2]
                           + r[1] * v[1 + triangles[2 + i * 3] * 2];
        }

        return s;
    }

    public static void polygon_side_data(int n, double side, ref double area, ref double radin,
            ref double radout)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_SIDE_DATA determines polygonal data from its side length.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of sides of the polygon.
        //    N must be at least 3.
        //
        //    Input, double SIDE, the length of one side of the polygon.
        //
        //    Output, double[] AREA, the area of the regular polygon.
        //
        //    Output, double[] RADIN, the inner radius of the polygon, that is,
        //    the radius of the largest circle that can be inscribed within
        //    the polygon.
        //
        //    Output, double[] RADOUT, the outer radius of the polygon, that is,
        //    the radius of the smallest circle that can be described about
        //    the polygon.
        //
    {
        switch (n)
        {
            case < 3:
                Console.WriteLine("");
                Console.WriteLine("POLYGON_SIDE_DATA - Fatal error!");
                Console.WriteLine("  Input value of N must be at least 3,");
                Console.WriteLine("  but your input value was N = " + n + "");
                return;
        }

        double angle = Math.PI / n;
        area = 0.25 * n * side * side / Math.Tan(angle);
        radin = 0.5 * side / Math.Tan(angle);
        radout = 0.5 * side / Math.Sin(angle);

    }


    public static double segment_point_dist(double[] p1, double[] p2, double[] p, int p1Index = 0, int p2Index = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SEGMENT_POINT_DIST: distance ( line segment, point ).
        //
        //  Discussion:
        //
        //    A line segment is the finite portion of a line that lies between
        //    two points.
        //
        //    The nearest point will satisfy the condition
        //
        //      PN = (1-T) * P1 + T * P2.
        //
        //    T will always be between 0 and 1.
        //
        //    Thanks to Kirill Speransky for pointing out that a previous version
        //    of this routine was incorrect, 02 May 2006.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[2], P2[2], the endpoints of the line segment.
        //
        //    Input, double P[2], the point whose nearest neighbor on the line
        //    segment is to be determined.
        //
        //    Output, double SEGMENT_POINT_DIST, the distance from the point 
        //    to the line segment.
        //
    {
        int i;
        double t;
        double[] pn = new double[2];
        //
        //  If the line segment is actually a point, then the answer is easy.
        //
        if (Math.Abs(p1[p1Index + 0] - p2[p2Index + 0]) <= double.Epsilon && Math.Abs(p1[p1Index + 1] - p2[p2Index + 1]) <= double.Epsilon)
        {
            t = 0.0;
        }
        else
        {
            double bot = 0.0;
            for (i = 0; i < 2; i++)
            {
                bot += Math.Pow(p2[p2Index + i] - p1[p1Index + i], 2);
            }

            t = 0.0;
            for (i = 0; i < 2; i++)
            {
                t += (p[i] - p1[p1Index + i]) * (p2[p2Index + i] - p1[p1Index + i]);
            }

            t /= bot;

            t = t switch
            {
                > 1.0 => 1.0,
                _ => t switch
                {
                    < 0.0 => 0.0,
                    _ => t
                }
            };
        }

        for (i = 0; i < 2; i++)
        {
            pn[i] = p1[p1Index + i] + t * (p2[p2Index + i] - p1[p1Index + i]);
        }

        double dist = 0.0;
        for (i = 0; i < 2; i++)
        {
            dist += Math.Pow(p[i] - pn[i], 2);
        }

        dist = Math.Sqrt(dist);

        return dist;
    }

    public static double[] segment_point_near(double[] p1, double[] p2, double[] p, int p1Index = 0, int p2Index = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SEGMENT_POINT_NEAR finds the point on a line segment nearest a point.
        //
        //  Discussion:
        //
        //    A line segment is the finite portion of a line that lies between
        //    two points.
        //
        //    The nearest point will satisfy the condition:
        //
        //      PN = (1-T) * P1 + T * P2.
        //
        //    and T will always be between 0 and 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[2], P2[2], the two endpoints of the line segment.
        //
        //    Input, double P[2], the point whose nearest neighbor
        //    on the line segment is to be determined.
        //
        //    Output, double SEGMENT_POINT_NEAR[2], the point on the line segment 
        //    which is nearest P.
        //
    {
        int i;
        double t;
        //
        //  If the line segment is actually a point, then the answer is easy.
        //
        if (Math.Abs(p1[p1Index + 0] - p2[p2Index + 0]) <= double.Epsilon && Math.Abs(p1[p1Index + 1] - p2[p2Index + 1]) <= double.Epsilon)
        {
            t = 0.0;
        }
        else
        {
            double bot = 0.0;
            for (i = 0; i < 2; i++)
            {
                bot += Math.Pow(p2[p2Index + i] - p1[p1Index + i], 2);
            }

            t = 0.0;
            for (i = 0; i < 2; i++)
            {
                t += (p[i] - p1[p1Index + i]) * (p2[p2Index + i] - p1[p1Index + i]);
            }

            t /= bot;

            t = t switch
            {
                > 1.0 => 1.0,
                _ => t switch
                {
                    < 0.0 => 0.0,
                    _ => t
                }
            };
        }

        double[] pn = new double[2];

        for (i = 0; i < 2; i++)
        {
            pn[i] = p1[p1Index + i] + t * (p2[p2Index + i] - p1[p1Index + i]);
        }

        return pn;
    }

    public static double triangle_area(double xa, double ya, double xb, double yb, double xc,
            double yc)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_AREA computes the signed area of a triangle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double XA, YA, XB, YB, XC, YC, the coordinates of
        //    the vertices of the triangle, given in counterclockwise order.
        //
        //    Output, double TRIANGLE_AREA, the signed area of the triangle.
        //
    {
        double value = 0;

        value = 0.5 * (
            (xb - xa) * (yc - ya)
            - (xc - xa) * (yb - ya));

        return value;
    }

    public static double[] triangle_barycentric(double[] t, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_BARYCENTRIC finds the barycentric coordinates of a point.
        //
        //  Discussion:
        //
        //    The barycentric coordinate of point X related to vertex A can be
        //    interpreted as the ratio of the area of the triangle with 
        //    vertex A replaced by vertex X to the area of the original 
        //    triangle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T[2*3], the vertices of the triangle.
        //
        //    Input, double P[2], the point to be checked.
        //
        //    Output, double C[3], the barycentric coordinates of the point with respect
        //    to the triangle.
        //
    {
        const int N = 2;
        const int RHS_NUM = 1;

        double[] a = new double[N * (N + RHS_NUM)];
        //
        //  Set up the linear system
        //
        //    ( X2-X1  X3-X1 ) C1  = X-X1
        //    ( Y2-Y1  Y3-Y1 ) C2    Y-Y1
        //
        //  which is satisfied by the barycentric coordinates.
        //
        a[0 + 0 * N] = t[0 + 1 * 2] - t[0 + 0 * 2];
        a[1 + 0 * N] = t[1 + 1 * 2] - t[1 + 0 * 2];

        a[0 + 1 * N] = t[0 + 2 * 2] - t[0 + 0 * 2];
        a[1 + 1 * N] = t[1 + 2 * 2] - t[1 + 0 * 2];

        a[0 + 2 * N] = p[0] - t[0 + 0 * 2];
        a[1 + 2 * N] = p[1] - t[1 + 0 * 2];
        //
        //  Solve the linear system.
        //
        int info = typeMethods.r8mat_solve(N, RHS_NUM, ref a);

        if (info != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_BARYCENTRIC - Fatal error!");
            Console.WriteLine("  The linear system is singular.");
            Console.WriteLine("  The input data does not form a proper triangle.");
            return null;
        }

        double[] c = new double[3];

        c[0] = a[0 + 2 * N];
        c[1] = a[1 + 2 * N];
        c[2] = 1.0 - c[0] - c[1];

        return c;
    }

    public static bool triangle_contains_point_1(double[] t, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_CONTAINS_POINT_1 finds if a point is inside a triangle.
        //
        //  Discussion:
        //
        //    It is conventional to list the triangle vertices in counter clockwise
        //    order.  However, this routine does not require a particular order
        //    for the vertices.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T[2*3], the triangle vertices.
        //
        //    Input, double P[2], the point to be checked.
        //
        //    Output, bool TRIANGLE_CONTAINS_POINT_1, is TRUE if the points 
        //    is inside the triangle or on its boundary, and FALSE otherwise.
        //
    {
        int i;

        double[] c = triangle_barycentric(t, p);

        bool value = true;

        for (i = 0; i < 3; i++)
        {
            value = c[i] switch
            {
                < 0.0 => false,
                _ => value
            };
        }

        return value;
    }
}