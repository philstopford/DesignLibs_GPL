using System;
using Burkardt.Geometry;
using Burkardt.Types;

namespace Burkardt.Polygon
{
    public static class Geometry
    {
        public static double polygon_1_2d(int n, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_1_2D integrates the function 1 over a polygon in 2D.
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
            //    Output, double POLYGON_1_2D, the value of the integral.
            //
        {
            int i;
            int im1;
            double result;

            result = 0.0;

            if (n < 3)
            {
                Console.WriteLine("");
                Console.WriteLine("POLYGON_1_2D - Warning!");
                Console.WriteLine("  The number of vertices must be at least 3.");
                Console.WriteLine("  The input value of N = " + n + "");
                return result;
            }

            for (i = 0; i < n; i++)
            {
                if (i == 0)
                {
                    im1 = n - 1;
                }
                else
                {
                    im1 = i - 1;
                }

                result = result + 0.5 * (v[0 + i * 2] + v[0 + im1 * 2])
                                      * (v[1 + i * 2] - v[1 + im1 * 2]);
            }

            return result;
        }

        public static double[] polygon_angles_2d(int n, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_ANGLES_2D computes the interior angles of a polygon in 2D.
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
            //    Output, double POLYGON_ANGLES_2D[N], the angles of the polygon,
            //    in radians.
            //
        {
            double[] angle;
            int i;
            int im1;
            int ip1;

            if (n < 1)
            {
                return null;
            }

            angle = new double[n];

            if (n <= 2)
            {
                for (i = 0; i < n; i++)
                {
                    angle[i] = 0.0;
                }

                return angle;
            }

            for (i = 0; i < n; i++)
            {
                im1 = typeMethods.i4_wrap(i - 1, 0, n - 1);
                ip1 = typeMethods.i4_wrap(i + 1, 0, n - 1);

                angle[i] = Angle.angle_rad_2d(v, v, v, p1Index: +im1 * 2, p2Index: +i * 2, p3Index: +ip1 * 2);
            }

            return angle;
        }

        public static double polygon_area_2d(int n, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_AREA_2D computes the area of a polygon in 2D.
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
            //    Output, double POLYGON_AREA_2D, the area of the polygon.
            //
        {
            double area;
            int i;
            int im1;
            int ip1;

            area = 0.0;

            for (i = 0; i < n; i++)
            {
                im1 = i - 1;
                ip1 = i + 1;

                if (im1 < 0)
                {
                    im1 = n - 1;
                }

                if (n <= ip1)
                {
                    ip1 = 0;
                }

                area = area + v[0 + i * 2] * (v[1 + ip1 * 2] - v[1 + im1 * 2]);
            }

            area = 0.5 * area;

            return area;
        }

        public static double polygon_area_2d_2(int n, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_AREA_2D_2 computes the area of a polygon in 2D.
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
            //    Output, double POLYGON_AREA_2D_2, the area of the polygon.
            //
        {
            double area;
            double area_triangle;
            int i;
            double[] t = new double[2 * 3];

            area = 0.0;

            for (i = 0; i < n - 2; i++)
            {
                t[0 + 0 * 2] = v[0 + i * 2];
                t[1 + 0 * 2] = v[1 + i * 2];
                t[0 + 1 * 2] = v[0 + (i + 1) * 2];
                t[1 + 1 * 2] = v[1 + (i + 1) * 2];
                t[0 + 2 * 2] = v[0 + (n - 1) * 2];
                t[1 + 2 * 2] = v[1 + (n - 1) * 2];

                area_triangle = typeMethods.triangle_area_2d(t);

                area = area + area_triangle;
            }

            return area;
        }

        public static double polygon_area_3d(int n, double[] v, ref double[] normal)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_AREA_3D computes the area of a polygon in 3D.
            //
            //  Discussion:
            //
            //    The computation is not valid unless the vertices really do lie
            //    in a plane, so that the polygon that is defined is "flat".
            //
            //    The polygon does not have to be "regular", that is, neither its
            //    sides nor its angles need to be equal.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 October 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Allen Van Gelder,
            //    Efficient Computation of Polygon Area and Polyhedron Volume,
            //    Graphics Gems V, edited by Alan Paeth,
            //    AP Professional, 1995, T385.G6975.
            //
            //  Parameters:
            //
            //    Input, int N, the number of vertices.
            //
            //    Input, double V[3*N], the coordinates of the vertices.
            //    The vertices should be listed in neighboring order.
            //
            //    Output, double NORMAL[3], the unit normal vector to the polygon.
            //
            //    Output, double POLYGON_AREA_3D, the area of the polygon.
            //
        {
            int DIM_NUM = 3;

            double area;
            int i;
            int ip1;

            normal[0] = 0.0;
            normal[1] = 0.0;
            normal[2] = 0.0;

            for (i = 0; i < n; i++)
            {
                if (i < n - 1)
                {
                    ip1 = i + 1;
                }
                else
                {
                    ip1 = 0;
                }

                //
                //  Compute the cross product and add it to NORMAL.
                //
                normal[0] = normal[0] + v[1 + i * 3] * v[2 + ip1 * 3] - v[2 + i * 3] * v[1 + ip1 * 3];
                normal[1] = normal[1] + v[2 + i * 3] * v[0 + ip1 * 3] - v[0 + i * 3] * v[2 + ip1 * 3];
                normal[2] = normal[2] + v[0 + i * 3] * v[1 + ip1 * 3] - v[1 + i * 3] * v[0 + ip1 * 3];
            }

            area = typeMethods.r8vec_norm(DIM_NUM, normal);

            if (area != 0.0)
            {
                normal[0] = normal[0] / area;
                normal[1] = normal[1] / area;
                normal[2] = normal[2] / area;
            }
            else
            {
                normal[0] = 1.0;
                normal[1] = 0.0;
                normal[2] = 0.0;
            }

            area = 0.5 * area;

            return area;
        }

        public static double polygon_area_3d_2(int n, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_AREA_3D_2 computes the area of a polygon in 3D.
            //
            //  Discussion:
            //
            //    The computation is not valid unless the vertices of the polygon
            //    lie in a plane, so that the polygon that is defined is "flat".
            //
            //    The polygon does not have to be "regular", that is, neither its
            //    sides nor its angles need to be equal.
            //
            //    The area is computed as the sum of the areas of the triangles
            //    formed by the last node with consecutive pairs of nodes (1,2),
            //    (2,3), ..., and (N-2,N-1).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 October 2005
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
            //    Input, double V[3*N], the coordinates of the vertices.
            //
            //    Output, double POLYGON_AREA_3D_2, the area of the polygon.
            //
        {
            int DIM_NUM = 3;

            double area;
            double[] area_vector = new double[DIM_NUM];
            double[] area_vector_triangle;
            int i;
            int j;
            double[] t = new double[DIM_NUM * 3];

            for (i = 0; i < DIM_NUM; i++)
            {
                area_vector[i] = 0.0;
            }

            for (j = 0; j < n - 2; j++)
            {
                t[0 + 0 * 3] = v[0 + j * 3];
                t[1 + 0 * 3] = v[1 + j * 3];
                t[2 + 0 * 3] = v[2 + j * 3];

                t[0 + 1 * 3] = v[0 + (j + 1) * 3];
                t[1 + 1 * 3] = v[1 + (j + 1) * 3];
                t[2 + 1 * 3] = v[2 + (j + 1) * 3];

                t[0 + 2 * 3] = v[0 + (n - 1) * 3];
                t[1 + 2 * 3] = v[1 + (n - 1) * 3];
                t[2 + 2 * 3] = v[2 + (n - 1) * 3];

                area_vector_triangle = Burkardt.TriangleNS.Geometry.triangle_area_vector_3d(t);

                for (i = 0; i < DIM_NUM; i++)
                {
                    area_vector[i] = area_vector[i] + area_vector_triangle[i];
                }

            }

            area = 0.5 * typeMethods.r8vec_norm(DIM_NUM, area_vector);

            return area;
        }

        public static double[] polygon_centroid_2d(int n, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_CENTROID_2D computes the centroid of a polygon in 2D.
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
            //    03 July 2005
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
            //    Output, double *POLYGON_CENTROID_2D[2], the coordinates of the centroid.
            //
        {
            int DIM_NUM = 2;

            double area;
            double[] centroid;
            int i;
            int ip1;
            int j;
            double temp;

            area = 0.0;
            centroid = new double[DIM_NUM];

            for (j = 0; j < DIM_NUM; j++)
            {
                centroid[j] = 0.0;
            }

            for (i = 0; i < n; i++)
            {
                if (i < n - 1)
                {
                    ip1 = i + 1;
                }
                else
                {
                    ip1 = 0;
                }

                temp = (v[0 + i * 2] * v[1 + ip1 * 2] - v[0 + ip1 * 2] * v[1 + i * 2]);

                area = area + temp;

                centroid[0] = centroid[0] + (v[0 + ip1 * 2] + v[0 + i * 2]) * temp;
                centroid[1] = centroid[1] + (v[1 + ip1 * 2] + v[1 + i * 2]) * temp;
            }

            area = area / 2.0;

            for (j = 0; j < DIM_NUM; j++)
            {
                centroid[j] = centroid[j] / (6.0 * area);
            }

            return centroid;
        }

        public static double[] polygon_centroid_2d_2(int n, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_CENTROID_2D_2 computes the centroid of a polygon in 2D.
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
            //    03 July 2005
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
            //    Output, double POLYGON_CENTROID_2D_2[2], the coordinates of the centroid.
            //
        {
            int DIM_NUM = 2;

            double area;
            double area_triangle;
            double[] centroid;
            int i;
            int j;
            double[] t = new double[DIM_NUM * 3];

            area = 0.0;
            centroid = new double[DIM_NUM];

            for (j = 0; j < DIM_NUM; j++)
            {
                centroid[j] = 0.0;
            }

            for (i = 0; i < n - 2; i++)
            {
                t[0 + 0 * 2] = v[0 + i * 2];
                t[1 + 0 * 2] = v[1 + i * 2];
                t[0 + 1 * 2] = v[0 + (i + 1) * 2];
                t[1 + 1 * 2] = v[1 + (i + 1) * 2];
                t[0 + 2 * 2] = v[0 + (n - 1) * 2];
                t[1 + 2 * 2] = v[1 + (n - 1) * 2];

                area_triangle = typeMethods.triangle_area_2d(t);

                area = area + area_triangle;
                centroid[0] = centroid[0]
                              + area_triangle * (v[0 + i * 2] + v[0 + (i + 1) * 2] + v[0 + (n - 1) * 2]) / 3.0;
                centroid[1] = centroid[1]
                              + area_triangle * (v[1 + i * 2] + v[1 + (i + 1) * 2] + v[1 + (n - 1) * 2]) / 3.0;
            }

            for (j = 0; j < DIM_NUM; j++)
            {
                centroid[j] = centroid[j] / area;
            }

            return centroid;

        }

        public static double[] polygon_centroid_3d(int n, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_CENTROID_3D computes the centroid of a polygon in 3D.
            //
            //  Discussion:
            //
            //    The centroid is the area-weighted sum of the centroids of
            //    disjoint triangles that make up the polygon.
            //
            //    Thanks to Jeremy Jarrett for pointing out a typographical error
            //    in an earlier version of this code.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 October 2005
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
            //    Input, double V[3*N], the coordinates of the vertices.
            //
            //    Output, double POLYGON_CENTROID_3D[3], the coordinates of the centroid.
            //
        {
            int DIM_NUM = 3;

            double area;
            double area_triangle;
            double[] centroid;
            int i;
            int j;
            double[] t = new double[DIM_NUM * 3];

            area = 0.0;
            centroid = new double[DIM_NUM];

            for (j = 0; j < DIM_NUM; j++)
            {
                centroid[j] = 0.0;
            }

            for (i = 0; i < n - 2; i++)
            {
                t[0 + 0 * 3] = v[0 + i * 3];
                t[1 + 0 * 3] = v[1 + i * 3];
                t[2 + 0 * 3] = v[2 + i * 3];

                t[0 + 1 * 3] = v[0 + (i + 1) * 3];
                t[1 + 1 * 3] = v[1 + (i + 1) * 3];
                t[2 + 1 * 3] = v[2 + (i + 1) * 3];

                t[0 + 2 * 3] = v[0 + (n - 1) * 3];
                t[1 + 2 * 3] = v[1 + (n - 1) * 3];
                t[2 + 2 * 3] = v[2 + (n - 1) * 3];

                area_triangle = typeMethods.triangle_area_3d(t);

                area = area + area_triangle;

                centroid[0] = centroid[0]
                              + area_triangle * (v[0 + i * 3] + v[0 + (i + 1) * 3] + v[0 + (n - 1) * 3]) / 3.0;
                centroid[1] = centroid[1]
                              + area_triangle * (v[1 + i * 3] + v[1 + (i + 1) * 3] + v[1 + (n - 1) * 3]) / 3.0;
                centroid[2] = centroid[2]
                              + area_triangle * (v[2 + i * 3] + v[2 + (i + 1) * 3] + v[2 + (n - 1) * 3]) / 3.0;
            }

            for (j = 0; j < DIM_NUM; j++)
            {
                centroid[j] = centroid[j] / area;
            }

            return centroid;

        }

        public static bool polygon_contains_point_2d(int n, double[] v, double[] p)

            //****************************************************************************9-
            //
            //  Purpose:
            //
            //    POLYGON_CONTAINS_POINT_2D finds if a point is inside a polygon in 2D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 November 2016
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
            //    Output, bool POLYGON_CONTAINS_POINT_2D, is TRUE if the point
            //    is inside the polygon or on its boundary, and FALSE otherwise.
            //
        {
            int i;
            int ip1;
            double px1;
            double px2;
            double py1;
            double py2;
            bool value;
            double xints;

            value = false;

            px1 = v[0 + 2 * 0];
            py1 = v[1 + 2 * 0];
            xints = p[0] - 1.0;

            for (i = 0; i < n; i++)
            {
                ip1 = ((i + 1) % n);
                px2 = v[0 + 2 * ip1];
                py2 = v[1 + 2 * ip1];

                if (Math.Min(py1, py2) < p[1])
                {
                    if (p[1] <= Math.Max(py1, py2))
                    {
                        if (p[0] <= Math.Max(px1, px2))
                        {
                            if (py1 != py2)
                            {
                                xints = (p[1] - py1) * (px2 - px1) / (py2 - py1) + px1;
                            }

                            if (px1 == px2 || p[0] <= xints)
                            {
                                value = !value;
                            }
                        }
                    }
                }

                px1 = px2;
                py1 = py2;
            }

            return value;
        }

        public static bool polygon_contains_point_2d_2(int n, double[] v, double[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_CONTAINS_POINT_2D_2 finds if a point is inside a convex polygon in 2D.
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
            //  Parameters:
            //
            //    Input, int N, the number of nodes or vertices in the polygon.
            //    N must be at least 3.
            //
            //    Input, double V[2*N], the coordinates of the vertices.
            //
            //    Input, double P[2], the coordinates of the point to be tested.
            //
            //    Output, bool POLYGON_CONTAINS_POINT_2D_2, is TRUE if the point
            //    is inside the polygon or on its boundary, and FALSE otherwise.
            //
        {
            int DIM_NUM = 2;

            int i;
            double[] t = new double[DIM_NUM * 3];
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

                if (Burkardt.TriangleNS.Geometry.triangle_contains_point_2d_1(t, p))
                {
                    return true;
                }
            }

            return false;

        }

        public static bool polygon_contains_point_2d_3(int n, double[] v, double[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_CONTAINS_POINT_2D_3 finds if a point is inside a simple polygon in 2D.
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
            //    06 August 2005
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
            //    Output, bool POLYGON_CONTAINS_POINT_2D_3, is TRUE if the point
            //    is inside the polygon or on its boundary, and FALSE otherwise.
            //
        {
            int i;
            bool value;
            double x1;
            double x2;
            double y1;
            double y2;

            value = false;

            for (i = 0; i < n; i++)
            {
                x1 = v[0 + i * 2];
                y1 = v[1 + i * 2];

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

                if ((y1 < p[1] && p[1] <= y2) ||
                    (p[1] <= y1 && y2 < p[1]))
                {
                    if ((p[0] - x1) - (p[1] - y1) * (x2 - x1) / (y2 - y1) < 0)
                    {
                        value = !value;
                    }
                }
            }

            return value;
        }

        public static double polygon_diameter_2d(int n, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_DIAMETER_2D computes the diameter of a polygon in 2D.
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
            //    03 July 2005
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
            //    Output, double POLYGON_DIAMETER_2D, the diameter of the polygon.
            //
        {
            double diameter;
            int i;
            int j;

            diameter = 0.0;

            for (i = 0; i < n; i++)
            {
                for (j = i + 1; j < n; j++)
                {
                    diameter = Math.Max(diameter,
                        Math.Sqrt((v[0 + i * 2] - v[0 + j * 2]) * (v[0 + i * 2] - v[0 + j * 2])
                                  + (v[1 + i * 2] - v[1 + j * 2]) * (v[1 + i * 2] - v[1 + j * 2])));
                }
            }

            return diameter;
        }

        public static double[] polygon_expand_2d(int n, double[] v, double h)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_EXPAND_2D expands a polygon in 2D.
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
            //    19 July 2005
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
            //    Output, double POLYGON_EXPAND_2D[2*N], the "expanded" coordinates.
            //
        {
            int DIM_NUM = 2;

            double angle;
            double h2;
            int i;
            int j;
            int jm1;
            int jp1;
            double[] p4;
            double[] w;

            w = new double[DIM_NUM * n];
            //
            //  Consider each angle, formed by the nodes P(I-1), P(I), P(I+1).
            //
            for (j = 0; j < n; j++)
            {
                jm1 = typeMethods.i4_wrap(j - 1, 0, n - 1);
                jp1 = typeMethods.i4_wrap(j + 1, 0, n - 1);
                //
                //        P1
                //        /
                //       /   P4
                //      /  .
                //     / .
                //    P2--------->P3
                //
                p4 = Angle.angle_half_2d(v, v, v, p1Index: +jm1 * DIM_NUM, p2Index: +j * DIM_NUM,
                    p3Index: +jp1 * DIM_NUM);
                //
                //  Compute the value of the half angle.
                //
                angle = Angle.angle_rad_2d(v, v, p4, p1Index: +jm1 * DIM_NUM, p2Index: +j * DIM_NUM);
                //
                //  The stepsize along the ray must be adjusted so that the sides
                //  move out by H.
                //
                h2 = h / Math.Sin(angle);

                for (i = 0; i < DIM_NUM; i++)
                {
                    w[i + j * DIM_NUM] = v[i + j * DIM_NUM] - h2 * (p4[i] - v[i + j * DIM_NUM]);
                }

            }

            return w;

        }

        public static void polygon_inrad_data_2d(int n, double radin, ref double area, ref double radout,
                ref double side)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_INRAD_DATA_2D determines polygonal data from its inner radius in 2D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    09 September 2003
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
            //    Output, double *AREA, the area of the regular polygon.
            //
            //    Output, double *RADOUT, the outer radius of the polygon, that is,
            //    the radius of the smallest circle that can be described about
            //    the polygon.
            //
            //    Output, double *SIDE, the length of one side of the polygon.
            //
        {
            double angle;

            if (n < 3)
            {
                Console.WriteLine("");
                Console.WriteLine("POLYGON_INRAD_DATA_2D - Fatal error!");
                Console.WriteLine("  Input value of N must be at least 3,");
                Console.WriteLine("  but your input value was N = " + n + "");
                return;
            }

            angle = Math.PI / ((double) n);
            area = ((double) n) * radin * radin * Math.Tan(angle);
            side = 2.0 * radin * Math.Tan(angle);
            radout = 0.5 * (side) / Math.Sin(angle);

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
            int NOT_CONVEX = -1;
            int DEGENERATE_CONVEX = 0;
            int CONVEX_CCW = 1;
            int CONVEX_CW = 2;

            double angle;
            double cross;
            double dot;
            double exterior_total;
            int i;
            int ip1;
            int ip2;
            double sense;
            double TOL = 1.0;
            int value = 0;

            exterior_total = 0.0;
            //
            //  If there are not at least 3 distinct vertices, we are done.
            //
            if (n < 3)
            {
                return DEGENERATE_CONVEX;
            }

            sense = 0.0;
            //
            //  Consider each polygonal vertex I.
            //
            for (i = 0; i < n; i++)
            {
                ip1 = typeMethods.i4_wrap(i + 1, 0, n - 1);
                ip2 = typeMethods.i4_wrap(i + 2, 0, n - 1);

                dot = (v[0 + ip2 * 2] - v[0 + ip1 * 2]) * (v[0 + i * 2] - v[0 + ip1 * 2])
                      + (v[1 + ip2 * 2] - v[1 + ip1 * 2]) * (v[1 + i * 2] - v[1 + ip1 * 2]);

                cross = (v[0 + ip2 * 2] - v[0 + ip1 * 2]) * (v[1 + i * 2] - v[1 + ip1 * 2])
                        - (v[0 + i * 2] - v[0 + ip1 * 2]) * (v[1 + ip2 * 2] - v[1 + ip1 * 2]);

                angle = Math.Atan2(cross, dot);
                //
                //  See if the turn defined by this vertex is our first indication of
                //  the "sense" of the polygon, or if it disagrees with the previously
                //  defined sense.
                //
                if (sense == 0.0)
                {
                    if (angle < 0.0)
                    {
                        sense = -1.0;
                    }
                    else if (0.0 < angle)
                    {
                        sense = +1.0;
                    }
                }
                else if (sense == 1.0)
                {
                    if (angle < 0.0)
                    {
                        return NOT_CONVEX;
                    }
                }
                else if (sense == -1.0)
                {
                    if (0.0 < angle)
                    {
                        return NOT_CONVEX;
                    }
                }

                //
                //  If the exterior total is greater than 360, then the polygon is
                //  going around again.
                //
                angle = Math.Atan2(-cross, -dot);

                exterior_total = exterior_total + angle;

                if (360.0 + TOL < Helpers.radians_to_degrees(Math.Abs(exterior_total)))
                {
                    return NOT_CONVEX;
                }

            }

            if (sense == +1.0)
            {
                value = CONVEX_CCW;
            }
            else if (sense == -1.0)
            {
                value = CONVEX_CW;
            }

            return value;
        }

        public static double polygon_lattice_area_2d(int i, int b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_LATTICE_AREA_2D computes the area of a lattice polygon in 2D.
            //
            //  Discussion:
            //
            //    We define a lattice to be the 2D plane, in which the points
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
            //    09 September 2003
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
            //    Output, double POLYGON_LATTICE_AREA_2D, the area of the lattice polygon.
            //
        {
            double value;

            value = (double) i + ((double) b) / 2.0 - 1.0;

            return value;
        }

        public static double[] polygon_normal_3d(int n, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_NORMAL_3D computes the normal vector to a polygon in 3D.
            //
            //  Discussion:
            //
            //    If the polygon is planar, then this calculation is correct.
            //
            //    Otherwise, the normal vector calculated is the simple average
            //    of the normals defined by the planes of successive triples
            //    of vertices.
            //
            //    If the polygon is "almost" planar, this is still acceptable.
            //    But as the polygon is less and less planar, so this averaged normal
            //    vector becomes more and more meaningless.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 August 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Paulo Cezar Pinto Carvalho, Paulo Roma Cavalcanti,
            //    Point in Polyhedron Testing Using Spherical Polygons,
            //    in Graphics Gems V,
            //    edited by Alan Paeth,
            //    Academic Press, 1995, T385.G6975.
            //
            //  Parameters:
            //
            //    Input, int N, the number of vertices.
            //
            //    Input, double V[3*N], the coordinates of the vertices.
            //
            //    Output, double POLYGON_NORMAL_3D[3], the averaged normal vector
            //    to the polygon.
            //
        {
            int DIM_NUM = 3;

            int i;
            int j;
            double[] normal;
            double normal_norm;
            double[] p;
            double[] v1;
            double[] v2;

            normal = new double[DIM_NUM];
            v1 = new double[DIM_NUM];
            v2 = new double[DIM_NUM];

            typeMethods.r8vec_zero(DIM_NUM, ref normal);

            for (i = 0; i < DIM_NUM; i++)
            {
                v1[i] = v[i + 1 * DIM_NUM] - v[i + 0 * DIM_NUM];
            }

            for (j = 2; j < n; j++)
            {
                for (i = 0; i < DIM_NUM; i++)
                {
                    v2[i] = v[i + j * DIM_NUM] - v[i + 0 * DIM_NUM];
                }

                p = typeMethods.r8vec_cross_product_3d(v1, v2);
                for (i = 0; i < DIM_NUM; i++)
                {
                    normal[i] = normal[i] + p[i];
                }

                typeMethods.r8vec_copy(DIM_NUM, v2, ref v1);

            }

            //
            //  Normalize.
            //
            normal_norm = typeMethods.r8vec_norm(DIM_NUM, normal);

            if (normal_norm != 0.0)
            {
                for (i = 0; i < DIM_NUM; i++)
                {
                    normal[i] = normal[i] / normal_norm;
                }
            }

            return normal;

        }

        public static void polygon_outrad_data_2d(int n, double radout, ref double area, ref double radin,
                ref double side)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_OUTRAD_DATA_2D determines polygonal data from its outer radius in 2D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    09 September 2003
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
            //    Output, double *AREA, the area of the regular polygon.
            //
            //    Output, double *RADIN, the inner radius of the polygon, that is,
            //    the radius of the largest circle that can be inscribed
            //    within the polygon.
            //
            //    Output, double *SIDE, the length of one side of the polygon.
            //
        {
            double angle;

            if (n < 3)
            {
                Console.WriteLine("");
                Console.WriteLine("POLYGON_OUTRAD_DATA_2D - Fatal error!");
                Console.WriteLine("  Input value of N must be at least 3,");
                Console.WriteLine("  but your input value was N = " + n + "");
                return;
            }

            angle = Math.PI / ((double) n);
            area = 0.5 * ((double) n) * radout * radout * Math.Sin(2.0 * angle);
            side = 2.0 * radout * Math.Sin(angle);
            radin = 0.5 * (side) / Math.Tan(angle);

        }

        public static void polygon_side_data_2d(int n, double side, ref double area, ref double radin,
                ref double radout)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_SIDE_DATA_2D determines polygonal data from its side length in 2D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 June 2005
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
            //    Output, double *AREA, the area of the regular polygon.
            //
            //    Output, double *RADIN, the inner radius of the polygon, that is,
            //    the radius of the largest circle that can be inscribed within
            //    the polygon.
            //
            //    Output, double *RADOUT, the outer radius of the polygon, that is,
            //    the radius of the smallest circle that can be described about
            //    the polygon.
            //
        {
            double angle;

            if (n < 3)
            {
                Console.WriteLine("");
                Console.WriteLine("POLYGON_SIDE_DATA_2D - Fatal error!");
                Console.WriteLine("  Input value of N must be at least 3,");
                Console.WriteLine("  but your input value was N = " + n + "");
                return;
            }

            angle = Math.PI / ((double) n);
            area = 0.25 * n * side * side / Math.Tan(angle);
            radin = 0.5 * side / Math.Tan(angle);
            radout = 0.5 * side / Math.Sin(angle);

        }

        public static double polygon_solid_angle_3d(int n, double[] v, double[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_SOLID_ANGLE_3D calculates the projected solid angle of a 3D plane polygon.
            //
            //  Discussion:
            //
            //    A point P is at the center of a unit sphere.  A planar polygon
            //    is to be projected onto the surface of this sphere, by drawing
            //    the ray from P to each polygonal vertex, and noting where this ray
            //    intersects the sphere.
            //
            //    We compute the area on the sphere of the projected polygon.
            //
            //    Since we are projecting the polygon onto a unit sphere, the area
            //    of the projected polygon is equal to the solid angle subtended by
            //    the polygon.
            //
            //    The value returned by this routine will include a sign.  The
            //    angle subtended will be NEGATIVE if the normal vector defined by
            //    the polygon points AWAY from the viewing point, and will be
            //    POSITIVE if the normal vector points towards the viewing point.
            //
            //    If the orientation of the polygon is of no interest to you,
            //    then you can probably simply take the absolute value of the
            //    solid angle as the information you want.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 October 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Paulo Cezar Pinto Carvalho, Paulo Roma Cavalcanti,
            //    Point in Polyhedron Testing Using Spherical Polygons,
            //    in Graphics Gems V,
            //    edited by Alan Paeth,
            //    Academic Press, 1995,
            //    ISBN: 0125434553,
            //    LC: T385.G6975.
            //
            //  Parameters:
            //
            //    Input, int N, the number of vertices.
            //
            //    Input, double V[3*N], the coordinates of the vertices.
            //
            //    Input, double P[3], the point at the center of the unit sphere.
            //
            //    Output, double POLYGON_SOLID_ANGLE_3D, the solid angle subtended
            //    by the polygon, as projected onto the unit sphere around the point P.
            //
        {
            int DIM_NUM = 3;

            double[] a = new double[DIM_NUM];
            double angle;
            double area = 0.0;
            double[] b = new double[DIM_NUM];
            int j;
            int jp1;
            double[] normal1;
            double normal1_norm;
            double[] normal2;
            double normal2_norm;
            double[] plane;
            double[] r1 = new double[DIM_NUM];
            double s;
            double value;

            if (n < 3)
            {
                return 0.0;
            }

            plane = polygon_normal_3d(n, v);

            a[0] = v[0 + (n - 1) * DIM_NUM] - v[0 + 0 * DIM_NUM];
            a[1] = v[1 + (n - 1) * DIM_NUM] - v[1 + 0 * DIM_NUM];
            a[2] = v[2 + (n - 1) * DIM_NUM] - v[2 + 0 * DIM_NUM];

            for (j = 0; j < n; j++)
            {
                r1[0] = v[0 + j * DIM_NUM] - p[0];
                r1[1] = v[1 + j * DIM_NUM] - p[1];
                r1[2] = v[2 + j * DIM_NUM] - p[2];

                jp1 = typeMethods.i4_wrap(j + 1, 0, n - 1);

                b[0] = v[0 + jp1 * DIM_NUM] - v[0 + j * DIM_NUM];
                b[1] = v[1 + jp1 * DIM_NUM] - v[1 + j * DIM_NUM];
                b[2] = v[2 + jp1 * DIM_NUM] - v[2 + j * DIM_NUM];

                normal1 = typeMethods.r8vec_cross_product_3d(a, r1);

                normal1_norm = typeMethods.r8vec_norm(DIM_NUM, normal1);

                normal2 = typeMethods.r8vec_cross_product_3d(r1, b);

                normal2_norm = typeMethods.r8vec_norm(DIM_NUM, normal2);

                s = typeMethods.r8vec_dot_product(DIM_NUM, normal1, normal2)
                    / (normal1_norm * normal2_norm);

                angle = typeMethods.r8_acos(s);

                s = typeMethods.r8vec_scalar_triple_product(b, a, plane);

                if (0.0 < s)
                {
                    area = area + Math.PI - angle;
                }
                else
                {
                    area = area + Math.PI + angle;
                }

                a[0] = -b[0];
                a[1] = -b[1];
                a[2] = -b[2];

            }

            area = area - Math.PI * (double) (n - 2);

            if (0.0 < typeMethods.r8vec_dot_product(DIM_NUM, plane, r1))
            {
                value = -area;
            }
            else
            {
                value = area;
            }

            return value;

        }

        public static double polygon_x_2d(int n, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_X_2D integrates the function X over a polygon in 2D.
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
            //    Input, double V[2*N], the coordinates of the vertices.
            //
            //    Output, double POLYGON_X_2D, the value of the integral.
            //
        {
            int i;
            int im1;
            double result;

            result = 0.0;

            if (n < 3)
            {
                Console.WriteLine("");
                Console.WriteLine("POLYGON_X_2D - Warning!");
                Console.WriteLine("  The number of vertices must be at least 3.");
                Console.WriteLine("  The input value of N = " + n + "");
                return result;
            }

            for (i = 0; i < n; i++)
            {
                if (i == 0)
                {
                    im1 = n - 1;
                }
                else
                {
                    im1 = i - 1;
                }

                result = result + (v[0 + i * 2] * v[0 + i * 2]
                                   + v[0 + i * 2] * v[0 + im1 * 2]
                                   + v[0 + im1 * 2] * v[0 + im1 * 2])
                    * (v[1 + i * 2] - v[1 + im1 * 2]);
            }

            result = result / 6.0;

            return result;
        }

        public static double polygon_y_2d(int n, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_Y_2D integrates the function Y over a polygon in 2D.
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
            //    Input, double V[2*N], the coordinates of the vertices.
            //
            //    Output, double POLYGON_Y_2D, the value of the integral.
            //
        {
            int i;
            int im1;
            double result;

            result = 0.0;

            if (n < 3)
            {
                Console.WriteLine("");
                Console.WriteLine("POLYGON_Y_2D - Warning!");
                Console.WriteLine("  The number of vertices must be at least 3.");
                Console.WriteLine("  The input value of N = " + n + "");
                return result;
            }

            for (i = 0; i < n; i++)
            {
                if (i == 0)
                {
                    im1 = n - 1;
                }
                else
                {
                    im1 = i - 1;
                }

                result = result - (v[1 + i * 2] * v[1 + i * 2]
                                   + v[1 + i * 2] * v[1 + im1 * 2]
                                   + v[1 + im1 * 2] * v[1 + im1 * 2])
                    * (v[0 + i * 2] - v[0 + im1 * 2]);
            }

            result = result / 6.0;

            return result;
        }

        public static double polygon_xx_2d(int n, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_XX_2D integrates the function X*X over a polygon in 2D.
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
            //    Input, double V[2*N], the coordinates of the vertices.
            //
            //    Output, double POLYGON_XX_2D, the value of the integral.
            //
        {
            int i;
            int im1;
            double result;

            result = 0.0;

            if (n < 3)
            {
                Console.WriteLine("");
                Console.WriteLine("POLYGON_XX_2D - Warning!");
                Console.WriteLine("  The number of vertices must be at least 3.");
                Console.WriteLine("  The input value of N = " + n + "");
                return result;
            }

            for (i = 0; i < n; i++)
            {
                if (i == 0)
                {
                    im1 = n - 1;
                }
                else
                {
                    im1 = i - 1;
                }

                result = result + (
                    v[0 + i * 2] * v[0 + i * 2] * v[0 + i * 2]
                    + v[0 + i * 2] * v[0 + i * 2] * v[0 + im1 * 2]
                    + v[0 + i * 2] * v[0 + im1 * 2] * v[0 + im1 * 2]
                    + v[0 + im1 * 2] * v[0 + im1 * 2] * v[0 + im1 * 2]) * (v[1 + i * 2] - v[1 + im1 * 2]);
            }

            result = result / 12.0;

            return result;
        }

        public static double polygon_xy_2d(int n, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_XY_2D integrates the function X*Y over a polygon in 2D.
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
            //    Input, double V[2*N], the coordinates of the vertices.
            //
            //    Output, double POLYGON_XY_2D, the value of the integral.
            //
        {
            int i;
            int im1;
            double result;

            result = 0.0;

            if (n < 3)
            {
                Console.WriteLine("");
                Console.WriteLine("POLYGON_XY_2D - Warning!");
                Console.WriteLine("  The number of vertices must be at least 3.");
                Console.WriteLine("  The input value of N = " + n + "");
                return result;
            }

            for (i = 0; i < n; i++)
            {
                if (i == 0)
                {
                    im1 = n - 1;
                }
                else
                {
                    im1 = i - 1;
                }

                result = result + (
                    v[1 + i * 2] * (3.0 * v[0 + i * 2] * v[0 + i * 2]
                                    + 2.0 * v[0 + i * 2] * v[0 + im1 * 2]
                                    + v[0 + im1 * 2] * v[0 + im1 * 2])
                    + v[1 + im1 * 2] * (v[0 + i * 2] * v[0 + i * 2]
                                        + 2.0 * v[0 + i * 2] * v[0 + im1 * 2]
                                        + 3.0 * v[0 + im1 * 2] * v[0 + im1 * 2])
                ) * (v[1 + i * 2] - v[1 + im1 * 2]);
            }

            result = result / 24.0;

            return result;
        }

        public static double polygon_yy_2d(int n, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_YY_2D integrates the function Y*Y over a polygon in 2D.
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
            //    Input, double V[2*N], the coordinates of the vertices.
            //
            //    Output, double POLYGON_YY_2D, the value of the integral.
            //
            //
        {
            int i;
            int im1;
            double result;

            result = 0.0;

            if (n < 3)
            {
                Console.WriteLine("");
                Console.WriteLine("POLYGON_YY_2D - Warning!");
                Console.WriteLine("  The number of vertices must be at least 3.");
                Console.WriteLine("  The input value of N = " + n + "");
                return result;
            }

            for (i = 0; i < n; i++)
            {

                if (i == 0)
                {
                    im1 = n - 1;
                }
                else
                {
                    im1 = i - 1;
                }

                result = result -
                         (v[1 + i * 2] * v[1 + i * 2] * v[1 + i * 2]
                          + v[1 + i * 2] * v[1 + i * 2] * v[1 + im1 * 2]
                          + v[1 + i * 2] * v[1 + im1 * 2] * v[1 + im1 * 2]
                          + v[1 + im1 * 2] * v[1 + im1 * 2] * v[1 + im1 * 2])
                         * (v[0 + i * 2] - v[0 + im1 * 2]);
            }

            result = result / 12.0;

            return result;
        }
    }
}