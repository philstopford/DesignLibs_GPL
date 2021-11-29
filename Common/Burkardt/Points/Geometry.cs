using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Geometry;
using Burkardt.Types;

namespace Burkardt.PointsNS;

public static class Geometry
{
    public static bool points_avoid_point_naive_2d(int n, double[] pset, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINTS_AVOID_POINT_NAIVE_2D finds if a point is "far enough" from a set of points in 2D.
        //
        //  Discussion:
        //
        //    The routine discards points that are too close to other points.
        //    The method used to check this is quadratic in the number of points,
        //    and may take an inordinate amount of time if there are a large
        //    number of points.  But in that case, what do you want?  If you want
        //    lots of points, you don't want to delete any because it won't matter.
        //
        //    The test point is "far enough" from an accepted point if
        //    the Euclidean distance is at least 100 times EPSILON.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of accepted points.
        //
        //    Input, double PSET[2*N], the accepted points.  The points are stored
        //    in a one dimensional array, beginning with the X and Y coordinates of
        //    the first point, and so on.
        //
        //    Input, double P[2], a point to be tested.
        //
        //    Output, bool POINTS_AVOID_POINT_NAIVE_2D, is TRUE if P is
        //    "far enough" from all the accepted points.
        //
    {
        int j;

        double tolsq = 100.0 * typeMethods.r8_epsilon();
        tolsq *= tolsq;

        for (j = 0; j < n; j++)
        {
            double normsq = (pset[0 + j * 2] - p[0]) * (pset[0 + j * 2] - p[0])
                            + (pset[1 + j * 2] - p[1]) * (pset[1 + j * 2] - p[1]);

            if (normsq < tolsq)
            {
                return false;
            }
        }

        return true;
    }

    public static void points_bisect_line_imp_2d(double[] p1, double[] p2, ref double a,
            ref double b, ref double c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINTS_BISECT_LINE_IMP_2D finds the implicit line bisecting the line between two points in 2D.
        //
        //  Discussion:
        //
        //    The implicit form of a line in 2D is:
        //
        //      A * X + B * Y + C = 0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 July 2005
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
        //    Input, double P1[2], P2[2], the coordinates of two points.
        //
        //    Output, double *A, *B, *C, the parameters of the implicit line
        //    equidistant from both points.
        //
    {
        a = p1[0] - p2[0];
        b = p1[1] - p2[1];
        c = -0.5 * (p1[0] * p1[0] + p1[1] * p1[1]
                    - (p2[0] * p2[0] + p2[1] * p2[1]));
    }

    public static void points_bisect_line_par_2d(double[] p1, double[] p2, ref double f,
            ref double g, ref double x, ref double y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINTS_BISECT_LINE_PAR_2D finds the parametric line bisecting the line between two points in 2D.
        //
        //  Discussion:
        //
        //    The parametric form of a line in 2D is:
        //
        //      X = X0 + F * T
        //      Y = Y0 + G * T
        //
        //    For normalization, we choose F*F+G*G = 1 and 0 <= F.
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
        //  Reference:
        //
        //    Adrian Bowyer, John Woodwark,
        //    A Programmer's Geometry,
        //    Butterworths, 1983.
        //
        //  Parameters:
        //
        //    Input, double P1[2], P2[2], the coordinates of two points.
        //
        //    Output, double *F, *G, *X, *Y, the parameters of the parametric line
        //    equidistant from both points.
        //
    {
        f = 0.5 * (p1[0] + p2[0]);
        g = 0.5 * (p1[1] + p2[1]);

        double norm = Math.Sqrt(Math.Pow(f, 2) + Math.Pow(g, 2));

        switch (norm)
        {
            case > 0.0:
                f /= norm;
                g /= norm;
                break;
        }

        switch (f)
        {
            case < 0.0:
                f = -f;
                g = -g;
                break;
        }

        x = -(p2[1] - p1[1]);
        y = p2[0] - p1[0];

    }

    public static int points_centroid_2d(int n, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINTS_CENTROID_2D computes the discrete centroid of a point set in 2D.
        //
        //  Discussion:
        //
        //    Given a discrete set of points S, the discrete centroid z is defined by
        //
        //                           Sum ( x in S ) ( x - z )^2
        //        = min ( y in S ) { Sum ( x in S ) ( x - y )^2
        //
        //    In other words, the discrete centroid is a point in the set whose distance
        //    to the other points is minimized.  The discrete centroid of a point set
        //    need not be unique.  Consider a point set that comprises the
        //    vertices of an equilateral triangle.
        //
        //    This discrete centroid may also be referred to as the K-means cluster.
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
        //    Input, int N, the number of points.
        //
        //    Input, double P[2*N], the coordinates of the points.
        //
        //    Output, int POINTS_CENTROID_2D, the index of a discrete
        //    centroid of the set, between 0 and N-1.
        //
    {
        int i;

        double dist_min = 0.0;
        int cent = -1;

        for (i = 0; i < n; i++)
        {
            double dist = 0.0;
            int j;
            for (j = 0; j < n; j++)
            {
                dist = dist + (p[0 + i * 2] - p[0 + j * 2]) * (p[0 + i * 2] - p[0 + j * 2])
                            + (p[1 + i * 2] - p[1 + j * 2]) * (p[1 + i * 2] - p[1 + j * 2]);
            }

            switch (i)
            {
                case 0:
                    dist_min = dist;
                    cent = i;
                    break;
                default:
                {
                    if (dist < dist_min)
                    {
                        dist_min = dist;
                        cent = i;
                    }

                    break;
                }
            }
        }

        return cent;
    }

    public static double points_colin_2d(double[] p1, double[] p2, double[] p3)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINTS_COLIN_2D estimates the colinearity of 3 points in 2D.
        //
        //  Discussion:
        //
        //    The estimate of collinearity is the ratio of the area of the triangle
        //    spanned by the points to the area of the equilateral triangle with the
        //    same perimeter.
        //
        //    This is 1.0 if the points are maximally noncolinear, 0.0 if the
        //    points are exactly colinear, and otherwise is closer to 1 or 0 depending
        //    on whether the points are far or close to colinearity.
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
        //    Input, double P1[2], P2[2], P3[2], the coordinates of the points.
        //
        //    Output, double POINTS_COLIN_2D, an estimate of colinearity,
        //
    {
        const int DIM_NUM = 2;

        double colin;
        double[] t = new double[DIM_NUM * 3];

        t[0 + 0 * 2] = p1[0];
        t[1 + 0 * 2] = p1[1];
        t[0 + 1 * 2] = p2[0];
        t[1 + 1 * 2] = p2[1];
        t[0 + 2 * 2] = p3[0];
        t[1 + 2 * 2] = p3[1];

        double area_triangle = typeMethods.triangle_area_2d(t);

        switch (area_triangle)
        {
            case 0.0:
                colin = 0.0;
                break;
            default:
                double s12 = Math.Sqrt(Math.Pow(p2[0] - p1[0], 2) + Math.Pow(p2[1] - p1[1], 2));
                double s23 = Math.Sqrt(Math.Pow(p3[0] - p2[0], 2) + Math.Pow(p3[1] - p2[1], 2));
                double s31 = Math.Sqrt(Math.Pow(p1[0] - p3[0], 2) + Math.Pow(p1[1] - p3[1], 2));

                double perim = s12 + s23 + s31;

                double side = perim / 3.0;

                double area2 = 0.25 * Math.Sqrt(3.0) * side * side;

                colin = Math.Abs(area_triangle) / area2;
                break;
        }

        return colin;
    }

    public static double points_colin_3d(double[] p1, double[] p2, double[] p3)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINTS_COLIN_3D estimates the colinearity of 3 points in 3D.
        //
        //  Discussion:
        //
        //    The estimate of collinearity is the ratio of the area of the triangle
        //    spanned by the points to the area of the equilateral triangle with the
        //    same perimeter.
        //
        //    This is 1.0 if the points are maximally noncolinear, 0.0 if the
        //    points are exactly colinear, and otherwise is closer to 1 or 0 depending
        //    on whether the points are far or close to colinearity.
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
        //    Input, double P1[3], P2[3], P3[3], the points.
        //
        //    Output, double POINTS_COLIN_3D, an estimate of colinearity.
        //
    {
        const int DIM_NUM = 3;

        double colin;
        double[] t = new double[DIM_NUM * 3];

        t[0 + 0 * 3] = p1[0];
        t[1 + 0 * 3] = p1[1];
        t[2 + 0 * 3] = p1[2];
        t[0 + 1 * 3] = p2[0];
        t[1 + 1 * 3] = p2[1];
        t[2 + 1 * 3] = p2[2];
        t[0 + 2 * 3] = p3[0];
        t[1 + 2 * 3] = p3[1];
        t[2 + 2 * 3] = p3[2];

        double area_triangle = typeMethods.triangle_area_3d(t);

        switch (area_triangle)
        {
            case 0.0:
                colin = 0.0;
                break;
            default:
                double s12 = Math.Sqrt(Math.Pow(p2[0] - p1[0], 2)
                                       + Math.Pow(p2[1] - p1[1], 2)
                                       + Math.Pow(p2[2] - p1[2], 2));
                double s23 = Math.Sqrt(Math.Pow(p3[0] - p2[0], 2)
                                       + Math.Pow(p3[1] - p2[1], 2)
                                       + Math.Pow(p3[2] - p2[2], 2));
                double s31 = Math.Sqrt(Math.Pow(p1[0] - p3[0], 2)
                                       + Math.Pow(p1[1] - p3[1], 2)
                                       + Math.Pow(p1[2] - p3[2], 2));

                double perim = s12 + s23 + s31;

                double side = perim / 3.0;

                double area2 = 0.25 * Math.Sqrt(3.0) * side * side;

                colin = Math.Abs(area_triangle) / area2;
                break;
        }

        return colin;
    }

    public static double points_dist_2d(double[] p1, double[] p2, int p1Index = 0, int p2Index = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINTS_DIST_2D finds the distance between two points in 2D.
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
        //    Input, double P1[2], P2[2], two points.
        //
        //    Output, double POINTS_DIST_2D, the distance between the points.
        //
    {
        double dist = Math.Sqrt(Math.Pow(p1[(0 + p1Index) % p1.Length] - p2[(0 + p2Index) % p2.Length], 2)
                                + Math.Pow(p1[(1 + p1Index) % p1.Length] - p2[(1 + p2Index) % p2.Length], 2));

        return dist;
    }

    public static double points_dist_3d(double[] p1, double[] p2, int p1Index = 0, int p2Index = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINTS_DIST_3D finds the distance between two points in 3D.
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
        //    Input, double P1[3], P2[3], two points.
        //
        //    Output, double POINTS_DIST_3D, the distance between the points.
        //
    {
        double dist = Math.Sqrt(Math.Pow(p1[(0 + p1Index) % p1.Length] - p2[(0 + p2Index) % p2.Length], 2)
                                + Math.Pow(p1[(1 + p1Index) % p1.Length] - p2[(1 + p2Index) % p2.Length], 2)
                                + Math.Pow(p1[(2 + p1Index) % p1.Length] - p2[(2 + p2Index) % p2.Length], 2));

        return dist;
    }

    public static double points_dist_nd(int dim_num, double[] p1, double[] p2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINTS_DIST_ND finds the distance between two points in ND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the dimension of the space.
        //
        //    Input, double P1[DIM_NUM], P2[DIM_NUM], the coordinates of two points.
        //
        //    Output, double POINTS_DIST_ND, the distance between the points.
        //
    {
        int i;

        double dist = 0.0;
        for (i = 0; i < dim_num; i++)
        {
            dist += (p1[i] - p2[i]) * (p1[i] - p2[i]);
        }

        dist = Math.Sqrt(dist);

        return dist;
    }

    public static void points_hull_2d(int node_num, double[] node_xy, ref int hull_num,
            ref int[] hull)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINTS_HULL_2D computes the convex hull of a set of nodes in 2D.
        //
        //  Discussion:
        //
        //    The work involved is N*log(H), where N is the number of points, and H is
        //    the number of points that are on the hull.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
        //
        //    Output, int *HULL_NUM, the number of nodes that lie on the convex hull.
        //
        //    Output, int HULL[NODE_NUM].  The first HULL_NUM entries contain
        //    the indices of the nodes that form the convex hull, in order.
        //    These indices are 1-based, not 0-based!
        //
    {
        int i;
        double[] p_xy = new double[2];
        double[] q_xy = new double[2];
        double[] r_xy = new double[2];

        hull_num = 0;

        switch (node_num)
        {
            case < 1:
                return;
            //
            //  If NODE_NUM = 1, the hull is the node.
            //
            case 1:
                hull[hull_num] = 1;
                hull_num += 1;
                return;
            //
            //  If NODE_NUM = 2, then the convex hull is either the two distinct nodes,
            //  or possibly a single (repeated) node.
            //
            case 2:
            {
                hull[hull_num] = 1;
                hull_num += 1;

                if (!(Math.Abs(node_xy[0 + 0 * 2] - node_xy[0 + 1 * 2]) > double.Epsilon) &&
                    !(Math.Abs(node_xy[1 + 0 * 2] - node_xy[1 + 1 * 2]) > double.Epsilon))
                {
                    return;
                }

                hull[hull_num] = 2;
                hull_num += 1;

                return;
            }
        }

        //
        //  Find the leftmost point, and take the bottom-most in a tie.
        //  Call it "Q".
        //
        int q = 1;
        for (i = 2; i <= node_num; i++)
        {
            if (node_xy[0 + (i - 1) * 2] < node_xy[0 + (q - 1) * 2] ||
                Math.Abs(node_xy[0 + (i - 1) * 2] - node_xy[0 + (q - 1) * 2]) <= double.Epsilon &&
                node_xy[1 + (i - 1) * 2] < node_xy[1 + (q - 1) * 2])
            {
                q = i;
            }
        }

        q_xy[0] = node_xy[0 + (q - 1) * 2];
        q_xy[1] = node_xy[1 + (q - 1) * 2];
        //
        //  Remember the starting point.
        //
        int first = q;
        hull[hull_num] = q;
        hull_num += 1;
        //
        //  For the first point, make a dummy previous point, 1 unit south,
        //  and call it "P".
        //
        p_xy[0] = q_xy[0];
        p_xy[1] = q_xy[1] - 1.0;
        //
        //  Now, having old point P, and current point Q, find the new point R
        //  so the angle PQR is maximal.
        //
        //  Watch out for the possibility that the two nodes are identical.
        //
        for (;;)
        {
            int r = 0;
            double angle_max = 0.0;

            for (i = 1; i <= node_num; i++)
            {
                if (i == q || (!(Math.Abs(node_xy[0 + (i - 1) * 2] - q_xy[0]) > double.Epsilon) &&
                               !(Math.Abs(node_xy[1 + (i - 1) * 2] - q_xy[1]) > double.Epsilon)))
                {
                    continue;
                }

                double angle = Angle.angle_rad_2d(p_xy, q_xy, node_xy, p3Index: +(i - 1) * 2);

                if (r == 0 || angle_max < angle)
                {
                    r = i;
                    r_xy[0] = node_xy[0 + (r - 1) * 2];
                    r_xy[1] = node_xy[1 + (r - 1) * 2];
                    angle_max = angle;
                }
                //
                //  In case of ties, choose the nearer point.
                //
                else if (Math.Abs(angle - angle_max) <= double.Epsilon)
                {
                    double di = Math.Sqrt(Math.Pow(node_xy[0 + (i - 1) * 2] - q_xy[0], 2)
                                          + Math.Pow(node_xy[1 + (i - 1) * 2] - q_xy[1], 2));

                    double dr = Math.Sqrt(Math.Pow(r_xy[0] - q_xy[0], 2)
                                          + Math.Pow(r_xy[1] - q_xy[1], 2));

                    if (!(di < dr))
                    {
                        continue;
                    }

                    r = i;
                    r_xy[0] = node_xy[0 + (r - 1) * 2];
                    r_xy[1] = node_xy[1 + (r - 1) * 2];
                    angle_max = angle;
                }
            }

            //
            //  If we've returned to our starting node, exit.
            //
            if (r == first)
            {
                break;
            }

            if (node_num < hull_num + 1)
            {
                Console.WriteLine("");
                Console.WriteLine("POINTS_HULL_2D - Fatal error!");
                Console.WriteLine("  The algorithm failed.");
                return;
            }

            //
            //  Add point R to the convex hull.
            //
            hull[hull_num] = r;
            hull_num += 1;
            //
            //  Set Q := P, P := R, and repeat.
            //
            q = r;

            p_xy[0] = q_xy[0];
            p_xy[1] = q_xy[1];

            q_xy[0] = r_xy[0];
            q_xy[1] = r_xy[1];
        }

    }

    public static void points_plot(string file_name, int node_num, double[] node_xy,
            bool node_label)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINTS_PLOT plots a pointset.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string FILE_NAME, the name of the file to create.
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the nodes.
        //
        //    Input, bool NODE_LABEL, is TRUE if the nodes are to be labeled.
        //
        //  Local parameters:
        //
        //    int CIRCLE_SIZE, controls the size of the circles depicting
        //    the nodes.  Currently set to 5.  3 is pretty small, and 1 is
        //    barely visible.
        //
    {
        const int circle_size = 3;
        int delta;
        List<string> file_unit = new();
        int node;
        int x_ps;
        int x_ps_max = 576;
        int x_ps_max_clip = 594;
        int x_ps_min = 36;
        int x_ps_min_clip = 18;
        int y_ps;
        int y_ps_max = 666;
        int y_ps_max_clip = 684;
        int y_ps_min = 126;
        int y_ps_min_clip = 108;
        //
        //  We need to do some figuring here, so that we can determine
        //  the range of the data, and hence the height and width
        //  of the piece of paper.
        //
        double x_max = -typeMethods.r8_huge();
        for (node = 0; node < node_num; node++)
        {
            if (x_max < node_xy[0 + node * 2])
            {
                x_max = node_xy[0 + node * 2];
            }
        }

        double x_min = typeMethods.r8_huge();
        for (node = 0; node < node_num; node++)
        {
            if (node_xy[0 + node * 2] < x_min)
            {
                x_min = node_xy[0 + node * 2];
            }
        }

        double x_scale = x_max - x_min;

        x_max += 0.05 * x_scale;
        x_min -= 0.05 * x_scale;
        x_scale = x_max - x_min;

        double y_max = -typeMethods.r8_huge();
        for (node = 0; node < node_num; node++)
        {
            if (y_max < node_xy[1 + node * 2])
            {
                y_max = node_xy[1 + node * 2];
            }
        }

        double y_min = typeMethods.r8_huge();
        for (node = 0; node < node_num; node++)
        {
            if (node_xy[1 + node * 2] < y_min)
            {
                y_min = node_xy[1 + node * 2];
            }
        }

        double y_scale = y_max - y_min;

        y_max += 0.05 * y_scale;
        y_min -= 0.05 * y_scale;
        y_scale = y_max - y_min;

        if (x_scale < y_scale)
        {
            delta = (int) typeMethods.r8_nint((x_ps_max - x_ps_min)
                * (y_scale - x_scale) / (2.0 * y_scale));

            x_ps_max -= delta;
            x_ps_min += delta;

            x_ps_max_clip -= delta;
            x_ps_min_clip += delta;

            x_scale = y_scale;
        }
        else if (y_scale < x_scale)
        {
            delta = (int) typeMethods.r8_nint((y_ps_max - y_ps_min)
                * (x_scale - y_scale) / (2.0 * x_scale));

            y_ps_max -= delta;
            y_ps_min += delta;

            y_ps_max_clip -= delta;
            y_ps_min_clip += delta;

            y_scale = x_scale;
        }

        file_unit.Add("%!PS-Adobe-3.0 EPSF-3.0");
        file_unit.Add("%%Creator: points_plot.C");
        file_unit.Add("%%Title: " + file_name + "");

        file_unit.Add("%%Pages: 1");
        file_unit.Add("%%BoundingBox:  "
                      + x_ps_min + "  "
                      + y_ps_min + "  "
                      + x_ps_max + "  "
                      + y_ps_max + "");
        file_unit.Add("%%Document-Fonts: Times-Roman");
        file_unit.Add("%%LanguageLevel: 1");
        file_unit.Add("%%EndComments");
        file_unit.Add("%%BeginProlog");
        file_unit.Add("/inch {72 mul} def");
        file_unit.Add("%%EndProlog");
        file_unit.Add("%%Page:      1     1");
        file_unit.Add("save");
        file_unit.Add("%");
        file_unit.Add("% Set the RGB line color to very light gray.");
        file_unit.Add("%");
        file_unit.Add(" 0.9000 0.9000 0.9000 setrgbcolor");
        file_unit.Add("%");
        file_unit.Add("% Draw a gray border around the page.");
        file_unit.Add("%");
        file_unit.Add("newpath");
        file_unit.Add(x_ps_min + "  "
                               + y_ps_min + "  moveto");
        file_unit.Add(x_ps_max + "  "
                               + y_ps_min + "  lineto");
        file_unit.Add(x_ps_max + "  "
                               + y_ps_max + "  lineto");
        file_unit.Add(x_ps_min + "  "
                               + y_ps_max + "  lineto");
        file_unit.Add(x_ps_min + "  "
                               + y_ps_min + "  lineto");
        file_unit.Add("stroke");
        file_unit.Add("%");
        file_unit.Add("% Set RGB line color to black.");
        file_unit.Add("%");
        file_unit.Add(" 0.0000 0.0000 0.0000 setrgbcolor");
        file_unit.Add("%");
        file_unit.Add("%  Set the font and its size:");
        file_unit.Add("%");
        file_unit.Add("/Times-Roman findfont");
        file_unit.Add("0.50 inch scalefont");
        file_unit.Add("setfont");
        file_unit.Add("%");
        file_unit.Add("%  Print a title:");
        file_unit.Add("%");
        file_unit.Add("%  210  702 moveto");
        file_unit.Add("%(Pointset) show");
        file_unit.Add("%");
        file_unit.Add("% Define a clipping polygon");
        file_unit.Add("%");
        file_unit.Add("newpath");
        file_unit.Add(x_ps_min_clip + "  "
                                    + y_ps_min_clip + "  moveto");
        file_unit.Add(x_ps_max_clip + "  "
                                    + y_ps_min_clip + "  lineto");
        file_unit.Add(x_ps_max_clip + "  "
                                    + y_ps_max_clip + "  lineto");
        file_unit.Add(x_ps_min_clip + "  "
                                    + y_ps_max_clip + "  lineto");
        file_unit.Add(x_ps_min_clip + "  "
                                    + y_ps_min_clip + "  lineto");
        file_unit.Add("clip newpath");
        //
        //  Draw the nodes.
        //
        file_unit.Add("%");
        file_unit.Add("%  Draw filled dots at each node:");
        file_unit.Add("%");
        file_unit.Add("%  Set the color to blue:");
        file_unit.Add("%");
        file_unit.Add("0.000  0.150  0.750  setrgbcolor");
        file_unit.Add("%");

        for (node = 0; node < node_num; node++)
        {
            x_ps = (int) (
                ((x_max - node_xy[0 + node * 2]) * x_ps_min
                 + (+node_xy[0 + node * 2] - x_min) * x_ps_max)
                / (x_max - x_min));

            y_ps = (int) (
                ((y_max - node_xy[1 + node * 2]) * y_ps_min
                 + (node_xy[1 + node * 2] - y_min) * y_ps_max)
                / (y_max - y_min));

            file_unit.Add("newpath  "
                          + x_ps + "  "
                          + y_ps + "  "
                          + circle_size + " 0 360 arc closepath fill");
        }

        //
        //  Label the nodes.
        //
        file_unit.Add("%");
        file_unit.Add("%  Label the nodes:");
        file_unit.Add("%");
        file_unit.Add("%  Set the color to darker blue:");
        file_unit.Add("%");
        file_unit.Add("0.000  0.250  0.850  setrgbcolor");
        file_unit.Add("/Times-Roman findfont");
        file_unit.Add("0.20 inch scalefont");
        file_unit.Add("setfont");

        file_unit.Add("%");

        for (node = 0; node < node_num; node++)
        {
            x_ps = (int) (
                ((x_max - node_xy[0 + node * 2]) * x_ps_min
                 + (+node_xy[0 + node * 2] - x_min) * x_ps_max)
                / (x_max - x_min));

            y_ps = (int) (
                ((y_max - node_xy[1 + node * 2]) * y_ps_min
                 + (node_xy[1 + node * 2] - y_min) * y_ps_max)
                / (y_max - y_min));

            file_unit.Add("newpath  "
                          + x_ps + "  "
                          + y_ps + 5 + "  moveto ("
                          + node + 1 + ") show");
        }

        file_unit.Add("%");
        file_unit.Add("restore showpage");
        file_unit.Add("%");
        file_unit.Add("% End of page");
        file_unit.Add("%");
        file_unit.Add("%%Trailer");
        file_unit.Add("%%EOF");


        try
        {
            File.WriteAllLines(file_name, file_unit);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("POINTS_PLOT - Fatal error!");
            Console.WriteLine("  Could not open the output EPS file.");
        }
    }

    public static int points_point_near_naive_2d(int nset, double[] pset, double[] ptest,
            ref double d_min)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINTS_POINT_NEAR_NAIVE_2D finds the nearest point to a given point in 2D.
        //
        //  Discussion:
        //
        //    A naive algorithm is used.  The distance to every point is calculated,
        //    in order to determine the smallest.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 October 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NSET, the number of points in the set.
        //
        //    Input, double PSET[2*NSET], the coordinates of the points in the set.
        //
        //    Input, double PTEST[2], the point whose nearest neighbor is sought.
        //
        //    Output, double *D_MIN, the distance between P and PSET(*,I_MIN).
        //
        //    Output, int POINTS_POINT_NEAR_NAIVE_2D, the index of the nearest
        //    point in PSET to P.
        //
    {
        const int DIM_NUM = 2;

        int j;

        d_min = typeMethods.r8_huge();
        int p_min = 0;

        for (j = 0; j < nset; j++)
        {
            double d = 0.0;
            int i;
            for (i = 0; i < DIM_NUM; i++)
            {
                d += Math.Pow(ptest[i] - pset[i + j * DIM_NUM], 2);
            }

            if (!(d < d_min))
            {
                continue;
            }

            d_min = d;
            p_min = j;
        }

        d_min = Math.Sqrt(d_min);

        return p_min;
    }

    public static int points_point_near_naive_3d(int nset, double[] pset, double[] ptest,
            ref double d_min)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINTS_POINT_NEAR_NAIVE_3D finds the nearest point to a given point in 3D.
        //
        //  Discussion:
        //
        //    A naive algorithm is used.  The distance to every point is calculated,
        //    in order to determine the smallest.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 October 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NSET, the number of points in the set.
        //
        //    Input, double PSET[3*NSET], the coordinates of the points in the set.
        //
        //    Input, double PTEST[3], the point whose nearest neighbor is sought.
        //
        //    Output, double *D_MIN, the distance between P and PSET(*,I_MIN).
        //
        //    Output, int POINTS_POINT_NEAR_NAIVE_3D, the index of the nearest
        //    point in PSET to P.
        //
    {
        const int DIM_NUM = 3;

        int j;

        d_min = typeMethods.r8_huge();
        int p_min = 0;

        for (j = 0; j < nset; j++)
        {
            double d = 0.0;
            int i;
            for (i = 0; i < DIM_NUM; i++)
            {
                d += Math.Pow(ptest[i] - pset[i + j * DIM_NUM], 2);
            }

            if (!(d < d_min))
            {
                continue;
            }

            d_min = d;
            p_min = j;
        }

        d_min = Math.Sqrt(d_min);

        return p_min;
    }

    public static int points_point_near_naive_nd(int dim_num, int nset, double[] pset,
            double[] ptest, ref double d_min)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINTS_POINT_NEAR_NAIVE_ND finds the nearest point to a given point in ND.
        //
        //  Discussion:
        //
        //    A naive algorithm is used.  The distance to every point is calculated,
        //    in order to determine the smallest.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 October 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NSET, the number of points in the set.
        //
        //    Input, double PSET[DIM_NUM*NSET], the coordinates of the points in the set.
        //
        //    Input, double PTEST[DIM_NUM], the point whose nearest neighbor is sought.
        //
        //    Output, double *D_MIN, the distance between P and PSET(*,I_MIN).
        //
        //    Output, int POINTS_POINT_NEAR_NAIVE_ND, the index of the nearest
        //    point in PSET to P.
        //
    {
        int j;

        d_min = typeMethods.r8_huge();
        int p_min = 0;

        for (j = 0; j < nset; j++)
        {
            double d = 0.0;
            int i;
            for (i = 0; i < dim_num; i++)
            {
                d += Math.Pow(ptest[i] - pset[i + j * dim_num], 2);
            }

            if (!(d < d_min))
            {
                continue;
            }

            d_min = d;
            p_min = j;
        }

        d_min = Math.Sqrt(d_min);

        return p_min;
    }

    public static int[] points_points_near_naive_2d(int nset, double[] pset, int ntest,
            double[] ptest)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINTS_POINTS_NEAR_NAIVE_2D finds the nearest point to given points in 2D.
        //
        //  Discussion:
        //
        //    A naive algorithm is used.  The distance to every point is calculated,
        //    in order to determine the smallest.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 January 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NSET, the number of points in the set.
        //
        //    Input, double PSET[2*NSET], the coordinates of the points in the set.
        //
        //    Input, int NTEST, the number of test points.
        //
        //    Input, double PTEST[2*NTEST], the coordinates of the test points.
        //
        //    Output, int POINTS_POINTS_NEAR_NAIVE_2D[NTEST], the index of the
        //    nearest point in PSET to each point in PTEST.
        //
    {
        const int DIM_NUM = 2;

        int test;

        int[] nearest = new int[ntest];

        for (test = 0; test < ntest; test++)
        {
            double d_min = typeMethods.r8_huge();
            nearest[test] = -1;

            int set;
            for (set = 0; set < nset; set++)
            {
                double d = 0.0;
                int i;
                for (i = 0; i < DIM_NUM; i++)
                {
                    d += Math.Pow(ptest[i + test * DIM_NUM] - pset[i + set * DIM_NUM], 2);
                }

                if (!(d < d_min))
                {
                    continue;
                }

                d_min = d;
                nearest[test] = set;
            }
        }

        return nearest;
    }

    public static int[] points_points_near_naive_3d(int nset, double[] pset, int ntest,
            double[] ptest)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINTS_POINTS_NEAR_NAIVE_3D finds the nearest point to given points in 3D.
        //
        //  Discussion:
        //
        //    A naive algorithm is used.  The distance to every point is calculated,
        //    in order to determine the smallest.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 January 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NSET, the number of points in the set.
        //
        //    Input, double PSET[3*NSET], the coordinates of the points in the set.
        //
        //    Input, int NTEST, the number of test points.
        //
        //    Input, double PTEST[3*NTEST], the coordinates of the test points.
        //
        //    Output, int POINTS_POINTS_NEAR_NAIVE_3D[NTEST], the index of the
        //    nearest point in PSET to each point in PTEST.
        //
    {
        const int DIM_NUM = 3;

        int test;

        int[] nearest = new int[ntest];

        for (test = 0; test < ntest; test++)
        {
            double d_min = typeMethods.r8_huge();
            nearest[test] = -1;

            int set;
            for (set = 0; set < nset; set++)
            {
                double d = 0.0;
                int i;
                for (i = 0; i < DIM_NUM; i++)
                {
                    d += Math.Pow(ptest[i + test * DIM_NUM] - pset[i + set * DIM_NUM], 2);
                }

                if (!(d < d_min))
                {
                    continue;
                }

                d_min = d;
                nearest[test] = set;
            }
        }

        return nearest;
    }
}