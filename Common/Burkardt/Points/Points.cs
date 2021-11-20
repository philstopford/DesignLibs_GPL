using System;

namespace Burkardt.PointsNS;

public static partial class Points
{
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

                if (Math.Abs(node_xy[0 + 0 * 2] - node_xy[0 + 1 * 2]) > double.Epsilon || Math.Abs(node_xy[1 + 0 * 2] - node_xy[1 + 1 * 2]) > double.Epsilon)
                {
                    hull[hull_num] = 2;
                    hull_num += 1;
                }

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
                if (i != q && (Math.Abs(node_xy[0 + (i - 1) * 2] - q_xy[0]) > double.Epsilon || Math.Abs(node_xy[1 + (i - 1) * 2] - q_xy[1]) > double.Epsilon))
                {
                    double angle = Helpers.angle_rad_2d(p_xy, q_xy, node_xy, p3Index: +(i - 1) * 2);

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

                        if (di < dr)
                        {
                            r = i;
                            r_xy[0] = node_xy[0 + (r - 1) * 2];
                            r_xy[1] = node_xy[1 + (r - 1) * 2];
                            angle_max = angle;
                        }
                    }
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
                Console.WriteLine("\n");
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

    public static int[] points_delaunay_naive_2d(int node_num, double[] node_xy,
            ref int triangle_num )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINTS_DELAUNAY_NAIVE_2D computes the Delaunay triangulation in 2D.
        //
        //  Discussion:
        //
        //    A naive and inefficient (but extremely simple) method is used.
        //
        //    This routine is only suitable as a demonstration code for small
        //    problems.  Its running time is of order NODE_NUM^4.  Much faster
        //    algorithms are available.
        //
        //    Given a set of nodes in the plane, a triangulation is a set of
        //    triples of distinct nodes, forming triangles, so that every
        //    point with the convex hull of the set of  nodes is either one
        //    of the nodes, or lies on an edge of one or more triangles,
        //    or lies within exactly one triangle.
        //
        //    The number of nodes must be at least 3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Joseph ORourke,
        //    Computational Geometry,
        //    Cambridge University Press,
        //    Second Edition, 1998, page 187.
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
        //
        //    Output, int *TRIANGLE_NUM, the number of triangles.
        //
        //    Output, int POINTS_DELAUNAY_NAIVE_2D[3*TRIANGLE_NUM], the indices of the
        //    nodes making each triangle.
        //
    {
        int i;
        int pass;
        int[] tri = new int[3];

        int count = 0;

        double[] z = new double [node_num];

        for (i = 0; i < node_num; i++)
        {
            z[i] = node_xy[0 + i * 2] * node_xy[0 + i * 2] + node_xy[1 + i * 2] * node_xy[1 + i * 2];
        }

        //
        //  First pass counts triangles,
        //  Second pass allocates triangles and sets them.
        //
        for (pass = 1; pass <= 2; pass++)
        {
            tri = pass switch
            {
                2 => new int[3 * count],
                _ => tri
            };

            count = 0;
            //
            //  For each triple (I,J,K):
            //
            for (i = 0; i < node_num - 2; i++)
            {
                int j;
                for (j = i + 1; j < node_num; j++)
                {
                    int k;
                    for (k = i + 1; k < node_num; k++)
                    {
                        if (j != k)
                        {
                            double xn = (node_xy[1 + j * 2] - node_xy[1 + i * 2]) * (z[k] - z[i])
                                        - (node_xy[1 + k * 2] - node_xy[1 + i * 2]) * (z[j] - z[i]);
                            double yn = (node_xy[0 + k * 2] - node_xy[0 + i * 2]) * (z[j] - z[i])
                                        - (node_xy[0 + j * 2] - node_xy[0 + i * 2]) * (z[k] - z[i]);
                            double zn = (node_xy[0 + j * 2] - node_xy[0 + i * 2])
                                        * (node_xy[1 + k * 2] - node_xy[1 + i * 2])
                                        - (node_xy[0 + k * 2] - node_xy[0 + i * 2])
                                        * (node_xy[1 + j * 2] - node_xy[1 + i * 2]);

                            bool flag = zn < 0;

                            switch (flag)
                            {
                                case true:
                                {
                                    int m;
                                    for (m = 0; m < node_num; m++)
                                    {
                                        flag = flag && (node_xy[0 + m * 2] - node_xy[0 + i * 2]) * xn
                                            + (node_xy[1 + m * 2] - node_xy[1 + i * 2]) * yn
                                            + (z[m] - z[i]) * zn <= 0;
                                    }

                                    break;
                                }
                            }

                            switch (flag)
                            {
                                case true:
                                {
                                    switch (pass)
                                    {
                                        case 2:
                                            tri[0 + count * 3] = i + 1;
                                            tri[1 + count * 3] = j + 1;
                                            tri[2 + count * 3] = k + 1;
                                            break;
                                    }

                                    count += 1;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }

        triangle_num = count;
 
        return tri;
    }
}