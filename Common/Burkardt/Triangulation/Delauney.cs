using System;
using Burkardt.Types;

namespace Burkardt.TriangulationNS
{
    public static class Delauney
    {
        public static int dtris2(int point_num, int base_, ref double[] point_xy, ref int tri_num,
                ref int[] tri_vert, ref int[] tri_nabe)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DTRIS2 constructs a Delaunay triangulation of 2D vertices.
            //
            //  Discussion:
            //
            //    The routine constructs the Delaunay triangulation of a set of 2D vertices
            //    using an incremental approach and diagonal edge swaps.  Vertices are
            //    first sorted in lexicographically increasing (X,Y) order, and
            //    then are inserted one at a time from outside the convex hull.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    03 June 2009
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Barry Joe.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Barry Joe,
            //    GEOMPACK - a software package for the generation of meshes
            //    using geometric algorithms,
            //    Advances in Engineering Software,
            //    Volume 13, pages 325-331, 1991.
            //
            //  Parameters:
            //
            //    Input, int POINT_NUM, the number of vertices.
            //
            //    Input, int BASE, the base for the indexing of TRI_VERT.
            //    0, use 0-based indexing.
            //    1, use 1-based indexing.
            //
            //    Input/output, double POINT_XY[POINT_NUM*2], the coordinates of the vertices.
            //    On output, the vertices have been sorted into dictionary order.
            //
            //    Output, int tri_num, the number of triangles in the triangulation;
            //    TRI_NUM is equal to 2*POINT_NUM - NB - 2, where NB is the number
            //    of boundary vertices.
            //
            //    Output, int TRI_VERT[TRI_NUM*3], the nodes that make up each triangle.
            //    The elements are indices of POINT_XY.  The vertices of the triangles are
            //    in counter clockwise order.
            //
            //    Output, int TRI_NABE[TRI_NUM*3], the triangle neighbor list.
            //    Positive elements are indices of TIL; negative elements are used for links
            //    of a counter clockwise linked list of boundary edges; LINK = -(3*I + J-1)
            //    where I, J = triangle, edge index; TRI_NABE[I,J] refers to
            //    the neighbor along edge from vertex J to J+1 (mod 3).
            //
            //    Output, int DTRIS2, is 0 for no error.
        {
            double cmax;
            int e;
            int error;
            int i;
            int[] indx;
            int j;
            int k;
            int l;
            int ledg;
            int lr;
            int ltri;
            int m;
            int m1;
            int m2;
            int n;
            int redg;
            int rtri;
            int[] stack;
            int t;
            double tol;
            int top;

            stack = new int[point_num];

            tol = 100.0 * typeMethods.r8_epsilon();
            //
            //  Sort the vertices by increasing (x,y).
            //
            indx = typeMethods.r82vec_sort_heap_index_a(point_num, base_, point_xy);

            typeMethods.r82vec_permute(point_num, indx, base_, ref point_xy);
            //
            //  Make sure that the data points are "reasonably" distinct.
            //
            m1 = 1;

            for (i = 2; i <= point_num; i++)
            {
                m = m1;
                m1 = i;

                k = -1;

                for (j = 0; j <= 1; j++)
                {
                    cmax = Math.Max(Math.Abs(point_xy[2 * (m - 1) + j]),
                        Math.Abs(point_xy[2 * (m1 - 1) + j]));

                    if (tol * (cmax + 1.0)
                        < Math.Abs(point_xy[2 * (m - 1) + j] - point_xy[2 * (m1 - 1) + j]))
                    {
                        k = j;
                        break;
                    }

                }

                if (k == -1)
                {
                    Console.WriteLine("");
                    Console.WriteLine("DTRIS2 - Fatal error!");
                    Console.WriteLine("  Fails for point number I = " + i + "");
                    Console.WriteLine("  M =  " + m + "");
                    Console.WriteLine("  M1 = " + m1 + "");
                    Console.WriteLine("  X,Y(M)  = " + point_xy[2 * (m - 1) + 0] + "  "
                                      + point_xy[2 * (m - 1) + 1] + "");
                    Console.WriteLine("  X,Y(M1) = " + point_xy[2 * (m1 - 1) + 0] + "  "
                                      + point_xy[2 * (m1 - 1) + 1] + "");
                    return 224;
                }

            }

            //
            //  Starting from points M1 and M2, search for a third point M that
            //  makes a "healthy" triangle (M1,M2,M)
            //
            m1 = 1;
            m2 = 2;
            j = 3;

            for (;;)
            {
                if (point_num < j)
                {
                    Console.WriteLine("");
                    Console.WriteLine("DTRIS2 - Fatal error!");
                    return 225;
                }

                m = j;

                lr = Helpers.lrline(point_xy[2 * (m - 1) + 0], point_xy[2 * (m - 1) + 1],
                    point_xy[2 * (m1 - 1) + 0], point_xy[2 * (m1 - 1) + 1],
                    point_xy[2 * (m2 - 1) + 0], point_xy[2 * (m2 - 1) + 1], 0.0);

                if (lr != 0)
                {
                    break;
                }

                j = j + 1;

            }

            //
            //  Set up the triangle information for (M1,M2,M), and for any other
            //  triangles you created because points were collinear with M1, M2.
            //
            tri_num = j - 2;

            if (lr == -1)
            {
                tri_vert[3 * 0 + 0] = m1;
                tri_vert[3 * 0 + 1] = m2;
                tri_vert[3 * 0 + 2] = m;
                tri_nabe[3 * 0 + 2] = -3;

                for (i = 2; i <= tri_num; i++)
                {
                    m1 = m2;
                    m2 = i + 1;
                    tri_vert[3 * (i - 1) + 0] = m1;
                    tri_vert[3 * (i - 1) + 1] = m2;
                    tri_vert[3 * (i - 1) + 2] = m;
                    tri_nabe[3 * (i - 1) + 0] = -3 * i;
                    tri_nabe[3 * (i - 1) + 1] = i;
                    tri_nabe[3 * (i - 1) + 2] = i - 1;

                }

                tri_nabe[3 * (tri_num - 1) + 0] = -3 * (tri_num) - 1;
                tri_nabe[3 * (tri_num - 1) + 1] = -5;
                ledg = 2;
                ltri = tri_num;
            }
            else
            {
                tri_vert[3 * 0 + 0] = m2;
                tri_vert[3 * 0 + 1] = m1;
                tri_vert[3 * 0 + 2] = m;
                tri_nabe[3 * 0 + 0] = -4;

                for (i = 2; i <= tri_num; i++)
                {
                    m1 = m2;
                    m2 = i + 1;
                    tri_vert[3 * (i - 1) + 0] = m2;
                    tri_vert[3 * (i - 1) + 1] = m1;
                    tri_vert[3 * (i - 1) + 2] = m;
                    tri_nabe[3 * (i - 2) + 2] = i;
                    tri_nabe[3 * (i - 1) + 0] = -3 * i - 3;
                    tri_nabe[3 * (i - 1) + 1] = i - 1;
                }

                tri_nabe[3 * (tri_num - 1) + 2] = -3 * (tri_num);
                tri_nabe[3 * 0 + 1] = -3 * (tri_num) - 2;
                ledg = 2;
                ltri = 1;
            }

            //
            //  Insert the vertices one at a time from outside the convex hull,
            //  determine visible boundary edges, and apply diagonal edge swaps until
            //  Delaunay triangulation of vertices (so far) is obtained.
            //
            top = 0;

            for (i = j + 1; i <= point_num; i++)
            {
                m = i;
                m1 = tri_vert[3 * (ltri - 1) + ledg - 1];

                if (ledg <= 2)
                {
                    m2 = tri_vert[3 * (ltri - 1) + ledg];
                }
                else
                {
                    m2 = tri_vert[3 * (ltri - 1) + 0];
                }

                lr = Helpers.lrline(point_xy[2 * (m - 1) + 0], point_xy[2 * (m - 1) + 1],
                    point_xy[2 * (m1 - 1) + 0], point_xy[2 * (m1 - 1) + 1],
                    point_xy[2 * (m2 - 1) + 0], point_xy[2 * (m2 - 1) + 1], 0.0);

                if (0 < lr)
                {
                    rtri = ltri;
                    redg = ledg;
                    ltri = 0;
                }
                else
                {
                    l = -tri_nabe[3 * (ltri - 1) + ledg - 1];
                    rtri = l / 3;
                    redg = (l % 3) + 1;
                }

                VBEDG.vbedg(point_xy[2 * (m - 1) + 0], point_xy[2 * (m - 1) + 1], point_num,
                    point_xy, tri_num, tri_vert, tri_nabe, ref ltri, ref ledg, ref rtri, ref redg);

                n = tri_num + 1;
                l = -tri_nabe[3 * (ltri - 1) + ledg - 1];

                for (;;)
                {
                    t = l / 3;
                    e = (l % 3) + 1;
                    l = -tri_nabe[3 * (t - 1) + e - 1];
                    m2 = tri_vert[3 * (t - 1) + e - 1];

                    if (e <= 2)
                    {
                        m1 = tri_vert[3 * (t - 1) + e];
                    }
                    else
                    {
                        m1 = tri_vert[3 * (t - 1) + 0];
                    }

                    tri_num = tri_num + 1;
                    tri_nabe[3 * (t - 1) + e - 1] = tri_num;
                    tri_vert[3 * (tri_num - 1) + 0] = m1;
                    tri_vert[3 * (tri_num - 1) + 1] = m2;
                    tri_vert[3 * (tri_num - 1) + 2] = m;
                    tri_nabe[3 * (tri_num - 1) + 0] = t;
                    tri_nabe[3 * (tri_num - 1) + 1] = tri_num - 1;
                    tri_nabe[3 * (tri_num - 1) + 2] = tri_num + 1;
                    top = top + 1;

                    if (point_num < top)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("DTRIS2 - Fatal error!");
                        Console.WriteLine("  Stack overflow.");
                        return 8;
                    }

                    stack[top - 1] = tri_num;

                    if (t == rtri && e == redg)
                    {
                        break;
                    }

                }

                tri_nabe[3 * (ltri - 1) + ledg - 1] = -3 * n - 1;
                tri_nabe[3 * (n - 1) + 1] = -3 * (tri_num) - 2;
                tri_nabe[3 * (tri_num - 1) + 2] = -l;
                ltri = n;
                ledg = 2;

                error = typeMethods.swapec(m, ref top, ref ltri, ref ledg, point_num, point_xy, tri_num,
                    ref tri_vert, ref tri_nabe, stack);

                if (error != 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("DTRIS2 - Fatal error!");
                    Console.WriteLine("  Error return from SWAPEC.");
                    return error;
                }

            }

            //
            //  Now account for the sorting that we did.
            //
            for (i = 0; i < 3; i++)
            {
                for (j = 0; j < tri_num; j++)
                {
                    tri_vert[i + j * 3] = indx[tri_vert[i + j * 3] - 1];
                }
            }

            typeMethods.perm_inverse(point_num, ref indx);

            typeMethods.r82vec_permute(point_num, indx, base_, ref point_xy);

            return 0;
        }


        public static double triangulation_delaunay_discrepancy_compute(int node_num,
            double[] node_xy, int triangle_order, int triangle_num, int[] triangle_node,
        int[] triangle_neighbor, ref double angle_min, ref int angle_min_triangle,
        ref double angle_max, ref int angle_max_triangle )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_DELAUNAY_DISCREPANCY_COMPUTE reports if a triangulation is Delaunay.
        //
        //  Discussion:
        //
        //    A (maximal) triangulation is Delaunay if and only if it is locally
        //    Delaunay.
        //
        //    A triangulation is Delaunay if the minimum angle over all triangles
        //    in the triangulation is maximized.  That is, there is no other
        //    triangulation of the points which has a larger minimum angle.
        //
        //    A triangulation is locally Delaunay if, for every pair of triangles
        //    that share an edge, the minimum angle in the two triangles is larger
        //    than the minimum angle in the two triangles formed by removing the
        //    common edge and joining the two opposing vertices.
        //
        //    This function examines the question of whether a given triangulation
        //    is locally Delaunay.  It does this by looking at every pair of
        //    neighboring triangles and comparing the minimum angle attained
        //    for the current triangle pair and the alternative triangle pair.
        //
        //    Let A(i,j) be the minimum angle formed by triangles T(i) and T(j),
        //    which are two triangles in the triangulation which share a common edge.
        //    Let B(I,J) be the minimum angle formed by triangles S(i) and S(j),
        //    where S(i) and S(j) are formed by removing the common edge of T(i)
        //    and T(j), and joining the opposing vertices.
        //
        //    Then the triangulation is Delaunay if B(i,j) <= A(i,j) for every
        //    pair of neighbors T(i) and T(j).
        //
        //    If A(i,j) < B(i,j) for at least one pair of neighbors, the triangulation
        //    is not a Delaunay triangulation.
        //
        //    This program returns VALUE = min ( A(i,j) - B(i,j) ) over all
        //    triangle neighbors.  VALUE is scaled to be in degrees, for
        //    comprehensibility.  If VALUE is negative, then at least one pair
        //    of triangles violates the Delaunay condition, and so the entire
        //    triangulation is not a Delaunay triangulation.  If VALUE is nonnegative,
        //    then the triangulation is a Delaunay triangulation.
        //
        //    It is useful to return VALUE, rather than a simple True/False value,
        //    because there can be cases where the Delaunay condition is only
        //    "slightly" violated.  A simple example is a triangulation formed
        //    by starting with a mesh of squares and dividing each square into
        //    two triangles by choosing one of the diagonals of the square.
        //    The Delaunay discrepancy for this mesh, if computed exactly, is 0,
        //    but roundoff could easily result in discrepancies that were very
        //    slightly negative.
        //
        //    If VALUE is positive, and not very small in magnitude, then every
        //    pair of triangles in the triangulation satisfies the local Delaunay
        //    condition, and so the triangulation is a Delaunay triangulation.
        //
        //    If VALUE is negative, and not very small in magnitude, then at least
        //    one pair of triangles violates the Delaunay condition, and to a
        //    significant degree.  The triangulation is not a Delaunay triangulation.
        //
        //    If the magnitude of VALUE is very close to zero, then the triangulation
        //    is numerically ambiguous.  At least one pair of triangles violates
        //    or almost violates the condition, but no triangle violates the
        //    condition to a great extent.  The user must judge whether the
        //    violation is significant or not.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 September 2009
        //
        //  Author:
        //
        //    John Burkardt.
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
        //
        //    Input, int TRIANGLE_ORDER, the order of the triangles.
        //
        //    Input, int TRIANGLE_NUM, the number of triangles in
        //    the triangulation.
        //
        //    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM],
        //    the nodes that make up each triangle.
        //
        //    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the
        //    triangle neighbor list.
        //
        //    Output, double *ANGLE_MIN, the minimum angle that occurred in
        //    the triangulation.
        //
        //    Output, int *ANGLE_MIN_TRIANGLE, the triangle in which
        //    the minimum angle occurred.
        //
        //    Output, double *ANGLE_MAX, the maximum angle that occurred in
        //    the triangulation.
        //
        //    Output, int *ANGLE_MAX_TRIANGLE, the triangle in which
        //    the maximum angle occurred.
        //
        //    Output, double TRIANGULATION_DELAUNAY_DISCREPANCY,
        //    the minimum value of ( A(i,j) - B(i,j) ).
        //    POSITIVE indicates the triangulation is Delaunay.
        //    VERY NEAR ZERO is a numerically ambiguous case.
        //    NEGATIVE indicates the triangulation is not Delaunay.
        //
        {
            double angle_min1;
            double angle_min2;
            double[] angles1;
            double[] angles2;
            int i;
            int i1;
            int i2;
            int i3;
            int i4;
            int n1;
            int n2;
            int n3;
            int n4;
            int neighbor;
            
            double[] t = new double[2 * 3];
            int triangle1;
            int triangle2;
            double value;

            angle_max = 0.0;
            angle_max_triangle = -1;
            angle_min = Math.PI;
            angle_min_triangle = -1;
            value = 0.0;
            //
            //  Consider triangle TRIANGLE1
            //
            for (triangle1 = 0; triangle1 < triangle_num; triangle1++)
            {
                //
                //  Consider the side opposite vertex NEIGHBOR.
                //
                for (neighbor = 0; neighbor < 3; neighbor++)
                {
                    triangle2 = triangle_neighbor[neighbor + triangle1 * 3];
                    //
                    //  There might be no neighbor on side NEIGHBOR.
                    //
                    if (triangle2 < 0)
                    {
                        continue;
                    }

                    //
                    //  We only need to check a pair of triangles once.
                    //
                    if (triangle2 < triangle1)
                    {
                        continue;
                    }

                    //
                    //  List the vertices of the quadrilateral in such a way
                    //  that the nodes of triangle 1 come first.
                    //
                    //  We rely on a property of the TRIANGLE_NEIGHBOR array, namely, that
                    //  neighbor #1 is on the side opposite to vertex #1, and so on.
                    //
                    i1 = typeMethods.i4_wrap(neighbor + 2, 0, 2);
                    i2 = typeMethods.i4_wrap(neighbor, 0, 2);
                    i3 = typeMethods.i4_wrap(neighbor + 1, 0, 2);

                    n1 = triangle_node[i1 + triangle1 * triangle_order];
                    n2 = triangle_node[i2 + triangle1 * triangle_order];
                    n3 = triangle_node[i3 + triangle1 * triangle_order];
                    //
                    //  The "odd" or "opposing" node of the neighboring triangle
                    //  is the one which follows common node I3.
                    //
                    n4 = -1;
                    for (i = 0; i < 3; i++)
                    {
                        if (triangle_node[i + triangle2 * triangle_order] == n3)
                        {
                            i4 = i + 1;
                            i4 = typeMethods.i4_wrap(i4, 0, 2);
                            n4 = triangle_node[i4 + triangle2 * triangle_order];
                            break;
                        }
                    }

                    if (n4 == -1)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("TRIANGULATION_DELAUNAY_DISCREPANCY_COMPUTE - Fatal error/!");
                        Console.WriteLine("  Could not identify the fourth node.");
                        Console.WriteLine("");
                        Console.WriteLine("  Triangle1 = " + triangle1 + "");
                        string cout = "  Nodes =     ";
                        for (i = 0; i < 3; i++)
                        {
                            cout += "  " + triangle_node[i + triangle1 * triangle_order];
                        }

                        Console.WriteLine(cout);
                        cout = "  Neighbors =     ";
                        for (i = 0; i < 3; i++)
                        {
                            cout += "  " + triangle_neighbor[i + triangle1 * 3];
                        }

                        Console.WriteLine(cout);
                        Console.WriteLine("");
                        Console.WriteLine("  Neighbor index = " + neighbor + "");
                        Console.WriteLine("");
                        Console.WriteLine("  Triangle2 = " + triangle2 + "");
                        cout = "  Nodes =     ";
                        for (i = 0; i < 3; i++)
                        {
                            cout += "  " + triangle_node[i + triangle2 * triangle_order];
                        }

                        Console.WriteLine(cout);
                        cout = "  Neighbors =     ";
                        for (i = 0; i < 3; i++)
                        {
                            cout += "  " + triangle_neighbor[i + triangle2 * 3];
                        }

                        Console.WriteLine(cout);
                        return 1;
                    }

                    //
                    //  Compute the minimum angle for (I1,I2,I3) and (I1,I3,I4).
                    //
                    t[0 + 0 * 2] = node_xy[0 + n1 * 2];
                    t[1 + 0 * 2] = node_xy[1 + n1 * 2];
                    t[0 + 1 * 2] = node_xy[0 + n2 * 2];
                    t[1 + 1 * 2] = node_xy[1 + n2 * 2];
                    t[0 + 2 * 2] = node_xy[0 + n3 * 2];
                    t[1 + 2 * 2] = node_xy[1 + n3 * 2];
                    angles1 = typeMethods.triangle_angles_2d_new(t);

                    t[0 + 0 * 2] = node_xy[0 + n1 * 2];
                    t[1 + 0 * 2] = node_xy[1 + n1 * 2];
                    t[0 + 1 * 2] = node_xy[0 + n3 * 2];
                    t[1 + 1 * 2] = node_xy[1 + n3 * 2];
                    t[0 + 2 * 2] = node_xy[0 + n4 * 2];
                    t[1 + 2 * 2] = node_xy[1 + n4 * 2];
                    angles2 = typeMethods.triangle_angles_2d_new(t);

                    angle_min1 =
                        Math.Min(typeMethods.r8vec_min(3, angles1), typeMethods.r8vec_min(3, angles2));

                    if (angle_max < typeMethods.r8vec_max(3, angles1))
                    {
                        angle_max = typeMethods.r8vec_max(3, angles1);
                        angle_max_triangle = triangle1;
                    }

                    if (angle_max < typeMethods.r8vec_max(3, angles2))
                    {
                        angle_max = typeMethods.r8vec_max(3, angles2);
                        angle_max_triangle = triangle2;
                    }

                    if (typeMethods.r8vec_min(3, angles1) < angle_min)
                    {
                        angle_min = typeMethods.r8vec_min(3, angles1);
                        angle_min_triangle = triangle1;
                    }

                    if (typeMethods.r8vec_min(3, angles2) < angle_min)
                    {
                        angle_min = typeMethods.r8vec_min(3, angles2);
                        angle_min_triangle = triangle2;
                    }

                    //
                    //  Compute the minimum angle for (I1,I2,I4) and (I2,I3,I4).
                    //
                    t[0 + 0 * 2] = node_xy[0 + n1 * 2];
                    t[1 + 0 * 2] = node_xy[1 + n1 * 2];
                    t[0 + 1 * 2] = node_xy[0 + n2 * 2];
                    t[1 + 1 * 2] = node_xy[1 + n2 * 2];
                    t[0 + 2 * 2] = node_xy[0 + n4 * 2];
                    t[1 + 2 * 2] = node_xy[1 + n4 * 2];
                    angles1 = typeMethods.triangle_angles_2d_new(t);

                    t[0 + 0 * 2] = node_xy[0 + n3 * 2];
                    t[1 + 0 * 2] = node_xy[1 + n3 * 2];
                    t[0 + 1 * 2] = node_xy[0 + n3 * 2];
                    t[1 + 1 * 2] = node_xy[1 + n3 * 2];
                    t[0 + 2 * 2] = node_xy[0 + n4 * 2];
                    t[1 + 2 * 2] = node_xy[1 + n4 * 2];
                    angles2 = typeMethods.triangle_angles_2d_new(t);

                    angle_min2 =
                        Math.Min(typeMethods.r8vec_min(3, angles1), typeMethods.r8vec_min(3, angles2));

                    //
                    //  Compare this value to the current minimum.
                    //
                    value = Math.Min(value, angle_min1 - angle_min2);
                }
            }

            //
            //  Scale the results to degrees.
            //
            value = value * 180.0 / Math.PI;
            angle_max = angle_max * 180.0 / Math.PI;
            angle_min = angle_min * 180.0 / Math.PI;

            return value;
        }

        public static int r8tris2(int node_num, ref double[] node_xy, ref int triangle_num,
                ref int[] triangle_node, ref int[] triangle_neighbor)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8TRIS2 constructs a Delaunay triangulation of 2D vertices.
            //
            //  Discussion:
            //
            //    The routine constructs the Delaunay triangulation of a set of 2D vertices
            //    using an incremental approach and diagonal edge swaps.  Vertices are
            //    first sorted in lexicographically increasing (X,Y) order, and
            //    then are inserted one at a time from outside the convex hull.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 January 2004
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Barry Joe.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Barry Joe,
            //    GEOMPACK - a software package for the generation of meshes
            //    using geometric algorithms,
            //    Advances in Engineering Software,
            //    Volume 13, pages 325-331, 1991.
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Input/output, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
            //    On output, the coordinates have been sorted into dictionary order.
            //
            //    Output, int *TRIANGLE_NUM, the number of triangles in the triangulation;
            //    TRIANGLE_NUM is equal to 2*node_num - NB - 2, where NB is the number
            //    of boundary vertices.
            //
            //    Output, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up each
            //    triangle.  The elements are indices of NODE_XY.  The vertices of the
            //    triangles are in counterclockwise order.
            //
            //    Output, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the triangle neighbor list.
            //    Positive elements are indices of TIL; negative elements are used for links
            //    of a counter clockwise linked list of boundary edges; LINK = -(3*I + J-1)
            //    where I, J = triangle, edge index; TRIANGLE_NEIGHBOR[I,J] refers to
            //    the neighbor along edge from vertex J to J+1 (mod 3).
            //
            //    Output, int R8TRIS2, is 0 for no error.
        {
            int base_;
            double cmax;
            int e;
            int error;
            int i;
            int[] indx;
            int j;
            int k;
            int l;
            int ledg;
            int lr;
            int ltri;
            int m;
            int m1;
            int m2;
            int n;
            int redg;
            int rtri;
            int[] stack;
            int t;
            double tol;
            int top;

            stack = new int[node_num];

            tol = 100.0 * typeMethods.r8_epsilon();
            //
            //  Sort the vertices by increasing (x,y).
            //
            base_ = 0;

            indx = typeMethods.r82vec_sort_heap_index_a(node_num, base_, node_xy);

            typeMethods.r82vec_permute(node_num, indx, base_, ref node_xy);
            //
            //  Make sure that the nodes are "reasonably" distinct.
            //
            m1 = 1;

            for (i = 2; i <= node_num; i++)
            {
                m = m1;
                m1 = i;

                k = -1;

                for (j = 0; j <= 1; j++)
                {
                    cmax = Math.Max(Math.Abs(node_xy[2 * (m - 1) + j]),
                        Math.Abs(node_xy[2 * (m1 - 1) + j]));

                    if (tol * (cmax + 1.0)
                        < Math.Abs(node_xy[2 * (m - 1) + j] - node_xy[2 * (m1 - 1) + j]))
                    {
                        k = j;
                        break;
                    }

                }

                if (k == -1)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8TRIS2 - Fatal error!");
                    Console.WriteLine("  Fails for point number I = " + i + "");
                    Console.WriteLine("  M =  " + m + "");
                    Console.WriteLine("  M1 = " + m1 + "");
                    Console.WriteLine("  X,Y(M)  = " + node_xy[2 * (m - 1) + 0] + "  "
                        + node_xy[2 * (m - 1) + 1] + "");
                    Console.WriteLine("  X,Y(M1) = " + node_xy[2 * (m1 - 1) + 0] + "  "
                        + node_xy[2 * (m1 - 1) + 1] + "");
                    return (1);
                }

            }

            //
            //  Starting from nodes M1 and M2, search for a third point M that
            //  makes a "healthy" triangle (M1,M2,M)
            //
            m1 = 1;
            m2 = 2;
            j = 3;

            for (;;)
            {
                if (node_num < j)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8TRIS2 - Fatal error!");
                    return 225;
                }

                m = j;

                lr = Helpers.lrline(node_xy[2 * (m - 1) + 0], node_xy[2 * (m - 1) + 1],
                    node_xy[2 * (m1 - 1) + 0], node_xy[2 * (m1 - 1) + 1],
                    node_xy[2 * (m2 - 1) + 0], node_xy[2 * (m2 - 1) + 1], 0.0);

                if (lr != 0)
                {
                    break;
                }

                j = j + 1;

            }

            //
            //  Set up the triangle information for (M1,M2,M), and for any other
            //  triangles you created because nodes were collinear with M1, M2.
            //
            triangle_num = j - 2;

            if (lr == -1)
            {
                triangle_node[3 * 0 + 0] = m1;
                triangle_node[3 * 0 + 1] = m2;
                triangle_node[3 * 0 + 2] = m;
                triangle_neighbor[3 * 0 + 2] = -3;

                for (i = 2; i <= triangle_num; i++)
                {
                    m1 = m2;
                    m2 = i + 1;
                    triangle_node[3 * (i - 1) + 0] = m1;
                    triangle_node[3 * (i - 1) + 1] = m2;
                    triangle_node[3 * (i - 1) + 2] = m;
                    triangle_neighbor[3 * (i - 1) + 0] = -3 * i;
                    triangle_neighbor[3 * (i - 1) + 1] = i;
                    triangle_neighbor[3 * (i - 1) + 2] = i - 1;

                }

                triangle_neighbor[3 * (triangle_num - 1) + 0] = -3 * (triangle_num) - 1;
                triangle_neighbor[3 * (triangle_num - 1) + 1] = -5;
                ledg = 2;
                ltri = triangle_num;
            }
            else
            {
                triangle_node[3 * 0 + 0] = m2;
                triangle_node[3 * 0 + 1] = m1;
                triangle_node[3 * 0 + 2] = m;
                triangle_neighbor[3 * 0 + 0] = -4;

                for (i = 2; i <= triangle_num; i++)
                {
                    m1 = m2;
                    m2 = i + 1;
                    triangle_node[3 * (i - 1) + 0] = m2;
                    triangle_node[3 * (i - 1) + 1] = m1;
                    triangle_node[3 * (i - 1) + 2] = m;
                    triangle_neighbor[3 * (i - 2) + 2] = i;
                    triangle_neighbor[3 * (i - 1) + 0] = -3 * i - 3;
                    triangle_neighbor[3 * (i - 1) + 1] = i - 1;
                }

                triangle_neighbor[3 * (triangle_num - 1) + 2] = -3 * (triangle_num);
                triangle_neighbor[3 * 0 + 1] = -3 * (triangle_num) - 2;
                ledg = 2;
                ltri = 1;

            }

            //
            //  Insert the vertices one at a time from outside the convex hull,
            //  determine visible boundary edges, and apply diagonal edge swaps until
            //  Delaunay triangulation of vertices (so far) is obtained.
            //
            top = 0;

            for (i = j + 1; i <= node_num; i++)
            {
                m = i;
                m1 = triangle_node[3 * (ltri - 1) + ledg - 1];

                if (ledg <= 2)
                {
                    m2 = triangle_node[3 * (ltri - 1) + ledg];
                }
                else
                {
                    m2 = triangle_node[3 * (ltri - 1) + 0];
                }

                lr = Helpers.lrline(node_xy[2 * (m - 1) + 0], node_xy[2 * (m - 1) + 1],
                    node_xy[2 * (m1 - 1) + 0], node_xy[2 * (m1 - 1) + 1],
                    node_xy[2 * (m2 - 1) + 0], node_xy[2 * (m2 - 1) + 1], 0.0);

                if (0 < lr)
                {
                    rtri = ltri;
                    redg = ledg;
                    ltri = 0;
                }
                else
                {
                    l = -triangle_neighbor[3 * (ltri - 1) + ledg - 1];
                    rtri = l / 3;
                    redg = (l % 3) + 1;
                }

                VBEDG.vbedg(node_xy[2 * (m - 1) + 0], node_xy[2 * (m - 1) + 1], node_num,
                    node_xy, triangle_num, triangle_node, triangle_neighbor,
                    ref ltri, ref ledg, ref rtri, ref redg);

                n = triangle_num + 1;
                l = -triangle_neighbor[3 * (ltri - 1) + ledg - 1];

                for (;;)
                {
                    try
                    {
                        t = l / 3;
                        e = (l % 3) + 1;
                        l = -triangle_neighbor[3 * (t - 1) + e - 1];
                        m2 = triangle_node[3 * (t - 1) + e - 1];

                        if (e <= 2)
                        {
                            m1 = triangle_node[3 * (t - 1) + e];
                        }
                        else
                        {
                            m1 = triangle_node[3 * (t - 1) + 0];
                        }

                        triangle_num = triangle_num + 1;
                        triangle_neighbor[3 * (t - 1) + e - 1] = triangle_num;
                        triangle_node[3 * (triangle_num - 1) + 0] = m1;
                        triangle_node[3 * (triangle_num - 1) + 1] = m2;
                        triangle_node[3 * (triangle_num - 1) + 2] = m;
                        triangle_neighbor[3 * (triangle_num - 1) + 0] = t;
                        triangle_neighbor[3 * (triangle_num - 1) + 1] = triangle_num - 1;
                        triangle_neighbor[3 * (triangle_num - 1) + 2] = triangle_num + 1;
                        top = top + 1;

                        if (node_num < top)
                        {
                            Console.WriteLine("");
                            Console.WriteLine("R8TRIS2 - Fatal error!");
                            Console.WriteLine("  Stack overflow.");
                            return 8;
                        }

                        stack[top - 1] = triangle_num;

                        if (t == rtri && e == redg)
                        {
                            break;
                        }
                    }
                    catch
                    {
                        break;
                    }

                }

                triangle_neighbor[3 * (ltri - 1) + ledg - 1] = -3 * n - 1;
                triangle_neighbor[3 * (n - 1) + 1] = -3 * (triangle_num) - 2;
                triangle_neighbor[3 * (triangle_num - 1) + 2] = -l;
                ltri = n;
                ledg = 2;

                error = typeMethods.swapec(m, ref top, ref ltri, ref ledg, node_num, node_xy, triangle_num,
                    ref triangle_node, ref triangle_neighbor, stack);

                if (error != 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8TRIS2 - Fatal error!");
                    Console.WriteLine("  Error return from SWAPEC.");
                    return error;
                }

            }

            //
            //  Now account for the sorting that we did.
            //
            for (i = 0; i < 3; i++)
            {
                for (j = 0; j < triangle_num; j++)
                {
                    triangle_node[i + j * 3] = indx[triangle_node[i + j * 3] - 1];
                }
            }

            typeMethods.perm_inverse(node_num, ref indx);

            typeMethods.r82vec_permute(node_num, indx, base_, ref node_xy);

            return 0;
        }

        public static bool delaunay_swap_test(double[] xy)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DELAUNAY_SWAP_TEST performs the Delaunay swap test.
            //
            //  Discussion:
            //
            //    The current triangles are formed by nodes [0+2,3) and [0+3,4).
            //    if a swap is recommended, the new triangles will be [0+2,4) and [1+3,4).
            //
            //      4     ?     4
            //     . .         .|.
            //    1---3  ==>  1 | 3
            //     . .         .|.
            //      2           2
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Graham Carey,
            //    Computational Grids:
            //    Generation, Adaptation and Solution Strategies,
            //    Taylor and Francis, 1997,
            //    ISBN13: 978-1560326359,
            //    LC: QA377.C32.
            //
            //  Parameters:
            //
            //    Input, double XY[2*4], the coordinates of four points.
            //
            //    Output, bool SWAP, is TRUE if the triangles [0+2,4) and [1+3,4)
            //    are to replace triangles [0+2,3) and [0+3,4).
            //
        {
            double a;
            double b;
            double c;
            double d;
            bool swap;
            double x13;
            double x14;
            double x23;
            double x24;
            double y13;
            double y14;
            double y23;
            double y24;

            x13 = xy[0 + 0 * 2] - xy[0 + 2 * 2];
            x14 = xy[0 + 0 * 2] - xy[0 + 3 * 2];
            x23 = xy[0 + 1 * 2] - xy[0 + 2 * 2];
            x24 = xy[0 + 1 * 2] - xy[0 + 3 * 2];

            y13 = xy[1 + 0 * 2] - xy[1 + 2 * 2];
            y14 = xy[1 + 0 * 2] - xy[1 + 3 * 2];
            y23 = xy[1 + 1 * 2] - xy[1 + 2 * 2];
            y24 = xy[1 + 1 * 2] - xy[1 + 3 * 2];

            a = x13 * x23 + y13 * y23;
            b = x24 * y14 - x14 * y24;
            c = x23 * y13 - x13 * y23;
            d = x24 * x14 + y14 * y24;
            //
            //  The reference gives two initial tests before the
            //  main one.  However, there seems to be an error
            //  in at least one of these tests.  Since they are
            //  intended to avoid error in borderline cases, but
            //  instead cause real error in common cases, they are
            //  omitted for now.
            //
            // if ( 0.0 <= a && 0.0 <= d )
            // {
            //   swap = true;
            // }
            // else if ( a < d && d < 0.0 )
            // {
            //   swap = true;
            //  }
            //  else if ...

            if (a * b < c * d)
            {
                swap = true;
            }
            else
            {
                swap = false;
            }

            return swap;
        }
        
    }
}