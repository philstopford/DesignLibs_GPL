using System;
using Burkardt.Types;

namespace Burkardt.TriangulationNS;

public static partial class Print
{
    public static void triangulation_order3_print(int node_num, int triangle_num,
            double[] node_xy, int[] triangle_node, int[] triangle_neighbor )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_ORDER3_PRINT prints information defining a triangulation.
        //
        //  Discussion:
        //
        //    Triangulations created by R8TRIS2 include extra information encoded
        //    in the negative values of TRIANGLE_NEIGHBOR.
        //
        //    Because some of the nodes counted in NODE_NUM may not actually be
        //    used in the triangulation, I needed to compute the true number
        //    of vertices.  I added this calculation on 13 October 2001.
        //
        //    Ernest Fasse pointed out an error in the indexing of VERTEX_LIST,
        //    which was corrected on 19 February 2004.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int TRIANGLE_NUM, the number of triangles.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
        //
        //    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up
        //    the triangles.
        //
        //    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the triangle neighbors
        //    on each side.  If there is no triangle neighbor on a particular side,
        //    the value of TRIANGLE_NEIGHBOR should be negative.  If the
        //    triangulation data was created by R8TRIS2, then there is more
        //    information encoded in the negative values.
        //
    {
        int DIM_NUM = 2;

        int boundary_num;
        int i;
        int j;
        int k;
        int n1;
        int n2;
        int s;
        int s1;
        int s2;
        bool skip;
        int t;
        int[] vertex_list;
        int vertex_num;

        Console.WriteLine("");
        Console.WriteLine("TRIANGULATION_ORDER3_PRINT");
        Console.WriteLine("  Information defining a triangulation.");
        Console.WriteLine("");
        Console.WriteLine("  The number of nodes is " + node_num + "");

        typeMethods.r8mat_transpose_print(DIM_NUM, node_num, node_xy, "  Node coordinates");

        Console.WriteLine("");
        Console.WriteLine("  The number of triangles is " + triangle_num + "");
        Console.WriteLine("");
        Console.WriteLine("  Sets of three nodes are used as vertices of");
        Console.WriteLine("  the triangles.  For each triangle, the nodes");
        Console.WriteLine("  are listed in counterclockwise order.");

        typeMethods.i4mat_transpose_print(3, triangle_num, triangle_node, "  Triangle nodes");

        Console.WriteLine("");
        Console.WriteLine("  On each side of a given triangle, there is either");
        Console.WriteLine("  another triangle, or a piece of the convex hull.");
        Console.WriteLine("  For each triangle, we list the indices of the three");
        Console.WriteLine("  neighbors, or (if negative) the codes of the");
        Console.WriteLine("  segments of the convex hull.");

        typeMethods.i4mat_transpose_print(3, triangle_num, triangle_neighbor,
            "  Triangle neighbors");
        //
        //  Determine VERTEX_NUM, the number of vertices.
        //
        vertex_list = new int[3 * triangle_num];

        k = 0;
        for (t = 0; t < triangle_num; t++)
        {
            for (s = 0; s < 3; s++)
            {
                vertex_list[k] = triangle_node[s + t * 3];
                k += 1;
            }
        }

        typeMethods.i4vec_sort_heap_a(3 * triangle_num, ref vertex_list);

        vertex_num = typeMethods.i4vec_sorted_unique(3 * triangle_num, ref vertex_list);

        //
        //  Determine the number of boundary points.
        //
        boundary_num = 2 * vertex_num - triangle_num - 2;

        Console.WriteLine("");
        Console.WriteLine("  The number of boundary points is " + boundary_num + "");
        Console.WriteLine("");
        Console.WriteLine("  The segments that make up the convex hull can be");
        Console.WriteLine("  determined from the negative entries of the triangle");
        Console.WriteLine("  neighbor list.");
        Console.WriteLine("");
        Console.WriteLine("     #   Tri  Side    N1    N2");
        Console.WriteLine("");

        skip = false;

        k = 0;

        for (i = 0; i < triangle_num; i++)
        {
            for (j = 0; j < 3; j++)
            {
                if (triangle_neighbor[j + i * 3] < 0)
                {
                    s = -triangle_neighbor[j + i * 3];
                    t = s / 3;

                    if (t < 1 || triangle_num < t)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("  Sorry, this data does not use the R8TRIS2");
                        Console.WriteLine("  convention for convex hull segments.");
                        skip = true;
                        break;
                    }

                    s1 = s % 3 + 1;
                    s2 = typeMethods.i4_wrap(s1 + 1, 1, 3);
                    k += 1;
                    n1 = triangle_node[s1 - 1 + (t - 1) * 3];
                    n2 = triangle_node[s2 - 1 + (t - 1) * 3];
                    Console.WriteLine("  "
                                      + k.ToString().PadLeft(4) + "  "
                                      + t.ToString().PadLeft(4) + "  "
                                      + s1.ToString().PadLeft(4) + "  "
                                      + n1.ToString().PadLeft(4) + "  "
                                      + n2.ToString().PadLeft(4) + "");
                }
            }

            if (skip)
            {
                break;
            }
        }
    }
}