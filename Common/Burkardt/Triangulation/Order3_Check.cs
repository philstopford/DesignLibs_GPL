using System;
using Burkardt.Types;

namespace Burkardt.TriangulationNS
{
    public static partial class Check
    {
       public static int triangulation_order3_check(int node_num, int triangle_num,
                int[] triangle_node)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGULATION_ORDER3_CHECK makes some simple checks on a triangulation.
            //
            //  Discussion:
            //
            //    Because this routine does not receive the physical coordinates of
            //    the nodes, it cannot ensure that the triangulation is maximal,
            //    that is, that no more triangles can be created.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 June 2005
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
            //    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up the
            //    triangles.  These should be listed in counterclockwise order.
            //
            //    Output, int TRIANGULATION_CHECK, error flag.
            //    0, no error occurred.
            //    nonzero, an error occurred, the triangulation is not valid.
            //
        {
            int boundary_num;
            int euler;
            int i;
            int j;
            int[] used;
            //
            //  Checks 1 and 2:
            //  node_num must be at least 3.
            //  TRIANGLE_NUM must be at least 1.
            //
            if (node_num < 3)
            {
                Console.WriteLine("");
                Console.WriteLine("TRIANGULATION_ORDER3_CHECK - Fatal error!");
                Console.WriteLine("  The number of nodes is less than 3!");
                return 1;
            }

            if (triangle_num < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("TRIANGULATION_ORDER3_CHECK - Fatal error!");
                Console.WriteLine("  The number of triangles is less than 1!");
                return 2;
            }

            //
            //  Checks 3 and 4:
            //  Verify that all node values are greater than or equal to 1
            //  and less than or equal to node_num.
            //
            for (j = 0; j < triangle_num; j++)
            {
                for (i = 0; i < 3; i++)
                {
                    if (triangle_node[i + j * 3] < 1)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("TRIANGULATION_ORDER3_CHECK - Fatal error!");
                        Console.WriteLine("  Some vertices are less than 1!");
                        return 3;
                    }
                }
            }

            for (j = 0; j < triangle_num; j++)
            {
                for (i = 0; i < 3; i++)
                {
                    if (node_num < triangle_node[i + j * 3])
                    {
                        Console.WriteLine("");
                        Console.WriteLine("TRIANGULATION_ORDER3_CHECK - Fatal error!");
                        Console.WriteLine("  Some vertices are greater than node_num!");
                        return 4;
                    }
                }
            }

            //
            //  Check 5:
            //  Verify that every node is used at least once.
            //
            used = new int[node_num];

            for (i = 0; i < node_num; i++)
            {
                used[i] = 0;
            }

            for (j = 0; j < triangle_num; j++)
            {
                for (i = 0; i < 3; i++)
                {
                    used[triangle_node[i + j * 3] - 1] = used[triangle_node[i + j * 3] - 1] + 1;
                }
            }

            for (i = 0; i < node_num; i++)
            {
                if (used[i] == 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("TRIANGULATION_ORDER3_CHECK - Fatal error!");
                    Console.WriteLine("  Some nodes are never used as triangle vertices!");
                    Console.WriteLine("  First example is node " + (i + 1) + "");
                    return 5;
                }
            }

            //
            //  Check 6:
            //  Verify that no node is repeated in a triangle.
            //
            for (j = 0; j < triangle_num; j++)
            {
                if (triangle_node[0 + j * 3] == triangle_node[1 + j * 3] ||
                    triangle_node[1 + j * 3] == triangle_node[2 + j * 3] ||
                    triangle_node[2 + j * 3] == triangle_node[0 + j * 3])
                {
                    Console.WriteLine("");
                    Console.WriteLine("TRIANGULATION_ORDER3_CHECK - Fatal error!");
                    Console.WriteLine("  A triangle contains a null edge!");
                    return 6;
                }
            }

            //
            //  Check 7:
            //  Verify that no edge is repeated, and that repeated edges occur in
            //  negated pairs.
            //
            boundary_num = triangulation_order3_edge_check(triangle_num,
                triangle_node);

            if (boundary_num < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("TRIANGULATION_ORDER3_CHECK - Fatal error!");
                Console.WriteLine("  Some edges are repeated or given in the wrong direction!");
                return 7;
            }

            //
            //  Check 8:
            //  Does the triangulation satisfy Euler's criterion?
            //  If not, then the triangulation is not proper.  (For instance, there
            //  might be a hole in the interior.)
            //
            euler = boundary_num + triangle_num + 2 - 2 * node_num;

            if (euler != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("TRIANGULATION_ORDER3_CHECK - Fatal error!");
                Console.WriteLine("  The triangulation does not satisfy Euler's criterion!");
                return 8;
            }

            return 0;
        }

        public static int triangulation_order3_edge_check(int triangle_num, int[] triangle_node)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGULATION_ORDER3_EDGE_CHECK checks the edges of a triangulation.
            //
            //  Discussion:
            //
            //    Converted from a row-based to a column-based calculation.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int TRIANGLE_NUM, the number of triangles.
            //
            //    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up
            //    each triangle.
            //
            //    Output, int TRIANGULATION_EDGE_CHECK is negative if an error was
            //    detected; otherwise, it is the number of edges that lie on the boundary.
            //
        {
            int boundary_num;
            int i;
            int j;
            int k;
            int[] col;
            int tri;
            int triangle_order = 3;
            //
            //  Step 1.
            //  From the list of nodes for triangle T, of the form: (I,J,K)
            //  construct the three neighbor relations:
            //
            //    (I,J,+1) or (J,I,-1),
            //    (J,K,+1) or (K,J,-1),
            //    (K,I,+1) or (I,K,-1)
            //
            //  where we choose (I,J,+1) if I < J, or else (J,I,-1) and so on.
            //
            col = new int[3 * (3 * triangle_num)];

            for (tri = 0; tri < triangle_num; tri++)
            {
                i = triangle_node[0 + tri * triangle_order];
                j = triangle_node[1 + tri * triangle_order];
                k = triangle_node[2 + tri * triangle_order];

                if (i < j)
                {
                    col[0 + (3 * tri + 0) * 3] = i;
                    col[1 + (3 * tri + 0) * 3] = j;
                    col[2 + (3 * tri + 0) * 3] = +1;
                }
                else
                {
                    col[0 + (3 * tri + 0) * 3] = j;
                    col[1 + (3 * tri + 0) * 3] = i;
                    col[2 + (3 * tri + 0) * 3] = -1;
                }

                if (j < k)
                {
                    col[0 + (3 * tri + 1) * 3] = j;
                    col[1 + (3 * tri + 1) * 3] = k;
                    col[2 + (3 * tri + 1) * 3] = +1;
                }
                else
                {
                    col[0 + (3 * tri + 1) * 3] = k;
                    col[1 + (3 * tri + 1) * 3] = j;
                    col[2 + (3 * tri + 1) * 3] = -1;
                }

                if (k < i)
                {
                    col[0 + (3 * tri + 2) * 3] = k;
                    col[1 + (3 * tri + 2) * 3] = i;
                    col[2 + (3 * tri + 2) * 3] = +1;
                }
                else
                {
                    col[0 + (3 * tri + 2) * 3] = i;
                    col[1 + (3 * tri + 2) * 3] = k;
                    col[2 + (3 * tri + 2) * 3] = -1;
                }
            }

            //
            //  Step 2. Perform an ascending dictionary sort on the neighbor relations.
            //
            typeMethods.i4col_sort_a(3, 3 * triangle_num, ref col);
            //
            //  Step 3.
            //
            //  If any record occurs twice, we have an error.
            //  Unpaired records lie on the convex hull.
            //
            i = 0;
            boundary_num = 0;

            while (i < 3 * triangle_num)
            {
                i = i + 1;

                if (i == 3 * triangle_num)
                {
                    boundary_num = boundary_num + 1;
                }
                else
                {
                    if (col[0 + (i - 1) * 3] == col[0 + i * 3] &&
                        col[1 + (i - 1) * 3] == col[1 + i * 3])
                    {
                        if (col[2 + (i - 1) * 3] == col[2 + i * 3])
                        {
                            Console.WriteLine("");
                            Console.WriteLine("TRIANGULATION_ORDER3_EDGE_CHECK - Warning!");
                            Console.WriteLine("  An edge occurs twice.");
                            boundary_num = -1;
                            return boundary_num;
                        }
                        else
                        {
                            i = i + 1;
                        }
                    }
                    else
                    {
                        boundary_num = boundary_num + 1;
                    }
                }
            }
            return boundary_num;
        }
    }
}