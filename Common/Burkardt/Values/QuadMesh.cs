using System;
using Burkardt.Types;

namespace Burkardt.Values
{
    public static class QuadMesh
    {
        public static void example1_q4_mesh(int node_num, int element_num, ref double[] node_xy,
                ref int[] element_node, ref int[] element_neighbor)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EXAMPLE1_Q4_MESH sets up example #1 Q4 mesh.
            //
            //  Discussion:
            //
            //    The appropriate values of NODE_NUM and ELEMENT_NUM can be found by
            //    calling EXAMPLE1_Q4_MESH_SIZE first.
            //
            //   23---24---25---26---27
            //    | 13 | 14 | 15 | 16 |
            //   17---18---19---20---21---22
            //    |  9 | -2 | 10 | 11 | 12 |
            //   11---12---13---14---15---16
            //    |  4 |  5 |  6 |  7 |  8 |
            //    5----6----7----8----9---10
            //    |  0 |  1 |  2 |  3 |
            //    0----1----2----3----4
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 April 2020
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of nodes.  
            //
            //    Input, int ELEMENT_NUM, the number of elements.  
            //
            //    Output, double NODE_XY[2*NODE_NUM], the coordinates of the
            //    nodes.
            //
            //    Output, int ELEMENT_NODE[4*ELEMENT_NUM], the nodes
            //    that make up the elements.
            //
            //    Output, int ELEMENT_NEIGHBOR[4*ELEMENT_NUM], the 
            //    element neighbors on each side.  Negative values indicate edges that 
            //    lie on the exterior.
            //
        {
            int ELEMENT_NUM_DATA = 17;
            int NODE_NUM_DATA = 28;

            int[] element_neighbor_data =
            {
                -1, 1, 4, -1,
                -1, 2, 5, 0,
                -1, 3, 6, 1,
                -1, -1, 7, 2,
                0, 5, 9, -1,
                1, 6, -2, 4,
                2, 7, 10, 5,
                3, 8, 11, 6,
                -1, -1, 12, 7,
                4, -2, 13, -1,
                6, 11, 15, -2,
                7, 12, 16, 10,
                8, -1, -1, 11,
                9, 14, -1, -1,
                -2, 15, -1, 13,
                10, 16, -1, 14,
                11, -1, -1, 15
            };

            int[] element_node_data =
            {
                0, 1, 6, 5,
                1, 2, 7, 6,
                2, 3, 8, 7,
                3, 4, 9, 8,
                5, 6, 12, 11,
                6, 7, 13, 12,
                7, 8, 14, 13,
                8, 9, 15, 14,
                9, 10, 16, 15,
                11, 12, 18, 17,
                13, 14, 20, 19,
                14, 15, 21, 20,
                15, 16, 22, 21,
                17, 18, 24, 23,
                18, 19, 25, 24,
                19, 20, 26, 25,
                20, 21, 27, 26
            };

            double[] node_xy_data =
            {
                0.0, 0.0,
                1.0, 0.0,
                2.0, 0.0,
                3.0, 0.0,
                4.0, 0.0,
                0.0, 1.0,
                1.0, 1.0,
                2.0, 1.0,
                3.0, 1.0,
                4.0, 1.0,
                5.0, 1.0,
                0.0, 2.0,
                1.0, 2.0,
                2.0, 2.0,
                3.0, 2.0,
                4.0, 2.0,
                5.0, 2.0,
                0.0, 3.0,
                1.0, 3.0,
                2.0, 3.0,
                3.0, 3.0,
                4.0, 3.0,
                5.0, 3.0,
                0.0, 4.0,
                1.0, 4.0,
                2.0, 4.0,
                3.0, 4.0,
                4.0, 4.0
            };

            typeMethods.i4mat_copy(4, element_num, element_neighbor_data, ref element_neighbor);

            typeMethods.i4mat_copy(4, element_num, element_node_data, ref element_node);

            typeMethods.r8mat_copy(2, node_num, node_xy_data, ref node_xy);

        }

        public static void example1_q4_mesh_size(ref int node_num, ref int element_num, ref int hole_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EXAMPLE1_Q4_MESH_SIZE sets sizes for example #1 Q4 mesh
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 February 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, int *NODE_NUM, the number of nodes.  
            //
            //    Output, int *ELEMENT_NUM, the number of elements. 
            //
            //    Output, int *HOLE_NUM, the number of holes.
            //
        {
            element_num = 17;
            hole_num = 1;
            node_num = 28;

        }

        public static void example2_q4_mesh(int node_num, int element_num, ref double[] node_xy,
                ref int[] element_node, ref int[] element_neighbor)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EXAMPLE2_Q4_MESH sets up example #2 Q4 mesh.
            //
            //  Discussion:
            //
            //    The region is a semicircle.  This example includes degenerate elements
            //    (the first layer of elements is touching the origin, and so has a side
            //    of length zero).  The elements are not parallelograms.  And the elements
            //    vary in size.
            //
            //    Because of the treatment of node 1, algorithms for counting boundary 
            //    edges may become "confused".
            //
            //    The appropriate values of NODE_NUM and ELEMENT_NUM can be found by
            //    calling EXAMPLE1_Q4_MESH_SIZE first.
            //
            //   28---29---30---31---32---33---34---35---36
            //    | 24 | 25 | 26 | 27 | 28 | 29 | 30 | 31 |
            //   19---20---21---22---23---24---25---26---27
            //    | 16 | 17 | 18 | 19 | 20 | 21 | 22 | 23 |
            //   10---11---12---13---14---15---16---17---18
            //    |  8 |  9 | 10 | 11 | 12 | 13 | 14 | 15 |
            //    1----2----3----4----5----6----7----8----9
            //    |  0 |  1 |  2 |  3 |  4 |  5 |  6 |  7 |
            //    0----0----0----0----0----0----0----0----0
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 April 2020
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of nodes.  
            //
            //    Input, int ELEMENT_NUM, the number of elements.  
            //
            //    Output, double NODE_XY[2*NODE_NUM], the coordinates of the
            //    nodes.
            //
            //    Output, int ELEMENT_NODE[4*ELEMENT_NUM], the nodes
            //    that make up the elements.
            //
            //    Output, int ELEMENT_NEIGHBOR[4*ELEMENT_NUM], the 
            //    element neighbors on each side.  Negative values indicate edges that 
            //    lie on the exterior.
            //
        {
            double a;
            int col;
            int element;
            int k;
            
            double r;
            int row;
            //
            //  Set x and y coordinates.
            //
            k = 0;
            node_xy[0 + k * 2] = 0.0;
            node_xy[1 + k * 2] = 0.0;

            for (row = 1; row <= 4; row++)
            {
                r = (double)(row);
                for (col = 0; col <= 8; col++)
                {
                    a = (double)(8 - col) * Math.PI / 8.0;
                    k = k + 1;
                    node_xy[0 + k * 2] = r * Math.Cos(a);
                    node_xy[1 + k * 2] = r * Math.Sin(a);
                }
            }

            //
            //  Set the nodes that form each element.
            //
            element = 0;
            for (row = 0; row <= 3; row++)
            {
                for (col = 0; col <= 7; col++)
                {
                    if (row == 0)
                    {
                        element_node[0 + element * 4] = 0;
                        element_node[1 + element * 4] = 0;
                        element_node[2 + element * 4] = col + 2;
                        element_node[3 + element * 4] = col + 1;
                    }
                    else
                    {
                        element_node[0 + element * 4] = element_node[3 + (element - 8) * 4];
                        element_node[1 + element * 4] = element_node[2 + (element - 8) * 4];
                        element_node[2 + element * 4] = element_node[1 + element * 4] + 9;
                        element_node[3 + element * 4] = element_node[0 + element * 4] + 9;
                    }

                    element = element + 1;
                }
            }

            element = 0;
            for (row = 0; row <= 3; row++)
            {
                for (col = 0; col <= 7; col++)
                {
                    if (row == 0)
                    {
                        element_neighbor[0 + element * 4] = -1;
                    }
                    else
                    {
                        element_neighbor[0 + element * 4] = element - 8;
                    }

                    if (col == 7)
                    {
                        element_neighbor[1 + element * 4] = -1;
                    }
                    else
                    {
                        element_neighbor[1 + element * 4] = element + 1;
                    }

                    if (row == 3)
                    {
                        element_neighbor[2 + element * 4] = -1;
                    }
                    else
                    {
                        element_neighbor[2 + element * 4] = element + 8;
                    }

                    if (col == 0)
                    {
                        element_neighbor[3 + element * 4] = -1;
                    }
                    else
                    {
                        element_neighbor[3 + element * 4] = element - 1;
                    }

                    element = element + 1;
                }
            }

            return;
        }

        public static void example2_q4_mesh_size(ref int node_num, ref int element_num, ref int hole_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EXAMPLE2_Q4_MESH_SIZE sets sizes for example #2 Q4 mesh
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 February 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, int *NODE_NUM, the number of nodes.  
            //
            //    Output, int *ELEMENT_NUM, the number of elements. 
            //
            //    Output, int *HOLE_NUM, the number of holes.
            //
        {
            element_num = 32;
            hole_num = 0;
            node_num = 37;

        }
    }
}