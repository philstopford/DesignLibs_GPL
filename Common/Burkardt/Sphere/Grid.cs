using System;
using Burkardt.Types;

namespace Burkardt.SphereNS;

public static class Grid
{
    public static int sphere_grid_element_num(string code, int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_GRID_ELEMENT_NUM returns the number of elements in a sphere grid.
        //
        //  Discussion:
        //
        //    The number of elements generated will be NELEMX * NELEMY for
        //    quadrilaterals, or 2 * NELEMX * ( NELEMY - 1 ) for triangles.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string CODE, identifies the element desired.
        //    Legal values include 'Q4', 'Q9', 'Q16', 'T3', 'T6'.
        //
        //    Input, int NELEMX, NELEMY, the number of quadrilaterals along the
        //    X and Y directions.  
        //
        //    Output, int SPHERE_GRID_ELEMENT_NUM, the number of elements in the grid.
        //
    {
        int element_num = 1;

        switch (code)
        {
            case "Q4":
                element_num = sphere_grid_q4_element_num(nelemx, nelemy);
                break;
            case "Q9":
                element_num = sphere_grid_q9_element_num(nelemx, nelemy);
                break;
            case "Q16":
                element_num = sphere_grid_q16_element_num(nelemx, nelemy);
                break;
            case "T3":
                element_num = sphere_grid_t3_element_num(nelemx, nelemy);
                break;
            case "T6":
                element_num = sphere_grid_t6_element_num(nelemx, nelemy);
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("SPHERE_GRID_ELEMENT_NUM - Fatal error!");
                Console.WriteLine("  Illegal value of CODE = \"" + code + "\".");
                element_num = -1;
                break;
        }

        return element_num;
    }

    public static int sphere_grid_node_num (string code, int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_GRID_node_num = returns the number of nodes in a sphere grid.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string CODE, identifies the element desired.
        //    Legal values include 'Q4', 'Q9', 'Q16', 'T3', 'T6'.
        //
        //    Input, int NELEMX, NELEMY, the number of quadrilaterals along the
        //    X and Y directions.  
        //
        //    Output, int SPHERE_GRID_node_num, the number of elements in the grid.
        //
    {
        int node_num;

        switch (code)
        {
            case "Q":
                node_num = sphere_grid_q4_node_num (nelemx, nelemy);
                break;
            case "Q9":
                node_num = sphere_grid_q9_node_num (nelemx, nelemy);
                break;
            case "Q16":
                node_num = sphere_grid_q16_node_num (nelemx, nelemy);
                break;
            case "T3":
                node_num = sphere_grid_t3_node_num (nelemx, nelemy);
                break;
            case "T6":
                node_num = sphere_grid_t6_node_num (nelemx, nelemy);
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("SPHERE_GRID_node_num = - Fatal error!");
                Console.WriteLine("  Illegal value of CODE = \"" + code + "\".");
                node_num = -1;
                break;
        }

        return node_num;
    }

    public static int[] sphere_grid_q4_element(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_GRID_Q4_ELEMENT produces a Q4 sphere grid.
        //
        //  Discussion:
        //
        //    This would be the same as the grid in a plane, except that all the
        //    nodes along the bottom edge are identified (replaced by a single node
        //    that is the south pole) and similarly for the top edge, and the
        //    nodes on the extreme right edge are identified pairwise with those 
        //    on the extreme left edge.
        //
        //  Example:
        //
        //    Input:
        //
        //      NELEMX = 3, NELEMY = 4
        //
        //    Output:
        //
        //      ELEMENT_NODE =
        //         1,  1,  3,  2;
        //         1,  1,  4,  3;
        //         1,  1,  2,  4;
        //         2,  3,  6,  5;
        //         3,  4,  7,  6;
        //         4,  2,  5,  7;
        //         5,  6,  9,  8;
        //         6,  7, 10,  9;
        //         7,  5,  8, 10;
        //         8,  9, 11, 11;
        //         9, 10, 11, 11;
        //        10,  8, 11, 11;
        //
        //  Grid:
        //
        //   11----11----11----11
        //    |     |     |     |
        //    | E10 | E11 | E12 |
        //    |     |     |     |
        //    8-----9----10-----8
        //    |     |     |     |
        //    | E7  | E8  | E9  |
        //    |     |     |     |
        //    5-----6-----7-----5
        //    |     |     |     |
        //    | E4  | E5  | E6  |
        //    |     |     |     |
        //    2-----3-----4-----2
        //    |     |     |     |
        //    | E1  | E2  | E3  |
        //    |     |     |     |
        //    1-----1-----1-----1
        //
        //  Reference Element Q4:
        //
        //    |
        //    1  4------3
        //    |  |      |
        //    S  |      |
        //    |  |      |
        //    |  |      |
        //    0  1------2
        //    |
        //    +--0--R---1-->
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NELEMX, NELEMY, the number of elements along the
        //    X and Y directions.  The number of elements generated will be
        //    NELEMX * NELEMY.
        //
        //    Output, int SPHERE_GRID_Q4_ELEMENT[4*NELEMX*NELEMY], the nodes 
        //    that form each element.
        //
    {
        int base1;
        int base2;
        int element;
        int[] element_node;
        int element_num;
        int element_order = 4;
        int i;
        int j;

        element_num = sphere_grid_q4_element_num(nelemx, nelemy);

        element_node = new int[element_order * element_num];
        element = 0;

        for (j = 1; j <= nelemy; j++)
        {
            base1 = (j - 1) * nelemx + 2 - nelemx;

            for (i = 1; i <= nelemx; i++)
            {
                base2 = base1 + i - 1;

                element_node[0 + element * element_order] = base2;
                if (i < nelemx)
                {
                    element_node[1 + element * element_order] = base2 + 1;
                }
                else
                {
                    element_node[1 + element * element_order] = base1;
                }

                element_node[2 + element * element_order] =
                    element_node[1 + element * element_order] + nelemx;
                element_node[3 + element * element_order] =
                    element_node[0 + element * element_order] + nelemx;

                switch (j)
                {
                    case 1:
                        element_node[0 + element * element_order] = 1;
                        element_node[1 + element * element_order] = 1;
                        break;
                    default:
                    {
                        if (j == nelemy)
                        {
                            element_node[2 + element * element_order] = base1 + nelemx;
                            element_node[3 + element * element_order] = base1 + nelemx;
                        }

                        break;
                    }
                }

                element += 1;
            }
        }

        return element_node;
    }

    public static int sphere_grid_q4_element_num(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_GRID_Q4_ELEMENT_NUM counts the elements in a Q4 sphere grid.
        //
        //  Example:
        //
        //    Input:
        //
        //      NELEMX = 3, NELEMY = 2
        //
        //    Output:
        //
        //      ELEMENT_NUM = NELEMX * NELEMY = 6
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NELEMX, NELEMY, the number of elements along the
        //    X and Y directions. 
        //
        //    Output, int SPHERE_GRID_Q4_ELEMENT_NUM, the number of elements in the grid.
        //
    {
        int element_num;

        element_num = nelemx * nelemy;

        return element_num;
    }

    public static int sphere_grid_q4_node_num (int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        // 
        //    SPHERE_GRID_Q4_node_num = counts nodes in a Q4 sphere grid.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NELEMX, NELEMY, the number of elements along the
        //    X and Y directions.
        //
        //    Output, int SPHERE_GRID_Q4_node_num, the number of nodes in the grid.
        //
    {
        int node_num;

        node_num = nelemx * (nelemy - 1) + 2;

        return node_num;
    }

    public static double[] sphere_grid_q4_node_xyz(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_GRID_Q4_NODE_XYZ produces node coordinates for a Q4 sphere grid.
        //
        //  Discussion:
        //
        //    The number of nodes to be generated is
        //
        //      node_num = = NELEMX * ( NELEMY - 1 ) + 2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NELEMX, NELEMY, the number of elements along the
        //    X and Y directions.
        //
        //    Output, double SPHERE_GRID_Q4_NODE_XYZ[3*NODE_NUM], 
        //    the node coordinates.
        //
    {
        int i;
        int j;
        int node;
        int node_num;
        double[] node_xyz;
        double phi;
        double theta;

        node_num = sphere_grid_t6_node_num (nelemx, nelemy);
        node_xyz = new double[3 * node_num];
        node = 0;

        node_xyz[0 + node * 3] = 0.0;
        node_xyz[1 + node * 3] = 0.0;
        node_xyz[2 + node * 3] = -1.0;
        node += 1;

        for (j = nelemy; 2 <= j; j--)
        {
            phi = (j - 1) * Math.PI / nelemy;

            for (i = 1; i <= nelemx; i++)
            {
                theta = (i - 1) * 2.0 * Math.PI / nelemx;

                node_xyz[0 + node * 3] = Math.Cos(theta) * Math.Sin(phi);
                node_xyz[1 + node * 3] = Math.Sin(theta) * Math.Sin(phi);
                node_xyz[2 + node * 3] = Math.Cos(phi);
                node += 1;
            }
        }

        node_xyz[0 + node * 3] = 0.0;
        node_xyz[1 + node * 3] = 0.0;
        node_xyz[2 + node * 3] = 1.0;
        node += 1;

        return node_xyz;
    }

    public static int[] sphere_grid_q9_element(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_GRID_Q9_ELEMENT produces a Q9 sphere grid.
        //
        //  Discussion:
        //
        //    This would be the same as the grid in a plane, except that all the
        //    nodes along the bottom edge are identified (replaced by a single node
        //    that is the south pole) and similarly for the top edge, and the
        //    nodes on the extreme right edge are identified pairwise with those 
        //    on the extreme left edge.
        //
        //  Example:
        //
        //    Input:
        //
        //      NELEMX = 3, NELEMY = 4
        //
        //    Output:
        //
        //      ELEMENT_NODE =
        //         1,  1, 10,  8,  1,  4,  9,  2,  3;
        //         1,  1, 12, 10,  1,  6, 11,  4,  5;
        //         1,  1,  8, 12,  1,  2, 13,  6,  7;
        //         8, 10, 22, 20,  9, 16, 21, 14, 15;
        //        10, 12, 24, 22, 11, 18, 23, 16, 17;
        //        12,  8, 20, 24, 13, 14, 25, 18, 19;
        //        20, 22, 34, 32, 21, 28, 33, 26, 27;
        //        22, 24, 36, 34, 23, 30, 35, 28, 29;
        //        24, 20, 32, 36, 25, 26, 37, 30, 31;
        //        32, 34, 44, 44, 33, 40, 44, 38, 39;
        //        34, 36, 44, 44, 35, 42, 44, 40, 41;
        //        36, 32, 44, 44, 37, 38, 44, 42, 43;
        //
        //  Grid:
        //
        //   44-44-44-44-44-44-44
        //    |     |     |     |
        //   38 39 40 41 42 43 38
        //    |     |     |     |
        //   32-33-34-35-36-37-32
        //    |     |     |     |
        //   26 27 28 29 30 31 26
        //    |     |     |     |
        //   20-21-22-23-24-25-20
        //    |     |     |     |
        //   14 15 16 17 18 19 14
        //    |     |     |     |
        //    8--9-10-11-12-13--8
        //    |     |     |     |
        //    2  3  4  5  6  7  2
        //    |     |     |     |
        //    1--1--1--1--1--1--1
        //
        //  Reference Element Q9:
        //
        //    |
        //    1  4--7--3
        //    |  |     |
        //    |  |     |
        //    S  8  9  6
        //    |  |     |
        //    |  |     |
        //    0  1--5--2
        //    |
        //    +--0--R--1-->
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NELEMX, NELEMY, the number of elements along the
        //    X and Y directions.
        //
        //    Output, int SPHERE_GRID_Q9_ELEMENT[9*NELEMX*NELEMY], 
        //    the nodes that form each element.
        //
    {
        int base1;
        int base2;
        int element;
        int[] element_node;
        int element_num;
        int element_order = 9;
        int i;
        int j;

        element_num = sphere_grid_q9_element_num(nelemx, nelemy);
        element_node = new int[element_order * element_num];
        element = 0;

        for (j = 1; j <= nelemy; j++)
        {
            base1 = (j - 1) * 2 * 2 * nelemx + 2 - 2 * nelemx;

            for (i = 1; i <= nelemx; i++)
            {
                base2 = base1 + 2 * (i - 1);

                element_node[0 + element * element_order] = base2;
                element_node[4 + element * element_order] = base2 + 1;

                if (i < nelemx)
                {
                    element_node[1 + element * element_order] = base2 + 2;
                }
                else
                {
                    element_node[1 + element * element_order] = base1;
                }

                element_node[7 + element * element_order] =
                    element_node[0 + element * element_order] + 2 * nelemx;
                element_node[8 + element * element_order] =
                    element_node[4 + element * element_order] + 2 * nelemx;
                element_node[5 + element * element_order] =
                    element_node[1 + element * element_order] + 2 * nelemx;

                element_node[3 + element * element_order] =
                    element_node[7 + element * element_order] + 2 * nelemx;
                element_node[6 + element * element_order] =
                    element_node[8 + element * element_order] + 2 * nelemx;
                element_node[2 + element * element_order] =
                    element_node[5 + element * element_order] + 2 * nelemx;

                switch (j)
                {
                    case 1:
                        element_node[0 + element * element_order] = 1;
                        element_node[4 + element * element_order] = 1;
                        element_node[1 + element * element_order] = 1;
                        break;
                    default:
                    {
                        if (j == nelemy)
                        {
                            element_node[3 + element * element_order] = base1 + 4 * nelemx;
                            element_node[6 + element * element_order] = base1 + 4 * nelemx;
                            element_node[2 + element * element_order] = base1 + 4 * nelemx;
                        }

                        break;
                    }
                }

                element += 1;
            }
        }

        return element_node;
    }

    public static  int sphere_grid_q9_element_num(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_GRID_Q9_ELEMENT_NUM counts the elements in a Q9 sphere grid.
        //
        //  Example:
        //
        //    Input:
        //
        //      NELEMX = 3, NELEMY = 2
        //
        //    Output:
        //
        //      ELEMENT_NUM = NELEMX * NELEMY = 6
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NELEMX, NELEMY, the number of elements along the
        //    X and Y directions. 
        //
        //    Output, int SPHERE_GRID_Q9_ELEMENT_NUM, the number of elements in the grid.
        //
    {
        int element_num;

        element_num = nelemx * nelemy;

        return element_num;
    }

    public static int sphere_grid_q9_node_num (int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_GRID_Q9_node_num = counts nodes in a Q9 sphere grid.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NELEMX, NELEMY, the number of elements along the
        //    X and Y directions.
        //
        //    Output, int SPHERE_GRID_Q9_node_num, the number of nodes in the grid.
        //
    {
        int node_num;

        node_num = 4 * nelemx * nelemy - 2 * nelemx + 2;

        return node_num;
    }

    public static double[] sphere_grid_q9_node_xyz(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_GRID_Q9_NODE_XYZ produces node coordinates for a Q9 sphere grid.
        //
        //  Discussion:
        //
        //    The number of nodes to be generated is
        //
        //      node_num = = 4 * NELEMX * NELEMY - 2 * NELEMX + 2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NELEMX, NELEMY, the number of elements along the
        //    X and Y directions.  
        //
        //    Output, double SPHERE_GRID_Q9_NODE_XYZ[3*NODE_NUM], 
        //    the node coordinates.
        //
    {
        int i;
        int j;
        int node;
        int node_num;
        double[] node_xyz;
        double phi;
        double theta;

        node_num = sphere_grid_q9_node_num (nelemx, nelemy);
        node_xyz = new double[3 * node_num];
        node = 0;

        node_xyz[0 + node * 3] = 0.0;
        node_xyz[1 + node * 3] = 0.0;
        node_xyz[2 + node * 3] = -1.0;
        node += 1;

        for (j = 2 * nelemy; 2 <= j; j--)
        {
            phi = (j - 1) * Math.PI / (2 * nelemy);

            for (i = 1; i <= 2 * nelemx; i++)
            {
                theta = (i - 1) * 2.0 * Math.PI / (2 * nelemx);

                node_xyz[0 + node * 3] = Math.Cos(theta) * Math.Sin(phi);
                node_xyz[1 + node * 3] = Math.Sin(theta) * Math.Sin(phi);
                node_xyz[2 + node * 3] = Math.Cos(phi);
                node += 1;
            }
        }

        node_xyz[0 + node * 3] = 0.0;
        node_xyz[1 + node * 3] = 0.0;
        node_xyz[2 + node * 3] = 1.0;
        node += 1;

        return node_xyz;
    }

    public static int[] sphere_grid_q16_element(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_GRID_Q16_ELEMENT produces a Q16 sphere grid.
        //
        //  Discussion:
        //
        //    This would be the same as the grid in a plane, except that all the
        //    nodes along the bottom edge are identified (replaced by a single node
        //    that is the south pole) and similarly for the top edge, and the
        //    nodes on the extreme right edge are identified pairwise with those 
        //    on the extreme left edge.
        //
        //  Example:
        //
        //    Input:
        //
        //      NELEMX = 2, NELEMY = 2
        //
        //    Output:
        //
        //      ELEMENT_NODE =
        //         1,  1,  1,  1,  2,  3,  4,  5,  8,  9, 10, 11, 14, 15, 16, 17;
        //         1,  1,  1,  1,  5,  6,  7,  2, 11, 12, 13,  8, 17, 18, 19, 14;
        //        14, 15, 16, 17, 20, 21, 22, 23, 26, 27, 28, 29, 32, 32, 32, 32;
        //        17, 18, 19, 14, 23, 24, 25, 20, 29, 30, 31, 26, 32, 32, 32, 32.
        //
        //  Grid:
        //
        //   32-32-32-32-32-32-32
        //    |        |        |
        //    |        |        |
        //   26 27 28 29 30 31 26
        //    |        |        |
        //    |        |        |
        //   20 21 22 23 24 25 20
        //    |        |        |
        //    | E3     | E4     |
        //   14-15-16-17-18-19-14
        //    |        |        |
        //    |        |        |
        //    8  9 10 11 12 13  8
        //    |        |        |
        //    |        |        |
        //    2  3  4  5  6  7  2
        //    |        |        |
        //    | E1     | E2     |
        //    1--1--1--1--1--1--1
        //
        //  Reference Element Q16:
        //
        //    |
        //    1 13--14--15--16
        //    |  |   :   :   |
        //    |  |   :   :   |
        //    |  9..10..11..12
        //    S  |   :   :   |
        //    |  |   :   :   |
        //    |  5...6...7...8
        //    |  |   :   :   |
        //    |  |   :   :   |
        //    0  1---2---3---4
        //    |
        //    +--0-----R-----1-->
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NELEMX, NELEMY, the number of elements along the
        //    X and Y directions.  The number of elements generated will be
        //    NELEMX * NELEMY.
        //
        //    Output, int SPHERE_GRID_Q16_ELEMENT[16*NELEMX*NELEMY], 
        //    the nodes that form each element.
        //
    {
        int base1;
        int base2;
        int element;
        int[] element_node;
        int element_order = 16;
        int i;
        int j;

        element_node = new int[element_order * nelemx * nelemy];

        element = 0;

        for (j = 1; j <= nelemy; j++)
        {
            base1 = (j - 1) * 3 * 3 * nelemx + 2 - 3 * nelemx;

            for (i = 1; i <= nelemx; i++)
            {
                base2 = base1 + 3 * (i - 1);

                element_node[0 + element * element_order] = base2;
                element_node[1 + element * element_order] = base2 + 1;
                element_node[2 + element * element_order] = base2 + 2;

                if (i < nelemx)
                {
                    element_node[3 + element * element_order] = base2 + 3;
                }
                else
                {
                    element_node[3 + element * element_order] = base1;
                }

                element_node[4 + element * element_order] =
                    element_node[0 + element * element_order] + 3 * nelemx;
                element_node[5 + element * element_order] =
                    element_node[1 + element * element_order] + 3 * nelemx;
                element_node[6 + element * element_order] =
                    element_node[2 + element * element_order] + 3 * nelemx;
                element_node[7 + element * element_order] =
                    element_node[3 + element * element_order] + 3 * nelemx;

                element_node[8 + element * element_order] =
                    element_node[4 + element * element_order] + 3 * nelemx;
                element_node[9 + element * element_order] =
                    element_node[5 + element * element_order] + 3 * nelemx;
                element_node[10 + element * element_order] =
                    element_node[6 + element * element_order] + 3 * nelemx;
                element_node[11 + element * element_order] =
                    element_node[7 + element * element_order] + 3 * nelemx;

                element_node[12 + element * element_order] =
                    element_node[8 + element * element_order] + 3 * nelemx;
                element_node[13 + element * element_order] =
                    element_node[9 + element * element_order] + 3 * nelemx;
                element_node[14 + element * element_order] =
                    element_node[10 + element * element_order] + 3 * nelemx;
                element_node[15 + element * element_order] =
                    element_node[11 + element * element_order] + 3 * nelemx;

                switch (j)
                {
                    case 1:
                        element_node[0 + element * element_order] = 1;
                        element_node[1 + element * element_order] = 1;
                        element_node[2 + element * element_order] = 1;
                        element_node[3 + element * element_order] = 1;
                        break;
                    default:
                    {
                        if (j == nelemy)
                        {
                            element_node[12 + element * element_order] = base1 + 9 * nelemx;
                            element_node[13 + element * element_order] = base1 + 9 * nelemx;
                            element_node[14 + element * element_order] = base1 + 9 * nelemx;
                            element_node[15 + element * element_order] = base1 + 9 * nelemx;
                        }

                        break;
                    }
                }

                element += 1;
            }
        }

        return element_node;
    }

    public static int sphere_grid_q16_element_num(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_GRID_Q16_ELEMENT_NUM counts the elements in a Q16 sphere grid.
        //
        //  Example:
        //
        //    Input:
        //
        //      NELEMX = 3, NELEMY = 2
        //
        //    Output:
        //
        //      ELEMENT_NUM = NELEMX * NELEMY = 6
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NELEMX, NELEMY, the number of elements along the
        //    X and Y directions. 
        //
        //    Output, int SPHERE_GRID_Q16_ELEMENT_NUM, the number of elements in the grid.
        //
    {
        int element_num;

        element_num = nelemx * nelemy;

        return element_num;
    }

    public static int sphere_grid_q16_node_num (int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_GRID_Q16_node_num = counts nodes in a Q16 sphere grid.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NELEMX, NELEMY, the number of elements along the
        //    X and Y directions.
        //
        //    Output, int SPHERE_GRID_Q16_node_num, the number of nodes in the grid.
        //
    {
        int node_num;

        node_num = 9 * nelemx * nelemy - 3 * nelemx + 2;

        return node_num;
    }

    public static double[] sphere_grid_q16_node_xyz(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_GRID_Q16_NODE_XYZ produces node coordinates for a Q16 sphere grid.
        //
        //  Discussion:
        //
        //    The number of nodes to be generated is
        //
        //      node_num = = 9 * NELEMX * NELEMY - 3 * NELEMX + 2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NELEMX, NELEMY, the number of elements along the
        //    X and Y directions.
        //
        //    Output, double SPHERE_GRID_Q16_NODE_XYZ[3*NODE_NUM], the node coordinates.
        //
    {
        int i;
        int j;
        int node;
        int node_num;
        double[] node_xyz;
        double phi;
        double theta;

        node_num = sphere_grid_q16_node_num (nelemx, nelemy);
        node_xyz = new double[3 * node_num];
        node = 0;

        for (j = 3 * nelemy + 1; 1 <= j; j--)
        {
            phi = (j - 1) * Math.PI / (3 * nelemy);

            switch (j)
            {
                case 1:
                    node_xyz[0 + node * 3] = 0.0;
                    node_xyz[1 + node * 3] = 0.0;
                    node_xyz[2 + node * 3] = 1.0;
                    node += 1;
                    break;
                default:
                {
                    if (j < 3 * nelemy + 1)
                    {
                        for (i = 1; i <= 3 * nelemx; i++)
                        {
                            theta = (i - 1) * 2.0 * Math.PI / (3 * nelemx);

                            node_xyz[0 + node * 3] = Math.Cos(theta) * Math.Sin(phi);
                            node_xyz[1 + node * 3] = Math.Sin(theta) * Math.Sin(phi);
                            node_xyz[2 + node * 3] = Math.Cos(phi);
                            node += 1;
                        }
                    }
                    else if (j == 3 * nelemy + 1)
                    {
                        node_xyz[0 + node * 3] = 0.0;
                        node_xyz[1 + node * 3] = 0.0;
                        node_xyz[2 + node * 3] = -1.0;
                        node += 1;
                    }

                    break;
                }
            }
        }

        return node_xyz;
    }

    public static int[] sphere_grid_t3_element(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_GRID_T3_ELEMENT produces a T3 sphere grid.
        //
        //  Discussion:
        //
        //    This would be the same as the grid in a plane, except that all the
        //    nodes along the bottom edge are identified (replaced by a single node
        //    that is the south pole) and similarly for the top edge, and the
        //    nodes on the extreme right edge are identified pairwise with those 
        //    on the extreme left edge.
        //
        //  Example:
        //
        //    Input:
        //
        //      NELEMX = 3, NELEMY = 4
        //
        //    Output:
        //
        //      ELEMENT_NODE =
        //         1,  3,  2;
        //         1,  4,  3;
        //         1,  2,  4;
        //         2,  3,  5
        //         6,  5,  3
        //         3,  4,  6
        //         7,  6,  4;
        //         4,  2,  7;
        //         5,  7,  2;
        //         5,  6,  8;
        //         9,  8,  6;
        //         6,  7,  9;
        //        10,  9,  7;
        //         7,  5, 10;
        //         8, 10,  5;
        //         8,  9, 11;
        //         9, 10, 11;
        //        10,  8, 11;
        //
        //  Grid:
        //
        //   11    11    11    11
        //    | .   | .   | .   |
        //    |  .  |  .  |  .  |
        //    |E16. |E17 .|E18. |
        //    8-----9----10-----8
        //    | .E11| .E13| .E15|
        //    |  .  |  .  |  .  |
        //    |E10. |E12. |E14. |
        //    5-----6-----7-----5
        //    | .E5 | .E7 | .E9 |
        //    |  .  |  .  |  .  |
        //    |E4 . |E6 . |E8 . |
        //    2-----3-----4-----2
        //      .E1 | .E2 | .E3 |
        //       .  |  .  |  .  |
        //        . |   . |   . |
        //          1     1     1
        //
        //  Reference Element T3:
        //
        //    |
        //    1  3
        //    |  ..
        //    |  . .
        //    S  .  .
        //    |  .   .
        //    |  .    .
        //    0  1-----2
        //    |
        //    +--0--R--1-->
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NELEMX, NELEMY, the number of elements along the
        //    X and Y directions.
        //
        //    Output, int SPHERE_GRID_T3_ELEMENT[3*(2*NELEMX*(NELEMY-1)], 
        //    the nodes that form each element.
        //
    {
        int base1;
        int base2;
        int element;
        int[] element_node;
        int element_num;
        int element_order = 3;
        int i;
        int j;
        int ne;
        int nw;
        int se;
        int sw;

        element_num = sphere_grid_t3_element_num(nelemx, nelemy);
        element_node = new int[element_order * element_num];

        element = 0;

        for (j = 1; j <= nelemy; j++)
        {
            base1 = (j - 1) * nelemx + 2 - nelemx;

            for (i = 1; i <= nelemx; i++)
            {
                base2 = base1 + i - 1;

                sw = base2;
                if (i < nelemx)
                {
                    se = base2 + 1;
                }
                else
                {
                    se = base1;
                }

                nw = sw + nelemx;
                ne = se + nelemx;

                switch (j)
                {
                    case 1:
                        sw = 1;
                        se = 1;
                        break;
                    default:
                    {
                        if (j == nelemx)
                        {
                            nw = base1 + nelemx;
                            ne = base1 + nelemx;
                        }

                        break;
                    }
                }

                switch (j)
                {
                    case > 1:
                        element_node[0 + element * element_order] = sw;
                        element_node[1 + element * element_order] = se;
                        element_node[2 + element * element_order] = nw;
                        element += 1;
                        break;
                }

                if (j < nelemy)
                {
                    element_node[0 + element * element_order] = ne;
                    element_node[1 + element * element_order] = nw;
                    element_node[2 + element * element_order] = se;
                    element += 1;
                }
            }
        }

        return element_node;
    }

    public static int sphere_grid_t3_element_num(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_GRID_T3_ELEMENT_NUM counts the elements in a T3 sphere grid.
        //
        //  Example:
        //
        //    Input:
        //
        //      NELEMX = 6, NELEMY = 4
        //
        //    Output:
        //
        //      ELEMENT_NUM = 2 * NELEMX * ( NELEMY - 1 ) = 36
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NELEMX, NELEMY, the number of elements along the
        //    X and Y directions. 
        //
        //    Output, int SPHERE_GRID_T3_ELEMENT_NUM, the number of elements in the grid.
        //
    {
        int element_num;

        element_num = 2 * nelemx * (nelemy - 1);

        return element_num;
    }

    public static int sphere_grid_t3_node_num (int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_GRID_T3_node_num = counts nodes in a T3 sphere grid.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NELEMX, NELEMY, the number of elements along the
        //    X and Y directions.
        //
        //    Output, int SPHERE_GRID_T3_node_num, the number of nodes in the grid.
        //
    {
        int node_num;

        node_num = nelemx * (nelemy - 1) + 2;

        return node_num;
    }

    public static double[] sphere_grid_t3_node_xyz(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_GRID_T3_NODE_XYZ produces node coordinates for a T3 sphere grid.
        //
        //  Discussion:
        //
        //    The number of nodes to be generated is
        //
        //      node_num = = NELEMX * ( NELEMY - 1 ) + 2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NELEMX, NELEMY, the number of elements along the
        //    X and Y directions. 
        //
        //    Output, double SPHERE_GRID_T3_NODE_XYZ[3*NODE_NUM], 
        //    the node coordinates.
        //
    {
        int i;
        int j;
        int node;
        int node_num;
        double[] node_xyz;
        double phi;
        double theta;

        node_num = sphere_grid_t3_node_num (nelemx, nelemy);
        node_xyz = new double[3 * node_num];

        node = 0;

        node_xyz[0 + node * 3] = 0.0;
        node_xyz[1 + node * 3] = 0.0;
        node_xyz[2 + node * 3] = -1.0;
        node += 1;

        for (j = nelemy; 2 <= j; j--)
        {
            phi = (j - 1) * Math.PI / nelemy;

            for (i = 1; i <= nelemx; i++)
            {
                theta = (i - 1) * 2.0 * Math.PI / nelemx;

                node_xyz[0 + node * 3] = Math.Cos(theta) * Math.Sin(phi);
                node_xyz[1 + node * 3] = Math.Sin(theta) * Math.Sin(phi);
                node_xyz[2 + node * 3] = Math.Cos(phi);
                node += 1;
            }
        }

        node_xyz[0 + node * 3] = 0.0;
        node_xyz[1 + node * 3] = 0.0;
        node_xyz[2 + node * 3] = 1.0;
        node += 1;

        return node_xyz;
    }

    public static int[] sphere_grid_t6_element(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_GRID_T6_ELEMENT produces a T6 sphere grid.
        //
        //  Discussion:
        //
        //    This would be the same as the grid in a plane, except that all the
        //    nodes along the bottom edge are identified (replaced by a single node
        //    that is the south pole) and similarly for the top edge, and the
        //    nodes on the extreme right edge are identified pairwise with those 
        //    on the extreme left edge.
        //
        //  Example:
        //
        //    Input:
        //
        //      NELEMX = 3, NELEMY = 4
        //
        //    Output:
        //
        //      ELEMENT_NODE =
        //        10,  8,  1,  9,  3,  4;
        //        12, 10,  1, 11,  5,  6;
        //         8, 12,  1, 13,  7,  2;
        //         8, 10, 20,  9, 15, 14;
        //        22, 20, 10, 21, 15, 16;
        //        10, 12, 22, 11, 17, 16;
        //        24, 22, 12, 23, 17, 18;
        //        12,  8, 24, 13, 19, 18;
        //        20, 24,  8, 25, 19, 14;
        //        20, 22, 32, 21, 27, 26;
        //        34, 32, 22, 33, 27, 28;
        //        22, 24, 34, 23, 29, 28;
        //        36, 34, 24, 35, 29, 30;
        //        24, 20, 36, 25, 31, 30;
        //        32, 36, 20, 37, 31, 26;
        //        32, 34, 44, 33, 39, 38;
        //        34, 36, 44, 35, 41, 40;
        //        36, 32, 44, 37, 43, 42;
        //
        //  Grid:
        //
        //   44    44    44
        //    |.    |.    |.
        //   38 39 40 41 42 43 
        //    |    .|    .|    .
        //   32-33-34-35-36-37-32
        //    |.    |.    |.    |
        //   26 27 28 29 30 31 26
        //    |    .|    .|    .|
        //   20-21-22-23-24-25-20
        //    |.    |.    |.    |
        //   14 15 16 17 18 19 14
        //    |    .|    .|    .|
        //    8--9-10-11-12-13--8
        //     .    |.    |.    |
        //       3  4  5  6  7  2
        //         .|    .|    .|
        //          1     1     1
        //
        //  Reference Element T6:
        //
        //    |
        //    1  3
        //    |  ..
        //    |  . .
        //    S  6  5
        //    |  .   .
        //    |  .    .
        //    0  1--4--2
        //    |
        //    +--0--R--1-->
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NELEMX, NELEMY, the number of elements along the
        //    X and Y directions.
        //
        //    Output, int SPHERE_GRID_T6_ELEMENT[6*2*NELEMX*(NELEMY-1)], 
        //    the nodes that form each element.
        //
    {
        int base1;
        int base2;
        int c;
        int e;
        int element;
        int[] element_node;
        int element_num;
        int element_order = 6;
        int i;
        int j;
        int n;
        int ne;
        int nw;
        int s;
        int se;
        int sw;
        int w;

        element_num = sphere_grid_t6_element_num(nelemx, nelemy);
        element_node = new int[element_order * element_num];
        element = 0;

        for (j = 1; j <= nelemy; j++)
        {
            base1 = (j - 1) * 2 * 2 * nelemx + 2 - 2 * nelemx;

            for (i = 1; i <= nelemx; i++)
            {
                base2 = base1 + 2 * (i - 1);

                sw = base2;
                s = base2 + 1;
                if (i < nelemx)
                {
                    se = base2 + 2;
                }
                else
                {
                    se = base1;
                }

                w = sw + 2 * nelemx;
                c = s + 2 * nelemx;
                e = se + 2 * nelemx;

                nw = w + 2 * nelemx;
                n = c + 2 * nelemx;
                ne = e + 2 * nelemx;

                switch (j)
                {
                    case 1:
                        sw = 1;
                        s = 1;
                        se = 1;
                        break;
                    default:
                    {
                        if (j == nelemy)
                        {
                            nw = base1 + 4 * nelemx;
                            n = base1 + 4 * nelemx;
                            ne = base1 + 4 * nelemx;
                        }

                        break;
                    }
                }

                switch (j)
                {
                    case > 1:
                        element_node[0 + element * element_order] = sw;
                        element_node[1 + element * element_order] = se;
                        element_node[2 + element * element_order] = nw;
                        element_node[3 + element * element_order] = s;
                        element_node[4 + element * element_order] = c;
                        element_node[5 + element * element_order] = w;
                        element += 1;
                        break;
                }

                if (j < nelemy)
                {
                    element_node[0 + element * element_order] = ne;
                    element_node[1 + element * element_order] = nw;
                    element_node[2 + element * element_order] = se;
                    element_node[3 + element * element_order] = n;
                    element_node[4 + element * element_order] = c;
                    element_node[5 + element * element_order] = e;
                    element += 1;
                }
            }
        }

        return element_node;
    }

    public static int sphere_grid_t6_element_num(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_GRID_T6_ELEMENT_NUM counts the elements in a T6 sphere grid.
        //
        //  Example:
        //
        //    Input:
        //
        //      NELEMX = 6, NELEMY = 4
        //
        //    Output:
        //
        //      ELEMENT_NUM = 2 * NELEMX * ( NELEMY - 1 ) = 36
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NELEMX, NELEMY, the number of elements along the
        //    X and Y directions. 
        //
        //    Output, int SPHERE_GRID_T6_ELEMENT_NUM, the number of elements in the grid.
        //
    {
        int element_num;

        element_num = 2 * nelemx * (nelemy - 1);

        return element_num;
    }

    public static int sphere_grid_t6_node_num (int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_GRID_T6_node_num = counts nodes in a T6 sphere grid.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NELEMX, NELEMY, the number of elements along the
        //    X and Y directions.
        //
        //    Output, int SPHERE_GRID_T6_node_num, the number of nodes in the grid.
        //
    {
        int node_num;

        node_num = 4 * nelemx * nelemy - 2 * nelemx + 2;

        return node_num;
    }

    public static double[] sphere_grid_t6_node_xyz(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_GRID_T6_NODE_XYZ produces node coordinates for a T6 sphere grid.
        //
        //  Discussion:
        //
        //    The number of nodes to be generated is
        //
        //      node_num = = 4 * NELEMX * NELEMY - 2 * NELEMX + 2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NELEMX, NELEMY, the number of elements along the
        //    X and Y directions.  
        //
        //    Output, double SPHERE_GRID_T6_NODE_XYZ[3*NODE_NUM], 
        //    the node coordinates.
        //
    {
        int i;
        int j;
        int node;
        int node_num;
        double[] node_xyz;
        double phi;
        double theta;

        node_num = sphere_grid_t6_node_num (nelemx, nelemy);
        node_xyz = new double[3 * node_num];
        node = 0;

        node_xyz[0 + node * 3] = 0.0;
        node_xyz[1 + node * 3] = 0.0;
        node_xyz[2 + node * 3] = -1.0;
        node += 1;

        for (j = 2 * nelemy; 2 <= j; j--)
        {
            phi = (j - 1) * Math.PI / (2 * nelemy);

            for (i = 1; i <= 2 * nelemx; i++)
            {
                theta = (i - 1) * 2.0 * Math.PI / (2 * nelemx);

                node_xyz[0 + node * 3] = Math.Cos(theta) * Math.Sin(phi);
                node_xyz[1 + node * 3] = Math.Sin(theta) * Math.Sin(phi);
                node_xyz[2 + node * 3] = Math.Cos(phi);
                node += 1;
            }
        }

        node_xyz[0 + node * 3] = 0.0;
        node_xyz[1 + node * 3] = 0.0;
        node_xyz[2 + node * 3] = 1.0;
        node += 1;

        return node_xyz;
    }

    public static int[] sphere_grid_q4(int lat_num, int long_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_GRID_Q4: rectangular grid on a sphere.
        //
        //  Discussion:
        //
        //    The point numbering system is the same used in SPHERE_GRIDPOINTS,
        //    and that routine may be used to compute the coordinates of the points.
        //
        //    A sphere in 3D satisfies the equation:
        //
        //      sum ( ( P(1:DIM_NUM) - PC(1:DIM_NUM) )^2 ) = R * R
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int LAT_NUM, the number of "rows" of rectangles to
        //    be created.  LAT_NUM must be at least 2. 
        //
        //    Input, int LONG_NUM, the number of "columns" of 
        //    rectangles to be created.
        //
        //    Output, int SPHERE_GRID_Q4[4*LAT_NUM*LONG_NUM], 
        //    the indices of the nodes that make up each rectangle.
        //
    {
        int i;
        int j;
        int n;
        int n_max;
        int n_min;
        int ne;
        int nw;
        int s;
        int s_max;
        int s_min;
        int se;
        int sw;
        int[] rectangle_node;
        int rectangle_num;

        rectangle_node = new int[4 * lat_num * long_num];
        rectangle_num = 0;
        //
        //  The first row.
        //
        n = 0;

        sw = 1;
        se = sw + 1;

        s_min = 1;
        s_max = long_num;

        for (j = 1; j <= long_num; j++)
        {
            rectangle_node[0 + rectangle_num * 4] = sw;
            rectangle_node[1 + rectangle_num * 4] = se;
            rectangle_node[2 + rectangle_num * 4] = n;
            rectangle_node[3 + rectangle_num * 4] = n;
            rectangle_num += 1;

            sw = se;

            if (se == s_max)
            {
                se = s_min;
            }
            else
            {
                se += 1;
            }
        }

        //
        //  The intermediate rows.
        //
        for (i = 2; i < lat_num; i++)
        {
            n_max = s_max;
            n_min = s_min;

            s_max += long_num;
            s_min += long_num;

            nw = n_min;
            ne = nw + 1;
            sw = s_min;
            se = sw + 1;

            for (j = 1; j <= long_num; j++)
            {

                rectangle_node[0 + rectangle_num * 4] = sw;
                rectangle_node[1 + rectangle_num * 4] = se;
                rectangle_node[2 + rectangle_num * 4] = ne;
                rectangle_node[3 + rectangle_num * 4] = nw;
                rectangle_num += 1;

                sw = se;
                nw = ne;
                if (se == s_max)
                {
                    se = s_min;
                }
                else
                {
                    se += 1;
                }

                if (ne == n_max)
                {
                    ne = n_min;
                }
                else
                {
                    ne += 1;
                }
            }
        }

        //
        //  The last row.
        //
        n_max = s_max;
        n_min = s_min;

        s = n_max + 1;

        nw = n_min;
        ne = nw + 1;

        for (j = 1; j <= long_num; j++)
        {
            rectangle_node[0 + rectangle_num * 4] = ne;
            rectangle_node[1 + rectangle_num * 4] = nw;
            rectangle_node[2 + rectangle_num * 4] = s;
            rectangle_node[3 + rectangle_num * 4] = s;
            rectangle_num += 1;

            nw = ne;
            if (ne == n_max)
            {
                ne = n_min;
            }
            else
            {
                ne += 1;
            }
        }

        return rectangle_node;
    }

    public static int[] sphere_grid_t3(int lat_num, int long_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_GRID_T3 produces a triangle grid on a sphere.
        //
        //  Discussion:
        //
        //    I want to convert this to 0-based indexing.  However, my first attempt
        //    to modify the internal coding to compute the 0-based indexes directly
        //    was a failure.  So I reverted to computing 1-based node indices,
        //    and then subtracting 1.  You would think this would be easy to fix.
        //
        //    The point numbering system is the same used in SPHERE_GRIDPOINTS,
        //    and that routine may be used to compute the coordinates of the points.
        //
        //    A sphere in 3D satisfies the equation:
        //
        //      sum ( ( P(1:DIM_NUM) - pc(1:DIM_NUM) )^2 ) = R * R
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int LAT_NUM, LONG_NUM, the number of latitude 
        //    and longitude lines to draw.  The latitudes do not include the North 
        //    and South poles, which will be included automatically, so LAT_NUM = 5, 
        //    for instance, will result in points along 7 lines of latitude.
        //
        //    Output, int SPHERE_GRID_T3[3*2*(LAT_NUM+1)*LONG_NUM], the
        //    triangle vertices.
        //
    {
        int i;
        int j;
        int n;
        int n_max;
        int n_min;
        int ne;
        int nw;
        int s;
        int s_max;
        int s_min;
        int se;
        int sw;
        int[] triangle_node;
        int triangle_num;

        triangle_node = new int[3 * 2 * (lat_num + 1) * long_num];
        triangle_num = 0;
        //
        //  The first row.
        //
        //n = 0;
        //sw = 1;
        //se = sw + 1;
        //s_min = 1;
        //s_max = long_num;

        n = 1;
        sw = 2;
        se = sw + 1;
        s_min = 2;
        s_max = long_num + 1;

        for (j = 0; j < long_num; j++)
        {
            triangle_node[0 + triangle_num * 3] = sw - 1;
            triangle_node[1 + triangle_num * 3] = se - 1;
            triangle_node[2 + triangle_num * 3] = n - 1;
            triangle_num += 1;

            sw = se;
            if (se == s_max)
            {
                se = s_min;
            }
            else
            {
                se += 1;
            }
        }

        //
        //  The intermediate rows.
        //
        for (i = 1; i <= lat_num; i++)
        {
            n_max = s_max;
            n_min = s_min;

            s_max += long_num;
            s_min += long_num;

            nw = n_min;
            ne = nw + 1;
            sw = s_min;
            se = sw + 1;

            for (j = 0; j < long_num; j++)
            {
                triangle_node[0 + triangle_num * 3] = sw - 1;
                triangle_node[1 + triangle_num * 3] = se - 1;
                triangle_node[2 + triangle_num * 3] = nw - 1;
                triangle_num += 1;

                triangle_node[0 + triangle_num * 3] = ne - 1;
                triangle_node[1 + triangle_num * 3] = nw - 1;
                triangle_node[2 + triangle_num * 3] = se - 1;
                triangle_num += 1;

                sw = se;
                nw = ne;
                if (se == s_max)
                {
                    se = s_min;
                }
                else
                {
                    se += 1;
                }

                if (ne == n_max)
                {
                    ne = n_min;
                }
                else
                {
                    ne += 1;
                }
            }
        }

        //
        //  The last row.
        //
        n_max = s_max;
        n_min = s_min;

        s = n_max + 1;

        nw = n_min;
        ne = nw + 1;

        for (j = 0; j < long_num; j++)
        {
            triangle_node[0 + triangle_num * 3] = ne - 1;
            triangle_node[1 + triangle_num * 3] = nw - 1;
            triangle_node[2 + triangle_num * 3] = s - 1;
            triangle_num += 1;

            nw = ne;

            if (ne == n_max)
            {
                ne = n_min;
            }
            else
            {
                ne += 1;
            }
        }

        return triangle_node;
    }

}