using System;
using System.Globalization;
using Burkardt.OrderNS;

namespace Burkardt.FEM;

public static class Grid
{
    public static int[] grid_element(string code, int element_order, int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRID_ELEMENT returns the element grid associated with any available element.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string CODE, identifies the element desired.
        //    Legal values include "Q4", "Q8", "Q9", "Q12", "Q16", "QL", "T3", 
        //    "T4", "T6" and "T10".
        //
        //    Input, int ELEMENT_ORDER, the number of nodes per element.
        //
        //    Input, int NELEMX, NELEMY, the number of quadrilaterals along the
        //    X and Y directions.  The number of elements generated will be
        //    NELEMX * NELEMY for quadrilaterals, or 2 * NELEMX * NELEMY for
        //    triangles.
        //
        //    Output, int GRID_ELEMENT[ELEMENT_ORDER*ELEMENT_NUM], the nodes 
        //    that form each element.  
        //
    {
        int[] element_node;

        switch (code)
        {
            case "Q4":
                element_node = grid_q4_element(nelemx, nelemy);
                break;
            case "Q8":
                element_node = grid_q8_element(nelemx, nelemy);
                break;
            case "Q9":
                element_node = grid_q9_element(nelemx, nelemy);
                break;
            case "Q12":
                element_node = grid_q12_element(nelemx, nelemy);
                break;
            case "Q16":
                element_node = grid_q16_element(nelemx, nelemy);
                break;
            case "QL":
                element_node = grid_ql_element(nelemx, nelemy);
                break;
            case "T3":
                element_node = grid_t3_element(nelemx, nelemy);
                break;
            case "T4":
                element_node = grid_t4_element(nelemx, nelemy);
                break;
            case "T6":
                element_node = grid_t6_element(nelemx, nelemy);
                break;
            case "T10":
                element_node = grid_t10_element(nelemx, nelemy);
                break;
            default:
                element_node = null;
                Console.WriteLine("");
                Console.WriteLine("GRID_ELEMENT - Fatal error!");
                Console.WriteLine("  Illegal value of CODE = \"" + code + "\".");
                break;
        }

        return element_node;
    }

    public static int grid_element_num(string code, int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRID_ELEMENT_NUM returns the number of elements in a grid.
        //
        //  Discussion:
        //
        //    The number of elements generated will be NELEMX * NELEMY for
        //    quadrilaterals, or 2 * NELEMX * NELEMY for triangles.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 March 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string CODE, identifies the element desired.
        //    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 'T3', 
        //    'T4', 'T6' and 'T10'.
        //
        //    Input, int NELEMX, NELEMY, the number of quadrilaterals along the
        //    X and Y directions.  
        //
        //    Output, int GRID_ELEMENT_NUM, the number of elements in the grid.
        //
    {
        int element_num;

        switch (code)
        {
            case "Q4":
                element_num = grid_q4_element_num(nelemx, nelemy);
                break;
            case "Q8":
                element_num = grid_q8_element_num(nelemx, nelemy);
                break;
            case "Q9":
                element_num = grid_q9_element_num(nelemx, nelemy);
                break;
            case "Q12":
                element_num = grid_q12_element_num(nelemx, nelemy);
                break;
            case "Q16":
                element_num = grid_q16_element_num(nelemx, nelemy);
                break;
            case "QL":
                element_num = grid_ql_element_num(nelemx, nelemy);
                break;
            case "T3":
                element_num = grid_t3_element_num(nelemx, nelemy);
                break;
            case "T4":
                element_num = grid_t4_element_num(nelemx, nelemy);
                break;
            case "T6":
                element_num = grid_t6_element_num(nelemx, nelemy);
                break;
            case "T10":
                element_num = grid_t10_element_num(nelemx, nelemy);
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("GRID_ELEMENT_NUM - Fatal error!");
                Console.WriteLine("  Illegal value of CODE = \"" + code + "\".");
                element_num = -1;
                break;
        }

        return element_num;
    }

    public static int grid_node_num (string code, int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRID_node_num = returns the number of nodes in a grid.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 March 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string CODE, identifies the element desired.
        //    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 'T3', 
        //    'T4', 'T6' and 'T10'.
        //
        //    Input, int NELEMX, NELEMY, the number of quadrilaterals along the
        //    X and Y directions.  
        //
        //    Output, int GRID_node_num, the number of elements in the grid.
        //
    {
        int node_num;

        switch (code)
        {
            case "Q4":
                node_num = grid_q4_node_num(nelemx, nelemy);
                break;
            case "Q8":
                node_num = grid_q8_node_num(nelemx, nelemy);
                break;
            case "Q9":
                node_num = grid_q9_node_num(nelemx, nelemy);
                break;
            case "Q12":
                node_num = grid_q12_node_num(nelemx, nelemy);
                break;
            case "Q16":
                node_num = grid_q16_node_num(nelemx, nelemy);
                break;
            case "QL":
                node_num = grid_ql_node_num(nelemx, nelemy);
                break;
            case "T3":
                node_num = grid_t3_node_num(nelemx, nelemy);
                break;
            case "T4":
                node_num = grid_t4_node_num(nelemx, nelemy);
                break;
            case "T6":
                node_num = grid_t6_node_num(nelemx, nelemy);
                break;
            case "T10":
                node_num = grid_t10_node_num(nelemx, nelemy);
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("GRID_node_num = - Fatal error!");
                Console.WriteLine("  Illegal value of CODE = \"" + code + "\".");
                node_num = -1;
                return 1;
        }

        return node_num;
    }

    public static double[] grid_nodes_01(int x_num, int y_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRID_NODES_01 returns an equally spaced rectangular grid of nodes in the unit square.
        //
        //  Example:
        //
        //    X_NUM = 5
        //    Y_NUM = 3
        //
        //    NODE_XY =
        //    ( 0, 0.25, 0.5, 0.75, 1, 0,   0.25, 0.5, 0.75, 1,   0, 0.25, 0.5, 0.75, 1;
        //      0, 0,    0,   0,    0, 0.5, 0.5,  0.5, 0.5,  0.5, 1, 1.0,  1.0, 1.0,  1 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 May 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X_NUM, Y_NUM, the number of nodes in the X and Y directions.
        //
        //    Output, double GRID_NODES_01[2*X_NUM*Y_NUM], the coordinates of
        //    the nodes.
        //
    {
        int i;
        int j;
        int k;
        double value;

        int node_num = x_num* y_num;

        double[] node_xy = new double[2 * node_num];

        switch (x_num)
        {
            case 1:
            {
                for (k = 0; k < node_num; k++)
                {
                    node_xy[0 + k * 2] = 0.5;
                }

                break;
            }
            default:
            {
                for (i = 0; i < x_num; i++)
                {
                    value = i / (double) (x_num - 1);
                    for (j = i; j < node_num; j += x_num)
                    {
                        node_xy[0 + j * 2] = value;
                    }
                }

                break;
            }
        }

        switch (y_num)
        {
            case 1:
            {
                for (k = 0; k < node_num; k++)
                {
                    node_xy[1 + k * 2] = 0.5;
                }

                break;
            }
            default:
            {
                for (j = 0; j < y_num; j++)
                {
                    value = j / (double) (y_num - 1);
                    for (i = j * x_num; i < (j + 1) * x_num; i++)
                    {
                        node_xy[1 + i * 2] = value;
                    }
                }

                break;
            }
        }

        return node_xy;
    }

    public static void grid_print(int element_order, int element_num, int[] element_node)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    GRID_PRINT prints the elements that form a grid.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ELEMENT_ORDER, the number of nodes per element.
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        // 
        //    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the nodes that form
        //    each element.
        //
    {
        int element;
        int i;

        Console.WriteLine("");
        Console.WriteLine("  GRID_PRINT: Element -> Node table.");
        Console.WriteLine("");
        Console.WriteLine("  Element order = " + element_order + "");
        Console.WriteLine("  Number of elements = " + element_num + "");
        Console.WriteLine("");
        string cout = "    #   ";
        for (i = 0; i < element_order; i++)
        {
            cout += i.ToString(CultureInfo.InvariantCulture).PadLeft(3);
        }

        Console.WriteLine(cout);
        Console.WriteLine("");

        for (element = 0; element < element_num; element++)
        {
            string cout1 = "  " +element.ToString(CultureInfo.InvariantCulture).PadLeft(3) + "   ";
            for (i = 0; i < element_order; i++)
            {
                cout1 += element_node[i + element * element_order].ToString(CultureInfo.InvariantCulture).PadLeft(3);
            }

            Console.WriteLine(cout1);
        }
    }

    public static int[] grid_q4_element(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    GRID_Q4_ELEMENT produces a grid of 4 node quadrilaterals.
        //
        //  Discussion:
        //
        //    For each element, the nodes are listed in counter-clockwise order.
        //
        //  Example:
        //
        //    Input:
        //
        //      NELEMX = 3, NELEMY = 2
        //
        //    Output:
        //
        //      ELEMENT_NODE = 
        //         1, 2,  6,  5;
        //         2, 3,  7,  6;
        //         3, 4,  8,  7;
        //         5, 6, 10,  9;
        //         6, 7, 11, 10;
        //         7, 8, 12, 11.
        //
        //  Grid:
        //
        //    9---10---11---12
        //    |    |    |    |
        //    |    |    |    |
        //    |  4 |  5 |  6 |
        //    |    |    |    |
        //    5----6----7----8
        //    |    |    |    |
        //    |    |    |    |
        //    |  1 |  2 |  3 |
        //    |    |    |    |
        //    1----2----3----4
        //
        //  Element Q4:
        //
        //    |
        //    1  4-----3
        //    |  |     |
        //    |  |     |
        //    S  |     |
        //    |  |     |
        //    |  |     |
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
        //    06 April 2005
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
        //    Output, int GRID_Q4[4*NELEMX*NELEMY], the nodes that form
        //    each element.
        //
    {
        const int element_order = 4;
        int j;

        int[] element_node = new int[element_order * nelemx * nelemy];
        //
        //  Node labeling:
        //
        //    NW---NE
        //     |    |
        //    SW---SE
        //
        int element = 0;

        for (j = 1; j <= nelemy; j++)
        {
            int i;
            for (i = 1; i <= nelemx; i++)
            {
                int sw = i + (j - 1) * (nelemx + 1);
                int se = i + 1 + (j - 1) * (nelemx + 1);
                int nw = i + j * (nelemx + 1);
                int ne = i + 1 + j * (nelemx + 1);

                element_node[0 + element * element_order] = sw;
                element_node[1 + element * element_order] = se;
                element_node[2 + element * element_order] = ne;
                element_node[3 + element * element_order] = nw;

                element += 1;
            }
        }

        return element_node;
    }

    public static int grid_q4_element_num(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRID_Q4_ELEMENT_NUM counts the elements in a grid of 4 node quadrilaterals.
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
        //    23 March 2006
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
        //    Output, int GRID_Q4_ELEMENT_NUM, the number of elements in the grid.
        //
    {
        int element_num = nelemx * nelemy;

        return element_num;
    }

    public static int grid_q4_node_num(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRID_Q4_node_num = counts the nodes in a grid of 4 node quadrilaterals.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 March 2006
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
        //    Output, int GRID_Q4_node_num, the number of nodes in the grid.
        //
    {
        int node_num = (nelemx + 1) * (nelemy + 1);

        return node_num;
    }

    public static int[] grid_q8_element(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    GRID_Q8_ELEMENT produces a grid of 8 node quadrilaterals.
        //
        //  Example:
        //
        //    Input:
        //
        //      NELEMX = 3, NELEMY = 2
        //
        //    Output:
        //
        //      ELEMENT_NODE =
        //         1,  3, 14, 12,  2,  9, 13,  8;
        //         3,  5, 16, 14,  4, 10, 15,  9;
        //         5,  7, 18, 16,  6, 11, 17, 10;
        //        12, 14, 25, 23, 13, 20, 24, 19;
        //        14, 16, 27, 25, 15, 21, 26, 20;
        //        16, 18, 29, 27, 17, 22, 28, 21.
        //
        //  Diagram:
        //
        //   23---24---25---26---27---28---29
        //    |         |         |         |
        //    |         |         |         |
        //   19        20        21        22
        //    |         |         |         |
        //    | 4       | 5       | 6       |
        //   12---13---14---15---16---17---18
        //    |         |         |         |
        //    |         |         |         |
        //    8         9        10        11
        //    |         |         |         |
        //    | 1       | 2       | 3       |
        //    1----2----3----4----5----6----7
        //
        //  Element Q8:
        //
        //    |
        //    1  4--7--3
        //    |  |     |
        //    |  |     |
        //    S  8     6
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
        //    06 April 2005
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
        //    Output, int GRID_Q8[8*NELEMX*NELEMY], the nodes that form
        //    each element.
        //
    {
        const int element_order = 8;
        int j;

        int[] element_node = new int[element_order * nelemx * nelemy];
        //
        //  Node labeling:
        //
        //    NW----N----NE
        //     |          |
        //     W   (C)    E
        //     |          |
        //    SW----S----SE
        //

        int element = 0;

        for (j = 1; j <= nelemy; j++)
        {
            int i;
            for (i = 1; i <= nelemx; i++)
            {
                int sw = (j - 1) * (3 * nelemx + 2) + 2 * i - 1;
                int w = sw + 2 * nelemx + 2 - i;
                int nw = sw + 3 * nelemx + 2;

                int s = sw + 1;
                int n = sw + 3 * nelemx + 2 + 1;

                int se = sw + 2;
                int e = sw + 2 * nelemx + 2 - i + 1;
                int ne = sw + 3 * nelemx + 2 + 2;

                element_node[0 + element * element_order] = sw;
                element_node[1 + element * element_order] = se;
                element_node[2 + element * element_order] = ne;
                element_node[3 + element * element_order] = nw;
                element_node[4 + element * element_order] = s;
                element_node[5 + element * element_order] = e;
                element_node[6 + element * element_order] = n;
                element_node[7 + element * element_order] = w;

                element += 1;
            }
        }

        return element_node;
    }

    public static int grid_q8_element_num(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRID_Q8_ELEMENT_NUM counts the elements in a grid of 8 node quadrilaterals.
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
        //    23 March 2006
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
        //    Output, int GRID_Q8_ELEMENT_NUM, the number of elements in the grid.
        //
    {
        int element_num = nelemx * nelemy;

        return element_num;
    }

    public static int grid_q8_node_num(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRID_Q8_node_num = counts the nodes in a grid of 8 node quadrilaterals.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 March 2006
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
        //    Output, int GRID_Q8_node_num, the number of nodes in the grid.
        //
    {
        int node_num = 3 * nelemx * nelemy + 2 * nelemx + 2 * nelemy + 1;

        return node_num;
    }

    public static int[] grid_q9_element(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    GRID_Q9_ELEMENT produces a grid of 9 node quadrilaterals.
        //
        //  Example:
        //
        //    Input:
        //
        //      NELEMX = 3, NELEMY = 2
        //
        //    Output:
        //
        //      ELEMENT_NODE = 
        //         1,  3, 17, 15,  2, 10, 16,  8,  9;
        //         3,  5, 19, 17,  4, 12, 18, 10, 11;
        //         5,  7, 21, 19,  6, 14, 20, 12, 13;
        //        15, 17, 31, 29, 16, 24, 30, 22, 23;
        //        17, 19, 33, 31, 18, 26, 32, 24, 25;
        //        19, 21, 35, 33, 20, 28, 34, 26, 27.
        //
        //  Grid:
        //
        //   29---30---31---32---33---34---35
        //    |    .    |    .    |    .    |
        //    |    .    |    .    |    .    |
        //   22 . 23 . 24 . 25 . 26 . 27 . 28
        //    |    .    |    .    |    .    |
        //    | 4  .    | 5  .    | 6  .    |
        //   15---16---17---18---19---20---21
        //    |    .    |    .    |    .    |
        //    |    .    |    .    |    .    |
        //    8 .  9 . 10 . 11 . 12 . 13 . 14
        //    |    .    |    .    |    .    |
        //    | 1  .    | 2  .    | 3  .    |
        //    1----2----3----4----5----6----7
        //
        //  Element Q9:
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
        //    06 April 2005
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
        //    Output, int GRID_Q9[9*NELEMX*NELEMY], the nodes that form
        //    each element.
        //
    {
        const int element_order = 9;
        int j;

        int[] element_node = new int[element_order * nelemx * nelemy];
        //
        //  Node labeling:
        //
        //    NW----N----NE
        //     |          |
        //     W    C     E
        //     |          |
        //    SW----S----SE
        //
        int element = 0;

        for (j = 1; j <= nelemy; j++)
        {
            int i;
            for (i = 1; i <= nelemx; i++)
            {
                int sw = 2 * (j - 1) * (2 * nelemx + 1) + 2 * (i - 1) + 1;
                int w = sw + 2 * nelemx + 1;
                int nw = sw + 2 * (2 * nelemx + 1);

                int s = sw + 1;
                int c = sw + 1 + 2 * nelemx + 1;
                int n = sw + 1 + 2 * (2 * nelemx + 1);

                int se = sw + 2;
                int e = sw + 2 + 2 * nelemx + 1;
                int ne = sw + 2 + 2 * (2 * nelemx + 1);

                element_node[0 + element * element_order] = sw;
                element_node[1 + element * element_order] = se;
                element_node[2 + element * element_order] = ne;
                element_node[3 + element * element_order] = nw;
                element_node[4 + element * element_order] = s;
                element_node[5 + element * element_order] = e;
                element_node[6 + element * element_order] = n;
                element_node[7 + element * element_order] = w;
                element_node[8 + element * element_order] = c;

                element += 1;
            }
        }

        return element_node;
    }

    public static int grid_q9_element_num(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRID_Q9_ELEMENT_NUM counts the elements in a grid of 9 node quadrilaterals.
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
        //    23 March 2006
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
        //    Output, int GRID_Q9_ELEMENT_NUM, the number of elements in the grid.
        //
    {
        int element_num = nelemx * nelemy;

        return element_num;
    }

    public static int grid_q9_node_num(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRID_Q9_node_num = counts the nodes in a grid of 9 node quadrilaterals.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 March 2006
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
        //    Output, int GRID_Q9_node_num, the number of nodes in the grid.
        //
    {
        int node_num = (2 * nelemx + 1) * (2 * nelemy + 1);

        return node_num;
    }

    public static int[] grid_q12_element(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    GRID_Q12_ELEMENT produces a grid of 12 node quadrilaterals.
        //
        //  Example:
        //
        //    Input:
        //
        //      NELEMX = 3, NELEMY = 2
        //
        //    Output:
        //
        //      ELEMENT_NODE =
        //         1,  2,  3,  4, 11, 12, 15, 16, 19, 20, 21, 22;
        //         4,  5,  6,  7, 12, 13, 16, 17, 22, 23, 24, 25;
        //         7,  8,  9, 10, 13, 14, 17, 18, 25, 26, 27, 28;
        //        19, 20, 21, 22, 29, 30, 33, 34, 37, 38, 39, 40;
        //        22, 23, 24, 25, 30, 31, 34, 35, 40, 41, 42, 43;
        //        25, 26, 27, 28, 31, 32, 35, 36, 43, 44, 45, 46.
        //
        //  Grid:
        //
        //   37-38-39-40-41-42-43-44-45-46
        //    |        |        |        |
        //   33       34       35       36
        //    |        |        |        |
        //   29       30       31       32
        //    | 4      | 5      | 6      |
        //   19-20-21-22-23-24-25-26-27-28
        //    |        |        |        |
        //   15       16       17       18
        //    |        |        |        |
        //   11       12       13       14
        //    | 1      | 2      | 3      |
        //    1--2--3--4--5--6--7--8--9-10
        //
        //  Element Q12:
        //
        //    |
        //    1  9-10-11-12
        //    |  |        |
        //    |  7        8
        //    S  |        |
        //    |  5        6
        //    |  |        |
        //    0  1--2--3--4
        //    |
        //    +--0---R---1-->
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 April 2005
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
        //    Output, int GRID_Q12[12*NELEMX*NELEMY], the nodes that form
        //    each element.
        //
    {
        const int element_order = 12;
        int j;

        int[] element_node = new int[element_order * nelemx * nelemy];

        int element = 0;

        for (j = 1; j <= nelemy; j++)
        {
            int i;
            for (i = 1; i <= nelemx; i++)
            {
                int base_ = (j - 1) * (5 * nelemx + 3) + 1;

                element_node[0 + element * element_order] = base_ + (i - 1) * 3;
                element_node[1 + element * element_order] = base_ + (i - 1) * 3 + 1;
                element_node[2 + element * element_order] = base_ + (i - 1) * 3 + 2;
                element_node[3 + element * element_order] = base_ + (i - 1) * 3 + 3;

                element_node[4 + element * element_order] = base_ + 3 * nelemx + i;
                element_node[5 + element * element_order] = base_ + 3 * nelemx + i + 1;

                element_node[6 + element * element_order] = base_ + 4 * nelemx + i + 1;
                element_node[7 + element * element_order] = base_ + 4 * nelemx + i + 2;

                element_node[8 + element * element_order] = base_ + 5 * nelemx + 3 * i;
                element_node[9 + element * element_order] = base_ + 5 * nelemx + 3 * i + 1;
                element_node[10 + element * element_order] = base_ + 5 * nelemx + 3 * i + 2;
                element_node[11 + element * element_order] = base_ + 5 * nelemx + 3 * i + 3;

                element += 1;
            }
        }

        return element_node;
    }

    public static int grid_q12_element_num(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRID_Q12_ELEMENT_NUM counts the elements in a grid of 12 node quadrilaterals.
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
        //    23 March 2006
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
        //    Output, int GRID_Q12_ELEMENT_NUM, the number of elements in the grid.
        //
    {
        int element_num = nelemx * nelemy;

        return element_num;
    }

    public static int grid_q12_node_num(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRID_Q12_node_num = counts the nodes in a grid of 12 node quadrilaterals.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 March 2006
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
        //    Output, int GRID_Q12_node_num, the number of nodes in the grid.
        //
    {
        int node_num = 5 * nelemx * nelemy + 3 * nelemx + 3 * nelemy + 1;

        return node_num;
    }

    public static int[] grid_q16_element(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    GRID_Q16_ELEMENT produces a grid of 16 node quadrilaterals.
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
        //         1,  2,  3,  4,  8,  9, 10, 11, 15, 16, 17, 18, 22, 23, 24, 25;
        //         4,  5,  6,  7, 11, 12, 13, 14, 18, 19, 20, 21, 25, 26, 27, 28;
        //        22, 23, 24, 25, 29, 30, 31, 32, 36, 37, 38, 39, 43, 44, 45, 46;
        //        25, 26, 27, 28, 32, 33, 34, 35, 39, 40, 41, 42, 46, 47, 48, 49. 
        //        
        //  Grid:
        //
        //   43-44-45-46-47-48-49
        //    |        |        |
        //    |        |        |
        //   36 37 38 39 40 41 42
        //    |        |        |
        //    |        |        |
        //   29 30 31 32 33 34 35
        //    |        |        |
        //    | 3      | 4      |
        //   22-23-24-25-26-27-28
        //    |        |        |
        //    |        |        |
        //   15 16 17 18 19 20 21
        //    |        |        |
        //    |        |        |
        //    8  9 10 11 12 13 14
        //    |        |        |
        //    | 1      | 2      |
        //    1--2--3--4--5--6--7
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 April 2005
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
        //    Output, int GRID_Q16[16*NELEMX*NELEMY], the nodes that form
        //    each element.
        //
    {
        const int element_order = 16;
        int j;

        int[] element_node = new int[element_order * nelemx * nelemy];

        int element = 0;

        for (j = 1; j <= nelemy; j++)
        {
            int i;
            for (i = 1; i <= nelemx; i++)
            {
                int base_ = (j - 1) * 3 * (3 * nelemx + 1) + 3 * i - 2;

                element_node[0 + element * element_order] = base_;
                element_node[1 + element * element_order] = base_ + 1;
                element_node[2 + element * element_order] = base_ + 2;
                element_node[3 + element * element_order] = base_ + 3;
                element_node[4 + element * element_order] = base_ + 3 * nelemx + 1;
                element_node[5 + element * element_order] = base_ + 3 * nelemx + 1 + 1;
                element_node[6 + element * element_order] = base_ + 3 * nelemx + 1 + 2;
                element_node[7 + element * element_order] = base_ + 3 * nelemx + 1 + 3;
                element_node[8 + element * element_order] = base_ + 2 * (3 * nelemx + 1);
                element_node[9 + element * element_order] = base_ + 2 * (3 * nelemx + 1) + 1;
                element_node[10 + element * element_order] = base_ + 2 * (3 * nelemx + 1) + 2;
                element_node[11 + element * element_order] = base_ + 2 * (3 * nelemx + 1) + 3;
                element_node[12 + element * element_order] = base_ + 3 * (3 * nelemx + 1);
                element_node[13 + element * element_order] = base_ + 3 * (3 * nelemx + 1) + 1;
                element_node[14 + element * element_order] = base_ + 3 * (3 * nelemx + 1) + 2;
                element_node[15 + element * element_order] = base_ + 3 * (3 * nelemx + 1) + 3;

                element += 1;
            }
        }

        return element_node;
    }

    public static int grid_q16_element_num(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRID_Q16_ELEMENT_NUM counts the elements in a grid of 16 node quadrilaterals.
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
        //    23 March 2006
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
        //    Output, int GRID_Q16_ELEMENT_NUM, the number of elements in the grid.
        //
    {
        int element_num = nelemx * nelemy;

        return element_num;
    }

    public static int grid_q16_node_num(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRID_Q16_node_num = counts the nodes in a grid of 16 node quadrilaterals.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 March 2006
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
        //    Output, int GRID_Q16_node_num, the number of nodes in the grid.
        //
    {
        int node_num = (3 * nelemx + 1) * (3 * nelemy + 1);

        return node_num;
    }

    public static int[] grid_ql_element(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    GRID_QL_ELEMENT produces a grid of 6 node quadratics/linears.
        //
        //  Example:
        //
        //    Input:
        //
        //      NELEMX = 3, NELEMY = 2
        //
        //    Output:
        //
        //      ELEMENT_NODE = 
        //         1,  2,  3,  8,  9, 10;
        //         3,  4,  5, 10, 11, 12;
        //         5,  6,  7, 12, 13, 14;
        //         8,  9, 10, 15, 16, 17;
        //        10, 11, 12, 17, 18, 19;
        //        12, 13, 14, 19, 20, 21.
        //
        //  Grid:
        //
        //   15---16---17---18---19---20---21
        //    |         |         |         |
        //    |         |         |         |
        //    |    4    |    5    |    6    |
        //    |         |         |         |
        //    |         |         |         |
        //    8----9---10---11---12---13---14
        //    |         |         |         |
        //    |         |         |         |
        //    |    1    |    2    |    3    |
        //    |         |         |         |
        //    |         |         |         |
        //    1----2----3----4----5----6----7
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 April 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NELEMX, NELEMY, the number of elements along the
        //    X and Y directions.  The number of elements generated will be
        //    NELEMX * NELEMY.  X will the the "quadratic direction", and
        //    Y will be the "linear direction".
        //
        //    Output, int GRID_QL[6*NELEMX*NELEMY], the nodes that form
        //    each element.
        //
    {
        const int element_order = 6;
        int j;

        int[] element_node = new int[element_order * nelemx * nelemy];

        int element = 0;

        for (j = 1; j <= nelemy; j++)
        {
            int i;
            for (i = 1; i <= nelemx; i++)
            {
                int base_ = (j - 1) * (2 * nelemx + 1) + 2 * i - 1;

                element_node[0 + element * element_order] = base_;
                element_node[1 + element * element_order] = base_ + 1;
                element_node[2 + element * element_order] = base_ + 2;
                element_node[3 + element * element_order] = base_ + 2 * nelemx + 1;
                element_node[4 + element * element_order] = base_ + 2 * nelemx + 1 + 1;
                element_node[5 + element * element_order] = base_ + 2 * nelemx + 1 + 2;

                element += 1;
            }
        }

        return element_node;
    }

    public static int grid_ql_element_num(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRID_QL_ELEMENT_NUM counts the elements in a grid of quadratic/linear quadrilaterals.
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
        //    23 March 2006
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
        //    Output, int GRID_QL_ELEMENT_NUM, the number of elements in the grid.
        //
    {
        int element_num = nelemx * nelemy;

        return element_num;
    }

    public static int grid_ql_node_num(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRID_QL_node_num = counts the nodes in a grid of quadratic/linear quadrilaterals.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 March 2006
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
        //    Output, int GRID_QL_node_num, the number of nodes in the grid.
        //
    {
        int node_num = 2 * nelemx * nelemy + 2 * nelemx + nelemy + 1;

        return node_num;
    }

    public static void grid_shape_2d(int n, double[] a, ref int n1, ref int n2 )

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    GRID_SHAPE_2D guesses the shape N1 by N2 of a vector of data.
        //
        //  Discussion:
        //
        //    The data vector A is assumed to contain N1 * N2 values, with
        //    where each of N2 values is repeated N1 times.
        //
        //  Example:
        //
        //    Input:
        //
        //      A = ( 2, 2, 2, 7, 7, 7 )
        //
        //    Output:
        //
        //      N1 = 3, N2 = 2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of data values.
        //
        //    Input, double A[N], the data, which should have the properties
        //    described above.
        //
        //    Output, int *N1, *N2, the "shape" of the data in the array.
        //
    {
        //
        //  Make a guess for N1.
        //
        int i = 1;
        n1 = 1;

        for (i = 1; i < n; i++)
        {
            if (Math.Abs(a[i] - a[0]) > double.Epsilon)
            {
                break;
            }

            n1 += 1;
        }

        //
        //  Guess that N2 = N / N1.
        //
        n2 = n / n1;
    }

    public static int[] grid_t3_element(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    GRID_T3_ELEMENT produces a grid of pairs of 3 node triangles.
        //
        //  Example:
        //
        //    Input:
        //
        //      NELEMX = 3, NELEMY = 2
        //
        //    Output:
        //
        //      ELEMENT_NODE = 
        //         1,  2,  5;
        //         6,  5,  2;
        //         2,  3,  6;
        //         7,  6,  3;
        //         3,  4,  7;
        //         8,  7,  4;
        //         5,  6,  9;
        //        10,  9,  6;
        //         6,  7, 10;
        //        11, 10,  7;
        //         7,  8, 11;
        //        12, 11,  8.
        //
        //  Grid:
        //
        //    9---10---11---12
        //    |. 8 |.10 |.12 |
        //    | .  | .  | .  |
        //    |  . |  . |  . |
        //    |  7.|  9.| 11.|
        //    5----6----7----8
        //    |. 2 |. 4 |. 6 |
        //    | .  | .  | .  |
        //    |  . |  . |  . |
        //    |  1.|  3.|  5.|
        //    1----2----3----4
        //
        //  Element T3:
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
        //    06 April 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NELEMX, NELEMY, the number of elements along the
        //    X and Y directions.  The number of elements generated will be
        //    2 * NELEMX * NELEMY.
        //
        //    Output, int GRID_T3[3*2*NELEMX*NELEMY], the nodes that form
        //    each element.
        //
    {
        const int element_order = 3;
        int j;

        int[] element_node = new int[element_order * 2 * nelemx * nelemy];
        //
        //  Node labeling:
        //
        //    NW--NE
        //     |\ |
        //     | \|
        //    SW--SE
        //
        int element = 0;

        for (j = 1; j <= nelemy; j++)
        {
            int i;
            for (i = 1; i <= nelemx; i++)
            {
                int sw = i + (j - 1) * (nelemx + 1);
                int se = i + 1 + (j - 1) * (nelemx + 1);
                int nw = i + j * (nelemx + 1);
                int ne = i + 1 + j * (nelemx + 1);

                element_node[0 + element * element_order] = sw;
                element_node[1 + element * element_order] = se;
                element_node[2 + element * element_order] = nw;
                element += 1;

                element_node[0 + element * element_order] = ne;
                element_node[1 + element * element_order] = nw;
                element_node[2 + element * element_order] = se;
                element += 1;
            }
        }

        return element_node;
    }

    public static int grid_t3_element_num(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRID_T3_ELEMENT_NUM counts the elements in a grid of 3 node triangles.
        //
        //  Example:
        //
        //    Input:
        //
        //      NELEMX = 3, NELEMY = 2
        //
        //    Output:
        //
        //      ELEMENT_NUM = 2 * NELEMX * NELEMY = 12
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 March 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NELEMX, NELEMY, the number of elements along the
        //    X and Y directions.  The number of elements generated will be
        //    2 * NELEMX * NELEMY.
        //
        //    Output, int GRID_T3_ELEMENT_NUM, the number of elements in the grid.
        //
    {
        int element_num = 2 * nelemx * nelemy;

        return element_num;
    }

    public static int grid_t3_node_num(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRID_T3_node_num = counts the nodes in a grid of 3 node triangles.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 March 2006
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
        //    Output, int GRID_T3_node_num, the number of nodes in the grid.
        //
    {
        int node_num = (nelemx + 1) * (nelemy + 1);

        return node_num;
    }

    public static int[] grid_t4_element(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    GRID_T4_ELEMENT produces a grid of pairs of 4 node triangles.
        //
        //  Example:
        //
        //    Input:
        //
        //      NELEMX = 3, NELEMY = 2
        //
        //    Output:
        //
        //      ELEMENT_NODE = 
        //         1,  2,  11,  5;
        //        12, 11,   2,  8;
        //         2,  3,  12,  6;
        //        13, 12,   3,  9;
        //         3   4   13,  7;
        //        14, 13,   4,  10;
        //        11, 12,  21,  15;
        //        22, 21,  12,  18;
        //        12, 13,  22,  16;
        //        23, 22,  13,  19;
        //        13  14   23,  17;
        //        24, 23,  14,  20;
        //
        //  Grid:
        //
        //   21---22---23---24
        //    |.18 |.19 |.20 |
        //    | .  | .  | .  |
        //    |  . |  . |  . |
        //    | 15.| 16.| 17.|
        //   11---12---13---14
        //    |. 8 |. 9 |.10 |
        //    | .  | .  | .  |
        //    |  . |  . |  . |
        //    | 5 .|  6.|  7.|
        //    1----2----3----4
        //
        //  Element T4:
        //
        //    |
        //    1  3
        //    |  ..
        //    |  . .
        //    S  .  .
        //    |  . 4 .
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
        //    23 March 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NELEMX, NELEMY, the number of elements along the
        //    X and Y directions.  The number of elements generated will be
        //    2 * NELEMX * NELEMY.
        //
        //    Output, int GRID_T4[4*2*NELEMX*NELEMY], the nodes that form
        //    each element.
        //
    {
        const int element_order = 4;
        int j;

        int[] element_node = new int[element_order * 2 * nelemx * nelemy];
        //
        //  Node labeling:
        //
        //    NW----NE
        //     |.   |
        //     | .NC|
        //     |SC. |
        //     |   .|
        //    SW---SE
        //
        int element = 0;

        for (j = 1; j <= nelemy; j++)
        {
            int i;
            for (i = 1; i <= nelemx; i++)
            {
                int sw = i + (j - 1) * (3 * nelemx + 1);
                int se = sw + 1;
                int sc = sw + nelemx + 1;
                int nc = sw + 2 * nelemx + 1;
                int nw = sw + 3 * nelemx + 1;
                int ne = sw + 3 * nelemx + 2;

                element_node[0 + element * element_order] = sw;
                element_node[1 + element * element_order] = se;
                element_node[2 + element * element_order] = nw;
                element_node[3 + element * element_order] = sc;
                element += 1;

                element_node[0 + element * element_order] = ne;
                element_node[1 + element * element_order] = nw;
                element_node[2 + element * element_order] = se;
                element_node[3 + element * element_order] = nc;
                element += 1;
            }
        }

        return element_node;
    }

    public static int grid_t4_element_num(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRID_T4_ELEMENT_NUM counts the elements in a grid of 4 node triangles.
        //
        //  Example:
        //
        //    Input:
        //
        //      NELEMX = 3, NELEMY = 2
        //
        //    Output:
        //
        //      ELEMENT_NUM = 2 * NELEMX * NELEMY = 12
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 March 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NELEMX, NELEMY, the number of elements along the
        //    X and Y directions.  The number of elements generated will be
        //    2 * NELEMX * NELEMY.
        //
        //    Output, int GRID_T4_ELEMENT_NUM, the number of elements in the grid.
        //
    {
        int element_num = 2 * nelemx * nelemy;

        return element_num;
    }

    public static int grid_t4_node_num(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRID_T4_node_num = counts the nodes in a grid of 4 node triangles.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 March 2006
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
        //    Output, int GRID_T4_node_num, the number of nodes in the grid.
        //
    {
        int node_num = (nelemx + 1) * (nelemy + 1) + 2 * nelemx * nelemy;

        return node_num;
    }

    public static int[] grid_t6_element(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    GRID_T6_ELEMENT produces a grid of pairs of 6 node triangles.
        //
        //  Example:
        //
        //    Input:
        //
        //      NELEMX = 3, NELEMY = 2
        //
        //    Output:
        //
        //      ELEMENT_NODE = 
        //         1,  3, 15,  2,  9,  8;
        //        17, 15,  3, 16,  9, 10;
        //         3,  5, 17,  4, 11, 10;
        //        19, 17,  5, 18, 11, 12;
        //         5,  7, 19,  6, 13, 12;
        //        21, 19,  7, 20, 13, 14;
        //        15, 17, 29, 16, 23, 22;
        //        31, 29, 17, 30, 23, 24;
        //        17, 19, 31, 18, 25, 24;
        //        33, 31, 19, 32, 25, 26;
        //        19, 21, 33, 20, 27, 26;
        //        35, 33, 21, 34, 27, 28.
        //
        //  Grid:
        //
        //   29-30-31-32-33-34-35
        //    |. 8  |.10  |.12  |
        //    | .   | .   | .   |
        //   22 23 24 25 26 27 28
        //    |   . |   . |   . |
        //    |  7 .|  9 .| 11 .|
        //   15-16-17-18-19-20-21
        //    |. 2  |. 4  |. 6  |
        //    | .   | .   | .   |
        //    8  9 10 11 12 13 14
        //    |   . |   . |   . |
        //    |  1 .|  3 .|  5 .|
        //    1--2--3--4--5--6--7
        //
        //  Element T6:
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
        //    06 April 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NELEMX, NELEMY, the number of elements along the
        //    X and Y directions.  The number of elements generated will be
        //    2 * NELEMX * NELEMY.
        //
        //    Output, int GRID_T6[6*2*NELEMX*NELEMY], the nodes that form
        //    each element.
        //
    {
        const int element_order = 6;
        int j;

        int[] element_node = new int[element_order * 2 * nelemx * nelemy];
        //
        //  Node labeling:
        //
        //    NW---N--NE
        //     | .     |
        //     W   C   E
        //     |    .  |
        //    SW---S--SE
        //
        int element = 0;

        for (j = 1; j <= nelemy; j++)
        {
            int i;
            for (i = 1; i <= nelemx; i++)
            {
                int sw = 2 * (j - 1) * (2 * nelemx + 1) + 2 * (i - 1) + 1;
                int w = sw + 2 * nelemx + 1;
                int nw = sw + 2 * (2 * nelemx + 1);

                int s = sw + 1;
                int c = sw + 1 + 2 * nelemx + 1;
                int n = sw + 1 + 2 * (2 * nelemx + 1);

                int se = sw + 2;
                int e = sw + 2 + 2 * nelemx + 1;
                int ne = sw + 2 + 2 * (2 * nelemx + 1);

                element_node[0 + element * element_order] = sw;
                element_node[1 + element * element_order] = se;
                element_node[2 + element * element_order] = nw;
                element_node[3 + element * element_order] = s;
                element_node[4 + element * element_order] = c;
                element_node[5 + element * element_order] = w;
                element += 1;

                element_node[0 + element * element_order] = ne;
                element_node[1 + element * element_order] = nw;
                element_node[2 + element * element_order] = se;
                element_node[3 + element * element_order] = n;
                element_node[4 + element * element_order] = c;
                element_node[5 + element * element_order] = e;
                element += 1;
            }
        }

        return element_node;
    }

    public static int grid_t6_element_num(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRID_T6_ELEMENT_NUM counts the elements in a grid of 6 node triangles.
        //
        //  Example:
        //
        //    Input:
        //
        //      NELEMX = 3, NELEMY = 2
        //
        //    Output:
        //
        //      ELEMENT_NUM = 2 * NELEMX * NELEMY = 12
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 March 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NELEMX, NELEMY, the number of elements along the
        //    X and Y directions.  The number of elements generated will be
        //    2 * NELEMX * NELEMY.
        //
        //    Output, int GRID_T6_ELEMENT_NUM, the number of elements in the grid.
        //
    {
        int element_num = 2 * nelemx * nelemy;

        return element_num;
    }

    public static int grid_t6_node_num(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRID_T6_node_num = counts the nodes in a grid of 6 node triangles.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 March 2006
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
        //    Output, int GRID_T6_node_num, the number of nodes in the grid.
        //
    {
        int node_num = (2 * nelemx + 1) * (2 * nelemy + 1);

        return node_num;
    }

    public static int[] grid_t10_element(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    GRID_T10_ELEMENT produces a grid of pairs of 10 node triangles.
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
        //         1,  2,  3,  4, 10, 16, 22, 15,  8,  9;
        //        25, 24, 23, 22, 16, 10,  4, 11, 18, 17;
        //         4,  5,  6,  7, 13, 19, 25, 18, 11, 12;
        //        28, 27, 26, 25, 19, 13,  7, 14, 21, 20;
        //        22, 23, 24, 25, 31, 37, 43, 36, 29, 30;
        //        46, 45, 44, 43, 37, 31, 25, 32, 39, 38;
        //        25, 26, 27, 28, 34, 40, 46, 39, 31, 33;
        //        49, 48, 47, 46, 40, 34, 28, 35, 42, 41.
        //        
        //  Grid:
        //
        //   43-44-45-46-47-48-49
        //    |\     6 |\     8 |
        //    | \      | \      |
        //   36 37 38 39 40 41 42
        //    |   \    |   \    |
        //    |    \   |    \   |
        //   29 30 31 32 33 34 35
        //    |      \ |      \ |
        //    | 5     \| 7     \|
        //   22-23-24-25-26-27-28
        //    |\     2 |\     4 |
        //    | \      | \      |
        //   15 16 17 18 19 20 21
        //    |   \    |   \    |
        //    |    \   |    \   |
        //    8  9 10 11 12 13 14
        //    |      \ |      \ |
        //    | 1     \| 3     \|
        //    1--2--3--4--5--6--7
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 April 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NELEMX, NELEMY, the number of elements along the
        //    X and Y directions.  The number of elements generated will be
        //    2 * NELEMX * NELEMY.
        //
        //    Output, int GRID_T10[10*2*NELEMX*NELEMY], the nodes that form
        //    each element.
        //
    {
        const int element_order = 10;
        int j;

        int[] element_node = new int[element_order * 2 * nelemx * nelemy];

        int element = 0;

        for (j = 1; j <= nelemy; j++)
        {
            int i;
            for (i = 1; i <= nelemx; i++)
            {
                int base_ = (j - 1) * 3 * (3 * nelemx + 1) + 3 * i - 2;

                element_node[0 + element * element_order] = base_;
                element_node[1 + element * element_order] = base_ + 1;
                element_node[2 + element * element_order] = base_ + 2;
                element_node[3 + element * element_order] = base_ + 3;
                element_node[4 + element * element_order] = base_ + 3 * nelemx + 1 + 2;
                element_node[5 + element * element_order] = base_ + 2 * (3 * nelemx + 1) + 1;
                element_node[6 + element * element_order] = base_ + 3 * (3 * nelemx + 1);
                element_node[7 + element * element_order] = base_ + 2 * (3 * nelemx + 1);
                element_node[8 + element * element_order] = base_ + 2 * nelemx + 1 + 2;
                element_node[9 + element * element_order] = base_ + 2 * nelemx + 1 + 3;
                element += 1;

                element_node[0 + element * element_order] = base_ + 3 * (3 * nelemx + 1) + 3;
                element_node[1 + element * element_order] = base_ + 3 * (3 * nelemx + 1) + 2;
                element_node[2 + element * element_order] = base_ + 3 * (3 * nelemx + 1) + 1;
                element_node[3 + element * element_order] = base_ + 3 * (3 * nelemx + 1);
                element_node[4 + element * element_order] = base_ + 2 * (3 * nelemx + 1) + 1;
                element_node[5 + element * element_order] = base_ + 3 * nelemx + 1 + 2;
                element_node[6 + element * element_order] = base_ + 3;
                element_node[7 + element * element_order] = base_ + 3 * nelemx + 1 + 3;
                element_node[8 + element * element_order] = base_ + 2 * (3 * nelemx + 1) + 3;
                element_node[9 + element * element_order] = base_ + 2 * (3 * nelemx + 1) + 2;
                element += 1;
            }
        }

        return element_node;
    }

    public static int grid_t10_element_num(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRID_T10_ELEMENT_NUM counts the elements in a grid of 10 node triangles.
        //
        //  Example:
        //
        //    Input:
        //
        //      NELEMX = 3, NELEMY = 2
        //
        //    Output:
        //
        //      ELEMENT_NUM = 2 * NELEMX * NELEMY = 12
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 March 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NELEMX, NELEMY, the number of elements along the
        //    X and Y directions.  The number of elements generated will be
        //    2 * NELEMX * NELEMY.
        //
        //    Output, int GRID_T10_ELEMENT_NUM, the number of elements in the grid.
        //
    {
        int element_num = 2 * nelemx * nelemy;

        return element_num;
    }

    public static int grid_t10_node_num(int nelemx, int nelemy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRID_T10_node_num = counts the nodes in a grid of 10 node triangles.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 March 2006
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
        //    Output, int GRID_T10_node_num, the number of nodes in the grid.
        //
    {
        int node_num = (3 * nelemx + 1) * (3 * nelemy + 1);

        return node_num;
    }

    public static void grid_t6(int nx, int ny, int nnodes, int element_num, ref int[] element_node)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRID_T6 produces a grid of pairs of 6 node triangles.
        //
        //  Example:
        //
        //    Input:
        //
        //      NX = 4, NY = 3
        //
        //    Output:
        //
        //      ELEMENT_NODE =
        //         1,  3, 15,  2,  9,  8;
        //        17, 15,  3, 16,  9, 10;
        //         3,  5, 17,  4, 11, 10;
        //        19, 17,  5, 18, 11, 12;
        //         5,  7, 19,  6, 13, 12;
        //        21, 19,  7, 20, 13, 14;
        //        15, 17, 29, 16, 23, 22;
        //        31, 29, 17, 30, 23, 24;
        //        17, 19, 31, 18, 25, 24;
        //        33, 31, 19, 32, 25, 26;
        //        19, 21, 33, 20, 27, 26;
        //        35, 33, 21, 34, 27, 28.
        //
        //  Diagram:
        //
        //   29-30-31-32-33-34-35
        //    |\ 8  |\10  |\12  |
        //    | \   | \   | \   |
        //   22 23 24 25 26 27 28
        //    |   \ |   \ |   \ |
        //    |  7 \|  9 \| 11 \|
        //   15-16-17-18-19-20-21
        //    |\ 2  |\ 4  |\ 6  |
        //    | \   | \   | \   |
        //    8  9 10 11 12 13 14
        //    |   \ |   \ |   \ |
        //    |  1 \|  3 \|  5 \|
        //    1--2--3--4--5--6--7
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 September 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NX, NY, controls the number of elements along the
        //    X and Y directions.  The number of elements will be
        //    2 * ( NX - 1 ) * ( NY - 1 ).
        //
        //    Input, int NNODES, the number of local nodes per element.
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        //
        //    Output, int ELEMENT_NODE[NNODES*ELEMENT_NUM];
        //    ELEMENT_NODE(I,J) is the index of the I-th node of the J-th element.
        //
    {
        int j;

        int element = 0;

        for (j = 1; j <= ny - 1; j++)
        {
            int i;
            for (i = 1; i <= nx - 1; i++)
            {
                int sw = (j - 1) * 2 * (2 * nx - 1) + 2 * i - 1;
                int w = sw + 1;
                int nw = sw + 2;

                int s = sw + 2 * nx - 1;
                int c = s + 1;
                int n = s + 2;

                int se = s + 2 * nx - 1;
                int e = se + 1;
                int ne = se + 2;

                element += 1;
                element_node[0 + (element - 1) * nnodes] = sw;
                element_node[1 + (element - 1) * nnodes] = se;
                element_node[2 + (element - 1) * nnodes] = nw;
                element_node[3 + (element - 1) * nnodes] = s;
                element_node[4 + (element - 1) * nnodes] = c;
                element_node[5 + (element - 1) * nnodes] = w;

                element += 1;
                element_node[0 + (element - 1) * nnodes] = ne;
                element_node[1 + (element - 1) * nnodes] = nw;
                element_node[2 + (element - 1) * nnodes] = se;
                element_node[3 + (element - 1) * nnodes] = n;
                element_node[4 + (element - 1) * nnodes] = c;
                element_node[5 + (element - 1) * nnodes] = e;
            }
        }
    }

    public static void grid_test(string code)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    GRID_TEST tests the grid routines.
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
        //    Input, string CODE, the code for the element.
        //    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
        //    'T3', 'T4', 'T6' and 'T10'.
        //
    {
        int element_num = 0;
        //
        //  NODE is defined as a vector rather than a two dimensional array,
        //  so that we can handle the various cases using a single array.
        //
        Console.WriteLine("");
        Console.WriteLine("  GRID_TEST: Test the grid routine for element " + code + "");

        const int nelemx = 3;
        const int nelemy = 2;

        switch (code)
        {
            case "Q4":
            case "Q8":
            case "Q9":
            case "Q12":
            case "Q16":
            case "QL":
                element_num = nelemx * nelemy;
                break;
            case "T3":
            case "T4":
            case "T6":
            case "T10":
                element_num = 2 * nelemx * nelemy;
                break;
        }

        int element_order = Order.order_code(code);

        int[] element_node = grid_element(code, element_order, nelemx, nelemy);

        grid_print(element_order, element_num, element_node);

        int width = grid_width(element_order, element_num, element_node);

        Console.WriteLine("");
        Console.WriteLine("  Grid width is " + width + "");
    }

    public static int grid_width(int element_order, int element_num, int[] element_node)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    GRID_WIDTH computes the width of a given grid.
        //
        //  Definition:
        //
        //    The grid width is defined to the maximum absolute
        //    difference of global indices of nodes in the same element.
        //
        //  Example:
        //
        //    For the following grid, the grid width is 13.
        //
        //   23---24---25---26---27---28---29
        //    |         |         |         |
        //    |         |         |         |
        //   19        20        21        22
        //    |         |         |         |
        //    | 4       | 5       | 6       |
        //   12---13---14---15---16---17---18
        //    |         |         |         |
        //    |         |         |         |
        //    8         9        10        11
        //    |         |         |         |
        //    | 1       | 2       | 3       |
        //    1----2----3----4----5----6----7
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ELEMENT_ORDER, the order of the elements.
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        // 
        //    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the nodes that form
        //    each element.
        //
        //    Output, int GRID_WIDTH, the grid width.
        //
    {
        int element;

        int width = 0;

        for (element = 0; element < element_num; element++)
        {
            int node1;
            for (node1 = 0; node1 < element_order; node1++)
            {
                int ip1 = element_node[node1 + element * element_order];
                int node2;
                for (node2 = 0; node2 < element_order; node2++)
                {
                    int ip2 = element_node[node2 + element * element_order];
                    width = Math.Max(width, Math.Abs(ip1 - ip2));
                }
            }
        }

        return width;
    }

}