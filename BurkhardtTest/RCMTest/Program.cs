using System;
using Burkardt.MatrixNS;
using Burkardt.TriangulationNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace RCMTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for RCM_TEST.
        //
        //  Discussion:
        //
        //    RCM_TEST tests the RCM library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("RCM_TEST");
        Console.WriteLine("  Test the RCM library.");

        test01();
        test02();
        test03();
        test04();
        test05();
        test06();
        test07();
        test08();
        test09();

        test10();
        test11();
        test12();

        Console.WriteLine("");
        Console.WriteLine("RCM_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests ADJ_SET.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int NODE_NUM = 10;
        int ADJ_MAX = NODE_NUM * (NODE_NUM - 1);

        int[] adj = new int[ADJ_MAX];
        int adj_num = 0;
        int[] adj_row = new int[NODE_NUM + 1];
        int i;
        int j;
        int k;
        int n_calls;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  ADJ_SET sets up an adjacency matrix incrementally.");

        n_calls = UniformRNG.i4_uniform(1, ADJ_MAX, ref seed);

        AdjacencyMatrix.adj_set(NODE_NUM, ADJ_MAX, ref adj_num, ref adj_row, ref adj, -1, -1);

        Console.WriteLine("");
        Console.WriteLine("  Creating and recording adjacency information:");
        Console.WriteLine("");

        for (k = 1; k <= n_calls; k++)
        {
            i = UniformRNG.i4_uniform(1, NODE_NUM, ref seed);
            j = UniformRNG.i4_uniform(1, NODE_NUM, ref seed);

            Console.WriteLine("  " + i.ToString().PadLeft(8)
                                   + "  " + j.ToString().PadLeft(8) + "");

            AdjacencyMatrix.adj_set(NODE_NUM, ADJ_MAX, ref adj_num, ref adj_row, ref adj, i, j);
        }

        AdjacencyMatrix.adj_print(NODE_NUM, adj_num, adj_row, adj,
            "  Random adjacency matrix:");

        AdjacencyMatrix.adj_show(NODE_NUM, adj_num, adj_row, adj);

    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests GENRCM;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] adj;
        int adj_num = 0;
        int[] adj_row;
        int bandwidth;
        int i;
        int node_num = 0;
        int[] perm;
        int[] perm_inv;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  GENRCM reorders the nodes in a graph using");
        Console.WriteLine("  the Reverse Cuthill McKee algorithm.");

        Burkardt.Graph.Adjacency.graph_01_size(ref node_num, ref adj_num);

        adj_row = new int[node_num + 1];
        adj = new int[adj_num];
        perm = new int[node_num];
        perm_inv = new int[node_num];

        Burkardt.Graph.Adjacency.graph_01_adj(node_num, adj_num, ref adj_row, ref adj);

        AdjacencyMatrix.adj_print(node_num, adj_num, adj_row, adj, "  Adjacency matrix:");

        AdjacencyMatrix.adj_show(node_num, adj_num, adj_row, adj);

        bandwidth = AdjacencyMatrix.adj_bandwidth(node_num, adj_num, adj_row, adj);

        Console.WriteLine("");
        Console.WriteLine("  ADJ bandwidth = " + bandwidth + "");

        perm = Burkardt.Graph.GenRCM.genrcm(node_num, adj_num, adj_row, adj);

        typeMethods.perm_inverse3(node_num, perm, ref perm_inv);

        Console.WriteLine("");
        Console.WriteLine("  The RCM permutation and inverse:");
        Console.WriteLine("");

        for (i = 0; i < node_num; i++)
        {
            Console.WriteLine("  " + (i + 1).ToString().PadLeft(8)
                                   + "  " + perm[i].ToString().PadLeft(8)
                                   + "  " + perm_inv[i].ToString().PadLeft(8) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Permuted adjacency matrix:");
        Console.WriteLine("");

        AdjacencyMatrix.adj_perm_show(node_num, adj_num, adj_row, adj, perm, perm_inv);

        bandwidth = AdjacencyMatrix.adj_perm_bandwidth(node_num, adj_num, adj_row, adj,
            perm, perm_inv);

        Console.WriteLine("");
        Console.WriteLine("  ADJ (permuted) bandwidth = " + bandwidth + "");
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests GENRCM
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] adj;
        int adj_num;
        int[] adj_row;
        int bandwidth;
        int hole_num = 0;
        int i;
        int j;
        int node;
        int node_num = 0;
        double[] node_xy;
        int[] perm;
        int[] perm_inv;
        int seed;
        int test;
        int triangle_num = 0;
        int[] triangle_neighbor;
        int[] triangle_node;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  GENRCM generates the Reverse Cuthill McKee ordering.");
        Console.WriteLine("");
        Console.WriteLine("  Do the test twice.  On the second test, randomly");
        Console.WriteLine("  permute the initial nodes.");

        Order3_Example.triangulation_order3_example2_size(ref node_num, ref triangle_num, ref hole_num);

        for (test = 1; test <= 2; test++)
        {
            node_xy = new double[2 * node_num];
            triangle_node = new int[3 * triangle_num];
            triangle_neighbor = new int[3 * triangle_num];

            Order3_Example.triangulation_order3_example2(node_num, triangle_num, ref node_xy,
                ref triangle_node, ref triangle_neighbor);
            switch (test)
            {
                //
                //  Randomly permute the nodes.
                //
                case 2:
                {
                    seed = 123456789;

                    perm = typeMethods.perm_uniform(node_num, 0, ref seed);

                    typeMethods.i4vec_print(node_num, perm, "  The random permutation:");

                    for (i = 0; i < 3; i++)
                    {
                        for (j = 0; j < triangle_num; j++)
                        {
                            node = triangle_node[i + j * 3];
                            triangle_node[i + j * 3] = perm[node - 1];
                        }
                    }

                    break;
                }
            }

            typeMethods.i4mat_transpose_print(3, triangle_num, triangle_node,
                "  TRIANGLE_NODE:");

            adj_row = new int[node_num + 1];

            adj_num = Adjacency.triangulation_order3_adj_count(node_num, triangle_num,
                triangle_node, triangle_neighbor, adj_row);

            adj = Adjacency.triangulation_order3_adj_set(node_num, triangle_num, triangle_node,
                triangle_neighbor, adj_num, adj_row);

            AdjacencyMatrix.adj_print(node_num, adj_num, adj_row, adj, "  ADJ array:");

            bandwidth = AdjacencyMatrix.adj_bandwidth(node_num, adj_num, adj_row, adj);

            Console.WriteLine("");
            Console.WriteLine("  ADJ bandwidth = " + bandwidth + "");

            perm = new int[node_num];

            perm = Burkardt.Graph.GenRCM.genrcm(node_num, adj_num, adj_row, adj);

            typeMethods.i4vec_print(node_num, perm, "  The RCM permutation:");

            perm_inv = new int[node_num];

            typeMethods.perm_inverse3(node_num, perm, ref perm_inv);

            bandwidth = AdjacencyMatrix.adj_perm_bandwidth(node_num, adj_num, adj_row, adj,
                perm, perm_inv);

            Console.WriteLine("");
            Console.WriteLine("  Permuted ADJ bandwidth = " + bandwidth + "");
        }
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests GENRCM
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] adj;
        int adj_num;
        int[] adj_row;
        int bandwidth;
        int hole_num = 0;
        int i;
        int j;
        int node;
        int node_num = 0;
        double[] node_xy;
        int[] perm;
        int[] perm_inv;
        int seed;
        int triangle_num = 0;
        int[] triangle_neighbor;
        int[] triangle_node;
        int triangle_order = 6;

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  GENRCM generates the Reverse Cuthill McKee ordering.");

        Order6_Example.triangulation_order6_example2_size(ref node_num, ref triangle_num, ref hole_num);

        node_xy = new double[2 * node_num];
        triangle_node = new int[triangle_order * triangle_num];
        triangle_neighbor = new int[3 * triangle_num];

        Order6_Example.triangulation_order6_example2(node_num, triangle_num, ref node_xy,
            ref triangle_node, ref triangle_neighbor);
        //
        //  Randomly permute the nodes.
        //
        seed = 123456789;

        perm = typeMethods.perm_uniform(node_num, 0, ref seed);

        typeMethods.i4vec_print(node_num, perm, "  The random permutation:");

        for (i = 0; i < triangle_order; i++)
        {
            for (j = 0; j < triangle_num; j++)
            {
                node = triangle_node[i + j * triangle_order];
                triangle_node[i + j * triangle_order] = perm[node - 1];
            }
        }

        typeMethods.i4mat_transpose_print(triangle_order, triangle_num, triangle_node,
            "  Permuted TRIANGLE_NODE");

        adj_row = new int[node_num + 1];

        adj_num = Adjacency.triangulation_order6_adj_count(node_num, triangle_num,
            triangle_node, ref triangle_neighbor, ref adj_row);

        adj = Adjacency.triangulation_order6_adj_set(node_num, triangle_num, triangle_node,
            triangle_neighbor, adj_num, adj_row);

        AdjacencyMatrix.adj_print(node_num, adj_num, adj_row, adj, "  ADJ array:");

        bandwidth = AdjacencyMatrix.adj_bandwidth(node_num, adj_num, adj_row, adj);

        Console.WriteLine("");
        Console.WriteLine("  ADJ bandwidth = " + bandwidth + "");

        perm = new int[node_num];

        perm = Burkardt.Graph.GenRCM.genrcm(node_num, adj_num, adj_row, adj);

        typeMethods.i4vec_print(node_num, perm, "  The RCM permutation:");

        perm_inv = new int[node_num];

        typeMethods.perm_inverse3(node_num, perm, ref perm_inv);

        bandwidth = AdjacencyMatrix.adj_perm_bandwidth(node_num, adj_num, adj_row, adj,
            perm, perm_inv);

        Console.WriteLine("");
        Console.WriteLine("  Permuted ADJ bandwidth = " + bandwidth + "");
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests GRAPH_01_ADJ and GRAPH_01_SIZE;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] adj;
        int adj_num = 0;
        int[] adj_row;
        int node_num = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  GRAPH_01_SIZE returns the sizes for graph 1.");
        Console.WriteLine("  GRAPH_01_ADJ returns the adjacency for graph 1.");
        Console.WriteLine("  ADJ_PRINT prints the adjacency information.");

        Burkardt.Graph.Adjacency.graph_01_size(ref node_num, ref adj_num);

        adj_row = new int[node_num + 1];
        adj = new int[adj_num];

        Burkardt.Graph.Adjacency.graph_01_adj(node_num, adj_num, ref adj_row, ref adj);

        AdjacencyMatrix.adj_print(node_num, adj_num, adj_row, adj,
            "  Adjacency for GRAPH_01:");

        AdjacencyMatrix.adj_show(node_num, adj_num, adj_row, adj);
    }

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tests LEVEL_SET;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] adj;
        int adj_num = 0;
        int[] adj_row;
        int i;
        int j;
        int[] level;
        int level_num = 0;
        int[] level_row;
        int[] mask;
        int node_num = 0;
        int root;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  LEVEL_SET computes the level sets of a graph,");
        Console.WriteLine("  given a root node (which defines level 1).");

        Burkardt.Graph.Adjacency.graph_01_size(ref node_num, ref adj_num);

        adj_row = new int[node_num + 1];
        adj = new int[adj_num];

        Burkardt.Graph.Adjacency.graph_01_adj(node_num, adj_num, ref adj_row, ref adj);

        AdjacencyMatrix.adj_print(node_num, adj_num, adj_row, adj, "  Adjacency matrix:");

        AdjacencyMatrix.adj_show(node_num, adj_num, adj_row, adj);
        //
        //  Choose different roots.
        //
        level = new int[node_num];
        level_row = new int[node_num + 1];
        mask = new int[node_num];

        for (i = 1; i <= 3; i++)
        {
            root = UniformRNG.i4_uniform(1, node_num, ref seed);

            for (j = 0; j < node_num; j++)
            {
                mask[j] = 1;
            }

            Burkardt.Graph.GenRCM.level_set(root, adj_num, adj_row, adj, ref mask, ref level_num,
                ref level_row, ref level, node_num);

            Burkardt.Graph.GenRCM.level_set_print(node_num, level_num, level_row, level);
        }
    }

    private static void test07()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 tests ROOT_FIND;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] adj;
        int adj_num = 0;
        int[] adj_row;
        int i;
        int j;
        int[] level;
        int level_num = 0;
        int[] level_row;
        int[] mask;
        int node_num = 0;
        int root;

        Console.WriteLine("");
        Console.WriteLine("TEST07");
        Console.WriteLine("  ROOT_FIND is given a node in the graph,");
        Console.WriteLine("  and returns a better node to use as a starting");
        Console.WriteLine("  point for reordering.");

        Burkardt.Graph.Adjacency.graph_01_size(ref node_num, ref adj_num);

        adj_row = new int[node_num + 1];
        adj = new int[adj_num];

        Burkardt.Graph.Adjacency.graph_01_adj(node_num, adj_num, ref adj_row, ref adj);

        AdjacencyMatrix.adj_print(node_num, adj_num, adj_row, adj, "  Adjacency matrix:");

        AdjacencyMatrix.adj_show(node_num, adj_num, adj_row, adj);

        level = new int[node_num];
        level_row = new int[node_num + 1];
        mask = new int[node_num];

        for (i = 1; i <= node_num; i++)
        {
            root = i;

            Console.WriteLine("");
            Console.WriteLine("  Starting root =    " + root + "");

            for (j = 0; j < node_num; j++)
            {
                mask[j] = 1;
            }

            Burkardt.Graph.GenRCM.root_find(ref root, adj_num, adj_row, adj, mask, ref level_num,
                ref level_row, ref level, node_num);

            Console.WriteLine("  Suggested root =   " + root + "");
            Console.WriteLine("  Number of levels = " + level_num + "");
        }
    }

    private static void test08()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 tests TRIANGULATION_ORDER3_ADJ_COUNT
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] adj_row;
        int hole_num = 0;
        int node_num = 0;
        double[] node_xy;
        int triangle_num = 0;
        int[] triangle_neighbor;
        int[] triangle_node;

        Console.WriteLine("");
        Console.WriteLine("TEST08");
        Console.WriteLine("  TRIANGULATION_ORDER3_ADJ_COUNT counts the (lower)");
        Console.WriteLine("  adjacencies defined by a triangulation.");

        Order3_Example.triangulation_order3_example2_size(ref node_num, ref triangle_num, ref hole_num);

        node_xy = new double[2 * node_num];
        triangle_node = new int[3 * triangle_num];
        triangle_neighbor = new int[3 * triangle_num];

        Order3_Example.triangulation_order3_example2(node_num, triangle_num, ref node_xy,
            ref triangle_node, ref triangle_neighbor);

        typeMethods.i4mat_transpose_print(3, triangle_num, triangle_node,
            "  TRIANGLE_NODE");

        adj_row = new int[node_num + 1];

        Adjacency.triangulation_order3_adj_count(node_num, triangle_num,
            triangle_node, triangle_neighbor, adj_row);

        typeMethods.i4vec_print(node_num + 1, adj_row, "  ADJ_ROW");
    }

    private static void test09()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST09 tests TRIANGULATION_ORDER3_ADJ_SET
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] adj;
        int adj_num;
        int[] adj_row;
        int bandwidth;
        int hole_num = 0;
        int node_num = 0;
        double[] node_xy;
        int triangle_num = 0;
        int[] triangle_neighbor;
        int[] triangle_node;

        Console.WriteLine("");
        Console.WriteLine("TEST09");
        Console.WriteLine("  TRIANGULATION_ORDER3_ADJ_SET sets the (lower)");
        Console.WriteLine("  adjacencies defined by a triangulation.");

        Order3_Example.triangulation_order3_example2_size(ref node_num, ref triangle_num, ref hole_num);

        node_xy = new double[2 * node_num];
        triangle_node = new int[3 * triangle_num];
        triangle_neighbor = new int[3 * triangle_num];

        Order3_Example.triangulation_order3_example2(node_num, triangle_num, ref node_xy,
            ref triangle_node, ref triangle_neighbor);

        typeMethods.i4mat_transpose_print(3, triangle_num, triangle_node,
            "  TRIANGLE_NODE");

        adj_row = new int[node_num + 1];

        adj_num = Adjacency.triangulation_order3_adj_count(node_num, triangle_num,
            triangle_node, triangle_neighbor, adj_row);

        adj = Adjacency.triangulation_order3_adj_set(node_num, triangle_num, triangle_node,
            triangle_neighbor, adj_num, adj_row);

        AdjacencyMatrix.adj_print(node_num, adj_num, adj_row, adj, "  ADJ array:");

        bandwidth = AdjacencyMatrix.adj_bandwidth(node_num, adj_num, adj_row, adj);

        Console.WriteLine("");
        Console.WriteLine("  ADJ bandwidth = " + bandwidth + "");
    }

    private static void test10()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST10 tests TRIANGULATION_NEIGHBOR_TRIANGLES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 September 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int TRIANGLE_ORDER = 3;
        int TRIANGLE_NUM = 16;

        int[] triangle_node =
        {
            3, 4, 1,
            3, 1, 2,
            3, 2, 8,
            2, 1, 5,
            8, 2, 13,
            8, 13, 9,
            3, 8, 9,
            13, 2, 5,
            9, 13, 7,
            7, 13, 5,
            6, 7, 5,
            9, 7, 6,
            10, 9, 6,
            6, 5, 12,
            11, 6, 12,
            10, 6, 11
        };
        int[] triangle_neighbor;

        Console.WriteLine("");
        Console.WriteLine("TEST10");
        Console.WriteLine("  For a triangulation of a set of nodes,");
        Console.WriteLine("  TRIANGULATION_NEIGHBOR_TRIANGLES determines the");
        Console.WriteLine("    adjacency relationships between triangles.");

        typeMethods.i4mat_transpose_print(TRIANGLE_ORDER, TRIANGLE_NUM, triangle_node,
            "  Triangles:");

        triangle_neighbor = NeighborElements.triangulation_neighbor_triangles(TRIANGLE_ORDER,
            TRIANGLE_NUM, triangle_node);

        typeMethods.i4mat_transpose_print(3, TRIANGLE_NUM, triangle_neighbor,
            "  Triangle neighbors:");
    }

    private static void test11()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST11 tests TRIANGULATION_ORDER6_ADJ_COUNT
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] adj_row;
        int hole_num = 0;
        int node_num = 0;
        double[] node_xy;
        int triangle_num = 0;
        int[] triangle_neighbor;
        int[] triangle_node;

        Console.WriteLine("");
        Console.WriteLine("TEST11");
        Console.WriteLine("  TRIANGULATION_ORDER6_ADJ_COUNT counts the (lower)");
        Console.WriteLine("  adjacencies defined by a triangulation.");

        Order6_Example.triangulation_order6_example2_size(ref node_num, ref triangle_num, ref hole_num);

        node_xy = new double[2 * node_num];
        triangle_node = new int[6 * triangle_num];
        triangle_neighbor = new int[3 * triangle_num];

        Order6_Example.triangulation_order6_example2(node_num, triangle_num, ref node_xy,
            ref triangle_node, ref triangle_neighbor);

        adj_row = new int[node_num + 1];

        Adjacency.triangulation_order6_adj_count(node_num, triangle_num,
            triangle_node, ref triangle_neighbor, ref adj_row);

        typeMethods.i4vec_print(node_num + 1, adj_row, "  ADJ_ROW");
    }

    private static void test12()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST12 tests TRIANGULATION_ORDER6_ADJ_SET
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] adj;
        int adj_num;
        int[] adj_row;
        int bandwidth;
        int hole_num = 0;
        int node_num = 0;
        double[] node_xy;
        int triangle_num = 0;
        int[] triangle_neighbor;
        int[] triangle_node;

        Console.WriteLine("");
        Console.WriteLine("TEST12");
        Console.WriteLine("  TRIANGULATION_ORDER6_ADJ_SET sets the (lower)");
        Console.WriteLine("  adjacencies defined by a triangulation.");

        Order6_Example.triangulation_order6_example2_size(ref node_num, ref triangle_num, ref hole_num);

        node_xy = new double[2 * node_num];
        triangle_node = new int[6 * triangle_num];
        triangle_neighbor = new int[3 * triangle_num];

        Order6_Example.triangulation_order6_example2(node_num, triangle_num, ref node_xy,
            ref triangle_node, ref triangle_neighbor);

        typeMethods.i4mat_transpose_print(6, triangle_num, triangle_node,
            "  TRIANGLE_NODE");

        adj_row = new int[node_num + 1];

        adj_num = Adjacency.triangulation_order6_adj_count(node_num, triangle_num, triangle_node,
            ref triangle_neighbor, ref adj_row);

        adj = new int [adj_num];

        adj = Adjacency.triangulation_order6_adj_set(node_num, triangle_num, triangle_node,
            triangle_neighbor, adj_num, adj_row);

        AdjacencyMatrix.adj_print(node_num, adj_num, adj_row, adj, "  ADJ array:");

        bandwidth = AdjacencyMatrix.adj_bandwidth(node_num, adj_num, adj_row, adj);

        Console.WriteLine("");
        Console.WriteLine("  ADJ bandwidth = " + bandwidth + "");
    }
}