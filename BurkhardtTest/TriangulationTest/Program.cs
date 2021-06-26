using System;
using Burkardt;
using Burkardt.TriangulationNS;
using Burkardt.Types;

namespace TriangulationTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for TRIANGULATION_TEST.
            //
            //  Discussion:
            //
            //    TRIANGULATION_TEST tests the TRIANGULATION library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    03 March 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGULATION_TEST");
            Console.WriteLine("  Test the TRIANGULATION library.");

            test01();
            test02();
            test025();
            test026();
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
            test125();
            test127();
            test13();
            test14();
            test15();
            test16();
            test17();
            test18();
            test19();

            test20();
            test21();
            test213();
            test215();
            test217();
            test219();
            test22();
            test23();
            test24();
            test25();
            test26();
            test265();
            test27();

            test31();
            test32();
            test33();

            Console.WriteLine("");
            Console.WriteLine("TRIANGULATION_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 tests ALPHA_MEASURE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    21 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double alpha_ave = 0;
            double alpha_area = 0;
            double alpha_min = 0;
            int hole_num = 0;
            int node_num = 0;
            double[] node_xy;
            int triangle_num = 0;
            int[] triangle_node;
            int[] triangle_neighbor;
            int triangle_order = 3;

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  ALPHA_MEASURE returns the ALPHA measure of");
            Console.WriteLine("  quality of a triangulation.");
            //
            //  Get the sizes.
            //
            TriangulationSampleData.triangulation_order3_example1_size(ref node_num, ref triangle_num, ref hole_num);
            //
            //  Allocate space.
            //
            node_xy = new double[2 * node_num];
            triangle_node = new int [triangle_order * triangle_num];
            triangle_neighbor = new int[3 * triangle_num];
            //
            //  Get the triangulation data.
            //
            TriangulationSampleData.triangulation_order3_example1(node_num, triangle_num, ref node_xy,
                ref triangle_node, ref triangle_neighbor);
            //
            //  Compute the triangulation quality.
            //
            Measure.alpha_measure(node_num, node_xy, triangle_order, triangle_num,
                triangle_node, ref alpha_min, ref alpha_ave, ref alpha_area);

            Console.WriteLine("");
            Console.WriteLine("  ALPHA_MIN  = " + alpha_min + "");
            Console.WriteLine("  ALPHA_AVE  = " + alpha_ave + "");
            Console.WriteLine("  ALPHA_AREA = " + alpha_area + "");
        }

        static void test02()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 tests AREA_MEASURE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    21 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double area_ave = 0;
            double area_max = 0;
            double area_min = 0;
            double area_ratio = 0;
            double area_std = 0;
            int hole_num = 0;
            int node_num = 0;
            double[] node_xy;
            int triangle_num = 0;
            int[] triangle_node;
            int[] triangle_neighbor;
            int triangle_order = 3;

            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  AREA_MEASURE returns the AREA measure of");
            Console.WriteLine("  quality of a triangulation.");
            //
            //  Get the sizes.
            //
            TriangulationSampleData.triangulation_order3_example1_size(ref node_num, ref triangle_num, ref hole_num);
            //
            //  Allocate space.
            //
            node_xy = new double[2 * node_num];
            triangle_node = new int [triangle_order * triangle_num];
            triangle_neighbor = new int[3 * triangle_num];
            //
            //  Get the triangulation data.
            //
            TriangulationSampleData.triangulation_order3_example1(node_num, triangle_num, ref node_xy,
                ref triangle_node, ref triangle_neighbor);
            //
            //  Compute the triangulation quality.
            //
            Measure.area_measure(node_num, node_xy, triangle_order, triangle_num,
                triangle_node, ref area_min, ref area_max, ref area_ratio, ref area_ave, ref area_std);

            Console.WriteLine("");
            Console.WriteLine("  AREA_MIN   = " + area_min + "");
            Console.WriteLine("  AREA_MAX   = " + area_max + "");
            Console.WriteLine("  AREA_RATIO = " + area_ratio + "");
            Console.WriteLine("  AREA_AVE   = " + area_ave + "");
            Console.WriteLine("  AREA_STD   = " + area_std + "");
        }

        static void test025()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST025 tests DELAUNAY_SWAP_TEST.
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
        {
            int node_num = 4;
            int triangle_num = 2;
            int triangle_order = 3;

            double alpha_area = 0;
            double alpha_ave = 0;
            double alpha_min_swapped = 0;
            double alpha_min_unswapped = 0;
            double[] node_xy = new double[2 * 4];
            int seed = 123456789;
            bool swap;
            int test;
            int test_num = 10;
            int[] triangle_node = new int[3 * 2];

            Console.WriteLine("");
            Console.WriteLine("TEST025");
            Console.WriteLine("  DELAUNAY_SWAP_TEST determines whether two triangles");
            Console.WriteLine("  with a common edge need to \"swap\" diagonals.");
            Console.WriteLine("  If swapping is indicated, then ALPHA_MIN should increase.");
            Console.WriteLine("");
            Console.WriteLine("  Swap   ALPHA_MIN   ALPHA_MIN");
            Console.WriteLine("         Unswapped   Swapped");
            Console.WriteLine("");

            for (test = 1; test <= test_num; test++)
            {
                //
                //  Generate a random quadrilateral (1,2,3,4).
                //
                TriangulationSampleData.quad_convex_random(ref seed, ref node_xy);
                //
                //  Does it need swapping?
                //
                swap = Delauney.delaunay_swap_test(node_xy);
                //
                //  Compute ALPHA_MIN unswapped.
                //
                triangle_node[0 + 0 * 3] = 1;
                triangle_node[1 + 0 * 3] = 2;
                triangle_node[2 + 0 * 3] = 3;
                triangle_node[0 + 1 * 3] = 1;
                triangle_node[1 + 1 * 3] = 3;
                triangle_node[2 + 1 * 3] = 4;

                Measure.alpha_measure(node_num, node_xy, triangle_order, triangle_num,
                    triangle_node, ref alpha_min_unswapped, ref alpha_ave, ref alpha_area);
                //
                //  Compute ALPHA_MIN swapped.
                //
                triangle_node[0 + 0 * 3] = 1;
                triangle_node[1 + 0 * 3] = 2;
                triangle_node[2 + 0 * 3] = 4;
                triangle_node[0 + 1 * 3] = 2;
                triangle_node[1 + 1 * 3] = 3;
                triangle_node[2 + 1 * 3] = 4;

                Measure.alpha_measure(node_num, node_xy, triangle_order, triangle_num,
                    triangle_node, ref alpha_min_swapped, ref alpha_ave, ref alpha_area);

                if (false)
                {
                    typeMethods.r8mat_transpose_print(2, node_num, node_xy, "  Quadrilateral");
                }

                Console.WriteLine("     " + swap
                    + "  " + alpha_min_unswapped.ToString().PadLeft(10)
                    + "  " + alpha_min_swapped.ToString().PadLeft(10) + "");
            }
        }

        static void test026()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST026 tests DIAEDG.
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
        {
            int node_num = 4;
            int triangle_num = 2;
            int triangle_order = 3;

            double alpha_area = 0;
            double alpha_ave = 0;
            double alpha_min_swapped = 0;
            double alpha_min_unswapped = 0;
            double[] node_xy = new double[2 * 4];
            int seed = 123456789;
            bool swap;
            int test;
            int test_num = 10;
            int[] triangle_node = new int[3 * 2];
            int value;

            Console.WriteLine("");
            Console.WriteLine("TEST026");
            Console.WriteLine("  DIAEDG determines whether two triangles");
            Console.WriteLine("  with a common edge need to \"swap\" diagonals.");
            Console.WriteLine("  If swapping is indicated, then ALPHA_MIN should increase.");
            Console.WriteLine("");
            Console.WriteLine("  Swap   ALPHA_MIN   ALPHA_MIN");
            Console.WriteLine("         Unswapped   Swapped");
            Console.WriteLine("");

            for (test = 1; test <= test_num; test++)
            {
                //
                //  Generate a random quadrilateral (1,2,3,4).
                //
                TriangulationSampleData.quad_convex_random(ref seed, ref node_xy);
                //
                //  Does it need swapping?
                //
                value = typeMethods.diaedg(
                    node_xy[0 + 0 * 2], node_xy[1 + 0 * 2],
                    node_xy[0 + 1 * 2], node_xy[1 + 1 * 2],
                    node_xy[0 + 2 * 2], node_xy[1 + 2 * 2],
                    node_xy[0 + 3 * 2], node_xy[1 + 3 * 2]);

                if (value == 1)
                {
                    swap = false;
                }
                else
                {
                    swap = true;
                }

                //
                //  Compute ALPHA_MIN unswapped.
                //
                triangle_node[0 + 0 * 3] = 1;
                triangle_node[1 + 0 * 3] = 2;
                triangle_node[2 + 0 * 3] = 3;
                triangle_node[0 + 1 * 3] = 1;
                triangle_node[1 + 1 * 3] = 3;
                triangle_node[2 + 1 * 3] = 4;

                Measure.alpha_measure(node_num, node_xy, triangle_order, triangle_num,
                    triangle_node, ref alpha_min_unswapped, ref alpha_ave, ref alpha_area);
                //
                //  Compute ALPHA_MIN swapped.
                //
                triangle_node[0 + 0 * 3] = 1;
                triangle_node[1 + 0 * 3] = 2;
                triangle_node[2 + 0 * 3] = 4;
                triangle_node[0 + 1 * 3] = 2;
                triangle_node[1 + 1 * 3] = 3;
                triangle_node[2 + 1 * 3] = 4;

                Measure.alpha_measure(node_num, node_xy, triangle_order, triangle_num,
                    triangle_node, ref alpha_min_swapped, ref alpha_ave, ref alpha_area);

                if (false)
                {
                    typeMethods.r8mat_transpose_print(2, node_num, node_xy, "  Quadrilateral");
                }

                Console.WriteLine("     " + swap
                    + "  " + alpha_min_unswapped.ToString().PadLeft(10)
                    + "  " + alpha_min_swapped.ToString().PadLeft(10) + "");
            }
        }

        static void test03()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST03 tests NODE_MERGE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 August 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;
            int NODE_NUM = 15;
            int TEST_NUM = 4;

            int node;
            int[] node_rep = new int[NODE_NUM];
            double[] node_xy =  {
                0.0, 0.0,
                1.0, 0.0,
                3.0, 0.0,
                4.0, 0.0,
                1.0, 1.0,
                4.0, 1.0,
                2.0, 2.0,
                3.0, 3.0,
                2.0, 3.5,
                0.5, 4.0,
                1.0, 4.0,
                1.5, 4.0,
                4.0, 4.0,
                1.0, 4.5,
                1.0, 4.5
            }
            ;
            int rep;
            int rep_num;
            int test;
            double tolerance;
            double[] tolerance_test =  {
                0.01, 0.75, 1.2, 1.5
            }
            ;

            Console.WriteLine("");
            Console.WriteLine("TEST03");
            Console.WriteLine("  NODE_MERGE identifies groups of nodes");
            Console.WriteLine("  that can be merged, with a given tolerance.");

            typeMethods.r8mat_transpose_print(DIM_NUM, NODE_NUM, node_xy, "  Node coordinates:");

            for (test = 0; test < TEST_NUM; test++)
            {
                tolerance = tolerance_test[test];

                Node.node_merge(DIM_NUM, NODE_NUM, node_xy, tolerance, ref node_rep);

                Console.WriteLine("");
                Console.WriteLine("  TOLERANCE = " + tolerance + "");
                Console.WriteLine("");
                Console.WriteLine("      Node  Representatives:");
                Console.WriteLine("");

                for (node = 0; node < NODE_NUM; node++)
                {
                    Console.WriteLine("  " + node.ToString().PadLeft(8)
                        + "  " + node_rep[node].ToString().PadLeft(8) + "");
                }

                //
                //  Make a list of the node representatives.
                //
                Console.WriteLine("");
                Console.WriteLine("      Rep   Coordinates:");
                Console.WriteLine("");

                typeMethods.i4vec_sort_heap_a(NODE_NUM, ref node_rep);

                rep_num = 0;

                for (node = 0; node < NODE_NUM; node++)
                {
                    if (1 <= node)
                    {
                        if (node_rep[node - 1] == node_rep[node])
                        {
                            continue;
                        }
                    }

                    rep = node_rep[node];

                    Console.WriteLine("  " + rep_num.ToString().PadLeft(8)
                        + "  " + node_xy[0 + rep * 2].ToString().PadLeft(12)
                        + "  " + node_xy[1 + rep * 2].ToString().PadLeft(12) + "");

                    rep_num = rep_num + 1;
                }
            }
        }

        static void test04()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST04 tests NS_ADJ_COL_SET, NS_ADJ_COUNT and NS_ADJ_ROW_SET.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    26 September 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int NODE_NUM = 15;
            int TRIANGLE_NUM = 4;
            int TRIANGLE_ORDER = 6;
            int VARIABLE_NUM = 36;

            int[] adj_col = new int[VARIABLE_NUM + 1];
            int adj_num;
            int[] adj_row;
            string file_name = "ns_triangulation.eps";
            int node_show;
            int[] node_p_variable =  {
                3, -1, 8, -1, 13,
                -1, -1, -1, -1,
                24, -1, 29,
                -1, -1,
                36
            }
            ;
            int[] node_u_variable =  {
                1, 4, 6, 9, 11,
                14, 16, 18, 20,
                22, 25, 27,
                30, 32,
                34
            }
            ;
            int[] node_v_variable =  {
                2, 5, 7, 10, 12,
                15, 17, 19, 21,
                23, 26, 28,
                31, 33,
                35
            }
            ;
            double[] node_xy =  {
                0.0, 0.0,
                0.0, 1.0,
                0.0, 2.0,
                0.0, 3.0,
                0.0, 4.0,
                1.0, 0.0,
                1.0, 1.0,
                1.0, 2.0,
                1.0, 3.0,
                2.0, 0.0,
                2.0, 1.0,
                2.0, 2.0,
                3.0, 0.0,
                3.0, 1.0,
                4.0, 0.0
            }
            ;
            int num;
            int r;
            int rhi;
            int rlo;
            int[] triangle_neighbor =  {
                -1, 2, -1,
                3, 1, 4,
                2, -1, -1,
                -1, -1, 2
            }
            ;
            int[] triangle_node =  {
                1, 10, 3, 6, 7, 2,
                12, 3, 10, 8, 7, 11,
                3, 12, 5, 8, 9, 4,
                10, 15, 12, 13, 14, 11
            }
            ;
            int triangle_show;
            int variable;

            Console.WriteLine("");
            Console.WriteLine("TEST04");
            Console.WriteLine("  For an order 3/order 6 Taylor Hood triangulation");
            Console.WriteLine("  for Navier Stokes velocity and pressure,");
            Console.WriteLine("  NS_ADJ_COUNT counts variable adjacencies");
            Console.WriteLine("  and sets up the sparse compressed column");
            Console.WriteLine("  column pointer array.");
            Console.WriteLine("  NS_ADJ_COL_SET sets up the sparse compressed column");
            Console.WriteLine("  COL vector.");
            Console.WriteLine("  NS_ADJ_ROW_SET sets up the sparse compressed column");
            Console.WriteLine("  ROW vector.");
            //
            //  Plot the example.
            //
            node_show = 2;
            triangle_show = 2;

            Plot.triangulation_order6_plot(file_name, NODE_NUM, node_xy,
                TRIANGLE_NUM, triangle_node, node_show, triangle_show);
            //
            //  Get the count of the variable adjacencies.
            //  We don't really need to make this call, since the next
            //  call does the calculation as part of getting ADJ_COL.
            //
            Console.WriteLine("");
            Console.WriteLine("  Number of variables is " + VARIABLE_NUM + "");

            adj_num = NavierStokes.ns_adj_count(NODE_NUM, TRIANGLE_NUM, VARIABLE_NUM, triangle_node,
                triangle_neighbor, node_u_variable, node_v_variable, node_p_variable);

            Console.WriteLine("");
            Console.WriteLine("  As computed by NS_ADJ_COUNT,");
            Console.WriteLine("  Number of variable adjacency entries is " + adj_num + "");
            //
            //  Get the count of the variable adjacencies and the COL vector.
            //
            adj_num = NavierStokes.ns_adj_col_set(NODE_NUM, TRIANGLE_NUM, VARIABLE_NUM, triangle_node,
                triangle_neighbor, node_u_variable, node_v_variable, node_p_variable,
                ref adj_col);

            Console.WriteLine("");
            Console.WriteLine("  As computed by NS_ADJ_COL_SET,");
            Console.WriteLine("  Number of variable adjacency entries is " + adj_num + "");

            Console.WriteLine("");
            Console.WriteLine("  Variable adjacency pointers:");
            Console.WriteLine("");
            Console.WriteLine("  Variable     First      Last    Number");
            Console.WriteLine("");

            for (variable = 0; variable < VARIABLE_NUM; variable++)
            {
                num = adj_col[variable + 1] - adj_col[variable];

                Console.WriteLine("  " + (variable + 1).ToString().PadLeft(8)
                    + "  " + adj_col[variable].ToString().PadLeft(8)
                    + "  " + (adj_col[variable + 1] - 1).ToString().PadLeft(8)
                    + "  " + num.ToString().PadLeft(8) + "");
            }

            //
            //  Get the variable adjacencies.
            //
            adj_row = new int[adj_num];

            NavierStokes.ns_adj_row_set(NODE_NUM, TRIANGLE_NUM, VARIABLE_NUM, triangle_node,
                triangle_neighbor, node_u_variable, node_v_variable, node_p_variable,
                adj_num, adj_col, adj_row);
            //
            //  This is a huge array.  We only print out the beginning and end.
            //
            Console.WriteLine("");
            Console.WriteLine("  Variable adjacency row entries:");
            Console.WriteLine("  (Partial printout only)");
            Console.WriteLine("");
            Console.WriteLine("     Entry     Row       Col");
            Console.WriteLine("");

            for (variable = 0; variable < VARIABLE_NUM; variable++)
            {
                rlo = adj_col[variable] - 1;
                rhi = adj_col[variable + 1] - 2;

                if (variable <= 2 || VARIABLE_NUM - 4 <= variable)
                {
                    Console.WriteLine("");

                    for (r = rlo; r <= rhi; r++)
                    {
                        Console.WriteLine("  " + (r + 1).ToString().PadLeft(8)
                            + "  " + adj_row[r].ToString().PadLeft(8)
                            + "  " + (variable + 1).ToString().PadLeft(8) + "");
                    }
                }

                if (variable == 2)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  (SKIPPING MANY MANY ENTRIES...)");
                    Console.WriteLine("");
                }

            }
        }

        static void test05()

            //****************************************************************************80
            //
            //  Purpose:  
            //
            //    TEST04 tests POINTS_DELAUNAY_NAIVE_2D.
            //
            //  Diagram:
            //
            //    !....3&11....
            //    !............
            //    !............
            //    X..9.........
            //    !.....5......
            //    !...........6
            //    !.4.2...10...
            //    !.....8...12.
            //    V............
            //    !..7.........
            //    !......1.....
            //    !............
            //    !............
            //    !----V----X--
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 August 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int NODE_NUM = 12;
            int DIM_NUM = 2;

            int triangle_num = 0;
            int[] triangle_node;
            double[] node_xy =  {
                7.0, 3.0,
                4.0, 7.0,
                5.0, 13.0,
                2.0, 7.0,
                6.0, 9.0,
                12.0, 10.0,
                3.0, 4.0,
                6.0, 6.0,
                3.0, 10.0,
                8.0, 7.0,
                5.0, 13.0,
                10.0, 6.0
            }
            ;

            Console.WriteLine("");
            Console.WriteLine("TEST05");
            Console.WriteLine("  POINTS_DELAUNAY_NAIVE_2D computes the Delaunay");
            Console.WriteLine("  triangulation of a set of nodes.");

            typeMethods.r8mat_transpose_print(DIM_NUM, NODE_NUM, node_xy, "  The nodes:");

            triangle_node = Points.points_delaunay_naive_2d(NODE_NUM, node_xy, ref triangle_num);

            Console.WriteLine("");
            Console.WriteLine("  Number of triangles is TRIANGLE_NUM = " + triangle_num + "");

            typeMethods.i4mat_transpose_print(3, triangle_num, triangle_node,
                "  The Delaunay triangles:");

        }

        static void test06()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST06 tests POINTS_HULL_2D.
            //
            //  Diagram:
            //
            //    !....3.......
            //    !............
            //    !..9.........
            //    !.....5......
            //    !...........6
            //    !.4.2...10...
            //    !.....8......
            //    !.........12.
            //    !..7.........
            //    !......1.....
            //    !............
            //    !............
            //    !-----------
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 August 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int NODE_NUM = 12;
            int DIM_NUM = 2;

            int i;
            int j;
            int[] ival = new int[NODE_NUM];
            int nval = 0;
            double[] node_xy =  {
                7.0, 3.0,
                4.0, 7.0,
                5.0, 13.0,
                2.0, 7.0,
                6.0, 9.0,
                12.0, 8.0,
                3.0, 4.0,
                6.0, 6.0,
                3.0, 10.0,
                8.0, 7.0,
                5.0, 13.0,
                10.0, 6.0
            }
            ;

            Console.WriteLine("");
            Console.WriteLine("TEST06");
            Console.WriteLine("  POINTS_HULL_2D computes the convex hull");
            Console.WriteLine("  of a set of nodes.");

            typeMethods.r8mat_transpose_print(DIM_NUM, NODE_NUM, node_xy, "  The nodes:");

            Points.points_hull_2d(NODE_NUM, node_xy, ref nval, ref ival);

            Console.WriteLine("");
            Console.WriteLine("  The convex hull is formed by connecting:");
            Console.WriteLine("");
            for (j = 0; j < nval; j++)
            {
                string cout = "  "
                     + j.ToString().PadLeft(3) + "  "
                     + ival[j].ToString().PadLeft(3) + "  ";
                for (i = 0; i < DIM_NUM; i++)
                {
                    cout += "  " + node_xy[i + (ival[j] - 1) * DIM_NUM].ToString().PadLeft(14);
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("  The correct sequence of nodes is:");
            Console.WriteLine("  4, 9, 3, 6, 12, 1, 7, (4).");

        }

        static void test07()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST07 tests Q_MEASURE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    21 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int hole_num = 0;
            int node_num = 0;
            double[] node_xy;
            double q_area = 0;
            double q_ave = 0;
            double q_max = 0;
            double q_min = 0;
            int triangle_num = 0;
            int[] triangle_node;
            int[] triangle_neighbor;
            int triangle_order = 3;

            Console.WriteLine("");
            Console.WriteLine("TEST07");
            Console.WriteLine("  Q_MEASURE returns the Q measure of");
            Console.WriteLine("  quality of a triangulation.");
            //
            //  Get the sizes.
            //
            TriangulationSampleData.triangulation_order3_example1_size(ref node_num, ref triangle_num, ref hole_num);
            //
            //  Allocate space.
            //
            node_xy = new double[2 * node_num];
            triangle_node = new int [triangle_order * triangle_num];
            triangle_neighbor = new int[3 * triangle_num];
            //
            //  Get the triangulation data.
            //
            TriangulationSampleData.triangulation_order3_example1(node_num, triangle_num, ref node_xy,
                ref triangle_node, ref triangle_neighbor);
            //
            //  Compute the triangulation quality.
            //
            Measure.q_measure(node_num, node_xy, triangle_order, triangle_num,
                triangle_node, ref q_min, ref q_max, ref q_ave, ref q_area);

            Console.WriteLine("");
            Console.WriteLine("  Q_MIN  = " + q_min + "");
            Console.WriteLine("  Q_MAX  = " + q_max + "");
            Console.WriteLine("  Q_AVE  = " + q_ave + "");
            Console.WriteLine("  Q_AREA = " + q_area + "");
        }

        static void test08()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST08 tests R8TRIS2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 August 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int NODE_NUM = 9;

            int error;
            double[] node_xy =  {
                0.0, 0.0,
                0.0, 1.0,
                0.2, 0.5,
                0.3, 0.6,
                0.4, 0.5,
                0.6, 0.4,
                0.6, 0.5,
                1.0, 0.0,
                1.0, 1.0
            }
            ;
            int[] triangle_node = new int[2 * NODE_NUM * 3];
            int[] triangle_neighbor = new int[2 * NODE_NUM * 3];
            int triangle_num = 0;

            Console.WriteLine("");
            Console.WriteLine("TEST08");
            Console.WriteLine("  R8TRIS2 computes the Delaunay triangulation of a");
            Console.WriteLine("  pointset in 2D.");
            //
            //  Set up the Delaunay triangulation.
            //
            error = Delauney.r8tris2(NODE_NUM, ref node_xy, ref triangle_num, ref triangle_node,
                ref triangle_neighbor);

            if (error == 0)
            {
                Console.WriteLine("");
                Console.WriteLine("  R8TRIS2 computed the Delaunay triangulation with no");
                Console.WriteLine("  errors detected.");
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("  R8TRIS2 detected an error condition of index " + error + "");
                return;
            }

            Print.triangulation_order3_print(NODE_NUM, triangle_num, node_xy,
                triangle_node, triangle_neighbor);

        }


        static void test09()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST09 tests TRIANGLE_ORDER3_PHYSICAL_TO_REFERENCE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 10;

            int j;
            double[] phy = new double[2 * N];
            double[] ref_ = new double [2 * N];
            double[] ref2 = new double[2 * N];
            int seed;
            double[] t =  {
                1.0, 1.0,
                3.0, 1.0,
                2.0, 5.0
            }
            ;

            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("TEST09");
            Console.WriteLine("  For an order 3 triangle,");
            Console.WriteLine("  TRIANGLE_ORDER3_PHYSICAL_TO_REFERENCE");
            Console.WriteLine("  maps a physical point to a reference point.");
            Console.WriteLine("  TRIANGLE_ORDER3_REFERENCE_TO_PHYSICAL ");
            Console.WriteLine("  maps a reference point to a physical point.");
            Console.WriteLine("");
            Console.WriteLine("   (  XSI    ETA ) ==> ( X      Y  ) ==> ( XSI2    ETA2 )");
            Console.WriteLine("");

            typeMethods.triangle_reference_sample(N, ref seed, ref ref_);

            Triangulation.triangle_order3_reference_to_physical(t, N, ref_, ref phy);

            Triangulation.triangle_order3_physical_to_reference(t, N, phy, ref ref2);

            for (j = 0; j < N; j++)
            {
                Console.WriteLine("  " + ref_[0 + j * 2].ToString().PadLeft(10)
                    + "  " + ref_[ 1 + j * 2].ToString().PadLeft(10)
                    + "  "
                    + "  " + phy[0 + j * 2].ToString().PadLeft(10)
                    + "  " + phy[1 + j * 2].ToString().PadLeft(10)
                    + "  "
                    + "  " + ref2[0 + j * 2].ToString().PadLeft(10)
                    + "  " + ref2[1 + j * 2].ToString().PadLeft(10) + "");
            }
        }

        static void test10()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST10 tests TRIANGLE_ORDER6_PHYSICAL_TO_REFERENCE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 10;

            int j;
            double[] phy = new double[2 * N];
            double[] ref_ = new double[2 * N];
            double[] ref2 = new double[2 * N];
            int seed;
            double[] t =  {
                7.0, 2.0,
                9.0, 2.0,
                7.0, 3.0,
                8.0, 2.0,
                8.0, 2.5,
                7.0, 2.5
            }
            ;

            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("TEST10");
            Console.WriteLine("  For an order 6 triangle,");
            Console.WriteLine("  TRIANGLE_ORDER6_PHYSICAL_TO_REFERENCE");
            Console.WriteLine("  maps a physical point to a reference point");
            Console.WriteLine("  TRIANGLE_ORDER6_REFERENCE_TO_PHYSICAL");
            Console.WriteLine("  maps a reference point to a physical point.");
            Console.WriteLine("");
            Console.WriteLine("   (  XSI    ETA ) ==> ( X      Y  ) ==> ( XSI2    ETA2 )");
            Console.WriteLine("");

            typeMethods.triangle_reference_sample(N, ref seed, ref ref_);

            Triangulation.triangle_order6_reference_to_physical(t, N, ref_, ref phy);

            Triangulation.triangle_order6_physical_to_reference(t, N, phy, ref ref2);

            for (j = 0; j < N; j++)
            {
                Console.WriteLine("  " + ref_[0 + j * 2].ToString().PadLeft(10)
                    + "  " + ref_[1 + j * 2].ToString().PadLeft(10)
                    + "  "
                    + "  " + phy[0 + j * 2].ToString().PadLeft(10)
                    + "  " + phy[1 + j * 2].ToString().PadLeft(10)
                    + "  "
                    + "  " + ref2[0 + j * 2].ToString().PadLeft(10)
                    + "  " + ref2[1 + j * 2].ToString().PadLeft(10) + "");
            }
        }

        static void test11()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST11 tests TRIANGULATION_NODE_ORDER.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 August 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int NODE_NUM = 36;
            int TRIANGLE_NUM = 41;
            int TRIANGLE_ORDER = 3;

            int[] node_order;
            int[] triangle_node =  {
                1, 8, 7,
                1, 2, 8,
                2, 9, 8,
                2, 3, 9,
                3, 10, 9,
                3, 4, 10,
                4, 11, 10,
                4, 5, 11,
                5, 12, 11,
                5, 6, 12,
                7, 14, 13,
                7, 8, 14,
                8, 15, 14,
                8, 9, 15,
                11, 18, 17,
                11, 12, 18,
                13, 20, 19,
                13, 14, 20,
                14, 21, 20,
                14, 15, 21,
                15, 22, 21,
                15, 16, 22,
                16, 23, 22,
                16, 17, 23,
                17, 24, 23,
                17, 18, 24,
                19, 26, 25,
                19, 20, 26,
                21, 28, 27,
                21, 22, 28,
                25, 30, 29,
                25, 26, 30,
                26, 31, 30,
                27, 32, 31,
                27, 28, 32,
                29, 34, 33,
                29, 30, 34,
                30, 35, 34,
                30, 31, 35,
                31, 36, 35,
                31, 32, 36
            }
            ;

            Console.WriteLine("");
            Console.WriteLine("TEST11");
            Console.WriteLine("  TRIANGULATION_NODE_ORDER computes the order");
            Console.WriteLine("  of the nodes in a triangulation.");

            node_order = NodeOrder.triangulation_node_order(TRIANGLE_ORDER, TRIANGLE_NUM,
                triangle_node, NODE_NUM);

            typeMethods.i4vec_print(NODE_NUM, node_order, "  NODE ORDER:");
        }

        static void test12()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST12 tests TRIANGULATION_ORDER3_ADJ_SET.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 August 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int[] adj;
            int adj_num;
            int[] adj_col;
            int hole_num = 0;
            int k;
            int node;
            int node_num = 0;
            double[] node_xy;
            int[] triangle_neighbor;
            int[] triangle_node;
            int triangle_order = 3;
            int triangle_num = 0;

            Console.WriteLine("");
            Console.WriteLine("TEST12");
            Console.WriteLine("  For an order3 triangulation:");
            Console.WriteLine("  TRIANGULATION_ORDER3_ADJ_COUNT counts adjacencies");
            Console.WriteLine("  TRIANGULATION_ORDER3_ADJ_SET sets adjacencies.");
            //
            //  Get the sizes.
            //
            TriangulationSampleData.triangulation_order3_example1_size(ref node_num, ref triangle_num, ref hole_num);

            adj_col = new int[node_num + 1];
            node_xy = new double[2 * node_num];
            triangle_neighbor = new int[3 * triangle_num];
            triangle_node = new int[triangle_order * triangle_num];
            //
            //  Get the example data.
            //
            TriangulationSampleData.triangulation_order3_example1(node_num, triangle_num, ref node_xy,
                ref triangle_node, ref triangle_neighbor);
            //
            //  Get the count of the adjacencies.
            //
            adj_num = Adjacency.triangulation_order3_adj_count(node_num, triangle_num,
                triangle_node, triangle_neighbor, adj_col);

            Console.WriteLine("");
            Console.WriteLine("  Number of adjacency entries is " + adj_num + "");

            Console.WriteLine("");
            Console.WriteLine("  Adjacency pointers:");
            Console.WriteLine("");
            for (node = 1; node <= node_num; node++)
            {
                Console.WriteLine("  " + node.ToString().PadLeft(8)
                    + "  " + adj_col[node - 1].ToString().PadLeft(8)
                    + "  " + (adj_col[node] - 1).ToString().PadLeft(8) + "");
            }

            //
            //  Get the adjacencies.
            //
            adj = Adjacency.triangulation_order3_adj_set(node_num, triangle_num, triangle_node,
                triangle_neighbor, adj_num, adj_col);
            //
            //  Print the adjacencies.
            //
            for (node = 1; node <= node_num; node++)
            {
                Console.WriteLine("");
                Console.WriteLine("  Nodes adjacent to node " + node + "");
                Console.WriteLine("");

                for (k = adj_col[node - 1]; k <= adj_col[node] - 1; k++)
                {
                    Console.WriteLine("  " + adj[k - 1].ToString().PadLeft(8) + "");
                }
            }
        }

        static void test125()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST125 tests TRIANGULATION_ORDER3_ADJ_SET2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 August 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int adj;
            int adj_num = 0;
            int[] adj_col;
            int hole_num = 0;
            int[] ia;
            int[] ja;
            int node;
            int node_num = 0;
            double[] node_xy;
            int[] triangle_neighbor;
            int[] triangle_node;
            int triangle_order = 3;
            int triangle_num = 0;

            Console.WriteLine("");
            Console.WriteLine("TEST125");
            Console.WriteLine("  For an order3 triangulation:");
            Console.WriteLine("  TRIANGULATION_ORDER3_ADJ_COUNT counts adjacencies");
            Console.WriteLine("  TRIANGULATION_ORDER3_ADJ_SET2 sets adjacencies");
            Console.WriteLine("  as a pair of vectors IA(*), JA(*).");
            //
            //  Get the sizes.
            //
            TriangulationSampleData.triangulation_order3_example1_size(ref node_num, ref triangle_num, ref hole_num);

            adj_col = new int[node_num + 1];
            node_xy = new double[2 * node_num];
            triangle_neighbor = new int[3 * triangle_num];
            triangle_node = new int[triangle_order * triangle_num];
            //
            //  Get the example data.
            //
            TriangulationSampleData.triangulation_order3_example1(node_num, triangle_num, ref node_xy,
                ref triangle_node, ref triangle_neighbor);
            //
            //  Get the count of the adjacencies.
            //
            adj_num = Adjacency.triangulation_order3_adj_count(node_num, triangle_num,
                triangle_node, triangle_neighbor, adj_col);

            Console.WriteLine("");
            Console.WriteLine("  Number of adjacency entries is " + adj_num + "");

            Console.WriteLine("");
            Console.WriteLine("  Adjacency pointers:");
            Console.WriteLine("");
            for (node = 1; node <= node_num; node++)
            {
                Console.WriteLine("  " + node.ToString().PadLeft(8)
                    + "  " + adj_col[node - 1].ToString().PadLeft(8)
                    + "  " + (adj_col[node] - 1).ToString().PadLeft(8) + "");
            }

            //
            //  Get the adjacencies.
            //
            ia = new int[adj_num];
            ja = new int[adj_num];

            Adjacency.triangulation_order3_adj_set2(node_num, triangle_num, triangle_node,
                triangle_neighbor, adj_num, adj_col, ia, ja);
            //
            //  Print the adjacencies.
            //
            Console.WriteLine("");
            Console.WriteLine("  Adjacency list:");
            Console.WriteLine("");
            for (adj = 0; adj < adj_num; adj++)
            {
                Console.WriteLine("  " + (adj + 1).ToString().PadLeft(8)
                    + "  (" + ia[adj].ToString().PadLeft(2)
                    + "," + ja[adj].ToString().PadLeft(2) + ")");
            }
        }

        static void test127()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST127 tests TRIANGULATION_ORDER3_ADJACENCY.
            //
            //  Discussion:
            //
            //    41--42--43--44  45--46--47--48
            //     | \ | \ | \ |   | \ | \ | \ |
            //    33--34--35--36  37--38--39--40
            //     | \ |                   | \ |
            //    29--30                  31--32
            //     | \ |                   | \ |
            //    25--26                  27--28
            //     | \ |                   | \ |
            //    21--22                  23--24
            //     | \ |                   | \ |
            //    17--18                  19--20
            //     | \ |                   | \ |
            //     9--10--11--12--13--14--15--16
            //     | \ | \ | \ | \ | \ | \ | \ |
            //     1---2---3---4---5---6---7---8
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 March 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int ELEMENT_NUM = 46;
            int NODE_NUM = 48;

            int[] adj;
            int[] element_node =  {
                1, 2, 9,
                2, 10, 9,
                2, 3, 10,
                3, 11, 10,
                3, 4, 11,
                4, 12, 11,
                4, 5, 12,
                5, 13, 12,
                5, 6, 13,
                6, 14, 13,
                6, 7, 14,
                7, 15, 14,
                7, 8, 15,
                8, 16, 15,
                9, 10, 17,
                10, 18, 17,
                15, 16, 19,
                16, 20, 19,
                17, 18, 21,
                18, 22, 21,
                19, 20, 23,
                20, 24, 23,
                21, 22, 25,
                22, 26, 25,
                23, 24, 27,
                24, 28, 27,
                25, 26, 29,
                26, 30, 29,
                27, 28, 31,
                28, 32, 31,
                29, 30, 33,
                30, 34, 33,
                31, 32, 39,
                32, 40, 39,
                33, 34, 41,
                34, 42, 41,
                34, 35, 42,
                35, 43, 42,
                35, 36, 43,
                36, 44, 43,
                37, 38, 45,
                38, 46, 45,
                38, 39, 46,
                39, 47, 46,
                39, 40, 47,
                40, 48, 47
            }
            ;
            int element_num = ELEMENT_NUM;
            int i;
            int j;
            int node_num = NODE_NUM;

            Console.WriteLine("");
            Console.WriteLine("TEST127");
            Console.WriteLine("  For an order3 triangulation:");
            Console.WriteLine("  TRIANGULATION_ORDER3_ADJACENCY sets the full");
            Console.WriteLine("  adjacency matrix.");

            for (j = 0; j < element_num; j++)
            {
                for (i = 0; i < 3; i++)
                {
                    element_node[i + j * 3] = element_node[i + j * 3] - 1;
                }
            }

            adj = Adjacency.triangulation_order3_adjacency(node_num, element_num, element_node);

            Console.WriteLine("");
            Console.WriteLine("  Adjacency matrix:");
            Console.WriteLine("");
            Console.WriteLine("                1         2         3         4       ");
            Console.WriteLine("      012345678901234567890123456789012345678901234567");
            Console.WriteLine("");
            for (i = 0; i < node_num; i++)
            {
                string cout = "  " + i.ToString().PadLeft(2) + "  ";
                for (j = 0; j < node_num; j++)
                {
                    cout += adj[i + j * node_num];
                }

                Console.WriteLine(cout);
            }
        }

        static void test13()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST13 tests TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 August 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;
            int NODE_NUM = 36;
            int TRIANGLE_NUM = 41;
            int TRIANGLE_ORDER = 3;

            int boundary_edge_num;
            string file_name = "triangulation_order3_plot2.eps";
            int node_show = 2;
            double[] node_xy =  {
                0.0, 0.0,
                1.0, 0.0,
                2.0, 0.0,
                3.0, 0.0,
                4.0, 0.0,
                5.0, 0.0,
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
                0.0, 5.0,
                1.0, 5.0,
                2.0, 5.0,
                3.0, 5.0,
                0.0, 6.0,
                1.0, 6.0,
                2.0, 6.0,
                3.0, 6.0
            }
            ;
            int[] triangle_node =  {
                1, 8, 7,
                1, 2, 8,
                2, 9, 8,
                2, 3, 9,
                3, 10, 9,
                3, 4, 10,
                4, 11, 10,
                4, 5, 11,
                5, 12, 11,
                5, 6, 12,
                7, 14, 13,
                7, 8, 14,
                8, 15, 14,
                8, 9, 15,
                11, 18, 17,
                11, 12, 18,
                13, 20, 19,
                13, 14, 20,
                14, 21, 20,
                14, 15, 21,
                15, 22, 21,
                15, 16, 22,
                16, 23, 22,
                16, 17, 23,
                17, 24, 23,
                17, 18, 24,
                19, 26, 25,
                19, 20, 26,
                21, 28, 27,
                21, 22, 28,
                25, 30, 29,
                25, 26, 30,
                26, 31, 30,
                27, 32, 31,
                27, 28, 32,
                29, 34, 33,
                29, 30, 34,
                30, 35, 34,
                30, 31, 35,
                31, 36, 35,
                31, 32, 36
            }
            ;
            int triangle_show = 2;

            Console.WriteLine("");
            Console.WriteLine("TEST13");
            Console.WriteLine("  For an order3 triangulation:");
            Console.WriteLine("  TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT counts the");
            Console.WriteLine("  boundary edges;");
            Console.WriteLine("  TRIANGULATION_ORDER3_PLOT plots the triangulation.");

            Plot.triangulation_order3_plot(file_name, NODE_NUM, node_xy,
                TRIANGLE_NUM, triangle_node, node_show, triangle_show);

            Console.WriteLine("");
            Console.WriteLine("  An Encapsulated PostScript image of this");
            Console.WriteLine("  triangulation is in \"" + file_name + "\".");

            boundary_edge_num = Boundary.triangulation_order3_boundary_edge_count(
                TRIANGLE_NUM, triangle_node);

            Console.WriteLine("");
            Console.WriteLine("  Number of boundary edges = " + boundary_edge_num + "");
            Console.WriteLine("  Correct number =           " + 33 + "");
        }

        static void test14()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST14 tests TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT_EULER.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 August 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int boundary_num;
            int hole_num = 2;
            int node_num = 36;
            int triangle_num = 41;

            Console.WriteLine("");
            Console.WriteLine("TEST14");
            Console.WriteLine("  For an order3 triangulation:");
            Console.WriteLine("  TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT_EULER");
            Console.WriteLine("  determines the number of edges that lie on the");
            Console.WriteLine("  boundary of a region that has been triangulated.");
            Console.WriteLine("");
            Console.WriteLine("  Number of points =         " + node_num + "");
            Console.WriteLine("  Number of triangles =      " + triangle_num + "");
            Console.WriteLine("  Number of holes =          " + hole_num + "");

            boundary_num = Boundary.triangulation_order3_boundary_edge_count_euler(node_num,
                triangle_num, hole_num);

            Console.WriteLine("  Number of boundary edges = " + boundary_num + "");
        }

        static void test15()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST15 tests TRIANGULATION_ORDER3_BOUNDARY_NODE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 August 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;
            int NODE_NUM = 36;
            int TRIANGLE_NUM = 41;
            int TRIANGLE_ORDER = 3;

            bool[] node_boundary;
            int[] triangle_node =  {
                1, 8, 7,
                1, 2, 8,
                2, 9, 8,
                2, 3, 9,
                3, 10, 9,
                3, 4, 10,
                4, 11, 10,
                4, 5, 11,
                5, 12, 11,
                5, 6, 12,
                7, 14, 13,
                7, 8, 14,
                8, 15, 14,
                8, 9, 15,
                11, 18, 17,
                11, 12, 18,
                13, 20, 19,
                13, 14, 20,
                14, 21, 20,
                14, 15, 21,
                15, 22, 21,
                15, 16, 22,
                16, 23, 22,
                16, 17, 23,
                17, 24, 23,
                17, 18, 24,
                19, 26, 25,
                19, 20, 26,
                21, 28, 27,
                21, 22, 28,
                25, 30, 29,
                25, 26, 30,
                26, 31, 30,
                27, 32, 31,
                27, 28, 32,
                29, 34, 33,
                29, 30, 34,
                30, 35, 34,
                30, 31, 35,
                31, 36, 35,
                31, 32, 36
            }
            ;

            Console.WriteLine("");
            Console.WriteLine("TEST15");
            Console.WriteLine("  For an order3 triangulation:");
            Console.WriteLine("  TRIANGULATION_ORDER3_BOUNDARY_NODE determines which");
            Console.WriteLine("  nodes lie on the boundary of a triangulation.");

            node_boundary = Boundary.triangulation_order3_boundary_node(NODE_NUM, TRIANGLE_NUM,
                triangle_node);

            typeMethods.lvec_print(NODE_NUM, node_boundary, "    Node  BN?");
        }

        static void test16()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST16 tests TRIANGULATION_ORDER3_CHECK.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 August 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int NODE_NUM = 13;
            int TRIANGLE_NUM = 16;
            int TRIANGLE_ORDER = 3;

            int ierror;
            int isave;
            int node_num2;
            int triangle_num2;
            int[] triangle_node =  {
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
            }
            ;

            Console.WriteLine("");
            Console.WriteLine("TEST16");
            Console.WriteLine("  For an order3 triangulation:");
            Console.WriteLine("  TRIANGULATION_ORDER3_CHECK checks the triangulation.");

            typeMethods.i4mat_transpose_print(TRIANGLE_ORDER, TRIANGLE_NUM, triangle_node,
                "  Triangles:");
            //
            //  Pass all tests.
            //
            ierror = Check.triangulation_order3_check(NODE_NUM, TRIANGLE_NUM, triangle_node);

            Console.WriteLine("  Error code = " + ierror + "");
            //
            //  Fail test 1.
            //
            node_num2 = 2;

            ierror = Check.triangulation_order3_check(node_num2, TRIANGLE_NUM,
                triangle_node);

            Console.WriteLine("  Error code = " + ierror + "");
            //
            //  Fail test 2.
            //
            triangle_num2 = 0;

            ierror = Check.triangulation_order3_check(NODE_NUM, triangle_num2,
                triangle_node);

            Console.WriteLine("  Error code = " + ierror + "");
            //
            //  Fail test 3.
            //
            isave = triangle_node[1 + 4 * 3];
            triangle_node[1 + 4 * 3] = 0;

            ierror = Check.triangulation_order3_check(NODE_NUM, TRIANGLE_NUM, triangle_node);

            Console.WriteLine("  Error code = " + ierror + "");
            triangle_node[1 + 4 * 3] = isave;
            //
            //  Fail test 4.
            //
            isave = triangle_node[2 + 9 * 3];
            triangle_node[2 + 9 * 3] = 2 * NODE_NUM + 1;

            ierror = Check.triangulation_order3_check(NODE_NUM, TRIANGLE_NUM, triangle_node);

            Console.WriteLine("  Error code = " + ierror + "");
            triangle_node[2 + 9 * 3] = isave;
            //
            //  Fail test 5.
            //
            triangle_node[2 + 3 * 3] = 3;
            triangle_node[2 + 7 * 3] = 3;
            triangle_node[2 + 9 * 3] = 3;
            triangle_node[2 + 10 * 3] = 3;
            triangle_node[1 + 13 * 3] = 3;

            ierror = Check.triangulation_order3_check(NODE_NUM, TRIANGLE_NUM, triangle_node);

            Console.WriteLine("  Error code = " + ierror + "");

            triangle_node[2 + 3 * 3] = 5;
            triangle_node[2 + 7 * 3] = 5;
            triangle_node[2 + 9 * 3] = 5;
            triangle_node[2 + 10 * 3] = 5;
            triangle_node[1 + 13 * 3] = 5;
            //
            //  Fail test 6.
            //
            triangle_node[0 + 8 * 3] = 7;

            ierror = Check.triangulation_order3_check(NODE_NUM, TRIANGLE_NUM, triangle_node);

            Console.WriteLine("  Error code = " + ierror + "");
            triangle_node[0 + 8 * 3] = 9;
            //
            //  Fail test 7.
            //
            triangle_node[2 + 6 * 3] = 2;

            ierror = Check.triangulation_order3_check(NODE_NUM, TRIANGLE_NUM, triangle_node);

            Console.WriteLine("  Error code = " + ierror + "");

            triangle_node[2 + 6 * 3] = 9;
        }

        static void test17()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST17 tests TRIANGULATION_ORDER3_EXAMPLE1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 August 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int hole_num = 0;
            int node_num = 0;
            double[] node_xy;
            int[] triangle_node;
            int[] triangle_neighbor;
            int triangle_order = 3;
            int triangle_num = 0;

            Console.WriteLine("");
            Console.WriteLine("TEST17");
            Console.WriteLine("  For an order3 triangulation:");
            Console.WriteLine("  TRIANGULATION_ORDER3_EXAMPLE1_SIZE gives the sizes");
            Console.WriteLine("  for an example triangulation;");
            Console.WriteLine("  TRIANGULATION_ORDER3_EXAMPLE1 returns the information");
            Console.WriteLine("  for an example triangulation;");
            Console.WriteLine("  TRIANGULATION_ORDER3_PRINT prints a triangulation.");
            //
            //  Get the sizes.
            //
            TriangulationSampleData.triangulation_order3_example1_size(ref node_num, ref triangle_num, ref hole_num);

            node_xy = new double[2 * node_num];
            triangle_node = new int[triangle_order * triangle_num];
            triangle_neighbor = new int[3 * triangle_num];
            //
            //  Get the data.
            //
            TriangulationSampleData.triangulation_order3_example1(node_num, triangle_num, ref node_xy,
                ref triangle_node, ref triangle_neighbor);

            Print.triangulation_order3_print(node_num, triangle_num, node_xy,
                triangle_node, triangle_neighbor);
        }

        static void test18()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST18 tests TRIANGULATION_ORDER3_NEIGHBOR.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 August 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int NODE_NUM = 13;
            int TRIANGLE_NUM = 16;
            int TRIANGLE_ORDER = 3;

            int s1;
            int s2 = 0;
            int t1;
            int t2 = 0;
            int[] triangle_node =  {
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
            }
            ;

            Console.WriteLine("");
            Console.WriteLine("TEST18");
            Console.WriteLine("  For an order3 triangulation:");
            Console.WriteLine("  TRIANGULATION_ORDER3_NEIGHBOR determines the");
            Console.WriteLine("  triangle neighbors.");
            Console.WriteLine("");
            Console.WriteLine("    T1    S1    T2    S2");
            Console.WriteLine("");

            for (t1 = 1; t1 <= TRIANGLE_NUM; t1++)
            {
                for (s1 = 1; s1 <= 3; s1++)
                {
                    Neighbor.triangulation_order3_neighbor(TRIANGLE_NUM, ref triangle_node,
                        t1, s1, ref t2, ref s2);

                    Console.WriteLine("  " + t1.ToString().PadLeft(4)
                        + "  " + s1.ToString().PadLeft(4)
                        + "  " + t2.ToString().PadLeft(4)
                        + "  " + s2.ToString().PadLeft(4) + "");
                }
            }
        }

        static void test19()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST19 tests TRIANGULATION_NEIGHBOR_ELEMENTS.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 September 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int TRIANGLE_NUM = 16;
            int TRIANGLE_ORDER = 3;

            int[] triangle_node =  {
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
            }
            ;
            int[] triangle_neighbor;

            Console.WriteLine("");
            Console.WriteLine("TEST19");
            Console.WriteLine("  TRIANGULATION_NEIGHBOR_ELEMENTS determines the");
            Console.WriteLine("  adjacency relationships between elements.");

            typeMethods.i4mat_transpose_print(TRIANGLE_ORDER, TRIANGLE_NUM, triangle_node,
                "  Elements:");

            triangle_neighbor = NeighborElements.triangulation_neighbor_elements(TRIANGLE_ORDER,
                TRIANGLE_NUM, triangle_node);

            typeMethods.i4mat_transpose_print(3, TRIANGLE_NUM, triangle_neighbor,
                "  Element neighbors:");
        }

        static void test20()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST20 tests TRIANGULATION_ORDER3_PLOT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 August 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            string file_name = "triangulation_order3_plot.eps";
            int hole_num = 0;
            int node_show = 0;
            int node_num = 0;
            double[] node_xy;
            int triangle_show = 2;
            int[] triangle_neighbor;
            int[] triangle_node;
            int triangle_num = 0;
            int triangle_order = 3;

            Console.WriteLine("");
            Console.WriteLine("TEST20");
            Console.WriteLine("  For an order3 triangulation:");
            Console.WriteLine("  TRIANGULATION_ORDER3_PLOT can plot a triangulation.");
            //
            //  Get the sizes.
            //
            TriangulationSampleData.triangulation_order3_example1_size(ref node_num, ref triangle_num, ref hole_num);

            Console.WriteLine("");
            Console.WriteLine("  Example data has " + node_num + " points,");
            Console.WriteLine("  organized into " + triangle_num + " triangles.");
            //
            //  Allocate space.
            //
            node_xy = new double[2 * node_num];
            triangle_node = new int[triangle_order * triangle_num];
            triangle_neighbor = new int[3 * triangle_num];
            //
            //  Get the example data.
            //
            TriangulationSampleData.triangulation_order3_example1(node_num, triangle_num, ref node_xy,
                ref triangle_node, ref triangle_neighbor);
            //
            //  Make the plot.
            //
            Plot.triangulation_order3_plot(file_name, node_num, node_xy, triangle_num,
                triangle_node, node_show, triangle_show);

            Console.WriteLine("");
            Console.WriteLine("  TRIANGULATION_ORDER3_PLOT has created an");
            Console.WriteLine("  Encapsulated PostScript file (EPS) containing");
            Console.WriteLine("  an image of the triangulation.");
            Console.WriteLine("");
            Console.WriteLine("  This file is called \"" + file_name + "\".");
        }

        static void test21()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST21 tests TRIANGULATION_ORDER3_PRINT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 August 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int NODE_NUM = 9;
            int TRIANGLE_NUM = 12;
            int TRIANGLE_ORDER = 3;

            double[] node_xy =  {
                0.0, 0.0,
                0.0, 1.0,
                0.2, 0.5,
                0.3, 0.6,
                0.4, 0.5,
                0.6, 0.4,
                0.6, 0.5,
                1.0, 0.0,
                1.0, 1.0
            }
            ;
            int[] triangle_node =  {
                2, 1, 3,
                3, 1, 6,
                2, 3, 4,
                4, 3, 5,
                7, 4, 5,
                5, 3, 6,
                7, 5, 6,
                9, 4, 7,
                6, 1, 8,
                7, 6, 8,
                7, 8, 9,
                2, 4, 9
            }
            ;
            int[] triangle_neighbor =  {
                -28, 2, 3,
                1, 9, 6,
                1, 4, 12,
                3, 6, 5,
                8, 4, 7,
                4, 2, 7,
                5, 6, 10,
                12, 5, 11,
                2, -34, 10,
                7, 9, 11,
                10, -38, 8,
                3, 8, -3
            }
            ;

            Console.WriteLine("");
            Console.WriteLine("TEST21");
            Console.WriteLine("  For an order3 triangulation:");
            Console.WriteLine("  TRIANGULATION_ORDER3_PRINT prints out a triangulation.");

            Print.triangulation_order3_print(NODE_NUM, TRIANGLE_NUM, node_xy,
                triangle_node, triangle_neighbor);
        }

        static void test213()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST213 tests TRIANGULATION_ORDER3_QUAD.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    22 January 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int QUAD_NUM = 6;

            int i;
            int j;
            int k;
            int n;
            int n11;
            int n12;
            int n21;
            int n22;
            int node_num;
            double[] node_xy;
            double quad_value = 0;
            double[] quad_w =  {
                0.1666666666666666,
                0.1666666666666666,
                0.1666666666666666,
                0.1666666666666666,
                0.1666666666666666,
                0.16666666666666660
            }
            ;
            double[] quad_xy =  {
                0.659027622374092, 0.231933368553031,
                0.659027622374092, 0.109039009072877,
                0.231933368553031, 0.659027622374092,
                0.231933368553031, 0.109039009072877,
                0.109039009072877, 0.659027622374092,
                0.109039009072877, 0.231933368553031
            }
            ;
            double region_area = 0;
            int test;
            int[] triangle_node;
            int test_num = 4;
            int triangle_order = 3;
            int triangle_num;
            double x;
            double y;

            Console.WriteLine("");
            Console.WriteLine("TEST213");
            Console.WriteLine("  TRIANGULATION_ORDER3_QUAD can apply a quadrature rule");
            Console.WriteLine("  to every triangle in a triangulated region,");
            Console.WriteLine("  and estimate the integral of a function over");
            Console.WriteLine("  that region.");
            Console.WriteLine("");
            Console.WriteLine("  NODE_NUM   TRI_NUM  Integral estim  Area of Region");
            Console.WriteLine("");

            for (test = 1; test <= test_num; test++)
            {
                //
                //  Set up the grid.
                //
                n = (int)Math.Pow(2, test - 1);
                node_num = (n + 1) * (n + 1);

                node_xy = new double[2 * node_num];

                k = 0;
                for (j = 1; j <= n + 1; j++)
                {
                    y = (double) (j - 1) / (double) (n + 1 - 1);
                    for (i = 1; i <= n + 1; i++)
                    {
                        x = (double) (i - 1) / (double) (n + 1 - 1);
                        node_xy[0 + k * 2] = x;
                        node_xy[1 + k * 2] = y;
                        k = k + 1;
                    }
                }

                //
                //  Set up the triangulation.
                //
                triangle_num = 2 * n * n;

                triangle_node = new int[triangle_order * triangle_num];

                k = 0;
                for (j = 1; j <= n; j++)
                {
                    for (i = 1; i <= n; i++)
                    {
                        n11 = i + (j - 1) * (n + 1);
                        n12 = i + (j + 1 - 1) * (n + 1);
                        n21 = i + 1 + (j - 1) * (n + 1);
                        n22 = i + 1 + (j + 1 - 1) * (n + 1);

                        triangle_node[0 + k * 3] = n11;
                        triangle_node[1 + k * 3] = n21;
                        triangle_node[2 + k * 3] = n12;
                        k = k + 1;

                        triangle_node[0 + k * 3] = n22;
                        triangle_node[1 + k * 3] = n12;
                        triangle_node[2 + k * 3] = n21;
                        k = k + 1;
                    }
                }

                //
                //  Estimate the integral.
                //
                Quad.triangulation_order3_quad(node_num, node_xy, triangle_order,
                    triangle_num, triangle_node, quad_fun, QUAD_NUM, quad_xy, quad_w,
                    ref quad_value, ref region_area);

                Console.WriteLine("  " + node_num.ToString().PadLeft(8)
                    + "  " + triangle_num.ToString().PadLeft(8)
                    + "  " + quad_value.ToString().PadLeft(14)
                    + "  " + region_area.ToString().PadLeft(14) + "");
            }
        }

        static double[] quad_fun(int n, double[] xy_vec, double[] f_vec )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    QUAD_FUN is a sample integrand function for TRIANGULATION_QUAD.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of evaluation points.
        //
        //    Input, double XY_VEC[2*N], the evaluation points.
        //
        //    Output, double F_VEC[N], the value of the integrand
        //    function at the evaluation points.
        //
        {
            int i;

            for (i = 0; i < n; i++)
            {
                f_vec[i] = Math.Exp(Math.Pow(xy_vec[0 + i * 2], 2)
                               + Math.Pow(xy_vec[1 + i * 2], 2));
            }

            return f_vec;
        }

        static void test215()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST215 tests TRIANGULATION_ORDER3_REFINE_COMPUTE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 January 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;
            int NODE_NUM1 = 5;
            int TRIANGLE_NUM1 = 3;
            int TRIANGLE_ORDER = 3;

            int[] edge_data;
            int node_num2 = 0;
            double[] node_xy1 =  {
                0.0, 0.0,
                1.0, 0.0,
                0.0, 1.0,
                1.0, 1.0,
                0.5, 1.5
            }
            ;
            double[] node_xy2;
            int[] triangle_node1 =  {
                1, 2, 3,
                4, 3, 2,
                3, 4, 5
            }
            ;
            int[] triangle_node2;
            int triangle_num2 = 0;

            Console.WriteLine("");
            Console.WriteLine("TEST215");
            Console.WriteLine("  For an order3 triangulation:");
            Console.WriteLine("  TRIANGULATION_ORDER3_REFINE_SIZE determines the");
            Console.WriteLine("  size of a refined triangulation.");
            Console.WriteLine("  TRIANGULATION_ORDER3_REFINE_COMPUTES computes the");
            Console.WriteLine("  refined triangulation.");

            Console.WriteLine("");
            Console.WriteLine("  The number of nodes is " + NODE_NUM1 + "");
            Console.WriteLine("  The number of triangles is " + TRIANGLE_NUM1 + "");

            typeMethods.r8mat_transpose_print(DIM_NUM, NODE_NUM1, node_xy1,
                "  The nodes");

            typeMethods.i4mat_transpose_print(TRIANGLE_ORDER, TRIANGLE_NUM1, triangle_node1,
                "  The triangles:");

            edge_data = new int[5 * (3 * TRIANGLE_NUM1)];

            Console.WriteLine("");
            Console.WriteLine("  Sizing the refined mesh:");

            Refine.triangulation_order3_refine_size(NODE_NUM1, TRIANGLE_NUM1,
                triangle_node1, ref node_num2, ref triangle_num2, ref edge_data);

            Console.WriteLine("");
            Console.WriteLine("  Information about the refined mesh:");
            Console.WriteLine("");
            Console.WriteLine("  The number of nodes is " + node_num2 + "");
            Console.WriteLine("  The number of triangles is " + triangle_num2 + "");

            Console.WriteLine("");
            Console.WriteLine("  Computing the refined mesh:");

            node_xy2 = new double[DIM_NUM * node_num2];
            triangle_node2 = new int[TRIANGLE_ORDER * triangle_num2];

            Refine.triangulation_order3_refine_compute(NODE_NUM1, TRIANGLE_NUM1,
                node_xy1, triangle_node1, node_num2, triangle_num2, edge_data, ref node_xy2,
                ref triangle_node2);

            typeMethods.r8mat_transpose_print(DIM_NUM, node_num2, node_xy2,
                "  The refined nodes");

            typeMethods.i4mat_transpose_print(TRIANGLE_ORDER, triangle_num2, triangle_node2,
                "  The refined triangles:");
        }

        static void test217()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST217 tests TRIANGULATION_SEARCH_DELAUNAY.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    18 August 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;
            int NODE_NUM = 13;
            int TEST_NUM = 10;
            int TRIANGLE_ORDER = 3;

            double alpha = 0;
            double beta = 0;
            double d1;
            double d2;
            double d3;
            double dist;
            double dnear = 0;
            int edge = 0;
            int error;
            double gamma = 0;
            int i;
            int i1;
            int i2;
            int i3;
            int nnear;
            double[] node_xy =  {
                0.0, 0.0,
                2.0, 2.0,
                -1.0, 3.0,
                -2.0, 2.0,
                8.0, 2.0,
                9.0, 5.0,
                7.0, 4.0,
                5.0, 6.0,
                6.0, 7.0,
                8.0, 8.0,
                11.0, 7.0,
                10.0, 4.0,
                6.0, 4.0
            }
            ;
            double[] p = new double[DIM_NUM];
            int seed;
            int step_num = 0;
            int[] td = new int[TEST_NUM];
            int test;
            int triangle_index = 0;
            int[] triangle_neighbor = new int[3 * 2 * NODE_NUM];
            int[] triangle_node = new int[TRIANGLE_ORDER * 2 * NODE_NUM];
            int triangle_num = 0;
            double[] xd = new double[DIM_NUM * TEST_NUM];

            Console.WriteLine("");
            Console.WriteLine("TEST217");
            Console.WriteLine("  Given a set of nodes NODE_XY, and a single point XD,");
            Console.WriteLine("  find the nearest node in NODE_XY to XD.");
            Console.WriteLine("");
            Console.WriteLine("  POINTS_POINT_NEAR_NAIVE_ND uses a naive method.");
            Console.WriteLine("  TRIANGULATION_SEARCH_DELAUNAY finds a triangle");
            Console.WriteLine("  containing the point.  Often, one of these vertices");
            Console.WriteLine("  is the closest point.");
            //
            //  Set up the Delaunay triangulation.
            //
            error = Delauney.r8tris2(NODE_NUM, ref node_xy, ref triangle_num, ref triangle_node,
                ref triangle_neighbor);

            if (error == 0)
            {
                Console.WriteLine("");
                Console.WriteLine("  R8TRIS2 computed the Delaunay triangulation.");
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("  R8TRIS2 returned an error condition.");
                return;
            }

            //
            //  Get the test points.
            //
            seed = 123456789;

            Sample.triangulation_order3_sample(NODE_NUM, node_xy, triangle_num,
                triangle_node, TEST_NUM, ref seed, ref xd, ref td);

            Console.WriteLine("");
            Console.WriteLine("              X         Y     Distance    Index     Steps");
            Console.WriteLine("");

            DelaunaySearchData data = new DelaunaySearchData();
            
            for (test = 0; test < TEST_NUM; test++)
            {
                for (i = 0; i < DIM_NUM; i++)
                {
                    p[i] = xd[i + test * DIM_NUM];
                }

                nnear = Helpers.points_point_near_naive_nd(DIM_NUM, NODE_NUM, node_xy,
                    p, ref dnear);

                Console.WriteLine("");
                Console.WriteLine("  XD       "
                     + p[0].ToString().PadLeft(8) + "  "
                     + p[1].ToString().PadLeft(8) + "");
                Console.WriteLine("  Naive    "
                     + node_xy[0 + (nnear - 1) * DIM_NUM].ToString().PadLeft(8) + "  "
                     + node_xy[1 + (nnear - 1) * DIM_NUM].ToString().PadLeft(8) + "  "
                     + dnear.ToString().PadLeft(8) + "  "
                     + nnear.ToString().PadLeft(8) + "");

                Search.triangulation_search_delaunay(ref data, NODE_NUM, node_xy, TRIANGLE_ORDER, triangle_num,
                    triangle_node, triangle_neighbor, p, ref triangle_index,
                    ref alpha, ref beta, ref gamma, ref edge, ref step_num);

                if (triangle_index < 1)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  Error: the search failed.");
                    continue;
                }

                i1 = triangle_node[0 + (triangle_index - 1) * TRIANGLE_ORDER];
                d1 = Math.Sqrt(Math.Pow(p[0] - node_xy[0 + i1 * 2], 2)
                          + Math.Pow(p[1] - node_xy[1 + i1 * 2], 2));

                dist = d1;
                nnear = i1;

                i2 = triangle_node[1 + (triangle_index - 1) * TRIANGLE_ORDER];
                d2 = Math.Sqrt(Math.Pow(p[0] - node_xy[0 + i2 * 2], 2)
                          + Math.Pow(p[1] - node_xy[1 + i2 * 2], 2));

                if (d2 < dist)
                {
                    dnear = d2;
                    nnear = i2;
                }

                i3 = triangle_node[2 + (triangle_index - 1) * TRIANGLE_ORDER];
                d3 = Math.Sqrt(Math.Pow(p[0] - node_xy[0 + i3 * 2], 2)
                          + Math.Pow(p[1] - node_xy[1 + i3 * 2], 2));

                if (d3 < dist)
                {
                    dnear = d3;
                    nnear = i3;
                }

                Console.WriteLine("  Delaunay "
                     + node_xy[0 + nnear * 2].ToString().PadLeft(9) + "  "
                     + node_xy[1 + nnear * 2].ToString().PadLeft(9) + "  "
                     + dnear.ToString().PadLeft(9) + "  "
                     + (nnear + 1).ToString().PadLeft(9)
                     + step_num.ToString().PadLeft(9) + "");
            }
        }

        static void test219()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST219 tests TRIANGULATION_SEARCH_DELAUNAY, TRIANGULATION_SEARCH_NAIVE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    01 August 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;
            int NODE_NUM = 13;
            int TEST_NUM = 10;
            int TRIANGLE_ORDER = 3;

            double alpha = 0;
            double beta = 0;
            int edge = 0;
            int error;
            double gamma = 0;
            double[] node_xy =  {
                0.0, 0.0,
                2.0, 2.0,
                -1.0, 3.0,
                -2.0, 2.0,
                8.0, 2.0,
                9.0, 5.0,
                7.0, 4.0,
                5.0, 6.0,
                6.0, 7.0,
                8.0, 8.0,
                11.0, 7.0,
                10.0, 4.0,
                6.0, 4.0
            }
            ;
            double[] p_test = new double[DIM_NUM * TEST_NUM];
            int seed = 123456789;
            int step_num = 0;
            int[] t_test = new int[TEST_NUM];
            int test;
            int triangle_index1 = 0;
            int triangle_index2 = 0;
            int[] triangle_neighbor = new int[3 * 2 * NODE_NUM];
            int[] triangle_node = new int[TRIANGLE_ORDER * 2 * NODE_NUM];
            int triangle_num = 0;

            Console.WriteLine("");
            Console.WriteLine("TEST219");
            Console.WriteLine("  Given a triangulation, and a point P,");
            Console.WriteLine("  find the triangle T containing to P.");
            Console.WriteLine("");
            Console.WriteLine("  TRIANGULATION_SEARCH_NAIVE uses a naive method.");
            Console.WriteLine("  TRIANGULATION_SEARCH_DELAUNAY uses a method that will work");
            Console.WriteLine("  fast if the triangulation is Delaunay.");
            //
            //  Set up the Delaunay triangulation.
            //
            error = Delauney.r8tris2(NODE_NUM, ref node_xy, ref triangle_num, ref triangle_node,
                ref triangle_neighbor);

            if (error == 0)
            {
                Console.WriteLine("");
                Console.WriteLine("  R8TRIS2 computed the Delaunay triangulation.");
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("  R8TRIS2 returned an error condition.");
                return;
            }

            //
            //  Get the test points.
            //
            Sample.triangulation_order3_sample(NODE_NUM, node_xy, triangle_num,
                triangle_node, TEST_NUM, ref seed, ref p_test, ref t_test);

            Console.WriteLine("");
            Console.WriteLine("         X           Y     Naive   Delaunay  Steps");
            Console.WriteLine("");

            DelaunaySearchData data = new DelaunaySearchData();
            
            for (test = 0; test < TEST_NUM; test++)
            {
                triangle_index1 = Search.triangulation_search_naive(NODE_NUM, node_xy,
                    TRIANGLE_ORDER, triangle_num, triangle_node, p_test, pIndex: + DIM_NUM * test);

                Search.triangulation_search_delaunay(ref data, NODE_NUM, node_xy, TRIANGLE_ORDER,
                    triangle_num, triangle_node, triangle_neighbor, p_test,
                    ref triangle_index2, ref alpha, ref beta, ref gamma, ref edge, ref step_num, pIndex:  + DIM_NUM * test);

                Console.WriteLine("  " + p_test[0 + test * DIM_NUM].ToString().PadLeft(10)
                    + "  " + p_test[1 + test * DIM_NUM].ToString().PadLeft(10)
                    + "  " + triangle_index1.ToString().PadLeft(8)
                    + "  " + triangle_index2.ToString().PadLeft(8)
                    + "  " + step_num.ToString().PadLeft(8) + "");

            }
        }

        static void test22()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST22 tests TRIANGULATION_ORDER6_ADJ_SET.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 August 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int[] adj;
            int adj_num;
            int[] adj_col;
            int hole_num = 0;
            int k;
            int node;
            int node_num = 0;
            double[] node_xy;
            int[] triangle_neighbor;
            int[] triangle_node;
            int triangle_order = 6;
            int triangle_num = 0;

            Console.WriteLine("");
            Console.WriteLine("TEST22");
            Console.WriteLine("  For an order6 triangulation:");
            Console.WriteLine("  TRIANGULATION_ORDER6_ADJ_COUNT counts adjacencies");
            Console.WriteLine("  TRIANGULATION_ORDER6_ADJ_SET sets adjacencies.");
            //
            //  Get the sizes.
            //
            TriangulationSampleData.triangulation_order6_example1_size(ref node_num, ref triangle_num, ref hole_num);

            adj_col = new int[node_num + 1];
            node_xy = new double[2 * node_num];
            triangle_neighbor = new int[3 * triangle_num];
            triangle_node = new int[triangle_order * triangle_num];
            //
            //  Get the example data.
            //
            TriangulationSampleData.triangulation_order6_example1(node_num, triangle_num, ref node_xy,
                ref triangle_node, ref triangle_neighbor);
            //
            //  Get the count of the adjacencies.
            //
            adj_num = Adjacency.triangulation_order6_adj_count(node_num, triangle_num,
                triangle_node, ref triangle_neighbor, ref adj_col);

            Console.WriteLine("");
            Console.WriteLine("  Number of adjacency entries is " + adj_num + "");

            Console.WriteLine("");
            Console.WriteLine("  Adjacency pointers:");
            Console.WriteLine("");
            for (node = 1; node <= node_num; node++)
            {
                Console.WriteLine("  " + node.ToString().PadLeft(8)
                    + "  " + adj_col[node - 1].ToString().PadLeft(8)
                    + "  " + (adj_col[node] - 1).ToString().PadLeft(8) + "");
            }

            //
            //  Get the adjacencies.
            //
            adj = Adjacency.triangulation_order6_adj_set(node_num, triangle_num, triangle_node,
                triangle_neighbor, adj_num, adj_col);
            //
            //  Print the adjacencies.
            //
            for (node = 1; node <= node_num; node++)
            {
                Console.WriteLine("");
                Console.WriteLine("  Nodes adjacent to node " + node + "");
                Console.WriteLine("");

                for (k = adj_col[node - 1]; k <= adj_col[node] - 1; k++)
                {
                    Console.WriteLine("  " + adj[k - 1].ToString().PadLeft(8) + "");
                }
            }
        }

        static void test23()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST23 tests TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 August 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int boundary_edge_num;
            int dim_num = 2;
            int hole_num = 0;
            int node_num = 0;
            double[] node_xy;
            int[] triangle_node;
            int[] triangle_neighbor;
            int triangle_num = 0;
            int triangle_order = 6;

            Console.WriteLine("");
            Console.WriteLine("TEST23");
            Console.WriteLine("  For an order6 triangulation:");
            Console.WriteLine("  TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT counts the");
            Console.WriteLine("  boundary edges.");

            TriangulationSampleData.triangulation_order6_example1_size(ref node_num, ref triangle_num, ref hole_num);

            node_xy = new double[dim_num * node_num];
            triangle_node = new int[triangle_order * triangle_num];
            triangle_neighbor = new int[3 * triangle_num];

            TriangulationSampleData.triangulation_order6_example1(node_num, triangle_num, ref node_xy,
                ref triangle_node, ref triangle_neighbor);

            boundary_edge_num = Boundary.triangulation_order6_boundary_edge_count(triangle_num,
                triangle_node);

            Console.WriteLine("");
            Console.WriteLine("  Number of boundary edges = " + boundary_edge_num + "");
            Console.WriteLine("  Correct number =           " + 16 + "");
        }

        static void test24()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST24 tests TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT_EULER.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 August 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int boundary_num = 0;
            int hole_num = 0;
            int node_num = 0;
            int triangle_num = 0;

            Console.WriteLine("");
            Console.WriteLine("TEST24");
            Console.WriteLine("  For an order6 triangulation:");
            Console.WriteLine("  TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT_EULER");
            Console.WriteLine("  determines the number of edges that lie on the");
            Console.WriteLine("  boundary of a region that has been triangulated.");

            TriangulationSampleData.triangulation_order6_example1_size(ref node_num, ref triangle_num, ref hole_num);

            Console.WriteLine("");
            Console.WriteLine("  Number of nodes =          " + node_num + "");
            Console.WriteLine("  Number of triangles =      " + triangle_num + "");
            Console.WriteLine("  Number of holes =          " + hole_num + "");

            boundary_num = Boundary.triangulation_order6_boundary_edge_count_euler(node_num,
                triangle_num, hole_num);

            Console.WriteLine("  Number of boundary edges = " + boundary_num + "");
            Console.WriteLine("  Correct number =           " + 16 + "");
        }

        static void test25()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST25 tests TRIANGULATION_ORDER6_BOUNDARY_NODE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 August 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            string file_name = "triangulation_order6_plot.eps";
            int i;
            int dim_num = 2;
            int hole_num = 0;
            bool[] node_boundary;
            int node_num = 0;
            int node_show = 2;
            double[] node_xy;
            int triangle_num = 0;
            int[] triangle_node;
            int[] triangle_neighbor;
            int triangle_order = 6;
            int triangle_show = 2;

            Console.WriteLine("");
            Console.WriteLine("TEST25");
            Console.WriteLine("  For an order6 triangulation:");
            Console.WriteLine("  TRIANGULATION_ORDER6_BOUNDARY_COUNT counts the boundary");
            Console.WriteLine("  edges.");
            Console.WriteLine("  TRIANGULATION_ORDER6_PLOT plots the triangulation.");

            TriangulationSampleData.triangulation_order6_example1_size(ref node_num, ref triangle_num, ref hole_num);

            node_xy = new double[dim_num * node_num];
            triangle_node = new int[triangle_order * triangle_num];
            triangle_neighbor = new int[3 * triangle_num];

            TriangulationSampleData.triangulation_order6_example1(node_num, triangle_num, ref node_xy,
                ref triangle_node, ref triangle_neighbor);
            //
            //  Make the plot.
            //
            Plot.triangulation_order6_plot(file_name, node_num, node_xy, triangle_num,
                triangle_node, node_show, triangle_show);

            Console.WriteLine("");
            Console.WriteLine("  An Encapsulated PostScript image of this");
            Console.WriteLine("  triangulation is in \"" + file_name + "\"");
            ;

            node_boundary = Boundary.triangulation_order6_boundary_node(node_num, triangle_num,
                triangle_node);

            Console.WriteLine("");
            Console.WriteLine("    Node  BN?");
            Console.WriteLine("");

            for (i = 1; i <= node_num; i++)
            {
                Console.WriteLine("  "
                     + i.ToString().PadLeft(6) + "  "
                     + node_boundary[i - 1] + "");
            }
        }

        static void test26()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST26 tests TRIANGULATION_ORDER6_PRINT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 August 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int hole_num = 0;
            int node_num= 0;
            double[] node_xy;
            int[] triangle_neighbor;
            int[] triangle_node;
            int triangle_num= 0;
            int triangle_order = 6;

            Console.WriteLine("");
            Console.WriteLine("TEST26");
            Console.WriteLine("  For an order6 triangulation:");
            Console.WriteLine("  TRIANGULATION_ORDER6_PRINT prints the data.");

            TriangulationSampleData.triangulation_order6_example1_size(ref node_num, ref triangle_num, ref hole_num);

            node_xy = new double[2 * node_num];
            triangle_neighbor = new int[3 * triangle_num];
            triangle_node = new int[triangle_order * triangle_num];

            TriangulationSampleData.triangulation_order6_example1(node_num, triangle_num, ref node_xy,
                ref triangle_node, ref triangle_neighbor);

            Print.triangulation_order6_print(node_num, triangle_num, node_xy,
                triangle_node, triangle_neighbor);
        }

        static void test265()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST265 tests TRIANGULATION_ORDER6_REFINE_COMPUTE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    10 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;
            int NODE_NUM1 = 12;
            int TRIANGLE_NUM1 = 3;
            int TRIANGLE_ORDER = 6;

            int[] edge_data = new int[1];
            int node_num2 = 0;
            double[] node_xy1 =  {
                0.0, 0.0,
                2.0, 0.0,
                0.0, 2.0,
                2.0, 2.0,
                1.0, 3.0,
                1.0, 0.0,
                0.0, 1.0,
                1.0, 1.0,
                2.0, 1.0,
                1.0, 2.0,
                0.5, 2.5,
                1.5, 2.5
            }
            ;
            double[] node_xy2;
            int[] triangle_node1 =  {
                1, 2, 3, 6, 8, 7,
                4, 3, 2, 9, 10, 8,
                3, 4, 5, 10, 12, 11
            }
            ;
            int[] triangle_node2;
            int triangle_num2 = 0;

            Console.WriteLine("");
            Console.WriteLine("TEST265");
            Console.WriteLine("  For an order6 triangulation:");
            Console.WriteLine("  TRIANGULATION_ORDER6_REFINE_SIZE determines the");
            Console.WriteLine("  size of a refined triangulation.");
            Console.WriteLine("  TRIANGULATION_ORDER6_REFINE_COMPUTES computes the");
            Console.WriteLine("  refined triangulation.");

            Console.WriteLine("");
            Console.WriteLine("  The number of nodes is " + NODE_NUM1 + "");
            Console.WriteLine("  The number of triangles is " + TRIANGLE_NUM1 + "");

            typeMethods.r8mat_transpose_print(DIM_NUM, NODE_NUM1, node_xy1,
                "  The nodes");

            typeMethods.i4mat_transpose_print(TRIANGLE_ORDER, TRIANGLE_NUM1, triangle_node1,
                "  The triangles:");

            edge_data = new int[5 * (3 * TRIANGLE_NUM1)];

            Console.WriteLine("");
            Console.WriteLine("  Sizing the refined mesh:");

            Refine.triangulation_order6_refine_size(NODE_NUM1, TRIANGLE_NUM1,
                triangle_node1, ref node_num2, ref triangle_num2, ref edge_data);

            Console.WriteLine("");
            Console.WriteLine("  Information about the refined mesh:");
            Console.WriteLine("");
            Console.WriteLine("  The number of nodes is " + node_num2 + "");
            Console.WriteLine("  The number of triangles is " + triangle_num2 + "");

            Console.WriteLine("");
            Console.WriteLine("  Computing the refined mesh:");

            node_xy2 = new double[DIM_NUM * node_num2];
            triangle_node2 = new int[TRIANGLE_ORDER * triangle_num2];

            Refine.triangulation_order6_refine_compute(NODE_NUM1, TRIANGLE_NUM1,
                node_xy1, triangle_node1, node_num2, triangle_num2, edge_data, ref node_xy2,
                ref triangle_node2);

            typeMethods.r8mat_transpose_print(DIM_NUM, node_num2, node_xy2,
                "  The refined nodes");

            typeMethods.i4mat_transpose_print(TRIANGLE_ORDER, triangle_num2, triangle_node2,
                "  The refined triangles:");

        }

        static void test27()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST27 tests TRIANGULATION_ORDER6_VERTEX_COUNT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 August 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int hole_num = 0;
            int midside_num;
            int node_num = 0;
            double[] node_xy;
            int[] triangle_neighbor;
            int[] triangle_node;
            int triangle_num = 0;
            int triangle_order = 6;
            int vertex_num;

            Console.WriteLine("");
            Console.WriteLine("TEST27");
            Console.WriteLine("  For an order6 triangulation:");
            Console.WriteLine("  TRIANGULATION_ORDER6_VERTEX_COUNT counts the ");
            Console.WriteLine("  vertex nodes and midside nodes.");

            TriangulationSampleData.triangulation_order6_example1_size(ref node_num, ref triangle_num, ref hole_num);

            node_xy = new double[2 * node_num];
            triangle_neighbor = new int[3 * triangle_num];
            triangle_node = new int[triangle_order * triangle_num];

            TriangulationSampleData.triangulation_order6_example1(node_num, triangle_num, ref node_xy,
                ref triangle_node, ref triangle_neighbor);

            vertex_num = VertexCount.triangulation_order6_vertex_count(triangle_num,
                triangle_node);

            midside_num = node_num - vertex_num;

            Console.WriteLine("");
            Console.WriteLine("  Number of nodes =         " + node_num + "");
            Console.WriteLine("  Number of vertex nodes =  " + vertex_num + "");
            Console.WriteLine("  Number of midside nodes = " + midside_num + "");
        }

        static void test31()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST31 tests VORONOI_POLYGON_AREA.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 August 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;
            int NEIGHBOR_NUM = 4;
            int NODE_NUM = 5;

            double area;
            double area_correct = 0.5;
            int[] neighbor_index =  {
                0, 1, 2, 3
            }
            ;
            int node = 4;
            double[] node_xy =  {
                0.0, 0.0,
                1.0, 0.0,
                1.0, 1.0,
                0.0, 1.0,
                0.5, 0.5
            }
            ;

            Console.WriteLine("");
            Console.WriteLine("TEST31");
            Console.WriteLine("  VORONOI_POLYGON_AREA computes the area of");
            Console.WriteLine("  a finite Voronoi polygon.");

            area = Voronoi.voronoi_polygon_area(node, NEIGHBOR_NUM, neighbor_index,
                NODE_NUM, node_xy);

            Console.WriteLine("");
            Console.WriteLine("  The computed area is " + area + "");
            Console.WriteLine("  The correct area is  " + area_correct + "");
        }

        static void test32()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST32 tests VORONOI_POLYGON_CENTROID.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 August 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;
            int NEIGHBOR_NUM = 4;
            int NODE_NUM = 5;

            double[] centroid;
            double[] centroid_exact =  {
                0.5, 0.5
            }
            ;
            int[] neighbor_index =  {
                0, 1, 2, 3
            }
            ;
            int node = 4;
            double[] node_xy =  {
                0.0, 0.0,
                1.0, 0.0,
                1.0, 1.0,
                0.0, 1.0,
                0.5, 0.5
            }
            ;

            Console.WriteLine("");
            Console.WriteLine("TEST32");
            Console.WriteLine("  VORONOI_POLYGON_CENTROID computes the centroid of");
            Console.WriteLine("  a finite Voronoi polygon.");
            Console.WriteLine();

            centroid = Voronoi.voronoi_polygon_centroid(node, NEIGHBOR_NUM,
                neighbor_index, NODE_NUM, node_xy);

            Console.WriteLine("");
            Console.WriteLine("  The computed centroid is "
                 + centroid[0].ToString().PadLeft(10) + "  "
                 + centroid[1].ToString().PadLeft(10) + "");
            Console.WriteLine("  The correct centroid is  "
                 + centroid_exact[0].ToString().PadLeft(10) + "  "
                 + centroid_exact[1].ToString().PadLeft(10) + "");
        }

        static void test33()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST33 tests VORONOI_POLYGON_VERTICES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 August 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;
            int NEIGHBOR_NUM = 4;
            int NODE_NUM = 5;

            int[] neighbor_index =  {
                0, 1, 2, 3
            }
            ;
            int node = 4;
            double[] v = new double[DIM_NUM * NEIGHBOR_NUM];
            double[] node_xy =  {
                0.0, 0.0,
                1.0, 0.0,
                1.0, 1.0,
                0.0, 1.0,
                0.5, 0.5
            }
            ;

            Console.WriteLine("");
            Console.WriteLine("TEST33");
            Console.WriteLine("  VORONOI_POLYGON_VERTICES computes the vertices of");
            Console.WriteLine("  a finite Voronoi polygon.");
            Console.WriteLine();

            Voronoi.voronoi_polygon_vertices(node, NEIGHBOR_NUM, neighbor_index,
                NODE_NUM, node_xy, ref v);

            typeMethods.r8mat_transpose_print(DIM_NUM, NEIGHBOR_NUM, v, "  Vertices:");
        }
    }
}