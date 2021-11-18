using System;
using Burkardt.RankingNS;
using Burkardt.Types;

namespace ComboTest;

internal partial class Program
{
    private static void edge_check_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EDGE_CHECK_TEST tests EDGE_CHECK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        bool check;
        int edge_num = 0;
        int[] edge_list = new int[1];
        int[] edge_list1 =  {
                1, 2,
                2, 3,
                3, 1
            }
            ;
        int[] edge_list3 =  {
                1, 2,
                2, 3,
                3, 4
            }
            ;
        int[] edge_list4 =  {
                1, 2,
                2, 2,
                3, 1
            }
            ;
        int[] edge_list5 =  {
                1, 2,
                2, 3,
                2, 1
            }
            ;
        int[] edge_list6 =  {
                1, 2,
                2, 3,
                3, 1
            }
            ;
        int node_num = 0;
        int test;

        Console.WriteLine("");
        Console.WriteLine("EDGE_CHECK TEST");
        Console.WriteLine("  EDGE_CHECK checks a graph described by edges.");
        Console.WriteLine("");
        Console.WriteLine("  Check?  Nodes  Edges    EdgeList");
        Console.WriteLine("");

        for (test = 1; test <= 6; test++)
        {
            switch (test)
            {
                case 1:
                    node_num = -5;
                    edge_num = 3;
                    edge_list = typeMethods.i4vec_copy_new(2 * edge_num, edge_list1);
                    break;
                case 2:
                    node_num = 3;
                    edge_num = -1;
                    edge_list = null;
                    break;
                case 3:
                    node_num = 3;
                    edge_num = 3;
                    edge_list = typeMethods.i4vec_copy_new(2 * edge_num, edge_list3);
                    break;
                case 4:
                    node_num = 3;
                    edge_num = 3;
                    edge_list = typeMethods.i4vec_copy_new(2 * edge_num, edge_list4);
                    break;
                case 5:
                    node_num = 3;
                    edge_num = 3;
                    edge_list = typeMethods.i4vec_copy_new(2 * edge_num, edge_list5);
                    break;
                case 6:
                    node_num = 3;
                    edge_num = 3;
                    edge_list = typeMethods.i4vec_copy_new(2 * edge_num, edge_list6);
                    break;
            }

            Console.WriteLine("");
            check = Ranking.edge_check(node_num, edge_num, edge_list);
            Console.WriteLine("      " + check.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                       + "     " + node_num.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                       + "     " + edge_num.ToString(CultureInfo.InvariantCulture).PadLeft(2) + "");
            typeMethods.i4mat_print(2, edge_num, edge_list, "  Edge list of graph:");

        }
    }

    private static void edge_degree_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EDGE_DEGREE_TEST tests EDGE_DEGREE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] d;
        int[] edge =  {
                1, 2, 2, 3, 4,
                2, 3, 4, 4, 5
            }
            ;
        int edge_num;
        int node_num;

        Console.WriteLine("");
        Console.WriteLine("EDGE_DEGREE_TEST");
        Console.WriteLine("  EDGE_DEGREE determines the degree of each node in a graph.");

        node_num = 5;
        edge_num = 5;

        typeMethods.i4mat_print(2, edge_num, edge, "  The edge array:");

        d = Ranking.edge_degree(node_num, edge_num, edge);

        typeMethods.i4vec_print(node_num, d, "  The degree vector:");

    }

    private static void edge_enum_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EDGE_ENUM_TEST tests EDGE_ENUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int edge_num;
        int node_num;

        Console.WriteLine("");
        Console.WriteLine("EDGE_ENUM_TEST");
        Console.WriteLine("  EDGE_ENUM enumerates the maximum number of edges");
        Console.WriteLine("  possible in a graph of NODE_NUM nodes.");
        Console.WriteLine("");
        Console.WriteLine("   NODE_NUM    EDGE_NUM(max)");
        Console.WriteLine("");

        for (node_num = 1; node_num <= 10; node_num++)
        {
            edge_num = Ranking.edge_enum(node_num);
            Console.WriteLine("       "
                              + "  " + node_num.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                              + "    "
                              + "  " + edge_num.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
        }
    }
}