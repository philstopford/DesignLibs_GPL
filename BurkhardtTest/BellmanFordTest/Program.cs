using System;
using Burkardt.MinDist;
using Burkardt.Types;

namespace BellmanFordTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for BELLMAN_FORD_TEST.
        //
        //  Discussion:
        //
        //    BELLMAN_FORD_TEST tests the BELLMAN_FORD library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 November 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("BELLMAN_FORD_TEST");
        Console.WriteLine("  Test the BELLMAN_FORD library.");

        test01();

        Console.WriteLine("");
        Console.WriteLine("BELLMAN_FORD_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 runs a simple test.
        //
        //  Discussion:
        //
        //    The correct distances are { 0, -6, -2, 3, 0, 0 }.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 November 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] e =
            {
                1, 0,
                4, 1,
                1, 2,
                2, 4,
                4, 0,
                2, 5,
                5, 0,
                3, 2,
                5, 3,
                3, 0,
            }
            ;
        int e_num = 10;
        double[] e_weight =
            {
                -3.0,
                6.0,
                -4.0,
                -1.0,
                4.0,
                -2.0,
                2.0,
                8.0,
                -3.0,
                3.0
            }
            ;
        int[] predecessor = new int[6];
        int source = 0;
        int v_num = 6;
        double[] v_weight = new double[6];

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  Bellman-Ford shortest path algorithm.");

        Console.WriteLine("");
        Console.WriteLine("  Number of vertices = " + v_num + "");
        Console.WriteLine("  Number of edges = " + e_num + "");
        Console.WriteLine("  The reference vertex is " + source + "");
        typeMethods.i4mat_transpose_print(2, e_num, e, "  The edge array:");
        typeMethods.r8vec_print(e_num, e_weight, "  The edge weights:");

        BellmanFord.bellman_ford(v_num, e_num, source, e, e_weight, ref v_weight, ref predecessor);

        typeMethods.r8vec_print(v_num, v_weight, "  The shortest distances:");

        typeMethods.i4vec_print(v_num, predecessor, "  The vertex predecessor parents for the shortest paths:");
    }
}