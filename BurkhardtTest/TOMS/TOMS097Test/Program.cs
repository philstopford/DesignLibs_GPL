using System;
using Burkardt.Types;

namespace TOMS097Test;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TOMS097_TEST.
        //
        //  Discussion:
        //
        //    TOMS097_TEST tests the TOMS097 library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 March 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {

        Console.WriteLine("");
        Console.WriteLine("TOMS097_TEST");
        Console.WriteLine("  Test the TOMS097 library.");

        i4mat_shortest_path_test();
        r8mat_shortest_path_test();

        Console.WriteLine("");
        Console.WriteLine("TOMS097_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void i4mat_shortest_path_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4MAT_SHORTEST_PATH_TEST tests I4MAT_SHORTEST_PATH.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 March 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 6;

        int[] a =
        {
            0, -1, -1, -1, -1, -1,
            2, 0, -1, -1, -1, 5,
            5, 7, 0, -1, 2, -1,
            -1, 1, 4, 0, -1, 2,
            -1, -1, -1, 3, 0, 4,
            -1, 8, -1, -1, 3, 0
        };
        int i;
        int j;
        int n = N;
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("I4MAT_SHORTEST_PATH_TEST");
        Console.WriteLine("  I4MAT_SHORTEST_PATH uses Floyd''s algorithm to find the");
        Console.WriteLine("  shortest distance between all pairs of nodes");
        Console.WriteLine("  in a directed graph, starting from the initial array");
        Console.WriteLine("  of direct node-to-node distances.");

        Console.WriteLine("");
        Console.WriteLine("  In the initial direct distance array, if");
        Console.WriteLine("    A(I,J) = HUGE,");
        Console.WriteLine("  this indicates there is NO directed link from");
        Console.WriteLine("  node I to node J.  In that case, the value of");
        Console.WriteLine("  of A(I,J) is essentially 'infinity'.");

        Console.WriteLine("");
        Console.WriteLine("  Initial direct-link distance matrix:");
        Console.WriteLine("");

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                cout += a[i + j * n].ToString().PadLeft(6);
            }

            Console.WriteLine(cout);
            cout = "";
        }

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                a[i + j * n] = a[i + j * n] switch
                {
                    -1 => typeMethods.i4_huge(),
                    _ => a[i + j * n]
                };
            }
        }

        typeMethods.i4mat_shortest_path(n, ref a);

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                if (a[i + j * n] == typeMethods.i4_huge())
                {
                    a[i + j * n] = -1;
                }
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  In the final shortest distance array, if");
        Console.WriteLine("    A(I,J) = -1,");
        Console.WriteLine("  this indicates there is NO directed path from");
        Console.WriteLine("  node I to node J.");

        Console.WriteLine("");
        Console.WriteLine("  Final distance matrix:");
        Console.WriteLine("");

        cout = "";
            
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                cout += a[i + j * n].ToString().PadLeft(6);
            }

            Console.WriteLine(cout);
            cout = "";
        }
    }

    private static void r8mat_shortest_path_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_SHORTEST_PATH_TEST tests R8MAT_SHORTEST_PATH.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 March 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 6;

        double[] a =
        {
            0.0, -1.0, -1.0, -1.0, -1.0, -1.0,
            2.0, 0.0, -1.0, -1.0, -1.0, 5.0,
            5.0, 7.0, 0.0, -1.0, 2.0, -1.0,
            -1.0, 1.0, 4.0, 0.0, -1.0, 2.0,
            -1.0, -1.0, -1.0, 3.0, 0.0, 4.0,
            -1.0, 8.0, -1.0, -1.0, 3.0, 0.0
        };
        int i;
        int j;
        int n = N;
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("R8MAT_SHORTEST_PATH_TEST");
        Console.WriteLine("  R8MAT_SHORTEST_PATH uses Floyd''s algorithm to find the");
        Console.WriteLine("  shortest distance between all pairs of nodes");
        Console.WriteLine("  in a directed graph, starting from the initial array");
        Console.WriteLine("  of direct node-to-node distances.");

        Console.WriteLine("");
        Console.WriteLine("  In the initial direct distance array, if");
        Console.WriteLine("    A(I,J) = -1,");
        Console.WriteLine("  this indicates there is NO directed link from");
        Console.WriteLine("  node I to node J.  In that case, the value of");
        Console.WriteLine("  of A(I,J) is essentially 'infinity'.");

        Console.WriteLine("");
        Console.WriteLine("  Initial direct-link distance matrix:");
        Console.WriteLine("");

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                cout += a[i + j * n].ToString().PadLeft(10);
            }

            Console.WriteLine(cout);
            cout = "";
        }

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                a[i + j * n] = a[i + j * n] switch
                {
                    -1.0 => typeMethods.r8_huge(),
                    _ => a[i + j * n]
                };
            }
        }

        typeMethods.r8mat_shortest_path(n, ref a);

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                if (a[i + j * n] == typeMethods.r8_huge())
                {
                    a[i + j * n] = -1.0;
                }
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  In the final shortest distance array, if");
        Console.WriteLine("    A(I,J) = -1,");
        Console.WriteLine("  this indicates there is NO directed path from");
        Console.WriteLine("  node I to node J.");

        Console.WriteLine("");
        Console.WriteLine("  Final distance matrix:");
        Console.WriteLine("");

        cout = "";
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                cout += a[i + j * n].ToString().PadLeft(10);
            }

            Console.WriteLine(cout);
            cout = "";
        }

    }
}