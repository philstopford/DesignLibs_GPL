﻿using System;
using System.Globalization;
using Burkardt.MinDist;
using Burkardt.Types;

namespace FloydTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for FLOYD_TEST.
        //
        //  Discussion:
        //
        //    FLOYD_TEST tests the FLOYD library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 March 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("FLOYD_TEST");
            
        Console.WriteLine("  Test the FLOYD library.");

        test01();
        test02();
        test03();
        test04();

        Console.WriteLine("");
        Console.WriteLine("FLOYD_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests I4MAT_FLOYD.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 November 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 6;

        int[] a =  {
                0, -1, -1, -1, -1, -1,
                2, 0, -1, -1, -1, 5,
                5, 7, 0, -1, 2, -1,
                -1, 1, 4, 0, -1, 2,
                -1, -1, -1, 3, 0, 4,
                -1, 8, -1, -1, 3, 0
            }
            ;
        int i;
        int j;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  I4MAT_FLOYO uses Floyd's algorithm to find the");
        Console.WriteLine("  shortest distance between all pairs of nodes");
        Console.WriteLine("  in a directed graph, starting from the initial array");
        Console.WriteLine("  of direct node-to-node distances.");

        Console.WriteLine("");
        Console.WriteLine("  In the initial direct distance array, if");
        Console.WriteLine("    A(I,J) = -1,");
        Console.WriteLine("  this indicates there is NO directed link from");
        Console.WriteLine("  node I to node J.  In that case, the value of");
        Console.WriteLine("  of A(I,J) is essentially \"infinity\".");

        typeMethods.i4mat_print(N, N, a, "  Initial direct distance array:");

        int huge = typeMethods.i4_huge() / 2;

        for (j = 0; j < N; j++)
        {
            for (i = 0; i < N; i++)
            {
                a[i + j * N] = a[i + j * N] switch
                {
                    -1 => huge,
                    _ => a[i + j * N]
                };
            }
        }

        Floyd.i4mat_floyd(N, ref a);

        for (j = 0; j < N; j++)
        {
            for (i = 0; i < N; i++)
            {
                if (a[i + j * N] == huge)
                {
                    a[i + j * N] = -1;
                }
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  In the final shortest distance array, if");
        Console.WriteLine("    A(I,J) = -1,");
        Console.WriteLine("  this indicates there is NO directed path from");
        Console.WriteLine("  node I to node J.");

        typeMethods.i4mat_print(N, N, a, "  Final shortest distance array:");
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests R8MAT_FLOYD.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 November 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 6;

        double[] a =  {
                0.0, -1.0, -1.0, -1.0, -1.0, -1.0,
                2.0, 0.0, -1.0, -1.0, -1.0, 5.0,
                5.0, 7.0, 0.0, -1.0, 2.0, -1.0,
                -1.0, 1.0, 4.0, 0.0, -1.0, 2.0,
                -1.0, -1.0, -1.0, 3.0, 0.0, 4.0,
                -1.0, 8.0, -1.0, -1.0, 3.0, 0.0
            }
            ;
        int i;
        int j;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  R8MAT_FLOYO uses Floyd's algorithm to find the");
        Console.WriteLine("  shortest distance between all pairs of nodes");
        Console.WriteLine("  in a directed graph, starting from the initial array");
        Console.WriteLine("  of direct node-to-node distances.");

        Console.WriteLine("");
        Console.WriteLine("  In the initial direct distance array, if");
        Console.WriteLine("    A(I,J) = -1,");
        Console.WriteLine("  this indicates there is NO directed link from");
        Console.WriteLine("  node I to node J.  In that case, the value of");
        Console.WriteLine("  of A(I,J) is essentially \"infinity\".");

        typeMethods.r8mat_print(N, N, a, "  Initial direct distance array:");

        double huge = typeMethods.r8_huge();

        for (j = 0; j < N; j++)
        {
            for (i = 0; i < N; i++)
            {
                a[i + j * N] = a[i + j * N] switch
                {
                    -1.0 => huge,
                    _ => a[i + j * N]
                };
            }
        }

        Floyd.r8mat_floyd(N, ref a);

        for (j = 0; j < N; j++)
        {
            for (i = 0; i < N; i++)
            {
                if (Math.Abs(a[i + j * N] - huge) <= typeMethods.r8_epsilon())
                {
                    a[i + j * N] = -1.0;
                }
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  In the final shortest distance array, if");
        Console.WriteLine("    A(I,J) = -1,");
        Console.WriteLine("  this indicates there is NO directed path from");
        Console.WriteLine("  node I to node J.");

        typeMethods.r8mat_print(N, N, a, "  Final shortest distance array:");
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 applies Floyd's method to problems of increasing size.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 March 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  Test I4MAT_FLOYD on the MOD(I,J) matrix.");
        Console.WriteLine("  The work is roughly N^3.");
        Console.WriteLine("");
        Console.WriteLine("         N   Time (seconds)  Time/N^3");
        Console.WriteLine("");

        int n = 1;
        while (n <= 2048)
        {
            double wtime = test03_sub(n);
            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + wtime.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + (1000000.0 * wtime / n
                                                               / n / n).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "");
            n *= 2;
        }
    }

    private static double test03_sub(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03_SUB tests I4MAT_FLOYD.
        //
        //  Discussion:
        //
        //    The matrix size is input by the user.
        //
        //    The matrix A has the property that
        //
        //      A(I,J) = 1 if I is divisible by J.
        //
        //  Example:
        //
        //    N = 6
        //
        //    1 0 0 0 0 0
        //    1 1 0 0 0 0
        //    1 0 1 0 0 0
        //    1 1 0 1 0 0
        //    1 0 0 0 1 0
        //    1 1 1 0 0 1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 July 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the size of the matrix.
        //
        //    Output, double TEST03, the CPU time required by I4MAT_FLOYD.
        //
    {
        int j;

        int[] a = new int[n * n];

        int huge = typeMethods.i4_huge() / 2;

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < n; i++)
            {
                a[i + j * n] = ((i + 1) % (j + 1)) switch
                {
                    0 => 1,
                    _ => huge
                };
            }
        }

        DateTime time1 = DateTime.Now;

        Floyd.i4mat_floyd(n, ref a);

        DateTime time2 = DateTime.Now;

        double wtime = (time2 - time1).TotalSeconds;

        return wtime;
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 uses Floyd's method for a triangulation.
        //
        //  Discussion:
        //
        //     8  41--42--43--44  45--46--47--48
        //     |   | \ | \ | \ |   | \ | \ | \ |
        //     7  33--34--35--36  37--38--39--40
        //     |   | \ |                   | \ |
        //     6  29--30                  31--32
        //     |   | \ |                   | \ |
        //     5  25--26                  27--28
        //     |   | \ |                   | \ |
        //     4  21--22                  23--24
        //     |   | \ |                   | \ |
        //     3  17--18                  19--20
        //     |   | \ |                   | \ |
        //     2   9--10--11--12--13--14--15--16
        //     |   | \ | \ | \ | \ | \ | \ | \ |
        //     1   1---2---3---4---5---6---7---8
        //     |    
        //     +---1---2---3---4---5---6---7---8
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 March 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int ELEMENT_NUM = 46;
        int NODE_NUM = 48;

        double[] d = new double[NODE_NUM * NODE_NUM];
        int element;
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
        int i;
        int j;
        double[] xy =  {
                1.0, 1.0,
                2.0, 1.0,
                3.0, 1.0,
                4.0, 1.0,
                5.0, 1.0,
                6.0, 1.0,
                7.0, 1.0,
                8.0, 1.0,
                1.0, 2.0,
                2.0, 2.0,
                3.0, 2.0,
                4.0, 2.0,
                5.0, 2.0,
                6.0, 2.0,
                7.0, 2.0,
                8.0, 2.0,
                1.0, 3.0,
                2.0, 3.0,
                7.0, 3.0,
                8.0, 3.0,
                1.0, 4.0,
                2.0, 4.0,
                7.0, 4.0,
                8.0, 4.0,
                1.0, 5.0,
                2.0, 5.0,
                7.0, 5.0,
                8.0, 5.0,
                1.0, 6.0,
                2.0, 6.0,
                7.0, 6.0,
                8.0, 6.0,
                1.0, 7.0,
                2.0, 7.0,
                3.0, 7.0,
                4.0, 7.0,
                5.0, 7.0,
                6.0, 7.0,
                7.0, 7.0,
                8.0, 7.0,
                1.0, 8.0,
                2.0, 8.0,
                3.0, 8.0,
                4.0, 8.0,
                5.0, 8.0,
                6.0, 8.0,
                7.0, 8.0,
                8.0, 8.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  Start with a triangulation, and use R8_FLOYD");
        Console.WriteLine("  to determine the pairwise distance matrix.");
        //
        //  Must initialize distances to -1
        //
        for (j = 0; j < NODE_NUM; j++)
        {
            for (i = 0; i < NODE_NUM; i++)
            {
                d[i + j * NODE_NUM] = -1.0;
            }
        }

        //
        //  Diagonals are 0.
        //
        for (i = 0; i < NODE_NUM; i++)
        {
            d[i + i * NODE_NUM] = 0.0;
        }

        //
        //  Initialize D to the one-step distance.
        //
        for (element = 0; element < ELEMENT_NUM; element++)
        {
            int n1 = element_node[2 + element * 3] - 1;
            for (i = 0; i < 3; i++)
            {
                int n2 = element_node[i + element * 3] - 1;
                d[n1 + n2 * NODE_NUM] = Math.Sqrt(Math.Pow(xy[0 + n1 * 2] - xy[0 + n2 * 2], 2)
                                                  + Math.Pow(xy[1 + n1 * 2] - xy[1 + n2 * 2], 2));
                d[n2 + n1 * NODE_NUM] = d[n1 + n2 * NODE_NUM];
                n1 = n2;
            }
        }

        //
        //  Reset -1 values to R8_HUGE.
        //
        for (j = 0; j < NODE_NUM; j++)
        {
            for (i = 0; i < NODE_NUM; i++)
            {
                d[i + j * NODE_NUM] = d[i + j * NODE_NUM] switch
                {
                    -1.0 => typeMethods.r8_huge(),
                    _ => d[i + j * NODE_NUM]
                };
            }
        }

        //
        //  Update D to the N-1 step distance.
        //
        Floyd.r8mat_floyd(NODE_NUM, ref d);

        typeMethods.r8mat_print(NODE_NUM, NODE_NUM, d, "  Distance matrix");
    }
}