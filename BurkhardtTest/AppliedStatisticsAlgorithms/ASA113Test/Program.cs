using System;
using System.Globalization;
using Burkardt;
using Burkardt.AppliedStatistics;

namespace ASA113Test;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for ASA113_TEST.
        //
        //  Discussion:
        //
        //    ASA136_TEST tests the ASA113 library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("ASA113_TEST:");
        Console.WriteLine("  Test the ASA113 library.");

        test01();

        Console.WriteLine("");
        Console.WriteLine("ASA113_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tries out the ASA113 routine.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int ci;
        int ifault = 0;
        const int k = 5;
        const int m = 100;
        const int n = 2;

        int[] c = new int[m];
        double[] c_center = new double[k * n];
        int[] c_size = new int[k];
        double[] wss = new double[k];

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Test the ASA113 classification algorithm.");
            
        /*
        //
        //  Read the data.
        */
        double[] a = Helpers.getExampleDoubleData();
            
        //
        //  Print a few data values.
        //
        Console.WriteLine("");
        Console.WriteLine("  First 5 data values:");
        Console.WriteLine("");

        for (int i = 1; i <= 5; i++)
        {
            string cout = "  " + i.ToString().PadLeft(8);
            for (int j = 1; j <= n; j++)
            {
                cout += "  " + a[i - 1 + (j - 1) * m].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);
        }

        //
        //  Assign points randomly to classes.
        //
        for (int i = 1; i <= m; i++)
        {
            c[i - 1] = i % k + 1;
        }

        //
        //  Define the critical value as the sum of the squares of the distances
        //  of the points to their cluster center.
        //
        for (int i = 1; i <= k; i++)
        {
            c_size[i - 1] = 0;
            for (int j = 1; j <= n; j++)
            {
                c_center[i - 1 + (j - 1) * k] = 0.0;
            }
        }

        for (int i = 1; i <= m; i++)
        {
            ci = c[i - 1];
            c_size[ci - 1] += 1;
            for (int j = 1; j <= n; j++)
            {
                c_center[ci - 1 + (j - 1) * k] += a[i - 1 + (j - 1) * m];
            }
        }

        for (int i = 1; i <= k; i++)
        {
            for (int j = 1; j <= n; j++)
            {
                c_center[i - 1 + (j - 1) * k] /= c_size[i - 1];
            }
        }

        for (int i = 1; i <= k; i++)
        {
            wss[i - 1] = 0.0;
        }

        for (int i = 1; i <= m; i++)
        {
            ci = c[i - 1];
            for (int j = 1; j <= n; j++)
            {
                wss[ci - 1] += Math.Pow(a[i - 1 + (j - 1) * m] - c_center[ci - 1 + (j - 1) * k], 2);
            }
        }

        double critvl = 0.0;
        for (int i = 1; i <= k; i++)
        {
            critvl += wss[i - 1];
        }

        Console.WriteLine("");
        Console.WriteLine("        Initial CRITVL = " + critvl + "");
        //
        //  Compute the clusters.
        //
        int ntrans1 = -1;
        int ntrans2 = -1;

        for (;;)
        {
            Algorithms.trnsfr(a, ref c, ref c_size, m, k, n, ref critvl, ref ntrans1, ref ifault);

            if (ntrans1 == 0 && ntrans2 == 0)
            {
                break;
            }

            Console.WriteLine("  After TRNSFR, CRITVL = " + critvl + "");

            Algorithms.swap(a, ref c, ref c_size, m, k, n, ref critvl, ref ntrans2, ref ifault);

            if (ntrans1 == 0 && ntrans2 == 0)
            {
                break;
            }

            Console.WriteLine("    After SWAP, CRITVL = " + critvl + "");
        }

        //
        //  Define the critical value as the sum of the squares of the distances
        //  of the points to their cluster center.
        //
        for (int i = 1; i <= k; i++)
        {
            for (int j = 1; j <= n; j++)
            {
                c_center[i - 1 + (j - 1) * k] = 0.0;
            }
        }

        for (int i = 1; i <= m; i++)
        {
            ci = c[i - 1];
            for (int j = 1; j <= n; j++)
            {
                c_center[ci - 1 + (j - 1) * k] += a[i - 1 + (j - 1) * m];
            }
        }

        for (int i = 1; i <= k; i++)
        {
            for (int j = 1; j <= n; j++)
            {
                c_center[i - 1 + (j - 1) * k] /= c_size[i - 1];
            }
        }

        for (int i = 1; i <= k; i++)
        {
            wss[i - 1] = 0.0;
        }

        for (int i = 1; i <= m; i++)
        {
            ci = c[i - 1];
            for (int j = 1; j <= n; j++)
            {
                wss[ci - 1] += Math.Pow(a[i - 1 + (j - 1) * m] - c_center[ci - 1 + (j - 1) * k], 2);
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  Cluster  Population  Energy");
        Console.WriteLine("");

        for (int i = 1; i <= k; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(8)
                                   + "  " + c_size[i - 1].ToString().PadLeft(8)
                                   + "  " + wss[i - 1].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("     Total"
                          + "  " + m.ToString().PadLeft(8)
                          + "  " + critvl.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

}