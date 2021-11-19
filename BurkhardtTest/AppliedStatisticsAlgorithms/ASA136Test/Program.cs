﻿using System;
using System.Globalization;
using Burkardt;
using Burkardt.AppliedStatistics;

namespace ASA136Test;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for ASA136_TEST.
        //
        //  Discussion:
        //
        //    ASA136_TEST tests the ASA136 library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("ASA136_TEST:");
        Console.WriteLine("  Test the ASA136 library.");

        test01();

        Console.WriteLine("");
        Console.WriteLine("ASA136_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tries out the ASA136 routine.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int ifault = 0;
        const int k = 5;
        const int m = 100;
        const int n = 2;

        double[] c = new double[k * n];
        int[] ic1 = new int[m];
        int[] nc = new int[k];
        double[] wss = new double[k];

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Test the KMNS algorithm,");
        Console.WriteLine("  Applied Statistics Algorithm #136.");

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
                cout += a[i - 1 + (j - 1) * m].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);
        }

        //
        //  Initialize the cluster centers.
        //  Here, we arbitrarily make the first K data points cluster centers.
        //
        for (int i = 1; i <= k; i++)
        {
            for (int j = 1; j <= n; j++)
            {
                c[i - 1 + (j - 1) * k] = a[i - 1 + (j - 1) * m];
            }
        }

        const int iter = 50;
        //
        //  Compute the clusters.
        //
        Algorithms.kmns(a, m, n, ref c, k, ref ic1, ref nc, iter, ref wss, ref ifault);

        if (ifault != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("TEST01 - Fatal error!");
            Console.WriteLine("  KMNS returned IFAULT = " + ifault + "");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("  Cluster  Population  Energy");
        Console.WriteLine("");

        int nc_sum = 0;
        double wss_sum = 0.0;

        for (int i = 1; i <= k; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(8)
                                   + "  " + nc[i - 1].ToString().PadLeft(8)
                                   + "  " + wss[i - 1].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            nc_sum += nc[i - 1];
            wss_sum += wss[i - 1];
        }

        Console.WriteLine("");
        Console.WriteLine("     Total"
                          + "  " + nc_sum.ToString().PadLeft(8)
                          + "  " + wss_sum.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

}