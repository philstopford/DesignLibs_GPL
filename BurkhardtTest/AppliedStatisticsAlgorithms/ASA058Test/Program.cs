using System;
using System.Globalization;
using Burkardt;
using Burkardt.AppliedStatistics;

namespace ASA058Test;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for ASA058_TEST.
        //
        //  Discussion:
        //
        //    ASA058_TEST tests the ASA058 library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("ASA058_TEST:");
        Console.WriteLine("  Test the ASA058 library.");

        test01 ( );

        Console.WriteLine("");
        Console.WriteLine("ASA058_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tries out the ASA058 routine.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int K = 5;
        const int M = 2;
        const int N = 100;

        int[] b = new int[N];
        double[] d = new double[K * M];
        double[] dev = new double[K];
        int[] e = new int[K];
        double[] f = new double[N];
        int i;
        int j;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Test the CLUSTR algorithm.");
        Console.WriteLine("  Applied Statistics Algorithm 58");
        /*
        //
        //  Read the data.
        */
        double[] x = Helpers.getExampleDoubleData();

        //
        //  Print a few data values.
        //
        Console.WriteLine("");
        Console.WriteLine("  First 5 data values:");
        Console.WriteLine("");

        for (i = 1; i <= 5; i++)
        {
            string cout = "  " + i.ToString().PadLeft(8);
            for (j = 1; j <= M; j++)
            {
                cout += "  " + x[i - 1 + (j - 1) * N].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);
        }

        //
        //  Initialize the cluster centers arbitrarily.
        //
        for (i = 1; i <= K; i++)
        {
            for (j = 1; j <= M; j++)
            {
                d[i - 1 + (j - 1) * K] = x[i - 1 + (j - 1) * N];
            }
        }

        //
        //  Compute the clusters.
        //
        int nz = 1;
        int k2 = K;

        Algorithms.clustr(x, ref d, ref dev, ref b, f, ref e, N, M, K, nz, k2);

        Console.WriteLine("");
        Console.WriteLine("  Cluster  Population  Energy");
        Console.WriteLine("");

        for (i = 1; i <= K; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(8)
                                   + "  " + e[i - 1].ToString().PadLeft(8)
                                   + "  " + dev[i - 1].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        int e_sum = 0;
        double dev_sum = 0.0;

        for (i = 1; i <= K; i++)
        {
            e_sum += e[i - 1];
            dev_sum += dev[i - 1];
        }

        Console.WriteLine("");
        Console.WriteLine("  " + "   Total"
                               + "  " + e_sum.ToString().PadLeft(8)
                               + "  " + dev_sum.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }


}