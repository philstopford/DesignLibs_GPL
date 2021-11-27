using System;
using Burkardt.AppliedStatistics;
using Burkardt.Types;

namespace ASA053Test;

internal static class Program
{
    private static void Main()
        //
        //  Purpose:
        //
        //    MAIN is the main program for ASA053_TEST.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("ASA053_TEST:");
        Console.WriteLine("  Test the ASA053 library.");

        test01 ( );
        test02 ( );

        Console.WriteLine("");
        Console.WriteLine("ASA053_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }


    private static void test01 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 generates a random Wishart variate.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int NP = 3;

        double[] d = {
            3.0, 
            2.0, 4.0, 
            1.0, 2.0, 5.0 };

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Generate a single Wishart deviate.");

        const int n = 1;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("  The number of variables is " + NP + "");
        Console.WriteLine("  The number of degrees of freedom is " + n + "");

        typeMethods.r8utp_print ( NP, d, "  The upper Cholesky factor:" );

        double[] sa = Algorithms.wshrt ( d, n, NP, ref seed );

        typeMethods.r8pp_print ( NP, sa, "  The sample matrix:" );
    }

    private static void test02()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 averages many Wishart samples.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int NP = 3;

        double[] d =  {
                3.0,
                2.0, 4.0,
                1.0, 2.0, 5.0
            }
            ;
        const int test_num = 100000;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Generate many Wishart deviates.");
        Console.WriteLine("  Compare to D' * D * np / n");
        const int n = 2;
        const int npp = NP * (NP + 1) / 2;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("  The number of variables is " + NP + "");
        Console.WriteLine("  The number of degrees of freedom is " + n + "");

        typeMethods.r8utp_print(NP, d, "  The upper Cholesky factor:");

        double[] s_average = new double[npp];

        for (int j = 0; j < npp; j++)
        {
            s_average[j] = 0.0;
        }

        for (int i = 1; i <= test_num; i++)
        {
            double[] sa = Algorithms.wshrt(d, n, NP, ref seed);
            for (int j = 0; j < npp; j++)
            {
                s_average[j] += sa[j];
            }
        }

        for (int j = 0; j < npp; j++)
        {
            s_average[j] /= test_num;
        }

        typeMethods.r8pp_print(NP, s_average, "  The averaged matrix:");
        //
        //  Compare the result to ( D' * D ) * np / n.
        //
        double[] sigma = new double[NP * NP];

        for (int i = 0; i < NP; i++)
        {
            for (int j = 0; j < NP; j++)
            {
                sigma[i + j * NP] = 0.0;
                int k;
                for (k = 0; k <= Math.Min(i, j); k++)
                {
                    int ki = k + i * (i + 1) / 2;
                    int kj = k + j * (j + 1) / 2;
                    sigma[i + j * NP] += d[ki] * d[kj];
                }

                sigma[i + j * NP] = sigma[i + j * NP] * NP / n;
            }
        }

        typeMethods.r8mat_print(NP, NP, sigma, "  Expected result:");

    }

}