using System;
using Burkardt.AppliedStatistics;
using Burkardt.Table;
using Burkardt.Types;

namespace Burkardt.ASA053Test
{
    class Program
    {
        static void Main(string[] args)
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
        
        
        static void test01 ( )
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
            int NP = 3;

            double[] d = {
                3.0, 
                2.0, 4.0, 
                1.0, 2.0, 5.0 };
            int n;
            int np = NP;
            double[] sa;
            int seed;

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  Generate a single Wishart deviate.");

            n = 1;
            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("  The number of variables is " + np + "");
            Console.WriteLine("  The number of degrees of freedom is " + n + "");

            typeMethods.r8utp_print ( np, d, "  The upper Cholesky factor:" );

            sa = Algorithms.wshrt ( d, n, np, ref seed );

            typeMethods.r8pp_print ( np, sa, "  The sample matrix:" );
        }

        static void test02()
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
            int NP = 3;

            double[] d =  {
                3.0,
                2.0, 4.0,
                1.0, 2.0, 5.0
            }
            ;
            int np = NP;
            int test_num = 100000;

            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  Generate many Wishart deviates.");
            Console.WriteLine("  Compare to D' * D * np / n");
            int n = 2;
            int npp = (np * (np + 1)) / 2;
            int seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("  The number of variables is " + np + "");
            Console.WriteLine("  The number of degrees of freedom is " + n + "");

            typeMethods.r8utp_print(np, d, "  The upper Cholesky factor:");

            double[] s_average = new double[npp];

            for (int j = 0; j < npp; j++)
            {
                s_average[j] = 0.0;
            }

            for (int i = 1; i <= test_num; i++)
            {
                double[] sa = Algorithms.wshrt(d, n, np, ref seed);
                for (int j = 0; j < npp; j++)
                {
                    s_average[j] = s_average[j] + sa[j];
                }
            }

            for (int j = 0; j < npp; j++)
            {
                s_average[j] = s_average[j] / (double) (test_num);
            }

            typeMethods.r8pp_print(np, s_average, "  The averaged matrix:");
            //
            //  Compare the result to ( D' * D ) * np / n.
            //
            double[] sigma = new double[np * np];

            for (int i = 0; i < np; i++)
            {
                for (int j = 0; j < np; j++)
                {
                    sigma[i + j * np] = 0.0;
                    int k;
                    for (k = 0; k <= Math.Min(i, j); k++)
                    {
                        int ki = k + (i * (i + 1)) / 2;
                        int kj = k + (j * (j + 1)) / 2;
                        sigma[i + j * np] = sigma[i + j * np] + d[ki] * d[kj];
                    }

                    sigma[i + j * np] = sigma[i + j * np] * (double) np / (double) n;
                }
            }

            TableMisc.r8mat_print(np, np, sigma, "  Expected result:");

        }

    }
}