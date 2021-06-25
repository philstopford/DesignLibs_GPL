using System;
using Burkardt;
using Burkardt.Types;

using MonteCarlo = Burkardt.Wedge.MonteCarlo;

namespace WedgeMonteCarloTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for WEDGE_MONTE_CARLO_TEST.
            //
            //  Discussion:
            //
            //    WEDGE_MONTE_CARLO_TEST tests the WEDGE_MONTE_CARLO library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 August 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("WEDGE_MONTE_CARLO_TEST");
            Console.WriteLine("  Test the WEDGE_MONTE_CARLO library.");

            test01();
            //
            //  Terminate.
            //
            Console.WriteLine("");
            Console.WriteLine("WEDGE_MONTE_CARLO_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 uses WEDGE01_SAMPLE with an increasing number of points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 August 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int[] e = new int[3];
            int[] e_test =  {
                0, 0, 0,
                1, 0, 0,
                0, 1, 0,
                0, 0, 1,
                2, 0, 0,
                1, 1, 0,
                0, 0, 2,
                3, 0, 0
            }
            ;
            int i;
            int j;
            int m = 3;
            int n;
            double result;
            int seed;
            double[] value;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  Use WEDGE01_SAMPLE for a Monte Carlo estimate of an");
            Console.WriteLine("  integral over the interior of the unit wedge in 3D.");

            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("         N        1               X               Y " + 
                              "              Z                X^2            XY              Z^2    " + 
                              "        X^3");
            Console.WriteLine("");

            n = 1;

            while (n <= 65536)
            {
                x = MonteCarlo.wedge01_sample(n, ref seed);

                string cout = "  " + n.ToString().PadLeft(8);

                for (j = 0; j < 8; j++)
                {
                    for (i = 0; i < m; i++)
                    {
                        e[i] = e_test[i + j * m];
                    }

                    value = Monomial.monomial_value(m, n, e, x);

                    result = MonteCarlo.wedge01_volume() * typeMethods.r8vec_sum(n, value) / (double) (n);
                    cout += "  " + result.ToString().PadLeft(14);
                }

                Console.WriteLine(cout);
                
                n = 2 * n;
            }

            string cout2 = "     Exact";

            for (j = 0; j < 8; j++)
            {
                for (i = 0; i < m; i++)
                {
                    e[i] = e_test[i + j * m];
                }

                result = MonteCarlo.wedge01_integral(e);
                cout2 += "  " + result.ToString().PadLeft(14);
            }

            Console.WriteLine(cout2);
        }
    }
}