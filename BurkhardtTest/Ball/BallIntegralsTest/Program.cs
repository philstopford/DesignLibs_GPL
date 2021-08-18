using System;
using Burkardt.Ball;
using Burkardt.MonomialNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace BallIntegralsTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for BALL_INTEGRALS_TEST.
            //
            //  Discussion:
            //
            //    BALL_INTEGRALS_TEST tests the BALL_INTEGRALS library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 January 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("BALL_INTEGRALS_TEST");
            Console.WriteLine("  Test the BALL_INTEGRALS library.");

            test01();

            Console.WriteLine("");
            Console.WriteLine("BALL_INTEGRALS_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 uses BALL01_SAMPLE to compare estimated and exact integrals.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 January 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int[] e;
            double error;
            double exact;
            int i;
            int m = 3;
            int n = 4192;
            double result;
            int seed;
            int test;
            int test_num = 20;
            double[] value;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  Estimate monomial integrals using Monte Carlo");
            Console.WriteLine("  over the interior of the unit ball in 3D.");
            //
            //  Get sample points.
            //
            seed = 123456789;
            x = Integrals.ball01_sample(n, ref seed);

            Console.WriteLine("");
            Console.WriteLine("  Number of sample points used is " + n + "");
            //
            //  Randomly choose X,Y exponents between 0 and 8.
            //
            Console.WriteLine("");
            Console.WriteLine("  If any exponent is odd, the integral is zero.");
            Console.WriteLine("  We will restrict this test to randomly chosen even exponents.");
            Console.WriteLine("");
            Console.WriteLine("  Ex  Ey  Ez     MC-Estimate           Exact      Error");
            Console.WriteLine("");

            for (test = 1; test <= test_num; test++)
            {
                e = UniformRNG.i4vec_uniform_ab_new(m, 0, 4, ref seed);
                for (i = 0; i < m; i++)
                {
                    e[i] = e[i] * 2;
                }

                value = Monomial.monomial_value(m, n, e, x);

                result = Integrals.ball01_volume() * typeMethods.r8vec_sum(n, value)
                         / (double)(n);
                exact = Integrals.ball01_monomial_integral(e);
                error = Math.Abs(result - exact);

                Console.WriteLine("  " + e[0].ToString().PadLeft(2)
                                       + "  " + e[1].ToString().PadLeft(2)
                                       + "  " + e[2].ToString().PadLeft(2)
                                       + "  " + result.ToString().PadLeft(14)
                                       + "  " + exact.ToString().PadLeft(14)
                                       + "  " + error.ToString().PadLeft(10) + "");

            }
        }
    }
}