using System;
using Burkardt.Types;

namespace LineMonteCarloTest
{
    using MonteCarlo = Burkardt.Line.MonteCarlo;
    using Monomial = Burkardt.Monomial;

    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for LINE_MONTE_CARLO_TEST.
            //
            //  Discussion:
            //
            //    LINE_MONTE_CARLO_TEST tests the LINE_MONTE_CARLO library.
            //    
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 June 2017
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("LINE_MONTE_CARLO_TEST");
            Console.WriteLine("  Test the LINE_MONTE_CARLO library.");

            line01_sample_random_test();
            line01_sample_ergodic_test();

            Console.WriteLine("");
            Console.WriteLine("LINE_MONTE_CARLO_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void line01_sample_random_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LINE01_SAMPLE_RANDOM_TEST compares exact and estimated monomial integrals.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 June 2017
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int e;
            double error;
            double exact;
            int n = 4192;
            double result;
            int seed;
            int test;
            int test_num = 11;
            double[] value;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("LINE01_SAMPLE_RANDOM_TEST");
            Console.WriteLine("  LINE01_SAMPLE_RANDOM randomly samples the unit line segment.");
            Console.WriteLine("  Use it to estimate integrals.");
            //
            //  Get sample points.
            //
            seed = 123456789;
            x = MonteCarlo.line01_sample_random(n, ref seed);

            Console.WriteLine("");
            Console.WriteLine("  Number of sample points used is " + n + "");
            Console.WriteLine("");
            Console.WriteLine("   E     MC-Estimate      Exact           Error");
            Console.WriteLine("");

            for (test = 1; test <= test_num; test++)
            {
                e = test - 1;

                value = Monomial.monomial_value_1d(n, e, x);

                result = MonteCarlo.line01_length() * typeMethods.r8vec_sum(n, value) / (double) (n);
                exact = MonteCarlo.line01_monomial_integral(e);
                error = Math.Abs(result - exact);

                Console.WriteLine("  " + e.ToString().PadLeft(2)
                                       + "  " + result.ToString().PadLeft(14)
                                       + "  " + exact.ToString().PadLeft(14)
                                       + "  " + error.ToString().PadLeft(10) + "");
            }

        }

        static void line01_sample_ergodic_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LINE01_SAMPLE_ERGODIC_TEST compares exact and estimated monomial integrals.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 June 2017
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int e;
            double error;
            double exact;
            int n = 4192;
            double result;
            double shift;
            int test;
            int test_num = 11;
            double[] value;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("LINE01_SAMPLE_ERGODIC_TEST");
            Console.WriteLine("  LINE01_SAMPLE_ERGODIC ergodically samples the unit line segment.");
            Console.WriteLine("  Use it to estimate integrals.");
            //
            //  Get sample points.
            //
            shift = 0.0;
            x = MonteCarlo.line01_sample_ergodic(n, ref shift);

            Console.WriteLine("");
            Console.WriteLine("  Number of sample points used is " + n + "");
            Console.WriteLine("");
            Console.WriteLine("   E     MC-Estimate      Exact           Error");
            Console.WriteLine("");

            for (test = 1; test <= test_num; test++)
            {
                e = test - 1;

                value = Monomial.monomial_value_1d(n, e, x);

                result = MonteCarlo.line01_length() * typeMethods.r8vec_sum(n, value) / (double) (n);
                exact = MonteCarlo.line01_monomial_integral(e);
                error = Math.Abs(result - exact);

                Console.WriteLine("  " + e.ToString().PadLeft(2)
                                       + "  " + result.ToString().PadLeft(14)
                                       + "  " + exact.ToString().PadLeft(14)
                                       + "  " + error.ToString().PadLeft(6) + "");

            }

        }
    }
}