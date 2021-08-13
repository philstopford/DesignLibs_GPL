using System;
using Burkardt;
using Burkardt.MonomialNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace HyperballIntegralsTest
{
    using Integral = Burkardt.HyperGeometry.Hyperball.Integrals;
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for HYPERBALL_INTEGRALS_TEST.
            //
            //  Discussion:
            //
            //    HYPERBALL_INTEGRALS_TEST tests the HYPERBALL_INTEGRALS library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 January 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("HYPERBALL_INTEGRALS_TEST");
            Console.WriteLine("  Test the HYPERBALL_INTEGRALS library.");

            test01();
            test02();

            Console.WriteLine("");
            Console.WriteLine("HYPERBALL_INTEGRALS_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 uses HYPERBALL01_SAMPLE to compare exact and estimated integrals in 3D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 January 2014
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
            Console.WriteLine("  Use the Monte Carlo method to estimate integrals over");
            Console.WriteLine("  the interior of the unit hyperball in M dimensions.");
            Console.WriteLine("");
            Console.WriteLine("  Spatial dimension M = " + m + "");
            //
            //  Get sample points.
            //
            seed = 123456789;
            typeMethods.r8vecNormalData data = new typeMethods.r8vecNormalData();
            x = Integral.hyperball01_sample(m, n, ref data, ref seed);
            Console.WriteLine("");
            Console.WriteLine("  Number of sample points used is " + n + "");
            //
            //  Randomly choose exponents between 0 and 8.
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

                result = Integral.hyperball01_volume(m) * typeMethods.r8vec_sum(n, value)
                         / (double) (n);
                exact = Integral.hyperball01_monomial_integral(m, e);
                error = Math.Abs(result - exact);

                string cout = "";
                
                for (i = 0; i < m; i++)
                {
                    cout += "  " + e[i].ToString().PadLeft(2);
                }

                Console.WriteLine(cout + "  " + result.ToString().PadLeft(14)
                    + "  " + exact.ToString().PadLeft(14)
                    + "  " + error.ToString().PadLeft(10) + "");
            }
        }

        static void test02()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 uses HYPERBALL01_SAMPLE to compare exact and estimated integrals in 6D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 January 2014
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
            int m = 6;
            int n = 4192;
            double result;
            int seed;
            int test;
            int test_num = 20;
            double[] value;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  Use the Monte Carlo method to estimate integrals over");
            Console.WriteLine("  the interior of the unit hyperball in M dimensions.");
            Console.WriteLine("");
            Console.WriteLine("  Spatial dimension M = " + m + "");
            //
            //  Get sample points.
            //
            seed = 123456789;
            typeMethods.r8vecNormalData data = new typeMethods.r8vecNormalData();
            x = Integral.hyperball01_sample(m, n, ref data, ref seed);
            Console.WriteLine("");
            Console.WriteLine("  Number of sample points used is " + n + "");
            //
            //  Randomly choose exponents between 0 and 8.
            //
            Console.WriteLine("");
            Console.WriteLine("  If any exponent is odd, the integral is zero.");
            Console.WriteLine("  We will restrict this test to randomly chosen even exponents.");
            Console.WriteLine("");
            Console.WriteLine("  E1  E2  E3  E4  E5  E6     MC-Estimate           Exact      Error");
            Console.WriteLine("");   

            for (test = 1; test <= test_num; test++)
            {
                e = UniformRNG.i4vec_uniform_ab_new(m, 0, 4, ref seed);

                for (i = 0; i < m; i++)
                {
                    e[i] = e[i] * 2;
                }

                value = Monomial.monomial_value(m, n, e, x);

                result = Integral.hyperball01_volume(m) * typeMethods.r8vec_sum(n, value)
                         / (double) (n);
                exact = Integral.hyperball01_monomial_integral(m, e);
                error = Math.Abs(result - exact);

                string cout = "";
                
                for (i = 0; i < m; i++)
                {
                    cout += "  " + e[i].ToString().PadLeft(2);
                }

                Console.WriteLine(cout + "  " + result.ToString().PadLeft(14)
                    + "  " + exact.ToString().PadLeft(14)
                    + "  " + error.ToString().PadLeft(10) + "");
            }
        }
    }
}