﻿using System;
using Burkardt.Cube;
using Burkardt.Types;
using Burkardt.Uniform;

namespace CubeIntegralsTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for CUBE_INTEGRALS_TEST.
            //
            //  Discussion:
            //
            //    CUBE_INTEGRALS_TEST tests the CUBE_INTEGRALS library.
            //    
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 January 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("CUBE_INTEGRALS_TEST");
            Console.WriteLine("  C++ version");
            Console.WriteLine("  Test the CUBE_INTEGRALS library.");

            test01();

            Console.WriteLine("");
            Console.WriteLine("CUBE_INTEGRALS_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 estimates integrals over the unit cube in 3D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 January 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int[] e;
            double error;
            double exact;
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
            Console.WriteLine("  Compare exact and estimated integrals");
            Console.WriteLine("  over the interior of the unit cube in 3D.");
            //
            //  Get sample points.
            //
            seed = 123456789;
            x = Integrals.cube01_sample(n, ref seed);
            Console.WriteLine("");
            Console.WriteLine("  Number of sample points is " + n + "");
            //
            //  Randomly choose exponents.
            //
            Console.WriteLine("");
            Console.WriteLine("  Ex  Ey  Ez     MC-Estimate           Exact      Error");
            Console.WriteLine("");

            for (test = 1; test <= test_num; test++)
            {
                e = UniformRNG.i4vec_uniform_ab_new(m, 0, 7, ref seed);

                value = typeMethods.monomial_value(m, n, e, x);

                result = Integrals.cube01_volume() * typeMethods.r8vec_sum(n, value) / (double) (n);
                exact = Integrals.cube01_monomial_integral(e);
                error = Math.Abs(result - exact);

                Console.WriteLine("  " + e[0].ToString().PadLeft(2)
                    + "  " + e[1].ToString().PadLeft(2)
                    + "  " + e[2].ToString().PadLeft(2)
                    + "  " + result.ToString().PadLeft(14)
                    + "  " + exact.ToString().PadLeft(14)
                    + "  " + error.ToString().PadLeft(14) + "");
            }
        }
    }
}