using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace ProbabilityTest
{
    partial class Program
    {
        static void i4_choose_test()

//****************************************************************************80
//
//  Purpose:
//
//    I4_CHOOSE_TEST tests I4_CHOOSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 October 2014
//
//  Author:
//
//    John Burkardt
//
        {
            int cnk;
            int k;
            int n;

            Console.WriteLine("");
            Console.WriteLine("I4_CHOOSE_TEST");
            Console.WriteLine("  I4_CHOOSE evaluates C(N,K).");
            Console.WriteLine("");
            Console.WriteLine("       N       K     CNK");

            for (n = 0; n <= 4; n++)
            {
                Console.WriteLine("");
                for (k = 0; k <= n; k++)
                {
                    cnk = typeMethods.i4_choose(n, k);

                    Console.WriteLine( "  "
                         + n.ToString().PadLeft(6) + "  "
                         + k.ToString().PadLeft(6) + "  "
                         + cnk.ToString().PadLeft(6) + "");
                }
            }

            return;
        }

        static void i4_choose_log_test()

//****************************************************************************80
//
//  Purpose:
//
//    I4_CHOOSE_LOG_TEST tests I4_CHOOSE_LOG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 March 2016
//
//  Author:
//
//    John Burkardt
//
        {
            int cnk;
            double elcnk;
            int k;
            double lcnk;
            int n;

            Console.WriteLine("");
            Console.WriteLine("I4_CHOOSE_LOG_TEST");
            Console.WriteLine("  I4_CHOOSE_LOG evaluates log(C(N,K)).");
            Console.WriteLine("");
            Console.WriteLine("       N       K            lCNK           elCNK     CNK");

            for (n = 0; n <= 4; n++)
            {
                Console.WriteLine("");
                for (k = 0; k <= n; k++)
                {
                    lcnk = typeMethods.i4_choose_log(n, k);
                    elcnk = Math.Exp(lcnk);
                    cnk = typeMethods.i4_choose(n, k);

                    Console.WriteLine("  " + n.ToString().PadLeft(6)
                                      + "  " + k.ToString().PadLeft(6)
                                      + "  " + lcnk.ToString().PadLeft(14)
                                      + "  " + elcnk.ToString().PadLeft(14)
                                      + "  " + cnk.ToString().PadLeft(6) + "");
                }
            }

            return;
        }

        static void i4_is_power_of_10_test()

//****************************************************************************80
//
//  Purpose:
//
//    I4_IS_POWER_OF_10_TEST tests I4_IS_POWER_OF_10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2016
//
//  Author:
//
//    John Burkardt
//
        {
            int i;

            Console.WriteLine("");
            Console.WriteLine("I4_IS_POWER_OF_10_TEST");
            Console.WriteLine("  I4_IS_POWER_OF_10 reports whether an I4 is a power of 10.");
            Console.WriteLine("");
            Console.WriteLine("  I     I4_IS_POWER_OF_10(I)");
            Console.WriteLine("");

            for (i = 97; i <= 103; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(6)
                                  + "  " + typeMethods.i4_is_power_of_10(i) + "");
            }

            return;
        }

        static void i4_uniform_ab_test()

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM_AB_TEST tests I4_UNIFORM_AB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 October 2014
//
//  Author:
//
//    John Burkardt
//
        {
            int a = -100;
            int b = 200;
            int i;
            int j;
            int seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("I4_UNIFORM_AB_TEST");
            Console.WriteLine("  I4_UNIFORM_AB computes pseudorandom values");
            Console.WriteLine("  in an interval [A,B].");

            Console.WriteLine("");
            Console.WriteLine("  The lower endpoint A = " + a + "");
            Console.WriteLine("  The upper endpoint B = " + b + "");
            Console.WriteLine("  The initial seed is " + seed + "");
            Console.WriteLine("");

            for (i = 1; i <= 20; i++)
            {
                j = UniformRNG.i4_uniform_ab(a, b, ref seed);

                Console.WriteLine("  " + i.ToString().PadLeft(8)
                                  + "  " + j.ToString().PadLeft(8) + "");
            }

            return;
        }

        static void i4vec_uniform_ab_new_test()

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_UNIFORM_AB_NEW_TEST tests I4VEC_UNIFORM_AB_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 October 2014
//
//  Author:
//
//    John Burkardt
//
        {
            int a = -100;
            int b = 200;
            int n = 20;
            int seed = 123456789;
            int[] v;

            Console.WriteLine("");
            Console.WriteLine("I4VEC_UNIFORM_AB_NEW_TEST");
            Console.WriteLine("  I4VEC_UNIFORM_AB_NEW computes pseudorandom values");
            Console.WriteLine("  in an interval [A,B].");

            Console.WriteLine("");
            Console.WriteLine("  The lower endpoint A = " + a + "");
            Console.WriteLine("  The upper endpoint B = " + b + "");
            Console.WriteLine("  The initial seed is " + seed + "");
            Console.WriteLine("");

            v = UniformRNG.i4vec_uniform_ab_new(n, a, b, ref seed);

            typeMethods.i4vec_print(n, v, "  The random vector:");

        }

        static void i4vec_unique_count_test()

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_UNIQUE_COUNT_TEST tests I4VEC_UNIQUE_COUNT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 March 2016
//
//  Author:
//
//    John Burkardt
//
        {
            int[] a;
            int a_hi;
            int a_lo;
            int a_unique;
            int n;
            int seed;

            Console.WriteLine("");
            Console.WriteLine("I4VEC_UNIQUE_COUNT_TEST");
            Console.WriteLine("  I4VEC_UNIQUE_COUNT counts unique entries in an I4VEC.");

            n = 20;
            a_lo = 0;
            a_hi = n;
            seed = 123456789;

            a = UniformRNG.i4vec_uniform_ab_new(n, a_lo, a_hi, ref seed);

            typeMethods.i4vec_print(n, a, "  Array:");

            a_unique = typeMethods.i4vec_unique_count(n, a);

            Console.WriteLine("");
            Console.WriteLine("  Number of unique entries is " + a_unique + "");
            
        }

    }
}