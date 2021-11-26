using System;
using System.Globalization;
using Burkardt.Types;
using Burkardt.Uniform;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void i4_choose_test()

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
        int n;

        Console.WriteLine("");
        Console.WriteLine("I4_CHOOSE_TEST");
        Console.WriteLine("  I4_CHOOSE evaluates C(N,K).");
        Console.WriteLine("");
        Console.WriteLine("       N       K     CNK");

        for (n = 0; n <= 4; n++)
        {
            Console.WriteLine("");
            int k;
            for (k = 0; k <= n; k++)
            {
                int cnk = typeMethods.i4_choose(n, k);

                Console.WriteLine( "  "
                                   + n.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                   + k.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                   + cnk.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
            }
        }
    }

    private static void i4_choose_log_test()

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
        int n;

        Console.WriteLine("");
        Console.WriteLine("I4_CHOOSE_LOG_TEST");
        Console.WriteLine("  I4_CHOOSE_LOG evaluates log(C(N,K)).");
        Console.WriteLine("");
        Console.WriteLine("       N       K            lCNK           elCNK     CNK");

        for (n = 0; n <= 4; n++)
        {
            Console.WriteLine("");
            int k;
            for (k = 0; k <= n; k++)
            {
                double lcnk = typeMethods.i4_choose_log(n, k);
                double elcnk = Math.Exp(lcnk);
                int cnk = typeMethods.i4_choose(n, k);

                Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                       + "  " + k.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                       + "  " + lcnk.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + elcnk.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + cnk.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
            }
        }
    }

    private static void i4_is_power_of_10_test()

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
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + typeMethods.i4_is_power_of_10(i) + "");
        }
    }

    private static void i4_uniform_ab_test()

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
        const int a = -100;
        const int b = 200;
        int i;
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
            int j = UniformRNG.i4_uniform_ab(a, b, ref seed);

            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
        }
    }

    private static void i4vec_uniform_ab_new_test()

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
        const int a = -100;
        const int b = 200;
        const int n = 20;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("I4VEC_UNIFORM_AB_NEW_TEST");
        Console.WriteLine("  I4VEC_UNIFORM_AB_NEW computes pseudorandom values");
        Console.WriteLine("  in an interval [A,B].");

        Console.WriteLine("");
        Console.WriteLine("  The lower endpoint A = " + a + "");
        Console.WriteLine("  The upper endpoint B = " + b + "");
        Console.WriteLine("  The initial seed is " + seed + "");
        Console.WriteLine("");

        int[] v = UniformRNG.i4vec_uniform_ab_new(n, a, b, ref seed);

        typeMethods.i4vec_print(n, v, "  The random vector:");

    }

    private static void i4vec_unique_count_test()

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
        Console.WriteLine("");
        Console.WriteLine("I4VEC_UNIQUE_COUNT_TEST");
        Console.WriteLine("  I4VEC_UNIQUE_COUNT counts unique entries in an I4VEC.");

        const int n = 20;
        const int a_lo = 0;
        int seed = 123456789;

        int[] a = UniformRNG.i4vec_uniform_ab_new(n, a_lo, n, ref seed);

        typeMethods.i4vec_print(n, a, "  Array:");

        int a_unique = typeMethods.i4vec_unique_count(n, a);

        Console.WriteLine("");
        Console.WriteLine("  Number of unique entries is " + a_unique + "");
            
    }

}