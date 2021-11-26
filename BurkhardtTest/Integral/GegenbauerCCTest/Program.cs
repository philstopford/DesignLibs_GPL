using System;
using System.Globalization;
using Burkardt;
using Burkardt.ChebyshevNS;
using Burkardt.Function;
using Burkardt.IntegralNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace GegenbauerCCTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for GEGENBAUER_CC_TEST.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("GEGENBAUER_CC_TEST:");
        Console.WriteLine("  Test the GEGENBAUER_CC library.");

        chebyshev_even1_test();
        chebyshev_even2_test();
        gegenbauer_cc1_test();
        gegenbauer_cc2_test();
        i4_uniform_ab_test();
        r8_mop_test();
        r8vec_print_test();
        r8vec2_print_test();

        Console.WriteLine("");
        Console.WriteLine("GEGENBAUER_CC_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static double f(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F is the function to be integrated.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument.
        //
        //    Output, double F, the value.
        //
    {
        const double a = 2.0;
        double value = Math.Cos(a * x);

        return value;
    }

    private static void chebyshev_even1_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBYSHEV_EVEN1_TEST tests CHEBYSHEV_EVEN1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a2_exact =  {
                0.4477815660,
                -0.7056685603,
                0.0680357987,
                -0.0048097159
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("CHEBYSHEV_EVEN1_TEST:");
        Console.WriteLine("  CHEBYSHEV_EVEN1 computes the even Chebyshev coefficients");
        Console.WriteLine("  of a function F, using the extreme points of Tn(x).");

        const int n = 6;

        double[] a2 = Chebyshev.chebyshev_even1(n, f);

        const int s = n / 2;
        typeMethods.r8vec2_print(s + 1, a2, a2_exact, "  Computed and Exact Coefficients:");
    }

    private static void chebyshev_even2_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBYSHEV_EVEN2_TEST tests CHEBYSHEV_EVEN2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("CHEBYSHEV_EVEN2_TEST:");
        Console.WriteLine("  CHEBYSHEV_EVEN2 computes the even Chebyshev coefficients");
        Console.WriteLine("  of a function F, using the zeros of Tn(x).");

        const int n = 6;

        double[] b2 = Chebyshev.chebyshev_even2(n, f);

        const int s = n / 2;
        typeMethods.r8vec_print(s + 1, b2, "  Computed Coefficients:");
    }

    private static void gegenbauer_cc1_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEGENBAUER_CC1_TEST tests GEGENBAUER_CC1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("GEGENBAUER_CC1_TEST:");
        Console.WriteLine("  GEGENBAUER_CC1 estimates the Gegenbauer integral of");
        Console.WriteLine("  a function f(x) using a Clenshaw-Curtis type approach");
        Console.WriteLine("  based on the extreme points of Tn(x).");

        const double lambda = 0.75;
        const double a = 2.0;
        const int n = 6;

        double value = Integral.gegenbauer_cc1(n, lambda, f);

        Console.WriteLine("");
        Console.WriteLine("  Value = " + value + "");
        double exact = Helpers.Gamma(lambda + 0.5) * Math.Sqrt(Math.PI) * BesselJ.besselj(lambda, a)
                       / Math.Pow(0.5 * a, lambda);
        Console.WriteLine("  Exact = " + exact + "");
    }

    private static void gegenbauer_cc2_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEGENBAUER_CC2_TEST tests GEGENBAUER_CC2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("GEGENBAUER_CC2_TEST:");
        Console.WriteLine("  GEGENBAUER_CC2 estimates the Gegenbauer integral of");
        Console.WriteLine("  a function f(x) using a Clenshaw-Curtis type approach");
        Console.WriteLine("  based on the zeros of Tn(x).");

        const double lambda = 0.75;
        const double a = 2.0;
        const int n = 6;

        double value = Integral.gegenbauer_cc2(n, lambda, f);

        Console.WriteLine("");
        Console.WriteLine("  Value = " + value + "");
        double exact = Helpers.Gamma(lambda + 0.5) * Math.Sqrt(Math.PI) * BesselJ.besselj(lambda, a)
                       / Math.Pow(0.5 * a, lambda);
        Console.WriteLine("  Exact = " + exact + "");
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

    private static void r8_mop_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_MOP_TEST tests R8_MOP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 December 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int seed = 123456789;
        int test;

        Console.WriteLine("");
        Console.WriteLine("R8_MOP_TEST");
        Console.WriteLine("  R8_MOP evaluates (-1.0)^I4 as an R8.");
        Console.WriteLine("");
        Console.WriteLine("    I4  R8_MOP(I4)");
        Console.WriteLine("");

        const int i4_min = -100;
        const int i4_max = +100;

        for (test = 1; test <= 10; test++)
        {
            int i4 = UniformRNG.i4_uniform_ab(i4_min, i4_max, ref seed);
            double r8 = typeMethods.r8_mop(i4);
            Console.WriteLine("  "
                              + i4.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                              + r8.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "");
        }
    }

    private static void r8vec_print_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_PRINT_TEST tests R8VEC_PRINT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a =  {
                123.456, 0.000005, -1.0E+06, 3.14159265
            }
            ;
        const int n = 4;

        Console.WriteLine("");
        Console.WriteLine("R8VEC_PRINT_TEST");
        Console.WriteLine("  R8VEC_PRINT prints an R8VEC.");

        typeMethods.r8vec_print(n, a, "  The R8VEC:");
    }

    private static void r8vec2_print_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC2_PRINT_TEST tests R8VEC2_PRINT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a =  {
                1.0, 2.0, 3.0, 4.0, 5.0
            }
            ;
        double[] b = new double[5];
        double[] c = new double[5];
        int i;
        const int n = 5;

        Console.WriteLine("");
        Console.WriteLine("R8VEC2_PRINT_TEST");
        Console.WriteLine("  R8VEC2_PRINT prints a pair of R8VEC's.");

        for (i = 0; i < n; i++)
        {
            b[i] = a[i] * a[i];
            c[i] = Math.Sqrt(a[i]);
        }

        typeMethods.r8vec2_print(n, b, c, "  Squares and square roots:");
    }
}