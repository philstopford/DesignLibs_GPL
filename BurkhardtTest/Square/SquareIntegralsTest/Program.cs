using System;
using Burkardt.MonomialNS;
using Burkardt.Square;
using Burkardt.Types;
using Burkardt.Uniform;

namespace SquareIntegralsTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SQUARE_INTEGRALS_TEST.
        //
        //  Discussion:
        //
        //    SQUARE_INTEGRALS_TEST tests the SQUARE_INTEGRALS library.
        //    
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 February 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("SQUARE_INTEGRALS_TEST");
            
        Console.WriteLine("  Test the SQUARE_INTEGRALS library.");

        square01_monomial_integral_test();
        squaresym_monomial_integral_test();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("SQUARE_INTEGRALS_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void square01_monomial_integral_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SQUARE01_MONOMIAL_INTEGRAL_TEST tests SQUARE01_MONOMIAL_INTEGRAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 February 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] e;
        double error;
        double exact;
        int m = 2;
        int n = 4192;
        double result;
        int seed;
        int test;
        int test_num = 20;
        double[] value;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("SQUARE01_MONOMIAL_INTEGRAL_TEST");
        Console.WriteLine("  SQUARE01_MONOMIAL_INTEGRAL returns the exact integral");
        Console.WriteLine("  of a monomial over the interior of the unit square in 2D.");
        Console.WriteLine("  Compare exact and estimated values.");
        //
        //  Get sample points.
        //
        seed = 123456789;
        x = Integrals.square01_sample(n, ref seed);
        Console.WriteLine("");
        Console.WriteLine("  Number of sample points is " + n + "");
        //
        //  Randomly choose exponents.
        //
        Console.WriteLine("");
        Console.WriteLine("  Ex  Ey     MC-Estimate           Exact      Error");
        Console.WriteLine("");

        for (test = 1; test <= test_num; test++)
        {
            e = UniformRNG.i4vec_uniform_ab_new(m, 0, 7, ref seed);

            value = Monomial.monomial_value(m, n, e, x);

            result = Integrals.square01_area() * typeMethods.r8vec_sum(n, value) / n;
            exact = Integrals.square01_monomial_integral(e);
            error = Math.Abs(result - exact);

            Console.WriteLine("  " + e[0].ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + e[1].ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + result.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + exact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + error.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

        }


    }

    private static void squaresym_monomial_integral_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SQUARESYM_MONOMIAL_INTEGRAL_TEST tests SQUARESYM_MONOMIAL_INTEGRAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 February 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] e;
        double error;
        double exact;
        int m = 2;
        int n = 4192;
        double result;
        int seed;
        int test;
        int test_num = 20;
        double[] value;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("SQUARESYM_MONOMIAL_INTEGRAL_TEST");
        Console.WriteLine("  SQUARESYM_MONOMIAL_INTEGRAL returns the exact integral");
        Console.WriteLine("  of a monomial over the interior of the symmetric unit square in 2D.");
        Console.WriteLine("  Compare exact and estimated values.");
        //
        //  Get sample points.
        //
        seed = 123456789;
        x = Integrals.squaresym_sample(n, ref seed);
        Console.WriteLine("");
        Console.WriteLine("  Number of sample points is " + n + "");
        //
        //  Randomly choose exponents.
        //
        Console.WriteLine("");
        Console.WriteLine("  Ex  Ey     MC-Estimate           Exact      Error");
        Console.WriteLine("");

        for (test = 1; test <= test_num; test++)
        {
            e = UniformRNG.i4vec_uniform_ab_new(m, 0, 7, ref seed);

            value = Monomial.monomial_value(m, n, e, x);

            result = Integrals.squaresym_area() * typeMethods.r8vec_sum(n, value) / n;
            exact = Integrals.squaresym_monomial_integral(e);
            error = Math.Abs(result - exact);

            Console.WriteLine("  " + e[0].ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + e[1].ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + result.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + exact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + error.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

        }

    }
}