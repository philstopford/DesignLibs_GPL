﻿using System;
using System.Globalization;
using Burkardt.LineNS;
using Burkardt.Types;
using Monomial = Burkardt.MonomialNS.Monomial;

namespace LineIntegralsTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for LINE_INTEGRALS_TEST.
        //
        //  Discussion:
        //
        //    LINE_INTEGRALS_TEST tests the LINE_INTEGRALS library.
        //    
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("LINE_INTEGRALS_TEST");
        Console.WriteLine("  Test the LINE_INTEGRALS library.");

        test01();

        Console.WriteLine("");
        Console.WriteLine("LINE_INTEGRALS_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 compares exact and estimated monomial integrals.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 4192;
        int test;
        const int test_num = 11;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Compare exact and estimated integrals ");
        Console.WriteLine("  over the length of the unit line in 1D.");
        //
        //  Get sample points.
        //
        int seed = 123456789;
        double[] x = Integrals.line01_sample(n, ref seed);

        Console.WriteLine("");
        Console.WriteLine("  Number of sample points used is " + n + "");
        Console.WriteLine("");
        Console.WriteLine("   E     MC-Estimate      Exact           Error");
        Console.WriteLine("");

        for (test = 1; test <= test_num; test++)
        {
            int e = test - 1;

            double[] value = Monomial.monomial_value_1d(n, e, x);

            double result = Integrals.line01_length() * typeMethods.r8vec_sum(n, value) / n;
            double exact = Integrals.line01_monomial_integral(e);
            double error = Math.Abs(result - exact);

            Console.WriteLine("  " + e.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + result.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + exact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + error.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");

        }
    }
}