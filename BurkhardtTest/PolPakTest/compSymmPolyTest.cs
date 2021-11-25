﻿using System;
using System.Globalization;
using Burkardt.PolynomialNS;
using Burkardt.Types;

namespace PolPakTest;

public static class compSymmPolyTest
{
    public static void complete_symmetric_poly_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COMPLETE_SYMMETRIC_POLY_TEST tests COMPLETE_SYMMETRIC_POLY.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 November 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 5;
        int nn;
        double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0 };

        Console.WriteLine("");
        Console.WriteLine("COMPLETE_SYMMETRIC_POLY_TEST");
        Console.WriteLine("  COMPLETE_SYMMETRIC_POLY evaluates a complete symmetric.");
        Console.WriteLine("  polynomial in a given set of variables X.");

        typeMethods.r8vec_print ( n, x, "  Variable vector X:" );

        Console.WriteLine("");
        Console.WriteLine("   N\\R     0       1       2       3       4       5");
        Console.WriteLine("");

        for ( nn = 0; nn <= n; nn++ )
        {
            string cout = "  " + nn.ToString().PadLeft(2);
            int rr;
            for ( rr = 0; rr <= 5; rr++ )
            {
                double value = CompleteSymmetric.complete_symmetric_poly ( nn, rr, x );
                cout += "  " + value.ToString(CultureInfo.InvariantCulture).PadLeft(6);
            }
            Console.WriteLine(cout);
        }

    }

}