using System;
using System.Globalization;
using Burkardt.PolynomialNS;

namespace PolPakTest;

public static class bernsteinTest
{
    public static void bernstein_poly_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNSTEIN_POLY_TEST tests BERNSTEIN_POLY.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 February 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double b = 0;
        double[] bvec = new double[11];
        int k = 0;
        int n = 0;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("BERNSTEIN_POLY_TEST:");
        Console.WriteLine("  BERNSTEIN_POLY evaluates the Bernstein polynomials.");
        Console.WriteLine("");
        Console.WriteLine("   N   K   X   Exact   B(N,K)(X)");
        Console.WriteLine("");

        int n_data = 0;

        for ( ; ; )
        {
            Burkardt.Values.Bernstein.bernstein_poly_values ( ref n_data, ref n, ref k, ref x, ref b );

            if ( n_data == 0 )
            {
                break;
            }

            BernsteinPolynomial.bernstein_poly ( n, x, ref bvec );

            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + k.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(7)
                                   + "  " + b.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + bvec[k].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

    }

    public static void bpab_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BPAB_TEST tests BPAB.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 10;

        double[] bern = new double[N+1];
        int i;

        Console.WriteLine("");
        Console.WriteLine("BPAB_TEST");
        Console.WriteLine("  BPAB evaluates Bernstein polynomials.");
        Console.WriteLine("");

        const double x = 0.3;
        const double a = 0.0;
        const double b = 1.0;

        BernsteinPolynomial.bpab ( N, x, a, b, ref bern );

        Console.WriteLine("  The Bernstein polynomials of degree " + N + "");
        Console.WriteLine("  based on the interval from " + a + "");
        Console.WriteLine("  to " + b + "");
        Console.WriteLine("  evaluated at X = " + x + "");
        Console.WriteLine("");

        for ( i = 0; i <= N; i++ )
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)      
                                   + "  " + bern[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

    }

}