using System;
using System.Globalization;
using Burkardt.Types;

namespace Burkardt.IntegralNS;

public static partial class Integral
{
    public static void legendre_2d_exactness ( double[] a, double[] b, int n, ref double[] x, 
            ref double[] y, ref double[] w, int t )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGENDRE_2D_EXACTNESS: monomial exactness for the 2D Legendre integral.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A[2], the lower limits of integration.
        //
        //    Input, double B[2], the upper limits of integration.
        //
        //    Input, int N, the number of points in the rule.
        //
        //    Input, double X[N], Y[N], the quadrature points.
        //
        //    Input, double W[N], the quadrature weights.
        //
        //    Input, int T, the maximum total degree.
        //    0 <= T.
        //
    {
        int[] p = new int[2];
        int tt;

        Console.WriteLine("");
        Console.WriteLine("  Quadrature rule for the 2D Legendre integral.");
        Console.WriteLine("  Number of points in rule is " + n + "");
        Console.WriteLine("");
        Console.WriteLine("   D   I       J          Relative Error");

        double[] v = new double[n];

        for ( tt = 0; tt <= t; tt++ )
        {
            Console.WriteLine("  " + tt + "");
            int j;
            for ( j = 0; j <= tt; j++ )
            {
                int i = tt - j;

                p[0] = i;
                p[1] = j;

                double s = legendre_2d_monomial_integral ( a, b, p );

                for ( i = 0; i < n; i++ )
                {
                    v[i] = Math.Pow( x[i], p[0] ) * Math.Pow( y[i], p[1] );
                }
                double q = typeMethods.r8vec_dot_product ( n, w, v );

                double e = s switch
                {
                    0.0 => Math.Abs(q),
                    _ => Math.Abs(q - s) / Math.Abs(s)
                };
                Console.WriteLine(p[0].ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                                             + p[1].ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                                             + e.ToString("0.################").PadLeft(24) + "");
            }
        }
    }
        
    public static double legendre_2d_monomial_integral ( double[] a, double[] b, int[] p )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGENDRE_2D_MONOMIAL_INTEGRAL the Legendre integral of a monomial.
        //
        //  Discussion:
        //
        //    The Legendre integral to be evaluated has the form
        //
        //      I(f) = integral ( y1 <= y <= y2 ) 
        //             integral ( x1 <= x <= x2 ) x^i y^j dx dy
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A[2], the lower limits of integration.
        //
        //    Input, double B[2], the upper limits of integration.
        //
        //    Input, int P[2], the exponents of X and Y.
        //
        //    Output, double LEGENDRE_2D_MONOMIAL_INTEGRAL, the value of the 
        //    exact integral.
        //
    {
        double exact = ( Math.Pow( b[0], p[0] + 1 ) - Math.Pow( a[0], p[0] + 1 ) ) 
                       / (p[0] + 1) 
                       * ( Math.Pow( b[1], p[1] + 1 ) - Math.Pow( a[1], p[1] + 1 ) ) 
                       / (p[1] + 1);

        return exact;
    }
    public static double legendre_integral ( int p )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGENDRE_INTEGRAL evaluates a monomial Legendre integral.
        //
        //  Discussion:
        //
        //    Integral ( -1 <= x <= +1 ) x^p dx
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int P, the exponent.
        //    0 <= P.
        //
        //    Output, double LEGENDRE_INTEGRAL, the value of the exact integral.
        //
    {
        double s = (p % 2) switch
        {
            0 => 2.0 / (p + 1),
            _ => 0.0
        };

        return s;
    }
}