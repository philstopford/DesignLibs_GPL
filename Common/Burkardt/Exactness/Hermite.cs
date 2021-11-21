﻿using System;
using System.Globalization;
using Burkardt.IntegralNS;
using Burkardt.Types;

namespace Burkardt.ExactnessNS;

public static partial class Exactness
{
    public static void hermite_exactness ( int n, double[] x, double[] w, int p_max )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_EXACTNESS investigates exactness of Hermite quadrature.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of points in the rule.
        //
        //    Input, double X[N], the quadrature points.
        //
        //    Input, double W[N], the quadrature weights.
        //
        //    Input, int P_MAX, the maximum exponent.
        //    0 <= P_MAX.
        //
    {
        int p;

        Console.WriteLine("");
        Console.WriteLine("  Quadrature rule for the Hermite integral.");
        Console.WriteLine("  Rule of order N = " + n + "");
        Console.WriteLine("  Degree          Relative Error");
        Console.WriteLine("");

        double[] v = new double[n];

        for ( p = 0; p <= p_max; p++ )
        {
            double s = Integral.hermite_integral ( p );

            int i;
            for ( i = 0; i < n; i++ )
            {
                v[i] = Math.Pow ( x[i], p );
            }

            double q = typeMethods.r8vec_dot_product ( n, w, v );

            double e = s switch
            {
                0.0 => Math.Abs(q),
                _ => Math.Abs(q - s) / Math.Abs(s)
            };

            Console.WriteLine(p.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                                      + e.ToString(CultureInfo.InvariantCulture).PadLeft(24) + "");
        }
    }
}