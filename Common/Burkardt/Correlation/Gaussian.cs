﻿using System;

namespace Burkardt.CorrelationNS
{
    public static partial class Correlation
    {
        public static double[] correlation_gaussian ( int n, double[] rho, double rho0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CORRELATION_GAUSSIAN evaluates the Gaussian correlation function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Petter Abrahamsen,
        //    A Review of Gaussian Random Fields and Correlation Functions,
        //    Norwegian Computing Center, 1997.
        //
        //  Parameters:
        //
        //    Input, int N, the number of arguments.
        //
        //    Input, double RHO[N], the arguments.
        //
        //    Input, double RHO0, the correlation length.
        //
        //    Output, double C[N], the correlations.
        //
        {
            double[] c;
            int i;

            c = new double[n];

            for ( i = 0; i < n; i++ )
            {
                c[i] = Math.Exp ( - Math.Pow ( rho[i] / rho0, 2 ) );
            }
            return c;
        }
    }
}