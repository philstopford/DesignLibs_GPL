﻿namespace Burkardt.CorrelationNS
{
    public static partial class Correlation
    {
        public static double[] correlation_constant ( int n, double[] rho, double rho0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CORRELATION_CONSTANT evaluates the constant correlation function.
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
                c[i] = 1.0;
            }
            return c;
        }
    }
}