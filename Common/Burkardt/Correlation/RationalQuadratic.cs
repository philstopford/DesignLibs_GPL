using System;
using Burkardt.FullertonFnLib;

namespace Burkardt.CorrelationNS;

public static partial class Correlation
{
    public static CorrelationResult correlation_rational_quadratic (FullertonLib.BesselData globaldata, FullertonLib.r8BESJ0Data data, int n, double[] rho, double rho0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CORRELATION_RATIONAL_QUADRATIC: rational quadratic correlation function.
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
        int i;

        double[] c = new double[n];

        for ( i = 0; i < n; i++ )
        {
            c[i] = 1.0 / ( 1.0 + Math.Pow ( rho[i] / rho0, 2 ) );
        }

        return new CorrelationResult {result = c, data = globaldata, j0data = data};
    }
}