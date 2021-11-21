using System;
using Burkardt.FullertonFnLib;

namespace Burkardt.CorrelationNS;

public static partial class Correlation
{
    public static CorrelationResult correlation_hole (FullertonLib.BesselData globaldata, FullertonLib.r8BESK1Data data, int n, double[] rho, double rho0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CORRELATION_HOLE evaluates the hole correlation function.
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
            c[i] = ( 1.0 - Math.Abs ( rho[i] ) / rho0 ) 
                   * Math.Exp ( - Math.Abs ( rho[i] ) / rho0 );
        }
        return new CorrelationResult {result = c, data = globaldata, k1data = data};
    }
}