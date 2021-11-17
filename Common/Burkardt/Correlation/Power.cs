using System;
using Burkardt.FullertonFnLib;

namespace Burkardt.CorrelationNS;

public static partial class Correlation
{
    public static CorrelationResult correlation_power (FullertonLib.BesselData globaldata, FullertonLib.r8BESK1Data data, int n, double[] rho, double rho0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CORRELATION_POWER evaluates the power correlation function.
        //
        //  Discussion:
        //
        //    In order to be able to call this routine under a dummy name, I had
        //    to drop E from the argument list.
        //
        //    The power correlation is
        //
        //      C(rho) = ( 1 - |rho| )^e  if 0 <= |rho| <= 1
        //             = 0                otherwise
        //
        //      The constraint on the exponent is 2 <= e.
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
        //    0.0 <= RHO.
        //
        //    Input, double RHO0, the correlation length.
        //    0.0 < RHO0.
        //
        //    Input, double E, the exponent.
        //    E has a default value of 2.0;
        //    2.0 <= E.
        //
        //    Output, double C[N], the correlations.
        //
    {
        double[] c;
        double e;
        int i;
        double rhohat;

        e = 2.0;

        c = new double[n];

        for ( i = 0; i < n; i++ )
        {
            rhohat = Math.Abs ( rho[i] ) / rho0;
            c[i] = rhohat switch
            {
                <= 1.0 => Math.Pow(1.0 - rhohat, e),
                _ => 0.0
            };
        }
        return new CorrelationResult(){result = c, data = globaldata};
    }
}