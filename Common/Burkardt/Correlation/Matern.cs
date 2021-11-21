using System;
using Burkardt.FullertonFnLib;

namespace Burkardt.CorrelationNS;

public static partial class Correlation
{
    public static CorrelationResult correlation_matern (FullertonLib.BesselData globaldata, FullertonLib.r8BESKData data, int n, double[] rho, double rho0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CORRELATION_MATERN evaluates the Matern correlation function.
        //
        //  Discussion:
        //
        //    In order to call this routine under a dummy name, I had to drop NU from
        //    the parameter list.
        //
        //    The Matern correlation is
        //
        //      rho1 = 2 * sqrt ( nu ) * rho / rho0
        //
        //      c(rho) = ( rho1 )^nu * BesselK ( nu, rho1 ) 
        //               / gamma ( nu ) / 2 ^ ( nu - 1 )
        //
        //    The Matern covariance has the form:
        //
        //      K(rho) = sigma^2 * c(rho)
        //
        //    A Gaussian process with Matern covariance has sample paths that are
        //    differentiable (nu - 1) times.
        //
        //    When nu = 0.5, the Matern covariance is the exponential covariance.
        //
        //    As nu goes to +oo, the correlation converges to exp ( - (rho/rho0)^2 ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 November 2012
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
        //    Output, double C[N], the correlations.
        //
    {
        int i;

        double nu = 2.5;

        double[] c = new double[n];

        for ( i = 0; i < n; i++ )
        {
            double rho1 = 2.0 * Math.Sqrt ( nu ) * Math.Abs ( rho[i] ) / rho0;

            c[i] = rho1 switch
            {
                0.0 => 1.0,
                _ => Math.Pow(rho1, nu) * FullertonLib.r8_besk(ref globaldata, ref data, nu, rho1) / r8_gamma(nu) /
                     Math.Pow(2.0, nu - 1.0)
            };
        }
        return new CorrelationResult {result = c, data = globaldata, kdata = data};
    }
}