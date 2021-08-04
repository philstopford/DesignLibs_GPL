using System;
using Burkardt.FullertonFnLib;

namespace Burkardt.CorrelationNS
{
    public static partial class Correlation
    {
        public static CorrelationResult correlation_matern (FullertonLib.BesselData data, int n, double[] rho, double rho0 )

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
            double[] c;
            int i;
            double nu;
            double rho1;

            nu = 2.5;

            c = new double[n];

            for ( i = 0; i < n; i++ )
            {
                rho1 = 2.0 * Math.Sqrt ( nu ) * Math.Abs ( rho[i] ) / rho0;

                if ( rho1 == 0.0 )
                {
                    c[i] = 1.0;
                }
                else
                {
                    c[i] = Math.Pow ( rho1, nu ) * FullertonLib.r8_besk (ref data, nu, rho1 ) / r8_gamma ( nu ) 
                                                                   / Math.Pow ( 2.0, nu - 1.0 );
                }
            }
            return new CorrelationResult(){result = c, data = data};
        }
    }
}