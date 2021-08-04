using System;
using Burkardt.FullertonFnLib;

namespace Burkardt.CorrelationNS
{
    public static partial class Correlation
    {
        public static CorrelationResult correlation_pentaspherical (FullertonLib.BesselData data,  int n, double[] rho, double rho0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CORRELATION_PENTASPHERICAL evaluates the pentaspherical correlation function.
        //
        //  Discussion:
        //
        //    This correlation is based on the volume of overlap of two spheres
        //    of radius RHO0 and separation RHO.
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
            double rhohat;

            c = new double[n];

            for ( i = 0; i < n; i++ )
            {
                rhohat = Math.Min ( Math.Abs ( rho[i] ) / rho0, 1.0 );

                c[i] = 1.0 - 1.875 * rhohat + 1.25 * Math.Pow ( rhohat, 3 )
                       - 0.375 * Math.Pow ( rhohat, 5 );
            }

            return new CorrelationResult(){result = c, data = data};
        }
    }
}