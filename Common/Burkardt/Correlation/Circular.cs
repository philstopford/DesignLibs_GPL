using System;
using Burkardt.FullertonFnLib;

namespace Burkardt.CorrelationNS
{
    public static partial class Correlation
    {
        public static CorrelationResult correlation_circular (FullertonLib.BesselData globaldata, FullertonLib.r8BESK1Data data,  int n, double[] rho, double rho0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CORRELATION_CIRCULAR evaluates the circular correlation function.
        //
        //  Discussion:
        //
        //    This correlation is based on the area of overlap of two circles
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
            double pi = 3.141592653589793;
            double rhohat;

            c = new double[n];

            for ( i = 0; i < n; i++ )
            {
                rhohat = Math.Min ( Math.Abs ( rho[i] ) / rho0, 1.0 );

                c[i] = ( 1.0 - ( 2.0 / pi ) 
                    * ( rhohat * Math.Sqrt ( 1.0 - rhohat * rhohat ) + Math.Asin ( rhohat ) ) );
            }

            return new CorrelationResult() {result = c, data = globaldata, k1data = data};
        }
    }
}