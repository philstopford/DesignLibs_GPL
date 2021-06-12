using System;
using Burkardt.Types;

namespace Burkardt.CorrelationNS
{
    public static partial class Correlation
    {
        public static double[] correlation_besselj(int n, double[] rho, double rho0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CORRELATION_BESSELJ evaluates the Bessel J correlation function.
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

            for (i = 0; i < n; i++)
            {
                rhohat = Math.Abs(rho[i]) / rho0;
                c[i] = r8_besj0(rhohat);
            }

            return c;
        }

        public static double[] correlation_besselk(int n, double[] rho, double rho0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CORRELATION_BESSELK evaluates the Bessel K correlation function.
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

            for (i = 0; i < n; i++)
            {
                if (rho[i] == 0.0)
                {
                    c[i] = 1.0;
                }
                else
                {
                    rhohat = Math.Abs(rho[i]) / rho0;
                    c[i] = rhohat * r8_besk1(rhohat);
                }
            }

            return c;
        }
    }
}