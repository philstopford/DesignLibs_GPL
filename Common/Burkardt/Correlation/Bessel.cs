using System;
using Burkardt.FullertonFnLib;

namespace Burkardt.CorrelationNS;

public static partial class Correlation
{
    public class CorrelationResult
    {
        public double[] result;
        public FullertonLib.BesselData data;
        public FullertonLib.r8BESJ0Data j0data;
        public FullertonLib.r8BESKData kdata;
        public FullertonLib.r8BESK1Data k1data;
    }

    public static CorrelationResult correlation_besselj(FullertonLib.BesselData globaldata, FullertonLib.r8BESJ0Data data, int n, double[] rho, double rho0 )

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
            c[i] = FullertonLib.r8_besj0(ref globaldata, ref data, rhohat);
        }

        CorrelationResult result = new() { result = c, data = globaldata, j0data = data };
            
        return result;
    }

    public static CorrelationResult correlation_besselk(FullertonLib.BesselData globaldata, FullertonLib.r8BESK1Data data, int n, double[] rho, double rho0 )

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
            switch (rho[i])
            {
                case 0.0:
                    c[i] = 1.0;
                    break;
                default:
                    rhohat = Math.Abs(rho[i]) / rho0;
                    c[i] = rhohat * FullertonLib.r8_besk1(ref globaldata, ref data, rhohat);
                    break;
            }
        }

        CorrelationResult res = new() { result = c, data = globaldata, k1data = data};
            
        return res;
    }
}