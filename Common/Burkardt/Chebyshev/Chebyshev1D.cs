using System;
using Burkardt.SolveNS;
using Burkardt.Types;

namespace Burkardt.ChebyshevNS
{
    public static class Chebyshev1D
    {
        public static double[] chebyshev_coef_1d(int nd, double[] xd, double[] yd, ref double xmin,
        ref double xmax )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBYSHEV_COEF_1D determines the Chebyshev interpolant coefficients.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ND, the number of data points.
        //    ND must be at least 1.
        //
        //    Input, double XD[ND], the data locations.
        //
        //    Input, double YD[ND], the data values.
        //
        //    Output, double &XMIN, &XMAX, the interpolation interval.
        //
        //    Output, double CHEBYSHEV_COEF_1D[ND], the Chebyshev coefficients.
        //
        {
            double[] a;
            double[] c;
            int i;
            int j;
            double[] x;

            if (nd == 1)
            {
                xmin = xd[0];
                xmax = xd[0];
                c = new double[nd];
                c[0] = 1.0;
                return c;
            }

            xmin = typeMethods.r8vec_min(nd, xd);
            xmax = typeMethods.r8vec_max(nd, xd);
            //
            //  Map XD to [-1,+1].
            //
            x = new double[nd];
            for (i = 0; i < nd; i++)
            {
                x[i] = (2.0 * xd[i] - xmin - xmax) / (xmax - xmin);
            }

            //
            //  Form the Chebyshev Vandermonde matrix.
            //
            a = new double[nd * nd];
            for (j = 0; j < nd; j++)
            {
                for (i = 0; i < nd; i++)
                {
                    a[i + j * nd] = Math.Cos(Math.Acos(x[i]) * (double) (j));
                }
            }

            //
            //  Solve for the expansion coefficients.
            //
            c = QRSolve.qr_solve(nd, nd, a, yd);

            return c;
        }

        public static double[] chebyshev_interp_1d(int nd, double[] xd, double[] yd, int ni,
        double[] xi )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBYSHEV_INTERP_1D determines and evaluates the Chebyshev interpolant.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ND, the number of data points.
        //    ND must be at least 1.
        //
        //    Input, double XD[ND], the data locations.
        //
        //    Input, double YD[ND], the data values.
        //
        //    Input, int NI, the number of interpolation points.
        //
        //    Input, double XI[NI], the interpolation points, which
        //    must be each be in the interval [ min(XD), max(XD)].
        //
        //    Output, double YI[NI], the interpolated values.
        //
        {
            double xmax = 0;
            double xmin = 0;

            double[] c = chebyshev_coef_1d(nd, xd, yd, ref xmin, ref xmax);

            double[] yi = chebyshev_value_1d(nd, c, xmin, xmax, ni, xi);
            
            return yi;
        }

        public static double[] chebyshev_value_1d(int nd, double[] c, double xmin, double xmax,
        int ni, double[] xi )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBYSHEV_VALUE_1D evaluates a Chebyshev interpolant, given its coefficients.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ND, the number of data points.
        //    ND must be at least 1.
        //
        //    Input, double C[ND], the Chebyshev coefficients.
        //
        //    Input, double XMIN, XMAX, the interpolation interval.
        //
        //    Input, int NI, the number of interpolation points.
        //
        //    Input, double XI[NI], the interpolation points, which
        //    must be each be in the interval [XMIN,XMAX].
        //
        //    Output, double YI[NI], the interpolated values.
        //
        {
            double[] yi;

            if (nd == 1)
            {
                yi = new double[nd];
                yi[0] = c[0];
                return yi;
            }

            //
            //  Map XI to [-1,+1].
            //
            double[] x = new double[ni];
            for (int i = 0; i < ni; i++)
            {
                x[i] = (2.0 * xi[i] - xmin - xmax) / (xmax - xmin);
            }

            double[] a = new double[ni * nd];
            for (int j = 0; j < nd; j++)
            {
                for (int i = 0; i < ni; i++)
                {
                    a[i + j * ni] = Math.Cos(Math.Acos(x[i]) * (double) (j));
                }
            }

            yi = typeMethods.r8mat_mv_new(ni, nd, a, c);

            return yi;
        }
    }
}