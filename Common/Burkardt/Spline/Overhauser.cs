using System;
using Burkardt.Function;
using Burkardt.Types;

namespace Burkardt.Spline
{
    public static class Overhauser
    {
        public static double spline_overhauser_nonuni_val(int ndata, double[] tdata,
                double[] ydata, double tval)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPLINE_OVERHAUSER_NONUNI_VAL evaluates the nonuniform Overhauser spline.
            //
            //  Discussion:
            //
            //    The nonuniformity refers to the fact that the abscissas values
            //    need not be uniformly spaced.
            //
            //    Thanks to Doug Fortune for pointing out that the point distances
            //    used to define ALPHA and BETA should be the Euclidean distances
            //    between the points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 May 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NDATA, the number of data points.
            //    NDATA must be at least 3.
            //
            //    Input, double TDATA[NDATA], the abscissas of the data points.
            //    The values of TDATA are assumed to be distinct and increasing.
            //
            //    Input, double YDATA[NDATA], the data values.
            //
            //    Input, double TVAL, the value where the spline is to
            //    be evaluated.
            //
            //    Output, double SPLINE_OVERHAUSER_NONUNI_VAL, the value of the 
            //    spline at TVAL.
            //
        {
            double alpha;
            double beta;
            double d21;
            double d32;
            double d43;
            int left = 0;
            double[] mbasis;
            int right = 0;
            double yval;
            //
            //  Check NDATA.
            //
            if (ndata < 3)
            {
                Console.WriteLine("");
                Console.WriteLine("SPLINE_OVERHAUSER_NONUNI_VAL - Fatal error!");
                Console.WriteLine("  NDATA < 3.");
                return (1);
            }

            //
            //  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
            //
            typeMethods.r8vec_bracket(ndata, tdata, tval, ref left, ref right);
            //
            //  Evaluate the spline in the given interval.
            //
            if (left == 1)
            {
                d21 = Math.Sqrt(Math.Pow(tdata[1] - tdata[0], 2)
                                + Math.Pow(ydata[1] - ydata[0], 2));

                d32 = Math.Sqrt(Math.Pow(tdata[2] - tdata[1], 2)
                                + Math.Pow(ydata[2] - ydata[1], 2));

                alpha = d21 / (d32 + d21);

                mbasis = Basis.basis_matrix_overhauser_nul(alpha);

                yval = Basis.basis_matrix_tmp(left, 3, mbasis, ndata, tdata, ydata, tval);
            }
            else if (left < ndata - 1)
            {
                d21 = Math.Sqrt(Math.Pow(tdata[left - 1] - tdata[left - 2], 2)
                                + Math.Pow(ydata[left - 1] - ydata[left - 2], 2));

                d32 = Math.Sqrt(Math.Pow(tdata[left] - tdata[left - 1], 2)
                                + Math.Pow(ydata[left] - ydata[left - 1], 2));

                d43 = Math.Sqrt(Math.Pow(tdata[left + 1] - tdata[left], 2)
                                + Math.Pow(ydata[left + 1] - ydata[left], 2));

                alpha = d21 / (d32 + d21);
                beta = d32 / (d43 + d32);

                mbasis = Basis.basis_matrix_overhauser_nonuni(alpha, beta);

                yval = Basis.basis_matrix_tmp(left, 4, mbasis, ndata, tdata, ydata, tval);
            }
            else if (left == ndata - 1)
            {
                d32 = Math.Sqrt(Math.Pow(tdata[ndata - 2] - tdata[ndata - 3], 2)
                                + Math.Pow(ydata[ndata - 2] - ydata[ndata - 3], 2));

                d43 = Math.Sqrt(Math.Pow(tdata[ndata - 1] - tdata[ndata - 2], 2)
                                + Math.Pow(ydata[ndata - 1] - ydata[ndata - 2], 2));

                beta = d32 / (d43 + d32);

                mbasis = Basis.basis_matrix_overhauser_nur(beta);

                yval = Basis.basis_matrix_tmp(left, 3, mbasis, ndata, tdata, ydata, tval);
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("SPLINE_OVERHAUSER_NONUNI_VAL - Fatal error!");
                Console.WriteLine("  Nonsensical value of LEFT = " + left + "");
                Console.WriteLine("  but 0 < LEFT < NDATA = " + ndata + "");
                Console.WriteLine("  is required.");
                return (1);
            }

            return yval;
        }

        public static double spline_overhauser_uni_val(int ndata, double[] tdata, double[] ydata,
                double tval)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPLINE_OVERHAUSER_UNI_VAL evaluates the uniform Overhauser spline.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 February 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NDATA, the number of data points.
            //    NDATA must be at least 3.
            //
            //    Input, double TDATA[NDATA], the abscissas of the data points.
            //    The values of TDATA are assumed to be distinct and increasing.
            //    This routine also assumes that the values of TDATA are uniformly
            //    spaced; for instance, TDATA(1) = 10, TDATA(2) = 11, TDATA(3) = 12...
            //
            //    Input, double YDATA[NDATA], the data values.
            //
            //    Input, double TVAL, the value where the spline is to
            //    be evaluated.
            //
            //    Output, double SPLINE_OVERHAUSER_UNI_VAL, the value of the spline at TVAL.
            //
        {
            int left = 0;
            double[] mbasis;
            int right = 0;
            double yval = 0;
            //
            //  Check NDATA.
            //
            if (ndata < 3)
            {
                Console.WriteLine("");
                Console.WriteLine("SPLINE_OVERHAUSER_UNI_VAL - Fatal error!");
                Console.WriteLine("  NDATA < 3.");
                return (1);
            }

            //
            //  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
            //
            typeMethods.r8vec_bracket(ndata, tdata, tval, ref left, ref right);
            //
            //  Evaluate the spline in the given interval.
            //
            if (left == 1)
            {

                mbasis = Basis.basis_matrix_overhauser_uni_l();

                yval = Basis.basis_matrix_tmp(left, 3, mbasis, ndata, tdata, ydata, tval);
            }
            else if (left < ndata - 1)
            {
                mbasis = Basis.basis_matrix_overhauser_uni();

                yval = Basis.basis_matrix_tmp(left, 4, mbasis, ndata, tdata, ydata, tval);
            }
            else if (left == ndata - 1)
            {
                mbasis = Basis.basis_matrix_overhauser_uni_r();

                yval = Basis.basis_matrix_tmp(left, 3, mbasis, ndata, tdata, ydata, tval);

            }

            return yval;
        }

        public static void spline_overhauser_val(int ndim, int ndata, double[] tdata,
                double[] ydata, double tval, double[] yval)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPLINE_OVERHAUSER_VAL evaluates an Overhauser spline.
            //
            //  Discussion:
            //
            //    Over the first and last intervals, the Overhauser spline is a 
            //    quadratic.  In the intermediate intervals, it is a piecewise cubic.
            //    The Overhauser spline is also known as the Catmull-Rom spline.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    11 December 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    JA Brewer, DC Anderson,
            //    Visual Interaction with Overhauser Curves and Surfaces,
            //    SIGGRAPH 77,
            //    in Proceedings of the 4th Annual Conference on Computer Graphics
            //    and Interactive Techniques,
            //    ASME, July 1977, pages 132-137.
            //
            //    Edwin Catmull, Raphael Rom,
            //    A Class of Local Interpolating Splines,
            //    in Computer Aided Geometric Design,
            //    edited by Robert Barnhill, Richard Reisenfeld,
            //    Academic Press, 1974,
            //    ISBN: 0120790505.
            //
            //    David Rogers, Alan Adams,
            //    Mathematical Elements of Computer Graphics,
            //    Second Edition,
            //    McGraw Hill, 1989,
            //    ISBN: 0070535299.
            //
            //  Parameters:
            //
            //    Input, int NDIM, the dimension of a single data point.
            //    NDIM must be at least 1.
            //
            //    Input, int NDATA, the number of data points.
            //    NDATA must be at least 3.
            //
            //    Input, double TDATA[NDATA], the abscissas of the data points.  The
            //    values in TDATA must be in strictly ascending order.
            //
            //    Input, double YDATA[NDIM*NDATA], the data points corresponding to
            //    the abscissas.
            //
            //    Input, double TVAL, the abscissa value at which the spline
            //    is to be evaluated.  Normally, TDATA[0] <= TVAL <= T[NDATA-1], and 
            //    the data will be interpolated.  For TVAL outside this range, 
            //    extrapolation will be used.
            //
            //    Output, double YVAL[NDIM], the value of the spline at TVAL.
            //
        {
            int i;
            int left = 0;
            int order;
            int right = 0;
            double[] yl;
            double[] yr;
            //
            //  Check.
            //
            order = typeMethods.r8vec_order_type(ndata, tdata);

            if (order != 2)
            {
                Console.WriteLine("");
                Console.WriteLine("SPLINE_OVERHAUSER_VAL - Fatal error!");
                Console.WriteLine("  The data abscissas are not strictly ascending.");
                return;
            }

            if (ndata < 3)
            {
                Console.WriteLine("");
                Console.WriteLine("SPLINE_OVERHAUSER_VAL - Fatal error!");
                Console.WriteLine("  NDATA < 3.");
                return;
            }

            // 
            //  Locate the abscissa interval T[LEFT], T[LEFT+1] nearest to or 
            //  containing TVAL. 
            //
            typeMethods.r8vec_bracket(ndata, tdata, tval, ref left, ref right);
            // 
            //  Evaluate the "left hand" quadratic defined at 
            //  T[LEFT-1], T[LEFT], T[RIGHT]. 
            //
            yl = new double[ndim];
            yr = new double[ndim];

            if (0 < left - 1)
            {
                Parabola.parabola_val2(ndim, ndata, tdata, ydata, left - 1, tval, ref yl);
            }

            // 
            //  Evaluate the "right hand" quadratic defined at 
            //  T[LEFT], T[RIGHT], T[RIGHT+1]. 
            //
            if (right + 1 <= ndata)
            {
                Parabola.parabola_val2(ndim, ndata, tdata, ydata, left, tval, ref yr);
            }

            //
            //  Blend the quadratics. 
            //
            if (left == 1)
            {
                for (i = 0; i < ndim; i++)
                {
                    yval[i] = yr[i];
                }
            }
            else if (right < ndata)
            {
                for (i = 0; i < ndim; i++)
                {
                    yval[i] = (
                                  (tdata[right - 1] - tval) * yl[i]
                                  + (tval - tdata[left - 1]) * yr[i])
                              / (tdata[right - 1] - tdata[left - 1]);
                }
            }
            else
            {
                for (i = 0; i < ndim; i++)
                {
                    yval[i] = yl[i];
                }
            }
        }
    }
}