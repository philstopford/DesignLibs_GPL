using System;
using Burkardt.Types;

namespace Burkardt.Spline;

public static class Beta
{
    public static double spline_beta_val(double beta1, double beta2, int ndata, double[] tdata,
            double[] ydata, double tval )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPLINE_BETA_VAL evaluates a cubic beta spline approximant.
        //
        //  Discussion:
        //
        //    The cubic beta spline will approximate the data, but is not
        //    designed to interpolate it.
        //
        //    If BETA1 = 1 and BETA2 = 0, the cubic beta spline will be the
        //    same as the cubic B spline approximant.
        //
        //    With BETA1 = 1 and BETA2 large, the beta spline becomes more like
        //    a linear spline.
        //
        //    In effect, two "phantom" data values are appended to the data,
        //    so that the spline will interpolate the first and last data values.
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
        //    Input, double BETA1, the skew or bias parameter.
        //    BETA1 = 1 for no skew or bias.
        //
        //    Input, double BETA2, the tension parameter.
        //    BETA2 = 0 for no tension.
        //
        //    Input, int NDATA, the number of data values.
        //
        //    Input, double TDATA[NDATA], the abscissas of the data.
        //
        //    Input, double YDATA[NDATA], the data values.
        //
        //    Input, double TVAL, a point at which the spline is to be evaluated.
        //
        //    Output, double SPLINE_BETA_VAL, the value of the function at TVAL.
        //
    {
        double a;
        double b;
        double bval;
        double c;
        double d;
        double delta;
        int left = 0;
        int right = 0;
        double u;
        double yval;
        //
        //  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
        //
        typeMethods.r8vec_bracket(ndata, tdata, tval, ref left, ref right);
        //
        //  Evaluate the 5 nonzero beta spline basis functions in the interval,
        //  weighted by their corresponding data values.
        //
        u = (tval - tdata[left - 1]) / (tdata[right - 1] - tdata[left - 1]);

        delta = ((2.0
                    * beta1 + 4.0)
                * beta1 + 4.0)
            * beta1 + 2.0 + beta2;

        yval = 0.0;
        //
        //  Beta function associated with node LEFT - 1, (or "phantom node"),
        //  evaluated in its 4th interval.
        //
        bval = 2.0 * Math.Pow(beta1 * (1.0 - u), 3) / delta;

        yval += (left - 1) switch
        {
            > 0 => ydata[left - 2] * bval,
            _ => (2.0 * ydata[0] - ydata[1]) * bval
        };

        //
        //  Beta function associated with node LEFT,
        //  evaluated in its third interval.
        //
        a = beta2 + (4.0 + 4.0 * beta1) * beta1;

        b = -6.0 * beta1 * (1.0 - beta1) * (1.0 + beta1);

        c = ((-6.0
                    * beta1 - 6.0)
                * beta1 + 0.0)
            * beta1 - 3.0 * beta2;

        d = ((+2.0
                    * beta1 + 2.0)
                * beta1 + 2.0)
            * beta1 + 2.0 * beta2;

        bval = (a + u * (b + u * (c + u * d))) / delta;

        yval += ydata[left - 1] * bval;
        //
        //  Beta function associated with node RIGHT,
        //  evaluated in its second interval.
        //
        a = 2.0;

        b = 6.0 * beta1;

        c = 3.0 * beta2 + 6.0 * beta1 * beta1;

        d = -2.0 * (1.0 + beta2 + beta1 + beta1 * beta1);

        bval = (a + u * (b + u * (c + u * d))) / delta;

        yval += ydata[right - 1] * bval;
        //
        //  Beta function associated with node RIGHT+1, (or "phantom node"),
        //  evaluated in its first interval.
        //
        bval = 2.0 * Math.Pow(u, 3) / delta;

        if (right + 1 <= ndata)
        {
            yval += ydata[right] * bval;
        }
        else
        {
            yval += (2.0 * ydata[ndata - 1] - ydata[ndata - 2]) * bval;
        }

        return yval;
    }
}