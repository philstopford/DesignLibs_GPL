using System;
using Burkardt.Types;

namespace Burkardt.Spline;

public static class CubicB
{
    public static double spline_b_val(int ndata, double[] tdata, double[] ydata, double tval)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPLINE_B_VAL evaluates a cubic B spline approximant.
        //
        //  Discussion:
        //
        //    The cubic B spline will approximate the data, but is not
        //    designed to interpolate it.
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
        //  Reference:
        //
        //    Carl deBoor,
        //    A Practical Guide to Splines,
        //    Springer, 2001,
        //    ISBN: 0387953663.
        //
        //  Parameters:
        //
        //    Input, int NDATA, the number of data values.
        //
        //    Input, double TDATA[NDATA], the abscissas of the data.
        //
        //    Input, double YDATA[NDATA], the data values.
        //
        //    Input, double TVAL, a point at which the spline is to be evaluated.
        //
        //    Output, double SPLINE_B_VAL, the value of the function at TVAL.
        //
    {
        int left = 0;
        int right = 0;
        //
        //  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
        //
        typeMethods.r8vec_bracket(ndata, tdata, tval, ref left, ref right);
        //
        //  Evaluate the 5 nonzero B spline basis functions in the interval,
        //  weighted by their corresponding data values.
        //
        double u = (tval - tdata[left - 1]) / (tdata[right - 1] - tdata[left - 1]);
        double yval = 0.0;
        //
        //  B function associated with node LEFT - 1, (or "phantom node"),
        //  evaluated in its 4th interval.
        //
        double bval = (((-1.0
                    * u + 3.0)
                * u - 3.0)
            * u + 1.0) / 6.0;

        yval += (left - 1) switch
        {
            > 0 => ydata[left - 2] * bval,
            _ => (2.0 * ydata[0] - ydata[1]) * bval
        };

        //
        //  B function associated with node LEFT,
        //  evaluated in its third interval.
        //
        bval = (((3.0
                    * u - 6.0)
                * u + 0.0)
            * u + 4.0) / 6.0;

        yval += ydata[left - 1] * bval;
        //
        //  B function associated with node RIGHT,
        //  evaluated in its second interval.
        //
        bval = (((-3.0
                    * u + 3.0)
                * u + 3.0)
            * u + 1.0) / 6.0;

        yval += ydata[right - 1] * bval;
        //
        //  B function associated with node RIGHT+1, (or "phantom node"),
        //  evaluated in its first interval.
        //
        bval = Math.Pow(u, 3) / 6.0;

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