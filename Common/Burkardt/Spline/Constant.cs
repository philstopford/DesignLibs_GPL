namespace Burkardt.Spline;

public static class Constant
{
    public static double spline_constant_val ( int ndata, double[] tdata, double[] ydata, 
            double tval )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPLINE_CONSTANT_VAL evaluates a piecewise constant spline at a point.
        //
        //  Discussion:
        //
        //    NDATA-1 points TDATA define NDATA intervals, with the first
        //    and last being semi-infinite.
        //
        //    The value of the spline anywhere in interval I is YDATA(I).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 February 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NDATA, the number of data points defining the spline.
        //
        //    Input, double TDATA[NDATA-1], the breakpoints.  The values of TDATA should
        //    be distinct and increasing.
        //
        //    Input, double YDATA[NDATA], the values of the spline in the intervals
        //    defined by the breakpoints.
        //
        //    Input, double TVAL, the point at which the spline is to be evaluated.
        //
        //    Output, double *SPLINE_CONSTANT_VAL, the value of the spline at TVAL.
        //
    {
        int i;

        for ( i = 0; i < ndata-1; i++ )
        {
            if ( tval <= tdata[i] )
            {
                return ydata[i];
            }
        }

        return ydata[ndata-1];
    }
}