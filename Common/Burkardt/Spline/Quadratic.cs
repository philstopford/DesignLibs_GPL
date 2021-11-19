﻿using System;
using Burkardt.Types;

namespace Burkardt.Spline;

public static class Quadratic
{
    public static void spline_quadratic_val(int ndata, double[] tdata, double[] ydata,
            double tval, ref double yval, ref double ypval )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPLINE_QUADRATIC_VAL evaluates a piecewise quadratic spline at a point.
        //
        //  Discussion:
        //
        //    Because of the simple form of a piecewise quadratic spline,
        //    the raw data points ( TDATA(I), YDATA(I)) can be used directly to
        //    evaluate the spline at any point.  No processing of the data
        //    is required.
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
        //    Input, int NDATA, the number of data points defining the spline.
        //    NDATA should be odd and at least 3.
        //
        //    Input, double TDATA[NDATA], YDATA[NDATA], the values of the independent
        //    and dependent variables at the data points.  The values of TDATA should
        //    be distinct and increasing.
        //
        //    Input, double TVAL, the point at which the spline is to be evaluated.
        //
        //    Output, double *YVAL, *YPVAL, the value of the spline and its first
        //    derivative dYdT at TVAL.  YPVAL is not reliable if TVAL is exactly
        //    equal to TDATA(I) for some I.
        //
    {
        int left = 0;
        int right = 0;

        switch (ndata)
        {
            case < 3:
                Console.WriteLine("");
                Console.WriteLine("SPLINE_QUADRATIC_VAL - Fatal error!");
                Console.WriteLine("  NDATA < 3.");
                return;
        }

        switch (ndata % 2)
        {
            case 0:
                Console.WriteLine("");
                Console.WriteLine("SPLINE_QUADRATIC_VAL - Fatal error!");
                Console.WriteLine("  NDATA must be odd.");
                return;
        }

        //
        //  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
        //  nearest to, TVAL.
        //
        typeMethods.r8vec_bracket(ndata, tdata, tval, ref left, ref right);
        switch (left % 2)
        {
            //
            //  Force LEFT to be odd.
            //
            case 0:
                left -= 1;
                break;
        }

        //
        //  Copy out the three abscissas.
        //
        double t1 = tdata[left - 1];
        double t2 = tdata[left];
        double t3 = tdata[left + 1];

        if (t2 <= t1 || t3 <= t2)
        {
            Console.WriteLine("");
            Console.WriteLine("SPLINE_QUADRATIC_VAL - Fatal error!");
            Console.WriteLine("  T2 <= T1 or T3 <= T2.");
            return;
        }

        //
        //  Construct and evaluate a parabolic interpolant for the data
        //  in each dimension.
        //
        double y1 = ydata[left - 1];
        double y2 = ydata[left];
        double y3 = ydata[left + 1];

        double dif1 = (y2 - y1) / (t2 - t1);

        double dif2 = ((y3 - y1) / (t3 - t1)
                       - (y2 - y1) / (t2 - t1)) / (t3 - t2);

        yval = y1 + (tval - t1) * (dif1 + (tval - t2) * dif2);
        ypval = dif1 + dif2 * (2.0 * tval - t1 - t2);

    }
}