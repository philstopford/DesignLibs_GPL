using System;
using Burkardt.MatrixNS;
using Burkardt.Types;

namespace Burkardt.Spline;

public static class Linear
{
    public static double spline_linear_int(int ndata, double[] tdata, double[] ydata,
            double a, double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPLINE_LINEAR_INT evaluates the integral of a piecewise linear spline.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 January 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NDATA, the number of data points defining the spline.
        //
        //    Input, double TDATA[NDATA], YDATA[NDATA], the values of the independent
        //    and dependent variables at the data points.  The values of TDATA should
        //    be distinct and increasing.
        //
        //    Input, double A, B, the interval over which the integral is desired.
        //
        //    Output, double SPLINE_LINEAR_INT, the value of the integral.
        //
    {
        int a_left = 0;
        int a_right = 0;
        int b_left = 0;
        int b_right = 0;
        int i_left;
        double tval;
        double yp;
        double yval;

        double int_val = 0.0;

        if (Math.Abs(a - b) <= typeMethods.r8_epsilon())
        {
            return int_val;
        }

        double a_copy = Math.Min(a, b);
        double b_copy = Math.Max(a, b);
        //
        //  Find the interval [ TDATA(A_LEFT), TDATA(A_RIGHT) ] that contains, or is
        //  nearest to, A.
        //
        typeMethods.r8vec_bracket(ndata, tdata, a_copy, ref a_left, ref a_right);
        //
        //  Find the interval [ TDATA(B_LEFT), TDATA(B_RIGHT) ] that contains, or is
        //  nearest to, B.
        //
        typeMethods.r8vec_bracket(ndata, tdata, b_copy, ref b_left, ref b_right);
        //
        //  If A and B are in the same interval...
        //
        if (a_left == b_left)
        {
            tval = (a_copy + b_copy) / 2.0;

            yp = (ydata[a_right - 1] - ydata[a_left - 1]) /
                 (tdata[a_right - 1] - tdata[a_left - 1]);

            yval = ydata[a_left - 1] + (tval - tdata[a_left - 1]) * yp;

            int_val = yval * (b_copy - a_copy);

            return int_val;
        }

        //
        //  Otherwise, integrate from:
        //
        //  A               to TDATA(A_RIGHT),
        //  TDATA(A_RIGHT)  to TDATA(A_RIGHT+1),...
        //  TDATA(B_LEFT-1) to TDATA(B_LEFT),
        //  TDATA(B_LEFT)   to B.
        //
        //  Use the fact that the integral of a linear function is the
        //  value of the function at the midpoint times the width of the interval.
        //
        tval = (a_copy + tdata[a_right - 1]) / 2.0;

        yp = (ydata[a_right - 1] - ydata[a_left - 1]) /
             (tdata[a_right - 1] - tdata[a_left - 1]);

        yval = ydata[a_left - 1] + (tval - tdata[a_left - 1]) * yp;

        int_val += yval * (tdata[a_right - 1] - a_copy);

        for (i_left = a_right; i_left <= b_left - 1; i_left++)
        {
            tval = (tdata[i_left] + tdata[i_left - 1]) / 2.0;

            yp = (ydata[i_left - 1] - ydata[i_left - 2]) /
                 (tdata[i_left - 1] - tdata[i_left - 2]);

            yval = ydata[i_left - 2] + (tval - tdata[i_left - 2]) * yp;

            int_val += yval * (tdata[i_left - 1] - tdata[i_left - 2]);
        }

        tval = (tdata[b_left - 1] + b_copy) / 2.0E+0;

        yp = (ydata[b_right - 1] - ydata[b_left - 1]) /
             (tdata[b_right - 1] - tdata[b_left - 1]);

        yval = ydata[b_left - 1] + (tval - tdata[b_left - 1]) * yp;

        int_val += yval * (b_copy - tdata[b_left - 1]);

        if (b < a)
        {
            int_val = -int_val;
        }

        return int_val;
    }

    public static void spline_linear_intset(int n, double[] int_x, double[] int_v,
            ref double[] data_x, ref double[] data_y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPLINE_LINEAR_INTSET sets a piecewise linear spline with given integral properties.
        //
        //  Discussion:
        //
        //    The user has in mind an interval, divided by N+1 points into
        //    N intervals.  A linear spline is to be constructed,
        //    with breakpoints at the centers of each interval, and extending
        //    continuously to the left of the first and right of the last
        //    breakpoints.  The constraint on the linear spline is that it is
        //    required that it have a given integral value over each interval.
        //
        //    A tridiagonal linear system of equations is solved for the
        //    values of the spline at the breakpoints.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 February 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of intervals.
        //
        //    Input, double INT_X[N+1], the points that define the intervals.
        //    Interval I lies between INT_X(I) and INT_X(I+1).
        //
        //    Input, double INT_V[N], the desired value of the integral of the
        //    linear spline over each interval.
        //
        //    Output, double DATA_X[N], DATA_Y[N], the values of the independent
        //    and dependent variables at the data points.  The values of DATA_X are
        //    the interval midpoints.  The values of DATA_Y are determined in such
        //    a way that the exact integral of the linear spline over interval I
        //    is equal to INT_V(I).
        //
    {
        int i;

        double[] a = new double[3 * n];
        double[] b = new double[n];
        //
        //  Set up the easy stuff.
        //
        for (i = 1; i <= n; i++)
        {
            data_x[i - 1] = 0.5 * (int_x[i - 1] + int_x[i]);
        }

        //
        //  Set up the coefficients of the linear system.
        //
        for (i = 0; i < n - 2; i++)
        {
            a[2 + i * 3] = 1.0 - (0.5 * (data_x[i + 1] + int_x[i + 1])
                                  - data_x[i]) / (data_x[i + 1] - data_x[i]);
        }

        a[2 + (n - 2) * 3] = 0.0;
        a[2 + (n - 1) * 3] = 0.0;

        a[1 + 0 * 3] = int_x[1] - int_x[0];

        for (i = 1; i < n - 1; i++)
        {
            a[1 + i * 3] = 1.0 + (0.5 * (data_x[i] + int_x[i])
                                  - data_x[i - 1]) / (data_x[i] - data_x[i - 1])
                           - (0.5 * (data_x[i] + int_x[i + 1]) - data_x[i])
                           / (data_x[i + 1] - data_x[i]);
        }

        a[1 + (n - 1) * 3] = int_x[n] - int_x[n - 1];

        a[0 + 0 * 3] = 0.0;
        a[0 + 1 * 3] = 0.0;
        for (i = 2; i < n; i++)
        {
            a[0 + i * 3] = (0.5 * (data_x[i - 1] + int_x[i])
                            - data_x[i - 1]) / (data_x[i] - data_x[i - 1]);
        }

        //
        //  Set up DATA_Y, which begins as the right hand side of the linear system.
        //
        b[0] = int_v[0];
        for (i = 2; i <= n - 1; i++)
        {
            b[i - 1] = 2.0 * int_v[i - 1] / (int_x[i] - int_x[i - 1]);
        }

        b[n - 1] = int_v[n - 1];
        //
        //  Solve the linear system.
        //
        double[] c = D3.d3_np_fs(n, a, b);

        for (i = 0; i < n; i++)
        {
            data_y[i] = c[i];
        }
    }

    public static void spline_linear_val(int ndata, double[] tdata, double[] ydata,
            double tval, ref double yval, ref double ypval)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPLINE_LINEAR_VAL evaluates a piecewise linear spline at a point.
        //
        //  Discussion:
        //
        //    Because of the extremely simple form of the linear spline,
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
        //
        //  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
        //  nearest to, TVAL.
        //
        typeMethods.r8vec_bracket(ndata, tdata, tval, ref left, ref right);
        //
        //  Now evaluate the piecewise linear function.
        //
        ypval = (ydata[right - 1] - ydata[left - 1])
                / (tdata[right - 1] - tdata[left - 1]);

        yval = ydata[left - 1] + (tval - tdata[left - 1]) * ypval;

    }

}