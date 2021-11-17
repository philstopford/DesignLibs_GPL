using System;

namespace Burkardt.Function;

public static class Parabola
{
    public static void parabola_val2(int ndim, int ndata, double[] tdata, double[] ydata,
            int left, double tval, ref double[] yval)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARABOLA_VAL2 evaluates a parabolic function through 3 points in a table.
        //
        //  Discussion:
        //
        //    This routine is a utility routine used by OVERHAUSER_SPLINE_VAL.
        //    It constructs the parabolic interpolant through the data in
        //    3 consecutive entries of a table and evaluates this interpolant
        //    at a given abscissa value.
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
        //  Parameters:
        //
        //    Input, integer NDIM, the dimension of a single data point.
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
        //    Input, int LEFT, the location of the first of the three
        //    consecutive data points through which the parabolic interpolant
        //    must pass.  0 <= LEFT <= NDATA - 3.
        //
        //    Input, double TVAL, the value of T at which the parabolic interpolant
        //    is to be evaluated.  Normally, TDATA[0] <= TVAL <= T[NDATA-1], and 
        //    the data will be interpolated.  For TVAL outside this range, 
        //    extrapolation will be used.
        //
        //    Output, double YVAL[NDIM], the value of the parabolic interpolant 
        //    at TVAL.
        //
    {
        double dif1;
        double dif2;
        int i;
        double t1;
        double t2;
        double t3;
        double y1;
        double y2;
        double y3;
        switch (left)
        {
            //
            //  Check. 
            //
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("PARABOLA_VAL2 - Fatal error!");
                Console.WriteLine("  LEFT < 0.");
                return;
        }

        if (ndata - 2 < left)
        {
            Console.WriteLine("");
            Console.WriteLine("PARABOLA_VAL2 - Fatal error!");
            Console.WriteLine(" NDATA-2 < LEFT.");
            return;
        }

        switch (ndim)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("PARABOLA_VAL2 - Fatal error!");
                Console.WriteLine(" NDIM < 1.");
                return;
        }

        // 
        //  Copy out the three abscissas. 
        //
        t1 = tdata[left - 1];
        t2 = tdata[left];
        t3 = tdata[left + 1];

        if (t2 <= t1 || t3 <= t2)
        {
            Console.WriteLine("");
            Console.WriteLine("PARABOLA_VAL2 - Fatal error!");
            Console.WriteLine("  T2 <= T1 or T3 <= T2.");
            Console.WriteLine("  T1 = " + t1 + "");
            Console.WriteLine("  T2 = " + t2 + "");
            Console.WriteLine("  T3 = " + t3 + "");
            return;
        }

        // 
        //  Construct and evaluate a parabolic interpolant for the data. 
        //
        for (i = 0; i < ndim; i++)
        {
            y1 = ydata[i + (left - 1) * ndim];
            y2 = ydata[i + left * ndim];
            y3 = ydata[i + (left + 1) * ndim];

            dif1 = (y2 - y1) / (t2 - t1);
            dif2 =
                ((y3 - y1) / (t3 - t1)
                 - (y2 - y1) / (t2 - t1)) / (t3 - t2);

            yval[i] = y1 + (tval - t1) * (dif1 + (tval - t2) * dif2);
        }

    }
}