using System;

namespace Burkardt.PolynomialNS;

public static class Cubic
{
    public static int chfev(double x1, double x2, double f1, double f2, double d1, double d2,
            int ne, double[] xe, ref double[] fe, ref int[] next, int xeIndex = 0, int feIndex = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHFEV evaluates a cubic polynomial given in Hermite form.
        //
        //  Discussion:
        //
        //    This routine evaluates a cubic polynomial given in Hermite form at an
        //    array of points.  While designed for use by SPLINE_PCHIP_VAL, it may
        //    be useful directly as an evaluator for a piecewise cubic
        //    Hermite function in applications, such as graphing, where
        //    the interval is known in advance.
        //
        //    The cubic polynomial is determined by function values
        //    F1, F2 and derivatives D1, D2 on the interval [X1,X2].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 August 2005
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Fred Fritsch, Lawrence Livermore National Laboratory.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Fred Fritsch, Ralph Carlson, 
        //    Monotone Piecewise Cubic Interpolation,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 17, Number 2, April 1980, pages 238-246.
        //
        //    David Kahaner, Cleve Moler, Steven Nash,
        //    Numerical Methods and Software,
        //    Prentice Hall, 1989,
        //    ISBN: 0-13-627258-4,
        //    LC: TA345.K34.
        //
        //  Parameters:
        //
        //    Input, double X1, X2, the endpoints of the interval of
        //    definition of the cubic.  X1 and X2 must be distinct.
        //
        //    Input, double F1, F2, the values of the function at X1 and
        //    X2, respectively.
        //
        //    Input, double D1, D2, the derivative values at X1 and
        //    X2, respectively.
        //
        //    Input, int NE, the number of evaluation points.
        //
        //    Input, double XE[NE], the points at which the function is to
        //    be evaluated.  If any of the XE are outside the interval
        //    [X1,X2], a warning error is returned in NEXT.
        //
        //    Output, double FE[NE], the value of the cubic function
        //    at the points XE.
        //
        //    Output, int NEXT[2], indicates the number of extrapolation points:
        //    NEXT[0] = number of evaluation points to the left of interval.
        //    NEXT[1] = number of evaluation points to the right of interval.
        //
        //    Output, int CHFEV, error flag.
        //    0, no errors.
        //    -1, NE < 1.
        //    -2, X1 == X2.
        //
    {
        int i;
        int ierr;

        switch (ne)
        {
            case < 1:
                ierr = -1;
                Console.WriteLine("");
                Console.WriteLine("CHFEV - Fatal error!");
                Console.WriteLine("  Number of evaluation points is less than 1.");
                Console.WriteLine("  NE = " + ne + "");
                return ierr;
        }

        double h = x2 - x1;

        switch (h)
        {
            case 0.0:
                ierr = -2;
                Console.WriteLine("");
                Console.WriteLine("CHFEV - Fatal error!");
                Console.WriteLine("  The interval [X1,X2] is of zero length.");
                return ierr;
        }

        //
        //  Initialize.
        //
        ierr = 0;
        next[0] = 0;
        next[1] = 0;
        double xmi = Math.Min(0.0, h);
        double xma = Math.Max(0.0, h);
        //
        //  Compute cubic coefficients expanded about X1.
        //
        double delta = (f2 - f1) / h;
        double del1 = (d1 - delta) / h;
        double del2 = (d2 - delta) / h;
        double c2 = -(del1 + del1 + del2);
        double c3 = (del1 + del2) / h;
        //
        //  Evaluation loop.
        //
        for (i = 0; i < ne; i++)
        {
            double x = xe[xeIndex + i] - x1;
            fe[feIndex + i] = f1 + x * (d1 + x * (c2 + x * c3));
            //
            //  Count the extrapolation points.
            //
            if (x < xmi)
            {
                next[0] += 1;
            }

            if (xma < x)
            {
                next[1] += 1;
            }

        }

        return 0;
    }
}