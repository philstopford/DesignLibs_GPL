using Burkardt.Types;

namespace Burkardt.Spline;

public static class Hermite
{
    public static double[] spline_hermite_set(int ndata, double[] tdata, double[] ydata,
            double[] ypdata)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPLINE_HERMITE_SET sets up a piecewise cubic Hermite interpolant.
        //
        //  Discussion:
        //
        //    Once the array C is computed, then in the interval
        //    (TDATA(I), TDATA(I+1)), the interpolating Hermite polynomial
        //    is given by
        //
        //      SVAL(TVAL) =                 C(1,I)
        //         + ( TVAL - TDATA(I) ) * ( C(2,I)
        //         + ( TVAL - TDATA(I) ) * ( C(3,I)
        //         + ( TVAL - TDATA(I) ) *   C(4,I) ) )
        //
        //    This is algorithm CALCCF from Conte and deBoor.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 February 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Samuel Conte, Carl deBoor,
        //    Elementary Numerical Analysis,
        //    Second Edition,
        //    McGraw Hill, 1972,
        //    ISBN: 07-012446-4.
        //
        //  Parameters:
        //
        //    Input, int NDATA, the number of data points.
        //    NDATA must be at least 2.
        //
        //    Input, double TDATA[NDATA], the abscissas of the data points.
        //    The entries of TDATA are assumed to be strictly increasing.
        //
        //    Input, double Y[NDATA], YP[NDATA], the value of the
        //    function and its derivative at TDATA(1:NDATA).
        //
        //    Output, double SPLINE_HERMITE_SET[4*NDATA], the coefficients of 
        //    the Hermite polynomial.  We will refer to this array as "C".
        //    C(1,1:NDATA) = Y(1:NDATA) and C(2,1:NDATA) = YP(1:NDATA).
        //    C(3,1:NDATA-1) and C(4,1:NDATA-1) are the quadratic and cubic
        //    coefficients.
        //
    {
        int i;
        int j;

        double[] c = new double[4 * ndata];

        for (j = 0; j < ndata; j++)
        {
            c[0 + j * 4] = ydata[j];
        }

        for (j = 0; j < ndata; j++)
        {
            c[1 + j * 4] = ypdata[j];
        }

        for (i = 1; i <= ndata - 1; i++)
        {
            double dt = tdata[i] - tdata[i - 1];
            double divdif1 = (c[0 + i * 4] - c[0 + (i - 1) * 4]) / dt;
            double divdif3 = c[1 + (i - 1) * 4] + c[1 + i * 4] - 2.0 * divdif1;
            c[2 + (i - 1) * 4] = (divdif1 - c[1 + (i - 1) * 4] - divdif3) / dt;
            c[3 + (i - 1) * 4] = divdif3 / (dt * dt);
        }

        c[2 + (ndata - 1) * 4] = 0.0;
        c[3 + (ndata - 1) * 4] = 0.0;

        return c;
    }

    public static void spline_hermite_val(int ndata, double[] tdata, double[] c, double tval,
            ref double sval, ref double spval)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPLINE_HERMITE_VAL evaluates a piecewise cubic Hermite interpolant.
        //
        //  Discussion:
        //
        //    SPLINE_HERMITE_SET must be called first, to set up the
        //    spline data from the raw function and derivative data.
        //
        //    In the interval (TDATA(I), TDATA(I+1)), the interpolating
        //    Hermite polynomial is given by
        //
        //      SVAL(TVAL) =                 C(1,I)
        //         + ( TVAL - TDATA(I) ) * ( C(2,I)
        //         + ( TVAL - TDATA(I) ) * ( C(3,I)
        //         + ( TVAL - TDATA(I) ) *   C(4,I) ) )
        //
        //    and
        //
        //      SVAL'(TVAL) =                    C(2,I)
        //         + ( TVAL - TDATA(I) ) * ( 2 * C(3,I)
        //         + ( TVAL - TDATA(I) ) *   3 * C(4,I) )
        //
        //    This is algorithm PCUBIC from Conte and deBoor.
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
        //    Samuel Conte, Carl deBoor,
        //    Elementary Numerical Analysis,
        //    Second Edition,
        //    McGraw Hill, 1972,
        //    ISBN: 07-012446-4.
        //
        //  Parameters:
        //
        //    Input, int NDATA, the number of data points.
        //    NDATA is assumed to be at least 2.
        //
        //    Input, double TDATA[NDATA], the abscissas of the data points.
        //    The entries of TDATA are assumed to be strictly increasing.
        //
        //    Input, double C[4*NDATA], the coefficient data computed by
        //    SPLINE_HERMITE_SET.
        //
        //    Input, double TVAL, the point where the interpolant is to
        //    be evaluated.
        //
        //    Output, double *SVAL, *SPVAL, the value of the interpolant
        //    and its derivative at TVAL.
        //
    {
        int left = 0;
        int right = 0;
        //
        //  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains
        //  or is nearest to TVAL.
        //
        typeMethods.r8vec_bracket(ndata, tdata, tval, ref left, ref right);
        //
        //  Evaluate the cubic polynomial.
        //
        double dt = tval - tdata[left - 1];

        sval = c[0 + (left - 1) * 4]
               + dt * (c[1 + (left - 1) * 4]
                       + dt * (c[2 + (left - 1) * 4]
                               + dt * c[3 + (left - 1) * 4]));

        spval = c[1 + (left - 1) * 4]
                + dt * (2.0 * c[2 + (left - 1) * 4]
                        + dt * 3.0 * c[3 + (left - 1) * 4]);

    }
}