using System;
using Burkardt.PolynomialNS;

namespace Burkardt.Spline;

public static class Bezier
{
    public static void bc_val(int n, double t, double[] xcon, double[] ycon, ref double xval,
            ref double yval)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BC_VAL evaluates a parameterized Bezier curve.
        //
        //  Discussion:
        //
        //    BC_VAL(T) is the value of a vector function of the form
        //
        //      BC_VAL(T) = ( X(T), Y(T) )
        //
        //    where
        //
        //      X(T) = Sum ( 0 <= I <= N ) XCON(I) * BERN(I,N)(T)
        //      Y(T) = Sum ( 0 <= I <= N ) YCON(I) * BERN(I,N)(T)
        //
        //    BERN(I,N)(T) is the I-th Bernstein polynomial of order N
        //    defined on the interval [0,1],
        //
        //    XCON(0:N) and YCON(0:N) are the coordinates of N+1 "control points".
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 February 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    David Kahaner, Cleve Moler, Steven Nash,
        //    Numerical Methods and Software,
        //    Prentice Hall, 1989,
        //    ISBN: 0-13-627258-4,
        //    LC: TA345.K34.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the Bezier curve, which
        //    must be at least 0.
        //
        //    Input, double T, the point at which the Bezier curve should
        //    be evaluated.  The best results are obtained within the interval
        //    [0,1] but T may be anywhere.
        //
        //    Input, double XCON[0:N], YCON[0:N], the X and Y coordinates
        //    of the control points.  The Bezier curve will pass through
        //    the points ( XCON(0), YCON(0) ) and ( XCON(N), YCON(N) ), but
        //    generally NOT through the other control points.
        //
        //    Output, double *XVAL, *YVAL, the X and Y coordinates of the point
        //    on the Bezier curve corresponding to the given T value.
        //
    {
        int i;

        double[] bval = BernsteinPolynomial.bernstein_poly_01(n, t);

        xval = 0.0;
        for (i = 0; i <= n; i++)
        {
            xval += xcon[i] * bval[i];
        }

        yval = 0.0;
        for (i = 0; i <= n; i++)
        {
            yval += ycon[i] * bval[i];
        }
    }

    public static double bez_val(int n, double x, double a, double b, double[] y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BEZ_VAL evaluates a Bezier function at a point.
        //
        //  Discussion:
        //
        //    The Bezier function has the form:
        //
        //      BEZ(X) = Sum ( 0 <= I <= N ) Y(I) * BERN(N,I)( (X-A)/(B-A) )
        //
        //    BERN(N,I)(X) is the I-th Bernstein polynomial of order N
        //    defined on the interval [0,1],
        //
        //    Y(0:N) is a set of coefficients,
        //
        //    and if, for I = 0 to N, we define the N+1 points
        //
        //      X(I) = ( (N-I)*A + I*B) / N,
        //
        //    equally spaced in [A,B], the pairs ( X(I), Y(I) ) can be regarded as
        //    "control points".
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 February 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    David Kahaner, Cleve Moler, Steven Nash,
        //    Numerical Methods and Software,
        //    Prentice Hall, 1989,
        //    ISBN: 0-13-627258-4,
        //    LC: TA345.K34.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the Bezier function, which
        //    must be at least 0.
        //
        //    Input, double X, the point at which the Bezier function should
        //    be evaluated.  The best results are obtained within the interval
        //    [A,B] but X may be anywhere.
        //
        //    Input, double A, B, the interval over which the Bezier function
        //    has been defined.  This is the interval in which the control
        //    points have been set up.  Note BEZ(A) = Y(0) and BEZ(B) = Y(N),
        //    although BEZ will not, in general pass through the other
        //    control points.  A and B must not be equal.
        //
        //    Input, double Y[0:N], a set of data defining the Y coordinates
        //    of the control points.
        //
        //    Output, double BEZ_VAL, the value of the Bezier function at X.
        //
    {
        int i;

        switch (b - a)
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("BEZ_VAL - Fatal error!");
                Console.WriteLine("  Null interval, A = B = " + a + "");
                return 1;
        }

        //
        //  X01 lies in [0,1], in the same relative position as X in [A,B].
        //
        double x01 = (x - a) / (b - a);

        double[] bval = BernsteinPolynomial.bernstein_poly_01(n, x01);

        double value = 0.0;
        for (i = 0; i <= n; i++)
        {
            value += y[i] * bval[i];
        }

        return value;
    }
}