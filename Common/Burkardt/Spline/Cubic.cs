using System;
using Burkardt.SolveNS;
using Burkardt.Types;

namespace Burkardt.Spline
{
    public static class Cubic
    {
        public static double[] spline_cubic_set(int n, double[] t, double[] y, int ibcbeg,
                double ybcbeg, int ibcend, double ybcend)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPLINE_CUBIC_SET computes the second derivatives of a piecewise cubic spline.
            //
            //  Discussion:
            //
            //    For data interpolation, the user must call SPLINE_SET to determine
            //    the second derivative data, passing in the data to be interpolated,
            //    and the desired boundary conditions.
            //
            //    The data to be interpolated, plus the SPLINE_SET output, defines
            //    the spline.  The user may then call SPLINE_VAL to evaluate the
            //    spline at any point.
            //
            //    The cubic spline is a piecewise cubic polynomial.  The intervals
            //    are determined by the "knots" or abscissas of the data to be
            //    interpolated.  The cubic spline has continous first and second
            //    derivatives over the entire interval of interpolation.
            //
            //    For any point T in the interval T(IVAL), T(IVAL+1), the form of
            //    the spline is
            //
            //      SPL(T) = A(IVAL)
            //             + B(IVAL) * ( T - T(IVAL) )
            //             + C(IVAL) * ( T - T(IVAL) )^2
            //             + D(IVAL) * ( T - T(IVAL) )^3
            //
            //    If we assume that we know the values Y(*) and YPP(*), which represent
            //    the values and second derivatives of the spline at each knot, then
            //    the coefficients can be computed as:
            //
            //      A(IVAL) = Y(IVAL)
            //      B(IVAL) = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
            //        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
            //      C(IVAL) = YPP(IVAL) / 2
            //      D(IVAL) = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
            //
            //    Since the first derivative of the spline is
            //
            //      SPL'(T) =     B(IVAL)
            //              + 2 * C(IVAL) * ( T - T(IVAL) )
            //              + 3 * D(IVAL) * ( T - T(IVAL) )^2,
            //
            //    the requirement that the first derivative be continuous at interior
            //    knot I results in a total of N-2 equations, of the form:
            //
            //      B(IVAL-1) + 2 C(IVAL-1) * (T(IVAL)-T(IVAL-1))
            //      + 3 * D(IVAL-1) * (T(IVAL) - T(IVAL-1))^2 = B(IVAL)
            //
            //    or, setting H(IVAL) = T(IVAL+1) - T(IVAL)
            //
            //      ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
            //      - ( YPP(IVAL) + 2 * YPP(IVAL-1) ) * H(IVAL-1) / 6
            //      + YPP(IVAL-1) * H(IVAL-1)
            //      + ( YPP(IVAL) - YPP(IVAL-1) ) * H(IVAL-1) / 2
            //      =
            //      ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
            //      - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * H(IVAL) / 6
            //
            //    or
            //
            //      YPP(IVAL-1) * H(IVAL-1) + 2 * YPP(IVAL) * ( H(IVAL-1) + H(IVAL) )
            //      + YPP(IVAL) * H(IVAL)
            //      =
            //      6 * ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
            //      - 6 * ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
            //
            //    Boundary conditions must be applied at the first and last knots.  
            //    The resulting tridiagonal system can be solved for the YPP values.
            //
            //    Anton Reinhard corrected the assignments:
            //      a2[i] = ( t[i+1] - t[i]   ) / 6.0;
            //      a3[i] = ( t[i+1] - t[i-1] ) / 3.0;
            //      a4[i] = ( t[i]   - t[i-1] ) / 6.0;
            //    to
            //      a2[i] = ( t[i]   - t[i-1] ) / 6.0;
            //      a3[i] = ( t[i+1] - t[i-1] ) / 3.0;
            //      a4[i] = ( t[i+1] - t[i]   ) / 6.0;
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 August 2020
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
            //  Input:
            //
            //    int N, the number of data points.  N must be at least 2.
            //    In the special case where N = 2 and IBCBEG = IBCEND = 0, the
            //    spline will actually be linear.
            //
            //    double T[N], the knot values, that is, the points were data is
            //    specified.  The knot values should be distinct, and increasing.
            //
            //    double Y[N], the data values to be interpolated.
            //
            //    int IBCBEG, left boundary condition flag:
            //    0: the cubic spline should be a quadratic over the first interval;
            //    1: the first derivative at the left endpoint should be YBCBEG;
            //    2: the second derivative at the left endpoint should be YBCBEG;
            //    3: Not-a-knot: the third derivative is continuous at T(2).
            //
            //    double YBCBEG, the values to be used in the boundary
            //    conditions if IBCBEG is equal to 1 or 2.
            //
            //    int IBCEND, right boundary condition flag:
            //    0: the cubic spline should be a quadratic over the last interval;
            //    1: the first derivative at the right endpoint should be YBCEND;
            //    2: the second derivative at the right endpoint should be YBCEND;
            //    3: Not-a-knot: the third derivative is continuous at T(N-1).
            //
            //    double YBCEND, the values to be used in the boundary
            //    conditions if IBCEND is equal to 1 or 2.
            //
            //  Output:
            //
            //    double SPLINE_CUBIC_SET[N], the second derivatives 
            //    of the cubic spline.
            //
        {
            double[] a1;
            double[] a2;
            double[] a3;
            double[] a4;
            double[] a5;
            double[] b;
            int i;
            double[] ypp;
            //
            //  Check.
            //
            if (n <= 1)
            {
                Console.WriteLine("");
                Console.WriteLine("SPLINE_CUBIC_SET - Fatal error!");
                Console.WriteLine("  The number of data points N must be at least 2.");
                Console.WriteLine("  The input value is " + n + ".");
                return null;
            }

            for (i = 0; i < n - 1; i++)
            {
                if (t[i + 1] <= t[i])
                {
                    Console.WriteLine("");
                    Console.WriteLine("SPLINE_CUBIC_SET - Fatal error!");
                    Console.WriteLine("  The knots must be strictly increasing, but");
                    Console.WriteLine("  T(" + i + ") = " + t[i] + "");
                    Console.WriteLine("  T(" + i + 1 + ") = " + t[i + 1] + "");
                    return null;
                }
            }

            a1 = new double[n];
            a2 = new double[n];
            a3 = new double[n];
            a4 = new double[n];
            a5 = new double[n];
            b = new double[n];

            for (i = 0; i < n; i++)
            {
                a1[i] = 0.0;
                a2[i] = 0.0;
                a3[i] = 0.0;
                a4[i] = 0.0;
                a5[i] = 0.0;
            }

            //
            //  Set up the first equation.
            //
            if (ibcbeg == 0)
            {
                b[0] = 0.0;
                a3[0] = 1.0;
                a4[0] = -1.0;
            }
            else if (ibcbeg == 1)
            {
                b[0] = (y[1] - y[0]) / (t[1] - t[0]) - ybcbeg;
                a3[0] = (t[1] - t[0]) / 3.0;
                a4[0] = (t[1] - t[0]) / 6.0;
            }
            else if (ibcbeg == 2)
            {
                b[0] = ybcbeg;
                a3[0] = 1.0;
                a4[0] = 0.0;
            }
            else if (ibcbeg == 3)
            {
                b[0] = 0.0;
                a3[0] = -(t[2] - t[1]);
                a4[0] = (t[2] - t[0]);
                a5[0] = -(t[1] - t[0]);
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("SPLINE_CUBIC_SET - Fatal error!");
                Console.WriteLine("  IBCBEG must be 0, 1, 2, or 3.");
                Console.WriteLine("  The input value is " + ibcbeg + ".");
                return null;
            }

            //
            //  Set up the intermediate equations.
            //  Note that these lines have been corrected as suggested
            //  by Anton Reinhard, 17 August 2020.
            //
            for (i = 1; i < n - 1; i++)
            {
                b[i] = (y[i + 1] - y[i]) / (t[i + 1] - t[i])
                       - (y[i] - y[i - 1]) / (t[i] - t[i - 1]);
                a2[i] = (t[i] - t[i - 1]) / 6.0;
                a3[i] = (t[i + 1] - t[i - 1]) / 3.0;
                a4[i] = (t[i + 1] - t[i]) / 6.0;
            }

            //
            //  Set up the last equation.
            //
            if (ibcend == 0)
            {
                b[n - 1] = 0.0;
                a2[n - 1] = -1.0;
                a3[n - 1] = 1.0;
            }
            else if (ibcend == 1)
            {
                b[n - 1] = ybcend - (y[n - 1] - y[n - 2]) / (t[n - 1] - t[n - 2]);
                a2[n - 1] = (t[n - 1] - t[n - 2]) / 6.0;
                a3[n - 1] = (t[n - 1] - t[n - 2]) / 3.0;
            }
            else if (ibcend == 2)
            {
                b[n - 1] = ybcend;
                a2[n - 1] = 0.0;
                a3[n - 1] = 1.0;
            }
            else if (ibcbeg == 3)
            {
                b[n - 1] = 0.0;
                a1[n - 1] = -(t[n - 1] - t[n - 2]);
                a2[n - 1] = (t[n - 1] - t[n - 3]);
                a3[n - 1] = -(t[n - 2] - t[n - 3]);
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("SPLINE_CUBIC_SET - Fatal error!");
                Console.WriteLine("  IBCEND must be 0, 1, 2, or 3.");
                Console.WriteLine("  The input value is " + ibcend + ".");
                return null;
            }

            //
            //  Solve the linear system.
            //
            if (n == 2 && ibcbeg == 0 && ibcend == 0)
            {
                ypp = new double[2];

                ypp[0] = 0.0;
                ypp[1] = 0.0;
            }
            else
            {
                ypp = Pentadiagonal.penta(n, a1, a2, a3, a4, a5, b);
            }

            return ypp;
        }

        public static double spline_cubic_val(int n, double[] t, double[] y, double[] ypp,
                double tval, ref double ypval, ref double yppval)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPLINE_CUBIC_VAL evaluates a piecewise cubic spline at a point.
            //
            //  Discussion:
            //
            //    SPLINE_CUBIC_SET must have already been called to define the values of YPP.
            //
            //    For any point T in the interval T(IVAL), T(IVAL+1), the form of
            //    the spline is
            //
            //      SPL(T) = A
            //             + B * ( T - T(IVAL) )
            //             + C * ( T - T(IVAL) )^2
            //             + D * ( T - T(IVAL) )^3
            //
            //    Here:
            //      A = Y(IVAL)
            //      B = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
            //        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
            //      C = YPP(IVAL) / 2
            //      D = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    25 August 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Input:
            //
            //    int N, the number of knots.
            //
            //    double T[N], the knot values.
            //
            //    double Y[N], the data values at the knots.
            //
            //    double YPP[N], the second derivatives of the spline at
            //    the knots.
            //
            //    double TVAL, a point, typically between T[0] and T[N-1], at
            //    which the spline is to be evalulated.  If TVAL lies outside
            //    this range, extrapolation is used.
            //
            //  Output:
            //
            //    double[] YPVAL, the derivative of the spline at TVAL.
            //
            //    double[] YPPVAL, the second derivative of the spline at TVAL.
            //
            //    double SPLINE_VAL, the value of the spline at TVAL.
            //
        {
            double dt;
            double h;
            int i;
            int ival;
            double yval;
            //
            //  Determine the interval [ T(I), T(I+1) ] that contains TVAL.
            //  Values below T[0] or above T[N-1] use extrapolation.
            //
            ival = n - 2;

            for (i = 0; i < n - 1; i++)
            {
                if (tval < t[i + 1])
                {
                    ival = i;
                    break;
                }
            }

            //
            //  In the interval I, the polynomial is in terms of a normalized
            //  coordinate between 0 and 1.
            //
            dt = tval - t[ival];
            h = t[ival + 1] - t[ival];

            yval = y[ival]
                   + dt * ((y[ival + 1] - y[ival]) / h
                           - (ypp[ival + 1] / 6.0 + ypp[ival] / 3.0) * h
                           + dt * (0.5 * ypp[ival]
                                   + dt * ((ypp[ival + 1] - ypp[ival]) / (6.0 * h))));

            ypval = (y[ival + 1] - y[ival]) / h
                    - (ypp[ival + 1] / 6.0 + ypp[ival] / 3.0) * h
                    + dt * (ypp[ival]
                            + dt * (0.5 * (ypp[ival + 1] - ypp[ival]) / h));

            yppval = ypp[ival] + dt * (ypp[ival + 1] - ypp[ival]) / h;

            return yval;
        }

        public static void spline_cubic_val2(int n, double[] t, double tval, ref int left, double[] y,
                double[] ypp, ref double yval, ref double ypval, ref double yppval)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPLINE_CUBIC_VAL2 evaluates a piecewise cubic spline at a point.
            //
            //  Discussion:
            //
            //    This routine is a modification of SPLINE_CUBIC_VAL; it allows the
            //    user to speed up the code by suggesting the appropriate T interval
            //    to search first.
            //
            //    SPLINE_CUBIC_SET must have already been called to define the
            //    values of YPP.
            //
            //    In the LEFT interval, let RIGHT = LEFT+1.  The form of the spline is
            //
            //    SPL(T) =
            //      A
            //    + B * ( T - T[LEFT] )
            //    + C * ( T - T[LEFT] )^2
            //    + D * ( T - T[LEFT] )^3
            //
            //    Here:
            //      A = Y[LEFT]
            //      B = ( Y[RIGHT] - Y[LEFT] ) / ( T[RIGHT] - T[LEFT] )
            //        - ( YPP[RIGHT] + 2 * YPP[LEFT] ) * ( T[RIGHT] - T[LEFT] ) / 6
            //      C = YPP[LEFT] / 2
            //      D = ( YPP[RIGHT] - YPP[LEFT] ) / ( 6 * ( T[RIGHT] - T[LEFT] ) )
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
            //    Input, int N, the number of knots.
            //
            //    Input, double T[N], the knot values.
            //
            //    Input, double TVAL, a point, typically between T[0] and T[N-1], at
            //    which the spline is to be evalulated.  If TVAL lies outside
            //    this range, extrapolation is used.
            //
            //    Input/output, int *LEFT, the suggested T interval to search.
            //    LEFT should be between 1 and N-1.  If LEFT is not in this range,
            //    then its value will be ignored.  On output, LEFT is set to the
            //    actual interval in which TVAL lies.
            //
            //    Input, double Y[N], the data values at the knots.
            //
            //    Input, double YPP[N], the second derivatives of the spline at
            //    the knots.
            //
            //    Output, double[] YVAL, *YPVAL, *YPPVAL, the value of the spline, and
            //    its first two derivatives at TVAL.
            //
        {
            double dt;
            double h;
            int right;
            //
            //  Determine the interval [T[LEFT], T[RIGHT]] that contains TVAL.  
            //  
            //  What you want from R8VEC_BRACKET3 is that TVAL is to be computed
            //  by the data in interval [T[LEFT-1], T[RIGHT-1]].  
            //
            typeMethods.r8vec_bracket3(n, t, tval, ref left);
            //
            // In the interval LEFT, the polynomial is in terms of a normalized
            // coordinate  ( DT / H ) between 0 and 1.
            //
            right = left + 1;

            int tlIndex = left % t.Length;
            if (tlIndex < 0)
            {
                tlIndex += t.Length;
            }
            
            int trIndex = right % t.Length;
            if (trIndex < 0)
            {
                trIndex += t.Length;
            }

            int ylIndex = left % y.Length;
            if (ylIndex < 0)
            {
                ylIndex += y.Length;
            }
            
            int yrIndex = right % y.Length;
            if (yrIndex < 0)
            {
                yrIndex += y.Length;
            }

            int ypplIndex = left % ypp.Length;
            if (ypplIndex < 0)
            {
                ypplIndex += ypp.Length;
            }
            
            int ypprIndex = right % ypp.Length;
            if (ypprIndex < 0)
            {
                ypprIndex += ypp.Length;
            }

            dt = tval - t[tlIndex];
            h = t[trIndex] - t[tlIndex];

            yval = y[ylIndex]
                   + dt * ((y[yrIndex] - y[ylIndex]) / h
                           - (ypp[ypprIndex] / 6.0 + ypp[ypplIndex] / 3.0) * h
                           + dt * (0.5 * ypp[ypplIndex]
                                   + dt * ((ypp[ypprIndex] - ypp[ypplIndex]) / (6.0 * h))));

            ypval = (y[yrIndex] - y[ylIndex]) / h
                    - (ypp[ypprIndex] / 6.0 + ypp[ypplIndex] / 3.0) * h
                    + dt * (ypp[ypplIndex]
                            + dt * (0.5 * (ypp[ypprIndex] - ypp[ypplIndex]) / h));

            yppval = ypp[ypplIndex] + dt * (ypp[ypprIndex] - ypp[ypplIndex]) / h;

        }

    }
}