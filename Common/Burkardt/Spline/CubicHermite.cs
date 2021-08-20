using System;

namespace Burkardt.Spline
{
    public static class CubicHermite
    {
        public static double pchst(double arg1, double arg2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PCHST: PCHIP sign-testing routine.
            //
            //  Discussion:
            //
            //    This routine essentially computes the sign of ARG1 * ARG2.
            //
            //    The object is to do this without multiplying ARG1 * ARG2, to avoid
            //    possible over/underflow problems.
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
            //  Parameters:
            //
            //    Input, double ARG1, ARG2, two values to check.
            //
            //    Output, double PCHST,
            //    -1.0, if ARG1 and ARG2 are of opposite sign.
            //     0.0, if either argument is zero.
            //    +1.0, if ARG1 and ARG2 are of the same sign.
            //
        {
            double value = 0;

            if (arg1 == 0.0)
            {
                value = 0.0;
            }
            else if (arg1 < 0.0)
            {
                if (arg2 < 0.0)
                {
                    value = 1.0;
                }
                else if (arg2 == 0.0)
                {
                    value = 0.0;
                }
                else if (0.0 < arg2)
                {
                    value = -1.0;
                }
            }
            else if (0.0 < arg1)
            {
                if (arg2 < 0.0)
                {
                    value = -1.0;
                }
                else if (arg2 == 0.0)
                {
                    value = 0.0;
                }
                else if (0.0 < arg2)
                {
                    value = 1.0;
                }
            }

            return value;
        }

        public static void spline_pchip_set(int n, double[] x, double[] f, ref double[] d)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPLINE_PCHIP_SET sets derivatives for a piecewise cubic Hermite interpolant.
            //
            //  Discussion:
            //
            //    This routine computes what would normally be called a Hermite 
            //    interpolant.  However, the user is only required to supply function
            //    values, not derivative values as well.  This routine computes
            //    "suitable" derivative values, so that the resulting Hermite interpolant
            //    has desirable shape and monotonicity properties.
            //
            //    The interpolant will have an extremum at each point where
            //    monotonicity switches direction.
            //
            //    The resulting piecewise cubic Hermite function may be evaluated
            //    by SPLINE_PCHIP_VAL..
            //
            //    This routine was originally called "PCHIM".
            //
            //    An "abs" was corrected to a "Math.Abs" on the report of Thomas Beutlich,
            //    10 October 2012.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 August 2005
            //
            //  Author:
            //
            //    FORTRAN77 original version by Fred Fritsch, Lawrence Livermore National Laboratory.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Fred Fritsch, Ralph Carlson,
            //    Monotone Piecewise Cubic Interpolation,
            //    SIAM Journal on Numerical Analysis,
            //    Volume 17, Number 2, April 1980, pages 238-246.
            //
            //    Fred Fritsch, Judy Butland,
            //    A Method for Constructing Local Monotone Piecewise 
            //    Cubic Interpolants,
            //    SIAM Journal on Scientific and Statistical Computing,
            //    Volume 5, Number 2, 1984, pages 300-304.
            //
            //  Parameters:
            //
            //    Input, int N, the number of data points.  N must be at least 2.
            //
            //    Input, double X[N], the strictly increasing independent
            //    variable values.
            //
            //    Input, double F[N], dependent variable values to be interpolated.  This 
            //    routine is designed for monotonic data, but it will work for any F-array.
            //    It will force extrema at points where monotonicity switches direction.
            //
            //    Output, double D[N], the derivative values at the
            //    data points.  If the data are monotonic, these values will determine
            //    a monotone cubic Hermite function.  
            //
        {
            double del1;
            double del2;
            double dmax;
            double dmin;
            double drat1;
            double drat2;
            double dsave;
            double h1;
            double h2;
            double hsum;
            double hsumt3;
            int i;
            int ierr;
            int nless1;
            double temp;
            double w1;
            double w2;
            //
            //  Check the arguments.
            //
            if (n < 2)
            {
                ierr = -1;
                Console.WriteLine("");
                Console.WriteLine("SPLINE_PCHIP_SET - Fatal error!");
                Console.WriteLine("  Number of data points less than 2.");
                return;
            }

            for (i = 1; i < n; i++)
            {
                if (x[i] <= x[i - 1])
                {
                    ierr = -3;
                    Console.WriteLine("");
                    Console.WriteLine("SPLINE_PCHIP_SET - Fatal error!");
                    Console.WriteLine("  X array not strictly increasing.");
                    return;
                }
            }

            ierr = 0;
            nless1 = n - 1;
            h1 = x[1] - x[0];
            del1 = (f[1] - f[0]) / h1;
            dsave = del1;
            //
            //  Special case N=2, use linear interpolation.
            //
            if (n == 2)
            {
                d[0] = del1;
                d[n - 1] = del1;
                return;
            }

            //
            //  Normal case, 3 <= N.
            //
            h2 = x[2] - x[1];
            del2 = (f[2] - f[1]) / h2;
            //
            //  Set D(1) via non-centered three point formula, adjusted to be
            //  shape preserving.
            //
            hsum = h1 + h2;
            w1 = (h1 + hsum) / hsum;
            w2 = -h1 / hsum;
            d[0] = w1 * del1 + w2 * del2;

            if (pchst(d[0], del1) <= 0.0)
            {
                d[0] = 0.0;
            }
            //
            //  Need do this check only if monotonicity switches.
            //
            else if (pchst(del1, del2) < 0.0)
            {
                dmax = 3.0 * del1;

                if (Math.Abs(dmax) < Math.Abs(d[0]))
                {
                    d[0] = dmax;
                }

            }

            //
            //  Loop through interior points.
            //
            for (i = 2; i <= nless1; i++)
            {
                if (2 < i)
                {
                    h1 = h2;
                    h2 = x[i] - x[i - 1];
                    hsum = h1 + h2;
                    del1 = del2;
                    del2 = (f[i] - f[i - 1]) / h2;
                }

                //
                //  Set D(I)=0 unless data are strictly monotonic.
                //
                d[i - 1] = 0.0;

                temp = pchst(del1, del2);

                if (temp < 0.0)
                {
                    ierr = ierr + 1;
                    dsave = del2;
                }
                //
                //  Count number of changes in direction of monotonicity.
                //
                else if (temp == 0.0)
                {
                    if (del2 != 0.0)
                    {
                        if (pchst(dsave, del2) < 0.0)
                        {
                            ierr = ierr + 1;
                        }

                        dsave = del2;
                    }
                }
                //
                //  Use Brodlie modification of Butland formula.
                //
                else
                {
                    hsumt3 = 3.0 * hsum;
                    w1 = (hsum + h1) / hsumt3;
                    w2 = (hsum + h2) / hsumt3;
                    dmax = Math.Max(Math.Abs(del1), Math.Abs(del2));
                    dmin = Math.Min(Math.Abs(del1), Math.Abs(del2));
                    drat1 = del1 / dmax;
                    drat2 = del2 / dmax;
                    d[i - 1] = dmin / (w1 * drat1 + w2 * drat2);
                }
            }

            //
            //  Set D(N) via non-centered three point formula, adjusted to be
            //  shape preserving.
            //
            w1 = -h2 / hsum;
            w2 = (h2 + hsum) / hsum;
            d[n - 1] = w1 * del1 + w2 * del2;

            if (pchst(d[n - 1], del2) <= 0.0)
            {
                d[n - 1] = 0.0;
            }
            else if (pchst(del1, del2) < 0.0)
            {
                //
                //  Need do this check only if monotonicity switches.
                //
                dmax = 3.0 * del2;

                if (Math.Abs(dmax) < Math.Abs(d[n - 1]))
                {
                    d[n - 1] = dmax;
                }

            }

            return;
        }

        public static void spline_pchip_val(int n, double[] x, double[] f, double[] d,
                int ne, double[] xe, ref double[] fe)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPLINE_PCHIP_VAL evaluates a piecewise cubic Hermite function.
            //
            //  Description:
            //
            //    This routine may be used by itself for Hermite interpolation, or as an
            //    evaluator for SPLINE_PCHIP_SET.
            //
            //    This routine evaluates the cubic Hermite function at the points XE.
            //
            //    Most of the coding between the call to CHFEV and the end of
            //    the IR loop could be eliminated if it were permissible to
            //    assume that XE is ordered relative to X.
            //
            //    CHFEV does not assume that X1 is less than X2.  Thus, it would
            //    be possible to write a version of SPLINE_PCHIP_VAL that assumes a strictly
            //    decreasing X array by simply running the IR loop backwards
            //    and reversing the order of appropriate tests.
            //
            //    The present code has a minor bug, which I have decided is not
            //    worth the effort that would be required to fix it.
            //    If XE contains points in [X(N-1),X(N)], followed by points less than
            //    X(N-1), followed by points greater than X(N), the extrapolation points
            //    will be counted (at least) twice in the total returned in IERR.
            //
            //    The evaluation will be most efficient if the elements of XE are
            //    increasing relative to X; that is, for all J <= K,
            //      X(I) <= XE(J)
            //    implies
            //      X(I) <= XE(K).
            //
            //    If any of the XE are outside the interval [X(1),X(N)],
            //    values are extrapolated from the nearest extreme cubic,
            //    and a warning error is returned.
            //
            //    This routine was originally named "PCHFE".
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 August 2005
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
            //  Parameters:
            //
            //    Input, int N, the number of data points.  N must be at least 2.
            //
            //    Input, double X[N], the strictly increasing independent
            //    variable values.
            //
            //    Input, double F[N], the function values.
            //
            //    Input, double D[N], the derivative values.
            //
            //    Input, int NE, the number of evaluation points.
            //
            //    Input, double XE[NE], points at which the function is to
            //    be evaluated.
            //
            //    Output, double FE[NE], the values of the cubic Hermite
            //    function at XE.
            //
        {
            int i;
            int ierc;
            int ierr;
            int ir;
            int j;
            int j_first;
            int j_new;
            int j_save;
            int[] next = new int[2];
            int nj;
            //
            //  Check arguments.
            //
            if (n < 2)
            {
                ierr = -1;
                Console.WriteLine("");
                Console.WriteLine("SPLINE_PCHIP_VAL - Fatal error!");
                Console.WriteLine("  Number of data points less than 2.");
                return;
            }

            for (i = 1; i < n; i++)
            {
                if (x[i] <= x[i - 1])
                {
                    ierr = -3;
                    Console.WriteLine("");
                    Console.WriteLine("SPLINE_PCHIP_VAL - Fatal error!");
                    Console.WriteLine("  X array not strictly increasing.");
                    return;
                }
            }

            if (ne < 1)
            {
                ierr = -4;
                Console.WriteLine("");
                Console.WriteLine("SPLINE_PCHIP_VAL - Fatal error!");
                Console.WriteLine("  Number of evaluation points less than 1.");
                return;
            }

            ierr = 0;
            //
            //  Loop over intervals.
            //  The interval index is IL = IR-1.
            //  The interval is X(IL) <= X < X(IR).
            //
            j_first = 1;
            ir = 2;

            for (;;)
            {
                //
                //  Skip out of the loop if have processed all evaluation points.
                //
                if (ne < j_first)
                {
                    break;
                }

                //
                //  Locate all points in the interval.
                //
                j_save = ne + 1;

                for (j = j_first; j <= ne; j++)
                {
                    if (x[ir - 1] <= xe[j - 1])
                    {
                        j_save = j;
                        if (ir == n)
                        {
                            j_save = ne + 1;
                        }

                        break;
                    }
                }

                //
                //  Have located first point beyond interval.
                //
                j = j_save;

                nj = j - j_first;
                //
                //  Skip evaluation if no points in interval.
                //
                if (nj != 0)
                {
                    //
                    //  Evaluate cubic at XE(J_FIRST:J-1).
                    //
                    ierc = Burkardt.PolynomialNS.Cubic.chfev(x[ir - 2], x[ir - 1], f[ir - 2], f[ir - 1], d[ir - 2],
                        d[ir - 1],
                        nj, xe, ref fe, ref next, xeIndex: +j_first - 1, feIndex: +j_first - 1);

                    if (ierc < 0)
                    {
                        ierr = -5;
                        Console.WriteLine("");
                        Console.WriteLine("SPLINE_PCHIP_VAL - Fatal error!");
                        Console.WriteLine("  Error return from CHFEV.");
                        return;
                    }

                    //
                    //  In the current set of XE points, there are NEXT(2) to the right of X(IR).
                    //
                    if (next[1] != 0)
                    {
                        if (ir < n)
                        {
                            ierr = -5;
                            Console.WriteLine("");
                            Console.WriteLine("SPLINE_PCHIP_VAL - Fatal error!");
                            Console.WriteLine("  IR < N.");
                            return;
                        }

                        //
                        //  These are actually extrapolation points.
                        //
                        ierr = ierr + next[1];

                    }

                    //
                    //  In the current set of XE points, there are NEXT(1) to the left of X(IR-1).
                    //
                    if (next[0] != 0)
                    {
                        //
                        //  These are actually extrapolation points.
                        //
                        if (ir <= 2)
                        {
                            ierr = ierr + next[0];
                        }
                        else
                        {
                            j_new = -1;

                            for (i = j_first; i <= j - 1; i++)
                            {
                                if (xe[i - 1] < x[ir - 2])
                                {
                                    j_new = i;
                                    break;
                                }
                            }

                            if (j_new == -1)
                            {
                                ierr = -5;
                                Console.WriteLine("");
                                Console.WriteLine("SPLINE_PCHIP_VAL - Fatal error!");
                                Console.WriteLine("  Could not bracket the data point.");
                                return;
                            }

                            //
                            //  Reset J.  This will be the new J_FIRST.
                            //
                            j = j_new;
                            //
                            //  Now find out how far to back up in the X array.
                            //
                            for (i = 1; i <= ir - 1; i++)
                            {
                                if (xe[j - 1] < x[i - 1])
                                {
                                    break;
                                }
                            }

                            //
                            //  At this point, either XE(J) < X(1) or X(i-1) <= XE(J) < X(I) .
                            //
                            //  Reset IR, recognizing that it will be incremented before cycling.
                            //
                            ir = Math.Max(1, i - 1);
                        }
                    }

                    j_first = j;
                }

                ir = ir + 1;

                if (n < ir)
                {
                    break;
                }

            }
        }
    }
}