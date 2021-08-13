using System;

namespace Burkardt.SolveNS
{
    public static class LocalMinimum
    {
        public static double local_min(double a, double b, double t, Func<double, double> f,
                ref double x, ref int calls)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    local_min() seeks a local minimum of a function F(X) in an interval [A,B].
            //
            //  Discussion:
            //
            //    If the function F is defined on the interval (A,B), then local_min
            //    finds an approximation X to the point at which F attatains its minimum
            //    (or the appropriate limit point), and returns the value of F at X.
            //
            //    T and EPS define a tolerance TOL = EPS * abs ( X ) + T.
            //    F is never evaluated at two points closer than TOL.  
            //
            //    If F is delta-unimodal for some delta less than TOL, the X approximates
            //    the global minimum of F with an error less than 3*TOL.
            //
            //    If F is not delta-unimodal, then X may approximate a local, but 
            //    perhaps non-global, minimum.
            //
            //    The method used is a combination of golden section search and
            //    successive parabolic interpolation.  Convergence is never much slower
            //    than that for a Fibonacci search.  If F has a continuous second
            //    derivative which is positive at the minimum (which is not at A or
            //    B), then, ignoring rounding errors, convergence is superlinear, and 
            //    usually of the order of about 1.3247.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 May 2021
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Richard Brent.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Richard Brent,
            //    Algorithms for Minimization Without Derivatives,
            //    Dover, 2002,
            //    ISBN: 0-486-41998-3,
            //    LC: QA402.5.B74.
            //
            //  Input:
            //
            //    double A, B, the endpoints of the interval.
            //
            //    double T, a positive absolute error tolerance.
            //
            //    double F(double x), the name of a user-supplied function.
            //
            //  Output:
            //
            //    double &X, the estimated value of an abscissa
            //    for which F attains a local minimum value in [A,B].
            //
            //    int &CALLS: the number of calls to F.
            //
            //    double LOCAL_MIN, the value F(X).
            //
        {
            double c;
            double d = 0;
            double e;
            double eps;
            double fu;
            double fv;
            double fw;
            double fx;
            double m;
            double p;
            double q;
            double r;
            double sa;
            double sb;
            double t2;
            double tol;
            double u;
            double v;
            double w;

            calls = 0;
            //
            //  C is the square of the inverse of the golden ratio.
            //
            c = 0.5 * (3.0 - Math.Sqrt(5.0));

            eps = Math.Sqrt(double.Epsilon);

            sa = a;
            sb = b;
            x = sa + c * (b - a);
            w = x;
            v = w;
            e = 0.0;
            fx = f(x);
            calls = calls + 1;
            fw = fx;
            fv = fw;

            for (;;)
            {
                m = 0.5 * (sa + sb);
                tol = eps * Math.Abs(x) + t;
                t2 = 2.0 * tol;
                //
                //  Check the stopping criterion.
                //
                if (Math.Abs(x - m) <= t2 - 0.5 * (sb - sa))
                {
                    break;
                }

                //
                //  Fit a parabola.
                //
                r = 0.0;
                q = r;
                p = q;

                if (tol < Math.Abs(e))
                {
                    r = (x - w) * (fx - fv);
                    q = (x - v) * (fx - fw);
                    p = (x - v) * q - (x - w) * r;
                    q = 2.0 * (q - r);
                    if (0.0 < q)
                    {
                        p = -p;
                    }

                    q = Math.Abs(q);
                    r = e;
                    e = d;
                }

                if (Math.Abs(p) < Math.Abs(0.5 * q * r) &&
                    q * (sa - x) < p &&
                    p < q * (sb - x))
                {
                    //
                    //  Take the parabolic interpolation step.
                    //
                    d = p / q;
                    u = x + d;
                    //
                    //  F must not be evaluated too close to A or B.
                    //
                    if ((u - sa) < t2 || (sb - u) < t2)
                    {
                        if (x < m)
                        {
                            d = tol;
                        }
                        else
                        {
                            d = -tol;
                        }
                    }
                }
                //
                //  A golden-section step.
                //
                else
                {
                    if (x < m)
                    {
                        e = sb - x;
                    }
                    else
                    {
                        e = sa - x;
                    }

                    d = c * e;
                }

                //
                //  F must not be evaluated too close to X.
                //
                if (tol <= Math.Abs(d))
                {
                    u = x + d;
                }
                else if (0.0 < d)
                {
                    u = x + tol;
                }
                else
                {
                    u = x - tol;
                }

                fu = f(u);
                calls = calls + 1;
                //
                //  Update A, B, V, W, and X.
                //
                if (fu <= fx)
                {
                    if (u < x)
                    {
                        sb = x;
                    }
                    else
                    {
                        sa = x;
                    }

                    v = w;
                    fv = fw;
                    w = x;
                    fw = fx;
                    x = u;
                    fx = fu;
                }
                else
                {
                    if (u < x)
                    {
                        sa = u;
                    }
                    else
                    {
                        sb = u;
                    }

                    if (fu <= fw || w == x)
                    {
                        v = w;
                        fv = fw;
                        w = u;
                        fw = fu;
                    }
                    else if (fu <= fv || v == x || v == w)
                    {
                        v = u;
                        fv = fu;
                    }
                }
            }

            return fx;
        }

        public class LocalMinimumData
        {
            public double arg;
            public double c;
            public double d;
            public double e;
            public double eps;
            public double fu;
            public double fv;
            public double fw;
            public double fx;
            public double midpoint;
            public double p;
            public double q;
            public double r;
            public double tol;
            public double tol1;
            public double tol2;
            public double u;
            public double v;
            public double w;
            public double x;
            
        }

        public static double local_min_rc (ref LocalMinimumData data, ref double a, ref double b, ref int status, double value )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    local_min_rc() seeks a minimizer of a scalar function of a scalar variable.
        //
        //  Discussion:
        //
        //    This routine seeks an approximation to the point where a function
        //    F attains a minimum on the interval (A,B).
        //
        //    The method used is a combination of golden section search and
        //    successive parabolic interpolation.  Convergence is never much
        //    slower than that for a Fibonacci search.  If F has a continuous
        //    second derivative which is positive at the minimum (which is not
        //    at A or B), then convergence is superlinear, and usually of the
        //    order of about 1.324...
        //
        //    The routine is a revised version of the Brent local minimization
        //    algorithm, using reverse communication.
        //
        //    It is worth stating explicitly that this routine will NOT be
        //    able to detect a minimizer that occurs at either initial endpoint
        //    A or B.  If this is a concern to the user, then the user must
        //    either ensure that the initial interval is larger, or to check
        //    the function value at the returned minimizer against the values
        //    at either endpoint.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 May 2021
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Richard Brent,
        //    Algorithms for Minimization Without Derivatives,
        //    Dover, 2002,
        //    ISBN: 0-486-41998-3,
        //    LC: QA402.5.B74.
        //
        //    David Kahaner, Cleve Moler, Steven Nash,
        //    Numerical Methods and Software,
        //    Prentice Hall, 1989,
        //    ISBN: 0-13-627258-4,
        //    LC: TA345.K34.
        //
        //  Parameters
        //
        //    Input/output, double &A, &B.  On input, the left and right
        //    endpoints of the initial interval.  On output, the lower and upper
        //    bounds for an interval containing the minimizer.  It is required
        //    that A < B.
        //
        //    Input/output, int &STATUS, used to communicate between
        //    the user and the routine.  The user only sets STATUS to zero on the first
        //    call, to indicate that this is a startup call.  The routine returns STATUS
        //    positive to request that the function be evaluated at ARG, or returns
        //    STATUS as 0, to indicate that the iteration is complete and that
        //    ARG is the estimated minimizer.
        //
        //    Input, double VALUE, the function value at ARG, as requested
        //    by the routine on the previous call.
        //
        //    Output, double LOCAL_MIN_RC, the currently considered point.
        //    On return with STATUS positive, the user is requested to evaluate the
        //    function at this point, and return the value in VALUE.  On return with
        //    STATUS zero, this is the routine's estimate for the function minimizer.
        //
        //  Local:
        //
        //    double C: the squared inverse of the golden ratio.
        //
        //    double EPS: the square root of the relative machine precision.
        //
        {
        //
        //  STATUS (INPUT) = 0, startup.
        //
        if ( status == 0 )
        {
        if ( b <= a )
        {
        Console.WriteLine("");
        Console.WriteLine("local_min_rc(): Fatal error!");
        Console.WriteLine("  A < B is required, but");
        Console.WriteLine("  A = " + a + "");
        Console.WriteLine("  B = " + b + "");
        status = -1;
        return ( 1 );
        }
        data.c = 0.5 * ( 3.0 - Math.Sqrt ( 5.0 ) );

        data.eps = Math.Sqrt ( double.Epsilon );
        data.tol = double.Epsilon;

        data.v = a + data.c * ( b - a );
        data.w = data.v;
        data.x = data.v;
        data.e = 0.0;

        status = 1;
        data.arg = data.x;

        return data.arg;
        }
        //
        //  STATUS (INPUT) = 1, return with initial function value of FX.
        //
        else if ( status == 1 )
        {
            data.fx = value;
            data.fv = data.fx;
        data.fw = data.fx;
        }
        //
        //  STATUS (INPUT) = 2 or more, update the data.
        //
        else if ( 2 <= status )
        {
            data.fu = value;

        if ( data.fu <= data.fx )
        {
        if ( data.x <= data.u )
        {
        a = data.x;
        }
        else
        {
        b = data.x;
        }
        data.v = data.w;
        data.fv = data.fw;
        data.w = data.x;
        data.fw = data.fx;
        data.x = data.u;
        data.fx = data.fu;
        }
        else
        {
        if ( data.u < data.x )
        {
        a = data.u;
        }
        else
        {
        b = data.u;
        }

        if ( data.fu <= data.fw || data.w == data.x )
        {
            data.v = data.w;
            data.fv = data.fw;
            data.w = data.u;
            data.fw = data.fu;
        }
        else if ( data.fu <= data.fv || data.v == data.x || data.v == data.w )
        {
            data.v = data.u;
            data.fv = data.fu;
        }
        }
        }
        //
        //  Take the next step.
        //
        data.midpoint = 0.5 * ( a + b );
        data.tol1 = data.eps * Math.Abs ( data.x ) + data.tol / 3.0;
        data.tol2 = 2.0 * data.tol1;
        //
        //  If the stopping criterion is satisfied, we can exit.
        //
        if ( Math.Abs ( data.x - data.midpoint ) <= ( data.tol2 - 0.5 * ( b - a ) ) )
        {
        status = 0;
        return data.arg;
        }
        //
        //  Is golden-section necessary?
        //
        if ( Math.Abs ( data.e ) <= data.tol1 )
        {
        if ( data.midpoint <= data.x )
        {
            data.e = a - data.x;
        }
        else
        {
            data.e = b - data.x;
        }
        data.d = data.c * data.e;
        }
        //
        //  Consider fitting a parabola.
        //
        else
        {
            data.r = ( data.x - data.w ) * ( data.fx - data.fv );
            data.q = ( data.x - data.v ) * ( data.fx - data.fw );
            data.p = ( data.x - data.v ) * data.q - ( data.x - data.w ) * data.r;
            data.q = 2.0 * ( data.q - data.r );
        if ( 0.0 < data.q )
        {
            data.p = - data.p;
        }
        data.q = Math.Abs ( data.q );
        data.r = data.e;
        data.e = data.d;
        //
        //  Choose a golden-section step if the parabola is not advised.
        //
        if (
        ( Math.Abs ( 0.5 * data.q * data.r ) <= Math.Abs ( data.p ) ) ||
        ( data.p <= data.q * ( a - data.x ) ) ||
        ( data.q * ( b - data.x ) <= data.p ) )
        {
        if ( data.midpoint <= data.x )
        {
            data.e = a - data.x;
        }
        else
        {
            data.e = b - data.x;
        }
        data.d = data.c * data.e;
        }
        //
        //  Choose a parabolic interpolation step.
        //
        else
        {
            data.d = data.p / data.q;
            data.u = data.x + data.d;

        if ( ( data.u - a ) < data.tol2 )
        {
            data.d = Math.CopySign ( data.tol1, data.midpoint - data.x );
        }

        if ( ( b - data.u ) < data.tol2 )
        {
            data.d = Math.CopySign ( data.tol1, data.midpoint - data.x );
        }
        }
        }
        //
        //  F must not be evaluated too close to X.
        //
        if ( data.tol1 <= Math.Abs ( data.d ) )
        {
            data.u = data.x + data.d;
        }
        if ( Math.Abs ( data.d ) < data.tol1 )
        {
            data.u = data.x + Math.CopySign ( data.tol1, data.d );
        }
        //
        //  Request value of F(U).
        //
        data.arg = data.u;
        status = status + 1;

        return data.arg;
        }

    }
}