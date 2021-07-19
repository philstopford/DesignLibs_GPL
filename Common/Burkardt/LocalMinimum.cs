using System;

namespace Burkardt
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

    }
}