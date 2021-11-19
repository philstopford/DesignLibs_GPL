using System;
using Burkardt.Types;

namespace Burkardt.SolveNS;

public static class Zero
{
    public static double zero(double a, double b, double t, Func<double, double> f, ref int calls)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    zero() seeks a root of a function F(X) in an interval [A,B].
        //
        //  Discussion:
        //
        //    The interval [A,B] must be a change of sign interval for F.
        //    That is, F(A) and F(B) must be of opposite signs.  Then
        //    assuming that F is continuous implies the existence of at least
        //    one value C between A and B for which F(C) = 0.
        //
        //    The location of the zero is determined to within an accuracy
        //    of 4 * EPSILON * abs ( C ) + 2 * T.
        //
        //    Thanks to Thomas Secretin for pointing out a transcription error in the
        //    setting of the value of P, 11 February 2013.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 July 2021
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
        //    double A, B, the endpoints of the change of sign interval.
        //
        //    double T, a positive error tolerance.
        //
        //    double F(double x), the name of a user-supplied function.
        //
        //  Output:
        //
        //    int &CALLS, the number of calls to the function F.
        //
        //    double ZERO, the estimated value of a zero of the function F.
        //
    {
        calls = 0;
        //
        //  Make local copies of A and B.
        //
        double sa = a;
        double fa = f(sa);
        calls += 1;

        double sb = b;
        double fb = f(sb);
        calls += 1;

        double c = sa;
        double fc = fa;
        double e = sb - sa;
        double d = e;

        for (;;)
        {
            if (Math.Abs(fc) < Math.Abs(fb))
            {
                sa = sb;
                sb = c;
                c = sa;
                fa = fb;
                fb = fc;
                fc = fa;
            }

            double tol = 2.0 * typeMethods.r8_epsilon() * Math.Abs(sb) + t;
            double m = 0.5 * (c - sb);

            if (Math.Abs(m) <= tol || fb == 0.0)
            {
                break;
            }

            if (Math.Abs(e) < tol || Math.Abs(fa) <= Math.Abs(fb))
            {
                e = m;
                d = e;
            }
            else
            {
                double s = fb / fa;

                double q;
                double p;
                if (Math.Abs(sa - c) <= double.Epsilon)
                {
                    p = 2.0 * m * s;
                    q = 1.0 - s;
                }
                else
                {
                    q = fa / fc;
                    double r = fb / fc;
                    p = s * (2.0 * m * q * (q - r) - (sb - sa) * (r - 1.0));
                    q = (q - 1.0) * (r - 1.0) * (s - 1.0);
                }

                switch (p)
                {
                    case > 0.0:
                        q = -q;
                        break;
                    default:
                        p = -p;
                        break;
                }

                s = e;
                e = d;

                if (2.0 * p < 3.0 * m * q - Math.Abs(tol * q) &&
                    p < Math.Abs(0.5 * s * q))
                {
                    d = p / q;
                }
                else
                {
                    e = m;
                    d = e;
                }
            }

            sa = sb;
            fa = fb;

            if (tol < Math.Abs(d))
            {
                sb += d;
            }
            else
            {
                switch (m)
                {
                    case > 0.0:
                        sb += tol;
                        break;
                    default:
                        sb -= tol;
                        break;
                }
            }

            fb = f(sb);
            calls += 1;

            switch (fb)
            {
                case > 0.0 when 0.0 < fc:
                case <= 0.0 when fc <= 0.0:
                    c = sa;
                    fc = fa;
                    e = sb - sa;
                    d = e;
                    break;
            }
        }

        return sb;
    }

}