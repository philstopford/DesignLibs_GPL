using System;

namespace Burkardt.ZeroRCNS
{
    public class ZeroRC_data
    {
        public double c;
        public double d;
        public double e;
        public double fa;
        public double fb;
        public double fc;
        public double sa;
        public double sb;
    }

    public static class ZeroRC
    {
        public static void zero_rc(ref ZeroRC_data data, double a, double b, double t, ref double arg, ref int status,
        double value )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    zero_rc() seeks a root of a function F(X) using reverse communication.
        //
        //  Discussion:
        //
        //    The interval [A,B] must be a change of sign interval for F.
        //    That is, F(A) and F(B) must be of opposite signs.  Then
        //    assuming that F is continuous implies the existence of at least
        //    one value C between A and B for which F(C) = 0.
        //
        //    The location of the zero is determined to within an accuracy
        //    of 6 * EPSILON * r8_abs ( C ) + 2 * T.
        //
        //    The routine is a revised version of the Brent zero finder
        //    algorithm, using reverse communication.
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
        //  Parameters:
        //
        //    Input, double A, B, the endpoints of the change of sign interval.
        //
        //    Input, double T, a positive error tolerance.
        //
        //    Output, double &ARG, the currently considered point.  The user
        //    does not need to initialize this value.  On return with STATUS positive,
        //    the user is requested to evaluate the function at ARG, and return
        //    the value in VALUE.  On return with STATUS zero, ARG is the routine's
        //    estimate for the function's zero.
        //
        //    Input/output, int &STATUS, used to communicate between
        //    the user and the routine.  The user only sets STATUS to zero on the first
        //    call, to indicate that this is a startup call.  The routine returns STATUS
        //    positive to request that the function be evaluated at ARG, or returns
        //    STATUS as 0, to indicate that the iteration is complete and that
        //    ARG is the estimated zero
        //
        //    Input, double VALUE, the function value at ARG, as requested
        //    by the routine on the previous call.
        //
        {
            double m;
            double p;
            double q;
            double r;
            double s;
            double tol;
            //
            //  Input STATUS = 0.
            //  Initialize, request F(A).
            //
            if (status == 0)
            {
                data.sa = a;
                data.sb = b;
                data.e = data.sb - data.sa;
                data.d = data.e;

                status = 1;
                arg = a;
                return;
            }
            //
            //  Input STATUS = 1.
            //  Receive F(A), request F(B).
            //
            else if (status == 1)
            {
                data.fa = value;
                status = 2;
                arg = data.sb;
                return;
            }
            //
            //  Input STATUS = 2
            //  Receive F(B).
            //
            else if (status == 2)
            {
                data.fb = value;

                if (0.0 < data.fa * data.fb)
                {
                    status = -1;
                    return;
                }

                data.c = data.sa;
                data.fc = data.fa;
            }
            else
            {
                data.fb = value;

                if ((0.0 < data.fb && 0.0 < data.fc) || (data.fb <= 0.0 && data.fc <= 0.0))
                {
                    data.c = data.sa;
                    data.fc = data.fa;
                    data.e = data.sb - data.sa;
                    data.d = data.e;
                }
            }

            //
            //  Compute the next point at which a function value is requested.
            //
            if (Math.Abs(data.fc) < Math.Abs(data.fb))
            {
                data.sa = data.sb;
                data.sb = data.c;
                data.c = data.sa;
                data.fa = data.fb;
                data.fb = data.fc;
                data.fc = data.fa;
            }

            tol = 2.0 * double.Epsilon * Math.Abs(data.sb) + t;
            m = 0.5 * (data.c - data.sb);

            if (Math.Abs(m) <= tol || data.fb == 0.0)
            {
                status = 0;
                arg = data.sb;
                return;
            }

            if (Math.Abs(data.e) < tol || Math.Abs(data.fa) <= Math.Abs(data.fb))
            {
                data.e = m;
                data.d = data.e;
            }
            else
            {
                s = data.fb / data.fa;

                if (data.sa == data.c)
                {
                    p = 2.0 * m * s;
                    q = 1.0 - s;
                }
                else
                {
                    q = data.fa / data.fc;
                    r = data.fb / data.fc;
                    p = s * (2.0 * m * q * (q - r) - (data.sb - data.sa) * (r - 1.0));
                    q = (q - 1.0) * (r - 1.0) * (s - 1.0);
                }

                if (0.0 < p)
                {
                    q = -q;
                }
                else
                {
                    p = -p;
                }

                s = data.e;
                data.e = data.d;

                if (2.0 * p < 3.0 * m * q - Math.Abs(tol * q) &&
                    p < Math.Abs(0.5 * s * q))
                {
                    data.d = p / q;
                }
                else
                {
                    data.e = m;
                    data.d = data.e;
                }
            }

            data.sa = data.sb;
            data.fa = data.fb;

            if (tol < Math.Abs(data.d))
            {
                data.sb = data.sb + data.d;
            }
            else if (0.0 < m)
            {
                data.sb = data.sb + tol;
            }
            else
            {
                data.sb = data.sb - tol;
            }

            arg = data.sb;
            status = status + 1;
        }
    }
}