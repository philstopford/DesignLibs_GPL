using System;

namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        public static double apser(double a, double b, double x, double eps)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    APSER computes the incomplete beta ratio I(SUB(1-X))(B,A).
            //
            //  Discussion:
            //
            //    APSER is used only for cases where
            //
            //      A <= min ( EPS, EPS * B ),
            //      B * X <= 1, and
            //      X <= 0.5.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 February 2021
            //
            //  Author:
            //
            //    Barry Brown, James Lovato, Kathy Russell.
            //
            //  Input:
            //
            //    double *A, *B, *X, the parameters of the incomplete beta ratio.
            //
            //    double *EPS, a tolerance.
            //
            //  Output:
            //
            //    double APSER, the value of the incomplete beta ratio.
            //
        {
            double aj;
            double apser = 0;
            double bx;
            double c;
            double g = 0.577215664901533e0;
            double j;
            double s;
            double t;
            double tol;

            bx = b * x;
            t = x - bx;

            if (b * eps <= 2e-2)
            {
                c = Math.Log(x) + psi(b) + g + t;
            }
            else
            {
                c = Math.Log(bx) + g + t;
            }

            tol = 5.0e0 * eps * Math.Abs(c);
            j = 1.0e0;
            s = 0.0e0;

            while (true)
            {
                j = j + 1.0e0;
                t = t * (x - bx / j);
                aj = t / j;
                s = s + aj;
                if (Math.Abs(aj) <= tol)
                {
                    break;
                }

                apser = -(a * (c + s));
            }

            return apser;
        }
    }
}