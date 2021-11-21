using System;

namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static double apser ( double a, double b, double x, double eps )

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
        double apser = 0;
        const double g = 0.577215664901533e0;

        double bx = b * x;
        double t = x - bx;

        double c = (b * eps) switch
        {
            <= 2e-2 => Math.Log(x) + psi(b) + g + t,
            _ => Math.Log(bx) + g + t
        };

        double tol = 5.0e0 * eps * Math.Abs ( c );
        double j = 1.0e0;
        double s = 0.0e0;

        while ( true )
        {
            j += 1.0e0;
            t *= ( x - bx / j );
            double aj = t / j;
            s += aj;
            if ( Math.Abs ( aj ) <= tol )
            {
                break;
            }
            apser = - ( a * ( c + s ) );
        }

        return apser;
    }
}