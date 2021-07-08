using System;

namespace Burkardt.PolynomialNS
{
    public static class Hermite
    {
        public static void hermite_recur ( ref double p2, ref double dp2, ref double p1, double x, int order )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_RECUR finds the value and derivative of a Hermite polynomial.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 April 2006
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Arthur Stroud, Don Secrest,
            //    Gaussian Quadrature Formulas,
            //    Prentice Hall, 1966,
            //    LC: QA299.4G3S7.
            //
            //  Parameters:
            //
            //    Output, double *P2, the value of H(ORDER)(X).
            //
            //    Output, double *DP2, the value of H'(ORDER)(X).
            //
            //    Output, double *P1, the value of H(ORDER-1)(X).
            //
            //    Input, double X, the point at which polynomials are evaluated.
            //
            //    Input, int ORDER, the order of the polynomial to be computed.
            //
        {
            int i;
            double dq0;
            double dq1;
            double dq2;
            double q0;
            double q1;
            double q2;

            q1 = 1.0;
            dq1 = 0.0;

            q2 = x;
            dq2 = 1.0;

            for ( i = 2; i <= order; i++ )
            {
                q0 = q1;
                dq0 = dq1;

                q1 = q2;
                dq1 = dq2;

                q2  = x * q1 - 0.5 * ( ( double ) ( i ) - 1.0 ) * q0;
                dq2 = x * dq1 + q1 - 0.5 * ( ( double ) ( i ) - 1.0 ) * dq0;
            }

            p2 = q2;
            dp2 = dq2;
            p1 = q1;
        }
        
        public static void hermite_root ( ref double x, int order, ref double dp2, ref double p1 )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_ROOT improves an approximate root of a Hermite polynomial.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 May 2006
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Arthur Stroud, Don Secrest,
            //    Gaussian Quadrature Formulas,
            //    Prentice Hall, 1966,
            //    LC: QA299.4G3S7.
            //
            //  Parameters:
            //
            //    Input/output, double *X, the approximate root, which
            //    should be improved on output.
            //
            //    Input, int ORDER, the order of the Hermite polynomial.
            //
            //    Output, double *DP2, the value of H'(ORDER)(X).
            //
            //    Output, double *P1, the value of H(ORDER-1)(X).
            //
        {
            double d;
            double eps = 1.0E-12;
            double p2 = 0;
            int step;
            int step_max = 10;

            for ( step = 1; step <= step_max; step++ )
            {
                hermite_recur ( ref p2, ref dp2, ref p1, x, order );

                d = p2 / ( dp2 );
                x = x - d;

                if ( Math.Abs ( d ) <= eps * ( Math.Abs ( x ) + 1.0 ) )
                {
                    return;
                }
            }
        }
    }
}