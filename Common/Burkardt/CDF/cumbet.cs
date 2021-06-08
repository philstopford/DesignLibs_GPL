namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        public static void cumbet ( double x, double y, double a, double b, ref double cum,
                ref double ccum )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CUMBET evaluates the cumulative incomplete beta distribution.
            //
            //  Discussion:
            //
            //    This routine calculates the CDF to X of the incomplete beta distribution
            //    with parameters A and B.  This is the integral from 0 to x
            //    of (1/B(a,b))*f(t)) where f(t) = t^(a-1) * (1-t)^(b-1)
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
            //  Reference:
            //
            //    A R Didonato and Alfred Morris,
            //    Algorithm 708:
            //    Significant Digit Computation of the Incomplete Beta Function Ratios.
            //    ACM Transactions on Mathematical Software,
            //    Volume 18, Number 3, September 1992, pages 360-373.
            //
            //  Input:
            //
            //    double *X, the upper limit of integration.
            //
            //    double *Y, the value of 1-X.
            //
            //    double *A, *B, the parameters of the distribution.
            //
            //  Output:
            //
            //    double *CUM, *CCUM, the values of the cumulative
            //    density function and complementary cumulative density function.
            //
        {
            int ierr = 0;

            if ( x <= 0.0 )
            {
                cum = 0.0;
                ccum = 1.0;
            }
            else if ( y <= 0.0 )
            {
                cum = 1.0;
                ccum = 0.0;
            }
            else
            {
                beta_inc ( a, b, x, y, ref cum, ref ccum, ref ierr );
            }
        }
    }
}