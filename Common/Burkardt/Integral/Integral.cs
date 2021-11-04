using System;

namespace Burkardt.IntegralNS
{
    public static partial class Integral
    {
        public static double fu_integral ( int d )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FU_INTEGRAL is the integral of the test function for the [0,1]^D interval.
            //
            //  Discussion:
            //
            //    The same function, integrated over [-1,+1]^D, has an integral
            //    that is 2^D times larger.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    08 December 2012
            //
            //  Author:
            //
            //    John Burkardt.
            //
            //  Parameters:
            //
            //    Input, int D, the spatial dimension.
            //
            //    Output, double FU_INTEGRAL, the integral value.
            //
        {
            double value;

            value = Math.Pow ( 0.5 * Helpers.Erf ( 0.5 / Math.Sqrt ( 2.0 ) ), d );

            return value;
        }

        public static double[] fu_value ( int d, int n, double[] x )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FU_VALUE is a sample function for the [0,1]^D interval.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    08 December 2012
            //
            //  Author:
            //
            //    John Burkardt.
            //
            //  Parameters:
            //
            //    Input, int D, the spatial dimension.
            //
            //    Input, int N, the number of points.
            //
            //    Input, double X[D*N], the points.
            //
            //    Output, double FU_VALUE[N], the function values.
            //
        {
            double[] fx;
            int i;
            int j;
            

            fx = new double[n];

            for ( j = 0; j < n; j++ )
            {
                fx[j] = 1.0;
                for ( i = 0; i < d; i++ )
                {
                    fx[j] = fx[j] * Math.Exp ( - Math.Pow ( x[i+j*d] / 2.0, 2 ) / 2.0 ) 
                            / 2.0 / Math.Sqrt ( 2.0 * Math.PI );
                }
            }

            return fx;
        }
    }
}