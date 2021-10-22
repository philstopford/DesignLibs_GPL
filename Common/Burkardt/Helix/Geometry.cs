using System;

namespace Burkardt.Helix
{
    public static class Geometry
    {
        public static void helix_shape_3d ( double a, int n, double r, double theta1, double theta2,
                ref double[] p )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HELIX_SHAPE_3D computes points on a helix in 3D.
            //
            //  Discussion:
            //
            //    The user specifies the parameters A and R, the first and last
            //    THETA values, and the number of equally spaced THETA values
            //    at which values are to be computed.
            //
            //      ( R * COS ( THETA ), R * SIN ( THETA ), A * THETA )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 September 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double A, the rate at which Z advances with THETA.
            //
            //    Input, int N, the number of points to compute on the helix.
            //
            //    Input, double R, the radius of the helix.
            //
            //    Input, double THETA1, THETA2, the first and last THETA values at
            //    which to compute points on the helix.  THETA is measured in
            //    radians.
            //
            //    Output, double P[3*N], the coordinates of points on the helix.
            //
        {
            int i;
            double theta;

            for ( i = 0; i < n; i++ )
            {
                if ( n == 1 )
                {
                    theta = 0.5 * ( theta1 + theta2 );
                }
                else
                {
                    theta = ( ( ( double ) ( n - i     ) ) * theta1
                              + ( ( double ) (     i - 1 ) ) * theta2 )
                            / ( ( double ) ( n     - 1 ) );
                }

                p[0+i*3] = r * Math.Cos ( theta );
                p[1+i*3] = r * Math.Sin ( theta );
                p[2+i*3] = a * theta;
            }
        }
    }
}