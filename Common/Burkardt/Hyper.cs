using System;

namespace Burkardt
{
    public static class Hyper
    {
        public static double hypersphere_unit_volume ( int m )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HYPERSPHERE_UNIT_VOLUME: volume of a unit hypersphere in M dimensions.
            //
            //  Discussion:
            //
            //    The unit sphere in M dimensions satisfies the equation:
            //
            //      Sum ( 1 <= I <= M ) X(I) * X(I) = 1
            //
            //     M  Volume
            //
            //     1    2
            //     2    1        * PI
            //     3  ( 4 /   3) * PI
            //     4  ( 1 /   2) * PI^2
            //     5  ( 8 /  15) * PI^2
            //     6  ( 1 /   6) * PI^3
            //     7  (16 / 105) * PI^3
            //     8  ( 1 /  24) * PI^4
            //     9  (32 / 945) * PI^4
            //    10  ( 1 / 120) * PI^5
            //
            //    For the unit sphere, Volume(M) = 2 * PI * Volume(M-2)/ M
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 August 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the dimension of the space.
            //
            //    Output, double HYPERSPHERE_UNIT_VOLUME, the volume of the sphere.
            //
        {
            int i;
            int m2;
            const double r8_pi = 3.141592653589793;
            double volume;

            if ( m % 2== 0 )
            {
                m2 = m / 2;
                volume = 1.0;
                for ( i = 1; i <= m2; i++ )
                {
                    volume = volume * r8_pi / ( ( double ) i );
                }
            }
            else
            {
                m2 = ( m - 1 ) / 2;
                volume = Math.Pow ( r8_pi, m2 ) * Math.Pow ( 2.0, m );
                for ( i = m2 + 1; i <= 2 * m2 + 1; i++ )
                {
                    volume = volume / ( ( double ) i );
                }
            }

            return volume;
        }
    }
}