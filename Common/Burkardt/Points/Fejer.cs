using System;

namespace Burkardt.PointsNS
{
    public static partial class Points
    {
        public static double[] fejer1 ( int n )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FEJER1 returns the Type 1 Fejer points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 March 2018
            //
            //  Author:
            //
            //    John Burkardt.
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, double FEJER1[N], the points.
            //
        {
            int i;
            
            double theta;
            double[] x;

            x = new double[n];

            for ( i = 0; i < n; i++ )
            {
                theta = Math.PI * ( double ) ( 2 * n - 1 - 2 * i ) 
                        / ( double ) ( 2 * n );
                x[i] = Math.Cos ( theta );
            }
            return x;
        }

        public static double[] fejer2 ( int n )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FEJER2 returns the Type 2 Fejer points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 March 2018
            //
            //  Author:
            //
            //    John Burkardt.
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, double FEJER2[N], the points.
            //
        {
            int i;
            
            double theta;
            double[] x;

            x = new double[n];

            for ( i = 0; i < n; i++ )
            {
                theta = Math.PI * ( double ) ( n - i ) 
                        / ( double ) ( n + 1 );
                x[i] = Math.Cos ( theta );
            }

            return x;
        }
    }
}