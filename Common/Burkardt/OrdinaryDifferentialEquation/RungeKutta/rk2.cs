using System;

namespace Burkardt.ODENS
{
    public static partial class RungeKutta
    {
        public static double rk2_leg ( double t1, double t2, double x, int n )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    RK2_LEG advances the value of X(T) using a Runge-Kutta method.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    22 October 2009
            //
            //  Author:
            //
            //    Original C++ version by Nick Hale.
            //    This C++ version by John Burkardt.
            //
            //  Parameters:
            //
            //    Input, double T1, T2, the range of the integration interval.
            //
            //    Input, double X, the value of X at T1.
            //
            //    Input, int N, the number of steps to take.
            //
            //    Output, double RK2_LEG, the value of X at T2.
            //
        {
            double f;
            double h;
            int j;
            double k1;
            double k2;
            int m = 10;
            double snn1;
            double t;

            h = ( t2 - t1 ) / ( double ) m;
            snn1 = Math.Sqrt ( ( double ) ( n * ( n + 1 ) ) );
            t = t1;

            for ( j = 0; j < m; j++ )
            {
                f = ( 1.0 - x ) * ( 1.0 + x );
                k1 = - h * f / ( snn1 * Math.Sqrt ( f ) - 0.5 * x * Math.Sin ( 2.0 * t ) );
                x = x + k1;

                t = t + h;

                f = ( 1.0 - x ) * ( 1.0 + x );
                k2 = - h * f / ( snn1 * Math.Sqrt ( f ) - 0.5 * x * Math.Sin ( 2.0 * t ) );
                x = x + 0.5 * ( k2 - k1 );
            }
            return x;
        }
    }
}