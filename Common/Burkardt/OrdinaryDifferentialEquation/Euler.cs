using System;

namespace Burkardt
{
    public static class Euler
    {
        public static void euler ( Func < double, double[], int, double[] > dydt, double[] tspan, 
        double[] y0, int n, int m, double[] t, double[] y )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    euler approximates the solution to an ODE using Euler's method.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 April 2020
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Input:
        //
        //    dydt: a function that evaluates the right hand side of the ODE.
        //
        //    double tspan[2]: contains the initial and final times.
        //
        //    double y0[m]: a column vector containing the initial condition.
        //
        //    int n: the number of steps to take.
        //
        //    int m: the number of variables.
        //
        //  Output:
        //
        //    double t[n+1], y[m*(n+1)]: the times and solution values.
        //
        {
            double dt;
            double[] dy;
            int i;
            int j;
            double tfirst;
            double tlast;

            tfirst = tspan[0];
            tlast = tspan[1];
            dt = ( tlast - tfirst ) / ( double ) ( n );
            j = 0;
            t[j] = tspan[0];
            for ( i = 0; i < m; i++ )
            {
                y[i+j*m] = y0[i];
            }

            for ( j = 0; j < n; j++ )
            {
                dy = dydt ( t[j], y, +j*m );
                t[j+1] = t[j] + dt;
                for ( i = 0; i < m; i++ )
                {
                    y[i+(j+1)*m] = y[i+j*m] + dt * dy[i];
                }
            }
        }
    }
}