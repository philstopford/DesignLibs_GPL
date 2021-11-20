using System;

namespace Burkardt.ODENS;

public static class MidpointFixed
{
    public static void midpoint_fixed ( Func < double, double[], int, double[] > dydt, 
            double[] tspan, double[] y0, int n, int m, int index, ref double[] t, ref double[] y )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    midpoint_fixed uses a fixed-point midpoint method to solve an ODE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 April 2021
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
        int i;

        double[] ym = new double[m];

        double dt = ( tspan[1] - tspan[0] ) / n;

        int it_max = 10;
        double theta = 0.5;

        t[0] = tspan[0];
        int j = 0;
        for ( i = 0; i < m; i++ )
        {
            y[i+j*m] = y0[i];
        }

        for ( j = 0; j < n; j++ )
        {
            double tm = t[j] + theta * dt;
            for ( i = 0; i < m; i++ )
            {
                ym[i] = y[i+j*m];
            }

            int k;
            for ( k = 0; k < it_max; k++ )
            {
                double[] f = dydt ( tm, ym, index );
                for ( i = 0; i < m; i++ )
                {
                    ym[i] = y[i+j*m] + theta * dt * f[i];
                }
            }
            t[j+1] = t[j] + dt;
            for ( i = 0; i < m; i++ )
            {
                y[i+(j+1)*m] = 1.0 / theta * ym[i]
                               + ( 1.0 - 1.0 / theta ) * y[i+j*m];
            }
        }
    }
}