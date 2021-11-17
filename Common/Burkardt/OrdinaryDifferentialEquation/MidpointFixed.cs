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
        double dt;
        double[] f;
        int i;
        int it_max;
        int j;
        int k;
        double theta;
        double tm;
        double[] ym;

        ym = new double[m];

        dt = ( tspan[1] - tspan[0] ) / n;

        it_max = 10;
        theta = 0.5;

        t[0] = tspan[0];
        j = 0;
        for ( i = 0; i < m; i++ )
        {
            y[i+j*m] = y0[i];
        }

        for ( j = 0; j < n; j++ )
        {
            tm = t[j] + theta * dt;
            for ( i = 0; i < m; i++ )
            {
                ym[i] = y[i+j*m];
            }
            for ( k = 0; k < it_max; k++ )
            {
                f = dydt ( tm, ym, index );
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