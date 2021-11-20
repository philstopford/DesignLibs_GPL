using System;

namespace Burkardt.ODENS;

public static class MidpointExplicit
{
    public static void midpoint_explicit ( Func < double, double[], int, double[] > dydt, 
            double[] tspan, double[] y0, int n, int m, ref double[] t, ref double[] y )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    midpoint_explicit uses the explicit midpoint method to solve an ODE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 April 2021
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

        t[0] = tspan[0];
        int j = 0;
        for ( i = 0; i < m; i++ )
        {
            y[i+j*m] = y0[i];
        }

        for ( j = 0; j < n; j++ )
        {
            double[] f = dydt ( t[j], y, +j*m );
            double tm = t[j] + 0.5 * dt;
            for ( i = 0; i < m; i++ )
            {
                ym[i] = y[i+j*m] + 0.5 * dt * f[i];
            }

            f = dydt ( tm, ym, 0 );
            t[j+1] = t[j] + dt;
            for ( i = 0; i < m; i++ )
            {
                y[i+(j+1)*m] = y[i+j*m] + dt * f[i];
            }
        }
    }
}