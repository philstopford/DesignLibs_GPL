using System;
using Burkardt.Types;

namespace Burkardt.ODENS;

public static partial class RungeKutta
{
    public static double rk2_tv_step(double x, double t, double h, double q,
            Func<double, double, double> fv, Func<double, double, double> gv,
            ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RK2_TV_STEP takes one step of a stochastic Runge Kutta scheme.
        //
        //  Discussion:
        //
        //    The Runge-Kutta scheme is second-order, and suitable for time-varying
        //    systems.
        //
        //    d/dx X(t,xsi) = F ( X(t,xsi), t ) + G ( X(t,xsi), t ) * w(t,xsi)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Jeremy Kasdin,
        //    Runge-Kutta algorithm for the numerical integration of
        //    stochastic differential equations,
        //    Journal of Guidance, Control, and Dynamics,
        //    Volume 18, Number 1, January-February 1995, pages 114-120.
        //
        //    Jeremy Kasdin,
        //    Discrete Simulation of Colored Noise and Stochastic Processes
        //    and 1/f^a Power Law Noise Generation,
        //    Proceedings of the IEEE,
        //    Volume 83, Number 5, 1995, pages 802-827.
        //
        //  Parameters:
        //
        //    Input, double X, the value at the current time.
        //
        //    Input, double T, the current time.
        //
        //    Input, double H, the time step.
        //
        //    Input, double Q, the spectral density of the input white noise.
        //
        //    Input, double FV ( double T, double X ), the name of the deterministic
        //    right hand side function.
        //
        //    Input, double GV ( double T, double X ), the name of the stochastic
        //    right hand side function.
        //
        //    Input/output, int *SEED, a seed for the random 
        //    number generator.
        //
        //    Output, double RK2_TV_STEP, the value at time T+H.
        //
    {
        typeMethods.r8NormalData data = new();

        const double a21 = 1.0;
        const double a31 = 0.5;
        const double a32 = 0.5;

        const double q1 = 2.0;
        const double q2 = 2.0;

        double w1 = typeMethods.r8_normal_01(ref data, ref seed) * Math.Sqrt(q1 * q / h);
        double k1 = h * fv(t, x) + h * gv(t, x) * w1;

        double t2 = t + a21 * h;
        double x2 = x + a21 * k1;
        double w2 = typeMethods.r8_normal_01(ref data, ref seed) * Math.Sqrt(q2 * q / h);
        double k2 = h * fv(t2, x2) + h * gv(t2, x2) * w2;

        double xstar = x + a31 * k1 + a32 * k2;

        return xstar;
    }

    public static double rk2_ti_step ( double x, double t, double h, double q, 
            Func< double, double > fi, Func< double, double > gi, ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RK2_TI_STEP takes one step of a stochastic Runge Kutta scheme.
        //
        //  Discussion:
        //
        //    The Runge-Kutta scheme is second-order, and suitable for time-invariant
        //    systems.
        //
        //    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Jeremy Kasdin,
        //    Runge-Kutta algorithm for the numerical integration of
        //    stochastic differential equations,
        //    Journal of Guidance, Control, and Dynamics,
        //    Volume 18, Number 1, January-February 1995, pages 114-120.
        //
        //    Jeremy Kasdin,
        //    Discrete Simulation of Colored Noise and Stochastic Processes
        //    and 1/f^a Power Law Noise Generation,
        //    Proceedings of the IEEE,
        //    Volume 83, Number 5, 1995, pages 802-827.
        //
        //  Parameters:
        //
        //    Input, double X, the value at the current time.
        //
        //    Input, double T, the current time.
        //
        //    Input, double H, the time step.
        //
        //    Input, double Q, the spectral density of the input white noise.
        //
        //    Input, double FI ( double X ), the name of the deterministic
        //    right hand side function.
        //
        //    Input, double GI ( double X ), the name of the stochastic
        //    right hand side function.
        //
        //    Input/output, int *SEED, a seed for the random 
        //    number generator.
        //
        //    Output, double RK2_TI_STEP, the value at time T+H.
        //
    {
        typeMethods.r8NormalData data = new();

        const double a21 = 1.0;
        const double a31 = 0.5;
        const double a32 = 0.5;

        const double q1 = 2.0;
        const double q2 = 2.0;

        double w1 = typeMethods.r8_normal_01(ref data, ref seed) * Math.Sqrt(q1 * q / h);
        double k1 = h * fi(x) + h * gi(x) * w1;

        double x2 = x + a21 * k1;
        double w2 = typeMethods.r8_normal_01(ref data, ref seed) * Math.Sqrt(q2 * q / h);
        double k2 = h * fi(x2) + h * gi(x2) * w2;

        double xstar = x + a31 * k1 + a32 * k2;

        return xstar;
    }

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
        //    C++ version by John Burkardt.
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
        int j;
        const int m = 10;

        double h = ( t2 - t1 ) / m;
        double snn1 = Math.Sqrt ( n * ( n + 1 ) );
        double t = t1;

        for ( j = 0; j < m; j++ )
        {
            double f = ( 1.0 - x ) * ( 1.0 + x );
            double k1 = - h * f / ( snn1 * Math.Sqrt ( f ) - 0.5 * x * Math.Sin ( 2.0 * t ) );
            x += k1;

            t += h;

            f = ( 1.0 - x ) * ( 1.0 + x );
            double k2 = - h * f / ( snn1 * Math.Sqrt ( f ) - 0.5 * x * Math.Sin ( 2.0 * t ) );
            x += 0.5 * ( k2 - k1 );
        }
        return x;
    }
}