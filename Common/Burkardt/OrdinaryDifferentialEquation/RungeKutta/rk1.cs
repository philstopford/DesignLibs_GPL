using System;
using Burkardt.Types;

namespace Burkardt.ODENS;

public static partial class RungeKutta
{
    public static double rk1_tv_step ( double x, double t, double h, double q, 
            Func < double, double, double > fv, Func < double, double, double > gv, 
            ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RK1_TV_STEP takes one step of a stochastic Runge Kutta scheme.
        //
        //  Discussion:
        //
        //    The Runge-Kutta scheme is first-order, and suitable for time-varying
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
        //    Output, double RK1_TV_STEP the value at time T+H.
        //
    {
        typeMethods.r8NormalData data = new();

        const double a21 = 1.0;

        const double q1 = 1.0;

        double w1 = typeMethods.r8_normal_01(ref data, ref seed) * Math.Sqrt(q1 * q / h);
        double k1 = h * fv(t, x) + h * gv(t, x) * w1;

        double xstar = x + a21 * k1;

        return xstar;
    }

    public static double rk1_ti_step ( double x, double t, double h, double q, 
            Func< double, double > fi, Func< double, double > gi, ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RK1_TI_STEP takes one step of a stochastic Runge Kutta scheme.
        //
        //  Discussion:
        //
        //    The Runge-Kutta scheme is first-order, and suitable for time-invariant
        //    systems in which F and G do not depend explicitly on time.
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
        //    Output, double RK1_TI_STEP, the value at time T+H.
        //
    {
        typeMethods.r8NormalData data = new();

        const double a21 = 1.0;

        const double q1 = 1.0;

        double w1 = typeMethods.r8_normal_01 (ref data, ref seed ) * Math.Sqrt ( q1 * q / h );
        double k1 = h * fi ( x ) + h * gi ( x ) * w1;

        double xstar = x + a21 * k1;

        return xstar;
    }
}