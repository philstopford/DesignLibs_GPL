using System;
using Burkardt.Types;

namespace Burkardt.ODENS;

public static partial class RungeKutta
{
    public static double rk3_ti_step(double x, double t, double h, double q,
            Func< double, double > fi, Func< double, double > gi, ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RK3_TI_STEP takes one step of a stochastic Runge Kutta scheme.
        //
        //  Discussion:
        //
        //    The Runge-Kutta scheme is third-order, and suitable for time-invariant
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
        //    Output, double RK3_TI_STEP, the value at time T+H.
        //
    {
        //double t1;
        typeMethods.r8NormalData data = new();

        const double a21 = 1.52880952525675;
        const double a31 = 0.0;
        const double a32 = 0.51578733443615;
        const double a41 = 0.53289582961739;
        const double a42 = 0.25574324768195;
        const double a43 = 0.21136092270067;

        const double q1 = 1.87653936176981;
        const double q2 = 3.91017166264989;
        const double q3 = 4.73124353935667;

        //t1 = t;
        double w1 = typeMethods.r8_normal_01(ref data, ref seed) * Math.Sqrt(q1 * q / h);
        double k1 = h * fi(x) + h * gi(x) * w1;

        //t2 = t1 + a21 * h;
        double x2 = x + a21 * k1;
        double w2 = typeMethods.r8_normal_01(ref data, ref seed) * Math.Sqrt(q2 * q / h);
        double k2 = h * fi(x2) + h * gi(x2) * w2;

        //t3 = t1 + a31 * h  + a32 * h;
        double x3 = x + a31 * k1 + a32 * k2;
        double w3 = typeMethods.r8_normal_01(ref data, ref seed) * Math.Sqrt(q3 * q / h);
        double k3 = h * fi(x3) + h * gi(x3) * w3;

        double xstar = x + a41 * k1 + a42 * k2 + a43 * k3;

        return xstar;
    }
}