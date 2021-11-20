using System;
using Burkardt.Types;

namespace Burkardt.ODENS;

public static partial class RungeKutta
{
    public static double rk4_tv_step(double x, double t, double h, double q,
            Func<double, double, double> fv, Func<double, double, double> gv,
            ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RK4_TV_STEP takes one step of a stochastic Runge Kutta scheme.
        //
        //  Discussion:
        //
        //    The Runge-Kutta scheme is fourth-order, and suitable for time-varying
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
        //    Output, double RK4_TV_STEP, the value at time T+H.
        //
    {
        typeMethods.r8NormalData data = new();

        const double a21 = 0.66667754298442;
        const double a31 = 0.63493935027993;
        const double a32 = 0.00342761715422;
        const double a41 = -2.32428921184321;
        const double a42 = 2.69723745129487;
        const double a43 = 0.29093673271592;
        const double a51 = 0.25001351164789;
        const double a52 = 0.67428574806272;
        const double a53 = -0.00831795169360;
        const double a54 = 0.08401868181222;

        const double q1 = 3.99956364361748;
        const double q2 = 1.64524970733585;
        const double q3 = 1.59330355118722;
        const double q4 = 0.26330006501868;

        double w1 = typeMethods.r8_normal_01(ref data, ref seed) * Math.Sqrt(q1 * q / h);
        double k1 = h * fv(t, x) + h * gv(t, x) * w1;

        double t2 = t + a21 * h;
        double x2 = x + a21 * k1;
        double w2 = typeMethods.r8_normal_01(ref data, ref seed) * Math.Sqrt(q2 * q / h);
        double k2 = h * fv(t2, x2) + h * gv(t2, x2) * w2;

        double t3 = t + a31 * h + a32 * h;
        double x3 = x + a31 * k1 + a32 * k2;
        double w3 = typeMethods.r8_normal_01(ref data, ref seed) * Math.Sqrt(q3 * q / h);
        double k3 = h * fv(t3, x3) + h * gv(t3, x3) * w3;

        double t4 = t + a41 * h + a42 * h + a43 * h;
        double x4 = x + a41 * k1 + a42 * k2 + a43 * k3;
        double w4 = typeMethods.r8_normal_01(ref data, ref seed) * Math.Sqrt(q4 * q / h);
        double k4 = h * fv(t4, x4) + h * gv(t4, x4) * w4;

        double xstar = x + a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4;

        return xstar;
    }

    public static double rk4_ti_step ( double x, double t, double h, double q, 
            Func< double, double > fi, Func< double, double > gi, ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RK4_TI_STEP takes one step of a stochastic Runge Kutta scheme.
        //
        //  Discussion:
        //
        //    The Runge-Kutta scheme is fourth-order, and suitable for time-invariant
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
        //    Output, double RK4_TI_STEP, the value at time T+H.
        //
    {
        //double a43;
        //double t1;
        typeMethods.r8NormalData data = new();

        const double a21 = 2.71644396264860;
        const double a31 = -6.95653259006152;
        const double a32 = 0.78313689457981;
        const double a41 = 0.0;
        const double a42 = 0.48257353309214;
        //a43 =   0.26171080165848;
        const double a51 = 0.47012396888046;
        const double a52 = 0.36597075368373;
        const double a53 = 0.08906615686702;
        const double a54 = 0.07483912056879;

        const double q1 = 2.12709852335625;
        const double q2 = 2.73245878238737;
        const double q3 = 11.22760917474960;
        const double q4 = 13.36199560336697;

        double w1 = typeMethods.r8_normal_01(ref data, ref seed) * Math.Sqrt(q1 * q / h);
        double k1 = h * fi(x) + h * gi(x) * w1;

        double x2 = x + a21 * k1;
        double w2 = typeMethods.r8_normal_01(ref data, ref seed) * Math.Sqrt(q2 * q / h);
        double k2 = h * fi(x2) + h * gi(x2) * w2;

        double x3 = x + a31 * k1 + a32 * k2;
        double w3 = typeMethods.r8_normal_01(ref data, ref seed) * Math.Sqrt(q3 * q / h);
        double k3 = h * fi(x3) + h * gi(x3) * w3;

        double x4 = x + a41 * k1 + a42 * k2;
        double w4 = typeMethods.r8_normal_01(ref data, ref seed) * Math.Sqrt(q4 * q / h);
        double k4 = h * fi(x4) + h * gi(x4) * w4;

        double xstar = x + a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4;

        return xstar;
    }

    public static double[] rk4vec(double t0, int m, double[] u0, double dt,
            Func<double, int, double[], int, double[]> f, int index = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RK4VEC takes one Runge-Kutta step for a vector ODE.
        //
        //  Discussion:
        //
        //    It is assumed that an initial value problem, of the form
        //
        //      du/dt = f ( t, u )
        //      u(t0) = u0
        //
        //    is being solved.
        //
        //    If the user can supply current values of t, u, a stepsize dt, and a
        //    function to evaluate the derivative, this function can compute the
        //    fourth-order Runge Kutta estimate to the solution at time t+dt.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T0, the current time.
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, double U0[M], the solution estimate at the current time.
        //
        //    Input, double DT, the time step.
        //
        //    Input, double[] F ( double T, int M, double U[] ), a function which evaluates
        //    the derivative, or right hand side of the problem.
        //
        //    Output, double RK4VEC[M], the fourth-order Runge-Kutta solution estimate
        //    at time T0+DT.
        //
    {
        int i;
        //
        //  Get four sample values of the derivative.
        //
        double[] f0 = f(t0, m, u0, index);

        double t1 = t0 + dt / 2.0;
        double[] u1 = new double[m];
        for (i = 0; i < m; i++)
        {
            u1[i] = u0[index + i] + dt * f0[i] / 2.0;
        }

        double[] f1 = f(t1, m, u1, 0);

        double t2 = t0 + dt / 2.0;
        double[] u2 = new double[m];
        for (i = 0; i < m; i++)
        {
            u2[i] = u0[index + i] + dt * f1[i] / 2.0;
        }

        double[] f2 = f(t2, m, u2, 0);

        double t3 = t0 + dt;
        double[] u3 = new double[m];
        for (i = 0; i < m; i++)
        {
            u3[i] = u0[index + i] + dt * f2[i];
        }

        double[] f3 = f(t3, m, u3, 0);
        //
        //  Combine them to estimate the solution.
        //
        double[] u = new double[m];
        for (i = 0; i < m; i++)
        {
            u[i] = u0[index + i] + dt * (f0[i] + 2.0 * f1[i] + 2.0 * f2[i] + f3[i]) / 6.0;
        }

        return u;
    }
}