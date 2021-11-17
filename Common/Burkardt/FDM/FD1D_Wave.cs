using System;

namespace Burkardt.FDM;

public static class FD1D_Wave
{
    public static double fd1d_wave_alpha(int x_num, double x1, double x2, int t_num, double t1,
            double t2, double c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FD1D_WAVE_ALPHA computes ALPHA for the 1D wave equation.
        //
        //  Discussion:
        //
        //    The explicit timestepping procedure uses the quantity ALPHA, which
        //    is determined by this function.
        //
        //    If the spatial region bounds are X1 <= X <= X2, containing X_NUM equally
        //    spaced nodes, including the endpoints, and the time domain similarly
        //    extends from T1 <= T <= T2 containing T_NUM equally spaced time values,
        //    then
        //
        //      ALPHA = C * DT / DX
        //            = C * ( ( T2 - T1 ) / ( T_NUM - 1 ) )
        //                / ( ( X2 - X1 ) / ( X_NUM - 1 ) ).
        //
        //    For a stable computation, it must be the case that ALPHA < 1.
        //
        //    If ALPHA is greater than 1, then the middle coefficient 1-C^2 DT^2 / DX^2 
        //    is negative, and the sum of the magnitudes of the three coefficients 
        //    becomes unbounded.  In such a case, the user must reduce the time step 
        //    size appropriately.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 January 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    George Lindfield, John Penny,
        //    Numerical Methods Using MATLAB,
        //    Second Edition,
        //    Prentice Hall, 1999,
        //    ISBN: 0-13-012641-1,
        //    LC: QA297.P45.
        //
        //  Parameters:
        //
        //    Input, int X_NUM, the number of nodes in the X direction.
        //
        //    Input, double X1, X2, the first and last X coordinates.
        //
        //    Input, int T_NUM, the number of time steps, including the 
        //    initial condition.
        //
        //    Input, double T1, T2, the first and last T coordinates.
        //
        //    Input, double C, a parameter which gives the speed of waves.
        //
        //    Output, double FD1D_WAVE_ALPHA, the stability coefficient.
        //
    {
        double alpha;
        double t_delta;
        double x_delta;

        t_delta = (t2 - t1) / (t_num - 1);
        x_delta = (x2 - x1) / (x_num - 1);
        alpha = c * t_delta / x_delta;

        Console.WriteLine("");
        Console.WriteLine("  Stability condition ALPHA = C * DT / DX = " + alpha + "");

        switch (Math.Abs(alpha))
        {
            case > 1.0:
                Console.WriteLine("");
                Console.WriteLine("FD1D_WAVE_ALPHA - Warning!");
                Console.WriteLine("  The stability condition |ALPHA| <= 1 fails.");
                Console.WriteLine("  Computed results are liable to be inaccurate.");
                break;
        }

        return alpha;
    }

    public static double[] fd1d_wave_start(int x_num, double[] x_vec, double t, double t_delta,
            double alpha, Func<double, double> u_x1, Func<double, double> u_x2,
            Func <int, double[], double[]> ut_t1, double[] u1 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FD1D_WAVE_START takes the first step for the wave equation.
        //
        //  Discussion:
        //
        //    This program solves the 1D wave equation of the form:
        //
        //      Utt = c^2 Uxx
        //
        //    over the spatial interval [X1,X2] and time interval [T1,T2],
        //    with initial conditions:
        //
        //      U(T1,X)  = U_T1(X),
        //      Ut(T1,X) = UT_T1(X),
        //
        //    and boundary conditions of Dirichlet type:
        //
        //      U(T,X1) = U_X1(T),
        //      U(T,X2) = U_X2(T).
        //
        //    The value C represents the propagation speed of waves.
        //
        //    The program uses the finite difference method, and marches
        //    forward in time, solving for all the values of U at the next
        //    time step by using the values known at the previous two time steps.
        //
        //    Central differences may be used to approximate both the time
        //    and space derivatives in the original differential equation.
        //
        //    Thus, assuming we have available the approximated values of U
        //    at the current and previous times, we may write a discretized
        //    version of the wave equation as follows:
        //
        //      Uxx(T,X) = ( U(T,   X+dX) - 2 U(T,X) + U(T,   X-dX) ) / dX^2
        //      Utt(T,X) = ( U(T+dt,X   ) - 2 U(T,X) + U(T-dt,X   ) ) / dT^2
        //
        //    If we multiply the first term by C^2 and solve for the single
        //    unknown value U(T+dt,X), we have:
        //
        //      U(T+dT,X) =        (     C^2 * dT^2 / dX^2 ) * U(T,   X+dX)
        //                  +  2 * ( 1 - C^2 * dT^2 / dX^2 ) * U(T,   X   )
        //                  +      (     C^2 * dT^2 / dX^2 ) * U(T,   X-dX)
        //                  -                                  U(T-dT,X   )
        //
        //    (Equation to advance from time T to time T+dT, except for FIRST step!)
        //
        //    However, on the very first step, we only have the values of U
        //    for the initial time, but not for the previous time step.
        //    In that case, we use the initial condition information for dUdT
        //    which can be approximated by a central difference that involves
        //    U(T+dT,X) and U(T-dT,X):
        //
        //      dU/dT(T,X) = ( U(T+dT,X) - U(T-dT,X) ) / ( 2 * dT )
        //
        //    and so we can estimate U(T-dT,X) as
        //
        //      U(T-dT,X) = U(T+dT,X) - 2 * dT * dU/dT(T,X)
        //
        //    If we replace the "missing" value of U(T-dT,X) by the known values
        //    on the right hand side, we now have U(T+dT,X) on both sides of the
        //    equation, so we have to rearrange to get the formula we use
        //    for just the first time step:
        //
        //      U(T+dT,X) =   1/2 * (     C^2 * dT^2 / dX^2 ) * U(T,   X+dX)
        //                  +       ( 1 - C^2 * dT^2 / dX^2 ) * U(T,   X   )
        //                  + 1/2 * (     C^2 * dT^2 / dX^2 ) * U(T,   X-dX)
        //                  +  dT *                         dU/dT(T,   X   )
        //
        //    (Equation to advance from time T to time T+dT for FIRST step.)
        //
        //    It should be clear now that the quantity ALPHA = C * DT / DX will affect
        //    the stability of the calculation.  If it is greater than 1, then
        //    the middle coefficient 1-C^2 DT^2 / DX^2 is negative, and the
        //    sum of the magnitudes of the three coefficients becomes unbounded.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 January 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    George Lindfield, John Penny,
        //    Numerical Methods Using MATLAB,
        //    Second Edition,
        //    Prentice Hall, 1999,
        //    ISBN: 0-13-012641-1,
        //    LC: QA297.P45.
        //
        //  Parameters:
        //
        //    Input, int X_NUM, the number of nodes in the X direction.
        //
        //    Input, double X_VEC[X_NUM], the coordinates of the nodes.
        //
        //    Input, double T, the time after the first step has been taken.
        //    In other words, T = T1 + T_DELTA.
        //
        //    Input, double T_DELTA, the time step.
        //
        //    Input, double ALPHA, the stability coefficient, computed 
        //    by FD1D_WAVE_ALPHA.
        //
        //    Input, double U_X1 ( double T ), U_X2 ( double T ), functions for the left and 
        //    right boundary conditions.
        //
        //    Input, double *UT_T1 ( int X_NUM, double X_VEC[] ), the function that 
        //    evaluates dUdT at the initial time.
        //
        //    Input, real U1[X_NUM], the initial condition.
        //
        //    Output, real FD1D_WAVE_START[X_NUM], the solution at the first time step.
        //
    {
        int j;
        double[] u2;
        double[] ut;

        ut = ut_t1(x_num, x_vec);

        u2 = new double[x_num];

        u2[0] = u_x1(t);

        for (j = 1; j < x_num - 1; j++)
        {
            u2[j] = alpha * alpha * u1[j + 1] / 2.0
                    + (1.0 - alpha * alpha) * u1[j]
                    + alpha * alpha * u1[j - 1] / 2.0
                    + t_delta * ut[j];
        }

        u2[x_num - 1] = u_x2(t);

        return u2;
    }

    public static double[] fd1d_wave_step(int x_num, double t, double alpha,
            Func < double, double > u_x1, Func <double, double > u_x2, double[] u1, double[] u2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FD1D_WAVE_STEP computes a step of the 1D wave equation.
        //
        //  Discussion:
        //
        //    This program solves the 1D wave equation of the form:
        //
        //      Utt = c^2 Uxx
        //
        //    over the spatial interval [X1,X2] and time interval [T1,T2],
        //    with initial conditions:
        //
        //      U(T1,X)  = U_T1(X),
        //      Ut(T1,X) = UT_T1(X),
        //
        //    and boundary conditions of Dirichlet type:
        //
        //      U(T,X1) = U_X1(T),
        //      U(T,X2) = U_X2(T).
        //
        //    The value C represents the propagation speed of waves.
        //
        //    The program uses the finite difference method, and marches
        //    forward in time, solving for all the values of U at the next
        //    time step by using the values known at the previous two time steps.
        //
        //    Central differences may be used to approximate both the time
        //    and space derivatives in the original differential equation.
        //
        //    Thus, assuming we have available the approximated values of U
        //    at the current and previous times, we may write a discretized
        //    version of the wave equation as follows:
        //
        //      Uxx(T,X) = ( U(T,   X+dX) - 2 U(T,X) + U(T,   X-dX) ) / dX^2
        //      Utt(T,X) = ( U(T+dt,X   ) - 2 U(T,X) + U(T-dt,X   ) ) / dT^2
        //
        //    If we multiply the first term by C^2 and solve for the single
        //    unknown value U(T+dt,X), we have:
        //
        //      U(T+dT,X) =        (     C^2 * dT^2 / dX^2 ) * U(T,   X+dX)
        //                  +  2 * ( 1 - C^2 * dT^2 / dX^2 ) * U(T,   X   )
        //                  +      (     C^2 * dT^2 / dX^2 ) * U(T,   X-dX)
        //                  -                                  U(T-dT,X   )
        //
        //    (Equation to advance from time T to time T+dT, except for FIRST step!)
        //
        //    However, on the very first step, we only have the values of U
        //    for the initial time, but not for the previous time step.
        //    In that case, we use the initial condition information for dUdT
        //    which can be approximated by a central difference that involves
        //    U(T+dT,X) and U(T-dT,X):
        //
        //      dU/dT(T,X) = ( U(T+dT,X) - U(T-dT,X) ) / ( 2 * dT )
        //
        //    and so we can estimate U(T-dT,X) as
        //
        //      U(T-dT,X) = U(T+dT,X) - 2 * dT * dU/dT(T,X)
        //
        //    If we replace the "missing" value of U(T-dT,X) by the known values
        //    on the right hand side, we now have U(T+dT,X) on both sides of the
        //    equation, so we have to rearrange to get the formula we use
        //    for just the first time step:
        //
        //      U(T+dT,X) =   1/2 * (     C^2 * dT^2 / dX^2 ) * U(T,   X+dX)
        //                  +       ( 1 - C^2 * dT^2 / dX^2 ) * U(T,   X   )
        //                  + 1/2 * (     C^2 * dT^2 / dX^2 ) * U(T,   X-dX)
        //                  +  dT *                         dU/dT(T,   X   )
        //
        //    (Equation to advance from time T to time T+dT for FIRST step.)
        //
        //    It should be clear now that the quantity ALPHA = C * DT / DX will affect
        //    the stability of the calculation.  If it is greater than 1, then
        //    the middle coefficient 1-C^2 DT^2 / DX^2 is negative, and the
        //    sum of the magnitudes of the three coefficients becomes unbounded.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 January 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    George Lindfield, John Penny,
        //    Numerical Methods Using MATLAB,
        //    Second Edition,
        //    Prentice Hall, 1999,
        //    ISBN: 0-13-012641-1,
        //    LC: QA297.P45.
        //
        //  Parameters:
        //
        //    Input, int X_NUM, the number of nodes in the X direction.
        //
        //    Input, double T, the new time, that is, the current 
        //    time + T_DELTA.
        //
        //    Input, double ALPHA, the stability coefficient, computed 
        //    by FD1D_WAVE_ALPHA.
        //
        //    Input, double U_X1 ( double T ), U_X2 ( double T ), functions for the left and 
        //    right boundary conditions.
        //
        //    Input, double U1[X_NUM], the solution at the old time.
        //
        //    Input, double U2[X_NUM], the solution at the current time.
        //
        //    Output, double FD1D_WAVE_STEP[X_NUM], the solution at the new time.
        //
    {
        int j;
        double[] u3;

        u3 = new double[x_num];

        u3[0] = u_x1(t);

        for (j = 1; j < x_num - 1; j++)
        {
            u3[j] = alpha * alpha * u2[j + 1]
                    + 2.0 * (1.0 - alpha * alpha) * u2[j]
                    + alpha * alpha * u2[j - 1]
                    - u1[j];
        }

        u3[x_num - 1] = u_x2(t);

        return u3;
    }
}