using System;

namespace Burkardt.FDM;

public static class FD1D_Heat_Explicit
{
    public static double[] fd1d_heat_explicit(int x_num, double[] x, double t, double dt,
            double cfl, Func<int, double[], double, double[] > rhs,
            Func <int, double[], double, double[], double[] > bc, double[] h )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FD1D_HEAT_EXPLICIT: Finite difference solution of 1D heat equation.
        //
        //  Discussion:
        //
        //    This program takes one time step to solve the 1D heat equation 
        //    with an explicit method.
        //
        //    This program solves
        //
        //      dUdT - k * d2UdX2 = F(X,T)
        //
        //    over the interval [A,B] with boundary conditions
        //
        //      U(A,T) = UA(T),
        //      U(B,T) = UB(T),
        //
        //    over the time interval [T0,T1] with initial conditions
        //
        //      U(X,T0) = U0(X)
        //
        //    The code uses the finite difference method to approximate the
        //    second derivative in space, and an explicit forward Euler approximation
        //    to the first derivative in time.
        //
        //    The finite difference form can be written as
        //
        //      U(X,T+dt) - U(X,T)                  ( U(X-dx,T) - 2 U(X,T) + U(X+dx,T) )
        //      ------------------  = F(X,T) + k *  ------------------------------------
        //               dt                                   dx * dx
        //
        //    or, assuming we have solved for all values of U at time T, we have
        //
        //      U(X,T+dt) = U(X,T) 
        //        + cfl * ( U(X-dx,T) - 2 U(X,T) + U(X+dx,T) ) + dt * F(X,T) 
        //
        //    Here "cfl" is the Courant-Friedrichs-Loewy coefficient:
        //
        //      cfl = k * dt / dx / dx
        //
        //    In order for accurate results to be computed by this explicit method,
        //    the CFL coefficient must be less than 0.5.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 January 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X_NUM, the number of points to use in the 
        //    spatial dimension.
        //
        //    Input, double X(X_NUM), the coordinates of the nodes.
        //
        //    Input, double T, the current time.
        //
        //    Input, double DT, the size of the time step.
        //
        //    Input, double CFL, the Courant-Friedrichs-Loewy coefficient,
        //    computed by FD1D_HEAT_EXPLICIT_CFL.
        //
        //    Input, double H[X_NUM], the solution at the current time.
        //
        //    Input, double *RHS ( int x_num, double x[], double t ), the function 
        //    which evaluates the right hand side.
        //
        //    Input, void BC ( int x_num, double x[], double t, double h[] ), 
        //    the function which evaluates the boundary conditions.
        //
        //    Output, double FD1D_HEAT_EXPLICIT[X_NUM)], the solution at time T+DT.
        //
    {
        int j;

        double[] f = rhs(x_num, x, t);

        double[] h_new = new double[x_num];

        h_new[0] = 0.0;

        for (j = 1; j < x_num - 1; j++)
        {
            h_new[j] = h[j] + dt * f[j]
                            + cfl * (h[j - 1]
                                     - 2.0 * h[j]
                                     + h[j + 1]);
        }

        h_new[x_num - 1] = 0.0;

        h_new = bc(x_num, x, t + dt, h_new);

        return h_new;
    }

    public static double fd1d_heat_explicit_cfl(double k, int t_num, double t_min, double t_max,
            int x_num, double x_min, double x_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FD1D_HEAT_EXPLICIT_CFL: compute the Courant-Friedrichs-Loewy coefficient.
        //
        //  Discussion:
        //
        //    The equation to be solved has the form:
        //
        //      dUdT - k * d2UdX2 = F(X,T)
        //
        //    over the interval [X_MIN,X_MAX] with boundary conditions
        //
        //      U(X_MIN,T) = U_X_MIN(T),
        //      U(X_MIN,T) = U_X_MAX(T),
        //
        //    over the time interval [T_MIN,T_MAX] with initial conditions
        //
        //      U(X,T_MIN) = U_T_MIN(X)
        //
        //    The code uses the finite difference method to approximate the
        //    second derivative in space, and an explicit forward Euler approximation
        //    to the first derivative in time.
        //
        //    The finite difference form can be written as
        //
        //      U(X,T+dt) - U(X,T)                  ( U(X-dx,T) - 2 U(X,T) + U(X+dx,T) )
        //      ------------------  = F(X,T) + k *  ------------------------------------
        //               dt                                   dx * dx
        //
        //    or, assuming we have solved for all values of U at time T, we have
        //
        //      U(X,T+dt) = U(X,T) 
        //        + cfl * ( U(X-dx,T) - 2 U(X,T) + U(X+dx,T) ) + dt * F(X,T) 
        //
        //    Here "cfl" is the Courant-Friedrichs-Loewy coefficient:
        //
        //      cfl = k * dt / dx / dx
        //
        //    In order for accurate results to be computed by this explicit method,
        //    the CFL coefficient must be less than 0.5!
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 January 2012
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
        //    Input, double K, the heat conductivity coefficient.
        //
        //    Input, int T_NUM, the number of time values, including 
        //    the initial value.
        //
        //    Input, double T_MIN, T_MAX, the minimum and maximum times.
        //
        //    Input, int X_NUM, the number of nodes.
        //
        //    Input, double X_MIN, X_MAX, the minimum and maximum spatial 
        //    coordinates.
        //
        //    Output, double FD1D_HEAT_EXPLICIT_CFL, the Courant-Friedrichs-Loewy coefficient.
        //
    {
        double dx = (x_max - x_min) / (x_num - 1);
        double dt = (t_max - t_min) / (t_num - 1);
        //
        //  Check the CFL condition, print out its value, and quit if it is too large.
        //
        double cfl = k * dt / dx / dx;

        Console.WriteLine("");
        Console.WriteLine("  CFL stability criterion value = " + cfl + "");

        switch (cfl)
        {
            case >= 0.5:
                Console.WriteLine("");
                Console.WriteLine("FD1D_HEAT_EXPLICIT_CFL - Fatal error!");
                Console.WriteLine("  CFL condition failed.");
                Console.WriteLine("  0.5 <= K * dT / dX / dX = CFL.");
                return 1;
            default:
                return cfl;
        }
    }
}