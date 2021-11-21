using System;
using Burkardt.Types;

namespace Burkardt.FDM;

public static class FD1D_Heat_Steady
{
    public static double[] fd1d_heat_steady(int n, double a, double b, double ua, double ub,
            Func < double, double > k, Func<double, double> f, double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FD1D_HEAT_STEADY solves the steady 1D heat equation.
        //
        //  Discussion:
        //
        //    This program seeks a solution of the steady heat equation:
        //
        //      - d/dx ( K(X) dUdx ) = F(X)
        //
        //    over the interval [A,B] with boundary conditions
        //
        //      U(A) = UA,
        //      U(B) = UB.
        //
        //    The code uses the finite difference method to approximate the
        //    second derivative in space.  This results in a sparse linear system.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of grid points.
        //
        //    Input, double A, B, the interval endpoints.
        //
        //    Input, double UA, UB, the values prescribed for U at the endpoints.
        //
        //    Input, double K ( double X ), evaluates the thermal conductance at the N
        //    points X.  Set K(X) = 1 if you don't care about this coefficient.
        //
        //    Input, double F ( double X ), evaluates the heat source term at the N
        //    points X.  Set F(X) = 0 if you don't want any heat sources.
        //
        //    Input, double X[N], the grid points.
        //
        //    Output, double FD1D_HEAT_STEADY[N], the approximation to the solution
        //    at the grid points.
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("FD1D_HEAT_STEADY");
            
        Console.WriteLine("");
        Console.WriteLine("  Finite difference solution of");
        Console.WriteLine("  the steady 1D heat equation");
        Console.WriteLine("");
        Console.WriteLine("    - d/dx ( k(x) dUdx ) = F(x)");
        Console.WriteLine("");
        Console.WriteLine("  for space interval A <= X <= B with boundary conditions");
        Console.WriteLine("");
        Console.WriteLine("    U(A) = UA");
        Console.WriteLine("    U(B) = UB");
        Console.WriteLine("");
        Console.WriteLine("  A second order difference approximation is used.");
        //
        //  Set the spacing.
        //
        double dx = (b - a) / (n - 1);
        //
        //  Set up the tridiagonal matrix.
        //
        double[] tri = new double[3 * n];
        double[] rhs = new double[n];

        tri[0 + 0 * 3] = 0.0;
        tri[1 + 0 * 3] = 1.0;
        tri[2 + 0 * 3] = 0.0;
        rhs[0] = ua;

        for (i = 1; i < n - 1; i++)
        {
            double xm = (x[i - 1] + x[i]) / 2.0;
            double xp = (x[i] + x[i + 1]) / 2.0;

            tri[0 + i * 3] = -k(xm) / dx / dx;
            tri[1 + i * 3] = (k(xm) + k(xp)) / dx / dx;
            tri[2 + i * 3] = -k(xp) / dx / dx;

            rhs[i] = f(x[i]);
        }

        tri[0 + (n - 1) * 3] = 0.0;
        tri[1 + (n - 1) * 3] = 1.0;
        tri[2 + (n - 1) * 3] = 0.0;
        rhs[n - 1] = ub;
        //
        //  Solve the linear system.
        //
        double[] u = typeMethods.r83np_fs(n, ref tri, rhs);

        return u;
    }
}