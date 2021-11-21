using System;
using Burkardt.Types;

namespace Burkardt.FDM;

public static class FDM_1D_BVP
{
    public static double[] fd1d_bvp(int n, Func<double,double> a , Func<double,double> aprime,
        Func<double,double> c, Func<double,double> f, double[] x )
//****************************************************************************80
//
//  Purpose:
//
//    FD1D_BVP solves a two point boundary value problem.
//
//  Discussion:
//
//    The program uses the finite difference method to solve a BVP
//    (boundary value problem) in one dimension.
//
//    The problem is defined on the region X[0] <= x <= X[N-1].
//
//    The following differential equation is imposed in the region:
//
//      - d/dx a(x) du/dx + c(x) * u(x) = f(x)
//
//    where a(x), c(x), and f(x) are given functions.  We write out
//    the equation in full as
//
//      - a(x) * u''(x) - a'(x) * u'(x) + c(x) * u(x) = f(x)
//
//    At the boundaries, the following conditions are applied:
//
//      u(X[0]) = 0.0
//      u(X[N-1]) = 0.0
//
//    We replace the function U(X) by a vector of N values U associated
//    with the nodes.
//
//    The first and last values of U are determined by the boundary conditions.
//
//    At each interior node I, we write an equation to help us determine
//    U(I).  We do this by approximating the derivatives of U(X) by
//    finite differences.  Let us write XL, XM, and XR for X(I-1), X(I) and X(I+1).
//    Similarly we have UL, UM, and UR.  Other quantities to be evaluated at
//    X(I) = XM will also be labeled with an M:
//
//      - AM * ( UL - 2 UM + UR ) / DX^2 
//      - A'M * ( UL - UR ) / ( 2 * DX ) 
//      + CM * UM = FM
//
//    These N-2 linear equations for the unknown coefficients complete the
//    linear system and allow us to compute the finite difference approximation
//    to the solution of the BVP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Input, double A ( double x ), evaluates a(x);
//
//    Input, double APRIME ( double x ), evaluates a'(x);
//
//    Input, double C ( double x ), evaluates c(x);
//
//    Input, double F ( double x ), evaluates f(x);
//
//    Input, double X[N], the mesh points, which may be nonuniformly spaced.
//
//    Output, double FD1D_BVP[N], the value of the finite difference
//    approximation to the solution.
//
    {
        int i;
        //
//  Equation 1 is the left boundary condition, U(X[0]) = 0.0;
//
        double[] tri = new double[3 * n];
        double[] rhs = new double[n];

        tri[0 + 0 * 3] = 0.0;
        tri[1 + 0 * 3] = 1.0;
        tri[2 + 0 * 3] = 0.0;
        rhs[0] = 0.0;
//
//  Now gather the multipliers of U(I-1) to get the matrix entry A(I,I-1),
//  and so on.
//
        for (i = 1; i < n - 1; i++)
        {
            double xm = x[i];
            double am = a(xm);
            double apm = aprime(xm);
            double cm = c(xm);
            double fm = f(xm);

            tri[0 + i * 3] = -2.0 * am / (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1])
                             + apm / (x[i + 1] - x[i - 1]);

            tri[1 + i * 3] = +2.0 * am / (x[i] - x[i - 1]) / (x[i + 1] - x[i])
                             + cm;

            tri[2 + i * 3] = -2.0 * am / (x[i + 1] - x[i]) / (x[i + 1] - x[i - 1])
                             - apm / (x[i + 1] - x[i - 1]);

            rhs[i] = fm;
        }

//
//  Equation N is the right boundary condition, U(X[N-1]) = 0.0;
//
        tri[0 + (n - 1) * 3] = 0.0;
        tri[1 + (n - 1) * 3] = 1.0;
        tri[2 + (n - 1) * 3] = 0.0;
        rhs[n - 1] = 0.0;
//
//  Solve the linear system.
//
        double[] u = typeMethods.r83np_fs(n, ref tri, rhs);

        return u;
    }
}