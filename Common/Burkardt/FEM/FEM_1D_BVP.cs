using System;
using Burkardt.Types;

namespace Burkardt.FEM;

public static class FEM_1D_BVP
{
    public static double[] fem1d_bvp_linear(int n, Func<double, double> a, Func<double, double> c,
        Func<double, double> f, double[] x)
//****************************************************************************80
//
//  Purpose:
//
//    FEM1D_BVP_LINEAR solves a two point boundary value problem.
//
//  Location:
//
//    http://people.sc.fsu.edu/~jburkardt/cpp_src/fem1d_bvp_linear/fem1d_bvp_linear.cpp
//
//  Discussion:
//
//    The program uses the finite element method, with piecewise linear basis
//    functions to solve a boundary value problem in one dimension.
//
//    The problem is defined on the region 0 <= x <= 1.
//
//    The following differential equation is imposed between 0 and 1:
//
//      - d/dx a(x) du/dx + c(x) * u(x) = f(x)
//
//    where a(x), c(x), and f(x) are given functions.
//
//    At the boundaries, the following conditions are applied:
//
//      u(0.0) = 0.0
//      u(1.0) = 0.0
//
//    A set of N equally spaced nodes is defined on this
//    interval, with 0 = X(1) < X(2) < ... < X(N) = 1.0.
//
//    At each node I, we associate a piecewise linear basis function V(I,X),
//    which is 0 at all nodes except node I.  This implies that V(I,X) is
//    everywhere 0 except that
//
//    for X(I-1) <= X <= X(I):
//
//      V(I,X) = ( X - X(I-1) ) / ( X(I) - X(I-1) ) 
//
//    for X(I) <= X <= X(I+1):
//
//      V(I,X) = ( X(I+1) - X ) / ( X(I+1) - X(I) )
//
//    We now assume that the solution U(X) can be written as a linear
//    sum of these basis functions:
//
//      U(X) = sum ( 1 <= J <= N ) U(J) * V(J,X)
//
//    where U(X) on the left is the function of X, but on the right,
//    is meant to indicate the coefficients of the basis functions.
//
//    To determine the coefficient U(J), we multiply the original
//    differential equation by the basis function V(J,X), and use
//    integration by parts, to arrive at the I-th finite element equation:
//
//        Integral A(X) * U'(X) * V'(I,X) + C(X) * U(X) * V(I,X) dx 
//      = Integral F(X) * V(I,X) dx
//
//    We note that the functions U(X) and U'(X) can be replaced by
//    the finite element form involving the linear sum of basis functions,
//    but we also note that the resulting integrand will only be nonzero
//    for terms where J = I - 1, I, or I + 1.
//
//    By writing this equation for basis functions I = 2 through N - 1,
//    and using the boundary conditions, we have N linear equations
//    for the N unknown coefficients U(1) through U(N), which can
//    be easily solved.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Input, double A ( double X ), evaluates a(x);
//
//    Input, double C ( double X ), evaluates c(x);
//
//    Input, double F ( double X ), evaluates f(x);
//
//    Input, double X[N], the mesh points.
//
//    Output, double FEM1D_BVP_LINEAR[N], the finite element coefficients, 
//    which are also the value of the computed solution at the mesh points.
//
    {
        const int QUAD_NUM = 2;

        double[] abscissa =
            {
                -0.577350269189625764509148780502,
                +0.577350269189625764509148780502
            }
            ;
        int e;
        int i;
        int ierror = 0;
        int j;
        double[] weight =
            {
                1.0, 1.0
            }
            ;
        //
//  Zero out the matrix and right hand side.
//
        double[] amat = typeMethods.r8mat_zero_new(n, n);
        double[] b = typeMethods.r8vec_zero_new(n);

        int e_num = n - 1;

        for (e = 0; e < e_num; e++)
        {
            int r = e + 1;

            double xl = x[e];
            double xr = x[r];

            int q;
            for (q = 0; q < QUAD_NUM; q++)
            {
                double xq = ((1.0 - abscissa[q]) * xl
                             + (1.0 + abscissa[q]) * xr)
                            / 2.0;

                double wq = weight[q] * (xr - xl) / 2.0;

                double vl = (xr - xq) / (xr - xl);
                double vlp = -1.0 / (xr - xl);

                double vr = (xq - xl) / (xr - xl);
                double vrp = +1.0 / (xr - xl);

                double axq = a(xq);
                double cxq = c(xq);
                double fxq = f(xq);

                amat[e + e * n] += wq * (vlp * axq * vlp + vl * cxq * vl);
                amat[e + r * n] += wq * (vlp * axq * vrp + vl * cxq * vr);
                b[e] += wq * (vl * fxq);

                amat[r + e * n] += wq * (vrp * axq * vlp + vr * cxq * vl);
                amat[r + r * n] += wq * (vrp * axq * vrp + vr * cxq * vr);
                b[r] += wq * (vr * fxq);
            }
        }

//
//  Equation 1 is the left boundary condition, U(0.0) = 0.0;
//
        for (j = 0; j < n; j++)
        {
            amat[0 + j * n] = 0.0;
        }

        b[0] = 0.0;
        for (i = 1; i < n; i++)
        {
            b[i] -= amat[i + 0 * n] * b[0];
        }

        for (i = 0; i < n; i++)
        {
            amat[i + 0 * n] = 0.0;
        }

        amat[0 + 0 * n] = 1.0;
//
//  Equation N is the right boundary condition, U(1.0) = 0.0;
//
        for (j = 0; j < n; j++)
        {
            amat[n - 1 + j * n] = 0.0;
        }

        b[n - 1] = 0.0;
        for (i = 0; i < n - 1; i++)
        {
            b[i] -= amat[i + (n - 1) * n] * b[n - 1];
        }

        for (i = 0; i < n; i++)
        {
            amat[i + (n - 1) * n] = 0.0;
        }

        amat[n - 1 + (n - 1) * n] = 1.0;
//
//  Solve the linear system.
//
        double[] u = typeMethods.r8mat_solve2(n, ref amat, ref b, ref ierror);

        return u;
    }

        
}