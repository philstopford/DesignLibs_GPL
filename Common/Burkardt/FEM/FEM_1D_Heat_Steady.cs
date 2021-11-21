using System;
using Burkardt.Types;

namespace Burkardt.FEM;

public static class FEM_1D_Heat_Steady
{
    public static double[] fem1d_heat_steady(int n, double a, double b, double ua, double ub,
        Func<double, double> k, Func<double, double> f, double[] x)
//****************************************************************************80
//
//  Purpose:
//
//    FEM1D_HEAT_STEADY solves the steady 1D heat equation with finite elements.
//
//  Discussion:
//
//    The program uses the finite element method, with piecewise linear basis
//    functions to solve the steady state heat equation in one dimension.
//
//    The problem is defined on the region A <= x <= B.
//
//    The following differential equation is imposed between A and B:
//
//      - d/dx k(x) du/dx = f(x)
//
//    where k(x) and f(x) are given functions.
//
//    At the boundaries, the following conditions are applied:
//
//      u(A) = UA
//      u(B) = UB
//
//    A set of N equally spaced nodes is defined on this
//    interval, with A = X(1) < X(2) < ... < X(N) = B.
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
//        Integral K(X) * U'(X) * V'(I,X) dx = Integral F(X) * V(I,X) dx
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
//    09 April 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Input, double A, B, the left and right endpoints.
//
//    Input, double UA, UB, the prescribed value of U at A and B.
//
//    Input, double K ( double X ), evaluates k(x);
//
//    Input, double F ( double X ), evaluates f(x);
//
//    Input, double X[N], the mesh points.
//
//    Output, double FEM1D_HEAT_STEADY[N], the finite element coefficients, 
//    which are also the value of the computed solution at the mesh points.
//
    {
        const int QUAD_NUM = 2;

        double[] abscissa =  {
                -0.577350269189625764509148780502,
                +0.577350269189625764509148780502
            }
            ;
        int i;
        int ierror = 0;
        double[] weight =  {
                1.0, 1.0
            }
            ;
        //
//  Zero out the matrix and right hand side.
//
        double[] amat = typeMethods.r8mat_zero_new(n, n);
        double[] bvec = typeMethods.r8vec_zero_new(n);
//
//  Equation 1 is the left boundary condition, U(A) = UA;
//
        amat[0 + 0 * n] = 1.0;
        bvec[0] = ua;
//
//  Equation I involves the basis function at node I.
//  This basis function is nonzero from X(I-1) to X(I+1).
//  Equation I looks like this:
//
//    Integral A(X) U'(X) V'(I,X) dx = Integral F(X) V(I,X) dx
//
//  Then, we realize that U(X) = sum ( 1 <= J <= N ) U(J) * V(J,X), 
//  (U(X) means the function; U(J) is the coefficient of V(J,X) ).
//
//  The only V functions that are nonzero when V(I,X) is nonzero are
//  V(I-1,X) and V(I+1,X). 
//
//  Let's use the shorthand 
//
//    VL(X) = V(I-1,X)
//    VM(X) = V(I,X)
//    VR(X) = V(I+1,X)
//
//  So our equation becomes
//
//    Integral A(X) [ VL'(X) U(I-1) + VM'(X) U(I) + VR'(X) U(I+1) ] * VM'(X) dx
//  = Integral F(X) VM(X) dx.
//
//  
//
//  This is actually a set of N-2 linear equations for the N coefficients U.
//
//  Now gather the multipliers of U(I-1) to get the matrix entry A(I,I-1), 
//  and so on.
//
        for (i = 1; i < n - 1; i++)
        {
//
//  Get the left, right and middle coordinates.
//
            double xl = x[i - 1];
            double xm = x[i];
            double xr = x[i + 1];
//
//  Make temporary variables for A(I,I-1), A(I,I), A(I,I+1) and B(I).
//
            double al = 0.0;
            double am = 0.0;
            double ar = 0.0;
            double bm = 0.0;
//
//  We approximate the integrals by using a weighted sum of
//  the integrand values at quadrature points.
//
int q;
for (q = 0; q < QUAD_NUM; q++)
            {
//
//  Integrate over the LEFT interval, between XL and XM, where:
//
//  VL(X) = ( XM - X       ) / ( XM - XL )
//  VM(X) = (      X  - XL ) / ( XM - XL )
//  VR(X) = 0
//
//  VL'(X) =             - 1 / ( XM - XL )
//  VM'(X) =             + 1 / ( XM - XL ) 
//  VR'(X) = 0
//
                double xq = ((1.0 - abscissa[q]) * xl
                             + (1.0 + abscissa[q]) * xm)
                            / 2.0;

                double wq = weight[q] * (xm - xl) / 2.0;

//    vl =  ( xm - xq ) / ( xm - xl );
                double vlp = -1.0 / (xm - xl);

                double vm = (xq - xl) / (xm - xl);
                double vmp = +1.0 / (xm - xl);

//    vr =  0.0;
                double vrp = 0.0;

                double kxq = k(xq);
                double fxq = f(xq);

                al += wq * (kxq * vlp * vmp);
                am += wq * (kxq * vmp * vmp);
                ar += wq * (kxq * vrp * vmp);
                bm += wq * (fxq * vm);
//
//  Integrate over the RIGHT interval, between XM and XR, where:
//
//  VL(X) = 0
//  VM(X) = ( XR - X       ) / ( XR - XM )
//  VR(X) = (      X  - XM ) / ( XR - XM )
//
//  VL'(X) = 0
//  VM'(X) =             - 1 / ( XR - XM )
//  VR'(X) =             + 1 / ( XR - XM ) 
//
                xq = ((1.0 - abscissa[q]) * xm
                      + (1.0 + abscissa[q]) * xr)
                     / 2.0;

                wq = weight[q] * (xr - xm) / 2.0;

//    vl = 0.0;
                vlp = 0.0;

                vm = (xr - xq) / (xr - xm);
                vmp = -1.0 / (xr - xm);

//    vr = ( xq - xm ) / ( xr - xm );
                vrp = 1.0 / (xr - xm);

                kxq = k(xq);
                fxq = f(xq);

                al += wq * (kxq * vlp * vmp);
                am += wq * (kxq * vmp * vmp);
                ar += wq * (kxq * vrp * vmp);
                bm += wq * (fxq * vm);
            }

            amat[i + (i - 1) * n] = al;
            amat[i + i * n] = am;
            amat[i + (i + 1) * n] = ar;

            bvec[i] = bm;
        }

//
//  Equation N is the right boundary condition, U(B) = UB;
//
        amat[n - 1 + (n - 1) * n] = 1.0;
        bvec[n - 1] = ub;
//
//  Solve the linear system.
//
        double[] u = typeMethods.r8mat_solve2(n, ref amat, ref bvec, ref ierror);

        return u;
    }
}