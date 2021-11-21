using System;
using Burkardt.Types;

namespace Burkardt.FEM;

public static class FEM_1D_BVP_Quadratic
{
    public static double[] fem1d_bvp_quadratic(int n, Func<double, double> a , Func<double, double> c,
        Func<double, double> f, double[] x )
//****************************************************************************80
//
//  Purpose:
//
//    FEM1D_BVP_QUADRATIC solves a two point boundary value problem.
//
//  Discussion:
//
//    The program uses the finite element method, with piecewise quadratic basis
//    functions to solve a boundary value problem in one dimension.
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
//    Output, double FEM1D_BVP_QUADRATIC[N], the finite element coefficients, 
//    which are also the value of the computed solution at the mesh points.
//
    {
        const int QUAD_NUM = 3;

        double[] abscissa =  {
                -0.774596669241483377035853079956,
                0.000000000000000000000000000000,
                0.774596669241483377035853079956
            }
            ;
        int ierror = 0;
        double[] weight =  {
                0.555555555555555555555555555556,
                0.888888888888888888888888888889,
                0.555555555555555555555555555556
            }
            ;
        //
//  Zero out the matrix and right hand side.
//
        double[] amat = typeMethods.r8mat_zero_new(n, n);
        double[] b = typeMethods.r8vec_zero_new(n);
//
//  Integrate over element E.
//
        int e_num = (n - 1) / 2;

        for (int e = 0; e < e_num; e++)
        {
//
//  Element E uses nodes
//    L = 2 * E
//    M = 2 * E + 1
//    R = 2 * E + 2
//
            int l = 2 * e;
            int m = 2 * e + 1;
            int r = 2 * e + 2;

            double xl = x[l];
            double xm = x[m];
            double xr = x[r];

            for (int q = 0; q < QUAD_NUM; q++)
            {

                double xq = ((1.0 - abscissa[q]) * xl
                             + (1.0 + abscissa[q]) * xr)
                            / 2.0;

                double wq = weight[q] * (xr - xl) / 2.0;

                double axq = a(xq);
                double cxq = c(xq);
                double fxq = f(xq);

                double vl = (xq - xm) / (xl - xm)
                            * ((xq - xr) / (xl - xr));

                double vm = (xq - xl) / (xm - xl)
                            * ((xq - xr) / (xm - xr));

                double vr = (xq - xl) / (xr - xl)
                            * ((xq - xm) / (xr - xm));

                double vlp = 1.0 / (xl - xm)
                             * ((xq - xr) / (xl - xr))
                             + (xq - xm) / (xl - xm)
                             * (1.0 / (xl - xr));

                double vmp = 1.0 / (xm - xl)
                             * ((xq - xr) / (xm - xr))
                             + (xq - xl) / (xm - xl)
                             * (1.0 / (xm - xr));

                double vrp = 1.0 / (xr - xl)
                             * ((xq - xm) / (xr - xm))
                             + (xq - xl) / (xr - xl)
                             * (1.0 / (xr - xm));

                amat[l + l * n] += wq * (vlp * axq * vlp + vl * cxq * vl);
                amat[l + m * n] += wq * (vlp * axq * vmp + vl * cxq * vm);
                amat[l + r * n] += wq * (vlp * axq * vrp + vl * cxq * vr);
                b[l] += wq * (vl * fxq);

                amat[m + l * n] += wq * (vmp * axq * vlp + vm * cxq * vl);
                amat[m + m * n] += wq * (vmp * axq * vmp + vm * cxq * vm);
                amat[m + r * n] += wq * (vmp * axq * vrp + vm * cxq * vr);
                b[m] += wq * (vm * fxq);

                amat[r + l * n] += wq * (vrp * axq * vlp + vr * cxq * vl);
                amat[r + m * n] += wq * (vrp * axq * vmp + vr * cxq * vm);
                amat[r + r * n] += wq * (vrp * axq * vrp + vr * cxq * vr);
                b[r] += wq * (vr * fxq);
            }
        }

//
//  Equation 0 is the left boundary condition, U(0.0) = 0.0;
//
        int i = 0;
        for (int j = 0; j < n; j++)
        {
            amat[i + j * n] = 0.0;
        }

        amat[0] = 1.0;
        b[i] = 0.0;
//
//  Equation N-1 is the right boundary condition, U(1.0) = 0.0;
//
        i = n - 1;
        for (int j = 0; j < n; j++)
        {
            amat[i + j * n] = 0.0;
        }

        amat[i + i * n] = 1.0;
        b[i] = 0.0;
//
//  Solve the linear system.
//
        double[] u = typeMethods.r8mat_solve2(n, ref amat, ref b, ref ierror);

        return u;
    }
}