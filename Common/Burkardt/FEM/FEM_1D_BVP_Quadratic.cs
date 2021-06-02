using System;
using Burkardt.Types;

namespace Burkardt.FEM
{
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
            int QUAD_NUM = 3;

            double[] abscissa =  {
                -0.774596669241483377035853079956,
                0.000000000000000000000000000000,
                0.774596669241483377035853079956
            }
            ;
            double[] amat;
            double axq;
            double[] b;
            double cxq;
            int e_num;
            double fxq;
            int i;
            int ierror = 0;
            int l;
            int m;
            int quad_num = QUAD_NUM;
            int r;
            double[] u;
            double[] weight =  {
                0.555555555555555555555555555556,
                0.888888888888888888888888888889,
                0.555555555555555555555555555556
            }
            ;
            double wq;
            double vl;
            double vlp;
            double vm;
            double vmp;
            double vr;
            double vrp;
            double xl;
            double xm;
            double xq;
            double xr;
//
//  Zero out the matrix and right hand side.
//
            amat = typeMethods.r8mat_zero_new(n, n);
            b = typeMethods.r8vec_zero_new(n);
//
//  Integrate over element E.
//
            e_num = (n - 1) / 2;

            for (int e = 0; e < e_num; e++)
            {
//
//  Element E uses nodes
//    L = 2 * E
//    M = 2 * E + 1
//    R = 2 * E + 2
//
                l = 2 * e;
                m = 2 * e + 1;
                r = 2 * e + 2;

                xl = x[l];
                xm = x[m];
                xr = x[r];

                for (int q = 0; q < quad_num; q++)
                {

                    xq = ((1.0 - abscissa[q]) * xl
                          + (1.0 + abscissa[q]) * xr)
                         / 2.0;

                    wq = weight[q] * (xr - xl) / 2.0;

                    axq = a(xq);
                    cxq = c(xq);
                    fxq = f(xq);

                    vl = ((xq - xm) / (xl - xm))
                         * ((xq - xr) / (xl - xr));

                    vm = ((xq - xl) / (xm - xl))
                         * ((xq - xr) / (xm - xr));

                    vr = ((xq - xl) / (xr - xl))
                         * ((xq - xm) / (xr - xm));

                    vlp = (1.0 / (xl - xm))
                          * ((xq - xr) / (xl - xr))
                          + ((xq - xm) / (xl - xm))
                          * (1.0 / (xl - xr));

                    vmp = (1.0 / (xm - xl))
                          * ((xq - xr) / (xm - xr))
                          + ((xq - xl) / (xm - xl))
                          * (1.0 / (xm - xr));

                    vrp = (1.0 / (xr - xl))
                          * ((xq - xm) / (xr - xm))
                          + ((xq - xl) / (xr - xl))
                          * (1.0 / (xr - xm));

                    amat[l + l * n] = amat[l + l * n] + wq * (vlp * axq * vlp + vl * cxq * vl);
                    amat[l + m * n] = amat[l + m * n] + wq * (vlp * axq * vmp + vl * cxq * vm);
                    amat[l + r * n] = amat[l + r * n] + wq * (vlp * axq * vrp + vl * cxq * vr);
                    b[l] = b[l] + wq * (vl * fxq);

                    amat[m + l * n] = amat[m + l * n] + wq * (vmp * axq * vlp + vm * cxq * vl);
                    amat[m + m * n] = amat[m + m * n] + wq * (vmp * axq * vmp + vm * cxq * vm);
                    amat[m + r * n] = amat[m + r * n] + wq * (vmp * axq * vrp + vm * cxq * vr);
                    b[m] = b[m] + wq * (vm * fxq);

                    amat[r + l * n] = amat[r + l * n] + wq * (vrp * axq * vlp + vr * cxq * vl);
                    amat[r + m * n] = amat[r + m * n] + wq * (vrp * axq * vmp + vr * cxq * vm);
                    amat[r + r * n] = amat[r + r * n] + wq * (vrp * axq * vrp + vr * cxq * vr);
                    b[r] = b[r] + wq * (vr * fxq);
                }
            }

//
//  Equation 0 is the left boundary condition, U(0.0) = 0.0;
//
            i = 0;
            for (int j = 0; j < n; j++)
            {
                amat[i + j * n] = 0.0;
            }

            amat[i + i * n] = 1.0;
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
            u = typeMethods.r8mat_solve2(n, ref amat, ref b, ref ierror);

            return u;
        }
    }
}