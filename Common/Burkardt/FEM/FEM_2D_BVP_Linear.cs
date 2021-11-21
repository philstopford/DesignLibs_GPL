using System;
using Burkardt.Types;

namespace Burkardt.FEM;

public static class FEM_2D_BVP_Linear
{
    public static double[] fem2d_bvp_linear(int nx, int ny, Func<double, double, double> a,
            Func<double, double, double> c, Func<double, double, double> f,
            double[] x, double[] y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEM2D_BVP_LINEAR solves a boundary value problem on a rectangle.
        //
        //  Discussion:
        //
        //    The procedure uses the finite element method, with piecewise linear basis
        //    functions to solve a 2D boundary value problem over a rectangle
        //
        //    The following differential equation is imposed inside the region:
        //
        //      - d/dx a(x,y) du/dx - d/dy a(x,y) du/dy + c(x,y) * u(x,y) = f(x,y)
        //
        //    where a(x,y), c(x,y), and f(x,y) are given functions.
        //
        //    On the boundary, the solution is constrained to have the value 0.
        //
        //    The finite element method will use a regular grid of NX nodes in X, and 
        //    NY nodes in Y.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NX, NY, the number of X and Y grid values.
        //
        //    Input, double A ( double X, double Y ), evaluates a(x,y);
        //
        //    Input, double C ( double X, double Y ), evaluates c(x,y);
        //
        //    Input, double F ( double X, double Y ), evaluates f(x,y);
        //
        //    Input, double X[NX], Y[NY], the grid coordinates.
        //
        //    Output, double FEM1D_BVP_LINEAR[NX*NY], the finite element coefficients, 
        //    which are also the value of the computed solution at the mesh points.
        //
    {
        int QUAD_NUM = 3;

        double[] abscissa =
        {
            -0.774596669241483377035853079956,
            0.000000000000000000000000000000,
            0.774596669241483377035853079956
        };
        int ex;
        int ierror = 0;
        int ii;
        int j;
        int jj;
        double[] weight =
        {
            0.555555555555555555555555555556,
            0.888888888888888888888888888889,
            0.555555555555555555555555555556
        };

        int mn = nx * ny;

        double[] amat = new double[mn * mn];
        double[] b = new double[mn];

        for (jj = 0; jj < mn; jj++)
        {
            for (ii = 0; ii < mn; ii++)
            {
                amat[ii + jj * mn] = 0.0;
            }
        }

        for (ii = 0; ii < mn; ii++)
        {
            b[ii] = 0.0;
        }

        for (ex = 0; ex < nx - 1; ex++)
        {
            int w = ex;
            int e = ex + 1;

            double xw = x[w];
            double xe = x[e];

            int ey;
            for (ey = 0; ey < ny - 1; ey++)
            {
                int s = ey;
                int n = ey + 1;

                double ys = y[s];
                double yn = y[n];

                int sw = ey * nx + ex;
                int se = ey * nx + ex + 1;
                int nw = (ey + 1) * nx + ex;
                int ne = (ey + 1) * nx + ex + 1;

                int qx;
                for (qx = 0; qx < QUAD_NUM; qx++)
                {
                    double xq = ((1.0 - abscissa[qx]) * xw
                                 + (1.0 + abscissa[qx]) * xe)
                                / 2.0;

                    int qy;
                    for (qy = 0; qy < QUAD_NUM; qy++)
                    {
                        double yq = ((1.0 - abscissa[qy]) * ys
                                     + (1.0 + abscissa[qy]) * yn)
                                    / 2.0;

                        double wq = weight[qx] * (xe - xw) / 2.0
                            * weight[qy] * (yn - ys) / 2.0;

                        double aq = a(xq, yq);
                        double cq = c(xq, yq);
                        double fq = f(xq, yq);

                        double vsw = (xe - xq) / (xe - xw) * (yn - yq) / (yn - ys);
                        double vswx = -1.0 / (xe - xw) * (yn - yq) / (yn - ys);
                        double vswy = (xe - xq) / (xe - xw) * -1.0 / (yn - ys);

                        double vse = (xq - xw) / (xe - xw) * (yn - yq) / (yn - ys);
                        double vsex = 1.0 / (xe - xw) * (yn - yq) / (yn - ys);
                        double vsey = (xq - xw) / (xe - xw) * -1.0 / (yn - ys);

                        double vnw = (xe - xq) / (xe - xw) * (yq - ys) / (yn - ys);
                        double vnwx = -1.0 / (xe - xw) * (yq - ys) / (yn - ys);
                        double vnwy = (xe - xq) / (xe - xw) * 1.0 / (yn - ys);

                        double vne = (xq - xw) / (xe - xw) * (yq - ys) / (yn - ys);
                        double vnex = 1.0 / (xe - xw) * (yq - ys) / (yn - ys);
                        double vney = (xq - xw) / (xe - xw) * 1.0 / (yn - ys);

                        amat[sw + sw * mn] += wq * (vswx * aq * vswx
                                                    + vswy * aq * vswy
                                                    + vsw * cq * vsw);
                        amat[sw + se * mn] += wq * (vswx * aq * vsex
                                                    + vswy * aq * vsey
                                                    + vsw * cq * vse);
                        amat[sw + nw * mn] += wq * (vswx * aq * vnwx
                                                    + vswy * aq * vnwy
                                                    + vsw * cq * vnw);
                        amat[sw + ne * mn] += wq * (vswx * aq * vnex
                                                    + vswy * aq * vney
                                                    + vsw * cq * vne);
                        b[sw] += wq * (vsw * fq);

                        amat[se + sw * mn] += wq * (vsex * aq * vswx
                                                    + vsey * aq * vswy
                                                    + vse * cq * vsw);
                        amat[se + se * mn] += wq * (vsex * aq * vsex
                                                    + vsey * aq * vsey
                                                    + vse * cq * vse);
                        amat[se + nw * mn] += wq * (vsex * aq * vnwx
                                                    + vsey * aq * vnwy
                                                    + vse * cq * vnw);
                        amat[se + ne * mn] += wq * (vsex * aq * vnex
                                                    + vsey * aq * vney
                                                    + vse * cq * vne);
                        b[se] += wq * (vse * fq);

                        amat[nw + sw * mn] += wq * (vnwx * aq * vswx
                                                    + vnwy * aq * vswy
                                                    + vnw * cq * vsw);
                        amat[nw + se * mn] += wq * (vnwx * aq * vsex
                                                    + vnwy * aq * vsey
                                                    + vnw * cq * vse);
                        amat[nw + nw * mn] += wq * (vnwx * aq * vnwx
                                                    + vnwy * aq * vnwy
                                                    + vnw * cq * vnw);
                        amat[nw + ne * mn] += wq * (vnwx * aq * vnex
                                                    + vnwy * aq * vney
                                                    + vnw * cq * vne);
                        b[nw] += wq * (vnw * fq);

                        amat[ne + sw * mn] += wq * (vnex * aq * vswx
                                                    + vney * aq * vswy
                                                    + vne * cq * vsw);
                        amat[ne + se * mn] += wq * (vnex * aq * vsex
                                                    + vney * aq * vsey
                                                    + vne * cq * vse);
                        amat[ne + nw * mn] += wq * (vnex * aq * vnwx
                                                    + vney * aq * vnwy
                                                    + vne * cq * vnw);
                        amat[ne + ne * mn] += wq * (vnex * aq * vnex
                                                    + vney * aq * vney
                                                    + vne * cq * vne);
                        b[ne] += wq * (vne * fq);
                    }
                }
            }
        }

        //
        //  Where a node is on the boundary, 
        //  replace the finite element equation by a boundary condition.
        //
        int k = 0;
        for (j = 0; j < ny; j++)
        {
            int i;
            for (i = 0; i < nx; i++)
            {
                if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1)
                {
                    for (jj = 0; jj < mn; jj++)
                    {
                        amat[k + jj * mn] = 0.0;
                    }

                    for (ii = 0; ii < mn; ii++)
                    {
                        amat[ii + k * mn] = 0.0;
                    }

                    amat[k + k * mn] = 1.0;
                    b[k] = 0.0;
                }

                k += 1;
            }
        }

        //
        //  Solve the linear system.
        //
        double[] u = typeMethods.r8mat_solve2(mn, ref amat, ref b, ref ierror);

        return u;
    }

    public static double fem2d_h1s_error_linear(int nx, int ny, double[] x, double[] y,
            double[] u, Func<double, double, double> exact_ux,
            Func<double, double, double> exact_uy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEM2D_H1S_ERROR_LINEAR: seminorm error of a finite element solution.
        //
        //  Discussion:
        //
        //    The finite element method has been used, over a rectangle,
        //    involving a grid of NX*NY nodes, with piecewise linear elements used 
        //    for the basis.
        //
        //    The finite element solution U(x,y) has been computed, and formulas for the
        //    exact derivatives Vx(x,y) and Vy(x,y) are known.
        //
        //    This function estimates the H1 seminorm of the error:
        //
        //      H1S = sqrt ( integral ( x, y )   ( Ux(x,y) - Vx(x,y) )^2 
        //                                     + ( Uy(x,y) - Vy(x,y) )^2 dx dy )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NX, NY, the number of X and Y grid values.
        //
        //    Input, double X[NX], Y[NY], the grid coordinates.
        //
        //    Input, double U[NX*NY], the finite element coefficients.
        //
        //    Input, function EXACT_UX(X,Y), EXACT_UY(X,Y) return the 
        //    value of the derivatives of the exact solution with respect to
        //    X and Y, respectively, at the point (X,Y).
        //
        //    Output, double FEM2D_H1S_ERROR_LINEAR, the estimated seminorm of 
        //    the error.
        //
    {
        const int QUAD_NUM = 3;

        double[] abscissa =
        {
            -0.774596669241483377035853079956,
            0.000000000000000000000000000000,
            0.774596669241483377035853079956
        };
        int ex;
        double[] weight =
        {
            0.555555555555555555555555555556,
            0.888888888888888888888888888889,
            0.555555555555555555555555555556
        };

        double h1s = 0.0;

        for (ex = 0; ex < nx - 1; ex++)
        {
            int w = ex;
            int e = ex + 1;

            double xw = x[w];
            double xe = x[e];

            int ey;
            for (ey = 0; ey < ny - 1; ey++)
            {
                int s = ey;
                int n = ey + 1;

                double ys = y[s];
                double yn = y[n];

                int qx;
                for (qx = 0; qx < QUAD_NUM; qx++)
                {
                    double xq = ((1.0 - abscissa[qx]) * xw
                                 + (1.0 + abscissa[qx]) * xe)
                                / 2.0;

                    int qy;
                    for (qy = 0; qy < QUAD_NUM; qy++)
                    {
                        double yq = ((1.0 - abscissa[qy]) * ys
                                     + (1.0 + abscissa[qy]) * yn)
                                    / 2.0;

                        double wq = weight[qx] * (xe - xw) / 2.0
                            * weight[qy] * (yn - ys) / 2.0;

                        double vswx = -1.0 / (xe - xw) * (yn - yq) / (yn - ys);
                        double vswy = (xe - xq) / (xe - xw) * -1.0 / (yn - ys);

                        double vsex = 1.0 / (xe - xw) * (yn - yq) / (yn - ys);
                        double vsey = (xq - xw) / (xe - xw) * -1.0 / (yn - ys);

                        double vnwx = -1.0 / (xe - xw) * (yq - ys) / (yn - ys);
                        double vnwy = (xe - xq) / (xe - xw) * 1.0 / (yn - ys);

                        double vnex = 1.0 / (xe - xw) * (yq - ys) / (yn - ys);
                        double vney = (xq - xw) / (xe - xw) * 1.0 / (yn - ys);
                        //
                        //  Note that the south-west component of U is stored in U(W,S), not U(S,W)!
                        //
                        double uxq = u[w + s * nx] * vswx + u[e + s * nx] * vsex
                                                          + u[w + n * nx] * vnwx + u[e + n * nx] * vnex;
                        double uyq = u[w + s * nx] * vswy + u[e + s * nx] * vsey
                                                          + u[w + n * nx] * vnwy + u[e + n * nx] * vney;

                        double exq = exact_ux(xq, yq);
                        double eyq = exact_uy(xq, yq);

                        h1s += wq * (Math.Pow(uxq - exq, 2) + Math.Pow(uyq - eyq, 2));
                    }
                }
            }
        }

        h1s = Math.Sqrt(h1s);

        return h1s;
    }

    public static double fem2d_l1_error(int nx, int ny, double[] x, double[] y, double[] u,
            Func<double, double, double> exact)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEM2D_L1_ERROR estimates the l1 error norm of a finite element solution.
        //
        //  Discussion:
        //
        //    The finite element method has been used, over a rectangle,
        //    involving a grid of NX*NY nodes, with piecewise linear elements used 
        //    for the basis.
        //
        //    The finite element coefficients have been computed, and a formula for the
        //    exact solution is known.
        //
        //    This function estimates the little l1 norm of the error:
        //      E1 = sum ( 1 <= I <= NX, 1 <= J <= NY ) 
        //        abs ( U(i,j) - EXACT(X(i),Y(j)) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NX, NY, the number of X and Y grid values.
        //
        //    Input, double X[NX], Y[NY], the grid coordinates.
        //
        //    Input, double U[NX*NY], the finite element coefficients.
        //
        //    Input, function EQ = EXACT(X,Y), returns the value of the exact
        //    solution at the point (X,Y).
        //
        //    Output, double FEM2D_L1_ERROR, the little l1 norm of the error.
        //
    {
        int j;

        double e1 = 0.0;

        for (j = 0; j < ny; j++)
        {
            int i;
            for (i = 0; i < nx; i++)
            {
                e1 += Math.Abs(u[i + j * nx] - exact(x[i], y[j]));
            }
        }

        e1 = e1 / nx / ny;

        return e1;
    }

    public static double fem2d_l2_error_linear(int nx, int ny, double[] x, double[] y,
            double[] u, Func<double, double, double> exact)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEM2D_L2_ERROR_LINEAR: L2 error norm of a finite element solution.
        //
        //  Discussion:
        //
        //    The finite element method has been used, over a rectangle,
        //    involving a grid of NX*NY nodes, with piecewise linear elements used 
        //    for the basis.
        //
        //    The finite element coefficients have been computed, and a formula for the
        //    exact solution is known.
        //
        //    This function estimates E2, the L2 norm of the error:
        //
        //      E2 = Integral ( X, Y ) ( U(X,Y) - EXACT(X,Y) )^2 dX dY
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NX, NY, the number of X and Y grid values.
        //
        //    Input, double X[NX], Y[NY], the grid coordinates.
        //
        //    Input, double U[NX*NY], the finite element coefficients.
        //
        //    Input, function EQ = EXACT(X,Y), returns the value of the exact
        //    solution at the point (X,Y).
        //
        //    Output, double FEM2D_L2_ERROR_LINEAR, the estimated L2 norm of the error.
        //
    {
        const int QUAD_NUM = 3;

        double[] abscissa =
        {
            -0.774596669241483377035853079956,
            0.000000000000000000000000000000,
            0.774596669241483377035853079956
        };
        int ex;
        double[] weight =
        {
            0.555555555555555555555555555556,
            0.888888888888888888888888888889,
            0.555555555555555555555555555556
        };

        double e2 = 0.0;
        //
        //  Integrate over each interval.
        //
        for (ex = 0; ex < nx - 1; ex++)
        {
            int w = ex;
            int e = ex + 1;

            double xw = x[w];
            double xe = x[e];

            int ey;
            for (ey = 0; ey < ny - 1; ey++)
            {
                int s = ey;
                int n = ey + 1;

                double ys = y[s];
                double yn = y[n];

                int qx;
                for (qx = 0; qx < QUAD_NUM; qx++)
                {
                    double xq = ((1.0 - abscissa[qx]) * xw
                                 + (1.0 + abscissa[qx]) * xe)
                                / 2.0;

                    int qy;
                    for (qy = 0; qy < QUAD_NUM; qy++)
                    {
                        double yq = ((1.0 - abscissa[qy]) * ys
                                     + (1.0 + abscissa[qy]) * yn)
                                    / 2.0;

                        double wq = weight[qx] * (xe - xw) / 2.0
                            * weight[qy] * (yn - ys) / 2.0;

                        double vsw = (xe - xq) / (xe - xw) * (yn - yq) / (yn - ys);
                        double vse = (xq - xw) / (xe - xw) * (yn - yq) / (yn - ys);
                        double vnw = (xe - xq) / (xe - xw) * (yq - ys) / (yn - ys);
                        double vne = (xq - xw) / (xe - xw) * (yq - ys) / (yn - ys);
                        //
                        //  Note that the south-west component of U is stored in U(W,S), not U(S,W)!
                        //
                        double uq = u[w + s * nx] * vsw + u[e + s * nx] * vse
                                                        + u[w + n * nx] * vnw + u[e + n * nx] * vne;
                        double eq = exact(xq, yq);

                        e2 += wq * Math.Pow(uq - eq, 2);
                    }
                }
            }
        }

        e2 = Math.Sqrt(e2);

        return e2;
    }
}