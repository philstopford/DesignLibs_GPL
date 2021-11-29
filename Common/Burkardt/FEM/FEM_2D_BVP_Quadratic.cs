using System;
using Burkardt.Types;

namespace Burkardt.FEM;

public static class FEM_2D_BVP_Quadratic
{
    public static double[] fem2d_bvp_quadratic(int nx, int ny, Func < double, double, double > a,
            Func <double, double, double> c, Func < double, double, double > f,
            double[] x, double[] y )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEM2D_BVP_LINEAR solves a boundary value problem on a rectangle.
        //
        //  Discussion:
        //
        //    The procedure uses the finite element method, with piecewise quadratic
        //    basis functions to solve a 2D boundary value problem over a rectangle
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
        //    Output, double FEM1D_BVP_QUADRATIC[NX*NY], the finite element coefficients, 
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
        int ex;
        int i;
        int ierror = 0;
        int ii;
        int j;
        int jj;
        int k;
        int[] node = new int[9];
        double[] v = new double[9];
        double[] vx = new double[9];
        double[] vy = new double[9];
        double[] weight =  {
                0.555555555555555555555555555556,
                0.888888888888888888888888888889,
                0.555555555555555555555555555556
            }
            ;
        double[] xx = new double[3];
        double[] yy = new double[3];

        int mn = nx * ny;

        double[] amat = typeMethods.r8mat_zero_new(mn, mn);
        double[] b = typeMethods.r8vec_zero_new(mn);

        int ex_num = (nx - 1) / 2;
        int ey_num = (ny - 1) / 2;

        for (ex = 0; ex < ex_num; ex++)
        {
            int w = 2 * ex;
            int cc = 2 * ex + 1;
            int e = 2 * ex + 2;

            xx[0] = x[w];
            xx[1] = x[cc];
            xx[2] = x[e];

            int ey;
            for (ey = 0; ey < ey_num; ey++)
            {
                int s = 2 * ey;
                int mm = 2 * ey + 1;
                int n = 2 * ey + 2;

                yy[0] = y[s];
                yy[1] = y[mm];
                yy[2] = y[n];
                //
                //  Node indices
                //
                //  7  8  9   wn cn en
                //  4  5  6   wm cm em
                //  1  2  3   ws cs es
                //
                node[0] = 2 * ey * nx + ex * 2 + 0;
                node[1] = 2 * ey * nx + ex * 2 + 1;
                node[2] = 2 * ey * nx + ex * 2 + 2;
                node[3] = (2 * ey + 1) * nx + ex * 2 + 0;
                node[4] = (2 * ey + 1) * nx + ex * 2 + 1;
                node[5] = (2 * ey + 1) * nx + ex * 2 + 2;
                node[6] = (2 * ey + 2) * nx + ex * 2 + 0;
                node[7] = (2 * ey + 2) * nx + ex * 2 + 1;
                node[8] = (2 * ey + 2) * nx + ex * 2 + 2;

                int qx;
                for (qx = 0; qx < QUAD_NUM; qx++)
                {
                    double xq = ((1.0 - abscissa[qx]) * xx[0]
                                 + (1.0 + abscissa[qx]) * xx[2])
                                / 2.0;

                    int qy;
                    for (qy = 0; qy < QUAD_NUM; qy++)
                    {
                        double yq = ((1.0 - abscissa[qy]) * yy[0]
                                     + (1.0 + abscissa[qy]) * yy[2])
                                    / 2.0;

                        double wq = weight[qx] * (xx[2] - xx[0]) / 2.0
                            * weight[qy] * (yy[2] - yy[0]) / 2.0;

                        k = 0;

                        int jl;
                        for (jl = 0; jl < 3; jl++)
                        {
                            int il;
                            for (il = 0; il < 3; il++)
                            {
                                v[k] = 1.0;
                                vx[k] = 0.0;
                                vy[k] = 0.0;
                                double t;
                                int il2;
                                int jl2;
                                for (il2 = 0; il2 < 3; il2++)
                                {
                                    if (il2 == il)
                                    {
                                        continue;
                                    }

                                    v[k] = v[k] * (xq - xx[il2]) / (xx[il] - xx[il2]);
                                    t = 1.0 / (xx[il] - xx[il2]);
                                    int il3;
                                    for (il3 = 0; il3 < 3; il3++)
                                    {
                                        if (il3 != il && il3 != il2)
                                        {
                                            t = t * (xq - xx[il3]) / (xx[il] - xx[il3]);
                                        }
                                    }

                                    for (jl2 = 0; jl2 < 3; jl2++)
                                    {
                                        if (jl2 != jl)
                                        {
                                            t = t * (yq - yy[jl2]) / (yy[jl] - yy[jl2]);
                                        }
                                    }

                                    vx[k] += t;
                                }

                                for (jl2 = 0; jl2 < 3; jl2++)
                                {
                                    if (jl2 == jl)
                                    {
                                        continue;
                                    }

                                    v[k] = v[k] * (yq - yy[jl2]) / (yy[jl] - yy[jl2]);
                                    t = 1.0 / (yy[jl] - yy[jl2]);
                                    for (il2 = 0; il2 < 3; il2++)
                                    {
                                        if (il2 != il)
                                        {
                                            t = t * (xq - xx[il2]) / (xx[il] - xx[il2]);
                                        }
                                    }

                                    int jl3;
                                    for (jl3 = 0; jl3 < 3; jl3++)
                                    {
                                        if (jl3 != jl && jl3 != jl2)
                                        {
                                            t = t * (yq - yy[jl3]) / (yy[jl] - yy[jl3]);
                                        }
                                    }

                                    vy[k] += t;
                                }

                                k += 1;
                            }
                        }

                        double aq = a(xq, yq);
                        double cq = c(xq, yq);
                        double 
                            fq = f(xq, yq);

                        for (i = 0; i < 9; i++)
                        {
                            ii = node[i];
                            for (j = 0; j < 9; j++)
                            {
                                jj = node[j];
                                amat[ii + jj * mn] += wq * (
                                    vx[i] * aq * vx[j]
                                    + vy[i] * aq * vy[j]
                                    + v[i] * cq * v[j]);
                            }

                            b[ii] += wq * (v[i] * fq);
                        }
                    }
                }
            }
        }

        //
        //  Where a node is on the boundary, 
        //  replace the finite element equation by a boundary condition.
        //
        k = 0;
        for (j = 0; j < ny; j++)
        {
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

    public static double fem2d_h1s_error_quadratic(int nx, int ny, double[] x, double[] y,
            double[] u, Func <double, double, double> exact_ux,
            Func <double, double, double> exact_uy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEM2D_H1S_ERROR_LINEAR: seminorm error of a finite element solution.
        //
        //  Discussion:
        //
        //    The finite element method has been used, over a rectangle,
        //    involving a grid of NX*NY nodes, with piecewise quadratic elements used 
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
        //    26 June 2014
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
        //    Output, double FEM2D_H1S_ERROR_QUADRATIC, the estimated seminorm of 
        //    the error.
        //
    {
        const int QUAD_NUM = 3;

        double[] abscissa =  {
                -0.774596669241483377035853079956,
                0.000000000000000000000000000000,
                0.774596669241483377035853079956
            }
            ;
        int ex;
        double h1s = 0;
        int[] node = new int[9];
        double[] vx = new double[9];
        double[] vy = new double[9];
        double[] weight =  {
                0.555555555555555555555555555556,
                0.888888888888888888888888888889,
                0.555555555555555555555555555556
            }
            ;
        double[] xx = new double[3];
        double[] yy = new double[3];

        int ex_num = (nx - 1) / 2;
        int ey_num = (ny - 1) / 2;

        for (ex = 0; ex < ex_num; ex++)
        {
            int w = 2 * ex;
            int cc = 2 * ex + 1;
            int e = 2 * ex + 2;

            xx[0] = x[w];
            xx[1] = x[cc];
            xx[2] = x[e];

            int ey;
            for (ey = 0; ey < ey_num; ey++)
            {
                int s = 2 * ey;
                int mm = 2 * ey + 1;
                int n = 2 * ey + 2;

                yy[0] = y[s];
                yy[1] = y[mm];
                yy[2] = y[n];
                //
                //  Node indices
                //
                //  7  8  9   wn cn en
                //  4  5  6   wm cm em
                //  1  2  3   ws cs es
                //
                node[0] = 2 * ey * nx + ex * 2 + 0;
                node[1] = 2 * ey * nx + ex * 2 + 1;
                node[2] = 2 * ey * nx + ex * 2 + 2;
                node[3] = (2 * ey + 1) * nx + ex * 2 + 0;
                node[4] = (2 * ey + 1) * nx + ex * 2 + 1;
                node[5] = (2 * ey + 1) * nx + ex * 2 + 2;
                node[6] = (2 * ey + 2) * nx + ex * 2 + 0;
                node[7] = (2 * ey + 2) * nx + ex * 2 + 1;
                node[8] = (2 * ey + 2) * nx + ex * 2 + 2;

                int qx;
                for (qx = 0; qx < QUAD_NUM; qx++)
                {
                    double xq = ((1.0 - abscissa[qx]) * xx[0]
                                 + (1.0 + abscissa[qx]) * xx[2])
                                / 2.0;

                    int qy;
                    for (qy = 0; qy < QUAD_NUM; qy++)
                    {
                        double yq = ((1.0 - abscissa[qy]) * yy[0]
                                     + (1.0 + abscissa[qy]) * yy[2])
                                    / 2.0;

                        double wq = weight[qx] * (xx[2] - xx[0]) / 2.0
                            * weight[qy] * (yy[2] - yy[0]) / 2.0;

                        double uxq = 0.0;
                        double 
                            uyq = 0.0;

                        int k = 0;

                        int jl;
                        for (jl = 0; jl < 3; jl++)
                        {
                            int il;
                            for (il = 0; il < 3; il++)
                            {
                                vx[k] = 0.0;
                                vy[k] = 0.0;
                                int jl2;
                                double t;
                                int il2;
                                for (il2 = 0; il2 < 3; il2++)
                                {
                                    if (il2 == il)
                                    {
                                        continue;
                                    }

                                    t = 1.0 / (xx[il] - xx[il2]);
                                    int il3;
                                    for (il3 = 0; il3 < 3; il3++)
                                    {
                                        if (il3 != il && il3 != il2)
                                        {
                                            t = t * (xq - xx[il3]) / (xx[il] - xx[il3]);
                                        }
                                    }

                                    for (jl2 = 0; jl2 < 3; jl2++)
                                    {
                                        if (jl2 != jl)
                                        {
                                            t = t * (yq - yy[jl2]) / (yy[jl] - yy[jl2]);
                                        }
                                    }

                                    vx[k] += t;
                                }

                                for (jl2 = 0; jl2 < 3; jl2++)
                                {
                                    if (jl2 == jl)
                                    {
                                        continue;
                                    }

                                    t = 1.0 / (yy[jl] - yy[jl2]);
                                    for (il2 = 0; il2 < 3; il2++)
                                    {
                                        if (il2 != il)
                                        {
                                            t = t * (xq - xx[il2]) / (xx[il] - xx[il2]);
                                        }
                                    }

                                    int jl3;
                                    for (jl3 = 0; jl3 < 3; jl3++)
                                    {
                                        if (jl3 != jl && jl3 != jl2)
                                        {
                                            t = t * (yq - yy[jl3]) / (yy[jl] - yy[jl3]);
                                        }
                                    }

                                    vy[k] += t;
                                }

                                uxq += u[node[k]] * vx[k];
                                uyq += u[node[k]] * vy[k];
                                k += 1;
                            }
                        }

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
        //    involving a grid of NX*NY nodes.
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

    public static double fem2d_l2_error_quadratic(int nx, int ny, double[] x, double[] y,
            double[] u, Func<double, double, double> exact )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEM2D_L2_ERROR_QUADRATIC: L2 error norm of a finite element solution.
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
        //    26 June 2014
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
        //    Output, double FEM2D_L2_ERROR_QUADRATIC, the estimated L2 norm of the error.
        //
    {
        const int QUAD_NUM = 3;

        double[] abscissa =  {
                -0.774596669241483377035853079956,
                0.000000000000000000000000000000,
                0.774596669241483377035853079956
            }
            ;
        double e2 = 0;
        int ex;
        int[] node = new int[9];
        double[] v = new double[9];
        double[] weight =  {
                0.555555555555555555555555555556,
                0.888888888888888888888888888889,
                0.555555555555555555555555555556
            }
            ;
        double[] xx = new double[3];
        double[] yy = new double[3];

        int ex_num = (nx - 1) / 2;
        int ey_num = (ny - 1) / 2;

        for (ex = 0; ex < ex_num; ex++)
        {
            int w = 2 * ex;
            int cc = 2 * ex + 1;
            int e = 2 * ex + 2;

            xx[0] = x[w];
            xx[1] = x[cc];
            xx[2] = x[e];

            int ey;
            for (ey = 0; ey < ey_num; ey++)
            {
                int s = 2 * ey;
                int mm = 2 * ey + 1;
                int n = 2 * ey + 2;

                yy[0] = y[s];
                yy[1] = y[mm];
                yy[2] = y[n];
                //
                //  Node indices
                //
                //  7  8  9   wn cn en
                //  4  5  6   wm cm em
                //  1  2  3   ws cs es
                //
                node[0] = 2 * ey * nx + ex * 2 + 0;
                node[1] = 2 * ey * nx + ex * 2 + 1;
                node[2] = 2 * ey * nx + ex * 2 + 2;
                node[3] = (2 * ey + 1) * nx + ex * 2 + 0;
                node[4] = (2 * ey + 1) * nx + ex * 2 + 1;
                node[5] = (2 * ey + 1) * nx + ex * 2 + 2;
                node[6] = (2 * ey + 2) * nx + ex * 2 + 0;
                node[7] = (2 * ey + 2) * nx + ex * 2 + 1;
                node[8] = (2 * ey + 2) * nx + ex * 2 + 2;

                int qx;
                for (qx = 0; qx < QUAD_NUM; qx++)
                {
                    double xq = ((1.0 - abscissa[qx]) * xx[0]
                                 + (1.0 + abscissa[qx]) * xx[2])
                                / 2.0;

                    int qy;
                    for (qy = 0; qy < QUAD_NUM; qy++)
                    {
                        double yq = ((1.0 - abscissa[qy]) * yy[0]
                                     + (1.0 + abscissa[qy]) * yy[2])
                                    / 2.0;

                        double wq = weight[qx] * (xx[2] - xx[0]) / 2.0
                            * weight[qy] * (yy[2] - yy[0]) / 2.0;

                        double uq = 0.0;
                        int k = 0;

                        int jl;
                        for (jl = 0; jl < 3; jl++)
                        {
                            int il;
                            for (il = 0; il < 3; il++)
                            {
                                v[k] = 1.0;
                                int il2;
                                for (il2 = 0; il2 < 3; il2++)
                                {
                                    if (il2 != il)
                                    {
                                        v[k] = v[k] * (xq - xx[il2]) / (xx[il] - xx[il2]);
                                    }
                                }

                                int jl2;
                                for (jl2 = 0; jl2 < 3; jl2++)
                                {
                                    if (jl2 != jl)
                                    {
                                        v[k] = v[k] * (yq - yy[jl2]) / (yy[jl] - yy[jl2]);
                                    }
                                }

                                uq += u[node[k]] * v[k];
                                k += 1;
                            }
                        }

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