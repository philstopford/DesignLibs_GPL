using System;
using Burkardt.Types;

namespace Burkardt.FEM
{
    public static class FEM_2D_BVP_Serene
    {
        public static double[] basis_serene(double xq, double yq, double xw, double ys, double xe,
            double yn, double[] xx, double[] yy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_SERENE evaluates the serendipity basis functions.
        //
        //  Discussion:
        //
        //    This procedure assumes that a serendipity element has been defined,
        //    whose sides are parallel to coordinate axes.
        //
        //    The local element numbering is
        //
        //      YN  3--2--1
        //       |  |     |
        //       |  4     8
        //       |  |     |
        //      YS  5--6--7
        //       |
        //       +--XW---XE--
        //
        //    We note that each basis function can be written as the product of
        //    three linear terms, which never result in an x^2y^2 term.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double XQ, YQ, the evaluation point.
        //
        //    Input, double XW, YS, the coordinates of the lower left corner.
        //
        //    Input, double XE, YN, the coordinates of the upper right corner.
        //
        //    Input, double XX[8], YY[8], the coordinates of the 8 nodes.
        //
        //    Output, double BASIS_SERENE[8], the value of the basis functions 
        //    at (XQ,YQ).
        //
        {
            double[] v;

            v = new double[8];

            v[0] = not1(xq, xw, xx[0])
                   * not1(yq, ys, yy[0])
                   * not2(xq, yq, xx[7], yy[7], xx[1], yy[1], xx[0], yy[0]);

            v[1] = not1(xq, xw, xx[1])
                   * not1(xq, xe, xx[1])
                   * not1(yq, ys, yy[1]);

            v[2] = not1(xq, xe, xx[2])
                   * not1(yq, ys, yy[2])
                   * not2(xq, yq, xx[1], yy[1], xx[3], yy[3], xx[2], yy[2]);

            v[3] = not1(xq, xe, xx[3])
                   * not1(yq, yn, yy[3])
                   * not1(yq, ys, yy[3]);

            v[4] = not1(xq, xe, xx[4])
                   * not1(yq, yn, yy[4])
                   * not2(xq, yq, xx[3], yy[3], xx[5], yy[5], xx[4], yy[4]);

            v[5] = not1(xq, xe, xx[5])
                   * not1(xq, xw, xx[5])
                   * not1(yq, yn, yy[5]);

            v[6] = not1(xq, xw, xx[6])
                   * not1(yq, yn, yy[6])
                   * not2(xq, yq, xx[5], yy[5], xx[7], yy[7], xx[6], yy[6]);

            v[7] = not1(yq, ys, yy[7])
                   * not1(yq, yn, yy[7])
                   * not1(xq, xw, xx[7]);

            return v;
        }
        //****************************************************************************80

        public static double[] basis_dx_serene(double xq, double yq, double xw, double ys,
            double xe, double yn, double[] xx, double[] yy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_DX_SERENE differentiates the serendipity basis functions for X.
        //
        //  Discussion:
        //
        //    This procedure assumes that a serendipity element has been defined,
        //    whose sides are parallel to coordinate axes.
        //
        //    The local element numbering is
        //
        //      YN  3--2--1
        //       |  |     |
        //       |  4     8
        //       |  |     |
        //      YS  5--6--7
        //       |
        //       +--XW---XE--
        //
        //    We note that each basis function can be written as the product of
        //    three linear terms, which never result in an x^2y^2 term.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double XQ, YQ, the evaluation point.
        //
        //    Input, double XW, YS, the coordinates of the lower left corner.
        //
        //    Input, double XE, YN, the coordinates of the upper right corner.
        //
        //    Input, double XX[8], YY[8], the coordinates of the 8 nodes.
        //
        //    Output, double BASIS_DX_SERENE[8], the derivatives of the basis 
        //    functions at (XQ,YQ) with respect to X.
        //
        {
            double[] vx;

            vx = new double[8];

            vx[0] =
                not1d(xw, xx[0])
                * not1(yq, ys, yy[0])
                * not2(xq, yq, xx[7], yy[7], xx[1], yy[1], xx[0], yy[0])
                + not1(xq, xw, xx[0])
                * not1(yq, ys, yy[0])
                * not2dx(xx[7], yy[7], xx[1], yy[1], xx[0], yy[0]);

            vx[1] =
                not1d(xw, xx[1])
                * not1(xq, xe, xx[1])
                * not1(yq, ys, yy[1])
                + not1(xq, xw, xx[1])
                * not1d(xe, xx[1])
                * not1(yq, ys, yy[1]);

            vx[2] =
                not1d(xe, xx[2])
                * not1(yq, ys, yy[2])
                * not2(xq, yq, xx[1], yy[1], xx[3], yy[3], xx[2], yy[2])
                + not1(xq, xe, xx[2])
                * not1(yq, ys, yy[2])
                * not2dx(xx[1], yy[1], xx[3], yy[3], xx[2], yy[2]);

            vx[3] =
                not1d(xe, xx[3])
                * not1(yq, yn, yy[3])
                * not1(yq, ys, yy[3]);

            vx[4] =
                not1d(xe, xx[4])
                * not1(yq, yn, yy[4])
                * not2(xq, yq, xx[3], yy[3], xx[5], yy[5], xx[4], yy[4])
                + not1(xq, xe, xx[4])
                * not1(yq, yn, yy[4])
                * not2dx(xx[3], yy[3], xx[5], yy[5], xx[4], yy[4]);

            vx[5] =
                not1d(xe, xx[5])
                * not1(xq, xw, xx[5])
                * not1(yq, yn, yy[5])
                + not1(xq, xe, xx[5])
                * not1d(xw, xx[5])
                * not1(yq, yn, yy[5]);

            vx[6] =
                not1d(xw, xx[6])
                * not1(yq, yn, yy[6])
                * not2(xq, yq, xx[5], yy[5], xx[7], yy[7], xx[6], yy[6])
                + not1(xq, xw, xx[6])
                * not1(yq, yn, yy[6])
                * not2dx(xx[5], yy[5], xx[7], yy[7], xx[6], yy[6]);

            vx[7] =
                not1(yq, ys, yy[7])
                * not1(yq, yn, yy[7])
                * not1d(xw, xx[7]);

            return vx;
        }
        //****************************************************************************80

        public static double[] basis_dy_serene(double xq, double yq, double xw, double ys,
            double xe, double yn, double[] xx, double[] yy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_DY_SERENE differentiates the serendipity basis functions for Y.
        //
        //  Discussion:
        //
        //    This procedure assumes that a serendipity element has been defined,
        //    whose sides are parallel to coordinate axes.
        //
        //    The local element numbering is
        //
        //      YN  3--2--1
        //       |  |     |
        //       |  4     8
        //       |  |     |
        //      YS  5--6--7
        //       |
        //       +--XW---XE--
        //
        //    We note that each basis function can be written as the product of
        //    three linear terms, which never result in an x^2y^2 term.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double XQ, YQ, the evaluation point.
        //
        //    Input, double XW, YS, the coordinates of the lower left corner.
        //
        //    Input, double XE, YN, the coordinates of the upper right corner.
        //
        //    Input, double XX[8], YY[8], the coordinates of the 8 nodes.
        //
        //    Output, double BASIS_DY_SERENE[8], the derivatives of the basis 
        //    functions at (XQ,YQ) with respect to Y.
        //
        {
            double[] vy;

            vy = new double[8];

            vy[0] =
                not1(xq, xw, xx[0])
                * not1d(ys, yy[0])
                * not2(xq, yq, xx[7], yy[7], xx[1], yy[1], xx[0], yy[0])
                + not1(xq, xw, xx[0])
                * not1(yq, ys, yy[0])
                * not2dy(xx[7], yy[7], xx[1], yy[1], xx[0], yy[0]);

            vy[1] =
                not1(xq, xw, xx[1])
                * not1(xq, xe, xx[1])
                * not1d(ys, yy[1]);

            vy[2] = not1(xq, xe, xx[2])
                    * not1d(ys, yy[2])
                    * not2(xq, yq, xx[1], yy[1], xx[3], yy[3], xx[2], yy[2])
                    + not1(xq, xe, xx[2])
                    * not1(yq, ys, yy[2])
                    * not2dy(xx[1], yy[1], xx[3], yy[3], xx[2], yy[2]);

            vy[3] =
                not1(xq, xe, xx[3])
                * not1d(yn, yy[3])
                * not1(yq, ys, yy[3])
                + not1(xq, xe, xx[3])
                * not1(yq, yn, yy[3])
                * not1d(ys, yy[3]);

            vy[4] =
                not1(xq, xe, xx[4])
                * not1d(yn, yy[4])
                * not2(xq, yq, xx[3], yy[3], xx[5], yy[5], xx[4], yy[4])
                + not1(xq, xe, xx[4])
                * not1(yq, yn, yy[4])
                * not2dy(xx[3], yy[3], xx[5], yy[5], xx[4], yy[4]);

            vy[5] =
                not1(xq, xe, xx[5])
                * not1(xq, xw, xx[5])
                * not1d(yn, yy[5]);

            vy[6] =
                not1(xq, xw, xx[6])
                * not1d(yn, yy[6])
                * not2(xq, yq, xx[5], yy[5], xx[7], yy[7], xx[6], yy[6])
                + not1(xq, xw, xx[6])
                * not1(yq, yn, yy[6])
                * not2dy(xx[5], yy[5], xx[7], yy[7], xx[6], yy[6]);

            vy[7] =
                not1d(ys, yy[7])
                * not1(yq, yn, yy[7])
                * not1(xq, xw, xx[7])
                + not1(yq, ys, yy[7])
                * not1d(yn, yy[7])
                * not1(xq, xw, xx[7]);

            return vy;
        }

        public static double[] fem2d_bvp_serene(int nx, int ny, Func<double,double,double> a ,
        Func<double,double,double> c, Func<double,double,double> f,
        double[] x, double[] y, bool show11 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEM2D_BVP_SERENE solves boundary value problem on a rectangle.
        //
        //  Discussion:
        //
        //    The program uses the finite element method, with piecewise 
        //    serendipity basis functions to solve a 2D boundary value problem 
        //    over a rectangle.
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
        //    NY nodes in Y.  Both NX and NY must be odd.
        //
        //    The local element numbering is
        //
        //      3--2--1
        //      |     |
        //      4     8
        //      |     |
        //      5--6--7
        //
        //    The serendipity element mass matrix is a multiple of:
        //
        //       6.0, -6.0,  2.0, -8.0,  3.0, -8.0,  2.0, -6.0
        //      -6.0, 32.0, -6.0, 20.0, -8.0, 16.0, -8.0, 20.0
        //       2.0, -6.0,  6.0, -6.0,  2.0, -8.0,  3.0, -8.0
        //      -8.0, 20.0, -6.0, 32.0, -6.0, 20.0, -8.0, 16.0
        //       3.0, -8.0,  2.0, -6.0,  6.0, -6.0,  2.0, -8.0
        //      -8.0, 16.0, -8.0, 20.0, -6.0, 32.0, -6.0, 20.0
        //       2.0, -8.0,  3.0, -8.0,  2.0, -6.0,  6.0, -6.0
        //      -6.0, 20.0, -8.0, 16.0, -8.0, 20.0, -6.0, 32.0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NX, NY, the number of X and Y grid values.
        //    NX and NY must be odd and at least 3.
        //
        //    Input, function A(X,Y), evaluates a(x,y)
        //
        //    Input, function C(X,Y), evaluates c(x,y)
        //
        //    Input, function F(X,Y), evaluates f(x,y)
        //
        //    Input, double X[NX], Y[NY], the mesh points.
        //
        //    Input, bool SHOW11, is true to print out the element matrix
        //    for the element in row 1, column 1.
        //
        //    Output, double FEM2D_BVP_SERENE[MN], the finite element coefficients, which 
        //    are also the value of the computed solution at the mesh points.
        //
        {
            int QUAD_NUM = 3;

            double[] abscissa =  {
                -0.774596669241483377035853079956,
                0.000000000000000000000000000000,
                0.774596669241483377035853079956
            }
            ;
            double[] ae = new double[1];
            double[] amat;
            double aq;
            double[] b;
            double[] be = new double[1];
            int cc;
            double cq;
            int e;
            int ex;
            int ex_num;
            int ey;
            int ey_num;
            double fq;
            int i;
            int ierror = 0;
            int ii;
            int inc;
            int j;
            int jj;
            int k;
            int mm;
            int mn;
            int n;
            int[] node = new int[8];
            int quad_num = QUAD_NUM;
            int qx;
            int qy;
            int s;
            double scale;
            double[] u;
            double[] v;
            double[] vx;
            double[] vy;
            int w;
            double[] weight =  {
                0.555555555555555555555555555556,
                0.888888888888888888888888888889,
                0.555555555555555555555555555556
            }
            ;
            double wq;
            double xe;
            double xq;
            double xw;
            double[] xx = new double[8];
            double yn;
            double yq;
            double ys;
            double[] yy = new double[8];
            //
            //  Make room for the matrix A and right hand side b.
            //
            mn = fem2d_bvp_serene_node_num(nx, ny);

            amat = typeMethods.r8mat_zero_new(mn, mn);
            b = typeMethods.r8vec_zero_new(mn);
            //
            //  Compute the matrix entries by integrating over each element.
            //
            ex_num = (nx - 1) / 2;
            ey_num = (ny - 1) / 2;

            for (ey = 0; ey < ey_num; ey++)
            {
                s = 2 * ey;
                mm = 2 * ey + 1;
                n = 2 * ey + 2;

                ys = y[s];
                yn = y[n];

                yy[0] = y[n];
                yy[1] = y[n];
                yy[2] = y[n];
                yy[3] = y[mm];
                yy[4] = y[s];
                yy[5] = y[s];
                yy[6] = y[s];
                yy[7] = y[mm];

                for (ex = 0; ex < ex_num; ex++)
                {
                    w = 2 * ex;
                    cc = 2 * ex + 1;
                    e = 2 * ex + 2;

                    xe = x[e];
                    xw = x[w];

                    xx[0] = x[e];
                    xx[1] = x[cc];
                    xx[2] = x[w];
                    xx[3] = x[w];
                    xx[4] = x[w];
                    xx[5] = x[cc];
                    xx[6] = x[e];
                    xx[7] = x[e];
                    //
                    //  Node indices
                    //
                    //  3  2  1  wn  cn  en
                    //  4     8  wm      em
                    //  5  6  7  ws  cs  es
                    //
                    node[0] = (3 * ey + 3) * ey_num + 2 * ey + 2 * ex + 4;
                    node[1] = (3 * ey + 3) * ey_num + 2 * ey + 2 * ex + 3;
                    node[2] = (3 * ey + 3) * ey_num + 2 * ey + 2 * ex + 2;
                    node[3] = (3 * ey + 2) * ey_num + 2 * ey + ex + 1;
                    node[4] = (3 * ey) * ey_num + 2 * ey + 2 * ex;
                    node[5] = (3 * ey) * ey_num + 2 * ey + 2 * ex + 1;
                    node[6] = (3 * ey) * ey_num + 2 * ey + 2 * ex + 2;
                    node[7] = (3 * ey + 2) * ey_num + 2 * ey + ex + 2;

                    if (show11)
                    {
                        if (ey == 0 && ex == 0)
                        {
                            ae = typeMethods.r8mat_zero_new(8, 8);
                            be = typeMethods.r8vec_zero_new(8);
                        }
                    }

                    for (qx = 0; qx < quad_num; qx++)
                    {
                        xq = ((1.0 - abscissa[qx]) * x[e]
                              + (1.0 + abscissa[qx]) * x[w])
                             / 2.0;

                        for (qy = 0; qy < quad_num; qy++)
                        {
                            yq = ((1.0 - abscissa[qy]) * y[n]
                                  + (1.0 + abscissa[qy]) * y[s])
                                 / 2.0;

                            wq = weight[qx] * (x[e] - x[w]) / 2.0
                                * weight[qy] * (y[n] - y[s]) / 2.0;

                            v = basis_serene(xq, yq, xw, ys, xe, yn, xx, yy);
                            vx = basis_dx_serene(xq, yq, xw, ys, xe, yn, xx, yy);
                            vy = basis_dy_serene(xq, yq, xw, ys, xe, yn, xx, yy);

                            aq = a(xq, yq);
                            cq = c(xq, yq);
                            fq = f(xq, yq);
                            //
                            //  Build the element matrix.
                            //
                            if (show11)
                            {
                                if (ey == 0 && ex == 0)
                                {
                                    for (i = 0; i < 8; i++)
                                    {
                                        for (j = 0; j < 8; j++)
                                        {
                                            ae[i + j * 8] = ae[i + j * 8] + wq * (vx[i] * aq * vx[j]
                                                + vy[i] * aq * vy[j]
                                                + v[i] * cq * v[j]);
                                        }

                                        be[i] = be[i] + wq * (v[i] * fq);
                                    }
                                }
                            }

                            for (i = 0; i < 8; i++)
                            {
                                ii = node[i];
                                for (j = 0; j < 8; j++)
                                {
                                    jj = node[j];
                                    amat[ii + jj * mn] = amat[ii + jj * mn] + wq * (vx[i] * aq * vx[j]
                                        + vy[i] * aq * vy[j]
                                        + v[i] * cq * v[j]);
                                }

                                b[ii] = b[ii] + wq * (v[i] * fq);
                            }
                        }
                    }

                    //
                    //  Print a sample element matrix.
                    //
                    if (show11)
                    {
                        if (ey == 0 && ex == 0)
                        {
                            scale = 0.5 * ae[0 + 2 * 8];
                            for (j = 0; j < 8; j++)
                            {
                                for (i = 0; i < 8; i++)
                                {
                                    ae[i + j * 8] = ae[i + j * 8] / scale;
                                }
                            }

                            typeMethods.r8mat_print(8, 8, ae, "  Wathen elementary mass matrix:");
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
                if ((j % 2) == 0)
                {
                    inc = 1;
                }
                else
                {
                    inc = 2;
                }

                for (i = 0; i < nx; i = i + inc)
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

                    k = k + 1;
                }
            }

            //
            //  Solve the linear system.
            //
            u = typeMethods.r8mat_solve2(mn, ref amat, ref b, ref ierror);

            return u;
        }

        public static int fem2d_bvp_serene_node_num(int nx, int ny)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FEM2D_BVP_SERENE_NODE_NUM counts the number of nodes.
            //
            //  Discussion:
            //
            //    The program uses the finite element method, with piecewise serendipity 
            //    basis functions to solve a 2D boundary value problem over a rectangle.
            //
            //    The grid uses NX nodes in the X direction and NY nodes in the Y direction.
            //
            //    Both NX and NY must be odd.
            //
            //    Because of the peculiar shape of the serendipity elements, counting the
            //    number of nodes and variables is a little tricky.  Here is a grid for
            //    the case when NX = 7 and NY = 5, for which there are 29 nodes 
            //    and variables.
            //
            //     23 24 25 26 27 28 29
            //     19    20    21    22
            //     12 13 14 15 16 17 18
            //      8     9    10    11
            //      1  2  3  4  5  6  7
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 July 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NX, NY, the number of X and Y grid values.
            //    NX and NY must be odd and at least 3.
            //
            //    Output, int FEM2D_BVP_SERENE_NODE_NUM, the number of nodes and variables.
            //
        {
            int value;

            value =
                nx * (ny + 1) / 2
                + (nx + 1) / 2 * (ny - 1) / 2;

            return value;
        }

        public static double fem2d_h1s_error_serene(int nx, int ny, double[] x, double[] y,
        double[] u, Func<double, double, double> exact_ux,
        Func<double, double, double> exact_uy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEM2D_H1S_ERROR_SERENE: seminorm error of a finite element solution.
        //
        //  Discussion:
        //
        //    We assume the finite element method has been used, over a product region
        //    involving a grid of NX*NY nodes, with serendipity elements used 
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
        //    07 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NX, NY, the number of nodes.
        //
        //    Input, double X[NX], Y[NY], the grid coordinates.
        //
        //    Input, double U[*], the finite element coefficients.
        //
        //    Input, function EQX = EXACT_UX(X,Y), function EQY = EXACT_UY(X,Y)
        //    returns the exact derivatives with respect to X and Y.
        //
        //    Output, double FEM2D_BVP_H1S, the estimated seminorm of the error.
        //
        {
            int QUAD_NUM = 3;

            double[] abscissa =  {
                -0.774596669241483377035853079956,
                0.000000000000000000000000000000,
                0.774596669241483377035853079956
            }
            ;
            int cc;
            int e;
            int ex;
            int ex_num;
            int ey;
            int ey_num;
            double exq;
            double eyq;
            double h1s;
            int k;
            int mm;
            int n;
            int[] node = new int[8];
            int quad_num = QUAD_NUM;
            int qx;
            int qy;
            int s;
            double uxq;
            double uyq;
            double[] vx;
            double[] vy;
            int w;
            double[] weight =  {
                0.555555555555555555555555555556,
                0.888888888888888888888888888889,
                0.555555555555555555555555555556
            }
            ;
            double wq;
            double xe;
            double xq;
            double xw;
            double[] xx = new double[8];
            double yn;
            double yq;
            double ys;
            double[] yy = new double[8];

            h1s = 0.0;
            //
            //  Quadrature definitions.
            //
            ex_num = (nx - 1) / 2;
            ey_num = (ny - 1) / 2;

            for (ey = 0; ey < ey_num; ey++)
            {
                s = 2 * ey;
                mm = 2 * ey + 1;
                n = 2 * ey + 2;

                ys = y[s];
                yn = y[n];

                yy[0] = y[n];
                yy[1] = y[n];
                yy[2] = y[n];
                yy[3] = y[mm];
                yy[4] = y[s];
                yy[5] = y[s];
                yy[6] = y[s];
                yy[7] = y[mm];

                for (ex = 0; ex < ex_num; ex++)
                {
                    w = 2 * ex;
                    cc = 2 * ex + 1;
                    e = 2 * ex + 2;

                    xe = x[e];
                    xw = x[w];

                    xx[0] = x[e];
                    xx[1] = x[cc];
                    xx[2] = x[w];
                    xx[3] = x[w];
                    xx[4] = x[w];
                    xx[5] = x[cc];
                    xx[6] = x[e];
                    xx[7] = x[e];
                    //
                    //  Node indices
                    //
                    //  3  2  1  wn  cn  en
                    //  4     8  wm      em
                    //  5  6  7  ws  cs  es
                    //
                    node[0] = (3 * ey + 3) * ey_num + 2 * ey + 2 * ex + 4;
                    node[1] = (3 * ey + 3) * ey_num + 2 * ey + 2 * ex + 3;
                    node[2] = (3 * ey + 3) * ey_num + 2 * ey + 2 * ex + 2;
                    node[3] = (3 * ey + 2) * ey_num + 2 * ey + ex + 1;
                    node[4] = (3 * ey) * ey_num + 2 * ey + 2 * ex;
                    node[5] = (3 * ey) * ey_num + 2 * ey + 2 * ex + 1;
                    node[6] = (3 * ey) * ey_num + 2 * ey + 2 * ex + 2;
                    node[7] = (3 * ey + 2) * ey_num + 2 * ey + ex + 2;

                    for (qx = 0; qx < quad_num; qx++)
                    {
                        xq = ((1.0 - abscissa[qx]) * x[e]
                              + (1.0 + abscissa[qx]) * x[w])
                             / 2.0;

                        for (qy = 0; qy < quad_num; qy++)
                        {
                            yq = ((1.0 - abscissa[qy]) * y[n]
                                  + (1.0 + abscissa[qy]) * y[s])
                                 / 2.0;

                            wq = weight[qx] * (x[e] - x[w]) / 2.0
                                * weight[qy] * (y[n] - y[s]) / 2.0;

                            vx = basis_dx_serene(xq, yq, xw, ys, xe, yn, xx, yy);
                            vy = basis_dy_serene(xq, yq, xw, ys, xe, yn, xx, yy);

                            uxq = 0.0;
                            uyq = 0.0;
                            for (k = 0; k < 8; k++)
                            {
                                uxq = uxq + u[node[k]] * vx[k];
                                uyq = uyq + u[node[k]] * vy[k];
                            }

                            exq = exact_ux(xq, yq);
                            eyq = exact_uy(xq, yq);

                            h1s = h1s + wq * (Math.Pow(uxq - exq, 2) + Math.Pow(uyq - eyq, 2));
                        }
                    }
                }
            }

            h1s = Math.Sqrt(h1s);

            return h1s;
        }

        public static double fem2d_l1_error_serene(int nx, int ny, double[] x, double[] y,
        double[] u, Func<double, double, double> exact )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEM2D_L1_ERROR_SERENE: l1 error norm of a finite element solution.
        //
        //  Discussion:
        //
        //    We assume the finite element method has been used, over a product
        //    region with NX*NY nodes and the serendipity element.
        //
        //    The coefficients U(*) have been computed, and a formula for the
        //    exact solution is known.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 July 2014
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
        //    Input, double U[*], the finite element coefficients.
        //
        //    Input, function EQ = EXACT(X,Y), returns the value of the exact
        //    solution at the point (X,Y).
        //
        //    Output, double FEM2D_L1_ERROR_SERENE, the little l1 norm of the error.
        //
        {
            int i;
            int inc;
            int j;
            int k;
            double e1;

            e1 = 0.0;
            k = 0;

            for (j = 0; j < ny; j++)
            {
                if ((j % 2) == 0)
                {
                    inc = 1;
                }
                else
                {
                    inc = 2;
                }

                for (i = 0; i < nx; i = i + inc)
                {
                    e1 = e1 + Math.Abs(u[k] - exact(x[i], y[j]));
                    k = k + 1;
                }
            }

            e1 = e1 / (double) (k);

            return e1;
        }

        public static double fem2d_l2_error_serene(int nx, int ny, double[] x, double[] y,
        double[] u, Func<double, double, double> exact )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEM2D_L2_ERROR_SERENE: L2 error norm of a finite element solution.
        //
        //  Discussion:
        //
        //    The finite element method has been used, over a rectangle,
        //    involving a grid of NXxNY nodes, with serendipity elements used 
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
        //    07 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NX, NY, the number of nodes in the X 
        //    and Y directions.
        //
        //    Input, double X[NX], Y[NY], the grid coordinates.
        //
        //    Input, double U[*], the finite element coefficients.
        //
        //    Input, function EQ = EXACT(X,Y), returns the value of the exact
        //    solution at the point (X,Y).
        //
        //    Output, double E2, the estimated L2 norm of the error.
        //
        {
            int QUAD_NUM = 3;

            double[] abscissa =  {
                -0.774596669241483377035853079956,
                0.000000000000000000000000000000,
                0.774596669241483377035853079956
            }
            ;
            int cc;
            int e;
            double e2;
            double eq;
            int ex;
            int ex_num;
            int ey;
            int ey_num;
            int k;
            int mm;
            int n;
            int[] node = new int[8];
            int quad_num = QUAD_NUM;
            int qx;
            int qy;
            int s;
            double uq;
            double[] v;
            int w;
            double[] weight =  {
                0.555555555555555555555555555556,
                0.888888888888888888888888888889,
                0.555555555555555555555555555556
            }
            ;
            double wq;
            double xe;
            double xq;
            double xw;
            double[] xx = new double[8];
            double yn;
            double yq;
            double ys;
            double[] yy = new double[8];

            e2 = 0.0;
            //
            //  Compute the matrix entries by integrating over each element.
            //
            ex_num = (nx - 1) / 2;
            ey_num = (ny - 1) / 2;

            for (ey = 0; ey < ey_num; ey++)
            {
                s = 2 * ey;
                mm = 2 * ey + 1;
                n = 2 * ey + 2;

                ys = y[s];
                yn = y[n];

                yy[0] = y[n];
                yy[1] = y[n];
                yy[2] = y[n];
                yy[3] = y[mm];
                yy[4] = y[s];
                yy[5] = y[s];
                yy[6] = y[s];
                yy[7] = y[mm];

                for (ex = 0; ex < ex_num; ex++)
                {
                    w = 2 * ex;
                    cc = 2 * ex + 1;
                    e = 2 * ex + 2;

                    xe = x[e];
                    xw = x[w];

                    xx[0] = x[e];
                    xx[1] = x[cc];
                    xx[2] = x[w];
                    xx[3] = x[w];
                    xx[4] = x[w];
                    xx[5] = x[cc];
                    xx[6] = x[e];
                    xx[7] = x[e];
                    //
                    //  Node indices
                    //
                    //  3  2  1  wn  cn  en
                    //  4     8  wm      em
                    //  5  6  7  ws  cs  es
                    //
                    node[0] = (3 * ey + 3) * ey_num + 2 * ey + 2 * ex + 4;
                    node[1] = (3 * ey + 3) * ey_num + 2 * ey + 2 * ex + 3;
                    node[2] = (3 * ey + 3) * ey_num + 2 * ey + 2 * ex + 2;
                    node[3] = (3 * ey + 2) * ey_num + 2 * ey + ex + 1;
                    node[4] = (3 * ey) * ey_num + 2 * ey + 2 * ex;
                    node[5] = (3 * ey) * ey_num + 2 * ey + 2 * ex + 1;
                    node[6] = (3 * ey) * ey_num + 2 * ey + 2 * ex + 2;
                    node[7] = (3 * ey + 2) * ey_num + 2 * ey + ex + 2;

                    for (qx = 0; qx < quad_num; qx++)
                    {
                        xq = ((1.0 - abscissa[qx]) * x[e]
                              + (1.0 + abscissa[qx]) * x[w])
                             / 2.0;

                        for (qy = 0; qy < quad_num; qy++)
                        {
                            yq = ((1.0 - abscissa[qy]) * y[n]
                                  + (1.0 + abscissa[qy]) * y[s])
                                 / 2.0;

                            wq = weight[qx] * (x[e] - x[w]) / 2.0
                                * weight[qy] * (y[n] - y[s]) / 2.0;

                            v = basis_serene(xq, yq, xw, ys, xe, yn, xx, yy);

                            uq = 0.0;
                            for (k = 0; k < 8; k++)
                            {
                                uq = uq + u[node[k]] * v[k];
                            }

                            eq = exact(xq, yq);
                            e2 = e2 + wq * Math.Pow(uq - eq, 2);
                        }
                    }
                }
            }

            e2 = Math.Sqrt(e2);

            return e2;
        }

        public static double not1(double x1, double x2, double x3)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NOT1 evaluates a factor for serendipity basis functions.
            //
            //  Discussion:
            //
            //    not1(x1,x2,x3) evaluates at the point x1, the basis factor that
            //    is 0 at x2 and 1 at x3:
            //
            //      not1(x1,x2,x3) = ( x1 - x2 ) / ( x3 - x2 )  
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 July 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X1, the evaluation point.
            //
            //    Input, double X2, X3, values that define the factor.
            //
            //    Output, double NOT1, the value of the basis function factor.
            //
        {
            double value;

            value = (x1 - x2) / (x3 - x2);

            return value;
        }

        public static double not1d(double x2, double x3)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NOT1D differentiates a factor for serendipity basis functions.
            //
            //  Discussion:
            //
            //    not1(x1,x2,x3) evaluates at the point x1, the basis factor that
            //    is 0 at x2 and 1 at x3:
            //
            //      not1(x1,x2,x3) = ( x1 - x2 ) / ( x3 - x2 )  
            //
            //    This function returns the derivative of the factor with respect to x1:
            //
            //      not1d(x1,x2,x3) = 1 / ( x3 - x2 )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 July 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X2, X3, values that define the factor.
            //
            //    Output, double NOT1D, the derivative of the basis function 
            //    factor.
            //
        {
            double value;

            value = 1.0 / (x3 - x2);

            return value;
        }

        public static double not2(double x1, double y1, double x2, double y2, double x3, double y3,
                double x4, double y4)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NOT2 evaluates a factor for serendipity basis functions.
            //
            //  Discussion:
            //
            //    not2(x1,y1,x2,y2,x3,y3,x4,y4) evaluates at the point (x1,y1), the basis 
            //    factor that is 0 at (x2,y2) and (x3,y3) and 1 at (x4,y4):
            //
            //          ( ( x1 - x2 ) * ( y3 - y2 ) - ( x3 - x2 ) * ( y1 - y2 ) )
            //        / ( ( x4 - x2 ) * ( y3 - y2 ) - ( x3 - x2 ) * ( y4 - y2 ) )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 July 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X1, Y2, the evaluation point.
            //
            //    Input, double X2, Y2, X3, Y3, values that define the factor.
            //
            //    Output, double NOT2, the value of the basis function factor.
            //
        {
            double value;

            value = ((x1 - x2) * (y3 - y2) - (x3 - x2) * (y1 - y2))
                    / ((x4 - x2) * (y3 - y2) - (x3 - x2) * (y4 - y2));

            return value;
        }

        public static double not2dx(double x2, double y2, double x3, double y3, double x4, double y4)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NOT2DX evaluates a factor for serendipity basis functions.
            //
            //  Discussion:
            //
            //    not2(x1,y1,x2,y2,x3,y3,x4,y4) evaluates at the point (x1,y1), the basis 
            //    factor that is 0 at (x2,y2) and (x3,y3) and 1 at (x4,y4):
            //
            //          ( ( x1 - x2 ) * ( y3 - y2 ) - ( x3 - x2 ) * ( y1 - y2 ) )
            //        / ( ( x4 - x2 ) * ( y3 - y2 ) - ( x3 - x2 ) * ( y4 - y2 ) )
            //
            //    not2dx returns the derivative of this function with respect to X1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 July 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X2, Y2, X3, Y3, values that define the factor.
            //
            //    Output, double NOT2DX, the derivative of the basis function 
            //    factor with respect to X1.
            //
        {
            double value;

            value = (1.0 * (y3 - y2) + 0.0)
                    / ((x4 - x2) * (y3 - y2) - (x3 - x2) * (y4 - y2));

            return value;
        }

        public static double not2dy(double x2, double y2, double x3, double y3, double x4, double y4)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NOT2DY evaluates a factor for serendipity basis functions.
            //
            //  Discussion:
            //
            //    not2(x1,y1,x2,y2,x3,y3,x4,y4) evaluates at the point (x1,y1), the basis 
            //    factor that is 0 at (x2,y2) and (x3,y3) and 1 at (x4,y4):
            //
            //          ( ( x1 - x2 ) * ( y3 - y2 ) - ( x3 - x2 ) * ( y1 - y2 ) )
            //        / ( ( x4 - x2 ) * ( y3 - y2 ) - ( x3 - x2 ) * ( y4 - y2 ) )
            //
            //    not2dy returns the derivatives of this function with respect to Y1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 July 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X2, Y2, X3, Y3, values that define the factor.
            //
            //    Output, double NOT2DY, the derivative of the basis function 
            //    factor with respect to Y1.
            //
        {
            double value;

            value = (0.0 - (x3 - x2) * 1.0)
                    / ((x4 - x2) * (y3 - y2) - (x3 - x2) * (y4 - y2));

            return value;
        }



    }
}