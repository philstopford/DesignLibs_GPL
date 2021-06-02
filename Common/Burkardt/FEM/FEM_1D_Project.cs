using System;
using Burkardt.Types;

namespace Burkardt.FEM
{
    public static class FEM_1D_Project
    {

        public static double[] fem1d_approximate(int sample_node_num, int sample_value_dim,
            double[] sample_node_x, double[] sample_value, int fem_node_num,
        double[] fem_node_x, int fem_element_order, int fem_element_num,
        int fem_value_dim, int fem_value_num )
//****************************************************************************80
//
//  Purpose:
//
//    FEM1D_APPROXIMATE approximates data at sample points with an FEM function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 May 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int SAMPLE_NODE_NUM, the number of sample points.
//
//    Input, int SAMPLE_VALUE_DIM, the value dimension.
//
//    Input, double SAMPLE_NODE_X[SAMPLE_NODE_NUM], the sample nodes.
//
//    Input, double SAMPLE_VALUE[VALUE_DIM*SAMPLE_NODE_NUM],
//    the values at sample nodes.
//
//    Input, int FEM_NODE_NUM, the number of FEM nodes.
//
//    Input, double FEM_NODE_X[FEM_NODE_NUM], the FEM nodes.  
//
//    Input, int FEM_ELEMENT_ORDER, the element order.
//
//    Input, int FEM_ELEMENT_NUM, the number of elements.
//
//    Input, int FEM_VALUE_DIM, the FEM value dimension.
//
//    Input, int FEM_VALUE_NUM, the number of FEM values.
//
//    Output, double FEM1D_APPROXIMATE[FEM_VALUE_DIM*FEM_VALUE_NUM], 
//    the FEM values.
//
        {
            int QUAD_NUM = 2;

            int phi_num = 0;
            double[] phi_v = new double[3];
            double[] phi_x = new double[3];
            int quad_num = QUAD_NUM;
            double[] quad_x =  {
                -0.577350269189625764509148780502,
                0.577350269189625764509148780502
            }
            ;
            double[] quad_w =  {
                1.0, 1.0
            }
            ;
            //
//  Set up the matrix A.
//
            double[] a = typeMethods.r8mat_zero_new(3, fem_node_num);

            for (int l = 0; l < fem_node_num - 1; l++)
            {
                int r = l + 1;
                double xl = fem_node_x[l];
                double xr = fem_node_x[r];

                for (int quad = 0; quad < quad_num; quad++)
                {
                    double xq = ((1.0 - quad_x[quad]) * xl
                                 + (1.0 + quad_x[quad]) * xr)
                                / 2.0;

                    double wq = quad_w[quad] * (xr - xl) / 2.0;

                    double phil = (xq - xr)
                                  / (xl - xr);

                    double phir = (xl - xq)
                                  / (xl - xr);

                    a[1 + l * 3] = a[1 + l * 3] + wq * phil * phil;
                    a[2 + l * 3] = a[2 + l * 3] + wq * phil * phir;

                    a[0 + r * 3] = a[0 + r * 3] + wq * phir * phil;
                    a[1 + r * 3] = a[1 + r * 3] + wq * phir * phir;
                }
            }

            typeMethods.r83_np_fa(fem_node_num, a);
//
//  Set up the right hand side b.
//
            double[] b = new double[fem_node_num];
            double[] v = new double[sample_node_num];
            double[] fem_value = new double[fem_value_dim * fem_value_num];

            for (int dim = 0; dim < fem_value_dim; dim++)
            {
                for (int i = 0; i < fem_node_num; i++)
                {
                    if (i == 0)
                    {
                        phi_num = 2;
                        phi_x[0] = fem_node_x[0];
                        phi_x[1] = fem_node_x[1];
                        phi_v[0] = 1.0;
                        phi_v[1] = 0.0;
                    }
                    else if (i < fem_node_num - 1)
                    {
                        phi_num = 3;
                        phi_x[0] = fem_node_x[i - 1];
                        phi_x[1] = fem_node_x[i];
                        phi_x[2] = fem_node_x[i + 1];
                        phi_v[0] = 0.0;
                        phi_v[1] = 1.0;
                        phi_v[2] = 0.0;
                    }
                    else if (i == fem_node_num - 1)
                    {
                        phi_num = 2;
                        phi_x[0] = fem_node_x[fem_node_num - 2];
                        phi_x[1] = fem_node_x[fem_node_num - 1];
                        phi_v[0] = 0.0;
                        phi_v[1] = 1.0;
                    }

                    double a1 = phi_x[0];
                    double b1 = phi_x[phi_num - 1];

                    for (int j = 0; j < sample_node_num; j++)
                    {
                        v[j] = sample_value[dim + j * sample_value_dim];
                    }

                    double integral = piecewise_linear_product_quad(a1, b1, phi_num, phi_x,
                        phi_v, sample_node_num, sample_node_x, v);

                    b[i] = integral;
                }

                int job = 0;
                double[] x = typeMethods.r83_np_sl(fem_node_num, a, b, job);

                for (int i = 0; i < fem_node_num; i++)
                {
                    fem_value[dim + i * fem_value_dim] = x[i];
                }
            }

            return fem_value;
        }

        public static double piecewise_linear_product_quad(double a, double b, int f_num,
            double[] f_x, double[] f_v, int g_num, double[] g_x, double[] g_v )
//****************************************************************************80
//
//  Purpose:
//
//    PIECEWISE_LINEAR_PRODUCT_QUAD: piecewise linear product integral.
//
//  Discussion:
//
//    We are given two piecewise linear functions F(X) and G(X) and we wish
//    to compute the exact value of the integral
//
//      INTEGRAL = Integral ( A <= X <= B ) F(X) * G(X) dx
//
//    The functions F(X) and G(X) are defined as tables of coordinates X and
//    values V.  A piecewise linear function is evaluated at a point X by 
//    evaluating the interpolant to the data at the endpoints of the interval 
//    containing X.  
//
//    It must be the case that A <= B.
//
//    It must be the case that the node coordinates F_X(*) and G_X(*) are
//    given in ascending order.
//
//    It must be the case that:
//
//      F_X(1) <= A and B <= F_X(F_NUM)
//      G_X(1) <= A and B <= G_X(G_NUM)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the limits of integration.
//
//    Input, int F_NUM, the number of nodes for F.
//
//    Input, double F_X[F_NUM], the node coordinates for F.
//
//    Input, double F_V[F_NUM], the nodal values for F.
//
//    Input, int G_NUM, the number of nodes for G.
//
//    Input, double G_X[G_NUM], the node coordinates for G.
//
//    Input, double G_V[G_NUM], the nodal values for G.
//
//    Output, double INTEGRAL, the integral of F(X) * G(X)
//    from A to B.
//
        {
            double integral = 0.0;

            if (f_x[f_num - 1] <= a || g_x[g_num - 1] <= a)
            {
                return integral;
            }

            if (f_num < 2 || g_num < 2)
            {
                return integral;
            }

            double xr = a;

            int f_left = 0;
            typeMethods.r8vec_bracket3(f_num, f_x, xr, ref f_left);
            double fr = f_v[f_left] + (xr - f_x[f_left]) * (f_v[f_left + 1] - f_v[f_left])
                / (f_x[f_left + 1] - f_x[f_left]);

            int g_left = 0;
            typeMethods.r8vec_bracket3(g_num, g_x, xr, ref g_left);
            double gr = g_v[g_left] + (xr - g_x[g_left]) * (g_v[g_left + 1] - g_v[g_left])
                / (g_x[g_left + 1] - g_x[g_left]);

            double xr_max = b;
            xr_max = Math.Min(xr_max, f_x[f_num - 1]);
            xr_max = Math.Min(xr_max, g_x[g_num - 1]);

            while (xr < xr_max)
            {
//
//  Shift right values to left.
//
                double xl = xr;
                double fl = fr;
                double gl = gr;
//
//  Determine the new right values.
//  The hard part is figuring out how to advance XR some, but not too much.
//
                xr = xr_max;

                for (int i = 1; i <= 2; i++)
                {
                    if (f_left + i <= f_num - 1)
                    {
                        if (xl < f_x[f_left + i] && f_x[f_left + i] < xr)
                        {
                            xr = f_x[f_left + i];
                            break;
                        }
                    }
                }

                for (int i = 1; i <= 2; i++)
                {
                    if (g_left + i <= g_num - 1)
                    {
                        if (xl < g_x[g_left + i] && g_x[g_left + i] < xr)
                        {
                            xr = g_x[g_left + i];
                            break;
                        }
                    }
                }

                typeMethods.r8vec_bracket3(f_num, f_x, xr, ref f_left);
                fr = f_v[f_left] + (xr - f_x[f_left]) * (f_v[f_left + 1] - f_v[f_left])
                    / (f_x[f_left + 1] - f_x[f_left]);

                typeMethods.r8vec_bracket3(g_num, g_x, xr, ref g_left);
                gr = g_v[g_left] + (xr - g_x[g_left]) * (g_v[g_left + 1] - g_v[g_left])
                    / (g_x[g_left + 1] - g_x[g_left]);
//
//  Form the linear polynomials for F(X) and G(X) over [XL,XR],
//  then the product H(X), integrate H(X) and add to the running total.
//
                if (double.Epsilon <= Math.Abs(xr - xl))
                {
                    double f1 = fl - fr;
                    double f0 = fr * xl - fl * xr;

                    double g1 = gl - gr;
                    double g0 = gr * xl - gl * xr;

                    double h2 = f1 * g1;
                    double h1 = f1 * g0 + f0 * g1;
                    double h0 = f0 * g0;

                    h2 = h2 / 3.0;
                    h1 = h1 / 2.0;

                    double bit = ((h2 * xr + h1) * xr + h0) * xr
                                 - ((h2 * xl + h1) * xl + h0) * xl;

                    integral = integral + bit / (xr - xl) / (xr - xl);
                }
            }

            return integral;
        }
    }
}