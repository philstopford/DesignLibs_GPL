using System;
using Burkardt.Types;

namespace Burkardt.PiecewiseLinear
{
    public static class ProductIntegral
    {
        public static double pwl_product_integral(double a, double b, int f_num,
                double[] f_x, double[] f_v, int g_num, double[] g_x, double[] g_v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    pwl_PRODUCT_INTEGRAL: piecewise linear product integral.
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
            double bit;
            int f_left;
            double f0;
            double f1;
            double fl;
            double fr;
            int g_left;
            double g0;
            double g1;
            double gl;
            double gr;
            double h0;
            double h1;
            double h2;
            int i;
            double integral;
            double xl;
            double xr;
            double xr_max;

            integral = 0.0;

            if (f_x[f_num - 1] <= a || g_x[g_num - 1] <= a)
            {
                return integral;
            }

            if (f_num < 2 || g_num < 2)
            {
                return integral;
            }

            xr = a;

            f_left = 0;
            typeMethods.r8vec_bracket3(f_num, f_x, xr, ref f_left);
            fr = f_v[f_left] + (xr - f_x[f_left]) * (f_v[f_left + 1] - f_v[f_left])
                / (f_x[f_left + 1] - f_x[f_left]);

            g_left = 0;
            typeMethods.r8vec_bracket3(g_num, g_x, xr, ref g_left);
            gr = g_v[g_left] + (xr - g_x[g_left]) * (g_v[g_left + 1] - g_v[g_left])
                / (g_x[g_left + 1] - g_x[g_left]);

            xr_max = b;
            xr_max = Math.Min(xr_max, f_x[f_num - 1]);
            xr_max = Math.Min(xr_max, g_x[g_num - 1]);

            while (xr < xr_max)
            {
                //
                //  Shift right values to left.
                //
                xl = xr;
                fl = fr;
                gl = gr;
                //
                //  Determine the new right values.
                //  The hard part is figuring out how to advance XR some, but not too much.
                //
                xr = xr_max;

                for (i = 1; i <= 2; i++)
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

                for (i = 1; i <= 2; i++)
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
                    f1 = fl - fr;
                    f0 = fr * xl - fl * xr;

                    g1 = gl - gr;
                    g0 = gr * xl - gl * xr;

                    h2 = f1 * g1;
                    h1 = f1 * g0 + f0 * g1;
                    h0 = f0 * g0;

                    h2 = h2 / 3.0;
                    h1 = h1 / 2.0;

                    bit = ((h2 * xr + h1) * xr + h0) * xr
                          - ((h2 * xl + h1) * xl + h0) * xl;

                    integral = integral + bit / (xr - xl) / (xr - xl);
                }
            }

            return integral;
        }

        public static double pwl_product_quad(double a, double b, int f_num,
                double[] f_x, double[] f_v, int g_num, double[] g_x, double[] g_v,
                int quad_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    pwl_PRODUCT_QUAD: estimate piecewise linear product integral.
            //
            //  Discussion:
            //
            //    We are given two piecewise linear functions F(X) and G(X) and we wish
            //    to estimate the value of the integral
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
            //    29 April 2009
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
            //    Input, int QUAD_NUM, the number of quadrature points.
            //
            //    Output, double pwl_PRODUCT_QUAD, an estimate for the integral 
            //    of F(X) * G(X) from A to B.
            //
        {
            double a2;
            double b2;
            int f_left;
            double fq;
            int g_left;
            double gq;
            int i;
            double quad;
            double xq;

            quad = 0.0;

            f_left = 0;
            g_left = 0;

            a2 = a;
            a2 = Math.Max(a2, f_x[0]);
            a2 = Math.Max(a2, g_x[0]);

            b2 = b;
            b2 = Math.Min(b2, f_x[f_num - 1]);
            b2 = Math.Min(b2, g_x[g_num - 1]);

            for (i = 1; i <= quad_num; i++)
            {
                xq = ((double)(2 * i - 1) * b2
                      + (double)(2 * quad_num - 2 * i + 1) * a2)
                     / (double)(2 * quad_num);

                typeMethods.r8vec_bracket3(f_num, f_x, xq, ref f_left);

                fq = f_v[f_left] + (xq - f_x[f_left]) * (f_v[f_left + 1] - f_v[f_left])
                    / (f_x[f_left + 1] - f_x[f_left]);

                typeMethods.r8vec_bracket3(g_num, g_x, xq, ref g_left);

                gq = g_v[g_left] + (xq - g_x[g_left]) * (g_v[g_left + 1] - g_v[g_left])
                    / (g_x[g_left + 1] - g_x[g_left]);

                quad = quad + fq * gq;
            }

            quad = quad * (b - a) / (double)(quad_num);

            return quad;
        }
    }
}