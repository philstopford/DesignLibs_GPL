using System;
using Burkardt.IntegralNS;
using Burkardt.Quadrature;
using Burkardt.Types;
using Monomial = Burkardt.MonomialNS.Monomial;

namespace StroudTest
{
    using StroudIntegral = Burkardt.IntegralNS.Stroud;
    using StroudQuadrature = Burkardt.Quadrature.Stroud;

    public static class cn
    {
        public static void cn_geg_test(int n, double alpha, int[] expon)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CN_GEG_TEST tests the rules for CN with Gegenbauer weight on a monomial.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 March 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double c1;
            int d;
            double delta0;
            double err;
            double exact;
            double gamma0;
            int i;
            int o;
            int p;
            double pi = 3.141592653589793;
            double quad;
            double[] v;
            double volume_1d;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("  N = " + n + "");
            Console.WriteLine("  ALPHA = " + alpha + "");
            string cout = "  EXPON = ";
            for (i = 0; i < n; i++)
            {
                cout += expon[i].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
            d = typeMethods.i4vec_sum(n, expon);
            Console.WriteLine("  Degree = " + d + "");
            Console.WriteLine("");

            exact = CN_Geg.cn_geg_monomial_integral(n, alpha, expon);

            p = 0;

            if (d <= p)
            {
                o = StroudQuadrature.cn_geg_00_1_size(n, alpha);
                x = new double[n * o];
                w = new double[o];
                StroudQuadrature.cn_geg_00_1(n, alpha, o, ref x, ref w);
                v = Monomial.monomial_value(n, o, x, expon);
                quad = typeMethods.r8vec_dot_product(o, w, v);
                err = Math.Abs(quad - exact);
                Console.WriteLine("  CN_GEG_00_1:   "
                                  + "  " + o.ToString().PadLeft(6)
                                  + "  " + quad.ToString().PadLeft(14)
                                  + "  " + err.ToString().PadLeft(14) + "");
            }

            p = 1;

            if (d <= p)
            {
                o = StroudQuadrature.cn_geg_01_1_size(n, alpha);
                x = new double[n * o];
                w = new double[o];
                StroudQuadrature.cn_geg_01_1(n, alpha, o, ref x, ref w);
                v = Monomial.monomial_value(n, o, x, expon);
                quad = typeMethods.r8vec_dot_product(o, w, v);
                err = Math.Abs(quad - exact);
                Console.WriteLine("  CN_GEG_01_1:   "
                                  + "  " + o.ToString().PadLeft(6)
                                  + "  " + quad.ToString().PadLeft(14)
                                  + "  " + err.ToString().PadLeft(14) + "");
            }

            p = 2;

            if (d <= p)
            {
                o = Xiu.cn_geg_02_xiu_size(n, alpha);
                x = new double[n * o];
                w = new double[o];
                Xiu.cn_geg_02_xiu(n, alpha, o, ref x, ref w);
                v = Monomial.monomial_value(n, o, x, expon);
                quad = typeMethods.r8vec_dot_product(o, w, v);
                err = Math.Abs(quad - exact);
                Console.WriteLine("  CN_GEG_02_XIU: "
                                  + "  " + o.ToString().PadLeft(6)
                                  + "  " + quad.ToString().PadLeft(14)
                                  + "  " + err.ToString().PadLeft(14) + "");

                o = GolubWelsch_Xiu.gw_02_xiu_size(n);
                gamma0 = 1.0;
                delta0 = 0.0;
                c1 = 1.0 / (2.0 * alpha + 3.0);
                volume_1d = Math.Sqrt(pi) * typeMethods.r8_gamma(alpha + 1.0)
                            / typeMethods.r8_gamma(alpha + 1.5);
                x = new double[n * o];
                w = new double[o];
                GolubWelsch_Xiu.gw_02_xiu(n, o, gamma0, delta0, c1, volume_1d, ref x, ref w);
                v = Monomial.monomial_value(n, o, x, expon);
                quad = typeMethods.r8vec_dot_product(o, w, v);
                err = Math.Abs(quad - exact);
                Console.WriteLine("  GW_02_XIU:     "
                                  + "  " + o.ToString().PadLeft(6)
                                  + "  " + quad.ToString().PadLeft(14)
                                  + "  " + err.ToString().PadLeft(14) + "");
            }

            p = 3;

            if (d <= p)
            {
                o = Xiu.cn_geg_03_xiu_size(n, alpha);
                x = new double[n * o];
                w = new double[o];
                Xiu.cn_geg_03_xiu(n, alpha, o, ref x, ref w);
                v = Monomial.monomial_value(n, o, x, expon);
                quad = typeMethods.r8vec_dot_product(o, w, v);
                err = Math.Abs(quad - exact);
                Console.WriteLine("  CN_GEG_03_XIU: "
                                  + "  " + o.ToString().PadLeft(6)
                                  + "  " + quad.ToString().PadLeft(14)
                                  + "  " + err.ToString().PadLeft(14) + "");
            }

            Console.WriteLine("  EXACT                  "
                              + "  " + exact.ToString().PadLeft(14) + "");

        }

        public static void cn_jac_test(int n, double alpha, double beta, int[] expon)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CN_JAC_TEST tests the rules for CN with Jacobi weight on a monomial.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 March 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double c1;
            int d;
            double delta0;
            double err;
            double exact;
            double gamma0;
            int i;
            int o;
            int p;
            double quad;
            double[] v;
            double volume_1d;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("  N = " + n + "");
            Console.WriteLine("  ALPHA = " + alpha + "");
            Console.WriteLine("  BETA =  " + beta + "");
            string cout = "  EXPON = ";
            for (i = 0; i < n; i++)
            {
                cout += expon[i].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
            d = typeMethods.i4vec_sum(n, expon);
            Console.WriteLine("  Degree = " + d + "");
            Console.WriteLine("");

            exact = CN_Geg.cn_jac_monomial_integral(n, alpha, beta, expon);

            p = 0;

            if (d <= p)
            {
                o = CN_Jac.cn_jac_00_1_size(n, alpha, beta);
                x = new double[n * o];
                w = new double[o];
                CN_Jac.cn_jac_00_1(n, alpha, beta, o, ref x, ref w);
                v = Monomial.monomial_value(n, o, x, expon);
                quad = typeMethods.r8vec_dot_product(o, w, v);
                err = Math.Abs(quad - exact);
                Console.WriteLine("  CN_JAC_00_1:   "
                                  + "  " + o.ToString().PadLeft(6)
                                  + "  " + quad.ToString().PadLeft(14)
                                  + "  " + err.ToString().PadLeft(14) + "");
            }

            p = 1;

            if (d <= p)
            {
                o = CN_Jac.cn_jac_01_1_size(n, alpha, beta);
                x = new double[n * o];
                w = new double[o];
                CN_Jac.cn_jac_01_1(n, alpha, beta, o, ref x, ref w);
                v = Monomial.monomial_value(n, o, x, expon);
                quad = typeMethods.r8vec_dot_product(o, w, v);
                err = Math.Abs(quad - exact);
                Console.WriteLine("  CN_JAC_01_1:   "
                                  + "  " + o.ToString().PadLeft(6)
                                  + "  " + quad.ToString().PadLeft(14)
                                  + "  " + err.ToString().PadLeft(14) + "");
            }

            p = 2;

            if (d <= p)
            {
                o = CN_Jac_Xiu.cn_jac_02_xiu_size(n, alpha, beta);
                x = new double[n * o];
                w = new double[o];
                CN_Jac_Xiu.cn_jac_02_xiu(n, alpha, beta, o, ref x, ref w);
                v = Monomial.monomial_value(n, o, x, expon);
                quad = typeMethods.r8vec_dot_product(o, w, v);
                err = Math.Abs(quad - exact);
                Console.WriteLine("  CN_JAC_02_XIU: "
                                  + "  " + o.ToString().PadLeft(6)
                                  + "  " + quad.ToString().PadLeft(14)
                                  + "  " + err.ToString().PadLeft(14) + "");

                o = GolubWelsch_Xiu.gw_02_xiu_size(n);
                gamma0 = (alpha + beta + 2.0) / 2.0;
                delta0 = (alpha - beta) / 2.0;
                c1 = 2.0 * (alpha + 1.0) * (beta + 1.0) / (alpha + beta + 3.0)
                                                        / (alpha + beta + 2.0);
                volume_1d = Math.Pow(2.0, alpha + beta + 1.0) * typeMethods.r8_gamma(alpha + 1.0)
                                                         * typeMethods.r8_gamma(beta + 1.0) / (alpha + beta + 1.0) /
                            typeMethods.r8_gamma(alpha + beta + 1.0);
                x = new double[n * o];
                w = new double[o];
                GolubWelsch_Xiu.gw_02_xiu(n, o, gamma0, delta0, c1, volume_1d, ref x, ref w);
                v = Monomial.monomial_value(n, o, x, expon);
                quad = typeMethods.r8vec_dot_product(o, w, v);
                err = Math.Abs(quad - exact);
                Console.WriteLine("  GW_02_XIU:     "
                                       + "  " + o.ToString().PadLeft(6)
                                       + "  " + quad.ToString().PadLeft(14)
                                       + "  " + err.ToString().PadLeft(14) + "");
            }

            Console.WriteLine("  EXACT                  "
                              + "  " + exact.ToString().PadLeft(14) + "");

        }

        public static void cn_leg_test(int n, int[] expon)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CN_LEG_TEST tests the rules for CN with Legendre weight on a monomial.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 March 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double c1;
            int d;
            double delta0;
            double err;
            double exact;
            double gamma0;
            int i;
            int o;
            int option;
            int p;
            double quad;
            double[] v;
            double volume_1d;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("  N = " + n + "");
            string cout = "  EXPON = ";
            for (i = 0; i < n; i++)
            {
                cout += expon[i].ToString().PadLeft(4);
            }

            Console.WriteLine("");
            d = typeMethods.i4vec_sum(n, expon);
            Console.WriteLine("  Degree = " + d + "");
            Console.WriteLine("");

            exact = CN_Geg.cn_leg_monomial_integral(n, expon);

            p = 1;

            if (d <= p)
            {
                o = CN_Leg.cn_leg_01_1_size(n);
                x = new double[n * o];
                w = new double[o];
                CN_Leg.cn_leg_01_1(n, o, ref x, ref w);
                v = Monomial.monomial_value(n, o, x, expon);
                quad = typeMethods.r8vec_dot_product(o, w, v);
                err = Math.Abs(quad - exact);
                Console.WriteLine("  CN_LEG_01_1:   "
                                  + "  " + o.ToString().PadLeft(6)
                                  + "  " + quad.ToString().PadLeft(14)
                                  + "  " + err.ToString().PadLeft(14) + "");
            }

            p = 2;

            if (d <= p)
            {
                o = Xiu.cn_leg_02_xiu_size(n);
                x = new double[n * o];
                w = new double[o];
                Xiu.cn_leg_02_xiu(n, o, ref x, ref w);
                v = Monomial.monomial_value(n, o, x, expon);
                quad = typeMethods.r8vec_dot_product(o, w, v);
                err = Math.Abs(quad - exact);
                Console.WriteLine("  CN_LEG_02_XIU: "
                                  + "  " + o.ToString().PadLeft(6)
                                  + "  " + quad.ToString().PadLeft(14)
                                  + "  " + err.ToString().PadLeft(14) + "");

                o = GolubWelsch_Xiu.gw_02_xiu_size(n);
                gamma0 = 1.0;
                delta0 = 0.0;
                c1 = 1.0 / 3.0;
                volume_1d = 2.0;
                x = new double[n * o];
                w = new double[o];
                GolubWelsch_Xiu.gw_02_xiu(n, o, gamma0, delta0, c1, volume_1d, ref x, ref w);
                v = Monomial.monomial_value(n, o, x, expon);
                quad = typeMethods.r8vec_dot_product(o, w, v);
                err = Math.Abs(quad - exact);
                Console.WriteLine("  GW_02_XIU:     "
                                  + "  " + o.ToString().PadLeft(6)
                                  + "  " + quad.ToString().PadLeft(14)
                                  + "  " + err.ToString().PadLeft(14) + "");
            }

            p = 3;

            if (d <= p)
            {
                o = StroudQuadrature.cn_leg_03_1_size(n);
                x = new double[n * o];
                w = new double[o];
                StroudQuadrature.cn_leg_03_1(n, o, ref x, ref w);
                v = Monomial.monomial_value(n, o, x, expon);
                quad = typeMethods.r8vec_dot_product(o, w, v);
                err = Math.Abs(quad - exact);
                Console.WriteLine("  CN_LEG_03_1:   "
                                  + "  " + o.ToString().PadLeft(6)
                                  + "  " + quad.ToString().PadLeft(14)
                                  + "  " + err.ToString().PadLeft(14) + "");

                o = Xiu.cn_leg_03_xiu_size(n);
                x = new double[n * o];
                w = new double[o];
                Xiu.cn_leg_03_xiu(n, o, ref x, ref w); 
                v = Monomial.monomial_value(n, o, x, expon);
                quad = typeMethods.r8vec_dot_product(o, w, v);
                err = Math.Abs(quad - exact);
                Console.WriteLine("  CN_LEG_03_XIU: "
                                  + "  " + o.ToString().PadLeft(6)
                                  + "  " + quad.ToString().PadLeft(14)
                                  + "  " + err.ToString().PadLeft(14) + "");
            }

            p = 5;

            if (d <= p)
            {
                if (4 <= n && n <= 6)
                {
                    o = StroudQuadrature.cn_leg_05_1_size(n);
                    x = new double[n * o];
                    w = new double[o];
                    option = 1;
                    StroudQuadrature.cn_leg_05_1(n, option, o, ref x, ref w);
                    v = Monomial.monomial_value(n, o, x, expon);
                    quad = typeMethods.r8vec_dot_product(o, w, v);
                    err = Math.Abs(quad - exact);
                    Console.WriteLine("  CN_LEG_05_1(1):"
                                      + "  " + o.ToString().PadLeft(6)
                                      + "  " + quad.ToString().PadLeft(14)
                                      + "  " + err.ToString().PadLeft(14) + "");
                }

                if (4 <= n && n <= 5)
                {
                    o = StroudQuadrature.cn_leg_05_1_size(n);
                    x = new double[n * o];
                    w = new double[o];
                    option = 2;
                    StroudQuadrature.cn_leg_05_1(n, option, o, ref x, ref w);
                    v = Monomial.monomial_value(n, o, x, expon);
                    quad = typeMethods.r8vec_dot_product(o, w, v);
                    err = Math.Abs(quad - exact);
                    Console.WriteLine("  CN_LEG_05_1(2):"
                                      + "  " + o.ToString().PadLeft(6)
                                      + "  " + quad.ToString().PadLeft(14)
                                      + "  " + err.ToString().PadLeft(14) + "");
                }

                if (2 <= n)
                {
                    o = StroudQuadrature.cn_leg_05_2_size(n);
                    x = new double[n * o];
                    w = new double[o];
                    StroudQuadrature.cn_leg_05_2(n, o, ref x, ref w);
                    v = Monomial.monomial_value(n, o, x, expon);
                    quad = typeMethods.r8vec_dot_product(o, w, v);
                    err = Math.Abs(quad - exact);
                    Console.WriteLine("  CN_LEG_05_2:   "
                                           + "  " + o.ToString().PadLeft(6)
                                           + "  " + quad.ToString().PadLeft(14)
                                           + "  " + err.ToString().PadLeft(14) + "");
                }
            }

            Console.WriteLine("  EXACT                  "
                                   + "  " + exact.ToString().PadLeft(14) + "");

        }

    }
}