using System;
using Burkardt.MonomialNS;
using Burkardt.Quadrature;
using Burkardt.Types;

namespace StroudTest
{
    using EpnLagIntegral = Burkardt.IntegralNS.Epn_Lag;
    using EpnLagQuadrature = Burkardt.Quadrature.Epn_Lag;

    public static class epn
    {
        public static void epn_glg_test(int n, int[] expon, double alpha)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EPN_GLG_TEST tests the rules for EPN with GLG weight on a monomial.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    27 January 2010
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
            string cout = "  EXPON = ";
            for (i = 0; i < n; i++)
            {
                cout += expon[i].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
            d = typeMethods.i4vec_sum(n, expon);
            Console.WriteLine("  Degree = " + d + "");
            Console.WriteLine("");

            exact = EpnLagIntegral.epn_glg_monomial_integral(n, expon, alpha);

            p = 0;

            if (d <= p)
            {
                o = EpnLagQuadrature.epn_glg_00_1_size(n, alpha);
                x = new double[n * o];
                w = new double[o];
                EpnLagQuadrature.epn_glg_00_1(n, alpha, o, ref x, ref w);
                v = Monomial.monomial_value(n, o, x, expon);
                quad = typeMethods.r8vec_dot_product(o, w, v);
                err = Math.Abs(quad - exact);
                Console.WriteLine("  EPN_GLG_00_1:   "
                                  + "  " + o.ToString().PadLeft(6)
                                  + "  " + quad.ToString().PadLeft(14)
                                  + "  " + err.ToString().PadLeft(14) + "");
            }

            p = 1;

            if (d <= p)
            {
                o = EpnLagQuadrature.epn_glg_01_1_size(n, alpha);
                x = new double[n * o];
                w = new double[o];
                EpnLagQuadrature.epn_glg_01_1(n, alpha, o, ref x, ref w);
                v = Monomial.monomial_value(n, o, x, expon);
                quad = typeMethods.r8vec_dot_product(o, w, v);
                err = Math.Abs(quad - exact);
                Console.WriteLine("  EPN_GLG_01_1:   "
                                  + "  " + o.ToString().PadLeft(6)
                                  + "  " + quad.ToString().PadLeft(14)
                                  + "  " + err.ToString().PadLeft(14) + "");
            }

            p = 2;

            if (d <= p)
            {
                o = EpnLagQuadrature.epn_glg_02_xiu_size(n, alpha);
                x = new double[n * o];
                w = new double[o];
                EpnLagQuadrature.epn_glg_02_xiu(n, alpha, o, ref x, ref w);
                v = Monomial.monomial_value(n, o, x, expon);
                quad = typeMethods.r8vec_dot_product(o, w, v);
                err = Math.Abs(quad - exact);
                Console.WriteLine("  EPN_GLG_02_XIU: "
                                  + "  " + o.ToString().PadLeft(6)
                                  + "  " + quad.ToString().PadLeft(14)
                                  + "  " + err.ToString().PadLeft(14) + "");

                o = GolubWelsch_Xiu.gw_02_xiu_size(n);
                gamma0 = -1.0;
                delta0 = alpha + 1.0;
                c1 = -alpha - 1.0;
                volume_1d = typeMethods.r8_gamma(1.0 + alpha);
                x = new double[n * o];
                w = new double[o];
                GolubWelsch_Xiu.gw_02_xiu(n, o, gamma0, delta0, c1, volume_1d, ref x, ref w);
                v = Monomial.monomial_value(n, o, x, expon);
                quad = typeMethods.r8vec_dot_product(o, w, v);
                err = Math.Abs(quad - exact);
                Console.WriteLine("  GW_02_XIU:      "
                                  + "  " + o.ToString().PadLeft(6)
                                  + "  " + quad.ToString().PadLeft(14)
                                  + "  " + err.ToString().PadLeft(14) + "");
            }

            Console.WriteLine("  EXACT                   "
                              + "  " + exact.ToString().PadLeft(14) + "");

        }

        public static void epn_lag_test(int n, int[] expon)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EPN_LAG_TEST tests the rules for EPN with Laguerre weight on a monomial.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    27 January 2010
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
            string cout = "  EXPON = ";
            for (i = 0; i < n; i++)
            {
                cout += expon[i].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
            d = typeMethods.i4vec_sum(n, expon);
            Console.WriteLine("  Degree = " + d + "");
            Console.WriteLine("");

            exact = EpnLagIntegral.epn_lag_monomial_integral(n, expon);

            p = 0;

            if (d <= p)
            {
                o = EpnLagQuadrature.epn_lag_00_1_size(n);
                x = new double[n * o];
                w = new double[o];
                EpnLagQuadrature.epn_lag_00_1(n, o, ref x, ref w);
                v = Monomial.monomial_value(n, o, x, expon);
                quad = typeMethods.r8vec_dot_product(o, w, v);
                err = Math.Abs(quad - exact);
                Console.WriteLine("  EPN_LAG_00_1:   "
                                  + "  " + o.ToString().PadLeft(6)
                                  + "  " + quad.ToString().PadLeft(14)
                                  + "  " + err.ToString().PadLeft(14) + "");
            }

            p = 1;

            if (d <= p)
            {
                o = EpnLagQuadrature.epn_lag_01_1_size(n);
                x = new double[n * o];
                w = new double[o];
                EpnLagQuadrature.epn_lag_01_1(n, o, ref x, ref w);
                v = Monomial.monomial_value(n, o, x, expon);
                quad = typeMethods.r8vec_dot_product(o, w, v);
                err = Math.Abs(quad - exact);
                Console.WriteLine("  EPN_LAG_01_1:   "
                                  + "  " + o.ToString().PadLeft(6)
                                  + "  " + quad.ToString().PadLeft(14)
                                  + "  " + err.ToString().PadLeft(14) + "");
            }

            p = 2;

            if (d <= p)
            {
                o = EpnLagQuadrature.epn_lag_02_xiu_size(n);
                x = new double[n * o];
                w = new double[o];
                EpnLagQuadrature.epn_lag_02_xiu(n, o, ref x, ref w);
                v = Monomial.monomial_value(n, o, x, expon);
                quad = typeMethods.r8vec_dot_product(o, w, v);
                err = Math.Abs(quad - exact);
                Console.WriteLine("  EPN_LAG_02_XIU: "
                                  + "  " + o.ToString().PadLeft(6)
                                  + "  " + quad.ToString().PadLeft(14)
                                  + "  " + err.ToString().PadLeft(14) + "");

                o = GolubWelsch_Xiu.gw_02_xiu_size(n);
                gamma0 = -1.0;
                delta0 = 1.0;
                c1 = -1.0;
                volume_1d = 1.0;
                x = new double[n * o];
                w = new double[o];
                GolubWelsch_Xiu.gw_02_xiu(n, o, gamma0, delta0, c1, volume_1d, ref x, ref w);
                v = Monomial.monomial_value(n, o, x, expon);
                quad = typeMethods.r8vec_dot_product(o, w, v);
                err = Math.Abs(quad - exact);
                Console.WriteLine("  GW_02_XIU:      "
                                  + "  " + o.ToString().PadLeft(6)
                                  + "  " + quad.ToString().PadLeft(14)
                                  + "  " + err.ToString().PadLeft(14) + "");
            }

            Console.WriteLine("  EXACT                   "
                              + "  " + exact.ToString().PadLeft(14) + "");

        }


    }
}