using System;
using Burkardt.MonomialNS;
using Burkardt.Quadrature;
using Burkardt.Types;

namespace StroudTest
{
    using StroudQuadrature = Burkardt.Quadrature.Stroud;
    using StroudIntegral = Burkardt.IntegralNS.Stroud;

    public static class en
    {
        public static void en_r2_test(int n, int[] expon)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EN_R2_TEST tests the Stroud EN_R2 rules on a monomial.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 January 2010
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
            double pi = 3.141592653589793;
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
                cout += "  " + expon[i].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
            d = typeMethods.i4vec_sum(n, expon);
            Console.WriteLine("  Degree = " + d + "");
            Console.WriteLine("");

            exact = StroudIntegral.en_r2_monomial_integral(n, expon);

            p = 1;

            if (d <= p)
            {
                o = StroudQuadrature.en_r2_01_1_size(n);
                x = new double[n * o];
                w = new double[o];
                StroudQuadrature.en_r2_01_1(n, o, ref x, ref w);
                v = Monomial.monomial_value(n, o, x, expon);
                quad = typeMethods.r8vec_dot_product(o, w, v);
                err = Math.Abs(quad - exact);
                Console.WriteLine("  EN_R2_01_1:    "
                                  + "  " + o.ToString().PadLeft(6)
                                  + "  " + quad.ToString().PadLeft(14)
                                  + "  " + err.ToString().PadLeft(14) + "");
            }

            p = 2;

            if (d <= p)
            {
                o = Xiu.en_r2_02_xiu_size(n);
                x = new double[n * o];
                w = new double[o];
                Xiu.en_r2_02_xiu(n, o, ref x, ref w);
                v = Monomial.monomial_value(n, o, x, expon);
                quad = typeMethods.r8vec_dot_product(o, w, v);
                err = Math.Abs(quad - exact);
                Console.WriteLine("  EN_R2_02_XIU:  "
                                  + "  " + o.ToString().PadLeft(6)
                                  + "  " + quad.ToString().PadLeft(14)
                                  + "  " + err.ToString().PadLeft(14) + "");

                o = GolubWelsch_Xiu.gw_02_xiu_size(n);
                gamma0 = 2.0;
                delta0 = 0.0;
                c1 = 1.0;
                volume_1d = Math.Sqrt(pi);
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
                o = StroudQuadrature.en_r2_03_1_size(n);
                x = new double[n * o];
                w = new double[o];
                StroudQuadrature.en_r2_03_1(n, o, ref x, ref w);
                v = Monomial.monomial_value(n, o, x, expon);
                quad = typeMethods.r8vec_dot_product(o, w, v);
                err = Math.Abs(quad - exact);
                Console.WriteLine("  EN_R2_03_1:    "
                                  + "  " + o.ToString().PadLeft(6)
                                  + "  " + quad.ToString().PadLeft(14)
                                  + "  " + err.ToString().PadLeft(14) + "");

                o = StroudQuadrature.en_r2_03_2_size(n);
                x = new double[n * o];
                w = new double[o];
                StroudQuadrature.en_r2_03_2(n, o, ref x, ref w);
                v = Monomial.monomial_value(n, o, x, expon);
                quad = typeMethods.r8vec_dot_product(o, w, v);
                err = Math.Abs(quad - exact);
                Console.WriteLine("  EN_R2_03_2:    "
                                  + "  " + o.ToString().PadLeft(6)
                                  + "  " + quad.ToString().PadLeft(14)
                                  + "  " + err.ToString().PadLeft(14) + "");

                o = Xiu.en_r2_03_xiu_size(n);
                x = new double[n * o];
                w = new double[o];
                Xiu.en_r2_03_xiu(n, o, ref x, ref w);
                v = Monomial.monomial_value(n, o, x, expon);
                quad = typeMethods.r8vec_dot_product(o, w, v);
                err = Math.Abs(quad - exact);
                Console.WriteLine("  EN_R2_03_XIU:  "
                                  + "  " + o.ToString().PadLeft(6)
                                  + "  " + quad.ToString().PadLeft(14)
                                  + "  " + err.ToString().PadLeft(14) + "");
            }

            p = 5;

            if (d <= p)
            {
                if (2 <= n && n <= 7)
                {
                    option = 1;
                    o = StroudQuadrature.en_r2_05_1_size(n);
                    x = new double[n * o];
                    w = new double[o];
                    StroudQuadrature.en_r2_05_1(n, option, o, ref x, ref w);
                    v = Monomial.monomial_value(n, o, x, expon);
                    quad = typeMethods.r8vec_dot_product(o, w, v);
                    err = Math.Abs(quad - exact);
                    Console.WriteLine("  EN_R2_05_1(1): "
                                      + "  " + o.ToString().PadLeft(6)
                                      + "  " + quad.ToString().PadLeft(14)
                                      + "  " + err.ToString().PadLeft(14) + "");
                }

                if (n == 3 || n == 5 || n == 6)
                {
                    option = 2;
                    o = StroudQuadrature.en_r2_05_1_size(n);
                    x = new double[n * o];
                    w = new double[o];
                    StroudQuadrature.en_r2_05_1(n, option, o, ref x, ref w);
                    v = Monomial.monomial_value(n, o, x, expon);
                    quad = typeMethods.r8vec_dot_product(o, w, v);
                    err = Math.Abs(quad - exact);
                    Console.WriteLine("  EN_R2_05_1(2): "
                                      + "  " + o.ToString().PadLeft(6)
                                      + "  " + quad.ToString().PadLeft(14)
                                      + "  " + err.ToString().PadLeft(14) + "");
                }

                o = StroudQuadrature.en_r2_05_2_size(n);
                x = new double[n * o];
                w = new double[o];
                StroudQuadrature.en_r2_05_2(n, o, ref x, ref w);
                v = Monomial.monomial_value(n, o, x, expon);
                quad = typeMethods.r8vec_dot_product(o, w, v);
                err = Math.Abs(quad - exact);
                Console.WriteLine("  EN_R2_05_2:    "
                                  + "  " + o.ToString().PadLeft(6)
                                  + "  " + quad.ToString().PadLeft(14)
                                  + "  " + err.ToString().PadLeft(14) + "");

                if (3 <= n)
                {
                    o = StroudQuadrature.en_r2_05_3_size(n);
                    x = new double[n * o];
                    w = new double[o];
                    StroudQuadrature.en_r2_05_3(n, o, ref x, ref w);
                    v = Monomial.monomial_value(n, o, x, expon);
                    quad = typeMethods.r8vec_dot_product(o, w, v);
                    err = Math.Abs(quad - exact);
                    Console.WriteLine("  EN_R2_05_3:    "
                                      + "  " + o.ToString().PadLeft(6)
                                      + "  " + quad.ToString().PadLeft(14)
                                      + "  " + err.ToString().PadLeft(14) + "");
                }

                o = StroudQuadrature.en_r2_05_4_size(n);
                x = new double[n * o];
                w = new double[o];
                StroudQuadrature.en_r2_05_4(n, o, ref x, ref w);
                v = Monomial.monomial_value(n, o, x, expon);
                quad = typeMethods.r8vec_dot_product(o, w, v);
                err = Math.Abs(quad - exact);
                Console.WriteLine("  EN_R2_05_4:    "
                                  + "  " + o.ToString().PadLeft(6)
                                  + "  " + quad.ToString().PadLeft(14)
                                  + "  " + err.ToString().PadLeft(14) + "");

                o = StroudQuadrature.en_r2_05_5_size(n);
                x = new double[n * o];
                w = new double[o];
                StroudQuadrature.en_r2_05_5(n, o, ref x, ref w);
                v = Monomial.monomial_value(n, o, x, expon);
                quad = typeMethods.r8vec_dot_product(o, w, v);
                err = Math.Abs(quad - exact);
                Console.WriteLine("  EN_R2_05_5:    "
                                  + "  " + o.ToString().PadLeft(6)
                                  + "  " + quad.ToString().PadLeft(14)
                                  + "  " + err.ToString().PadLeft(14) + "");

                if (5 <= n)
                {
                    o = StroudQuadrature.en_r2_05_6_size(n);
                    x = new double[n * o];
                    w = new double[o];
                    StroudQuadrature.en_r2_05_6(n, o, ref x, ref w);
                    v = Monomial.monomial_value(n, o, x, expon);
                    quad = typeMethods.r8vec_dot_product(o, w, v);
                    err = Math.Abs(quad - exact);
                    Console.WriteLine("  EN_R2_05_6:    "
                                      + "  " + o.ToString().PadLeft(6)
                                      + "  " + quad.ToString().PadLeft(14)
                                      + "  " + err.ToString().PadLeft(14) + "");
                }
            }

            p = 7;

            if (d <= p)
            {
                if (n == 3 || n == 4 || n == 6 || n == 7)
                {
                    option = 1;
                    o = StroudQuadrature.en_r2_07_1_size(n);
                    x = new double[n * o];
                    w = new double[o];
                    StroudQuadrature.en_r2_07_1(n, option, o, ref x, ref w);
                    v = Monomial.monomial_value(n, o, x, expon);
                    quad = typeMethods.r8vec_dot_product(o, w, v);
                    err = Math.Abs(quad - exact);
                    Console.WriteLine("  EN_R2_07_1(1): "
                                      + "  " + o.ToString().PadLeft(6)
                                      + "  " + quad.ToString().PadLeft(14)
                                      + "  " + err.ToString().PadLeft(14) + "");
                }

                if (n == 3 || n == 4)
                {
                    option = 2;
                    o = StroudQuadrature.en_r2_07_1_size(n);
                    x = new double[n * o];
                    w = new double[o];
                    StroudQuadrature.en_r2_07_1(n, option, o, ref x, ref w);
                    v = Monomial.monomial_value(n, o, x, expon);
                    quad = typeMethods.r8vec_dot_product(o, w, v);
                    err = Math.Abs(quad - exact);
                    Console.WriteLine("  EN_R2_07_1(2): "
                                      + "  " + o.ToString().PadLeft(6)
                                      + "  " + quad.ToString().PadLeft(14)
                                      + "  " + err.ToString().PadLeft(14) + "");
                }

                if (3 <= n)
                {
                    o = StroudQuadrature.en_r2_07_2_size(n);
                    x = new double[n * o];
                    w = new double[o];
                    StroudQuadrature.en_r2_07_2(n, o, ref x, ref w);
                    v = Monomial.monomial_value(n, o, x, expon);
                    quad = typeMethods.r8vec_dot_product(o, w, v);
                    err = Math.Abs(quad - exact);
                    Console.WriteLine("  EN_R2_07_2:    "
                                      + "  " + o.ToString().PadLeft(6)
                                      + "  " + quad.ToString().PadLeft(14)
                                      + "  " + err.ToString().PadLeft(14) + "");
                }

                if (3 <= n && n <= 6)
                {
                    option = 1;
                    o = StroudQuadrature.en_r2_07_3_size(n);
                    x = new double[n * o];
                    w = new double[o];
                    StroudQuadrature.en_r2_07_3(n, option, o, ref x, ref w);
                    v = Monomial.monomial_value(n, o, x, expon);
                    quad = typeMethods.r8vec_dot_product(o, w, v);
                    err = Math.Abs(quad - exact);
                    Console.WriteLine("  EN_R2_07_3(1): "
                                      + "  " + o.ToString().PadLeft(6)
                                      + "  " + quad.ToString().PadLeft(14)
                                      + "  " + err.ToString().PadLeft(14) + "");
                }

                if (n == 3 || n == 4)
                {
                    option = 2;
                    o = StroudQuadrature.en_r2_07_3_size(n);
                    x = new double[n * o];
                    w = new double[o];
                    StroudQuadrature.en_r2_07_3(n, option, o, ref x, ref w);
                    v = Monomial.monomial_value(n, o, x, expon);
                    quad = typeMethods.r8vec_dot_product(o, w, v);
                    err = Math.Abs(quad - exact);
                    Console.WriteLine("  EN_R2_07_3(2): "
                                      + "  " + o.ToString().PadLeft(6)
                                      + "  " + quad.ToString().PadLeft(14)
                                      + "  " + err.ToString().PadLeft(14) + "");
                }
            }

            p = 9;

            if (d <= p)
            {
                if (3 <= n && n <= 6)
                {
                    option = 1;
                    o = StroudQuadrature.en_r2_09_1_size(n);
                    x = new double[n * o];
                    w = new double[o];
                    StroudQuadrature.en_r2_09_1(n, option, o, ref x, ref w);
                    v = Monomial.monomial_value(n, o, x, expon);
                    quad = typeMethods.r8vec_dot_product(o, w, v);
                    err = Math.Abs(quad - exact);
                    Console.WriteLine("  EN_R2_09_1(1): "
                                      + "  " + o.ToString().PadLeft(6)
                                      + "  " + quad.ToString().PadLeft(14)
                                      + "  " + err.ToString().PadLeft(14) + "");

                    option = 2;
                    o = StroudQuadrature.en_r2_09_1_size(n);
                    x = new double[n * o];
                    w = new double[o];
                    StroudQuadrature.en_r2_09_1(n, option, o, ref x, ref w);
                    v = Monomial.monomial_value(n, o, x, expon);
                    quad = typeMethods.r8vec_dot_product(o, w, v);
                    err = Math.Abs(quad - exact);
                    Console.WriteLine("  EN_R2_09_1(2): "
                                      + "  " + o.ToString().PadLeft(6)
                                      + "  " + quad.ToString().PadLeft(14)
                                      + "  " + err.ToString().PadLeft(14) + "");
                }
            }

            p = 11;

            if (d <= p)
            {
                if (3 <= n && n <= 5)
                {
                    option = 1;
                    o = StroudQuadrature.en_r2_11_1_size(n);
                    x = new double[n * o];
                    w = new double[o];
                    StroudQuadrature.en_r2_11_1(n, option, o, ref x, ref w);
                    v = Monomial.monomial_value(n, o, x, expon);
                    quad = typeMethods.r8vec_dot_product(o, w, v);
                    err = Math.Abs(quad - exact);
                    Console.WriteLine("  EN_R2_11_1(1): "
                                      + "  " + o.ToString().PadLeft(6)
                                      + "  " + quad.ToString().PadLeft(14)
                                      + "  " + err.ToString().PadLeft(14) + "");

                    option = 2;
                    o = StroudQuadrature.en_r2_11_1_size(n);
                    x = new double[n * o];
                    w = new double[o];
                    StroudQuadrature.en_r2_11_1(n, option, o, ref x, ref w);
                    v = Monomial.monomial_value(n, o, x, expon);
                    quad = typeMethods.r8vec_dot_product(o, w, v);
                    err = Math.Abs(quad - exact);
                    Console.WriteLine("  EN_R2_11_1(2): "
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