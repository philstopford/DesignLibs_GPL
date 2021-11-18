using System;
using Burkardt.IntegralNS;
using Burkardt.Quadrature;
using Burkardt.Types;
using Monomial = Burkardt.MonomialNS.Monomial;
using Stroud = Burkardt.Quadrature.Stroud;

namespace SandiaCubatureTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SANDIA_CUBATURE_TEST.
        //
        //  Discussion:
        //
        //    SANDIA_CUBATURE_TEST tests the SANDIA_CUBATURE library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 February 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("SANDIA_CUBATURE_TEST");
        Console.WriteLine("  Test the SANDIA_CUBATURE library.");

        cn_geg_tests();
        cn_jac_tests();
        cn_leg_tests();
        en_her_tests();
        epn_glg_tests();
        epn_lag_tests();
        gw_tests();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("SANDIA_CUBATURE_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void cn_geg_tests()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CN_GEG_TESTS tests the rules for CN with Gegenbauer weight on monomials.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int TEST_NUM = 5;

        double alpha;
        double[] alpha_test = {-0.5, 0.0, 0.5, 1.0, 1.5};
        int[] expon = null;
        int n;
        int test;

        Console.WriteLine("");
        Console.WriteLine("CN_GEG_TESTS");
        Console.WriteLine("  Demonstrate the use of quadrature rules for the region");
        Console.WriteLine("  CN_GEG, that is, the hypercube [-1,+1]^N, with the");
        Console.WriteLine("  weight W(ALPHA;X) = product ( 1 <= I <= N )");
        Console.WriteLine("    (1-X(I)^2)^ALPHA");
        Console.WriteLine("");
        Console.WriteLine("  We use the formulas to integrate various monomials of");
        Console.WriteLine("  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)");
        Console.WriteLine("  and compare to the exact integral.");
        Console.WriteLine("");
        Console.WriteLine("  The precision of each formula is known, and we only use");
        Console.WriteLine("  a formula if its precision indicates it should be able to");
        Console.WriteLine("  produce an exact result.");

        for (n = 1; n <= 6; n++)
        {
            expon = new int[n];

            for (test = 0; test < TEST_NUM; test++)
            {
                alpha = alpha_test[test];

                typeMethods.i4vec_zero(n, ref expon);
                cn_geg_test(n, alpha, expon);
            }

            for (test = 0; test < TEST_NUM; test++)
            {
                alpha = alpha_test[test];

                typeMethods.i4vec_zero(n, ref expon);
                expon[n - 1] = 1;
                cn_geg_test(n, alpha, expon);
            }

            switch (n)
            {
                case >= 2:
                {
                    for (test = 0; test < TEST_NUM; test++)
                    {
                        alpha = alpha_test[test];

                        typeMethods.i4vec_zero(n, ref expon);
                        expon[0] = 1;
                        expon[1] = 1;
                        cn_geg_test(n, alpha, expon);
                    }

                    break;
                }
            }

            for (test = 0; test < TEST_NUM; test++)
            {
                alpha = alpha_test[test];

                typeMethods.i4vec_zero(n, ref expon);
                expon[0] = 2;
                cn_geg_test(n, alpha, expon);

            }
        }
    }

    private static void cn_geg_test(int n, double alpha, int[] expon)

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
            cout += expon[i].ToString(CultureInfo.InvariantCulture).PadLeft(4);
        }

        Console.WriteLine(cout);
        cout = "";
        d = typeMethods.i4vec_sum(n, expon);
        Console.WriteLine("  Degree = " + d + "");
        Console.WriteLine("");

        exact = CN_Geg.cn_geg_monomial_integral(n, alpha, expon);

        p = 1;

        if (d <= p)
        {
            o = Stroud.cn_geg_01_1_size(n, alpha);
            x = new double[n * o];
            w = new double[o];
            Stroud.cn_geg_01_1(n, alpha, o, ref x, ref w);
            v = Monomial.monomial_value(n, o, x, expon);
            quad = typeMethods.r8vec_dot_product(o, w, v);
            err = typeMethods.r8_abs(quad - exact);
            Console.WriteLine("  CN_GEG_01_1:   "
                              + "  " + o.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
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
            err = typeMethods.r8_abs(quad - exact);
            Console.WriteLine("  CN_GEG_02_XIU: "
                              + "  " + o.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

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
            err = typeMethods.r8_abs(quad - exact);
            Console.WriteLine("  GW_02_XIU:     "
                              + "  " + o.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
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
            err = typeMethods.r8_abs(quad - exact);
            Console.WriteLine("  CN_GEG_03_XIU: "
                              + "  " + o.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        Console.WriteLine("  EXACT                  "
                          + "  " + exact.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
    }

    private static void cn_jac_tests()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CN_JAC_TESTS tests the rules for CN with Jacobi weight on monomials.
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
        int TEST_NUM = 4;

        double alpha;
        double[] alpha_test = {0.0, 1.0, 0.0, 0.5};
        double beta;
        double[] beta_test = {0.0, 0.0, 2.0, 1.5};
        int[] expon;
        int n;
        int test;

        Console.WriteLine("");
        Console.WriteLine("CN_JAC_TESTS");
        Console.WriteLine("  Demonstrate the use of quadrature rules for the region");
        Console.WriteLine("  CN_JAC, that is, the hypercube [-1,+1]^N, with the");
        Console.WriteLine("  weight W(ALPHA,BETA;X) = product ( 1 <= I <= N )");
        Console.WriteLine("    (1-X(I))^ALPHA (1+X(I))^BETA");
        Console.WriteLine("");
        Console.WriteLine("  We use the formulas to integrate various monomials of");
        Console.WriteLine("  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)");
        Console.WriteLine("  and compare to the exact integral.");
        Console.WriteLine("");
        Console.WriteLine("  The precision of each formula is known, and we only use");
        Console.WriteLine("  a formula if its precision indicates it should be able to");
        Console.WriteLine("  produce an exact result.");

        for (n = 1; n <= 6; n++)
        {
            expon = new int[n];

            for (test = 0; test < TEST_NUM; test++)
            {
                alpha = alpha_test[test];
                beta = beta_test[test];

                typeMethods.i4vec_zero(n, ref expon);
                cn_jac_test(n, alpha, beta, expon);
            }

            for (test = 0; test < TEST_NUM; test++)
            {
                alpha = alpha_test[test];
                beta = beta_test[test];

                typeMethods.i4vec_zero(n, ref expon);
                expon[n - 1] = 1;
                cn_jac_test(n, alpha, beta, expon);
            }

            switch (n)
            {
                case >= 2:
                {
                    for (test = 0; test < TEST_NUM; test++)
                    {
                        alpha = alpha_test[test];
                        beta = beta_test[test];

                        typeMethods.i4vec_zero(n, ref expon);
                        expon[0] = 1;
                        expon[1] = 1;
                        cn_jac_test(n, alpha, beta, expon);
                    }

                    break;
                }
            }

            for (test = 0; test < TEST_NUM; test++)
            {
                alpha = alpha_test[test];
                beta = beta_test[test];

                typeMethods.i4vec_zero(n, ref expon);
                expon[0] = 2;
                cn_jac_test(n, alpha, beta, expon);

            }
        }
    }

    private static void cn_jac_test(int n, double alpha, double beta, int[] expon)

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
            cout += expon[i].ToString(CultureInfo.InvariantCulture).PadLeft(4);
        }

        Console.WriteLine(cout);
        cout = "";
        d = typeMethods.i4vec_sum(n, expon);
        Console.WriteLine("  Degree = " + d + "");
        Console.WriteLine("");

        exact = CN_Geg.cn_jac_monomial_integral(n, alpha, beta, expon);

        p = 1;

        if (d <= p)
        {
            o = CN_Jac.cn_jac_01_1_size(n, alpha, beta);
            x = new double[n * o];
            w = new double[o];
            CN_Jac.cn_jac_01_1(n, alpha, beta, o, ref x, ref w);
            v = Monomial.monomial_value(n, o, x, expon);
            quad = typeMethods.r8vec_dot_product(o, w, v);
            err = typeMethods.r8_abs(quad - exact);
            Console.WriteLine("  CN_JAC_01_1:   "
                              + "  " + o.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
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
            err = typeMethods.r8_abs(quad - exact);
            Console.WriteLine("  CN_JAC_02_XIU: "
                              + "  " + o.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            o = GolubWelsch_Xiu.gw_02_xiu_size(n);
            gamma0 = (alpha + beta + 2.0) / 2.0;
            delta0 = (alpha - beta) / 2.0;
            c1 = 2.0 * (alpha + 1.0) * (beta + 1.0) / (alpha + beta + 3.0)
                                                    / (alpha + beta + 2.0);
            volume_1d = Math.Pow(2.0, alpha + beta + 1.0) * typeMethods.r8_gamma(alpha + 1.0)
                                                          * typeMethods.r8_gamma(beta + 1.0) /
                        (alpha + beta + 1.0) / typeMethods.r8_gamma(alpha + beta + 1.0);
            x = new double[n * o];
            w = new double[o];
            GolubWelsch_Xiu.gw_02_xiu(n, o, gamma0, delta0, c1, volume_1d, ref x, ref w);
            v = Monomial.monomial_value(n, o, x, expon);
            quad = typeMethods.r8vec_dot_product(o, w, v);
            err = typeMethods.r8_abs(quad - exact);
            Console.WriteLine("  GW_02_XIU:     "
                              + "  " + o.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        Console.WriteLine("  EXACT                  "
                          + "  " + exact.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void cn_leg_tests()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CN_LEG_TESTS tests the rules for CN with Legendre weight on monomials.
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
        int[] expon;
        int n;

        Console.WriteLine("");
        Console.WriteLine("CN_LEG_TESTS");
        Console.WriteLine("  Demonstrate the use of quadrature rules for the region");
        Console.WriteLine("  CN_LEG, that is, the hypercube [-1,+1]^N, with the");
        Console.WriteLine("  Legendre weight W(X) = 1.");
        Console.WriteLine("");
        Console.WriteLine("  We use the formulas to integrate various monomials of");
        Console.WriteLine("  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)");
        Console.WriteLine("  and compare to the exact integral.");
        Console.WriteLine("");
        Console.WriteLine("  The precision of each formula is known, and we only use");
        Console.WriteLine("  a formula if its precision indicates it should be able to");
        Console.WriteLine("  produce an exact result.");

        for (n = 1; n <= 6; n++)
        {
            expon = new int[n];

            typeMethods.i4vec_zero(n, ref expon);
            cn_leg_test(n, expon);

            typeMethods.i4vec_zero(n, ref expon);
            expon[n - 1] = 1;
            cn_leg_test(n, expon);

            switch (n)
            {
                case >= 2:
                    typeMethods.i4vec_zero(n, ref expon);
                    expon[0] = 1;
                    expon[1] = 1;
                    cn_leg_test(n, expon);
                    break;
            }

            typeMethods.i4vec_zero(n, ref expon);
            expon[0] = 2;
            cn_leg_test(n, expon);

            typeMethods.i4vec_zero(n, ref expon);
            expon[0] = 3;
            cn_leg_test(n, expon);

            typeMethods.i4vec_zero(n, ref expon);
            expon[n - 1] = 4;
            cn_leg_test(n, expon);

            switch (n)
            {
                case >= 2:
                    typeMethods.i4vec_zero(n, ref expon);
                    expon[0] = 3;
                    expon[1] = 2;
                    cn_leg_test(n, expon);
                    break;
            }
        }
    }

    private static void cn_leg_test(int n, int[] expon)

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
            cout += expon[i].ToString(CultureInfo.InvariantCulture).PadLeft(4);
        }

        Console.WriteLine(cout);
        cout = "";
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
            err = typeMethods.r8_abs(quad - exact);
            Console.WriteLine("  CN_LEG_01_1:   "
                              + "  " + o.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
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
            err = typeMethods.r8_abs(quad - exact);
            Console.WriteLine("  CN_LEG_02_XIU: "
                              + "  " + o.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

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
            err = typeMethods.r8_abs(quad - exact);
            Console.WriteLine("  GW_02_XIU:     "
                              + "  " + o.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        p = 3;

        if (d <= p)
        {
            o = Stroud.cn_leg_03_1_size(n);
            x = new double[n * o];
            w = new double[o];
            Stroud.cn_leg_03_1(n, o, ref x, ref w);
            v = Monomial.monomial_value(n, o, x, expon);
            quad = typeMethods.r8vec_dot_product(o, w, v);
            err = typeMethods.r8_abs(quad - exact);
            Console.WriteLine("  CN_LEG_03_1:   "
                              + "  " + o.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            o = Xiu.cn_leg_03_xiu_size(n);
            x = new double[n * o];
            w = new double[o];
            Xiu.cn_leg_03_xiu(n, o, ref x, ref w);
            v = Monomial.monomial_value(n, o, x, expon);
            quad = typeMethods.r8vec_dot_product(o, w, v);
            err = typeMethods.r8_abs(quad - exact);
            Console.WriteLine("  CN_LEG_03_XIU: "
                              + "  " + o.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        p = 5;

        if (d <= p)
        {
            switch (n)
            {
                case >= 4 and <= 6:
                    o = Stroud.cn_leg_05_1_size(n);
                    x = new double[n * o];
                    w = new double[o];
                    option = 1;
                    Stroud.cn_leg_05_1(n, option, o, ref x, ref w);
                    v = Monomial.monomial_value(n, o, x, expon);
                    quad = typeMethods.r8vec_dot_product(o, w, v);
                    err = typeMethods.r8_abs(quad - exact);
                    Console.WriteLine("  CN_LEG_05_1(1):"
                                      + "  " + o.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                      + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                      + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                    break;
            }

            switch (n)
            {
                case >= 4 and <= 5:
                    o = Stroud.cn_leg_05_1_size(n);
                    x = new double[n * o];
                    w = new double[o];
                    option = 2;
                    Stroud.cn_leg_05_1(n, option, o, ref x, ref w);
                    v = Monomial.monomial_value(n, o, x, expon);
                    quad = typeMethods.r8vec_dot_product(o, w, v);
                    err = typeMethods.r8_abs(quad - exact);
                    Console.WriteLine("  CN_LEG_05_1(2):"
                                      + "  " + o.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                      + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                      + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                    break;
            }

            switch (n)
            {
                case >= 2:
                    o = Stroud.cn_leg_05_2_size(n);
                    x = new double[n * o];
                    w = new double[o];
                    Stroud.cn_leg_05_2(n, o, ref x, ref w);
                    v = Monomial.monomial_value(n, o, x, expon);
                    quad = typeMethods.r8vec_dot_product(o, w, v);
                    err = typeMethods.r8_abs(quad - exact);
                    Console.WriteLine("  CN_LEG_05_2:   "
                                      + "  " + o.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                      + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                      + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                    break;
            }
        }

        Console.WriteLine("  EXACT                  "
                          + "  " + exact.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void en_her_tests()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EN_HER_TESTS tests the Stroud EN_HER rules on monomials.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int d;
        int[] expon;
        int i;
        int n;

        Console.WriteLine("");
        Console.WriteLine("EN_HER_TESTS");
        Console.WriteLine("  Demonstrate the use of Stroud rules for the region");
        Console.WriteLine("  EN_HER, that is, all of N-dimensional space, with the");
        Console.WriteLine("  weight function W(X) = exp ( - X1^2 - X2^2 ... -XN^2 )");
        Console.WriteLine("");
        Console.WriteLine("  We use the formulas to integrate various monomials of");
        Console.WriteLine("  the form X1^ALPHA1 * X2^ALPHA2 * ... XN^ALPHAN");
        Console.WriteLine("  and compare to the exact integral.");
        Console.WriteLine("");
        Console.WriteLine("  The precision of each formula is known, and we only use");
        Console.WriteLine("  a formula if its precision indicates it should be able to");
        Console.WriteLine("  produce an exact result.");

        for (n = 1; n <= 7; n++)
        {
            expon = new int[n];

            for (i = 0; i < n; i++)
            {
                expon[i] = 0;
            }

            en_her_test(n, expon);

            for (i = 0; i < n; i++)
            {
                expon[i] = 0;
            }

            expon[0] = 2;
            en_her_test(n, expon);

            for (i = 0; i < n; i++)
            {
                expon[i] = 0;
            }

            expon[1] = 4;
            en_her_test(n, expon);

            for (i = 0; i < n; i++)
            {
                expon[i] = i + 1;
            }

            d = typeMethods.i4vec_sum(n, expon);
            switch (d)
            {
                case <= 5:
                    en_her_test(n, expon);
                    break;
            }

            for (i = 0; i < n; i++)
            {
                expon[i] = 2;
            }

            d = typeMethods.i4vec_sum(n, expon);
            switch (d)
            {
                case <= 5:
                    en_her_test(n, expon);
                    break;
            }

        }
    }

    private static void en_her_test(int n, int[] expon)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EN_HER_TEST tests the Stroud EN_HER rules on a monomial.
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
            cout += "  " + expon[i].ToString(CultureInfo.InvariantCulture).PadLeft(4);
        }

        Console.WriteLine(cout);
        cout = "";
        d = typeMethods.i4vec_sum(n, expon);
        Console.WriteLine("  Degree = " + d + "");
        Console.WriteLine("");

        exact = Burkardt.IntegralNS.Stroud.en_her_monomial_integral(n, expon);

        p = 1;

        if (d <= p)
        {
            o = Stroud.en_her_01_1_size(n);
            x = new double[n * o];
            w = new double[o];
            Stroud.en_her_01_1(n, o, ref x, ref w);
            v = Monomial.monomial_value(n, o, x, expon);
            quad = typeMethods.r8vec_dot_product(o, w, v);
            err = typeMethods.r8_abs(quad - exact);
            Console.WriteLine("  EN_HER_01_1:   "
                              + "  " + o.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        p = 2;

        if (d <= p)
        {
            o = Xiu.en_her_02_xiu_size(n);
            x = new double[n * o];
            w = new double[o];
            Xiu.en_her_02_xiu(n, o, ref x, ref w);
            v = Monomial.monomial_value(n, o, x, expon);
            quad = typeMethods.r8vec_dot_product(o, w, v);
            err = typeMethods.r8_abs(quad - exact);
            Console.WriteLine("  EN_HER_02_XIU: "
                              + "  " + o.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

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
            err = typeMethods.r8_abs(quad - exact);
            Console.WriteLine("  GW_02_XIU:     "
                              + "  " + o.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        p = 3;

        if (d <= p)
        {
            o = Stroud.en_her_03_1_size(n);
            x = new double[n * o];
            w = new double[o];
            Stroud.en_her_03_1(n, o, ref x, ref w);
            v = Monomial.monomial_value(n, o, x, expon);
            quad = typeMethods.r8vec_dot_product(o, w, v);
            err = typeMethods.r8_abs(quad - exact);
            Console.WriteLine("  EN_HER_03_1:   "
                              + "  " + o.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            o = Xiu.en_her_03_xiu_size(n);
            x = new double[n * o];
            w = new double[o];
            Xiu.en_her_03_xiu(n, o, ref x, ref w);
            v = Monomial.monomial_value(n, o, x, expon);
            quad = typeMethods.r8vec_dot_product(o, w, v);
            err = typeMethods.r8_abs(quad - exact);
            Console.WriteLine("  EN_HER_03_XIU: "
                              + "  " + o.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        p = 5;

        if (d <= p)
        {
            switch (n)
            {
                case >= 2 and <= 7:
                    option = 1;
                    o = Stroud.en_her_05_1_size(n);
                    x = new double[n * o];
                    w = new double[o];
                    Stroud.en_her_05_1(n, option, o, ref x, ref w);
                    v = Monomial.monomial_value(n, o, x, expon);
                    quad = typeMethods.r8vec_dot_product(o, w, v);
                    err = typeMethods.r8_abs(quad - exact);
                    Console.WriteLine("  EN_HER_05_1(1):"
                                      + "  " + o.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                      + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                      + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                    break;
            }

            switch (n)
            {
                case 3:
                case 5:
                case 6:
                    option = 2;
                    o = Stroud.en_her_05_1_size(n);
                    x = new double[n * o];
                    w = new double[o];
                    Stroud.en_her_05_1(n, option, o, ref x, ref w);
                    v = Monomial.monomial_value(n, o, x, expon);
                    quad = typeMethods.r8vec_dot_product(o, w, v);
                    err = typeMethods.r8_abs(quad - exact);
                    Console.WriteLine("  EN_HER_05_1(2):"
                                      + "  " + o.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                      + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                      + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                    break;
            }

            o = Stroud.en_her_05_2_size(n);
            x = new double[n * o];
            w = new double[o];
            Stroud.en_her_05_2(n, o, ref x, ref w);
            v = Monomial.monomial_value(n, o, x, expon);
            quad = typeMethods.r8vec_dot_product(o, w, v);
            err = typeMethods.r8_abs(quad - exact);
            Console.WriteLine("  EN_HER_05_2:   "
                              + "  " + o.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        Console.WriteLine("  EXACT                  "
                          + "  " + exact.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void epn_glg_tests()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EPN_GLG_TESTS tests the rules for EPN with GLG on monomials.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int TEST_NUM = 5;

        double alpha;
        double[] alpha_test = {-0.5, 0.0, 0.5, 1.0, 2.0};
        int[] expon;
        int n;
        int test;

        Console.WriteLine("");
        Console.WriteLine("EPN_GLG_TESTS");
        Console.WriteLine("  Demonstrate the use of quadrature rules for the region");
        Console.WriteLine("  EPN_GLG, that is, the positive half space [0,+oo)^N, with the");
        Console.WriteLine("  weight W(ALPHA;X) = product ( 1 <= I <= N ) X(I)^ALPHA exp ( -X(I) )");
        Console.WriteLine("");
        Console.WriteLine("  We use the formulas to integrate various monomials of");
        Console.WriteLine("  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)");
        Console.WriteLine("  and compare to the exact integral.");
        Console.WriteLine("");
        Console.WriteLine("  The precision of each formula is known, and we only use");
        Console.WriteLine("  a formula if its precision indicates it should be able to");
        Console.WriteLine("  produce an exact result.");

        for (n = 1; n <= 6; n++)
        {
            expon = new int[n];

            typeMethods.i4vec_zero(n, ref expon);
            for (test = 0; test < TEST_NUM; test++)
            {
                alpha = alpha_test[test];
                epn_glg_test(n, expon, alpha);
            }

            typeMethods.i4vec_zero(n, ref expon);
            expon[n - 1] = 1;
            for (test = 0; test < TEST_NUM; test++)
            {
                alpha = alpha_test[test];
                epn_glg_test(n, expon, alpha);
            }

            switch (n)
            {
                case >= 2:
                {
                    typeMethods.i4vec_zero(n, ref expon);
                    expon[0] = 1;
                    expon[1] = 1;
                    for (test = 0; test < TEST_NUM; test++)
                    {
                        alpha = alpha_test[test];
                        epn_glg_test(n, expon, alpha);
                    }

                    break;
                }
            }

            typeMethods.i4vec_zero(n, ref expon);
            expon[0] = 2;
            for (test = 0; test < TEST_NUM; test++)
            {
                alpha = alpha_test[test];
                epn_glg_test(n, expon, alpha);
            }
        }
    }

    private static void epn_glg_test(int n, int[] expon, double alpha)

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
        string cout = "  EXPON = ";
        for (i = 0; i < n; i++)
        {
            cout += expon[i].ToString(CultureInfo.InvariantCulture).PadLeft(4);
        }

        Console.WriteLine(cout);
        cout = "";
        d = typeMethods.i4vec_sum(n, expon);
        Console.WriteLine("  Degree = " + d + "");
        Console.WriteLine("");

        exact = Burkardt.IntegralNS.Epn_Lag.epn_glg_monomial_integral(n, expon, alpha);

        p = 1;

        if (d <= p)
        {
            o = Burkardt.Quadrature.Epn_Lag.epn_glg_01_1_size(n, alpha);
            x = new double[n * o];
            w = new double[o];
            Burkardt.Quadrature.Epn_Lag.epn_glg_01_1(n, alpha, o, ref x, ref w);
            v = Monomial.monomial_value(n, o, x, expon);
            quad = typeMethods.r8vec_dot_product(o, w, v);
            err = typeMethods.r8_abs(quad - exact);
            Console.WriteLine("  EPN_GLG_01_1:   "
                              + "  " + o.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        p = 2;

        if (d <= p)
        {
            o = Burkardt.Quadrature.Epn_Lag.epn_glg_02_xiu_size(n, alpha);
            x = new double[n * o];
            w = new double[o];
            Burkardt.Quadrature.Epn_Lag.epn_glg_02_xiu(n, alpha, o, ref x, ref w);
            v = Monomial.monomial_value(n, o, x, expon);
            quad = typeMethods.r8vec_dot_product(o, w, v);
            err = typeMethods.r8_abs(quad - exact);
            Console.WriteLine("  EPN_GLG_02_XIU: "
                              + "  " + o.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

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
            err = typeMethods.r8_abs(quad - exact);
            Console.WriteLine("  GW_02_XIU:      "
                              + "  " + o.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        Console.WriteLine("  EXACT                   "
                          + "  " + exact.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void epn_lag_tests()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EPN_LAG_TESTS tests the rules for EPN with Laguerre weight on monomials.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] expon;
        int n;

        Console.WriteLine("");
        Console.WriteLine("EPN_LAG_TESTS");
        Console.WriteLine("  Demonstrate the use of quadrature rules for the region");
        Console.WriteLine("  EPN_LAG, that is, the positive half space [0,+oo)^N, with the");
        Console.WriteLine("  weight W(X) = product ( 1 <= I <= N ) exp ( -X(I) )");
        Console.WriteLine("");
        Console.WriteLine("  We use the formulas to integrate various monomials of");
        Console.WriteLine("  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)");
        Console.WriteLine("  and compare to the exact integral.");
        Console.WriteLine("");
        Console.WriteLine("  The precision of each formula is known, and we only use");
        Console.WriteLine("  a formula if its precision indicates it should be able to");
        Console.WriteLine("  produce an exact result.");

        for (n = 1; n <= 6; n++)
        {
            expon = new int[n];

            typeMethods.i4vec_zero(n, ref expon);
            epn_lag_test(n, expon);

            typeMethods.i4vec_zero(n, ref expon);
            expon[n - 1] = 1;
            epn_lag_test(n, expon);

            switch (n)
            {
                case >= 2:
                    typeMethods.i4vec_zero(n, ref expon);
                    expon[0] = 1;
                    expon[1] = 1;
                    epn_lag_test(n, expon);
                    break;
            }

            typeMethods.i4vec_zero(n, ref expon);
            expon[0] = 2;
            epn_lag_test(n, expon);

        }
    }

    private static void epn_lag_test(int n, int[] expon)

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
        string cout = "  EXPON = ";
        for (i = 0; i < n; i++)
        {
            cout += expon[i].ToString(CultureInfo.InvariantCulture).PadLeft(4);
        }

        Console.WriteLine(cout);
        cout = "";
        d = typeMethods.i4vec_sum(n, expon);
        Console.WriteLine("  Degree = " + d + "");
        Console.WriteLine("");

        exact = Burkardt.IntegralNS.Epn_Lag.epn_lag_monomial_integral(n, expon);

        p = 1;

        if (d <= p)
        {
            o = Burkardt.Quadrature.Epn_Lag.epn_lag_01_1_size(n);
            x = new double[n * o];
            w = new double[o];
            Burkardt.Quadrature.Epn_Lag.epn_lag_01_1(n, o, ref x, ref w);
            v = Monomial.monomial_value(n, o, x, expon);
            quad = typeMethods.r8vec_dot_product(o, w, v);
            err = typeMethods.r8_abs(quad - exact);
            Console.WriteLine("  EPN_LAG_01_1:   "
                              + "  " + o.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        p = 2;

        if (d <= p)
        {
            o = Burkardt.Quadrature.Epn_Lag.epn_lag_02_xiu_size(n);
            x = new double[n * o];
            w = new double[o];
            Burkardt.Quadrature.Epn_Lag.epn_lag_02_xiu(n, o, ref x, ref w);
            v = Monomial.monomial_value(n, o, x, expon);
            quad = typeMethods.r8vec_dot_product(o, w, v);
            err = typeMethods.r8_abs(quad - exact);
            Console.WriteLine("  EPN_LAG_02_XIU: "
                              + "  " + o.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

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
            err = typeMethods.r8_abs(quad - exact);
            Console.WriteLine("  GW_02_XIU:      "
                              + "  " + o.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        Console.WriteLine("  EXACT                   "
                          + "  " + exact.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void gw_tests()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GW_TESTS tests the rules for GW on monomials.
        //
        //  Discussion:
        //
        //    Right now, this test simply calls the GW rule for each type of
        //    weight function for which the orthogonal polynomials are known.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 February 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double alpha;
        double beta;
        double c1;
        double delta0;
        double gamma0;
        int i;
        int j;
        int n;
        int o;
        double pi = 3.141592653589793;
        double volume_1d;
        double[] w;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("GW_TESTS");
        Console.WriteLine("  Demonstrate the use of quadrature rules for a Golub Welsch rule");
        Console.WriteLine("  defined over some interval and some weight function for which");
        Console.WriteLine("  the three term recursion of the orthogonal polynomials is known.");
        //
        //  For a given dimension N, the rule is always the same size.
        //
        n = 2;
        o = GolubWelsch_Xiu.gw_02_xiu_size(n);
        x = new double[n * o];
        w = new double[o];
        //
        //  Chebyshev Type 1.
        //
        gamma0 = 1.0;
        delta0 = 0.0;
        c1 = 0.5;
        volume_1d = pi;
        GolubWelsch_Xiu.gw_02_xiu(n, o, gamma0, delta0, c1, volume_1d, ref x, ref w);
        Console.WriteLine("");
        Console.WriteLine("  Chebyshev1:");
        Console.WriteLine("");
        for (j = 0; j < o; j++)
        {
            string cout = w[j].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            for (i = 0; i < n; i++)
            {
                cout += x[i + j * n].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);
        }

        //
        //  Chebyshev Type 2.
        //
        gamma0 = 2.0;
        delta0 = 0.0;
        c1 = 0.5;
        volume_1d = pi / 2.0;
        GolubWelsch_Xiu.gw_02_xiu(n, o, gamma0, delta0, c1, volume_1d, ref x, ref w);
        Console.WriteLine("");
        Console.WriteLine("  Chebyshev2:");
        Console.WriteLine("");
        for (j = 0; j < o; j++)
        {
            string cout = w[j].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            for (i = 0; i < n; i++)
            {
                cout += x[i + j * n].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);
        }

        //
        //  Gegenbauer
        //
        alpha = 1.0;
        gamma0 = 1.0;
        delta0 = 0.0;
        c1 = 1.0 / (2.0 * alpha + 3.0);
        volume_1d = Math.Sqrt(pi) * typeMethods.r8_gamma(alpha + 1.0) / typeMethods.r8_gamma(alpha + 1.5);
        GolubWelsch_Xiu.gw_02_xiu(n, o, gamma0, delta0, c1, volume_1d, ref x, ref w);
        Console.WriteLine("");
        Console.WriteLine("  Gegenbauer:");
        Console.WriteLine("  ALPHA = " + alpha + "");
        Console.WriteLine("");
        for (j = 0; j < o; j++)
        {
            string cout = w[j].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            for (i = 0; i < n; i++)
            {
                cout += x[i + j * n].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);
        }

        //
        //  Generalized Hermite.
        //
        alpha = 1.0;
        gamma0 = 2.0;
        delta0 = 0.0;
        c1 = 2.0 + 2.0 * alpha;
        volume_1d = typeMethods.r8_gamma((alpha + 1.0) / 2.0);
        GolubWelsch_Xiu.gw_02_xiu(n, o, gamma0, delta0, c1, volume_1d, ref x, ref w);
        Console.WriteLine("");
        Console.WriteLine("  Generalized Hermite:");
        Console.WriteLine("  ALPHA = " + alpha + "");
        Console.WriteLine("");
        for (j = 0; j < o; j++)
        {
            string cout = w[j].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            for (i = 0; i < n; i++)
            {
                cout += x[i + j * n].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);
        }

        //
        //  Generalized Laguerre.
        //
        alpha = 1.0;
        gamma0 = -1.0;
        delta0 = alpha + 1.0;
        c1 = -alpha - 1.0;
        volume_1d = typeMethods.r8_gamma(alpha + 1.0);
        GolubWelsch_Xiu.gw_02_xiu(n, o, gamma0, delta0, c1, volume_1d, ref x, ref w);
        Console.WriteLine("");
        Console.WriteLine("  Generalized Laguerre:");
        Console.WriteLine("  ALPHA = " + alpha + "");
        Console.WriteLine("");
        for (j = 0; j < o; j++)
        {
            string cout = w[j].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            for (i = 0; i < n; i++)
            {
                cout += x[i + j * n].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);
        }

        //
        //  Hermite (physicist)
        //
        gamma0 = 2.0;
        delta0 = 0.0;
        c1 = 1.0;
        volume_1d = Math.Sqrt(pi);
        GolubWelsch_Xiu.gw_02_xiu(n, o, gamma0, delta0, c1, volume_1d, ref x, ref w);
        Console.WriteLine("");
        Console.WriteLine("  Hermite (physicist):");
        Console.WriteLine("");
        for (j = 0; j < o; j++)
        {
            string cout = w[j].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            for (i = 0; i < n; i++)
            {
                cout += x[i + j * n].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);
        }

        //
        //  Hermite (probabilist)
        //
        gamma0 = 1.0;
        delta0 = 0.0;
        c1 = 1.0;
        volume_1d = Math.Sqrt(2.0 * pi);
        GolubWelsch_Xiu.gw_02_xiu(n, o, gamma0, delta0, c1, volume_1d, ref x, ref w);
        Console.WriteLine("");
        Console.WriteLine("  Hermite ( probabilist):");
        Console.WriteLine("");
        for (j = 0; j < o; j++)
        {
            string cout = w[j].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            for (i = 0; i < n; i++)
            {
                cout += x[i + j * n].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);
        }

        //
        //  Jacobi.
        //
        alpha = 0.5;
        beta = 1.5;
        gamma0 = (alpha + beta + 2.0) / 2.0;
        delta0 = (alpha - beta) / 2.0;
        c1 = 2.0 * (alpha + 1.0) * (beta + 1.0) / (alpha + beta + 3.0)
                                                / (alpha + beta + 2.0);
        volume_1d = Math.Pow(2.0, alpha + beta + 1.0) * typeMethods.r8_gamma(alpha + 1)
                                                      * typeMethods.r8_gamma(beta + 1.0) / (alpha + beta + 1.0)
            / typeMethods.r8_gamma(alpha + beta + 1.0);
        GolubWelsch_Xiu.gw_02_xiu(n, o, gamma0, delta0, c1, volume_1d, ref x, ref w);
        Console.WriteLine("");
        Console.WriteLine("  Jacobi:");
        Console.WriteLine("  ALPHA = " + alpha + "");
        Console.WriteLine("  BETA  = " + beta + "");
        Console.WriteLine("");
        for (j = 0; j < o; j++)
        {
            string cout = w[j].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            for (i = 0; i < n; i++)
            {
                cout += x[i + j * n].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);
        }

        //
        //  Laguerre.
        //
        gamma0 = -1.0;
        delta0 = 1.0;
        c1 = -1.0;
        volume_1d = 1.0;
        GolubWelsch_Xiu.gw_02_xiu(n, o, gamma0, delta0, c1, volume_1d, ref x, ref w);
        Console.WriteLine("");
        Console.WriteLine("  Laguerre:");
        Console.WriteLine("");
        for (j = 0; j < o; j++)
        {
            string cout = w[j].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            for (i = 0; i < n; i++)
            {
                cout += x[i + j * n].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);
        }

        //
        //  Legendre.
        //
        gamma0 = 1.0;
        delta0 = 0.0;
        c1 = 1.0 / 3.0;
        volume_1d = 2.0;
        GolubWelsch_Xiu.gw_02_xiu(n, o, gamma0, delta0, c1, volume_1d, ref x, ref w);
        Console.WriteLine("");
        Console.WriteLine("  Legendre:");
        Console.WriteLine("");
        for (j = 0; j < o; j++)
        {
            string cout = w[j].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            for (i = 0; i < n; i++)
            {
                cout += x[i + j * n].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);
        }
    }
}