using System;
using System.Globalization;
using Burkardt.MonomialNS;
using Burkardt.Quadrature;
using Burkardt.Types;

namespace StroudTest;

using EpnLagIntegral = Burkardt.IntegralNS.Epn_Lag;
using EpnLagQuadrature = Epn_Lag;

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
        double err;
        int i;
        int o;
        double quad;
        double[] v;
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
        int d = typeMethods.i4vec_sum(n, expon);
        Console.WriteLine("  Degree = " + d + "");
        Console.WriteLine("");

        double exact = EpnLagIntegral.epn_glg_monomial_integral(n, expon, alpha);

        int p = 0;

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
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
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
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
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
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            o = GolubWelsch_Xiu.gw_02_xiu_size(n);
            double gamma0 = -1.0;
            double delta0 = alpha + 1.0;
            double c1 = -alpha - 1.0;
            double volume_1d = typeMethods.r8_gamma(1.0 + alpha);
            x = new double[n * o];
            w = new double[o];
            GolubWelsch_Xiu.gw_02_xiu(n, o, gamma0, delta0, c1, volume_1d, ref x, ref w);
            v = Monomial.monomial_value(n, o, x, expon);
            quad = typeMethods.r8vec_dot_product(o, w, v);
            err = Math.Abs(quad - exact);
            Console.WriteLine("  GW_02_XIU:      "
                              + "  " + o.ToString().PadLeft(6)
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        Console.WriteLine("  EXACT                   "
                          + "  " + exact.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

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
        double err;
        int i;
        int o;
        double quad;
        double[] v;
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
        int d = typeMethods.i4vec_sum(n, expon);
        Console.WriteLine("  Degree = " + d + "");
        Console.WriteLine("");

        double exact = EpnLagIntegral.epn_lag_monomial_integral(n, expon);

        int p = 0;

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
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
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
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
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
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            o = GolubWelsch_Xiu.gw_02_xiu_size(n);
            const double gamma0 = -1.0;
            const double delta0 = 1.0;
            const double c1 = -1.0;
            const double volume_1d = 1.0;
            x = new double[n * o];
            w = new double[o];
            GolubWelsch_Xiu.gw_02_xiu(n, o, gamma0, delta0, c1, volume_1d, ref x, ref w);
            v = Monomial.monomial_value(n, o, x, expon);
            quad = typeMethods.r8vec_dot_product(o, w, v);
            err = Math.Abs(quad - exact);
            Console.WriteLine("  GW_02_XIU:      "
                              + "  " + o.ToString().PadLeft(6)
                              + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        Console.WriteLine("  EXACT                   "
                          + "  " + exact.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }


}