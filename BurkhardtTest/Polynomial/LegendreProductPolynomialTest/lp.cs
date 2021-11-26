﻿using System;
using System.Globalization;
using Burkardt.PolynomialNS;

namespace LegendreProductPolynomialTest;

public static class lpTest
{
    public static void lp_coefficients_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LP_COEFFICIENTS_TEST tests LP_COEFFICIENTS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int m = 1;
        int n;
        const int N_MAX = 10;
        int o = 0;

        Console.WriteLine("");
        Console.WriteLine("LP_COEFFICIENTS_TEST");
        Console.WriteLine("  LP_COEFFICIENTS: coefficients of Legendre polynomial P(n,x).");
        Console.WriteLine("");

        for (n = 0; n <= N_MAX; n++)
        {
            double[] c = new double[n + 1];
            int[] f = new int[n + 1];

            Legendre.lp_coefficients(n, ref o, ref c, ref f);

            int[] e = new int[o];
            int i;
            for (i = 0; i < o; i++)
            {
                e[i] = f[i] + 1;
            }

            string label = "  P(" + n + ",x) = ";
            Polynomial.polynomial_print(m, o, c, e, label);
        }
    }

    public static void lp_value_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LP_VALUE_TEST tests LP_VALUE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int o = 0;
        double x = 0;
        double[] xvec = new double[1];
        double fx1 = 0;

        const int n = 1;

        Console.WriteLine("");
        Console.WriteLine("LP_VALUE_TEST:");
        Console.WriteLine("  LP_VALUE evaluates a Legendre polynomial.");
        Console.WriteLine("");
        Console.WriteLine("                        Tabulated                 Computed");
        Console.WriteLine("     O        X           L(O,X)                    L(O,X)" +
                          "                   Error");
        Console.WriteLine("");

        int n_data = 0;

        for (;;)
        {
            Burkardt.Values.Legendre.lp_values(ref n_data, ref o, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            xvec[0] = x;

            double[] fx2 = Legendre.lp_value(n, o, xvec);

            double e = fx1 - fx2[0];

            Console.WriteLine(o.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                                      + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                      + fx1.ToString(CultureInfo.InvariantCulture).PadLeft(24) + "  "
                                                      + fx2[0].ToString(CultureInfo.InvariantCulture).PadLeft(24) + "  "
                                                      + e.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
        }

    }

    public static void lp_values_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LP_VALUES_TEST tests LP_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int o = 0;
        double x = 0;
        double fx = 0;

        Console.WriteLine("");
        Console.WriteLine("LP_VALUES_TEST:");
        Console.WriteLine("  LP_VALUES stores values of");
        Console.WriteLine("  the Legendre polynomial P(o,x).");
        Console.WriteLine("");
        Console.WriteLine("                        Tabulated");
        Console.WriteLine("     O        X           L(O,X)");
        Console.WriteLine("");

        int n_data = 0;

        for (;;)
        {
            Burkardt.Values.Legendre.lp_values(ref n_data, ref o, ref x, ref fx);

            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine(o.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                                      + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                      + fx.ToString(CultureInfo.InvariantCulture).PadLeft(24) + "");
        }
    }
}