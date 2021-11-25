﻿using System;
using System.Globalization;
using Burkardt.Function;

namespace PolPakTest;

public static class zetaTest
{
    public static void zeta_m1_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZETA_M1_TEST tests ZETA_M1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 January 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double p = 0;
        double z1 = 0;

        const double tol = 1.0E-10;

        Console.WriteLine("");
        Console.WriteLine("ZETA_M1_TEST");
        Console.WriteLine("  ZETA_M1 computes the Zeta Minus One function.");
        Console.WriteLine("  Requested relative accuracy= " + tol + "");
        Console.WriteLine("");
        Console.WriteLine("       P            Zeta_M1(P)         Zeta_M1(P)");
        Console.WriteLine("                    tabulated          computed.");
        Console.WriteLine("");

        int n_data = 0;

        for (;;)
        {
            Burkardt.Values.Zeta.zeta_m1_values(ref n_data, ref p, ref z1);

            if (n_data == 0)
            {
                break;
            }

            double z2 = Zeta.zeta_m1(p, tol);

            Console.WriteLine("  " + p.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + z1.ToString(CultureInfo.InvariantCulture).PadLeft(20)
                                   + "  " + z2.ToString(CultureInfo.InvariantCulture).PadLeft(20) + "");

        }

    }

    public static void zeta_naive_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZETA_NAIVE_TEST tests ZETA_NAIVE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 March 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n = 0;
        double z1 = 0;

        Console.WriteLine("");
        Console.WriteLine("ZETA_NAIVE_TEST");
        Console.WriteLine("  ZETA_NAIVE computes the Zeta function.");
        Console.WriteLine("");
        Console.WriteLine("       N            exact Zeta         computed Zeta");
        Console.WriteLine("");

        int n_data = 0;

        for (;;)
        {
            Burkardt.Values.Zeta.zeta_values(ref n_data, ref n, ref z1);

            if (n_data == 0)
            {
                break;
            }

            double n_real = n;

            double z2 = Zeta.zeta_naive(n_real);

            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + z1.ToString(CultureInfo.InvariantCulture).PadLeft(20)
                                   + "  " + z2.ToString(CultureInfo.InvariantCulture).PadLeft(20) + "");
        }
    }
}