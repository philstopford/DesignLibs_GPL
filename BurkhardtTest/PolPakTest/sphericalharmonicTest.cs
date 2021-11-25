﻿using System;
using System.Globalization;
using Burkardt.Function;

namespace PolPakTest;

public static class sphericalharmonicTest
{
    public static void spherical_harmonic_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERICAL_HARMONIC_TEST tests SPHERICAL_HARMONIC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N_MAX = 20;

        double[] c = new double[N_MAX + 1];
        int l = 0;
        int m = 0;
        double phi = 0;
        double[] s = new double[N_MAX + 1];
        double theta = 0;
        double yi = 0;
        double yr = 0;

        Console.WriteLine("");
        Console.WriteLine("SPHERICAL_HARMONIC_TEST:");
        Console.WriteLine("  SPHERICAL_HARMONIC evaluates spherical harmonic functions.");
        Console.WriteLine("");
        Console.WriteLine(
            "         N         M    THETA      PHI            YR            YI");
        Console.WriteLine("");

        int n_data = 0;

        for (;;)
        {
            Burkardt.Values.SphericalHarmonic.spherical_harmonic_values(ref n_data, ref l, ref m, ref theta,
                ref phi, ref yr, ref yi);

            if (n_data == 0)
            {
                break;
            }

            SphericalHarmonic.spherical_harmonic(l, m, theta, phi, ref c, ref s);

            double yr2 = c[l];
            double yi2 = s[l];

            Console.WriteLine("  "
                              + l.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + m.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + theta.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + phi.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + yr.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "  "
                              + yi.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            Console.WriteLine("  "
                              + "        " + "  "
                              + "        " + "  "
                              + "        " + "  "
                              + "        " + "  "
                              + yr2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "  "
                              + yi2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

    }

}