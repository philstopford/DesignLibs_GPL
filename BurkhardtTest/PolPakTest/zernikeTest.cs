using System;
using System.Globalization;
using Burkardt.PolynomialNS;
using Burkardt.Types;

namespace PolPakTest;

public static class zernikeTest
{
    public static void zernike_poly_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZERNIKE_POLY_TEST tests ZERNIKE_POLY.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 November 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] c;
        int m;
        int n;

        Console.WriteLine("");
        Console.WriteLine("ZERNIKE_POLY_TEST");
        Console.WriteLine("  ZERNIKE_POLY evaluates a Zernike polynomial directly.");
        Console.WriteLine("");
        Console.WriteLine("  Table of polynomial coefficients:");
        Console.WriteLine("");
        Console.WriteLine("   N   M");
        Console.WriteLine("");

        for (n = 0; n <= 5; n++)
        {
            Console.WriteLine("");
            for (m = 0; m <= n; m++)
            {
                c = Zernike.zernike_poly_coef(m, n);
                string cout = "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + m.ToString(CultureInfo.InvariantCulture).PadLeft(2);
                int i;
                for (i = 0; i <= n; i++)
                {
                    cout += c[i].ToString(CultureInfo.InvariantCulture).PadLeft(7);
                }

                Console.WriteLine(cout);
            }
        }

        double rho = 0.987654321;

        Console.WriteLine("");
        Console.WriteLine("  Z1: Compute polynomial coefficients,");
        Console.WriteLine("  then evaluate by Horner's method;");
        Console.WriteLine("  Z2: Evaluate directly by recursion.");
        Console.WriteLine("");
        Console.WriteLine("   N   M       Z1              Z2");
        Console.WriteLine("");

        for (n = 0; n <= 5; n++)
        {
            Console.WriteLine("");
            for (m = 0; m <= n; m++)
            {
                c = Zernike.zernike_poly_coef(m, n);
                double z1 = typeMethods.r8poly_value_horner(n, c, rho);

                double z2 = Zernike.zernike_poly(m, n, rho);
                Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                       + "  " + m.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                       + "  " + z1.ToString(CultureInfo.InvariantCulture).PadLeft(16)
                                       + "  " + z2.ToString(CultureInfo.InvariantCulture).PadLeft(16) + "");

            }
        }

    }

    public static void zernike_poly_coef_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZERNIKE_POLY_COEF_TEST tests ZERNIKE_POLY_COEF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 November 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int m;

        Console.WriteLine("");
        Console.WriteLine("ZERNIKE_POLY_COEF_TEST");
        Console.WriteLine("  ZERNIKE_POLY_COEF determines the Zernike");
        Console.WriteLine("  polynomial coefficients.");

        int n = 5;

        for (m = 0; m <= n; m++)
        {
            double[] c = Zernike.zernike_poly_coef(m, n);
            typeMethods.r8poly_print(n, c, "  Zernike polynomial");
        }

    }

}