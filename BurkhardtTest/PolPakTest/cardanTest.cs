using System;
using Burkardt.PolynomialNS;
using Burkardt.Types;

namespace PolPakTest;

public static class cardanTest
{
    public static void cardan_poly_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CARDAN_POLY_TEST tests CARDAN_POLY.
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
        const int N_MAX = 10;

        double[] c = new double[N_MAX + 1];
        double cx1;
        double[] cx2;
        int n;
        double s;
        double x;

        Console.WriteLine("");
        Console.WriteLine("CARDAN_POLY_TEST");
        Console.WriteLine("  CARDAN_POLY evaluates a Cardan polynomial directly.");
        Console.WriteLine("");

        n = N_MAX;
        x = 0.25;
        s = 0.5;

        Console.WriteLine("");
        Console.WriteLine("  Compare CARDAN_POLY_COEF + R8POLY_VALUE_HORNER");
        Console.WriteLine("  versus CARDAN_POLY alone.");
        Console.WriteLine("");
        Console.WriteLine("  Evaluate polynomials at X = " + x + "");
        Console.WriteLine("  We use the parameter S = " + s + "");
        Console.WriteLine("");
        Console.WriteLine("  Order       Horner          Direct");
        Console.WriteLine("");

        cx2 = Cardan.cardan_poly(n, x, s);

        for (n = 0; n <= N_MAX; n++)
        {
            Cardan.cardan_poly_coef(n, s, ref c);

            cx1 = typeMethods.r8poly_value_horner(n, c, x);

            Console.WriteLine("  " + n.ToString().PadLeft(2)
                                   + "  " + cx1.ToString().PadLeft(14)
                                   + "  " + cx2[n].ToString().PadLeft(14) + "");
        }

    }

    public static void cardan_poly_coef_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CARDAN_POLY_COEF_TEST tests CARDAN_POLY_COEF.
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
        const int N_MAX = 10;

        double[] c = new double[N_MAX + 1];
        int i;
        int n;
        double s;

        s = 1.0;

        Console.WriteLine("");
        Console.WriteLine("CARDAN_POLY_COEF_TEST");
        Console.WriteLine("  CARDAN_POLY_COEF returns the coefficients of a");
        Console.WriteLine("  Cardan polynomial.");
        Console.WriteLine("");
        Console.WriteLine("  We use the parameter S = " + s + "");
        Console.WriteLine("");
        Console.WriteLine("  Table of polynomial coefficients:");
        Console.WriteLine("");

        for (n = 0; n <= N_MAX; n++)
        {
            Cardan.cardan_poly_coef(n, s, ref c);
            string cout = "  " + n.ToString().PadLeft(2) + "  ";
            for (i = 0; i <= n; i++)
            {
                cout += c[i].ToString().PadLeft(3) + "  ";
            }

            Console.WriteLine(cout);
        }

    }

}