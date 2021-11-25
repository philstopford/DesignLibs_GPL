using System;
using System.Globalization;
using Burkardt.PolynomialNS;

namespace PolPakTest;

using Symbol = Burkardt.Symbol.Jacobi;

public static class jacobiTest
{
    public static void jacobi_poly_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    JACOBI_POLY_TEST tests JACOBI_POLY.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double b = 0;
        double fx = 0;
        int n = 0;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("JACOBI_POLY_TEST");
        Console.WriteLine("  JACOBI_POLY evaluates the Jacobi polynomials.");
        Console.WriteLine("  the Jacobi polynomials.");
        Console.WriteLine("");
        Console.WriteLine("        N       A       B      X       JPV      JACOBI");
        Console.WriteLine("");

        int n_data = 0;

        for (;;)
        {

            Burkardt.Values.Jacobi.jacobi_poly_values(ref n_data, ref n, ref a, ref b, ref x, ref fx);

            if (n_data == 0)
            {
                break;
            }

            double[] c = Jacobi.jacobi_poly(n, a, b, x);
            double fx2 = c[n];

            Console.WriteLine("  "
                              + n.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + a.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + b.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                              + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "  "
                              + fx2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

        }

    }

    public static void jacobi_symbol_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    JACOBI_SYMBOL_TEST tests JACOBI_SYMBOL.
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
        const int N_TEST = 4;

        int i;
        int[] ptest = { 3, 9, 10, 12 };

        Console.WriteLine("");
        Console.WriteLine("JACOBI_SYMBOL_TEST");
        Console.WriteLine("  JACOBI_SYMBOL computes the Jacobi symbol");
        Console.WriteLine("  (Q/P), which records if Q is a quadratic");
        Console.WriteLine("  residue modulo the number P.");

        for (i = 0; i < N_TEST; i++)
        {
            int p = ptest[i];
            Console.WriteLine("");
            Console.WriteLine("Jacobi Symbols for P = " + p + "");
            Console.WriteLine("");
            int q;
            for (q = 0; q <= p; q++)
            {
                Console.WriteLine("  " + p.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + q.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + Symbol.jacobi_symbol(q, p).ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
            }
        }
    }
}