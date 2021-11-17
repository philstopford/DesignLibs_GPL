using System;

namespace PolPakTest;

using Polynomial = Burkardt.PolynomialNS.Legendre;
using Function = Burkardt.Function.Legendre;
using TestValues = Burkardt.Values.Legendre;
using Symbol = Burkardt.Symbol.Legendre;

public static class legendreTest
{
    public static void legendre_poly_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGENDRE_POLY_TEST tests LEGENDRE_POLY.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N_MAX = 12;

        double fx = 0;
        double[] fp2 = new double[N_MAX + 1];
        double[] fx2 = new double[N_MAX + 1];
        int n = 0;
        int n_data = 0;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("LEGENDRE_POLY_TEST:");
        Console.WriteLine("  LEGENDRE_POLY evaluates the Legendre PN function.");
        Console.WriteLine("");
        Console.WriteLine("     N      X        Exact F       P(N)(X)");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            TestValues.legendre_polynomial_values(ref n_data, ref n, ref x, ref fx);

            if (n_data == 0)
            {
                break;
            }

            Polynomial.legendre_poly(n, x, ref fx2, ref fp2);

            Console.WriteLine("  "
                              + n.ToString().PadLeft(8) + "  "
                              + x.ToString().PadLeft(8) + "  "
                              + fx.ToString().PadLeft(14) + "  "
                              + fx2[n].ToString().PadLeft(14) + "");

        }

    }

    public static void legendre_poly_coef_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGENDRE_POLY_COEF_TEST tests LEGENDRE_POLY_COEF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 5;

        double[] c = new double[(N + 1) * (N + 1)];
        int i;
        int j;

        Console.WriteLine("");
        Console.WriteLine("LEGENDRE_POLY_COEF_TEST");
        Console.WriteLine("  LEGENDRE_POLY_COEF determines the Legendre P ");
        Console.WriteLine("  polynomial coefficients.");

        Polynomial.legendre_poly_coef(N, ref c);

        for (i = 0; i <= N; i++)
        {
            Console.WriteLine("");
            Console.WriteLine("  P(" + i + ")");
            Console.WriteLine("");
            for (j = i; 0 <= j; j--)
            {
                switch (j)
                {
                    case 0:
                        Console.WriteLine(c[i + j * (N + 1)].ToString().PadLeft(14) + "");
                        ;
                        break;
                    case 1:
                        Console.WriteLine(c[i + j * (N + 1)].ToString().PadLeft(14) + " * x");
                        break;
                    default:
                        Console.WriteLine(c[i + j * (N + 1)].ToString().PadLeft(14) + " * x^" + j + "");
                        break;
                }
            }
        }

    }

    public static void legendre_associated_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGENDRE_ASSOCIATED_TEST tests LEGENDRE_ASSOCIATED.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N_MAX = 20;

        double[] fx2 = new double[N_MAX + 1];
        double fx = 0;
        int m = 0;
        int n = 0;
        int n_data = 0;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("LEGENDRE_ASSOCIATED_TEST:");
        Console.WriteLine("  LEGENDRE_ASSOCIATED evaluates associated Legendre functions.");
        Console.WriteLine("");
        Console.WriteLine("      N       M    X     Exact F     PNM(X)");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Burkardt.Values.Legendre.legendre_associated_values(ref n_data, ref n, ref m, ref x, ref fx);

            if (n_data == 0)
            {
                break;
            }

            Polynomial.legendre_associated(n, m, x, ref fx2);

            Console.WriteLine("  "
                              + n.ToString().PadLeft(8) + "  "
                              + m.ToString().PadLeft(8) + "  "
                              + x.ToString().PadLeft(8) + "  "
                              + fx.ToString().PadLeft(14) + "  "
                              + fx2[n].ToString().PadLeft(14) + "");

        }

    }

    public static void legendre_associated_normalized_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGENDRE_ASSOCIATED_NORMALIZED_TEST tests LEGENDRE_ASSOCIATED_NORMALIZED.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 February 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N_MAX = 20;

        double[] fx2 = new double[N_MAX + 1];
        double fx = 0;
        int m = 0;
        int n = 0;
        int n_data = 0;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("LEGENDRE_ASSOCIATED_NORMALIZED_TEST:");
        Console.WriteLine("  LEGENDRE_ASSOCIATED_NORMALIZED evaluates ");
        Console.WriteLine("  normalized associated Legendre functions.");
        Console.WriteLine("");
        Console.WriteLine("      N       M    X     Exact F     PNM(X)");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Burkardt.Values.Legendre.legendre_associated_normalized_sphere_values(ref n_data, ref n, ref m,
                ref x, ref fx);

            if (n_data == 0)
            {
                break;
            }

            Polynomial.legendre_associated_normalized(n, m, x, ref fx2);

            Console.WriteLine("  "
                              + n.ToString().PadLeft(8) + "  "
                              + m.ToString().PadLeft(8) + "  "
                              + x.ToString().PadLeft(8) + "  "
                              + fx.ToString().PadLeft(14) + "  "
                              + fx2[n].ToString().PadLeft(14) + "");

        }

    }

    public static void legendre_function_q_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGENDRE_FUNCTION_Q_TEST tests LEGENDRE_FUNCTION_Q.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N_MAX = 12;

        double fx = 0;
        double[] fx2 = new double[N_MAX + 1];
        int n = 0;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("LEGENDRE_FUNCTION_Q_TEST:");
        Console.WriteLine("  LEGENDRE_FUNCTION_Q evaluates the Legendre Q function.");
        Console.WriteLine("");
        Console.WriteLine("     N      X        Exact F       Q(N)(X)");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            TestValues.legendre_function_q_values(ref n_data, ref n, ref x, ref fx);

            if (n_data == 0)
            {
                break;
            }

            Function.legendre_function_q(n, x, ref fx2);

            Console.WriteLine("  "
                              + n.ToString().PadLeft(8) + "  "
                              + x.ToString().PadLeft(8) + "  "
                              + fx.ToString().PadLeft(14) + "  "
                              + fx2[n].ToString().PadLeft(14) + "");

        }

    }

    public static void legendre_symbol_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGENDRE_SYMBOL_TEST tests LEGENDRE_SYMBOL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N_TEST = 4;

        int i;
        int p;
        int[] ptest = { 7, 11, 13, 17 };
        int q;

        Console.WriteLine("");
        Console.WriteLine("LEGENDRE_SYMBOL_TEST");
        Console.WriteLine("  LEGENDRE_SYMBOL computes the Legendre");
        Console.WriteLine("  symbol (Q/P) which records whether Q is ");
        Console.WriteLine("  a quadratic residue modulo the prime P.");

        for (i = 0; i < N_TEST; i++)
        {
            p = ptest[i];
            Console.WriteLine("");
            Console.WriteLine("  Legendre Symbols for P = " + p + "");
            Console.WriteLine("");
            for (q = 0; q <= p; q++)
            {
                Console.WriteLine("  "
                                  + p.ToString().PadLeft(8) + "  "
                                  + q.ToString().PadLeft(8) + "  "
                                  + Symbol.legendre_symbol(q, p).ToString().PadLeft(8) + "");
            }
        }
    }

}