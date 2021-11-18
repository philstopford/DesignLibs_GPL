using System;
using Burkardt.PolynomialNS;
using Burkardt.Types;

namespace PolPakTest;

public static class laguerreTest
{
    public static void gen_laguerre_poly_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEN_LAGUERRE_POLY_TEST tests GEN_LAGUERRE_POLY.
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
        int N = 10;
        int N_TEST = 6;

        double alpha;
        double[] alpha_test = { 0.0, 0.0, 0.1, 0.1, 0.5, 1.0 };
        double[] c = new double [N + 1];
        int i;
        int j;
        double x;
        double[] x_test = { 0.0, 1.0, 0.0, 0.5, 0.5, 0.5 };

        Console.WriteLine("");
        Console.WriteLine("GEN_LAGUERRE_POLY_TEST");
        Console.WriteLine("  GEN_LAGUERRE_POLY evaluates the generalized Laguerre");
        Console.WriteLine("  functions.");

        for (i = 0; i < N_TEST; i++)
        {

            x = x_test[i];
            alpha = alpha_test[i];

            Console.WriteLine("");
            Console.WriteLine("  Table of L(N,ALPHA)(X) for");
            Console.WriteLine("");
            Console.WriteLine("    N(max) = " + N + "");
            Console.WriteLine("    ALPHA =  " + alpha + "");
            Console.WriteLine("    X =      " + x + "");
            Console.WriteLine("");

            Laguerre.gen_laguerre_poly(N, alpha, x, ref c);

            for (j = 0; j <= N; j++)
            {
                Console.WriteLine("  "
                                  + j.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                  + c[j].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }
        }
    }

    public static void laguerre_associated_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGUERRE_ASSOCIATED_TEST tests LAGUERRE_ASSOCIATED.
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
        int N = 6;
        int N_TEST = 6;

        double[] c = new double[N + 1];
        int i;
        int j;
        int m;
        int[] m_test =  {
                0, 0, 1, 2, 3, 1
            }
            ;
        double x;
        double[] x_test =  {
                0.0, 1.0, 0.0, 0.5, 0.5, 0.5
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("LAGUERRE_ASSOCIATED_TEST");
        Console.WriteLine("  LAGUERRE_ASSOCIATED evaluates the associated Laguerre");
        Console.WriteLine("  polynomials.");

        for (i = 0; i < N_TEST; i++)
        {
            m = m_test[i];
            x = x_test[i];

            Console.WriteLine("");
            Console.WriteLine("  Table of L(N,M)(X) for");
            Console.WriteLine("");
            Console.WriteLine("  N(max) = " + N + "");
            Console.WriteLine("  M      = " + m + "");
            Console.WriteLine("  X =      " + x + "");
            Console.WriteLine("");

            Laguerre.laguerre_associated(N, m, x, ref c);

            for (j = 0; j <= N; j++)
            {
                Console.WriteLine("  "
                                  + j.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                  + c[j].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }

        }
    }

    public static void laguerre_poly_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGUERRE_POLY_TEST tests LAGUERRE_POLY.
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
        int n_data = 0;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("LAGUERRE_POLY_COEF:");
        Console.WriteLine("  LAGUERRE_POLY evaluates the Laguerre polynomial.");
        Console.WriteLine("");
        Console.WriteLine("     N      X        Exact F       L(N)(X)");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Burkardt.Values.Laguerre.laguerre_polynomial_values(ref n_data, ref n, ref x, ref fx);

            if (n_data == 0)
            {
                break;
            }

            Laguerre.laguerre_poly(n, x, ref fx2);

            Console.WriteLine("  "
                              + n.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "  "
                              + fx2[n].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

    }

    public static void laguerre_poly_coef_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST055 tests LAGUERRE_POLY_COEF.
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
        double fact;
        int i;
        int j;

        Console.WriteLine("");
        Console.WriteLine("LAGUERRE_POLY_COEF_TEST");
        Console.WriteLine("  LAGUERRE_POLY_COEF determines Laguerre ");
        Console.WriteLine("  polynomial coefficients.");

        Laguerre.laguerre_poly_coef(N, ref c);

        for (i = 0; i <= N; i++)
        {
            Console.WriteLine("");
            Console.WriteLine("  L(" + i + ")");
            Console.WriteLine("");
            for (j = i; 0 <= j; j--)
            {
                switch (j)
                {
                    case 0:
                        Console.WriteLine(c[i + j * (N + 1)].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                        ;
                        break;
                    case 1:
                        Console.WriteLine(c[i + j * (N + 1)].ToString(CultureInfo.InvariantCulture).PadLeft(14) + " * x");
                        break;
                    default:
                        Console.WriteLine(c[i + j * (N + 1)].ToString(CultureInfo.InvariantCulture).PadLeft(14) + " * x^" + j + "");
                        break;
                }
            }
        }

        for (i = 0; i <= N; i++)
        {
            fact = typeMethods.r8_factorial(i);
            Console.WriteLine("");
            Console.WriteLine("  Factorially scaled L(" + i + ")");
            Console.WriteLine("");
            for (j = i; 0 <= j; j--)
            {
                switch (j)
                {
                    case 0:
                        Console.WriteLine((fact * c[i + j * (N + 1)]).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                        ;
                        break;
                    case 1:
                        Console.WriteLine((fact * c[i + j * (N + 1)]).ToString(CultureInfo.InvariantCulture).PadLeft(14) + " * x");
                        break;
                    default:
                        Console.WriteLine((fact * c[i + j * (N + 1)]).ToString(CultureInfo.InvariantCulture).PadLeft(14) + " * x^" + j + "");
                        break;
                }
            }
        }
    }

}