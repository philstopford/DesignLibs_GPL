﻿using System;
using System.Globalization;

namespace PolPakTest;

using TestValues = Burkardt.Values.Chebyshev;
using Polynomial = Burkardt.PolynomialNS.Chebyshev;

public static class chebyshevTest
{
    public static void cheby_t_poly_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBY_T_POLY_TEST tests CHEBY_T_POLY.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int n = 0;
        double x = 0;
        double[] x_vec = new double[1];

        Console.WriteLine("");
        Console.WriteLine("CHEBY_T_POLY_TEST:");
        Console.WriteLine("  CHEBY_T_POLY evaluates the Chebyshev T polynomial.");
        Console.WriteLine("");
        Console.WriteLine("     N      X        Exact F       T(N)(X)");
        Console.WriteLine("");

        int n_data = 0;

        for (;;)
        {
            TestValues.cheby_t_poly_values(ref n_data, ref n, ref x, ref fx);

            if (n_data == 0)
            {
                break;
            }

            x_vec[0] = x;
            double[] fx2 = Polynomial.cheby_t_poly(1, n, x_vec);

            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + fx2[n].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

        }

    }

    public static void cheby_t_poly_zero_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBY_T_POLY_ZERO_TEST tests CHEBY_T_POLY_ZERO.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 March 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N_MAX = 4;

        int n;

        Console.WriteLine("");
        Console.WriteLine("CHEBY_T_POLY_ZERO_TEST:");
        Console.WriteLine("  CHEBY_T_POLY_ZERO returns zeroes of T(N,X).");
        Console.WriteLine("");
        Console.WriteLine("       N      X        T(N,X)");
        Console.WriteLine("");

        for (n = 1; n <= N_MAX; n++)
        {
            double[] z = Polynomial.cheby_t_poly_zero(n);
            double[] fx = Polynomial.cheby_t_poly(n, n, z);
            int i;
            for (i = 0; i < n; i++)
            {
                Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + z[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + fx[i + n * n].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }

            Console.WriteLine("");
        }

    }

    public static void cheby_t_poly_coef_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBY_T_POLY_COEF_TEST tests CHEBY_T_POLY_COEF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        const int n = 5;

        Console.WriteLine("");
        Console.WriteLine("CHEBY_T_POLY_COEF_TEST");
        Console.WriteLine("  CHEBY_T_POLY_COEF determines the  polynomial coefficients");
        Console.WriteLine("  of the Chebyshev polynomial T(n,x).");

        double[] c = Polynomial.cheby_t_poly_coef(n);

        for (i = 0; i <= n; i++)
        {
            Console.WriteLine("");
            Console.WriteLine("  T(" + i + ",x)");
            Console.WriteLine("");
            int j;
            for (j = i; 0 <= j; j--)
            {
                if (c[i + j * (n + 1)] != 0.0)
                {
                    switch (j)
                    {
                        case 0:
                            Console.WriteLine(c[i + j * (n + 1)].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                            break;
                        case 1:
                            Console.WriteLine(c[i + j * (n + 1)].ToString(CultureInfo.InvariantCulture).PadLeft(14) + " * x");
                            break;
                        default:
                            Console.WriteLine(c[i + j * (n + 1)].ToString(CultureInfo.InvariantCulture).PadLeft(14) + " * x^" + j + "");
                            break;
                    }
                }
            }
        }


    }

    public static void cheby_u_poly_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBY_U_POLY_TEST tests CHEBY_U_POLY.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 January 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int n = 0;
        double x = 0;
        double[] x_vec = new double[1];

        Console.WriteLine("");
        Console.WriteLine("CHEBY_U_POLY_TEST:");
        Console.WriteLine("  CHEBY_U_POLY evaluates the Chebyshev U polynomial.");
        Console.WriteLine("");
        Console.WriteLine("     N      X        Exact F       U(N)(X)");
        Console.WriteLine("");

        int n_data = 0;

        for (;;)
        {
            TestValues.cheby_u_poly_values(ref n_data, ref n, ref x, ref fx);

            if (n_data == 0)
            {
                break;
            }

            x_vec[0] = x;
            double[] fx2 = Polynomial.cheby_u_poly(1, n, x_vec);

            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + fx2[n].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");


        }

    }

    public static void cheby_u_poly_coef_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBY_U_POLY_COEF_TEST tests CHEBY_U_POLY_COEF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 5;

        double[] c = new double[(N + 1) * (N + 1)];
        int i;

        Console.WriteLine("");
        Console.WriteLine("CHEBY_U_POLY_COEF_TEST");
        Console.WriteLine("  CHEBY_U_POLY_COEF determines the polynomial coefficients");
        Console.WriteLine("  of the Chebyshev polynomial U(n,x).");

        Polynomial.cheby_u_poly_coef(N, ref c);

        for (i = 0; i <= N; i++)
        {
            Console.WriteLine("");
            Console.WriteLine("  U(" + i + ",x)");
            Console.WriteLine("");
            int j;
            for (j = i; 0 <= j; j--)
            {
                if (c[i + j * (N + 1)] != 0.0)
                {
                    switch (j)
                    {
                        case 0:
                            Console.WriteLine(c[i + j * (N + 1)].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
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
        }

    }

    public static void cheby_u_poly_zero_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBY_U_POLY_ZERO_TEST tests CHEBY_U_POLY_ZERO.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 January 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N_MAX = 4;

        int n;

        Console.WriteLine("");
        Console.WriteLine("CHEBY_U_POLY_ZERO_TEST:");
        Console.WriteLine("  CHEBY_U_POLY_ZERO returns zeroes of U(N,X).");
        Console.WriteLine("");
        Console.WriteLine("       N      X        U(N,X)");
        Console.WriteLine("");

        for (n = 1; n <= N_MAX; n++)
        {
            double[] z = Polynomial.cheby_u_poly_zero(n);
            double[] fx = Polynomial.cheby_u_poly(n, n, z);
            int i;
            for (i = 0; i < n; i++)
            {
                Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + z[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + fx[i + n * n].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }

            Console.WriteLine("");
        }

    }

    public static void chebyshev_discrete_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBYSHEV_DISCRETE_TEST tests CHEBYSHEV_DISCRETE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 March 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 5;

        int j;
        double[] value = new double[N + 1];

        Console.WriteLine("");
        Console.WriteLine("CHEBYSHEV_DISCRETE_TEST:");
        Console.WriteLine("  CHEBYSHEV_DISCRETE evaluates discrete Chebyshev polynomials.");
        Console.WriteLine("");
        Console.WriteLine("       N      M         X        T(N,M,X)");

        const int m = 5;

        for (j = 0; j <= 5; j++)
        {
            double x = j / 2.0;

            Polynomial.chebyshev_discrete(N, m, x, ref value);

            Console.WriteLine("");
            int i;
            for (i = 0; i <= 5; i++)
            {
                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                       + "  " + m.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + value[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }
        }

    }

}