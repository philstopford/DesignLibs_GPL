using System;
using Burkardt.PolynomialNS;

namespace LobattoPolynomialTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOBATTO_POLYNOMIAL_TEST tests the LOBATTO_POLYNOMIAL library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 November 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("LOBATTO_POLYNOMIAL_TEST:");
        Console.WriteLine("  Test the LOBATTO_POLYNOMIAL library.");

        lobatto_polynomial_value_test();
        lobatto_polynomial_derivative_test();
        lobatto_polynomial_plot_test();

        Console.WriteLine("");
        Console.WriteLine("LOBATTO_POLYNOMIAL_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void lobatto_polynomial_value_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOBATTO_POLYNOMIAL_VALUE_TEST tests LOBATTO_POLYNOMIAL_VALUE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 November 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double e;
        double fx1 = 0;
        double fx2;
        double[] l;
        int m;
        int n = 0;
        int n_data;
        double x = 0;
        double[] xvec = new double[1];

        m = 1;

        Console.WriteLine("");
        Console.WriteLine("LOBATTO_POLYNOMIAL_VALUE_TEST:");
        Console.WriteLine("  LOBATTO_POLYNOMIAL_VALUES stores values of");
        Console.WriteLine("  the completed Lobatto polynomial L(n,x).");
        Console.WriteLine("  LOBATTO_POLYNOMIAL_VALUE evaluates the completed Lobatto polynomial.");
        Console.WriteLine("");
        Console.WriteLine("                                       Tabulated                 Computed");
        Console.WriteLine("     N        X                        L(N,X)                    L(N,X)        Error");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Lobatto.lobatto_polynomial_values(ref n_data, ref n, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            xvec[0] = x;

            l = Lobatto.lobatto_polynomial_value(m, n, xvec);

            fx2 = l[0 + (n - 1) * m];

            e = fx1 - fx2;

            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + fx1.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + fx2.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + e.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");

        }
    }

    private static void lobatto_polynomial_derivative_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOBATTO_POLYNOMIAL_DERIVATIVE_TEST tests LOBATTO_POLYNOMIAL_DERIVATIVE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 November 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double e;
        double fx1 = 0;
        double fx2;
        double[] lp;
        int m;
        int n = 0;
        int n_data;
        double x = 0;
        double[] xvec = new double[1];

        m = 1;

        Console.WriteLine("");
        Console.WriteLine("LOBATTO_POLYNOMIAL_DERIVATIVE_TEST:");
        Console.WriteLine("  LOBATTO_POLYNOMIAL_DERIVATIVES stores derivatives of");
        Console.WriteLine("  the completed Lobatto polynomial L(n,x).");
        Console.WriteLine("  LOBATTO_POLYNOMIAL_DERIVATIVE evaluates the completed Lobatto polynomial.");
        Console.WriteLine("");
        Console.WriteLine("                                       Tabulated                 Computed");
        Console.WriteLine("     N        X                        L''(N,X)                   L''(N,X)       Error");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Lobatto.lobatto_polynomial_derivatives(ref n_data, ref n, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            xvec[0] = x;
            lp = Lobatto.lobatto_polynomial_derivative(m, n, xvec);
            fx2 = lp[0 + (n - 1) * m];

            e = fx1 - fx2;

            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + fx1.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + fx2.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + e.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");

        }

    }

    private static void lobatto_polynomial_plot_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOBATTO_POLYNOMIAL_PLOT_TEST tests LOBATTO_POLYNOMIAL_PLOT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 November 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] ndx = {1, 2, 3, 4, 5, 6, 7};
        int ndx_num = 7;
        string prefix = "test";

        Console.WriteLine("");
        Console.WriteLine("LOBATTO_POLYNOMIAL_PLOT_TEST:");
        Console.WriteLine("  LOBATTO_POLYNOMIAL_PLOT plots Lobatto polynomials.");

        Lobatto.lobatto_polynomial_plot(ndx_num, ndx, prefix);

    }
}