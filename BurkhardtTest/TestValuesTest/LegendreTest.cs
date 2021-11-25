using System;
using System.Globalization;
using Burkardt.Values;

namespace TestValuesTest;

public static class LegendreTest
{
    public static void legendre_associated_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGENDRE_ASSOCIATED_VALUES_TEST tests LEGENDRE_ASSOCIATED_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int m = 0;
        int n = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("LEGENDRE_ASSOCIATED_VALUES_TEST:");
        Console.WriteLine("  LEGENDRE_ASSOCIATED_VALUES stores values of");
        Console.WriteLine("  the associated Legendre polynomials.");
        Console.WriteLine("");
        Console.WriteLine("     N     M    X             P(N,M)(X)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Legendre.legendre_associated_values(ref n_data, ref n, ref m, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + n.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + m.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void legendre_associated_normalized_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGENDRE_ASSOCIATED_NORMALIZED_VALUES_TEST tests LEGENDRE_ASSOCIATED_NORMALIZED_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int m = 0;
        int n = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("LEGENDRE_ASSOCIATED_NORMALIZED_VALUES_TEST:");
        Console.WriteLine("  LEGENDRE_ASSOCIATED_NORMALIZED_VALUES stores values of");
        Console.WriteLine("  the normalized associated Legendre polynomials.");
        Console.WriteLine("");
        Console.WriteLine("     N     M    X             P(N,M)(X)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Legendre.legendre_associated_normalized_values(ref n_data, ref n, ref m, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + n.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + m.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void legendre_associated_normalized_sphere_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGENDRE_ASSOCIATED_NORMALIZED_SPHERE_VALUES_TEST tests LEGENDRE_ASSOCIATED_NORMALIZED_SPHERE_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 March 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int m = 0;
        int n = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("LEGENDRE_ASSOCIATED_NORMALIZED_SPHERE_VALUES_TEST:");
        Console.WriteLine("  LEGENDRE_ASSOCIATED_NORMALIZED_SPHERE_VALUES stores values of");
        Console.WriteLine("  the associated Legendre polynomials, ref normalized for the unit sphere.");
        Console.WriteLine("");
        Console.WriteLine("     N     M    X             P(N,M)(X)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Legendre.legendre_associated_normalized_sphere_values(ref n_data, ref n, ref m, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + n.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + m.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void legendre_normalized_polynomial_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGENDRE_NORMALIZED_POLYNOMIAL_VALUES_TEST tests LEGENDRE_NORMALIZED_POLYNOMIAL_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 March 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int n = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("LEGENDRE_NORMALIZED_POLYNOMIAL_VALUES_TEST:");
        Console.WriteLine("  LEGENDRE_NORMALIZED_POLYNOMIAL_VALUES stores values of ");
        Console.WriteLine("  the normalized Legendre polynomials.");
        Console.WriteLine("");
        Console.WriteLine("     N    X             Pn(N)(X)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Legendre.legendre_normalized_polynomial_values(ref n_data, ref n, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + n.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void legendre_polynomial_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGENDRE_POLYNOMIAL_VALUES_TEST tests LEGENDRE_POLYNOMIAL_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 March 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int n = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("LEGENDRE_POLYNOMIAL_VALUES_TEST:");
        Console.WriteLine("  LEGENDRE_POLYNOMIAL_VALUES stores values of ");
        Console.WriteLine("  the Legendre polynomials.");
        Console.WriteLine("");
        Console.WriteLine("     N    X             P(N)(X)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Legendre.legendre_polynomial_values(ref n_data, ref n, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + n.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void legendre_shifted_polynomial_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGENDRE_SHIFTED_POLYNOMIAL_VALUES_TEST tests LEGENDRE_SHIFTED_POLYNOMIAL_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 March 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int n = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("LEGENDRE_SHIFTED_POLYNOMIAL_VALUES_TEST:");
        Console.WriteLine("  LEGENDRE_SHIFTED_POLYNOMIAL_VALUES stores values of ");
        Console.WriteLine("  the shifted Legendre polynomials.");
        Console.WriteLine("");
        Console.WriteLine("     N    X             P(N)(X)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Legendre.legendre_shifted_polynomial_values(ref n_data, ref n, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + n.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void legendre_function_q_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGENDRE_FUNCTION_Q_VALUES_TEST tests LEGENDRE_FUNCTION_Q_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int n = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("LEGENDRE_FUNCTION_Q_VALUES_TEST:");
        Console.WriteLine("  LEGENDRE_FUNCTION_Q_VALUES stores values of");
        Console.WriteLine("  the Legendre Q function.");
        Console.WriteLine("");
        Console.WriteLine("     N    X             Q(N)(X)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Legendre.legendre_function_q_values(ref n_data, ref n, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + n.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }
}