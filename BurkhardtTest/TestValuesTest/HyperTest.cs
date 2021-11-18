using System;
using Burkardt.Values;

namespace TestValuesTest;

public static class HyperTest
{
    public static void hyper_1f1_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPER_1F1_VALUES_TEST tests HYPER_1F1_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 January 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double b = 0;
        double fx = 0;
        int n_data;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("HYPER_1F1_VALUES_TEST:");
        Console.WriteLine("  HYPER_1F1_VALUES stores values of");
        Console.WriteLine("  the hypergeometric function 1F1.");
        Console.WriteLine("");
        Console.WriteLine("      A      B      X   Hyper_1F1(A,B,X)");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Hypergeometric.hyper_1f1_values(ref n_data, ref a, ref b, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + a.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + b.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void hyper_2f1_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPER_2F1_VALUES_TEST tests HYPER_2F1_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 September 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double b = 0;
        double c = 0;
        double fx = 0;
        int n_data;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("HYPER_2F1_VALUES_TEST:");
        Console.WriteLine("  HYPER_2F1_VALUES stores values of");
        Console.WriteLine("  the hypergeometric function 2F1.");
        Console.WriteLine("");
        Console.WriteLine("      A      B     C      X   Hyper_2F1(A,B,C,X)");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Hypergeometric.hyper_2f1_values(ref n_data, ref a, ref b, ref c, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + a.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + b.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + c.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void hypergeometric_cdf_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERGEOMETRIC_CDF_VALUES_TEST tests HYPERGEOMETRIC_CDF_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int n_data;
        int pop = 0;
        int sam = 0;
        int succ = 0;
        int x = 0;
        Console.WriteLine("");
        Console.WriteLine("HYPERGEOMETRIC_CDF_VALUES_TEST:");
        Console.WriteLine("  HYPERGEOMETRIC_CDF_VALUES stores values of");
        Console.WriteLine("  the Hypergeometric CDF.");
        Console.WriteLine("");
        Console.WriteLine("     SAM    SUC   POP     X   HyperCDF(S,S,P)(X)");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Hypergeometric.hypergeometric_cdf_values(ref n_data, ref sam, ref succ, ref pop, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + sam.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + succ.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + pop.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void hypergeometric_pdf_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERGEOMETRIC_PDF_VALUES_TEST tests HYPERGEOMETRIC_PDF_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 January 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int n_data;
        int pop = 0;
        int sam = 0;
        int succ = 0;
        int x = 0;
        Console.WriteLine("");
        Console.WriteLine("HYPERGEOMETRIC_PDF_VALUES_TEST:");
        Console.WriteLine("  HYPERGEOMETRIC_PDF_VALUES stores values of");
        Console.WriteLine("  the Hypergeometric PDF.");
        Console.WriteLine("");
        Console.WriteLine("     SAM    SUC   POP     X   HyperPDF(S,S,P)(X)");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Hypergeometric.hypergeometric_pdf_values(ref n_data, ref sam, ref succ, ref pop, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + sam.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + succ.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + pop.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void hypergeometric_u_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERGEOMETRIC_U_VALUES_TEST tests HYPERGEOMETRIC_U_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 January 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double b = 0;
        double fx = 0;
        int n_data;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("HYPERGEOMETRIC_U_VALUES_TEST:");
        Console.WriteLine("  HYPERGEOMETRIC_U_VALUES stores values of");
        Console.WriteLine("  the hypergeometric function U.");
        Console.WriteLine("");
        Console.WriteLine("      A      B      X   HyperU(A,B,X)");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Hypergeometric.hypergeometric_u_values(ref n_data, ref a, ref b, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + a.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + b.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }
}