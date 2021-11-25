using System;
using System.Globalization;
using Burkardt.Values;

namespace TestValuesTest;

public static class NormalTest
{
    public static void normal_01_cdf_values_test()
        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    NORMAL_01_CDF_VALUES_TEST tests NORMAL_01_CDF_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("NORMAL_01_CDF_VALUES_TEST:");
        Console.WriteLine("  NORMAL_01_CDF_VALUES stores values of");
        Console.WriteLine("  the Normal 01 Cumulative Density Function.");
        Console.WriteLine("");
        Console.WriteLine("            X                   CDF(X)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Normal.normal_01_cdf_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString("0.################").PadLeft(24) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void normal_01_pdf_values_test()
        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    NORMAL_01_PDF_VALUES_TEST tests NORMAL_01_PDF_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("NORMAL_01_PDF_VALUES_TEST:");
        Console.WriteLine("  NORMAL_01_PDF_VALUES stores values of");
        Console.WriteLine("  the Normal 01 Probability Density Function.");
        Console.WriteLine("");
        Console.WriteLine("            X                   PDF(X)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Normal.normal_01_pdf_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString("0.################").PadLeft(24) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void normal_cdf_values_test()
        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    NORMAL_CDF_VALUES_TEST tests NORMAL_CDF_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        double mu = 0;
        double sigma = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("NORMAL_CDF_VALUES_TEST:");
        Console.WriteLine("  NORMAL_CDF_VALUES stores values of");
        Console.WriteLine("  the Normal Cumulative Density Function.");
        Console.WriteLine("");
        Console.WriteLine("            X                   CDF(X)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Normal.normal_cdf_values(ref n_data, ref mu, ref sigma, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + mu.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + sigma.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void normal_pdf_values_test()
        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    NORMAL_PDF_VALUES_TEST tests NORMAL_PDF_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        double mu = 0;
        double sigma = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("NORMAL_PDF_VALUES_TEST:");
        Console.WriteLine("  NORMAL_PDF_VALUES stores values of");
        Console.WriteLine("  the Normal Probability Density Function.");
        Console.WriteLine("");
        Console.WriteLine("            X                   PDF(X)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Normal.normal_pdf_values(ref n_data, ref mu, ref sigma, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + mu.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + sigma.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }
}