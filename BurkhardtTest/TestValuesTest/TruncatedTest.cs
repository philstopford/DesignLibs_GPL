using System;
using System.Globalization;
using Burkardt.Values;

namespace TestValuesTest;

public static class TruncatedTest
{
    public static void truncated_normal_ab_cdf_test()
        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    TRUNCATED_NORMAL_AB_CDF_VALUES_TEST tests TRUNCATED_NORMAL_AB_CDF_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 September 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double b = 0;
        double fx = 0;
        double mu = 0;
        double sigma = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("TRUNCATED_NORMAL_AB_CDF_VALUES_TEST:");
        Console.WriteLine("  TRUNCATED_NORMAL_AB_CDF_VALUES stores values of");
        Console.WriteLine("  the Truncated Normal Cumulative Density Function.");
        Console.WriteLine("");
        Console.WriteLine("        MU     SIGMA       A         B         X        CDF(X)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Truncated.truncated_normal_ab_cdf_values(ref n_data, ref mu, ref sigma, ref a, ref b, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  " + mu.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + sigma.ToString(CultureInfo.InvariantCulture).PadLeft(8) + sigma
                                   + "  " + a.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + b.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void truncated_normal_ab_pdf_test()
        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    TRUNCATED_NORMAL_AB_PDF_VALUES_TEST tests TRUNCATED_NORMAL_AB_PDF_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 September 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double b = 0;
        double fx = 0;
        double mu = 0;
        double sigma = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("TEST1547:");
        Console.WriteLine("  TRUNCATED_NORMAL_AB_PDF_VALUES stores values of");
        Console.WriteLine("  the Truncated Normal Probability Density Function.");
        Console.WriteLine("");
        Console.WriteLine("        MU     SIGMA       A         B         X        PDF(X)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Truncated.truncated_normal_ab_pdf_values(ref n_data, ref mu, ref sigma, ref a, ref b, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  " + mu.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + sigma.ToString(CultureInfo.InvariantCulture).PadLeft(8) + sigma
                                   + "  " + a.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + b.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void truncated_normal_a_cdf_test()
        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    TRUNCATED_NORMAL_A_CDF_VALUES_TEST tests TRUNCATED_NORMAL_A_CDF_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 September 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double fx = 0;
        double mu = 0;
        double sigma = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("TEST1548:");
        Console.WriteLine("  TRUNCATED_NORMAL_A_CDF_VALUES stores values of");
        Console.WriteLine("  the lower Truncated Normal Cumulative Density Function.");
        Console.WriteLine("");
        Console.WriteLine("        MU     SIGMA       A         X        CDF(X)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Truncated.truncated_normal_a_cdf_values(ref n_data, ref mu, ref sigma, ref a, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  " + mu.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + sigma.ToString(CultureInfo.InvariantCulture).PadLeft(8) + sigma
                                   + "  " + a.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void truncated_normal_a_pdf_test()
        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    TRUNCATED_NORMAL_A_PDF_VALUES_TEST tests TRUNCATED_NORMAL_A_PDF_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 September 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double fx = 0;
        double mu = 0;
        double sigma = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("TRUNCATED_NORMAL_A_PDF_VALUES_TEST:");
        Console.WriteLine("  TRUNCATED_NORMAL_A_PDF_VALUES stores values of");
        Console.WriteLine("  the lower Truncated Normal Probability Density Function.");
        Console.WriteLine("");
        Console.WriteLine("        MU     SIGMA       A         X        PDF(X)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Truncated.truncated_normal_a_pdf_values(ref n_data, ref mu, ref sigma, ref a, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  " + mu.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + sigma.ToString(CultureInfo.InvariantCulture).PadLeft(8) + sigma
                                   + "  " + a.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void truncated_normal_b_cdf_test()
        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    TRUNCATED_NORMAL_B_CDF_VALUES_TEST tests TRUNCATED_NORMAL_B_CDF_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 September 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double b = 0;
        double fx = 0;
        double mu = 0;
        double sigma = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("TRUNCATED_NORMAL_B_CDF_VALUES_TEST:");
        Console.WriteLine("  TRUNCATED_NORMAL_B_CDF_VALUES stores values of");
        Console.WriteLine("  the upper Truncated Normal Cumulative Density Function.");
        Console.WriteLine("");
        Console.WriteLine("        MU     SIGMA       B         X        CDF(X)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Truncated.truncated_normal_b_cdf_values(ref n_data, ref mu, ref sigma, ref b, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  " + mu.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + sigma.ToString(CultureInfo.InvariantCulture).PadLeft(8) + sigma
                                   + "  " + b.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void truncated_normal_b_pdf_test()
        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    TRUNCATED_NORMAL_B_PDF_VALUES_TEST tests TRUNCATED_NORMAL_B_PDF_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 September 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double b = 0;
        double fx = 0;
        double mu = 0;
        double sigma = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("TRUNCATED_NORMAL_B_PDF_VALUES_TEST:");
        Console.WriteLine("  TRUNCATED_NORMAL_B_PDF_VALUES stores values of");
        Console.WriteLine("  the upper Truncated Normal Probability Density Function.");
        Console.WriteLine("");
        Console.WriteLine("        MU     SIGMA       B         X        PDF(X)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Truncated.truncated_normal_b_pdf_values(ref n_data, ref mu, ref sigma, ref b, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  " + mu.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + sigma.ToString(CultureInfo.InvariantCulture).PadLeft(8) + sigma
                                   + "  " + b.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + fx.ToString("0.################").PadLeft(24) + "");
        }
    }
}