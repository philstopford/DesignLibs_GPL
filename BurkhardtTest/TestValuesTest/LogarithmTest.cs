using System;
using Burkardt.Values;

namespace TestValuesTest;

public static class LogarithmTest
{
    public static void log_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_VALUES_TEST tests LOG_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int n_data;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("LOG_VALUES_TEST:");
        Console.WriteLine("   LOG_VALUES stores values of the natural logarithm function.");
        Console.WriteLine("");
        Console.WriteLine("                X                     FX");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Log.log_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString().PadLeft(12) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void log_normal_cdf_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_NORMAL_CDF_VALUES_TEST tests LOG_NORMAL_CDF_VALUES.
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
        double mu = 0;
        int n_data;
        double sigma = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("LOG_NORMAL_CDF_VALUES_TEST:");
        Console.WriteLine("  LOG_NORMAL_CDF_VALUES returns values of ");
        Console.WriteLine("  the Log Normal Cumulative Density Function.");
        Console.WriteLine("");
        Console.WriteLine("     Mu      Sigma        X   CDF(X)");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Log.log_normal_cdf_values(ref n_data, ref mu, ref sigma, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + mu.ToString().PadLeft(8) + "  "
                              + sigma.ToString().PadLeft(8) + "  "
                              + x.ToString().PadLeft(8) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void log_series_cdf_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_SERIES_CDF_VALUES_TEST tests LOG_SERIES_CDF_VALUES.
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
        int n_data;
        double t = 0;
        Console.WriteLine("");
        Console.WriteLine("LOG_SERIES_CDF_VALUES_TEST:");
        Console.WriteLine("  LOG_SERIES_CDF_VALUES returns values of ");
        Console.WriteLine("  the Log Series Cumulative Density Function.");
        Console.WriteLine("");
        Console.WriteLine("     T      N   CDF(X)");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Log.log_series_cdf_values(ref n_data, ref t, ref n, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + t.ToString("0.################").PadLeft(24) + "  "
                              + n.ToString().PadLeft(6) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void log10_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG10_VALUES_TEST tests LOG10_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 January 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int n_data;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("LOG10_VALUES_TEST:");
        Console.WriteLine("   LOG10_VALUES stores values of the base 10 logarithm function.");
        Console.WriteLine("");
        Console.WriteLine("                X                     FX");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Log.log10_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString("0.################").PadLeft(24) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void logarithmic_integral_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOGARITHMIC_INTEGRAL_VALUES_TEST tests LOGARITHMIC_INTEGRAL_VALUES.
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
        int n_data;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("LOGARITHMIC_INTEGRAL_VALUES_TEST:");
        Console.WriteLine("  LOGARITHMIC_INTEGAL_VALUES stores values of");
        Console.WriteLine("  the logarithmic integral function.");
        Console.WriteLine("");
        Console.WriteLine("      X            LI(X)");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Log.logarithmic_integral_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString().PadLeft(12) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }
        
    public static void polylogarithm_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYLOGARITHM_VALUES_TEST tests POLYLOGARITHM_VALUES.
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
        int n = 0;
        int n_data;
        double z = 0;
        Console.WriteLine("");
        Console.WriteLine("POLYLOGARITHM_VALUES_TEST:");
        Console.WriteLine("  POLYLOGARITHM_VALUES returns values of ");
        Console.WriteLine("  the polylogarithm function.");
        Console.WriteLine("");
        Console.WriteLine("     N      Z          Fx");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Log.polylogarithm_values(ref n_data, ref n, ref z, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + n.ToString().PadLeft(6) + "  "
                              + z.ToString("0.################").PadLeft(24) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }
}