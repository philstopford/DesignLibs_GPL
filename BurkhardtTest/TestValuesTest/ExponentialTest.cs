using System;
using Burkardt.Values;

namespace TestValuesTest;

public class ExponentialTest
{
    public static void exp_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXP_VALUES_TEST tests EXP_VALUES.
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
        Console.WriteLine("EXP_VALUES_TEST:");
        Console.WriteLine("   EXP_VALUES stores values of the exponential function.");
        Console.WriteLine("");
        Console.WriteLine("                X                     FX");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Exponential.exp_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  " + x.ToString("0.################").PadLeft(24)
                                   + "  " + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void exp3_int_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXP3_INT_VALUES_TEST tests EXP3_INT_VALUES.
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
        int n_data;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("EXP3_INT_VALUES_TEST:");
        Console.WriteLine("  EXP3_INT_VALUES stores values of ");
        Console.WriteLine("  the exponential integral function.");
        Console.WriteLine("");
        Console.WriteLine("                X                     FX");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Exponential.exp3_int_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString("0.################").PadLeft(24) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void exponential_01_pdf_values_test()
        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    EXPONENTIAL_01_PDF_VALUES_TEST tests EXPONENTIAL_01_PDF_VALUES.
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
        int n_data;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("EXPONENTIAL_01_PDF_VALUES_TEST:");
        Console.WriteLine("  EXPONENTIAL_01_PDF_VALUES stores values of");
        Console.WriteLine("  the standard exponential Probability Density Function.");
        Console.WriteLine("");
        Console.WriteLine("            X                   PDF(X)");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Exponential.exponential_01_pdf_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString("0.################").PadLeft(24) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void exponential_cdf_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXPONENTIAL_CDF_VALUES_TEST tests EXPONENTIAL_CDF_VALUES.
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
        double lambda = 0;
        ;
        int n_data;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("EXPONENTIAL_CDF_VALUES_TEST:");
        Console.WriteLine("  EXPONENTIAL_CDF_VALUES stores values of ");
        Console.WriteLine("  the exponential CDF.");
        Console.WriteLine("");
        Console.WriteLine("       LAMBDA         X                     FX");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Exponential.exponential_cdf_values(ref n_data, ref lambda, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + lambda.ToString("0.########").PadLeft(24) + "  "
                              + x.ToString("0.########").PadLeft(24) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void exponential_pdf_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXPONENTIAL_PDF_VALUES_TEST tests EXPONENTIAL_PDF_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        double lambda = 0;
        ;
        int n_data;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("EXPONENTIAL_PDF_VALUES_TEST:");
        Console.WriteLine("  EXPONENTIAL_PDF_VALUES stores values of ");
        Console.WriteLine("  the exponential PDF.");
        Console.WriteLine("");
        Console.WriteLine("       LAMBDA         X                     FX");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Exponential.exponential_pdf_values(ref n_data, ref lambda, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + lambda.ToString("0.########").PadLeft(24) + lambda + "  "
                              + x.ToString("0.########").PadLeft(24) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }
}