using System;
using Burkardt.Values;

namespace TestValuesTest;

public static class FTest
{
    public static void f_cdf_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F_CDF_VALUES_TEST tests F_CDF_VALUES.
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
        int a = 0;
        int b = 0;
        double fx = 0;
        int n_data;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("F_CDF_VALUES_TEST:");
        Console.WriteLine("  F_CDF_VALUES stores values of");
        Console.WriteLine("  the F cumulative density function.");
        Console.WriteLine("");
        Console.WriteLine("     A       B            X            CDF(X)");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            F.f_cdf_values(ref n_data, ref a, ref b, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + a.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + b.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + fx.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    public static void f_noncentral_cdf_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F_NONCENTRAL_CDF_VALUES_TEST tests F_NONCENTRAL_CDF_VALUES.
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
        int a = 0;
        int b = 0;
        double fx = 0;
        double lambda = 0;
        ;
        int n_data;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("F_NONCENTRAL_CDF_VALUES_TEST:");
        Console.WriteLine("  F_NONCENTRAL_CDF_VALUES stores values of");
        Console.WriteLine("  the F cumulative density function.");
        Console.WriteLine("");
        Console.WriteLine("     A       B            LAMBDA    X            CDF");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            F.f_noncentral_cdf_values(ref n_data, ref a, ref b, ref lambda, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + a.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + b.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + lambda.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + fx.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

}