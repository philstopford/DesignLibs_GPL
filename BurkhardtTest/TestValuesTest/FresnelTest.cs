using System;
using Burkardt.Values;

namespace TestValuesTest;

public static class FresnelTest
{
    public static void fresnel_cos_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FRESNEL_COS_VALUES_TEST tests FRESNEL_COS_VALUES.
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
        Console.WriteLine("FRESNEL_COS_VALUES_TEST:");
        Console.WriteLine("  FRESNEL_COS_VALUES stores values of");
        Console.WriteLine("  the Fresnel cosine integral C(X).");
        Console.WriteLine("");
        Console.WriteLine("      X           C(X)");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Fresnel.fresnel_cos_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString().PadLeft(12) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void fresnel_sin_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FRESNEL_SIN_VALUES_TEST tests FRESNEL_SIN_VALUES.
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
        Console.WriteLine("FRESNEL_SIN_VALUES_TEST:");
        Console.WriteLine("  FRESNEL_SIN_VALUES stores values of");
        Console.WriteLine("  the Fresnel sine integral S(X).");
        Console.WriteLine("");
        Console.WriteLine("      X           S(X)");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Fresnel.fresnel_sin_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString().PadLeft(12) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

}