using System;
using Burkardt.Values;

namespace TestValuesTest;

public class DebyeTest
{
    public static void debye1_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DEBYE1_VALUES_TEST tests DEBYE1_VALUES.
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
        int n_data;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("DEBYE1_VALUES_TEST:");
        Console.WriteLine("  DEBYE1_VALUES stores values of ");
        Console.WriteLine("  the Debye function of order 1.");
        Console.WriteLine("");
        Console.WriteLine("                X                     FX");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Debye.debye1_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString("0.################").PadLeft(24) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void debye2_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DEBYE2_VALUES_TEST tests DEBYE2_VALUES.
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
        int n_data;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("DEBYE2_VALUES_TEST:");
        Console.WriteLine("  DEBYE2_VALUES stores values of ");
        Console.WriteLine("  the Debye function of order 2.");
        Console.WriteLine("");
        Console.WriteLine("                X                     FX");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Debye.debye2_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString("0.################").PadLeft(24) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void debye3_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DEBYE3_VALUES_TEST tests DEBYE3_VALUES.
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
        int n_data;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("DEBYE3_VALUES_TEST:");
        Console.WriteLine("  DEBYE3_VALUES stores values of ");
        Console.WriteLine("  the Debye function of order 3.");
        Console.WriteLine("");
        Console.WriteLine("                X                     FX");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Debye.debye3_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString("0.################").PadLeft(24) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void debye4_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DEBYE4_VALUES_TEST tests DEBYE4_VALUES.
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
        int n_data;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("DEBYE4_VALUES_TEST:");
        Console.WriteLine("  DEBYE4_VALUES stores values of ");
        Console.WriteLine("  the Debye function of order 4.");
        Console.WriteLine("");
        Console.WriteLine("                X                     FX");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Debye.debye4_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString("0.################").PadLeft(24) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }
}