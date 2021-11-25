﻿using System;
using System.Globalization;
using Burkardt.Values;

namespace TestValuesTest;

public static class StruveTest
{
    public static void struve_h0_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STRUVE_H0_VALUES_TEST tests STRUVE_H0_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("STRUVE_H0_VALUES_TEST:");
        Console.WriteLine("  STRUVE_H0_VALUES stores values of");
        Console.WriteLine("  the Struve H0 function.");
        Console.WriteLine("");
        Console.WriteLine("      X            H0(X)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Struve.struve_h0_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void struve_h1_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STRUVE_H1_VALUES_TEST tests STRUVE_H1_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("STRUVE_H1_VALUES_TEST:");
        Console.WriteLine("  STRUVE_H1_VALUES stores values of");
        Console.WriteLine("  the Struve H1 function.");
        Console.WriteLine("");
        Console.WriteLine("      X            H1(X)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Struve.struve_h1_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void struve_l0_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STRUVE_L0_VALUES_TEST tests STRUVE_L0_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("STRUVE_L0_VALUES_TEST:");
        Console.WriteLine("  STRUVE_L0_VALUES stores values of");
        Console.WriteLine("  the Struve L0 function.");
        Console.WriteLine("");
        Console.WriteLine("      X            L0(X)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Struve.struve_l0_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void struve_l1_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STRUVE_L1_VALUES_TEST tests STRUVE_L1_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("STRUVE_L1_VALUES_TEST:");
        Console.WriteLine("  STRUVE_L1_VALUES stores values of");
        Console.WriteLine("  the Struve L1 function.");
        Console.WriteLine("");
        Console.WriteLine("      X            L1(X)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Struve.struve_l1_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }
}