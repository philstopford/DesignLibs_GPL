using System;
using System.Globalization;
using Burkardt.Values;

namespace TestValuesTest;

public static class ViscosityTest
{
    public static void viscosity_values_test()
        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    VISCOSITY_VALUES_TEST tests VISCOSITY_VALUES.
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
        double eta = 0;
        double p = 0;
        double tc = 0;
        Console.WriteLine("");
        Console.WriteLine("VISCOSITY_VALUES_TEST:");
        Console.WriteLine("  VISCOSITY_VALUES stores values of");
        Console.WriteLine("  the viscosity of water");
        Console.WriteLine("  as a function of temperature and pressure.");
        Console.WriteLine("");
        Console.WriteLine("      T            P            ETA(T,P)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Viscosity.viscosity_values(ref n_data, ref tc, ref p, ref eta);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + tc.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + p.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + eta.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }
}