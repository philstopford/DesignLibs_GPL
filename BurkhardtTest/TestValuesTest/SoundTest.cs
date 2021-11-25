using System;
using System.Globalization;
using Burkardt.Values;

namespace TestValuesTest;

public static class SoundTest
{
    public static void sound_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SOUND_VALUES_TEST tests SOUND_VALUES.
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
        double c = 0;
        double p = 0;
        double tc = 0;
        Console.WriteLine("");
        Console.WriteLine("SOUND_VALUES_TEST:");
        Console.WriteLine("  SOUND_VALUES stores values of");
        Console.WriteLine("  the spead of sound in water");
        Console.WriteLine("  as a function of temperature and pressure.");
        Console.WriteLine("");
        Console.WriteLine("      T            P            C(T,P)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Sound.sound_values(ref n_data, ref tc, ref p, ref c);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + tc.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + p.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + c.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

}