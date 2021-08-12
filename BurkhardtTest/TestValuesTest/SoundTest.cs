using System;
using Burkardt.TestValues;

namespace TestValuesTest
{
    public class SoundTest
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
            int n_data;
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
            n_data = 0;
            for (;;)
            {
                Sound.sound_values(ref n_data, ref tc, ref p, ref c);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + tc.ToString().PadLeft(12) + "  "
                                  + p.ToString().PadLeft(12) + "  "
                                  + c.ToString().PadLeft(12) + "");
            }
        }

    }
}