﻿using System;
using TestValues;

namespace TestValuesTest
{
    public class TsatTest
    {
        public static void tsat_values_test()
            //****************************************************************************80
            //
            //  Purpose: 
            //
            //    TSAT_VALUES_TEST tests TSAT_VALUES.
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
            int n_data;
            double p = 0;
            double tc = 0;
            Console.WriteLine("");
            Console.WriteLine("TSAT_VALUES_TEST:");
            Console.WriteLine("  TSAT_VALUES stores values of");
            Console.WriteLine("  the saturation temperature");
            Console.WriteLine("  as a function of pressure.");
            Console.WriteLine("");
            Console.WriteLine("      P           Tsat(P)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Tsat.tsat_values(ref n_data, ref p, ref tc);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + p.ToString().PadLeft(12) + p + "  "
                                  + tc.ToString().PadLeft(12) + tc + "");
            }
        }
    }
}