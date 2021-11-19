using System;
using Burkardt.Values;

namespace TestValuesTest;

public class ClebschGordanTest
{
    public static void clebsch_gordan_values_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CLEBSCH_GORDAN_VALUES_TEST tests CLEBSCH_GORDAN_VALUES.
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
        double j1 = 0;
        double j2 = 0;
        double j3 = 0;
        double m1 = 0;
        double m2 = 0;
        double m3 = 0;
        int n_data;
        Console.WriteLine("");
        Console.WriteLine("CLEBSCH_GORDAN_VALUES_TEST:");
        Console.WriteLine("  CLEBSCH_GORDAN_VALUES returns values of");
        Console.WriteLine("  the Clebsch Gordan coefficient.");
        Console.WriteLine("");
        Console.WriteLine("      J1      J2      J3      M1      M2      M3        CG");
        Console.WriteLine("");
        n_data = 0;
        for ( ; ; )
        {
            ClebschGordan.clebsch_gordan_values ( ref n_data, ref j1, ref j2, ref j3, ref m1, ref m2, ref m3, ref fx );
            if ( n_data == 0 )
            {
                break;
            }
            Console.WriteLine("  " + j1.ToString().PadLeft(6)
                                   + "  " + j2.ToString().PadLeft(6)
                                   + "  " + j3.ToString().PadLeft(6)
                                   + "  " + m1.ToString().PadLeft(6)
                                   + "  " + m2.ToString().PadLeft(6)
                                   + "  " + m3.ToString().PadLeft(6) + m3
                                   + "  " + fx.ToString("0.################").PadLeft(24) + "");
        }
    }
}