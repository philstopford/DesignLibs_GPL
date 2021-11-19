using System;
using Burkardt.Function;

namespace PolPakTest;

public static class alignTest
{
    public static void align_enum_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ALIGN_ENUM_TEST tests ALIGN_ENUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int M_MAX = 10;
        const int N_MAX = 10;

        int i;
        int j;

        Console.WriteLine("");
        Console.WriteLine("ALIGN_ENUM_TEST");
        Console.WriteLine("  ALIGN_ENUM counts the number of possible");
        Console.WriteLine("  alignments of two biological sequences.");

        Console.WriteLine("");
        Console.WriteLine("  Alignment enumeration table:");
        Console.WriteLine("");

        string cout = "      ";
        for ( j = 0; j <= 5; j++ )
        {
            cout += j.ToString().PadLeft(8) + "  ";
        }
        Console.WriteLine(cout);
        Console.WriteLine("");

        for ( i = 0; i <= M_MAX; i++ )
        {
            cout = "  " + i.ToString().PadLeft(2) + "  ";
            for ( j = 0; j <= 5; j++ )
            {
                cout += Align.align_enum ( i, j ).ToString().PadLeft(8) + "  ";
            }
            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        cout = "      ";
        for ( j = 6; j <= N_MAX; j++ )
        {
            cout += j.ToString().PadLeft(8) + "  ";
        }
        Console.WriteLine(cout);
        Console.WriteLine("");

        for ( i = 0; i <= M_MAX; i++ )
        {
            cout = "  " + i.ToString().PadLeft(2) + "  ";
            for ( j = 6; j <= N_MAX; j++ )
            {
                cout += Align.align_enum ( i, j ).ToString().PadLeft(8) + "  ";
            }
            Console.WriteLine(cout);
        }
    }
}