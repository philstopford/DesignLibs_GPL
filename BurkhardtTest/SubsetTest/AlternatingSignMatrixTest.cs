using System;
using Burkardt.MatrixNS;

namespace SubsetTestNS;

public static class AlternatingSignMatrixTest
{
    public static void asm_enum_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ASM_ENUM_TEST tests ASM_ENUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;

        Console.WriteLine("");
        Console.WriteLine("ASM_ENUM_TEST");
        Console.WriteLine("  ASM_ENUM returns the number of alternating sign");
        Console.WriteLine("  matrices of a given order.");

        Console.WriteLine("");

        for ( n = 0; n <= 7; n++ )
        {
            Console.WriteLine(n.ToString().PadLeft(4) + "  "
                                                      + AlternatingSignMatrix.asm_enum ( n ).ToString().PadLeft(6) + "");
        }
    }

    public static void asm_triangle_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ASM_TRIANGLE_TEST tests ASM_TRIANGLE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N_MAX = 7;

        int[] a = new int[N_MAX+1];
        int n;

        Console.WriteLine("");
        Console.WriteLine("ASM_TRIANGLE_TEST");
        Console.WriteLine("  ASM_TRIANGLE returns a row of the alternating sign");
        Console.WriteLine("  matrix triangle.");
        Console.WriteLine("");

        for ( n = 0; n <= N_MAX; n++ )
        {
            AlternatingSignMatrix.asm_triangle ( n, ref a );
            string cout = n.ToString().PadLeft(4) + "  ";
            int i;
            for ( i = 0; i <= n; i++ )
            {
                cout += a[i].ToString().PadLeft(8) + "  ";
            }
            Console.WriteLine(cout);
        }

    }
}