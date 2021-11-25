using System;
using Burkardt.MatrixNS;
using Burkardt.Types;

namespace SubsetTestNS;

public static class PartialOrderingTest
{
    public static void pord_check_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PORD_CHECK_TEST tests PORD_CHECK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] a = {
            1,0,1,0,1,0,1,0,0,1, 
            0,1,0,0,1,0,0,0,0,0, 
            0,0,1,0,1,0,1,0,0,1, 
            0,1,1,1,1,1,1,1,0,1, 
            0,0,0,0,1,0,0,0,0,0, 
            0,1,0,0,1,1,1,0,0,0, 
            0,0,0,0,1,0,1,0,0,0, 
            0,1,0,0,1,1,1,1,0,1, 
            0,0,0,0,0,0,0,0,0,0, 
            0,0,0,0,1,0,1,0,0,1 };
        const int n = 10;

        Console.WriteLine("");
        Console.WriteLine("PORD_CHECK_TEST");
        Console.WriteLine("  PORD_CHECK checks a partial ordering.");

        typeMethods.i4mat_print ( n, n, a, "  The partial ordering matrix:" );
 
        bool ierror = PartialOrdering.pord_check ( n, a );
 
        Console.WriteLine("");
        Console.WriteLine("  CHECK FLAG = " + ierror + "");
        Console.WriteLine("  0 means no error.");
        Console.WriteLine("  1 means illegal value of N.");
        Console.WriteLine("  2 means some A(I,J) and A(J,I) are both nonzero.");
    }
}