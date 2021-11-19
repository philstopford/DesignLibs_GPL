using System;
using Burkardt.AppliedStatistics;
using Burkardt.Types;

namespace ASA314Test;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for ASA314_TEST.
        //
        //  Discussion:
        //
        //    ASA314_TEST tests ASA314.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Roger Payne,
        //    Inversion of matrices with contents subject to modulo arithmetic,
        //    Applied Statistics,
        //    Volume 46, Number 2, 1997, pages 295-298.
        //
    {
        Console.WriteLine("");
        Console.WriteLine("ASA314_TEST:");
        Console.WriteLine("  Test the ASA314 library.");

        test01 ( );

        Console.WriteLine("");
        Console.WriteLine("ASA314_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
            
    }

    private static void test01 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests INVMOD.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Roger Payne,
        //    Inversion of matrices with contents subject to modulo arithmetic,
        //    Applied Statistics,
        //    Volume 46, Number 2, 1997, pages 295-298.
        //
    {
        int ifault = 0;
        int[] jmat = { 1, 0, 0, 2, 1, 0, 1, 0, 1 };
        int[] mat = { 1, 0, 0, 1, 1, 0, 2, 0, 1 };

        const int nrow = 3;
        int[] cmod = new int[nrow];
        int[] rmod = new int[nrow];
        int[] imat = new int[nrow*nrow];

        for (int i = 0; i < nrow; i++ )
        {
            cmod[i] = 3;
            rmod[i] = 3;
        }

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  INVMOD computes the inverse of a matrix");
        Console.WriteLine("  whose elements are subject to modulo arithmetic.");

        typeMethods.i4mat_print ( nrow, nrow, mat, "  The matrix to be inverted:" );

        Algorithms.invmod ( ref mat, ref imat, rmod, cmod, nrow, ref ifault );

        typeMethods.i4mat_print ( nrow, nrow, imat, "  The computed inverse:" );

        typeMethods.i4mat_print ( nrow, nrow, jmat, "  The correct inverse:" );
    }
        
}