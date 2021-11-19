using System;
using Burkardt.Latin;
using Burkardt.Types;

namespace LatinCoverTest;

internal static class Program
{
    private static void Main()
    {
        Console.WriteLine();
        ;
        Console.WriteLine("LATIN_COVER_TEST:");
        Console.WriteLine("  Test the LATIN_COVER library.");

        test01();
        test02();
        test03();

        Console.WriteLine();
        Console.WriteLine("LATIN_COVER_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine();
    }

    private static void test01 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests LATIN_COVER.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 August 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine();
        Console.WriteLine("TEST01");
        Console.WriteLine("  LATIN_COVER:");

        for (int n = 3; n <= 9; n += 2 )
        {
            int seed = 123456789;

            for (int test = 1; test <= 3; test++ )
            {
                int[] p = typeMethods.perm_uniform_new ( n, ref seed );
 
                typeMethods.perm_print ( n, p, "  Permutation" );

                int[] a = LatinVariants.latin_cover ( n, p );

                typeMethods.i4mat_print ( n, n, a, "  Latin cover" );
            }
        }
    }

    private static void test02 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //   TEST02 tests LATIN_COVER_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 June 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine();
        Console.WriteLine("TEST02");
        Console.WriteLine("  LATIN_COVER_2D:");

        for ( int n = 3; n <= 9; n += 2 )
        {
            int seed = 123456789;
            for ( int test = 1; test <= 3; test++ )
            {
                int[] p1 = typeMethods.perm_uniform_new ( n, ref seed );
                typeMethods.perm_print ( n, p1, "  Permutation 1" );

                int[] p2 = typeMethods.perm_uniform_new ( n, ref seed ); 
                typeMethods.perm_print ( n, p2, "  Permutation 2" );

                int[] a = LatinVariants.latin_cover_2d ( n, p1, p2 );
                typeMethods.i4mat_print ( n, n, a, "  Latin cover" );
            }
        }
    }

    private static void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests LATIN_COVER_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 August 2012
//
//  Author:
//
//    John Burkardt
//
    {
        Console.WriteLine();
        Console.WriteLine("TEST03");
        Console.WriteLine("  LATIN_COVER_3D");

        for ( int n = 3; n <= 9; n += 2 )
        {
            int seed = 123456789;
            for ( int test = 1; test <= 3; test++ )
            {
                int[] p1 = typeMethods.perm_uniform_new ( n, ref seed );
                typeMethods.perm_print ( n, p1, "  Permutation 1" );

                int[] p2 = typeMethods.perm_uniform_new ( n, ref seed );
                typeMethods.perm_print ( n, p2, "  Permutation 2" );

                int[] p3 = typeMethods.perm_uniform_new ( n, ref seed );
                typeMethods.perm_print ( n, p3, "  Permutation 1" );

                int[] a = LatinVariants.latin_cover_3d ( n, p1, p2, p3 );

                typeMethods.i4block_print ( n, n, n, a, "  Latin cover" );
            }
        }
    }
}