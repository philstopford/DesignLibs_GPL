using System;
using Burkardt.Latin;
using Burkardt.Types;

namespace LatinEdgeTest;

internal class Program
{
    private static void Main(string[] args)
    {
        int seed = typeMethods.get_seed();
        int seed_save = seed;

        Console.WriteLine();
        Console.WriteLine("LATIN_EDGE_TEST:");
        Console.WriteLine("  Test the LATIN_EDGE library.");

        test01 ( ref seed );

        Console.WriteLine();
        Console.WriteLine("LATIN_EDGE_TEST:");
        Console.WriteLine("  Repeat TEST01, but with different seed from first run.");
            
        seed = typeMethods.get_seed();
            
        test01 ( ref seed );

        Console.WriteLine();
        Console.WriteLine("LATIN_EDGE_TEST:");
        Console.WriteLine("  Repeat TEST01 with same seed as first run.");

        seed = seed_save;
        test01 ( ref seed );

        Console.WriteLine();
        Console.WriteLine("LATIN_EDGE_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine();
    }


    private static void test01 ( ref int seed )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests LATIN_EDGE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 March 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int POINT_NUM = 10;

        Console.WriteLine();
        Console.WriteLine("TEST01");
        Console.WriteLine("  LATIN_EDGE chooses a Latin cell arrangement,");
        Console.WriteLine("  which includes the edge points.");
        Console.WriteLine();
        Console.WriteLine("  Spatial dimension = " + DIM_NUM);
        Console.WriteLine("  Number of points =  " + POINT_NUM);
        Console.WriteLine("  Initial seed for UNIFORM = " + seed);

        double[] x = LatinVariants.latin_edge( DIM_NUM, POINT_NUM, ref seed );

        Console.WriteLine();
        Console.WriteLine("  The Latin Edge Square points:");
        Console.WriteLine();

        int k = 0;
        for (int j = 0; j < POINT_NUM; j++ )
        {
            int kk = k;
            string cout = "";
            for (int i = 0; i < DIM_NUM; i++ )
            {
                cout += x[kk].ToString("0.########").PadLeft(10) + "  ";
                kk += POINT_NUM;
            }
            Console.WriteLine(cout);
            k += 1;
        }
    }
}