using System;
using Burkardt.Table;
using Burkardt.Types;
using Random = Burkardt.Latin.Random;

namespace Burkardt.LatinRandomTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for LATIN_RANDOM_TEST.
            //
            //  Discussion:
            //
            //    LATIN_RANDOM_TEST tests the LATIN_RANDOM library.
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
            Console.WriteLine();
            Console.WriteLine("LATIN_RANDOM_TEST:");
            Console.WriteLine("  Test the LATIN_RANDOM library.");

            int seed = 123456789;

            for (int test = 0; test < 3; test++)
            {
                test01(ref seed);
            }

            Console.WriteLine();
            Console.WriteLine("LATIN_RANDOM_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine();
        }

        static void test01(ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests LATIN_RANDOM_NEW.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 November 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        {
            int m = 2;
            int n = 10;

            Console.WriteLine();
            Console.WriteLine("TEST01");
            Console.WriteLine("  LATIN_RANDOM chooses a Latin Square cell arrangement,");
            Console.WriteLine("  and then chooses a random point from each cell.");
            Console.WriteLine();
            Console.WriteLine("  Spatial dimension = " + m);
            Console.WriteLine("  Number of points =  " + n);
            Console.WriteLine("  Initial seed for UNIFORM = " + seed);

            double[] x = Random.latin_random_new(m, n, ref seed);

            typeMethods.r8mat_transpose_print(m, n, x, "  Latin Random Square:");
        }
    }
}
