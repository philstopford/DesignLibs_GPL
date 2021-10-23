using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace GeometryTest
{
    public static class r8MatSolveTest
    {
        public static void test0234 ( )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST0234 tests R8MAT_SOLVE_2D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    16 November 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] a;
            double[] b = null;
            double det = 0;
            int i;
            int n = 2;
            int seed;
            int test;
            int test_num = 5;
            double[] x;
            double[] x2;

            Console.WriteLine("");
            Console.WriteLine("TEST0234");
            Console.WriteLine("  R8MAT_SOLVE_2D solves 2D linear systems.");

            seed = 123456789;

            for ( test = 1; test <= test_num; test++ )
            {
                a = UniformRNG.r8mat_uniform_01_new ( n, n, ref seed );
                x = UniformRNG.r8vec_uniform_01_new ( n, ref seed );
                typeMethods.r8mat_mv ( n, n, a, x, ref b );

                x2 = typeMethods.r8mat_solve_2d ( a, b, ref det );

                Console.WriteLine("");
                Console.WriteLine("  Solution / Computed:");
                Console.WriteLine("");

                for ( i = 0; i < n; i++ )
                {
                    Console.WriteLine("  " + x[i].ToString().PadLeft(14)
                                      + "  " + x2[i].ToString().PadLeft(14) + "");
                }

            }
        }

    }
}