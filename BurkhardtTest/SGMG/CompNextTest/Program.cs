using System;
using Burkardt.Composition;

namespace CompNextTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for COMP_NEXT_TEST.
            //
            //  Discussion:
            //
            //     COMP_NEXT_TEST tests COMP_NEXT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    21 December 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Fabio Nobile, Raul Tempone, Clayton Webster,
            //    A Sparse Grid Stochastic Collocation Method for Partial Differential
            //    Equations with Random Input Data,
            //    SIAM Journal on Numerical Analysis,
            //    Volume 46, Number 5, 2008, pages 2309-2345.
            //
        {
            int dim_num;
            int[] dim_num_array =
            {
                2, 2, 2, 2, 2,
                3, 3, 3, 3, 3,
                4, 4
            };
            int level_max;
            int[] level_max_array =
            {
                0, 1, 2, 3, 4,
                0, 1, 2, 3, 4,
                2, 3
            };
            int test;
            int test_num = 12;

            Console.WriteLine("");
            Console.WriteLine(" COMP_NEXT_TEST");
            Console.WriteLine("");
            Console.WriteLine("  Test COMP_NEXT.");
            //
            //  Check that COMP_NEXT generates compositions correctly.
            //
            for (test = 0; test < test_num; test++)
            {
                dim_num = dim_num_array[test];
                level_max = level_max_array[test];
                comp_next_test(dim_num, level_max);
            }

            //
            //  Terminate.
            //
            Console.WriteLine("");
            Console.WriteLine(" COMP_NEXT_TEST");
            Console.WriteLine("  Normal end of execution.");

            Console.WriteLine("");
        }

        static void comp_next_test(int dim_num, int level_max)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    COMP_NEXT_TEST tests COMP_NEXT, which computes 1D level vectors.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 August 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, int LEVEL_MAX, the maximum level.
            //
        {
            int dim;
            int h;
            int i;
            int level;
            int[] level_1d;
            int level_min;
            bool more_grids;
            int t;

            level_1d = new int[dim_num];
            level_min = Math.Max(0, level_max + 1 - dim_num);

            Console.WriteLine("");
            Console.WriteLine("COMP_NEXT_TEST");
            Console.WriteLine("  COMP_NEXT generates, one at a time, vectors");
            Console.WriteLine("  LEVEL_1D(1:DIM_NUM) whose components add up to LEVEL.");
            Console.WriteLine("");
            Console.WriteLine("  We call with:");
            Console.WriteLine("  DIM_NUM = " + dim_num + "");
            Console.WriteLine("  " + level_min + " = LEVEL_MIN <= LEVEL <= LEVEL_MAX = "
                              + level_max + "");
            Console.WriteLine("");
            Console.WriteLine("     LEVEL     INDEX  LEVEL_1D Vector");
            //
            //  The outer loop generates values of LEVEL from LEVEL_MIN to LEVEL_MAX.
            //
            for (level = level_min; level <= level_max; level++)
            {
                Console.WriteLine("");
                //
                //  The inner loop generates vectors LEVEL_1D(1:DIM_NUM) whose components 
                //  add up to LEVEL.
                //
                more_grids = false;
                h = 0;
                t = 0;
                i = 0;

                for (;;)
                {
                    Comp.comp_next(level, dim_num, ref level_1d, ref more_grids, ref h, ref t);

                    i = i + 1;
                    string cout = "  " + level.ToString().PadLeft(8)
                                       + "  " + i.ToString().PadLeft(8);
                    for (dim = 0; dim < dim_num; dim++)
                    {
                        cout += "  " + level_1d[dim].ToString().PadLeft(8);
                    }

                    Console.WriteLine(cout);

                    if (!more_grids)
                    {
                        break;
                    }
                }
            }
        }
    }
}