using System;
using System.Linq;
using Burkardt;
using Burkardt.Function;
using Burkardt.SortNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace BinsTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for BINS_TEST.
            //
            //  Discussion:
            //
            //    BINS_TEST tests the BINS library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 April 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("BINS_TEST");
            Console.WriteLine("  Test the BINS library.");

            test0081();
            test0835();
            test0836();
            test180();

            Console.WriteLine("");
            Console.WriteLine("BINS_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test0081()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST0081 tests BIN_TO_R8_EVEN2, R8_TO_BIN_EVEN2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 April 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a = 10.0;
            double b = 20.0;
            double cmax = 0;
            double cmin = 0;
            int nbin = 5;
            double rmax = 23.0;
            double rmin = 8.0;

            Console.WriteLine("");
            Console.WriteLine("TEST0676");
            Console.WriteLine("  R8_TO_BIN_EVEN2 puts a number into a bin.");
            Console.WriteLine("  BIN_TO_R8_EVEN2 returns the bin limits.");
            Console.WriteLine("  The bins are equally spaced between A and B.");
            Console.WriteLine("");
            Console.WriteLine("  A = " + a + "");
            Console.WriteLine("  B = " + b + "");
            Console.WriteLine("  Total number of bins = " + nbin + "");
            Console.WriteLine("");
            Console.WriteLine("  Generate some random values C and put them in bins.");
            Console.WriteLine("");
            Console.WriteLine("       C      Bin   Bin_Min  Bin_Max");
            Console.WriteLine("");

            int seed = typeMethods.get_seed();

            for (int i = 1; i <= 30; i++)
            {
                double c = UniformRNG.r8_uniform(rmin, rmax, ref seed);
                int bin = Bins.r8_to_bin_even2(nbin, a, b, c);
                Bins.bin_to_r8_even2(nbin, bin, a, b, ref cmin, ref cmax);
                Console.WriteLine("  " + c.ToString().PadLeft(10)
                                       + "  " + bin.ToString().PadLeft(6)
                                       + "  " + cmin.ToString().PadLeft(10)
                                       + "  " + cmax.ToString().PadLeft(10) + "");
            }
        }

        static void test0835()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST0835 tests R82VEC_PART_QUICK_A.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 April 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 12;

            double[] a = new double[N * 2];
            double[] ahi =
                {
                    10.0, 3.0
                }
                ;
            double[] alo =
                {
                    0.0, 2.0
                }
                ;
            int l = 0;
            int r = 0;
            int seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("TEST0835");
            Console.WriteLine("  R82VEC_PART_QUICK_A reorders an R2 vector");
            Console.WriteLine("    as part of a quick sort.");
            Console.WriteLine("  Using initial random number seed = " + seed + "");

            UniformRNG.r82vec_uniform(N, alo, ahi, ref seed, ref a);

            typeMethods.r82vec_print(N, a, "  Before rearrangment:");

            typeMethods.r82vec_part_quick_a(N, ref a, 0, ref l, ref r);

            Console.WriteLine("");
            Console.WriteLine("  Rearranged array");
            Console.WriteLine("  Left index =  " + l + "");
            Console.WriteLine("  Key index =   " + l + 1 + "");
            Console.WriteLine("  Right index = " + r + "");

            typeMethods.r82vec_print(l, a, "  Left half:");
            typeMethods.r82vec_print(1, a.Skip(2 * l).ToArray(), "  Key:");
            typeMethods.r82vec_print(N - l - 1, a.Skip(2 * (l + 1)).ToArray(), "  Right half:");

        }

        static void test0836()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST0836 tests R82VEC_SORT_QUICK_A.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 April 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 12;

            double[] a = new double[N * 2];
            double[] ahi =
                {
                    10.0, 3.0
                }
                ;
            double[] alo =
                {
                    0.0, 2.0
                }
                ;
            int seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("TEST0836");
            Console.WriteLine("  D2VEC_SORT_QUICK_A sorts an R2 vector");
            Console.WriteLine("    as part of a quick sort.");
            Console.WriteLine("  Using initial random number seed = " + seed + "");

            UniformRNG.r82vec_uniform(N, alo, ahi, ref seed, ref a);
            //
            //  For better testing, give a few elements the same first component.
            //
            a[2 * (3 - 1) + 0] = a[2 * (5 - 1) + 0];
            a[2 * (4 - 1) + 0] = a[2 * (12 - 1) + 0];
            //
            //  Make two entries equal.
            //
            a[2 * (7 - 1) + 0] = a[2 * (11 - 1) + 0];
            a[2 * (7 - 1) + 1] = a[2 * (11 - 1) + 1];

            typeMethods.r82vec_print(N, a, "  Before sorting:");

            typeMethods.r82vec_sort_quick_a(N, ref a);

            typeMethods.r82vec_print(N, a, "  Sorted array:");


        }

        static void test180()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST180 tests SORT_HEAP_EXTERNAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 April 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 20;

            int[] a = new int[N];
            int i;
            int indx;
            int isgn;
            int j;
            int seed;

            Console.WriteLine("");
            Console.WriteLine("TEST180");
            Console.WriteLine("  SORT_HEAP_EXTERNAL sorts objects externally.");

            indx = 0;
            i = 0;
            j = 0;
            isgn = 0;
            seed = 123456789;

            for (i = 0; i < N; i++)
            {
                a[i] = UniformRNG.i4_uniform(1, N, ref seed);
            }

            typeMethods.i4vec_print(N, a, "  Unsorted array:");
            //
            //  Call the sort routine over and over.
            //
            SortHeapExternalData data = new SortHeapExternalData();
            for (;;)
            {
                Sort.sort_heap_external(ref data, N, ref indx, ref i, ref j, isgn);
                //
                //  If the return value of INDX is negative, we're asked to compare
                //  array elements I and J;
                //
                if (indx < 0)
                {
                    if (a[i] <= a[j])
                    {
                        isgn = -1;
                    }
                    else
                    {
                        isgn = 1;
                    }

                }
                //
                //  ...and if the return value of INDX is positive, we're asked to switch
                //  array elements I and J;
                //
                else if (0 < indx)
                {
                    typeMethods.i4_swap(ref a[i], ref a[j]);
                    //
                    //  ...and if the return value of INDX is 0, we're done.
                    //
                }
                else
                {
                    break;
                }

            }

            typeMethods.i4vec_print(N, a, "  Sorted array:");
        }
    }
}