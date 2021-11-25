using System;
using Burkardt.SortNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace SortRCTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SORT_RC_TEST.
        //
        //  Discussion:
        //
        //    SORT_RC_TEST tests the SORT_RC library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 March 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("SORT_RC_TEST");
        Console.WriteLine("  Test the SORT_RC library.");

        sort_rc_test();
        sort_safe_rc_test();

        Console.WriteLine("");
        Console.WriteLine("SORT_RC_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void sort_rc_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SORT_RC_TEST tests SORT_RC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 March 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i = 0;
        int isgn = 0;
        int j = 0;
        const int n = 20;
        Sort.SortRCData data = new();

        Console.WriteLine("");
        Console.WriteLine("SORT_RC_TEST");
        Console.WriteLine("  SORT_RC sorts objects externally.");
        Console.WriteLine("  This function relies on the use of persistent");
        Console.WriteLine("  data stored internally.");
        //
        //  Generate some data to sort.
        //
        const int i4_lo = 1;
        int seed = 123456789;

        int[] a = UniformRNG.i4vec_uniform_ab_new(n, i4_lo, n, ref seed);

        typeMethods.i4vec_print(n, a, "  Unsorted array:");
        //
        //  Sort the data.
        // 
        int indx = 0;

        for (;;)
        {
            Sort.sort_rc(ref data, n, ref indx, ref i, ref j, isgn);

            if (indx < 0)
            {
                isgn = 1;
                if (a[i - 1] <= a[j - 1])
                {
                    isgn = -1;
                }
            }
            else if (0 < indx)
            {
                (a[i - 1], a[j - 1]) = (a[j - 1], a[i - 1]);
            }
            else
            {
                break;
            }
        }

        //
        //  Display the sorted data.
        //
        typeMethods.i4vec_print(n, a, "  Sorted array:");

    }

    private static void sort_safe_rc_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SORT_SAFE_RC_TEST tests SORT_SAFE_RC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 March 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i = 0;
        int i_save = 0;
        int isgn = 0;
        int j = 0;
        int j_save = 0;
        int k_save = 0;
        int l_save = 0;
        const int n = 20;
        int n_save = 0;

        Console.WriteLine("");
        Console.WriteLine("SORT_SAFE_RC_TEST");
        Console.WriteLine("  SORT_SAFE_RC sorts objects externally.");
        Console.WriteLine("  This version of the algorithm does not rely on");
        Console.WriteLine("  internally saved or 'persistent' or 'static' memory.");
        //
        //  Generate some data to sort.
        //
        const int i4_lo = 1;
        int seed = 123456789;

        int[] a = UniformRNG.i4vec_uniform_ab_new(n, i4_lo, n, ref seed);

        typeMethods.i4vec_print(n, a, "  Unsorted array:");
        //
        //  Sort the data.
        //
        int indx = 0;

        for (;;)
        {
            Sort.sort_safe_rc(n, ref indx, ref i, ref j, isgn, ref i_save, ref j_save, ref k_save, ref l_save,
                ref n_save);

            if (indx < 0)
            {
                isgn = 1;
                if (a[i - 1] <= a[j - 1])
                {
                    isgn = -1;
                }
            }
            else if (0 < indx)
            {
                (a[i - 1], a[j - 1]) = (a[j - 1], a[i - 1]);
            }
            else
            {
                break;
            }
        }

        //
        //  Display the sorted data.
        //
        typeMethods.i4vec_print(n, a, "  Sorted array:");
    }
}