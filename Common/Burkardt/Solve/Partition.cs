using System;
using Burkardt.SubsetNS;
using Burkardt.Types;

namespace Burkardt.SolveNS;

public static class Partition
{
    public static void partition_brute(int n, int[] w, ref int[] c, ref int discrepancy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARTITION_BRUTE approaches the partition problem using brute force.
        //
        //  Discussion:
        //
        //    We are given a set of N integers W.
        //
        //    We seek to partition W into subsets W0 and W1, such that the subsets
        //    have equal sums.
        //
        //    The "discrepancy" is the absolute value of the difference between the
        //    two sums, and will be zero if we have solved the problem.
        //
        //    For a given set of integers, there may be zero, one, or many solutions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 May 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the size of the set.
        //
        //    Input, int W[N], the integers.
        //
        //    Output, int C[N], indicates the proposed solution.
        //    C(I) is 0 for items in set W0 and 1 for items in set W1.
        //
        //    Output, int &DISCREPANCY, the discrepancy.
        //
    {
        int w_sum = typeMethods.i4vec_sum(n, w);
        discrepancy = w_sum;

        int rank = -1;
        int[] d = new int[n];

        while (true)
        {
            Subset.subset_next(n, ref d, ref rank);

            if (rank == -1)
            {
                break;
            }

            int d_discrepancy = Math.Abs(w_sum - 2 * typeMethods.i4vec_dot_product(n, d, w));

            if (d_discrepancy < discrepancy)
            {
                discrepancy = d_discrepancy;
                typeMethods.i4vec_copy(n, d, ref c);
            }

            if (discrepancy == 0)
            {
                break;
            }
        }
    }

    public static int partition_count(int n, int[] w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARTITION_COUNT counts the solutions to a partition problem.
        //
        //  Discussion:
        //
        //    We are given a set of N integers W.
        //
        //    We seek to partition W into subsets W0 and W1, such that the subsets
        //    have equal sums.
        //
        //    The "discrepancy" is the absolute value of the difference between the
        //    two sums, and will be zero if we have solved the problem.
        //
        //    For a given set of integers, there may be zero, one, or many solutions.
        //
        //    In the case where the weights are distinct, the count returned by this
        //    function may be regarded as twice as big as it should be, since the
        //    partition (W0,W1) is counted a second time as (W1,W0).  A more serious
        //    overcount can occur if the set W contains duplicate elements - in the
        //    extreme case, W might be entirely 1's, in which case there is really
        //    only one (interesting) solution, but this function will count many.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 May 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the size of the set.
        //
        //    Input, int W[N], the integers.
        //
        //    Output, int PARTITION_COUNT, the number of solutions.
        //
    {
        int w_sum = typeMethods.i4vec_sum(n, w);

        int[] c = new int[n];
        int rank = -1;
        int count = 0;

        while (true)
        {
            Subset.subset_next(n, ref c, ref rank);

            if (rank == -1)
            {
                break;
            }

            int discrepancy = Math.Abs(w_sum - 2 * typeMethods.i4vec_dot_product(n, c, w));

            switch (discrepancy)
            {
                case 0:
                    count += 1;
                    break;
            }
        }
            
        return count;
    }
}