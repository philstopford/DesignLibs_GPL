using Burkardt.SubsetNS;

namespace Burkardt
{
    public static class Knapsack
    {
        public static int[] knapsack_01 ( int n, int[] w, int c )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    KNAPSACK_01 seeks a solution of the 0/1 Knapsack problem.
        //
        //  Discussion:
        //
        //    In the 0/1 knapsack problem, a knapsack of capacity C is given,
        //    as well as N items, with the I-th item of weight W(I).
        //
        //    A selection is "acceptable" if the total weight is no greater than C.
        //
        //    It is desired to find an optimal acceptable selection, that is,
        //    an acceptable selection such that there is no acceptable selection
        //    of greater weight.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of weights.
        //
        //    Input, inte W[N], the weights.
        //
        //    Input, int C, the maximum weight.
        //
        //    Output, int KNAPSACK_01[N], is a binary vector which defines an 
        //    optimal selection.  It is 1 for the weights to be selected, and 
        //    0 otherwise.
        //
        {
            int i;
            int iadd = 0;
            bool more;
            int ncard;
            int[] s;
            int[] s_test;
            int t;
            int t_test;

            s = new int[n];
            s_test = new int[n];

            more = false;
            ncard = 0;

            for ( i = 0; i < n; i++ )
            {
                s_test[i] = 0;
            }
            t_test = 0;

            for ( i = 0; i < n; i++ )
            {
                s[i] = s_test[i];
            }
            t = 0;

            for ( ; ; )
            {
                Subset.subset_gray_next ( n, ref s_test, ref more, ref ncard, ref iadd );
                t_test = 0;
                for ( i = 0; i < n; i++ )
                {
                    t_test = t_test + s_test[i] * w[i];
                }

                if ( t < t_test && t_test <= c )
                {
                    t = t_test;
                    for ( i = 0; i < n; i++ )
                    {
                        s[i] = s_test[i];
                    }
                }

                if ( ! more )
                {
                    break;
                }
            }
            
            return s;
        }
    }
}