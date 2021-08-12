namespace Burkardt
{
    public static class RestrictedGrowth
    {
        public static void regro_next(ref bool done, int n, ref int[] v, ref int[] vmax )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    REGRO_NEXT computes restricted growth functions one at a time.
        //
        //  Discussion:
        //
        //    A restricted growth function on N is a vector (V(1), ..., V(N) )
        //    of values V(I) between 1 and N, satisfying the requirements:
        //      V(1) = 1;
        //      V(I) <= 1 + max ( V(1), V(2), ..., V(I-1) ).
        //
        //    The number of restricted growth functions on N is equal to
        //    the Bell number B(N).
        //
        //    There is a bijection between restricted growth functions on N
        //    and set partitions of N.
        //
        //  Example:
        //
        //    The 15 restricted growth functions for N = 4 are:
        //
        //    (1111), (1112), (1121), (1122), (1123),
        //    (1211), (1212), (1213), (1221), (1222),
        //    (1223), (1231), (1232), (1233), (1234).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 May 2003
        //
        //  Author:
        //
        //    John Burkardt.
        //
        //  Reference:
        //
        //    Dennis Stanton, Dennis White,
        //    Constructive Combinatorics,
        //    Springer, 1986,
        //    ISBN: 0387963472,
        //    LC: QA164.S79.
        //
        //  Parameters:
        //
        //    Input/output, bool &DONE.
        //    On first call, set DONE to TRUE, and then do not alter it.
        //    On output, DONE will be FALSE if the routine has computed another
        //    restricted growth function, or TRUE if all the restricted
        //    growth functions have been returned.
        //
        //    Input, int N, the number of components in the restricted growth
        //    function.
        //
        //    Input/output, int V[N].  The user need not set this quantity
        //    before the initial call, and should not alter it between successive
        //    calls.  On each return from the routine, with DONE = .FALSE.,
        //    V will contain the componentwise values of the next restricted
        //    growth function.
        //
        //    Input/output, int VMAX[N].  The user need not set this quantity
        //    before the initial call, and should not alter it between calls.
        //    VMAX(I) records the largest value that component V(I) could take,
        //    given the values of components 1 through I-1.
        //
        {
            int i;
            int j;
            //
            //  First call:
            //
            if (done)
            {
                for (i = 0; i < n; i++)
                {
                    v[i] = 1;
                }

                vmax[0] = 1;
                for (i = 1; i < n; i++)
                {
                    vmax[i] = 2;
                }

                done = false;
            }
            //
            //  Later calls.
            //
            else
            {
                j = n;

                for (;;)
                {
                    if (j == 1)
                    {
                        done = true;
                        return;
                    }

                    if (v[j - 1] != vmax[j - 1])
                    {
                        break;
                    }

                    j = j - 1;

                }

                v[j - 1] = v[j - 1] + 1;

                for (i = j + 1; i <= n; i++)
                {
                    v[i - 1] = 1;

                    if (v[j - 1] == vmax[j - 1])
                    {
                        vmax[i - 1] = vmax[j - 1] + 1;
                    }
                    else
                    {
                        vmax[i - 1] = vmax[j - 1];
                    }
                }
            }

        }
    }
}