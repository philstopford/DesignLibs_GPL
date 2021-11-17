using Burkardt.Types;

namespace Burkardt.RankingNS;

public static partial class Ranking
{
    public static int dist_enum(int k, int m)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    DIST_ENUM returns the number of distributions of indistinguishable objects.
        // 
        //  Licensing:
        // 
        //    This code is distributed under the GNU LGPL license.
        // 
        //  Modified:
        // 
        //    24 July 2011
        // 
        //  Author:
        // 
        //    John Burkardt
        // 
        //  Parameters:
        // 
        //    Input, int K, the number of distinguishable "slots".
        // 
        //    Input, int M, the number of indistinguishable objects.
        // 
        //    Output, int DIST_ENUM, the number of distributions of M
        //    indistinguishable objects about K distinguishable slots.
        // 
    {
        int value = typeMethods.i4_choose(m + k - 1, m);

        return value;
    }

    public static void dist_next(int k, int m, ref int[] q, ref int leftmost, ref bool more )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    DIST_NEXT returns the next distribution of indistinguishable objects.
        // 
        //  Discussion:
        // 
        //    A distribution of M objects into K parts is an ordered sequence
        //    of K nonnegative integers which sum to M.  This is similar to
        //    a partition of a set into K subsets, except that here the order
        //    matters.  That is, (1,1,2) and (1,2,1) are considered to be
        //    different distributions.
        // 
        //    On the first call to this routine, the user should set MORE = FALSE,
        //    to signal that this is a startup for the given computation.  The routine
        //    will return the first distribution, and set MORE = TRUE.
        // 
        //    If the user calls again, with MORE = TRUE, the next distribution
        //    is being requested.  If the routine returns with MORE = TRUE, then
        //    that distribution was found and returned.  However, if the routine
        //    returns with MORE = FALSE, then no more distributions were found;
        //    the enumeration of distributions has terminated.
        // 
        //    A "distribution of M indistinguishable objects into K slots" is
        //    sometimes called a "composition of the integer M into K parts".
        // 
        //  Example:
        // 
        //    K = 3, M = 5
        // 
        //    0           0           5
        //    0           1           4
        //    0           2           3
        //    0           3           2
        //    0           4           1
        //    0           5           0
        //    1           0           4
        //    1           1           3
        //    1           2           2
        //    1           3           1
        //    1           4           0
        //    2           0           3
        //    2           1           2
        //    2           2           1
        //    2           3           0
        //    3           0           2
        //    3           1           1
        //    3           2           0
        //    4           0           1
        //    4           1           0
        //    5           0           0
        // 
        //  Licensing:
        // 
        //    This code is distributed under the GNU LGPL license.
        // 
        //  Modified:
        // 
        //    24 December 2015
        // 
        //  Author:
        // 
        //    John Burkardt
        // 
        //  Reference:
        // 
        //    Robert Fenichel,
        //    Algorithm 329:
        //    Distribution of Indistinguishable Objects into
        //    Distinguishable Slots,
        //    Communications of the ACM,
        //    Volume 11, Number 6, June 1968, page 430.
        // 
        //  Parameters:
        // 
        //    Input, int K, the number of distinguishable "slots".
        // 
        //    Input, int M, the number of indistinguishable objects.
        // 
        //    Input/output, int Q[K], the number of objects in each
        //    slot.
        // 
        //    Input/output, int &LEFTMOST, used to speed up the calculation.
        //    Set this to 0 before the first call.
        //
        //    Input/output, bool &MORE, used by the user to start the computation,
        //    and by the routine to stop the computation.
        // 
    {
        int i;
        switch (more)
        {
            // 
            //  The startup call.
            // 
            case false:
            {
                more = true;
                for (i = 0; i < k - 1; i++)
                {
                    q[i] = 0;
                }

                q[k - 1] = m;

                leftmost = k + 1;
                break;
            }
            // 
            default:
            {
                if (q[0] == m)
                {
                    more = false;

                    for (i = 0; i < k - 1; i++)
                    {
                        q[i] = 0;
                    }

                    q[k - 1] = m;

                    leftmost = k + 1;
                }
                else if (leftmost < k + 1)
                {
                    leftmost -= 1;
                    q[k - 1] = q[leftmost - 1] - 1;
                    q[leftmost - 1] = 0;
                    q[leftmost - 2] += 1;
                    if (q[k - 1] != 0)
                    {
                        leftmost = k + 1;
                    }
                }
                else
                {
                    leftmost = q[k - 1] switch
                    {
                        1 => k,
                        _ => leftmost
                    };

                    q[k - 1] -= 1;
                    q[k - 2] += 1;
                }

                break;
            }
        }
    }
}