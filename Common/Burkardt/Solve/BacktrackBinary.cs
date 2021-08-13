namespace Burkardt.SolveNS
{
    public static class BacktrackBinary
    {
        public static void backbin_rc(int n, bool reject, ref int n2, ref int[] choice)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BACKBIN_RC uses reverse communication for binary backtracking.
            //
            //  Discussion:
            //
            //    If this procedure returns a solution with N2 = N, which is acceptable
            //    to the user, then a full solution has been found.
            //
            //    If this procedure returns N2 = -1, no more potential solutions are
            //    available to consider.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 January 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the length of the full solution.
            //
            //    Input, bool REJECT, is TRUE if the proposed partial solution
            //    in the first N2 entries of CHOICE must be rejected.
            //
            //    Input/output, int &N2, the length of the current
            //    partial solution.  On first call for a given problem, the user
            //    should set N2 to -1.  If the program has exhausted the search space,
            //    the value of N2 will be returned as -1.
            //
            //    Input/output, int CHOICE[N], indicates the current
            //    partial solution in entries 1 through N2, which will contain 0 or 1.
            //
        {
            //
            //  N2 = -1 means an initialization call.
            //
            if (n2 == -1)
            {
                for (int i = 0; i < n; i++)
                {
                    choice[i] = -1;
                }

                n2 = 1;
                choice[n2 - 1] = 1;
            }
            //
            //  1 <= FOCUS means we asked the user to evaluate CHOICE(1:N2).
            //
            //  N2 = N means we returned a full prospective solution
            //  so in any case we must increment CHOICE.
            //
            //  Returning REJECT = 1 means no solution begins this way
            //  so we must increment CHOICE.
            //
            else if (n2 == n || reject)
            {
                while (1 < n2)
                {
                    if (choice[n2 - 1] == 1)
                    {
                        choice[n2 - 1] = 0;
                        break;
                    }

                    choice[n2 - 1] = -1;
                    n2 = n2 - 1;
                }

                //
                //  Have we exhausted the solution space?
                //
                if (n2 == 1)
                {
                    if (choice[n2 - 1] == 1)
                    {
                        choice[n2 - 1] = 0;
                    }
                    else
                    {
                        choice[n2 - 1] = -1;
                        n2 = -1;
                    }
                }
            }
            //
            //  N2 < N and not REJECT means we can increment N2.
            //
            else
            {
                n2 = n2 + 1;
                choice[n2 - 1] = 1;
            }
        }
    }
}