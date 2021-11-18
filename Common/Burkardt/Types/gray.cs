using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static bool gray_code_check(int n, int[] t)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    GRAY_CODE_CHECK checks a Gray code element.
        // 
        //  Licensing:
        // 
        //    This code is distributed under the GNU LGPL license.
        // 
        //  Modified:
        // 
        //    25 November 2015
        // 
        //  Author:
        // 
        //    John Burkardt
        // 
        //  Parameters:
        // 
        //    Input, int N, the number of digits in each element.
        //    N must be positive.
        // 
        //    Input, int T[N], an element of the Gray code.
        //    Each entry T(I) is either 0 or 1.
        // 
        //    Output, bool GRAY_CODE_CHECK.
        //    TRUE, the data is legal.
        //    FALSE, the data is not legal.
        // 
    {
        int i;

        bool check = true;

        switch (n)
        {
            case < 1:
                check = false;
                return check;
        }

        for (i = 0; i < n; i++)
        {
            if (t[i] == 0 || t[i] == 1)
            {
                continue;
            }

            check = false;
            return check;
        }

        return check;
    }

    public static int gray_code_enum(int n)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    GRAY_CODE_ENUM enumerates the Gray codes on N digits.
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
        //    Input, int N, the number of digits in each element.
        //    N must be nonnegative.
        // 
        //    Output, int GRAY_CODE_ENUM, the number of distinct elements.
        // 
    {
        int value = (int) Math.Pow(2, n);

        return value;
    }

    public static void gray_code_successor(int n, ref int[] t, ref int rank)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    GRAY_CODE_SUCCESSOR computes the binary reflected Gray code successor.
        // 
        //  Example:
        // 
        //    000, 001, 011, 010, 110, 111, 101, 100,
        //    after which the sequence repeats.
        // 
        //  Discussion:
        // 
        //    In the original code, the successor of the element that has an
        //    initial 1 followed by N-1 zeroes is undefined.  In this version,
        //    the successor is the element with N zeroes.
        // 
        //  Licensing:
        // 
        //    This code is distributed under the GNU LGPL license.
        // 
        //  Modified:
        // 
        //    26 July 2011
        // 
        //  Author:
        // 
        //    John Burkardt
        // 
        //  Reference:
        // 
        //    Donald Kreher, Douglas Simpson,
        //    Combinatorial Algorithms,
        //    CRC Press, 1998,
        //    ISBN: 0-8493-3988-X,
        //    LC: QA164.K73.
        // 
        //  Parameters:
        // 
        //    Input, int N, the number of digits in each element.
        //    N must be positive.
        // 
        //    Input/output, int T[N].
        //    On input, T contains an element of the Gray code, that is,
        //    each entry T(I) is either 0 or 1.
        //    On output, T contains the successor to the input value; this
        //    is an element of the Gray code, which differs from the input
        //    value in a single position.
        // 
        //    Input/output, int &RANK, the rank.
        //    If RANK = -1 on input, then the routine understands that this is
        //    the first call, and that the user wishes the routine to supply
        //    the first element in the ordering, which has RANK = 0.
        //    In general, the input value of RANK is increased by 1 for output,
        //    unless the very last element of the ordering was input, in which
        //    case the output value of RANK is 0.
        // 
    {
        int i;
        switch (rank)
        {
            // 
            //  Return the first element.
            // 
            case -1:
            {
                for (i = 0; i < n; i++)
                {
                    t[i] = 0;
                }

                rank = 0;
                return;
            }
        }

        // 
        //  Check.
        // 
        bool check = gray_code_check(n, t);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("GRAY_CODE_SUCCESSOR - Fatal error!");
                Console.WriteLine("  The input array is illegal.");
                return;
        }

        int weight = i4vec_sum(n, t);

        switch (weight % 2)
        {
            case 0:
            {
                t[n - 1] = t[n - 1] switch
                {
                    0 => 1,
                    _ => 0
                };

                rank += 1;
                break;
            }
            default:
            {
                for (i = n - 1; 1 <= i; i--)
                {
                    switch (t[i])
                    {
                        case 1:
                        {
                            t[i - 1] = t[i - 1] switch
                            {
                                0 => 1,
                                _ => 0
                            };

                            rank += 1;
                            return;
                        }
                    }
                }

                // 
                //  The final element was input.
                //  Return the first element.
                // 
                for (i = 0; i < n; i++)
                {
                    t[i] = 0;
                }

                rank = 0;
                break;
            }
        }
    }

    public class GrayData
    {
        public int k;
        public int n_save = -1;
        public int[] a;
    }

    public static void gray_next_(ref GrayData data, int n, ref int change)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRAY_NEXT generates the next Gray code by switching one item at a time.
        //
        //  Discussion:
        //
        //    On the first call only, the user must set CHANGE = -N.
        //    This initializes the routine to the Gray code for N zeroes.
        //
        //    Each time it is called thereafter, it returns in CHANGE the index
        //    of the item to be switched in the Gray code.  The sign of CHANGE
        //    indicates whether the item is to be added or subtracted (or
        //    whether the corresponding bit should become 1 or 0).  When
        //    CHANGE is equal to N+1 on output, all the Gray codes have been
        //    generated.
        //
        //    The routine has internal memory that is set up on call with
        //    CHANGE = -N, and released on final return.
        //
        //  Example:
        //
        //    N  CHANGE         Subset in/out   Binary Number
        //                      Interpretation  Interpretation
        //                       1 2 4 8
        //   --  ---------      --------------  --------------
        //
        //    4   -4 / 0         0 0 0 0         0
        //
        //        +1             1 0 0 0         1
        //          +2           1 1 0 0         3
        //        -1             0 1 0 0         2
        //            +3         0 1 1 0         6
        //        +1             1 1 1 0         7
        //          -2           1 0 1 0         5
        //        -1             0 0 1 0         4
        //              +4       0 0 1 1        12
        //        +1             1 0 1 1        13
        //          +2           1 1 1 1        15
        //        -1             0 1 1 1        14
        //            -3         0 1 0 1        10
        //        +1             1 1 0 1        11
        //          -2           1 0 0 1         9
        //        -1             0 0 0 1         8
        //              -4       0 0 0 0         0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Albert Nijenhuis, Herbert Wilf,
        //    Combinatorial Algorithms for Computers and Calculators,
        //    Second Edition,
        //    Academic Press, 1978,
        //    ISBN: 0-12-519260-6,
        //    LC: QA164.N54.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the total set from which
        //    subsets will be drawn.
        //
        //    Input/output, int *CHANGE.  This item is used for input only
        //    on the first call for a particular sequence of Gray codes,
        //    at which time it must be set to -N.  This corresponds to
        //    all items being excluded, or all bits being 0, in the Gray code.
        //    On output, CHANGE indicates which of the N items must be "changed", 
        //    and the sign indicates whether the item is to be added or removed
        //    (or the bit is to become 1 or 0).  Note that on return from the 
        //    first call, CHANGE is set to 0, indicating that we begin with
        //    the empty set.
        //
    {
        int i;

        switch (n)
        {
            case <= 0:
                Console.WriteLine("");
                Console.WriteLine("GRAY_NEXT - Fatal error!");
                Console.WriteLine("  Input value of N <= 0.");
                return;
        }

        if (change == -n)
        {

            data.a = new int[n];
            for (i = 0; i < n; i++)
            {
                data.a[i] = 0;
            }

            data.n_save = n;
            data.k = 1;
            change = 0;

            return;
        }

        if (n != data.n_save)
        {
            Console.WriteLine("");
            Console.WriteLine("GRAY_NEXT - Fatal error!");
            Console.WriteLine("  Input value of N has changed from definition value.");
        }

        switch (data.k % 2)
        {
            //
            //  First determine WHICH item is to be changed.
            //
            case 1:
                change = 1;
                break;
            default:
            {
                for (i = 1; i <= data.n_save; i++)
                {
                    if (data.a[i - 1] != 1)
                    {
                        continue;
                    }

                    change = i + 1;
                    break;
                }

                break;
            }
        }

        //
        //  Take care of the terminal case CHANGE = N + 1.
        //
        if (change == n + 1)
        {
            change = n;
        }

        switch (data.a[change - 1])
        {
            //
            //  Now determine HOW the item is to be changed.
            //
            case 0:
                data.a[change - 1] = 1;
                break;
            case 1:
                data.a[change - 1] = 0;
                change = -change;
                break;
        }

        //
        //  Update the counter.
        //
        data.k += 1;
        //
        //  If the output CHANGE = -N_SAVE, then we're done.
        //
        if (change != -data.n_save)
        {
            return;
        }

        data.a = null;
        data.n_save = 0;
        data.k = 0;

    }

    public static void gray_next(ref GrayData data, int n, ref int change)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRAY_NEXT generates the next Gray code by switching one item at a time.
        //
        //  Discussion:
        //
        //    On the first call only, the user must set CHANGE = -N.
        //    This initializes the routine to the Gray code for N zeroes.
        //
        //    Each time it is called thereafter, it returns in CHANGE the index
        //    of the item to be switched in the Gray code.  The sign of CHANGE
        //    indicates whether the item is to be added or subtracted (or
        //    whether the corresponding bit should become 1 or 0).  When
        //    CHANGE is equal to N+1 on output, all the Gray codes have been
        //    generated.
        //
        //  Example:
        //
        //    N  CHANGE         Subset in/out   Binary Number
        //                      Interpretation  Interpretation
        //                       1 2 4 8
        //   --  ---------      --------------  --------------
        //
        //    4   -4 / 0         0 0 0 0         0
        //
        //        +1             1 0 0 0         1
        //          +2           1 1 0 0         3
        //        -1             0 1 0 0         2
        //            +3         0 1 1 0         6
        //        +1             1 1 1 0         7
        //          -2           1 0 1 0         5
        //        -1             0 0 1 0         4
        //              +4       0 0 1 1        12
        //        +1             1 0 1 1        13
        //          +2           1 1 1 1        15
        //        -1             0 1 1 1        14
        //            -3         0 1 0 1        10
        //        +1             1 1 0 1        11
        //          -2           1 0 0 1         9
        //        -1             0 0 0 1         8
        //              -4       0 0 0 0         0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 June 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Albert Nijenhuis, Herbert Wilf,
        //    Combinatorial Algorithms for Computers and Calculators,
        //    Second Edition,
        //    Academic Press, 1978,
        //    ISBN: 0-12-519260-6,
        //    LC: QA164.N54.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the total set from which
        //    subsets will be drawn.
        //
        //    Input/output, int &CHANGE.  This item is used for input only
        //    on the first call for a particular sequence of Gray codes,
        //    at which time it must be set to -N.  This corresponds to
        //    all items being excluded, or all bits being 0, in the Gray code.
        //    On output, CHANGE indicates which of the N items must be "changed", 
        //    and the sign indicates whether the item is to be added or removed
        //    (or the bit is to become 1 or 0).  Note that on return from the 
        //    first call, CHANGE is set to 0, indicating that we begin with
        //    the empty set.
        //
        //    Input/output, int &K, a bookkeeping variable.
        //    The user must declare this variable before the first call.
        //    The output value from one call should be the input value for the next.
        //    The user should not change this variable.
        //
        //    Input/output, int A[N], a bookkeeping variable.
        //    The user must declare this variable before the first call.
        //    The output value from one call should be the input value for the next.
        //    The user should not change this variable.
        //
    {
        int i;

        switch (n)
        {
            case <= 0:
                Console.WriteLine("");
                Console.WriteLine("GRAY_NEXT - Fatal error!");
                Console.WriteLine("  Input value of N <= 0.");
                return;
        }

        if (change == -n)
        {
            change = 0;
            data.k = 1;
            for (i = 0; i < n; i++)
            {
                data.a[i] = 0;
            }

            return;
        }

        switch (data.k % 2)
        {
            //
            //  First determine WHICH item is to be changed.
            //
            case 1:
                change = 1;
                break;
            default:
            {
                for (i = 1; i <= n; i++)
                {
                    if (data.a[i - 1] != 1)
                    {
                        continue;
                    }

                    change = i + 1;
                    break;
                }

                break;
            }
        }

        //
        //  Take care of the terminal case CHANGE = N + 1.
        //
        if (change == n + 1)
        {
            change = n;
        }

        switch (data.a[change - 1])
        {
            //
            //  Now determine HOW the item is to be changed.
            //
            case 0:
                data.a[change - 1] = 1;
                break;
            case 1:
                data.a[change - 1] = 0;
                change = -change;
                break;
        }

        //
        //  Update the counter.
        //
        data.k += 1;
        //
        //  If the output CHANGE = -N, then we're done.
        //
        if (change == -n)
        {
            data.k = 0;
        }

    }

}