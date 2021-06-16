using System;
using Burkardt.Types;

namespace Burkardt.RankingNS
{
    public static partial class Ranking
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

            if (n < 1)
            {
                check = false;
                return check;
            }

            for (i = 0; i < n; i++)
            {
                if (t[i] != 0 && t[i] != 1)
                {
                    check = false;
                    return check;
                }
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
            int value = (int)Math.Pow(2, n);

            return value;
        }

        public static int gray_code_rank(int n, int[] t)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    GRAY_CODE_RANK computes the rank of a Gray code element.
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
            //    Input, int T[N], an element of the Gray code.
            //    Each entry is either 0 or 1.
            // 
            //    Output, int GRAY_CODE_RANK, the rank of the element.
            // 
        {
            int b;
            bool check;
            int i;
            int rank;
            // 
            //  Check.
            // 
            check = gray_code_check(n, t);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("GRAY_CODE_RANK - Fatal error!");
                Console.WriteLine("  The input array is illegal.");
                return (1);
            }

            rank = 0;
            b = 0;

            for (i = n - 1; 0 <= i; i--)
            {
                if (t[n - i - 1] != 0)
                {
                    b = 1 - b;
                }

                if (b == 1)
                {
                    rank = rank + (int)Math.Pow(2, i);
                }
            }

            return rank;
        }

        public static void gray_code_successor(int n, ref int[] t, ref int rank )

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
            bool check;
            int i;
            int weight;
            // 
            //  Return the first element.
            // 
            if (rank == -1)
            {
                for (i = 0; i < n; i++)
                {
                    t[i] = 0;
                }

                rank = 0;
                return;
            }

            // 
            //  Check.
            // 
            check = gray_code_check(n, t);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("GRAY_CODE_SUCCESSOR - Fatal error!");
                Console.WriteLine("  The input array is illegal.");
                return;
            }

            weight = typeMethods.i4vec_sum(n, t);

            if ((weight % 2) == 0)
            {
                if (t[n - 1] == 0)
                {
                    t[n - 1] = 1;
                }
                else
                {
                    t[n - 1] = 0;
                }

                rank = rank + 1;
                return;
            }
            else
            {
                for (i = n - 1; 1 <= i; i--)
                {
                    if (t[i] == 1)
                    {
                        if (t[i - 1] == 0)
                        {
                            t[i - 1] = 1;
                        }
                        else
                        {
                            t[i - 1] = 0;
                        }

                        rank = rank + 1;
                        return;
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
            }
        }

        public static int[] gray_code_unrank(int rank, int n)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    GRAY_CODE_UNRANK computes the Gray code element of given rank.
            // 
            //  Licensing:
            // 
            //    This code is distributed under the GNU LGPL license.
            // 
            //  Modified:
            // 
            //    25 July 2011
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
            //    Input, int RANK, the rank of the element.
            //    0 <= RANK <= 2^N.
            // 
            //    Input, int N, the number of digits in each element.
            //    N must be positive.
            // 
            //    Output, int GRAY_CODE_UNRANK[N], the element of the Gray code which has
            //    the given rank.
            // 
        {
            int b;
            int bprime;
            int i;
            int ngray;
            int rank_copy;
            int[] t;
            // 
            //  Check.
            // 
            if (n < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("GRAY_CODE_UNRANK - Fatal error!");
                Console.WriteLine("  Input N is illegal.");
                return null;
            }

            ngray = gray_code_enum(n);

            if (rank < 0 || ngray < rank)
            {
                Console.WriteLine("");
                Console.WriteLine("GRAY_CODE_UNRANK - Fatal error!");
                Console.WriteLine("  The input rank is illegal.");
                return null;
            }

            t = new int[n];

            rank_copy = rank;
            for (i = 0; i < n; i++)
            {
                t[i] = 0;
            }

            bprime = 0;

            for (i = n - 1; 0 <= i; i--)
            {
                b = rank_copy / (int)Math.Pow(2, i);

                if (b != bprime)
                {
                    t[n - i - 1] = 1;
                }

                bprime = b;
                rank_copy = rank_copy - b * (int)Math.Pow(2, i);
            }

            return t;
        }
    }
}