using System;

namespace Burkardt.Sequence;

public static class Thue
{
    public static void thue_binary_next(ref int n, ref int[] thue )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    THUE_BINARY_NEXT returns the next element in a binary Thue sequence.
        //
        //  Discussion:
        //
        //    Thue demonstrated that arbitrarily long sequences of 0's and
        //    1's could be generated which had the "cubefree" property.  In
        //    other words, for a given string S, there was no substring W
        //    such that S contained "WWW".  In fact, a stronger result holds:
        //    if "a" is the first letter of W, it is never the case that S
        //    contains the substring "WWa".
        //
        //    In this example, the digits allowed are binary, that is, just
        //    "0" and "1".  The replacement rules are:
        //
        //    "0" -> "01"
        //    "1" -> "10"
        //
        //    This routine produces the next binary Thue sequence in a given series.
        //    However, the input sequence must be a Thue sequence in order for
        //    us to guarantee that the output sequence will also have the
        //    cubic nonrepetition property.
        //
        //    Also, enough space must be set aside in THUE to hold the
        //    output sequence.  This will always be twice the input
        //    value of N.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, int &N.  On input, the length of the input sequence.
        //    On output, the length of the output sequence.
        //
        //    Input, int THUE[N].  On input, the initial Thue sequence, and on
        //    output, the result of applying the substitution rules once.
        //
    {
        int i;

        int n_out = 0;
        int[] thue_out = new int[2 * n];

        for (i = 0; i < n; i++)
        {
            switch (thue[i])
            {
                case 0:
                    thue_out[n_out] = 0;
                    n_out += 1;
                    thue_out[n_out] = 1;
                    n_out += 1;
                    break;
                case 1:
                    thue_out[n_out] = 1;
                    n_out += 1;
                    thue_out[n_out] = 0;
                    n_out += 1;
                    break;
                default:
                    Console.WriteLine("");
                    Console.WriteLine("THUE_BINARY_NEXT - Fatal error!");
                    Console.WriteLine("  The input sequence contains a non-binary digit");
                    Console.WriteLine("  THUE[" + i + "] = " + thue[i] + "");
                    return;
            }
        }

        n = n_out;

        for (i = 0; i < n; i++)
        {
            thue[i] = thue_out[i];
        }
    }

    public static void thue_ternary_next(ref int n, ref int[] thue )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    THUE_TERNARY_NEXT returns the next element in a ternary Thue sequence.
        //
        //  Discussion:
        //
        //    Thue was interested in showing that there were arbitrarily long
        //    sequences of digits which never displayed a pair of contiguous
        //    repetitions of any length.  That is, there was no occurrence of
        //    "00" or "1010" or "121121", anywhere in the string.  This makes
        //    the string "squarefree".
        //
        //    To do this, he demonstrated a way to start with a single digit,
        //    and to repeatedly apply a series of transformation rules to each
        //    digit of the sequence, deriving nonrepeating sequences of ever
        //    greater length.
        //
        //    In this example, the digits allowed are ternary, that is, just
        //    "0", "1" and "2".  The replacement rules are:
        //
        //    "0" -> "12"
        //    "1" -> "102"
        //    "2" -> "0"
        //
        //    This routine produces the next Thue sequence in a given series.
        //    However, the input sequence must be a Thue sequence in order for
        //    us to guarantee that the output sequence will also have the
        //    nonrepetition property.
        //
        //    Also, enough space must be set aside in THUE to hold the
        //    output sequence.  This will never be more than 3 times the input
        //    value of N.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Brian Hayes,
        //    Third Base,
        //    American Scientist,
        //    Volume 89, Number 6, pages 490-494, November-December 2001.
        //
        //  Parameters:
        //
        //    Input/output, int &N.  On input, the length of the input sequence.
        //    On output, the length of the output sequence.
        //
        //    Input, int THUE[*N].  On input, the initial Thue sequence, and on
        //    output, the result of applying the substitution rules once.
        //
    {
        int i;

        int n_out = 0;
        int[] thue_out = new int[3 * n];

        for (i = 0; i < n; i++)
        {
            switch (thue[i])
            {
                case 0:
                    thue_out[n_out] = 1;
                    n_out += 1;
                    thue_out[n_out] = 2;
                    n_out += 1;
                    break;
                case 1:
                    thue_out[n_out] = 1;
                    n_out += 1;
                    thue_out[n_out] = 0;
                    n_out += 1;
                    thue_out[n_out] = 2;
                    n_out += 1;
                    break;
                case 2:
                    thue_out[n_out] = 0;
                    n_out += 1;
                    break;
                default:
                    Console.WriteLine("");
                    Console.WriteLine("THUE_TERNARY_NEXT - Fatal error!");
                    Console.WriteLine("  The input sequence contains a non-ternary digit");
                    Console.WriteLine("  THUE[" + i + "] = " + thue[i] + "");
                    return;
            }
        }

        n = n_out;
        for (i = 0; i < n_out; i++)
        {
            thue[i] = thue_out[i];
        }
    }
}