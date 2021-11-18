using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static i8 s_to_i8 ( string s )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    S_TO_I8 reads an I8 from a string.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string S, a string to be examined.
        //
        //    Output, int *LAST, the last character of S used to make IVAL.
        //
        //    Output, bool *ERROR is TRUE if an error occurred.
        //
        //    Output, int *S_TO_I4, the integer value read from the string.
        //    If the string is blank, then IVAL will be returned 0.
        //
    {
        i8 ival = new() {val = 0};

        int istate = 0;
        long isgn = 1;
        int i = 0;

        while ( i < s.Length )
        {
            char c = s[i];
            i += 1;
            switch (istate)
            {
                //
                //  Haven't read anything.
                //
                case 0 when c == ' ':
                    break;
                case 0 when c == '-':
                    istate = 1;
                    isgn = -1;
                    break;
                case 0 when c == '+':
                    istate = 1;
                    isgn = + 1;
                    break;
                case 0 when c is >= '0' and <= '9':
                    istate = 2;
                    ival.val = c - '0';
                    break;
                case 0:
                    ival.error = true;
                    return ival;
                //
                //  Have read the sign, expecting digits.
                //
                case 1 when c == ' ':
                    break;
                case 1 when c is >= '0' and <= '9':
                    istate = 2;
                    ival.val = c - '0';
                    break;
                case 1:
                    ival.error = true;
                    return ival;
                //
                //  Have read at least one digit, expecting more.
                //
                case 2 when c is >= '0' and <= '9':
                    ival.val = 10 * ival.val + c - '0';
                    break;
                case 2:
                    ival.val = isgn * ival.val;
                    ival.lchar = i - 1;
                    return ival;
            }
        }

        switch (istate)
        {
            //
            //  If we read all the characters in the string, see if we're OK.
            //
            case 2:
                ival.val = isgn * ival.val;
                ival.lchar = s_len_trim ( s );
                break;
            default:
                ival.error = true;
                ival.lchar = 0;
                break;
        }

        return ival;
    }

    public static int i8_bit_hi1(long n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I8_BIT_HI1 returns the position of the high 1 bit base 2 in an integer.
        //
        //  Example:
        //
        //       N    Binary    Hi 1
        //    ----    --------  ----
        //       0           0     0
        //       1           1     1
        //       2          10     2
        //       3          11     2 
        //       4         100     3
        //       5         101     3
        //       6         110     3
        //       7         111     3
        //       8        1000     4
        //       9        1001     4
        //      10        1010     4
        //      11        1011     4
        //      12        1100     4
        //      13        1101     4
        //      14        1110     4
        //      15        1111     4
        //      16       10000     5
        //      17       10001     5
        //    1023  1111111111    10
        //    1024 10000000000    11
        //    1025 10000000001    11
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, long long int N, the integer to be measured.
        //    N should be nonnegative.  If N is nonpositive, I8_BIT_HI1
        //    will always be 0.
        //
        //    Output, int I8_BIT_HI1, the number of bits base 2.
        //
    {
        int bit = 0;

        while (0 < n)
        {
            bit += 1;
            n /= 2;
        }

        return bit;
    }

    public static int i8_bit_lo0(long n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I8_BIT_LO0 returns the position of the low 0 bit base 2 in an integer.
        //
        //  Example:
        //
        //       N    Binary    Lo 0
        //    ----    --------  ----
        //       0           0     1
        //       1           1     2
        //       2          10     1
        //       3          11     3 
        //       4         100     1
        //       5         101     2
        //       6         110     1
        //       7         111     4
        //       8        1000     1
        //       9        1001     2
        //      10        1010     1
        //      11        1011     3
        //      12        1100     1
        //      13        1101     2
        //      14        1110     1
        //      15        1111     5
        //      16       10000     1
        //      17       10001     2
        //    1023  1111111111    11
        //    1024 10000000000     1
        //    1025 10000000001     2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, long long int N, the integer to be measured.
        //    N should be nonnegative.
        //
        //    Output, int I8_BIT_LO0, the position of the low 1 bit.
        //
    {
        int bit = 0;

        while (true)
        {
            bit += 1;
            long n2 = n / 2;

            if (n == 2 * n2)
            {
                break;
            }

            n = n2;

        }

        return bit;
    }

    public static long i8_choose(long n, long k )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I8_CHOOSE computes the binomial coefficient C(N,K) as an I8.
        //
        //  Discussion:
        //
        //    The value is calculated in such a way as to avoid overflow and
        //    roundoff.
        //
        //    The formula used is:
        //
        //      C(N,K) = N! / ( K! * (N-K)! )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    ML Wolfson, HV Wright,
        //    Algorithm 160:
        //    Combinatorial of M Things Taken N at a Time,
        //    Communications of the ACM,
        //    Volume 6, Number 4, April 1963, page 161.
        //
        //  Parameters:
        //
        //    Input, long long int N, K, the values of N and K.
        //
        //    Output, long long int I8_CHOOSE, the number of combinations of N
        //    things taken K at a time.
        //
    {
        long value;

        long mn = Math.Min(k, n - k);

        switch (mn)
        {
            case < 0:
                value = 0;
                break;
            case 0:
                value = 1;
                break;
            default:
            {
                long mx = Math.Max(k, n - k);
                value = mx + 1;

                for (long i = 2; i <= mn; i++)
                {
                    value = value * (mx + i) / i;
                }

                break;
            }
        }

        return value;
    }

    public static long i8_huge()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I8_HUGE returns a "huge" I8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, long long int I8_HUGE, a "huge" I8.
        //
    {
        const long value = 9223372036854775807;

        return value;
    }

    public static double i8_huge_normalizer()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I8_HUGE_NORMALIZER returns the "normalizer" for I8_HUGE.
        //
        //  Discussion:
        //
        //    The value returned is 1 / ( I8_HUGE + 1 ).
        //
        //    For any I8, it should be the case that
        //
        //     -1 < I8 * I8_HUGE_NORMALIZER < 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double I8_HUGE_NORMALIZER, the "normalizer"
        //    for I8_HUGE.
        //
    {
        const double value = 1.084202172485504434007E-19;

        return value;
    }
        
}