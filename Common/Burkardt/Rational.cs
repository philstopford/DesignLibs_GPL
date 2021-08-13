using System;
using Burkardt.Types;

namespace Burkardt.RationalNS
{
    public static class Rational
    {
        public static void rat_add(int itop1, int ibot1, int itop2, int ibot2, ref int itop, ref int ibot,
                ref bool error)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    RAT_ADD adds two rational values.
            //
            //  Discussion:
            //
            //    The routine computes
            //
            //      ITOP/IBOT = ITOP1/IBOT1 + ITOP2/IBOT2
            //
            //    while trying to avoid int overflow.
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
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int ITOP1, IBOT1, the first value to add.
            //
            //    Input, int ITOP2, IBOT2, the second value to add.
            //
            //    Output, ref int ITOP, &IBOT, the sum.
            //
            //    Output, ref bool ERROR, is TRUE if an error occurred.
            //
        {
            int i_max;
            const int i4_huge = 2147483647;
            int itemp;
            int ibot3;

            i_max = i4_huge;

            error = false;

            if (itop1 == 0)
            {
                itop = itop2;
                ibot = ibot2;
                return;
            }
            else if (itop2 == 0)
            {
                itop = itop1;
                ibot = ibot1;
                return;
            }

            //
            //  Divide out  the greatest common factor of the two denominators.
            //
            ibot3 = typeMethods.i4_gcd(ibot1, ibot2);
            ibot1 = ibot1 / ibot3;
            ibot2 = ibot2 / ibot3;
            //
            //  The fraction may now be formally written as:
            //
            //    (itop1*ibot2 + itop2*ibot1) / (ibot1*ibot2*ibot3)
            //
            //  Check the tops for overflow.
            //
            if ((double)(i_max) < Math.Abs((double)(itop1) * (double)(ibot2)))
            {
                error = true;
                Console.WriteLine("");
                Console.WriteLine("RAT_ADD - Fatal error!");
                Console.WriteLine("  Overflow of top of rational sum.");
                itop = 0;
                ibot = 1;
                return;
            }

            itop1 = itop1 * ibot2;

            if ((double)(i_max) < Math.Abs((double)(itop2) * (double)(ibot1)))
            {
                error = true;
                Console.WriteLine("");
                Console.WriteLine("RAT_ADD - Fatal error!");
                Console.WriteLine("  Overflow of top of rational sum.");
                itop = 0;
                ibot = 1;
                return;
            }

            itop2 = itop2 * ibot1;

            if ((double)(i_max) < Math.Abs((double)(itop1) + (double)(itop2)))
            {
                error = true;
                Console.WriteLine("");
                Console.WriteLine("RAT_ADD - Fatal error!");
                Console.WriteLine("  Overflow of top of rational sum.");
                itop = 0;
                ibot = 1;
                return;
            }

            itop = itop1 + itop2;
            //
            //  Check the bottom for overflow.
            //
            if ((double)(i_max) <
                Math.Abs((double)(ibot1) * (double)(ibot2) * (double)(ibot3)))
            {
                error = true;
                Console.WriteLine("");
                Console.WriteLine("RAT_ADD - Fatal error!");
                Console.WriteLine("  Overflow of bottom of rational sum.");
                itop = 0;
                ibot = 1;
                return;
            }

            ibot = ibot1 * ibot2 * ibot3;
            //
            //  Put the fraction in lowest terms.
            //
            itemp = typeMethods.i4_gcd(itop, ibot);
            itop = itop / itemp;
            ibot = ibot / itemp;
            //
            //  The bottom should be positive.
            //
            if (ibot < 0)
            {
                ibot = -ibot;
                itop = -itop;
            }

        }

        public static void rat_div(int itop1, int ibot1, int itop2, int ibot2, ref int itop,
                ref int ibot, ref bool error)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    RAT_DIV divides one rational value by another.
            //
            //  Discussion:
            //
            //    The routine computes
            //
            //      ITOP / IBOT = ( ITOP1 / IBOT1 ) / ( ITOP2 / IBOT2 ).
            //
            //    while avoiding integer overflow.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 July 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int ITOP1, IBOT1, the numerator.
            //
            //    Input, int ITOP2, IBOT2, the denominator.
            //
            //    Output, ref int ITOP, &IBOT, the result.
            //
            //    Output, ref bool ERROR, is TRUE if an error occurred.
            //
        {
            int i_max;
            const int i4_huge = 2147483647;
            int itemp;

            error = false;

            i_max = i4_huge;

            if (ibot1 == 0 || itop2 == 0 || ibot2 == 0)
            {
                error = true;
                return;
            }

            if (itop1 == 0)
            {
                itop = 0;
                ibot = 1;
                return;
            }

            //
            //  Get rid of all common factors in top and bottom.
            //
            itemp = typeMethods.i4_gcd(itop1, ibot1);
            itop1 = itop1 / itemp;
            ibot1 = ibot1 / itemp;
            itemp = typeMethods.i4_gcd(itop1, itop2);
            itop1 = itop1 / itemp;
            itop2 = itop2 / itemp;
            itemp = typeMethods.i4_gcd(ibot2, ibot1);
            ibot2 = ibot2 / itemp;
            ibot1 = ibot1 / itemp;
            itemp = typeMethods.i4_gcd(ibot2, itop2);
            ibot2 = ibot2 / itemp;
            itop2 = itop2 / itemp;
            //
            //  The fraction (ITOP1*IBOT2)/(IBOT1*ITOP2) is in lowest terms.
            //
            //  Check the top for overflow.
            //
            if ((double)i_max < Math.Abs((double)itop1 * (double)ibot2))
            {
                error = true;
                Console.WriteLine("");
                Console.WriteLine("RAT_DIV - Fatal error!");
                Console.WriteLine("  Overflow of top of rational fraction.");
                itop = 0;
                ibot = 0;
                return;
            }

            itop = itop1 * ibot2;
            //
            //  Check the bottom IBOT1*ITOP2 for overflow.
            //
            if ((double)i_max < Math.Abs((double)ibot1 * (double)itop2))
            {
                error = true;
                Console.WriteLine("");
                Console.WriteLine("RAT_DIV - Fatal error!");
                Console.WriteLine("  Overflow of bottom of rational fraction.");
                itop = 0;
                ibot = 1;
                return;
            }

            ibot = ibot1 * itop2;
            //
            //  The bottom should be positive.
            //
            if (ibot < 0)
            {
                ibot = -ibot;
                itop = -itop;
            }
            //
            //  The fraction is ITOP/IBOT with no loss of accuracy.
            //
        }

        public static void rat_farey(int n, int max_frac, ref int num_frac, ref int[] a, ref int[] b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    RAT_FAREY computes the N-th row of the Farey fraction table.
            //
            //  Example:
            //
            //    N = 5
            //
            //    NUM_FRAC = 11
            //    A =  0  1  1  1  2  1  3  2  3  4  1
            //    B =  1  5  4  3  5  2  5  3  4  5  1
            //
            //  Discussion:
            //
            //    In this form of the Farey fraction table, fractions in row N lie between
            //    0 and 1, are in lowest terms, and have a denominator that is no greater
            //    than N.  Row N is computed directly, and does not require the computation
            //    of previous rows.
            //
            //    The data satisfy the relationship:
            //
            //      A(K+1) * B(K) - A(K) * B(K+1) = 1
            //
            //    The number of items in the N-th row is roughly N^2 / PI^2.
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
            //    John Burkardt
            //
            //  Reference:
            //
            //    Donald Knuth,
            //    The Art of Computer Programming,
            //    Volume 1, Fundamental Algorithms,
            //    Addison Wesley, 1968, page 157.
            //
            //  Parameters:
            //
            //    Input, int N, the desired row number.  N must be positive.
            //
            //    Input, int MAX_FRAC, the maximum number of fractions to compute.
            //
            //    Output, ref int NUM_FRAC, the number of fractions computed.
            //
            //    Output, int A[MAX_FRAC], B[MAX_FRAC], contains the NUM_FRAC
            //    numerators and denominators of the N-th row of the Farey fraction table.
            //
        {
            int c;
            int k;

            if (n <= 0)
            {
                num_frac = 0;
                return;
            }

            k = 0;

            if (max_frac <= 0)
            {
                num_frac = k;
                return;
            }

            a[k] = 0;
            b[k] = 1;
            k = 1;

            if (max_frac <= 1)
            {
                num_frac = k;
                return;
            }

            a[k] = 1;
            b[k] = n;
            k = 2;

            while (k < max_frac)
            {
                if (a[k - 1] == 1 && b[k - 1] == 1)
                {
                    break;
                }

                c = (b[k - 2] + n) / b[k - 1];
                a[k] = c * a[k - 1] - a[k - 2];
                b[k] = c * b[k - 1] - b[k - 2];
                k = k + 1;
            }

            num_frac = k;

        }

        public static void rat_farey2(int n, ref int[] a, ref int[] b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    RAT_FAREY2 computes the next row of the Farey fraction table.
            //
            //  Example:
            //
            //    Input:
            //
            //      N = 3
            //      A =  0  1  1  2  1
            //      B =  1  3  2  3  1
            //
            //    Output:
            //
            //      A =  0  1  1  2  1  3  2  3  1
            //      B =  1  4  3  5  2  5  3  4  1
            //
            //  Discussion:
            //
            //    In this form of the Farey fraction table, fractions in row N lie between
            //    0 and 1, and are in lowest terms.  For every adjacent pair of input
            //    fractions, A1/B1 and A2/B2, the mediant (A1+A2)/(B1+B2) is computed
            //    and inserted between them.
            //
            //    The number of items in the N-th row is 1+2**(N-1).
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
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the input row number.  N must be nonnegative.
            //    If N is zero, then the input is ignored, and the entries of
            //    row 1 are computed directly.
            //
            //    Input/output, int A[1+2**N], B[1+2**N].
            //    On input, entries 1 through 1+2**(N-1) contain the entries of row N.
            //    On output, entries 1 through 1+2**N contain the entries of row N+1.
            //
        {
            int i;
            int ihi;

            if (n == 0)
            {
                a[0] = 0;
                b[0] = 1;
                a[1] = 1;
                b[1] = 1;
                return;
            }

            //
            //  Shift the current data.
            //
            ihi = (int)Math.Pow(2, n - 1);
            for (i = ihi; 0 <= i; i--)
            {
                a[2 * i] = a[i];
                b[2 * i] = b[i];
            }

            //
            //  Compute the mediants.
            //
            ihi = (int)Math.Pow(2, n) - 1;

            for (i = 1; i <= ihi; i = i + 2)
            {
                a[i] = a[i - 1] + a[i + 1];
                b[i] = b[i - 1] + b[i + 1];
            }

            return;
        }

        public static void rat_mul(int itop1, int ibot1, int itop2, int ibot2, ref int itop,
                ref int ibot, ref bool error)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    RAT_MUL multiplies two fractions.
            //
            //  Discussion:
            //
            //    The routine computes
            //
            //      ITOP / IBOT = ( ITOP1 / IBOT1 ) * ( ITOP2 / IBOT2 ).
            //
            //    while avoiding int overflow.
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
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int ITOP1, IBOT1, the first factor.
            //
            //    Input, int ITOP2, IBOT2, the second factor.
            //
            //    Output, ref int ITOP, &IBOT, the product.
            //
            //    Output, ref bool ERROR, is TRUE if an error occurred.
            //
        {
            int i_max;
            const int i4_huge = 2147483647;
            int itemp;

            error = false;

            i_max = i4_huge;

            if (itop1 == 0 || itop2 == 0)
            {
                itop = 0;
                ibot = 1;
                return;
            }

            //
            //  Get rid of all common factors in top and bottom.
            //
            itemp = typeMethods.i4_gcd(itop1, ibot1);
            itop1 = itop1 / itemp;
            ibot1 = ibot1 / itemp;
            itemp = typeMethods.i4_gcd(itop1, ibot2);
            itop1 = itop1 / itemp;
            ibot2 = ibot2 / itemp;
            itemp = typeMethods.i4_gcd(itop2, ibot1);
            itop2 = itop2 / itemp;
            ibot1 = ibot1 / itemp;
            itemp = typeMethods.i4_gcd(itop2, ibot2);
            itop2 = itop2 / itemp;
            ibot2 = ibot2 / itemp;
            //
            //  The fraction (ITOP1 * ITOP2) / (IBOT1*IBOT2) is in lowest terms.
            //
            //  Check the top ITOP1 * ITOP2 for overflow.
            //
            if ((double)(i_max) < Math.Abs((double)itop1 * (double)itop2))
            {
                error = true;
                Console.WriteLine("");
                Console.WriteLine("RAT_MUL - Fatal error!");
                Console.WriteLine("  Overflow of top of rational product.");
                itop = 0;
                return;
            }

            itop = itop1 * itop2;
            //
            //  Check the bottom IBOT1*IBOT2 for overflow.
            //
            if ((double)(i_max) < Math.Abs((double)ibot1 * (double)ibot2))
            {
                error = true;
                Console.WriteLine("");
                Console.WriteLine("RAT_MUL - Fatal error!");
                Console.WriteLine("  Overflow of bottom of rational product.");
                ibot = 1;
                return;
            }

            ibot = ibot1 * ibot2;
            //
            //  The bottom should be positive.
            //
            if (ibot < 0)
            {
                ibot = -ibot;
                itop = -itop;
            }

        }

        public static void rat_normalize(ref int a, ref int b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    RAT_NORMALIZE normalizes a rational number.
            //
            //  Discussion:
            //
            //    If A = B = 0, return.
            //
            //    If A = 0 (and B nonzero) set B => 1 and return.
            //
            //    If A nonzero, and B = 0, then A => 1 and return.
            //
            //    If A = B, then set A => 1, B => 1 and return.
            //
            //    If B < 0, then A => -A, B => -B.
            //
            //    If 1 < C = GCD(|A|,|B|), A => A/C, B => B/C.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 June 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input/output, ref int A, &B, the rational number.
            //
        {
            int c;
            //
            //  Cases where one or both is 0.
            //
            if (a == 0 && b == 0)
            {
                return;
            }
            else if (a == 0 && b != 0)
            {
                b = 1;
                return;
            }
            else if (a != 0 && b == 0)
            {
                a = 1;
                return;
            }

            if (a == b)
            {
                a = 1;
                b = 1;
                return;
            }

            if (b < 0)
            {
                a = -a;
                b = -b;
            }

            c = typeMethods.i4_gcd(Math.Abs(a), Math.Abs(b));

            if (1 < c)
            {
                a = a / c;
                b = b / c;
            }

            return;
        }

        public static void rat_sum_formula(int n, ref int[] a, ref int[] b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    RAT_SUM_FORMULA computes the formulas for sums of powers of integers.
            //
            //  Example:
            //
            //    N = 6
            //
            //        1    2    3    4    5    6    7
            //    -----------------------------------
            //    0 | 1    0    0    0    0    0    0
            //      |
            //    1 | 1    1    0    0    0    0    0
            //      | 2    2
            //      |
            //    2 | 1    1    1    0    0    0    0
            //      | 3    2    6
            //      |
            //    3 | 1    1    1    0    0    0    0
            //      | 4    2    4
            //      | 
            //    4 | 1    1    1    0   -1    0    0
            //      | 5    2    3        30
            //      |
            //    5 | 1    1    5    0   -1    0    0
            //      | 6    2   12        12
            //      |
            //    6 | 1    1    1    0   -1    0    1
            //      | 7    2    2         6        42
            //
            //    The interpretation of row 2, for instance, is:
            //
            //      sum ( 1 <= I <= N ) I^2 = 1/3 N^3 + 1/2 N^2 + 1/6 N
            //
            //    This suggests that a more sensible way to display the table
            //    is to reverse the order of the entries in the row, so that
            //    the entry in column J is the coeficient of N^J, which is
            //    not the case now.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    18 July 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Robert Owens,
            //    Sums of Powers of Integers,
            //    Mathematics Magazine,
            //    Volume 65, Number 1, February 1992, pages 38-40.
            //
            //  Parameters:
            //
            //    Input, int N, the number of rows of coefficients to compute.
            //
            //    Output, int A[(N+1)*(N+1)], B[(N+1)*(N+1)], the numerator and denominator
            //    of the coefficients.
            //
        {
            int asum = 0;
            int aval = 0;
            int bsum = 0;
            int bval = 0;
            bool error = false;
            int i;
            int j;

            a[0 + 0 * (n + 1)] = 1;
            for (j = 1; j < n + 1; j++)
            {
                a[0 + j * (n + 1)] = 0;
            }

            b[0 + 0 * (n + 1)] = 1;
            for (j = 1; j < n + 1; j++)
            {
                b[0 + j * (n + 1)] = 1;
            }

            for (i = 1; i <= n; i++)
            {
                asum = 0;
                bsum = 0;
                //
                //  Subdiagonal entries are multiples of entries above them.
                //
                for (j = 0; j < i; j++)
                {
                    aval = a[i - 1 + j * (n + 1)];
                    bval = b[i - 1 + j * (n + 1)];

                    rat_mul(aval, bval, i, i + 1 - j, ref aval, ref bval, ref error);

                    a[i + j * (n + 1)] = aval;
                    b[i + j * (n + 1)] = bval;

                    rat_add(asum, bsum, aval, bval, ref asum, ref bsum, ref error);
                }

                //
                //  Diagonal entry is 1 - sum of previous entries in row.
                //
                asum = -asum;

                rat_add(1, 1, asum, bsum, ref aval, ref bval, ref error);

                a[i + i * (n + 1)] = aval;
                b[i + i * (n + 1)] = bval;
                //
                //  Superdiagonal entries are zero.
                //
                for (j = i + 1; j < n + 1; j++)
                {
                    a[i + j * (n + 1)] = 0;
                }

                for (j = i + 1; j <= n; j++)
                {
                    b[i + j * (n + 1)] = 1;
                }
            }

            return;
        }

        public static void rat_to_cfrac(int ip, int iq, int m, ref int n, ref int[] a, ref bool error)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    RAT_TO_CFRAC converts a rational value to a continued fraction.
            //
            //  Discussion:
            //
            //    The routine is given a rational number represented by IP/IQ, and
            //    computes the monic or "simple" continued fraction representation
            //    with int coefficients of the number:
            //
            //      A(1) + 1/ (A(2) + 1/ (A(3) + ... + 1/A(N) ...))
            //
            //    The user must dimension A to a value M which is "large enough".
            //    The actual number of terms needed in the continued fraction
            //    representation cannot be known beforehand.
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
            //    FORTRAN77 version by John Hart, Ward Cheney, Charles Lawson, 
            //    Hans Maehly, Charles Mesztenyi, John Rice, Henry Thatcher, 
            //    Christoph Witzgall.
            //    C++ version by John Burkardt
            //
            //  Reference:
            //
            //    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi, 
            //    John Rice, Henry Thatcher, Christoph Witzgall,
            //    Computer Approximations,
            //    Wiley, 1968.
            //
            //  Parameters:
            //
            //    Input, int IP, IQ, the numerator and denominator of the
            //    rational value whose continued fraction representation is
            //    desired.
            //
            //    Input, int M, the dimension of A.  If M is not great
            //    enough, the algorithm may run out of space.
            //
            //    Output, ref int N, the actual number of entries used in A.
            //
            //    Output, int A[M], contains the continued fraction
            //    representation of the number.
            //
            //    Output, ref bool ERROR, is TRUE if an error occurred.
            //
        {
            error = false;

            n = 0;

            for (;;)
            {
                if (m <= n)
                {
                    error = true;
                    return;
                }

                a[n] = ip / iq;
                n = n + 1;

                ip = ip % iq;

                if (ip == 0)
                {
                    return;
                }

                if (m <= n)
                {
                    error = true;
                    return;
                }

                a[n] = iq / ip;
                n = n + 1;

                iq = iq % ip;

                if (iq == 0)
                {
                    break;
                }
            }

        }

        public static void rat_to_dec(int rat_top, int rat_bot, ref int mantissa, ref int exponent)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    RAT_TO_DEC converts a rational to a decimal representation.
            //
            //  Discussion:
            //
            //    A rational value is represented by RAT_TOP / RAT_BOT.
            //
            //    A decimal value is represented as MANTISSA * 10^EXPONENT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    22 November 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int RAT_TOP, RAT_BOT, the rational value.
            //
            //    Output, ref int MANTISSA, &EXPONENT, the decimal number.
            //
        {
            int gcd;
            const int i4_huge = 2147483647;
            double r;
            double r_max;
            int s;
            //
            //  Handle an input value of 0.
            //
            if (rat_top == 0)
            {
                mantissa = 0;
                exponent = 0;
                return;
            }

            //
            //  Take out any common factor.
            //
            gcd = typeMethods.i4_gcd(rat_top, rat_bot);
            rat_top = rat_top / gcd;
            rat_bot = rat_bot / gcd;
            //
            //  Force the bottom of the fraction to be positive.
            //
            if (rat_bot < 0)
            {
                rat_top = -rat_top;
                rat_bot = -rat_bot;
            }

            //
            //  If the fraction is a whole number, we can finish now.
            //
            if (rat_bot == 1)
            {
                mantissa = rat_top;
                exponent = 0;
                return;
            }

            //
            //  Whole factors of 10 in the bottom or top go into the decimal exponent.
            //
            exponent = 0;

            while ((rat_bot % 10) == 0)
            {
                exponent = exponent - 1;
                rat_bot = rat_bot / 10;
            }

            while ((rat_top % 10) == 0)
            {
                exponent = exponent + 1;
                rat_top = rat_top / 10;
            }

            //
            //  Convert to a real value.
            //
            r = (double)rat_top / (double)rat_bot;

            if (r < 0)
            {
                s = -1;
                r = -r;
            }
            else
            {
                s = 1;
            }

            r_max = ((double)i4_huge) / 10.0;

            while (r != (double)((int)r) && r < r_max)
            {
                r = r * 10.0;
                exponent = exponent - 1;
            }

            mantissa = s * (int)r;

            return;
        }

        public static double rat_to_r8(int top, int bot)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    RAT_TO_R8 converts rational values to R8's.
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
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int TOP, BOT, the rational quantity
            //    (TOP/BOT) that is to be converted.
            //
            //    Output, double RAT_TO_R8, the real value of the fraction.
            //
        {
            double value;

            value = (double)top / (double)bot;

            return value;
        }

        public static string rat_to_s(int a, int b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    RAT_TO_S converts a rational value to a string.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    01 June 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int A, B, the numerator and denominator.
            //
            //    Output, string S, a string representing the value.
            //
        {
            string s;
            string sa;
            string sb;

            if (b == 0)
            {
                if (a == 0)
                {
                    s = "NaN";
                }
                else if (a < 0)
                {
                    s = "-Inf";
                }
                else
                {
                    s = "Inf";
                }
            }
            else
            {
                sa = (a.ToString());
                sb = (b.ToString());
                s = sa + "/" + sb;
            }

            return s;
        }

        public static int rat_width(int a, int b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    RAT_WIDTH returns the "width" of a rational number.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 June 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int A, B, the rational number.
            //
            //    Output, int RAT_WIDTH, the "width" of the rational number.
            //
        {
            int abs_a;
            int abs_b;
            int ten_pow;
            int value;

            value = 1;
            ten_pow = 10;

            if (a == 0)
            {
                return value;
            }

            abs_a = Math.Abs(a);

            while (ten_pow <= abs_a)
            {
                value = value + 1;
                ten_pow = ten_pow * 10;
            }

            //
            //  If the fraction is negative, a minus sign will be prepended to the
            //  numerator.
            //
            if (a * b < 0)
            {
                value = value + 1;
                ten_pow = ten_pow * 10;
            }

            abs_b = Math.Abs(b);

            while (ten_pow <= abs_b)
            {
                value = value + 1;
                ten_pow = ten_pow * 10;
            }

            return value;
        }

        public static void ratmat_det(int n, int[] iatop, int[] iabot, ref int idtop, ref int idbot,
                ref bool error)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    RATMAT_DET finds the determinant of an N by N matrix of rational entries.
            //
            //  Discussion:
            //
            //    The brute force method is used.
            //
            //    This routine should only be used for small matrices, since this
            //    calculation requires the summation of N! products of N numbers.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 December 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of rows and columns of A.
            //
            //    Input, int IATOP[N*N], IABOT[N*N], the numerators
            //    and denominators of the entries of the matrix.
            //
            //    Output, ref int IDTOP, &IDBOT, the determinant of the matrix,
            //    expressed as IDTOP/IDBOT.
            //
            //    Output, ref bool ERROR, is TRUE if an error occurred.
            //
        {
            bool even = false;
            int i;
            int[] iarray;
            int ibot;
            int ibot2;
            int itop;
            int itop2;
            bool more;

            error = false;

            more = false;
            idtop = 0;
            idbot = 1;

            iarray = new int[n];

            for (;;)
            {
                Permutation.perm0_next(n, iarray, ref more, ref even);

                if (even)
                {
                    itop = 1;
                }
                else
                {
                    itop = -1;
                }

                ibot = 1;

                for (i = 0; i < n; i++)
                {
                    itop2 = iatop[i + iarray[i] * n];
                    ibot2 = iabot[i + iarray[i] * n];

                    rat_mul(itop, ibot, itop2, ibot2, ref itop, ref ibot, ref error);

                    if (error)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("RATMAT_DET - Fatal error!");
                        Console.WriteLine("  An overflow occurred.");
                        Console.WriteLine("  The determinant calculation cannot be done");
                        Console.WriteLine("  for this matrix.");
                        idtop = 0;
                        idbot = 1;
                        return;
                    }
                }

                rat_add(itop, ibot, idtop, idbot, ref idtop, ref idbot, ref error);

                if (error)
                {
                    Console.WriteLine("");
                    Console.WriteLine("RATMAT_DET - Fatal error!");
                    Console.WriteLine("  An overflow occurred.");
                    Console.WriteLine("  The determinant calculation cannot be done");
                    Console.WriteLine("  for this matrix.");
                    idtop = 0;
                    idbot = 1;
                    return;
                }

                if (!more)
                {
                    break;
                }
            }

            //
            //  The bottom should be positive.
            //
            if (idbot < 0)
            {
                idbot = -idbot;
                idtop = -idtop;
            }

        }

        public static void ratmat_print(int m, int n, int[] a, int[] b, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    RATMAT_PRINT prints out a rational matrix.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns in the matrix.
            //
            //    Input, int A[M*N], B[M*N], the rational matrix.
            //
            //    Input, string TITLE, a title.
            //
        {
            bool all_ones;
            int i;
            int j;
            int jmax;
            int jmin;
            int kmax;
            int ncolum = 80;
            int npline;
            string cout = "";
            //
            //  Figure out how many rationals we can get in NCOLUM columns.
            //
            //  Can I pass the VALUE of these entries this way?
            //
            kmax = 0;

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < n; j++)
                {
                    kmax = Math.Max(kmax, rat_width(a[i + j * m], b[i + j * m]));
                }
            }

            kmax = kmax + 2;
            npline = ncolum / kmax;
            //
            //  Now do the printing.
            //
            for (jmin = 0; jmin < n; jmin = jmin + npline)
            {
                jmax = Math.Min(jmin + npline - 1, n - 1);

                Console.WriteLine("");

                if (jmin == 0)
                {
                    if (0 < title.Length)
                    {
                        Console.WriteLine("");
                        Console.WriteLine(title + "");
                    }
                }

                if (0 < jmin || jmax < n - 1)
                {
                    Console.WriteLine("Columns " + jmin + " to " + jmax + "");
                    Console.WriteLine("");
                }

                for (i = 0; i < m; i++)
                {
                    for (j = jmin; j <= jmax; j++)
                    {
                        cout += a[i + j * m].ToString().PadLeft(kmax) + "  ";
                    }

                    Console.WriteLine(cout);
                    cout = "";
                    //
                    //  Delete each denominator that is 1.  If all are 1, don't
                    //  even print out the line.
                    //
                    all_ones = true;

                    for (j = jmin; j <= jmax; j++)
                    {
                        if (b[i + j * m] != 1)
                        {
                            all_ones = false;
                        }
                    }

                    if (!all_ones)
                    {
                        for (j = jmin; j <= jmax; j++)
                        {
                            cout += b[i + j * m].ToString().PadLeft(kmax) + "  ";
                        }

                        Console.WriteLine(cout);
                        cout = "";
                    }

                    if (jmax == n - 1 && i == m - 1)
                    {
                    }
                    else
                    {
                        Console.WriteLine(cout);
                        cout = "";
                    }
                }
            }

        }
    }
}