using System;
using System.Globalization;
using Burkardt.Function;

namespace Burkardt.Types;

public class i4
{
    public bool error { get; set; }
    public int val { get; set; }
    public int lchar { get; set; }
}

public class i4vec
{
    public bool error { get; set; }
    public int[] ivec { get; set; }
}

public static partial class typeMethods
{
    public static i4 s_to_i4(string s)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    S_TO_I4 reads an I4 from a string.
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
        i4 ival = new() {val = 0};

        int istate = 0;
        int isgn = 1;
        int i = 0;

        while (i < s.Length)
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
                    isgn = +1;
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
                ival.lchar = s_len_trim(s);
                break;
            default:
                ival.error = true;
                ival.lchar = 0;
                break;
        }

        return ival;
    }


    public static i4vec s_to_i4vec(string s, int n)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    S_TO_I4VEC reads an I4VEC from a string.
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
        //    Input, string S, the string to be read.
        //
        //    Input, int N, the number of values expected.
        //
        //    Output, int IVEC[N], the values read from the string.
        //
        //    Output, bool S_TO_I4VEC, is TRUE if an error occurred.
        //
    {
        i4vec ret = new() {ivec = new int[n]};

        string[] tokens = Helpers.splitStringByWhitespace(s);

        for (int i = 0; i < n; i++)
        {
            ret.ivec[i] = s_to_i4(tokens[i]).val;
        }

        return ret;
    }

    public static void i4_to_dvec ( int i4, int n, ref int[] dvec )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_TO_DVEC makes a signed decimal vector from an I4.
        //
        //  Discussion:
        //
        //    A DVEC is an integer vector of decimal digits, intended to
        //    represent an integer.  DVEC(1) is the units digit, DVEC(N-1)
        //    is the coefficient of 10^(N-2), and DVEC(N) contains sign
        //    information.  It is 0 if the number is positive, and 9 if
        //    the number is negative.
        //
        //    Negative values have a ten's complement operation applied.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I4, an integer to be represented.
        //
        //    Input, int N, the dimension of the vector.
        //
        //    Output, int DVEC[N], the signed decimal representation.
        //
    {
        const int base_ = 10;
        int i;

        int i4_copy = Math.Abs ( i4 );

        for ( i = 0; i < n-1; i++ )
        {
            dvec[i] = i4_copy % base_;

            i4_copy /= base_;
        }

        dvec[n-1] = 0;

        switch (i4)
        {
            case < 0:
                dvec_complementx ( n, dvec, ref dvec );
                break;
        }
    }

    public static void i4_to_i4poly ( int intval, int base_, int degree_max, ref int degree, ref int[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_TO_I4POLY converts an I4 to an integer polynomial in a given base.
        //
        //  Example:
        //
        //    INTVAL  BASE  Degree     A (in reverse order!)
        //
        //         1     2       0     1
        //         6     2       2     1  1  0
        //        23     2       4     1  0  1  1  1
        //        23     3       2     2  1  2
        //        23     4       2     1  1  3
        //        23     5       1     4  3
        //        23     6       1     3  5
        //        23    23       1     1  0
        //        23    24       0    23
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 July 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int INTVAL, an integer to be converted.
        //
        //    Input, int BASE, the base, which should be greater than 1.
        //
        //    Input, int DEGREE_MAX, the maximum degree.
        //
        //    Output, int &DEGREE, the degree of the polynomial.
        //
        //    Output, int A[DEGREE_MAX+1], contains the coefficients
        //    of the polynomial expansion of INTVAL in base BASE.
        //
    {
        int i;

        for ( i = 0; i <= degree_max; i++ )
        {
            a[i] = 0;
        }

        int j = Math.Abs ( intval );

        degree = 0;

        a[degree] = j % base_;

        j -= a[degree];
        j /= base_;

        while ( 0 < j )
        {
            degree += 1;

            if ( degree < degree_max )
            {
                a[degree] = j % base_;
            }

            j -= a[degree];
            j /= base_;
        }

        switch (intval)
        {
            case < 0:
            {
                for ( i = 0; i <= degree_max; i++ )
                {
                    a[i] = -a[i];
                }

                break;
            }
        }
    }

    public static void i4_to_triangle_lower(int k, ref int i, ref int j)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_TO_TRIANGLE_LOWER converts an integer to lower triangular coordinates.
        //
        //  Discussion:
        //
        //    Triangular coordinates are handy when storing a naturally triangular
        //    array (such as the lower half of a matrix) in a linear array.
        //
        //    Thus, for example, we might consider storing
        //
        //    (0,0)
        //    (1,0) (1,1)
        //    (2,0) (2,1) (2,2)
        //    (3,0) (3,1) (3,2) (3,3)
        //
        //    as the linear array
        //
        //    (0,0) (1,0) (1,1) (2,0) (2,1) (2,2) (3,0) (3,1) (3,2) (3,3)
        //
        //    Here, the quantities in parenthesis represent the natural row and
        //    column indices of a single number when stored in a rectangular array.
        //
        //    In this routine, we are given the location K of an item in the
        //    linear array, and wish to determine the row I and column J
        //    of the item when stored in the triangular array.
        //
        //  Example:
        //
        //    K  I  J
        //
        //    0  0  0
        //    1  1  0
        //    2  1  1
        //    3  2  0
        //    4  2  1
        //    5  2  2
        //    6  3  0
        //    7  3  1
        //    8  3  2
        //    9  3  3
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 January 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int K, the linear index of the (I,J) element, which
        //    must be nonnegative.
        //
        //    Output, int &I, &J, the row and column indices.
        //
    {
        switch (k)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("I4_TO_TRIANGLE_LOWER - Fatal error!");
                Console.WriteLine("  K < 0.");
                Console.WriteLine("  K = " + k + "");
                return;
            case 0:
                i = 0;
                j = 0;
                return;
        }

        //
        //   ( N - 1 )^2 + ( N - 1 ) < 2 * K <= N^2 + N
        //
        int r = (int)Math.Sqrt(2 * (k + 1));

        if (r * r + r < 2 * (k + 1))
        {
            r += 1;
        }

        r -= 1;

        int c = k - r * (r + 1) / 2;

        i = r;
        j = c;

    }

    public static void i4_to_triangle_upper(int k, ref int i, ref int j)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_TO_TRIANGLE_UPPER converts an integer to upper triangular coordinates.
        //
        //  Discussion:
        //
        //    Triangular coordinates are handy when storing a naturally triangular
        //    array (such as the upper half of a matrix) in a linear array.
        //
        //    Thus, for example, we might consider storing
        //
        //    (0,0) (0,1) (0,2) (0,3)
        //          (1,1) (1,2) (1,3)
        //                (2,2) (2,3)
        //                      (3,3)
        //
        //    as the linear array
        //
        //    (0,0) (0,1) (1,1) (0,2) (1,2) (2,2) (0,3) (1,3) (2,3) (3,3)
        //
        //    Here, the quantities in parenthesis represent the natural row and
        //    column indices of a single number when stored in a rectangular array.
        //
        //    In this routine, we are given the location K of an item in the
        //    linear array, and wish to determine the row I and column J
        //    of the item when stored in the triangular array.
        //
        //  Example:
        //
        //    K  I  J
        //
        //    0  0  0
        //    1  0  1
        //    2  1  1
        //    3  0  2
        //    4  1  2
        //    5  2  2
        //    6  0  3
        //    7  1  3
        //    8  2  3
        //    9  3  3
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 March 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int K, the linear index of the (I,J) element, which
        //    must be nonnegative.
        //
        //    Output, int &I, &J, the row and column indices.
        //
    {
        switch (k)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("I4_TO_TRIANGLE_UPPER - Fatal error!");
                Console.WriteLine("  K < 0.");
                Console.WriteLine("  K = " + k + "");
                return;
            case 0:
                i = 0;
                j = 0;
                return;
        }

        //
        //   ( N - 1 )^2 + ( N - 1 ) < 2 * K <= N^2 + N
        //
        int r = (int)Math.Sqrt(2 * (k + 1));

        if (r * r + r < 2 * (k + 1))
        {
            r += 1;
        }

        r -= 1;

        int c = k - r * (r + 1) / 2;

        j = r;
        i = c;

    }

    public static int i4_wrap(int ival, int ilo, int ihi)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_WRAP forces an I4 to lie between given limits by wrapping.
        //
        //  Example:
        //
        //    ILO = 4, IHI = 8
        //
        //    I   Value
        //
        //    -2     8
        //    -1     4
        //     0     5
        //     1     6
        //     2     7
        //     3     8
        //     4     4
        //     5     5
        //     6     6
        //     7     7
        //     8     8
        //     9     4
        //    10     5
        //    11     6
        //    12     7
        //    13     8
        //    14     4
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 August 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int IVAL, an integer value.
        //
        //    Input, int ILO, IHI, the desired bounds for the integer value.
        //
        //    Output, int I4_WRAP, a "wrapped" version of IVAL.
        //
    {
        int jlo = Math.Min(ilo, ihi);
        int jhi = Math.Max(ilo, ihi);

        int wide = jhi + 1 - jlo;

        int value = wide switch
        {
            1 => jlo,
            _ => jlo + i4_modp(ival - jlo, wide)
        };

        return value;
    }

    public static void i4block_print(int l, int m, int n, int[] a, string title)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4BLOCK_PRINT prints an I4BLOCK.
        //
        //  Discussion:
        //
        //    An I4BLOCK is a 3D array of I4 values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 June 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int L, M, N, the dimensions of the block.
        //
        //    Input, int A[L*M*N], the matrix to be printed.
        //
        //    Input, string TITLE, a title.
        //
    {
        Console.WriteLine();
        Console.WriteLine(title);

        for (int k = 0; k < n; k++)
        {
            Console.WriteLine();
            Console.WriteLine("  K = " + k);
            for (int jlo = 0; jlo < m; jlo += 10)
            {
                int jhi = Math.Min(jlo + 10, m);
                Console.WriteLine();
                string cout = "        J:";
                for (int j = jlo; j < jhi; j++)
                {
                    cout += "  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(6);
                }

                Console.WriteLine(cout);
                Console.WriteLine("       I:");
                for (int i = 0; i < l; i++)
                {
                    cout = "  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6) + ": ";
                    for (int j = jlo; j < jhi; j++)
                    {
                        cout += "  " + a[i + j * l + k * l * m].ToString(CultureInfo.InvariantCulture).PadLeft(6);
                    }

                    Console.WriteLine(cout);
                }
            }
        }
    }
        
    public static bool i4_is_fibonacci ( int i4 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_IS_FIBONACCI reports whether an integer is prime.
        //
        //  Discussion:
        //
        //    The positive integer i4 is a Fibonacci number if and only if
        //    5*I4^2+4 or 5*I4^2-4 is a perfect square.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 February 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I4, the integer to be tested.
        //
        //    Output, bool I4_IS_FIBONACCI, is TRUE if I4 is a Fibonacci number.
        //
    {
        switch (i4)
        {
            case <= 0:
                return false;
        }

        int t1 = 5 * i4 * i4 + 4;
        double t2 = Math.Sqrt ( t1 );
        int t3 = (int)Math.Floor ( t2 + 0.5 );
        if ( t3 * t3 == t1 )
        {
            return true;
        }

        t1 = 5 * i4 * i4 - 4;
        t2 = Math.Sqrt ( t1 );
        t3 = (int)Math.Floor ( t2 + 0.5 );
        return t3 * t3 == t1;
    }


    public static bool i4_is_power_of_10(int n)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_IS_POWER_OF_10 reports whether an I4 is a power of 10.
        //
        //  Discussion:
        //
        //    The powers of 10 are 1, 10, 100, 1000, 10000, and so on.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 March 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the integer to be tested.
        //
        //    Output, bool I4_IS_POWER_OF_10, is TRUE if N is a power of 10.
        //
    {
        switch (n)
        {
            case <= 0:
                return false;
        }

        while (1 < n)
        {
            if (n % 10 != 0)
            {
                return false;
            }

            n /= 10;
        }

        return true;
    }

    public static int i4_choose(int n, int k)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_CHOOSE computes the binomial coefficient C(N,K).
        //
        //  Discussion:
        //
        //    The value is calculated in such a way as to avoid overflow and
        //    roundoff.  The calculation is done in integer arithmetic.
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
        //    06 January 2013
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
        //    Input, int N, K, the values of N and K.
        //
        //    Output, int I4_CHOOSE, the number of combinations of N
        //    things taken K at a time.
        //
    {
        int value;

        int mn = k;
        int mx = n - k;
        if (mx < mn)
        {
            mn = n - k;
            mx = k;
        }

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
                value = mx + 1;

                for (int i = 2; i <= mn; i++)
                {
                    value = value * (mx + i) / i;
                }

                break;
            }
        }

        return value;
    }

    public static double i4_choose_log(int n, int k)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_CHOOSE_LOG computes the logarithm of the Binomial coefficient.
        //
        //  Discussion:
        //
        //    LOG ( C(N,K) ) = LOG ( N! / ( K! * (N-K)! ) ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 March 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, K, are the values of N and K.
        //
        //    Output, double I4_CHOOSE_LOG, the logarithm of C(N,K).
        //
    {
        double value = Helpers.LogGamma(n + 1)
                       - Helpers.LogGamma(k + 1)
                       - Helpers.LogGamma(n - k + 1);

        return value;
    }


    public static int i4_huge()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_HUGE returns a "huge" I4
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int I4_HUGE, a "huge" integer.
        //
    {
        const int value = 2147483647;

        return value;
    }


    public static int i4_sign(int i)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_SIGN returns the sign of an I4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 March 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, the integer whose sign is desired.
        //
        //    Output, int I4_SIGN, the sign of I.
    {
        return i switch
        {
            < 0 => -1,
            _ => 1
        };
    }

    public static void i4_sqrt ( int n, ref int q, ref int r )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_SQRT finds the integer square root of N by solving N = Q^2 + R.
        //
        //  Discussion:
        //
        //    The integer square root of N is an integer Q such that
        //    Q^2 <= N but N < (Q+1)^2.
        //
        //    A simpler calculation would be something like
        //
        //      Q = INT ( SQRT ( REAL ( N ) ) )
        //
        //    but this calculation has the virtue of using only integer arithmetic.
        //
        //    To avoid the tedium of worrying about negative arguments, the routine
        //    automatically considers the absolute value of the argument.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 June 2003
        //
        //  Author:
        //
        //   John Burkardt
        //
        //  Reference:
        //
        //    Mark Herkommer,
        //    Number Theory, A Programmer's Guide,
        //    McGraw Hill, 1999, pages 294-307.
        //
        //  Parameters:
        //
        //    Input, int N, the number whose integer square root is desired.
        //    Actually, only the absolute value of N is considered.
        //
        //    Output, int &Q, &R, the integer square root, and positive remainder,
        //    of N.
        //
    {
        n = Math.Abs(n);

        q = n;

        switch (n)
        {
            case > 0:
            {
                while (n / q < q)
                {
                    q = (q + n / q) / 2;
                }

                break;
            }
        }

        r = n - q * q;

    }

    public static void i4_sqrt_cf ( int n, int max_term, ref int n_term, ref int[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_SQRT_CF: continued fraction representation of a square root of an integer.
        //
        //  Discussion:
        //
        //    The continued fraction representation of the square root of an integer
        //    has the form
        //
        //      [ B0, (B1, B2, B3, ..., BM), ... ]
        //
        //    where
        //
        //      B0 = int ( sqrt ( real ( N ) ) )
        //      BM = 2 * B0
        //      the sequence ( B1, B2, B3, ..., BM ) repeats in the representation.
        //      the value M is termed the period of the representation.
        //
        //  Example:
        //
        //     N  Period  Continued Fraction
        //
        //     2       1  [ 1, 2, 2, 2, ... ]
        //     3       2  [ 1, 1, 2, 1, 2, 1, 2... ]
        //     4       0  [ 2 ]
        //     5       1  [ 2, 4, 4, 4, ... ]
        //     6       2  [ 2, 2, 4, 2, 4, 2, 4, ... ]
        //     7       4  [ 2, 1, 1, 1, 4, 1, 1, 4, 1, 1, 4... ]
        //     8       2  [ 2, 1, 4, 1, 4, 1, 4, 1, 4, ... ]
        //     9       0  [ 3 ]
        //    10       1  [ 3, 6, 6, 6, ... ]
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
        //   John Burkardt
        //
        //  Reference:
        //
        //    Mark Herkommer,
        //    Number Theory, A Programmer's Guide,
        //    McGraw Hill, 1999, pages 294-307.
        //
        //  Parameters:
        //
        //    Input, int N, the number whose continued fraction square root
        //    is desired.
        //
        //    Input, int MAX_TERM, the maximum number of terms that may
        //    be computed.
        //
        //    Output, int &N_TERM, the number of terms computed beyond the
        //    0 term.  The routine should stop if it detects that the period
        //    has been reached.
        //
        //    Output, int B[MAX_TERM+1], contains the continued fraction
        //    coefficients for indices 0 through N_TERM.
        //
    {
        int r = 0;
        int s = 0;

        n_term = 0;

        i4_sqrt(n, ref s, ref r);
        b[0] = s;

        switch (r)
        {
            case > 0:
            {
                int p = 0;
                int q = 1;

                for (;;)
                {
                    p = b[n_term] * q - p;
                    q = (n - p * p) / q;

                    if (max_term <= n_term)
                    {
                        return;
                    }

                    n_term += 1;
                    b[n_term] = (p + s) / q;

                    if (q == 1)
                    {
                        break;
                    }
                }

                break;
            }
        }
    }

    public static void i4_swap(ref int i, ref int j)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_SWAP switches two I4's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 January 2002
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, int *I, *J.  On output, the values of I and
        //    J have been interchanged.
        //
    {
        (i, j) = (j, i);
    }

    public static int i4_and(int i, int j)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_AND calculates the AND of two I4's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 January 2016
        //
        //  Author:
        //
        //   John Burkardt
        //
        //  Parameters:
        //
        //    Input, unsigned int I, J, two values whose AND is needed.
        //
        //    Output, unsigned int I4_AND, the AND of I and J.
        //
    {
        int k = 0;
        int l = 1;

        while (i != 0 || j != 0)
        {
            int i2 = i / 2;
            int j2 = j / 2;

            if (i != 2 * i2 && j != 2 * j2)
            {
                k += l;
            }

            i = i2;
            j = j2;
            l = 2 * l;
        }

        return k;
    }

    public static int i4_bclr(int i4, int pos)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_BCLR returns a copy of an I4 in which the POS-th bit is set to 0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Military Standard 1753,
        //    FORTRAN, DoD Supplement To American National Standard X3.9-1978,
        //    9 November 1978.
        //
        //  Parameters:
        //
        //    Input, int I4, the integer to be tested.
        //
        //    Input, int POS, the bit position, between 0 and 31.
        //
        //    Output, int I4_BCLR, a copy of I4, but with the POS-th bit
        //    set to 0.
        //
    {
        int value = i4;

        switch (pos)
        {
            case < 0:
                break;
            case < 31:
            {
                int sub = 1;

                int j = i4 switch
                {
                    >= 0 => i4,
                    _ => i4_huge() + i4 + 1
                };

                int k;
                for (k = 1; k <= pos; k++)
                {
                    j /= 2;
                    sub *= 2;
                }

                value = (j % 2) switch
                {
                    1 => i4 - sub,
                    _ => value
                };

                break;
            }
            case 31:
            {
                value = i4 switch
                {
                    < 0 => i4_huge() + i4 + 1,
                    _ => value
                };

                break;
            }
            case > 31:
                value = i4;
                break;
        }

        return value;
    }

    public static int i4_bit_hi1(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_BIT_HI1 returns the position of the high 1 bit base 2 in an I4.
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
        //    13 March 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the integer to be measured.
        //    N should be nonnegative.  If N is nonpositive, I4_BIT_HI1
        //    will always be 0.
        //
        //    Output, int I4_BIT_HI1, the location of the high order bit.
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

    public static int i4_bit_lo0(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_BIT_LO0 returns the position of the low 0 bit base 2 in an I4.
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
        //    08 February 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the integer to be measured.
        //    N should be nonnegative.
        //
        //    Output, int I4_BIT_LO0, the position of the low 1 bit.
        //
    {
        int bit = 0;

        while (true)
        {
            bit += 1;
            int n2 = n / 2;

            if (n == 2 * n2)
            {
                break;
            }

            n = n2;

        }

        return bit;
    }

    public static int i4_bit_lo1(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_BIT_LO1 returns the position of the low 1 bit base 2 in an I4.
        //
        //  Example:
        //
        //       N    Binary    Lo 1
        //    ----    --------  ----
        //       0           0     0
        //       1           1     1
        //       2          10     2
        //       3          11     1
        //       4         100     3
        //       5         101     1
        //       6         110     2
        //       7         111     1
        //       8        1000     4
        //       9        1001     1
        //      10        1010     2
        //      11        1011     1
        //      12        1100     3
        //      13        1101     1
        //      14        1110     2
        //      15        1111     1
        //      16       10000     5
        //      17       10001     1
        //    1023  1111111111     1
        //    1024 10000000000    11
        //    1025 10000000001     1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the integer to be measured.
        //    N should be nonnegative.
        //
        //    Output, int I4_BIT_LO1, the position of the low 1 bit.
        //
    {
        int bit = 0;
        int i = n;

        for (;;)
        {
            bit += 1;
            int i2 = i / 2;

            if (i != 2 * i2)
            {
                break;
            }

            i = i2;
        }

        return bit;
    }

    public static int i4_bit_reverse(int i, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_BIT_REVERSE reverses the bits in an I4.
        //
        //  Discussion:
        //
        //    An I4 is an int value.
        //
        //  Example:
        //
        //       I      N  2^N     I4_BIT_REVERSE ( I, N )
        //    ----    --------  -----------------------
        //       0      0    1     0
        //       1      0    1     1
        //
        //       0      3    8     0
        //       1      3    8     4
        //       2      3    8     2
        //       3      3    8     6
        //       4      3    8     1
        //       5      3    8     5
        //       6      3    8     3
        //       7      3    8     7
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, the integer to be bit reversed.
        //    I should be nonnegative.  Normally I < 2^N.
        //
        //    Input, int N, indicates the number of bits to
        //    be reverse (N+1) or the base with respect to which the integer is to
        //    be reversed (2^N).  N should be nonnegative.
        //
        //    Output, int I4_BIT_REVERSE, the bit reversed value.
        //
    {
        int value;

        switch (i)
        {
            case < 0:
                value = -1;
                break;
            default:
            {
                switch (n)
                {
                    case < 0:
                        value = -1;
                        break;
                    default:
                    {
                        int b = (int) Math.Pow(2, n);
                        int j = i % b;

                        value = 0;

                        for (;;)
                        {
                            if (b == 1)
                            {
                                value += j;
                                break;
                            }

                            switch (j % 2)
                            {
                                case 1:
                                    value += b / 2;
                                    j -= 1;
                                    break;
                            }

                            j /= 2;
                            b /= 2;
                        }

                        break;
                    }
                }

                break;
            }
        }

        return value;
    }

    public static int i4_bset(int i4, int pos)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_BSET returns a copy of an I4 in which the POS-th bit is set to 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Military Standard 1753,
        //    FORTRAN, DoD Supplement To American National Standard X3.9-1978,
        //    9 November 1978.
        //
        //  Parameters:
        //
        //    Input, int I4, the integer to be tested.
        //
        //    Input, int POS, the bit position, between 0 and 31.
        //
        //    Output, int I4_BSET, a copy of I4, but with the POS-th bit
        //    set to 1.
        //
    {
        int value = i4;

        switch (pos)
        {
            case < 0:
                break;
            case < 31:
            {
                int add = 1;

                int j = i4 switch
                {
                    >= 0 => i4,
                    _ => i4_huge() + i4 + 1
                };

                int k;
                for (k = 1; k <= pos; k++)
                {
                    j /= 2;
                    add *= 2;
                }

                value = (j % 2) switch
                {
                    0 => i4 + add,
                    _ => value
                };

                break;
            }
            case 31:
            {
                value = i4 switch
                {
                    > 0 => -(i4_huge() - i4) - 1,
                    _ => value
                };

                break;
            }
            case > 31:
                value = i4;
                break;
        }

        return value;
    }

    public static bool i4_btest(int i4, int pos)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_BTEST returns TRUE if the POS-th bit of an I4 is 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Military Standard 1753,
        //    FORTRAN, DoD Supplement To American National Standard X3.9-1978,
        //    9 November 1978.
        //
        //  Parameters:
        //
        //    Input, int I4, the integer to be tested.
        //
        //    Input, int POS, the bit position, between 0 and 31.
        //
        //    Output, bool I4_BTEST, is TRUE if the POS-th bit of I4 is 1.
        //
    {
        bool value;

        switch (pos)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("I4_BTEST - Fatal error!");
                Console.WriteLine("  POS < 0.");
                return false;
            case < 31:
            {
                int j = i4 switch
                {
                    >= 0 => i4,
                    _ => i4_huge() + i4 + 1
                };

                int k;
                for (k = 1; k <= pos; k++)
                {
                    j /= 2;
                }

                value = (j % 2) switch
                {
                    0 => false,
                    _ => true
                };

                break;
            }
            case 31 when i4 < 0:
                value = true;
                break;
            case 31:
                value = false;
                break;
            case > 31:
                Console.WriteLine("");
                Console.WriteLine("I4_BTEST - Fatal error!");
                Console.WriteLine("  31 < POS.");
                return false;
        }

        return value;
    }

    public static int i4_characteristic(int q)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_CHARACTERISTIC gives the characteristic for an I4.
        //
        //  Discussion:
        //
        //    For any positive integer Q, the characteristic is:
        //
        //    Q, if Q is a prime;
        //    P, if Q = P**N for some prime P and some integer N;
        //    0, otherwise, that is, if Q is negative, 0, 1, or the product
        //       of more than one distinct prime.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 March 2004
        //
        //  Reference:
        //
        //    Paul Bratley, Bennett Fox, Harald Niederreiter,
        //    Algorithm 738:
        //    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
        //    ACM Transactions on Mathematical Software,
        //    Volume 20, Number 4, pages 494-495, 1994.
        //
        //  Parameters:
        //
        //    Input, int Q, the value to be tested.
        //
        //    Output, int I4_CHARACTERISTIC, the characteristic of Q.
        //
    {
        int i;

        switch (q)
        {
            case <= 1:
                return 0;
        }

        //
        //  If Q is not prime, then there is at least one prime factor
        //  of Q no greater than SQRT(Q)+1.
        //
        //  A faster code would only consider prime values of I,
        //  but that entails storing a table of primes and limiting the
        //  size of Q.  Simplicity and flexibility for now//
        //
        int i_maximum = (int) Math.Sqrt(q) + 1;

        for (i = 2; i <= i_maximum; i++)
        {
            switch (q % i)
            {
                case 0:
                {
                    while (q % i == 0)
                    {
                        q /= i;
                    }

                    return q switch
                    {
                        1 => i,
                        _ => 0
                    };
                }
            }
        }

        //
        //  If no factor was found, then Q is prime.
        //
        return q;
    }

    public static bool i4_choose_check(int n, int k)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_CHOOSE_CHECK reports whether the binomial coefficient can be computed.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 March 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, K, the binomial parameters.
        //
        //    Output, bool I4_CHOOSE_CHECK is:
        //    TRUE, if C(N,K) < maximum integer.
        //    FALSE, otherwise.
        //
    {
        double i4_huge_log = Math.Log(i4_huge());

        double choose_nk_log = r8_gamma_log(n + 1)
                               - r8_gamma_log(k + 1)
                               - r8_gamma_log(n - k + 1);

        bool check = choose_nk_log < i4_huge_log;

        return check;
    }

    public static int i4_divp(int i, int j)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_DIVP returns the smallest multiple of J greater than or equal to an I4.
        //
        //  Example:
        //
        //    I  J  I4_DIVP(I,J)
        //
        //    0  4    0
        //    1  4    1
        //    2  4    1
        //    3  4    1
        //    4  4    1
        //    5  4    2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 March 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, the number to be analyzed.
        //
        //    Input, int J, the number, multiples of which will
        //    be compared against I.  J may not be zero.
        //
        //    Output, int I4_DIVP, the smallest multiple of J that
        //    is greater than or equal to I.
        //
    {
        switch (j)
        {
            case 0:
                Console.WriteLine("");
                Console.WriteLine("I4_DIVP - Fatal error!");
                Console.WriteLine("  The input value of J was zero.");
                return 0;
            default:
            {
                int value = 1 + (i - 1) / j;

                return value;
            }
        }
    }

    public static void i4_factor ( int n, int maxfactor, ref int nfactor, ref int[] factor, 
            ref int[] exponent, ref int nleft )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_FACTOR factors an I4 into prime factors.
        //
        //  Discussion:
        //
        //    N = NLEFT * Product ( 1 <= I <= NFACTOR ) FACTOR(I)^EXPONENT(I).
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
        //    Input, int N, the integer to be factored.  N may be positive,
        //    negative, or 0.
        //
        //    Input, int MAXFACTOR, the maximum number of prime factors for
        //    which storage has been allocated.
        //
        //    Output, int &NFACTOR, the number of prime factors of N discovered
        //    by the routine.
        //
        //    Output, int FACTOR[MAXFACTOR], the prime factors of N.
        //
        //    Output, int EXPONENT[MAXFACTOR].  EXPONENT(I) is the power of
        //    the FACTOR(I) in the representation of N.
        //
        //    Output, int &NLEFT, the factor of N that the routine could not
        //    divide out.  If NLEFT is 1, then N has been completely factored.
        //    Otherwise, NLEFT represents factors of N involving large primes.
        //
    {
        int i;

        nfactor = 0;

        for (i = 0; i < maxfactor; i++)
        {
            factor[i] = 0;
        }

        for (i = 0; i < maxfactor; i++)
        {
            exponent[i] = 0;
        }

        nleft = n;

        switch (n)
        {
            case 0:
                return;
        }

        switch (Math.Abs(n))
        {
            case 1:
                nfactor = 1;
                factor[0] = 1;
                exponent[0] = 1;
                return;
        }

        //
        //  Find out how many primes we stored.
        //
        int maxprime = Prime.prime(-1);
        //
        //  Try dividing the remainder by each prime.
        //
        for (i = 1; i <= maxprime; i++)
        {
            int p = Prime.prime(i);

            if (Math.Abs(nleft) % p != 0)
            {
                continue;
            }

            if (nfactor >= maxfactor)
            {
                continue;
            }

            nfactor += 1;
            factor[nfactor - 1] = p;
            exponent[nfactor - 1] = 0;

            for (;;)
            {
                exponent[nfactor - 1] += 1;
                nleft /= p;

                if (Math.Abs(nleft) % p != 0)
                {
                    break;
                }
            }

            if (Math.Abs(nleft) == 1)
            {
                break;
            }
        }
    }

    public static int i4_factorial(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_FACTORIAL computes the factorial of N.
        //
        //  Discussion:
        //
        //    factorial ( N ) = product ( 1 <= I <= N ) I
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 June 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the argument of the factorial function.
        //    If N is less than 1, the function value is returned as 1.
        //    0 <= N <= 13 is required.
        //
        //    Output, int I4_FACTORIAL, the factorial of N.
        //
    {
        int value = 1;

        switch (n)
        {
            case > 13:
                Console.WriteLine("I4_FACTORIAL - Fatal error!");
                Console.WriteLine("  I4_FACTORIAL(N) cannot be computed as an integer");
                Console.WriteLine("  for 13 < N.");
                Console.WriteLine("  Input value N = " + n + "");
                return value;
        }

        for (int i = 1; i <= n; i++)
        {
            value *= i;
        }

        return value;
    }

    public static double i4_factorial_log(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_FACTORIAL_LOG returns the logarithm of N factorial.
        //
        //  Discussion:
        //
        //    N! = Product ( 1 <= I <= N ) I
        //
        //    N! = Gamma(N+1).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 March 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the argument of the function.
        //    0 <= N.
        //
        //    Output, double I4_FACTORIAL_LOG, the logarithm of N factorial.
        //
    {
        switch (n)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("I4_FACTORIAL_LOG - Fatal error!");
                Console.WriteLine("  N < 0.");
                return 0;
        }

        double value = 0.0;

        for (int i = 2; i <= n; i++)
        {
            value += Math.Log(i);
        }

        return value;
    }

    public static void i4_factorial_values(ref int n_data, ref int n, ref int fn)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_FACTORIAL_VALUES returns values of the factorial function.
        //
        //  Discussion:
        //
        //    0! = 1
        //    I! = Product ( 1 <= J <= I ) I
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      n!
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 August 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //    Stephen Wolfram,
        //    The Mathematica Book,
        //    Fourth Edition,
        //    Cambridge University Press, 1999,
        //    ISBN: 0-521-64314-7,
        //    LC: QA76.95.W65.
        //
        //  Parameters:
        //
        //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, int &N, the argument of the function.
        //
        //    Output, int &FN, the value of the function.
        //
    {
        const int N_MAX = 13;

        int[] fn_vec =
            {
                1,
                1,
                2,
                6,
                24,
                120,
                720,
                5040,
                40320,
                362880,
                3628800,
                39916800,
                479001600
            }
            ;

        int[] n_vec =
            {
                0, 1, 2, 3,
                4, 5, 6, 7,
                8, 9, 10, 11,
                12
            }
            ;

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

        if (N_MAX < n_data)
        {
            n_data = 0;
            n = 0;
            fn = 0;
        }
        else
        {
            n = n_vec[n_data - 1];
            fn = fn_vec[n_data - 1];
        }
    }

    public static int i4_factorial2(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_FACTORIAL2 computes the double factorial function.
        //
        //  Discussion:
        //
        //    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
        //                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the argument of the double factorial function.
        //    If N is less than 1, I4_FACTORIAL2 is returned as 1.
        //
        //    Output, int I4_FACTORIAL2, the value of the double factorial function.
        //
    {
        int value;

        switch (n)
        {
            case < 1:
                value = 1;
                return value;
        }

        int n_copy = n;
        value = 1;

        while (1 < n_copy)
        {
            value *= n_copy;
            n_copy -= 2;
        }

        return value;
    }

    public static void i4_factorial2_values(ref int n_data, ref int n, ref int fn)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_FACTORIAL2_VALUES returns values of the double factorial function.
        //
        //  Formula:
        //
        //    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
        //                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      n!!
        //
        //  Example:
        //
        //     N    N!!
        //
        //     0     1
        //     1     1
        //     2     2
        //     3     3
        //     4     8
        //     5    15
        //     6    48
        //     7   105
        //     8   384
        //     9   945
        //    10  3840
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 August 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //    Stephen Wolfram,
        //    The Mathematica Book,
        //    Fourth Edition,
        //    Cambridge University Press, 1999,
        //    ISBN: 0-521-64314-7,
        //    LC: QA76.95.W65.
        //
        //    Daniel Zwillinger,
        //    CRC Standard Mathematical Tables and Formulae,
        //    30th Edition,
        //    CRC Press, 1996, page 16.
        //
        //  Parameters:
        //
        //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, int &N, the argument of the function.
        //
        //    Output, int &FN, the value of the function.
        //
    {
        const int N_MAX = 16;

        int[] fn_vec =
            {
                1,
                1,
                2,
                3,
                8,
                15,
                48,
                105,
                384,
                945,
                3840,
                10395,
                46080,
                135135,
                645120,
                2027025
            }
            ;

        int[] n_vec =
            {
                0,
                1, 2, 3, 4, 5,
                6, 7, 8, 9, 10,
                11, 12, 13, 14, 15
            }
            ;

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

        if (N_MAX < n_data)
        {
            n_data = 0;
            n = 0;
            fn = 0;
        }
        else
        {
            n = n_vec[n_data - 1];
            fn = fn_vec[n_data - 1];
        }
    }

    public static int i4_fall(int x, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_FALL computes the falling factorial function [X]_N.
        //
        //  Discussion:
        //
        //    Note that the number of "injections" or 1-to-1 mappings from
        //    a set of N elements to a set of M elements is [M]_N.
        //
        //    The number of permutations of N objects out of M is [M]_N.
        //
        //    Moreover, the Stirling numbers of the first kind can be used
        //    to convert a falling factorial into a polynomial, as follows:
        //
        //      [X]_N = S^0_N + S^1_N * X + S^2_N * X^2 + ... + S^N_N X^N.
        //
        //    The formula is:
        //
        //      [X]_N = X * ( X - 1 ) * ( X - 2 ) * ... * ( X - N + 1 ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 November 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X, the argument of the falling factorial function.
        //
        //    Input, int N, the order of the falling factorial function.
        //    If N = 0, FALL = 1, if N = 1, FALL = X.  Note that if N is
        //    negative, a "rising" factorial will be computed.
        //
        //    Output, int I4_FALL, the value of the falling factorial function.
        //
    {
        int i;

        int value = 1;

        switch (n)
        {
            case > 0:
            {
                for (i = 1; i <= n; i++)
                {
                    value *= x;
                    x -= 1;
                }

                break;
            }
            case < 0:
            {
                for (i = -1; n <= i; i--)
                {
                    value *= x;
                    x += 1;
                }

                break;
            }
        }

        return value;
    }

    public static void i4_fall_values(ref int n_data, ref int m, ref int n, ref int fmn)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_FALL_VALUES returns values of the integer falling factorial function.
        //
        //  Discussion:
        //
        //    The definition of the falling factorial function is
        //
        //      (m)_n = (m)! / (m-n)!
        //            = ( m ) * ( m - 1 ) * ( m - 2 ) ... * ( m - n + 1 )
        //            = Gamma ( m + 1 ) / Gamma ( m - n + 1 )
        //
        //    We assume 0 <= N <= M.
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      FactorialPower[m,n]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 December 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //    Stephen Wolfram,
        //    The Mathematica Book,
        //    Fourth Edition,
        //    Cambridge University Press, 1999,
        //    ISBN: 0-521-64314-7,
        //    LC: QA76.95.W65.
        //
        //  Parameters:
        //
        //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, int &M, &N, the arguments of the function.
        //
        //    Output, int &FMN, the value of the function.
        //
    {
        const int N_MAX = 15;

        int[] fmn_vec =
            {
                1, 5, 20, 60, 120,
                120, 0, 1, 10, 4000,
                90, 4896, 24, 912576, 0
            }
            ;

        int[] m_vec =
            {
                5, 5, 5, 5, 5,
                5, 5, 50, 10, 4000,
                10, 18, 4, 98, 1
            }
            ;

        int[] n_vec =
            {
                0, 1, 2, 3, 4,
                5, 6, 0, 1, 1,
                2, 3, 4, 3, 7
            }
            ;

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

        if (N_MAX < n_data)
        {
            n_data = 0;
            m = 0;
            n = 0;
            fmn = 0;
        }
        else
        {
            m = m_vec[n_data - 1];
            n = n_vec[n_data - 1];
            fmn = fmn_vec[n_data - 1];
        }
    }

    public static int i4_fraction(int i, int j)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_FRACTION computes a ratio and returns an integer result.
        //
        //  Discussion:
        //
        //    Given integer variables I and J, FORTRAN will evaluate the expression
        //    "I/J" using integer arithmetic.  This routine, which carries out the
        //    same operation, is thus not needed in FORTRAN.  It is provided simply
        //    to match the corresponding function in MATLAB, where the default
        //    result of "I/J" is a real number.
        //
        //  Example:
        //
        //       I     J     Real     K = I4_FRACTION ( I, J)
        //
        //       1     2     0.5      0
        //       8     4     2.00     2
        //       9     4     2.25     2
        //       7     4     1.75     1
        //      -7     4    -1.75    -1
        //       7    -4    -1.75    -1
        //      -7    -4     1.75     1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, J, the arguments.
        //
        //    Output, int K, the value of the ratio.
    {
        int k = i / j;

        return k;
    }

    public static int i4_gcd(int i, int j)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_GCD finds the greatest common divisor of two I4's.
        //
        //  Discussion:
        //
        //    Note that only the absolute values of I and J are
        //    considered, so that the result is always nonnegative.
        //
        //    If I or J is 0, I4_GCD is returned as max ( 1, abs ( I ), abs ( J ) ).
        //
        //    If I and J have no common factor, I4_GCD is returned as 1.
        //
        //    Otherwise, using the Euclidean algorithm, I4_GCD is the
        //    greatest common divisor of I and J.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 March 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, J, two numbers whose GCD is desired.
        //
        //    Output, int I4_GCD, the greatest common divisor of I and J.
        //
    {
        int q;
        switch (i)
        {
            //
            //  Return immediately if either I or J is zero.
            //
            case 0:
                q = Math.Max(1, Math.Abs(j));
                return q;
        }

        switch (j)
        {
            case 0:
                q = Math.Max(1, Math.Abs(i));
                return q;
        }

        //
        //  Set IP to the larger of I and J, IQ to the smaller.
        //  This way, we can alter IP and IQ as we go.
        //
        int p = Math.Max(Math.Abs(i), Math.Abs(j));
        q = Math.Min(Math.Abs(i), Math.Abs(j));
        //
        //  Carry out the Euclidean algorithm.
        //
        for (;;)
        {
            int r = p % q;

            if (r == 0)
            {
                break;
            }

            p = q;
            q = r;
        }

        return q;
    }

    public static int i4_gcdb(int i, int j, int k)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_GCDB finds the greatest common divisor of the form K**N of two I4's.
        //
        //  Discussion:
        //
        //    Note that if J is negative, I4_GCDB will also be negative.
        //    This is because it is likely that the caller is forming
        //    the fraction I/J, and so any minus sign should be
        //    factored out of J.
        //
        //    If I and J are both zero, I4_GCDB is returned as 1.
        //
        //    If I is zero and J is not, I4_GCDB is returned as J,
        //    and vice versa.
        //
        //    If I and J are nonzero, and have no common divisor of the
        //    form K**N, I4_GCDB is returned as 1.
        //
        //    Otherwise, I4_GCDB is returned as the largest common divisor
        //    of the form K**N shared by I and J.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 March 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, J, two numbers whose greatest common divisor K**N
        //    is desired.
        //
        //    Input, int K, the possible divisor of I and J.
        //
        //    Output, int I4_GCDB, the greatest common divisor of
        //    the form K^N shared by I and J.
        //
    {
        int value;
        switch (i)
        {
            //
            //  If both I and J are zero, I4_GCDB is 1.
            //
            case 0 when j == 0:
                value = 1;
                return value;
            //
            //  If just one of I and J is zero, I4_GCDB is the other one.
            //
            case 0:
                value = j;
                return value;
            default:
            {
                switch (j)
                {
                    case 0:
                        value = i;
                        return value;
                }

                break;
            }
        }

        value = j switch
        {
            //
            //  Divide out K as long as you can.
            //
            > 0 => 1,
            _ => -1
        };

        for (;;)
        {
            if (i % k != 0 || j % k != 0)
            {
                break;
            }

            value *= k;
            i /= k;
            j /= k;

        }

        return value;
    }

    public static double i4_huge_normalizer()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_HUGE_NORMALIZER returns the "normalizer" for I4_HUGE.
        //
        //  Discussion:
        //
        //    The value returned is 1 / ( I4_HUGE + 1 ).
        //
        //    For any I4, it should be the case that
        //
        //     -1 < I4 * I4_HUGE_NORMALIZER < 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double I4_HUGE_NORMALIZER, the "normalizer"
        //    for I4_HUGE.
        //
    {
        const double value = 4.656612873077392578125E-10;

        return value;
    }

    public static bool i4_is_power_of_2(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_IS_POWER_OF_2 reports whether an I4 is a power of 2.
        //
        //  Discussion:
        //
        //    The powers of 2 are 1, 2, 4, 8, 16, and so on.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 March 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the integer to be tested.
        //
        //    Output, bool I4_IS_POWER_OF_2, is TRUE if N is a power of 2.
        //
    {
        bool value;

        switch (n)
        {
            case <= 0:
                value = false;
                return value;
        }

        while (1 < n)
        {
            switch (n % 2)
            {
                case 1:
                    value = false;
                    return value;
                default:
                    n /= 2;
                    break;
            }
        }

        value = true;

        return value;
    }

    public static bool i4_is_prime(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_IS_PRIME reports whether an I4 is prime.
        //
        //  Discussion:
        //
        //    A simple, unoptimized sieve of Erasthosthenes is used to
        //    check whether N can be divided by any integer between 2
        //    and SQRT(N).
        //
        //    Note that negative numbers, 0 and 1 are not considered prime.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 March 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the integer to be tested.
        //
        //    Output, bool I4_IS_PRIME, is TRUE if N is prime, and FALSE
        //    otherwise.
        //
    {
        int i;

        switch (n)
        {
            case <= 0:
            case 1:
                return false;
            case <= 3:
                return true;
        }

        int nhi = (int) Math.Sqrt(n);

        for (i = 2; i <= nhi; i++)
        {
            switch (n % i)
            {
                case 0:
                    return false;
            }
        }

        return true;
    }
        
    public static bool i4_is_triangular ( int i )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_IS_TRIANGULAR determines whether an I4 is triangular.
        //
        //  Discussion:
        //
        //    The N-th triangular number is equal to the sum of the first
        //    N integers.
        //
        //  First Values:
        //
        //    Index  Value
        //     0      0
        //     1      1
        //     2      3
        //     3      6
        //     4     10
        //     5     15
        //     6     21
        //     7     28
        //     8     36
        //     9     45
        //    10     55
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 February 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, the integer to be checked.
        //
        //    Output, bool I4_IS_TRIANGULAR, is TRUE if I is triangular.
        //
    {
        int j = 0;
        int k = 0;

        switch (i)
        {
            case < 0:
                return false;
            case 0:
                return true;
            default:
            {
                i4_to_triangle_lower ( i, ref j, ref k );

                return j == k;
            }
        }

    }

    public static int i4_lcm(int i, int j)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_LCM computes the least common multiple of two I4's.
        //
        //  Discussion:
        //
        //    The least common multiple may be defined as
        //
        //      LCM(I,J) = ABS( I * J ) / GCF(I,J)
        //
        //    where GCF(I,J) is the greatest common factor of I and J.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 May 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, J, the integers whose LCM is desired.
        //
        //    Output, int I4_LCM, the least common multiple of I and J.
        //    I4_LCM is never negative.  I4_LCM is 0 if either I or J is zero.
        //
    {
        int value = Math.Abs(i * (j / i4_gcd(i, j)));

        return value;
    }

    public static int i4_lcm_12n(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_LCM_12N computes the least common multiple of the integers 1 through N.
        //
        //  Examples:
        //
        //    N    LCM_12N
        //
        //    1          1
        //    2          2
        //    3          3
        //    4         12
        //    5         60
        //    6         60
        //    7        420
        //    8        840
        //    9       2520
        //   10       2520
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 August 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the value of N.
        //
        //    Output, int I4_LCM_12N, the least common multiple of
        //    the integers 1 to N.
        //
    {
        int i;

        int value = 1;

        for (i = 2; i <= n; i++)
        {
            int mult = i;
            int j;
            for (j = 1; j < i; j++)
            {
                switch (mult % (i - j))
                {
                    case 0:
                        mult /= i - j;
                        break;
                }
            }

            value *= mult;
        }

        return value;
    }

    public static void i4_mant(double x, ref int s, ref int j, ref int k, ref int l)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_MANT computes the "mantissa" of an R4.
        //
        //  Discussion:
        //
        //    This function computes the "mantissa" or "fraction part" of a real
        //    number X, which it stores as a pair of integers, (J/K).
        //
        //    It also computes the sign, and the integer part of the logarithm
        //    (base 2) of X.
        //
        //    On return:
        //
        //      X = S * (J/K) * 2^L
        //
        //    where
        //
        //      S is +1 or -1,
        //      K is a power of 2,
        //      1 <= (J/K) < 2,
        //      L is an integer.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the real number to be decomposed.
        //
        //    Output, int &S, the "sign" of the number.
        //    S will be -1 if X is less than 0, and +1 if X is greater
        //    than or equal to zero.
        //
        //    Output, int &J, the top part of the mantissa fraction.
        //
        //    Output, int &K, the bottom part of the mantissa
        //    fraction.  K is a power of 2.
        //
        //    Output, int &L, the integer part of the logarithm (base 2) of X.
        //
    {
        switch (x)
        {
            //
            //  1: Handle the special case of 0.
            //
            case 0.0:
                s = 1;
                j = 0;
                k = 1;
                l = 0;
                return;
            //
            //  2: Determine the sign IS.
            //
            case > 0.0:
                s = 1;
                break;
            default:
                s = -1;
                x = -x;
                break;
        }

        //
        //  3: Force X to lie between 1 and 2, and compute the logarithm L.
        //
        l = 0;

        while (2.0 <= x)
        {
            x /= 2.0;
            l += 1;
        }

        while (x < 1.0)
        {
            x *= 2.0;
            l -= 1;
        }

        //
        //  4: Now strip out the mantissa as J/K.
        //
        j = 0;
        k = 1;

        for (;;)
        {
            j = 2 * j;

            switch (x)
            {
                case >= 1.0:
                    j += 1;
                    x -= 1.0;
                    break;
            }

            if (x == 0.0)
            {
                break;
            }

            k = 2 * k;
            x *= 2.0;

        }
    }

    public static int i4_mod_inv(int b, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_MOD_INV calculates the inverse of B mod N.
        //
        //  Discussion:
        //
        //    This function uses the extended Euclidean algorithm.
        //
        //    Unless the algorithm fails, the output value Y will satisfy
        //
        //      ( B * Y ) mod N = 1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 August 2011
        //
        //  Author:
        //
        //    Original MATLAB version by Wade Trappe, Lawrence Washington.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Wade Trappe, Lawrence Washington,
        //    Introduction to Cryptography with Coding Theory,
        //    Prentice Hall, 2005,
        //    ISBN13: 978-0131862395,
        //    LC: QA268.T73.
        //
        //  Parameters:
        //
        //    Input, int B, the value whose inverse is desired.
        //    B must not be 0, or a multiple of N.  However, B can be negative.
        //
        //    Input, int N, the value with respect to which the inverse
        //    is desired.  N must be 2 or greater.
        //
        //    Output, int I4_MOD_INV, the inverse of B mod N.  However, if the
        //    inverse does not exist, Y is returned as 0.
        //
    {
        int y = 0;

        int n0 = n;
        int b0 = Math.Abs(b);
        int t0 = 0;
        int t = 1;

        int q = n0 / b0;
        int r = n0 - q * b0;

        while (0 < r)
        {
            int temp = t0 - q * t;

            switch (temp)
            {
                case >= 0:
                    temp %= n;
                    break;
                default:
                    temp = n - -temp % n;
                    break;
            }

            n0 = b0;
            b0 = r;
            t0 = t;
            t = temp;

            q = n0 / b0;
            r = n0 - q * b0;
        }

        if (b0 != 1)
        {
            y = 0;
        }
        else
        {
            y = b switch
            {
                < 0 => -y,
                _ => t % n
            };
        }

        return y;
    }

    public static void i4_moddiv(int n, int d, ref int m, ref int r)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_MODDIV breaks an I4 into a multiple of a divisor and remainder.
        //
        //  Discussion:
        //
        //    N = M * D + R
        //
        //    0 <= || R || < || D ||
        //
        //    R has the sign of N.
        //
        //  Example:
        //
        //    N         D       M      R
        //
        //   107       50      2      7
        //   107      -50     -2      7
        //  -107       50     -2     -7
        //  -107      -50      2     -7
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number to be decomposed.
        //
        //    Input, int D, the divisor.  D may not be zero.
        //
        //    Output, int &M, the number of times N is evenly divided by D.
        //
        //    Output, int &R, a remainder.
        //
    {
        switch (d)
        {
            case 0:
                Console.WriteLine("");
                Console.WriteLine("I4_MODDIV - Fatal error!");
                Console.WriteLine("  Input divisor D = 0");
                return;
            default:
                m = n / d;
                r = n - d * m;
                break;
        }
    }

    public static int i4_modp(int i, int j)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_MODP returns the nonnegative remainder of I4 division.
        //
        //  Discussion:
        //
        //    If
        //      NREM = I4_MODP ( I, J )
        //      NMULT = ( I - NREM ) / J
        //    then
        //      I = J * NMULT + NREM
        //    where NREM is always nonnegative.
        //
        //    The MOD function computes a result with the same sign as the
        //    quantity being divided.  Thus, suppose you had an angle A,
        //    and you wanted to ensure that it was between 0 and 360.
        //    Then mod(A,360) would do, if A was positive, but if A
        //    was negative, your result would be between -360 and 0.
        //
        //    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
        //
        //        I         J     MOD  I4_MODP   I4_MODP Factorization
        //
        //      107        50       7       7    107 =  2 *  50 + 7
        //      107       -50       7       7    107 = -2 * -50 + 7
        //     -107        50      -7      43   -107 = -3 *  50 + 43
        //     -107       -50      -7      43   -107 =  3 * -50 + 43
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 May 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, the number to be divided.
        //
        //    Input, int J, the number that divides I.
        //
        //    Output, int I4_MODP, the nonnegative remainder when I is
        //    divided by J.
        //
    {
        switch (j)
        {
            case 0:
                Console.WriteLine("");
                Console.WriteLine("I4_MODP - Fatal error!");
                Console.WriteLine("  I4_MODP ( I, J ) called with J = " + j + "");
                return 0;
        }

        int value = i % j;

        switch (value)
        {
            case < 0:
                value += Math.Abs(j);
                break;
        }

        return value;
    }

    public static int i4_moebius(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_MOEBIUS returns the value of MU(N), the Moebius function of N.
        //
        //  Discussion:
        //
        //    MU(N) is defined as follows:
        //
        //      MU(N) = 1 if N = 1;
        //              0 if N is divisible by the square of a prime;
        //              (-1)^K, if N is the product of K distinct primes.
        //
        //    The Moebius function MU(D) is related to Euler's totient 
        //    function PHI(N):
        //
        //      PHI(N) = sum ( D divides N ) MU(D) * ( N / D ).
        //
        //    As special cases, MU(N) is -1 if N is a prime, and MU(N) is 0
        //    if N is a square, cube, etc.
        //
        //  Example:
        //
        //     N  MU(N)
        //
        //     1    1
        //     2   -1
        //     3   -1
        //     4    0
        //     5   -1
        //     6    1
        //     7   -1
        //     8    0
        //     9    0
        //    10    1
        //    11   -1
        //    12    0
        //    13   -1
        //    14    1
        //    15    1
        //    16    0
        //    17   -1
        //    18    0
        //    19   -1
        //    20    0
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
        //    Input, int N, the value to be analyzed.
        //
        //    Output, int I4_MOEBIUS, the value of MU(N).
        //    If N is less than or equal to 0, MU will be returned as -2.
        //    If there was not enough internal space for factoring, MU
        //    is returned as -3.
        //
    {
        const int FACTOR_MAX = 20;

        int[] exponent = new int[FACTOR_MAX];
        int[] factor = new int[FACTOR_MAX];
        int i;
        int mu;
        int nfactor = 0;
        int nleft = 0;

        switch (n)
        {
            case <= 0:
                mu = -2;
                return mu;
            case 1:
                mu = 1;
                return mu;
        }

        //
        //  Factor N.
        //
        i4_factor(n, FACTOR_MAX, ref nfactor, ref factor, ref exponent, ref nleft);

        if (nleft != 1)
        {
            Console.WriteLine("");
            Console.WriteLine("I4_MOEBIUS - Fatal error!");
            Console.WriteLine("  Not enough factorization space.");
            mu = -3;
            return mu;
        }

        mu = 1;

        for (i = 0; i < nfactor; i++)
        {
            mu = -mu;

            switch (exponent[i])
            {
                case > 1:
                    mu = 0;
                    return mu;
            }
        }

        return mu;
    }

    public static int i4_mop(int i)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_MOP returns the I-th power of -1 as an I4 value.
        //
        //  Discussion:
        //
        //    An I4 is an int value.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 November 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, the power of -1.
        //
        //    Output, int I4_MOP, the I-th power of -1.
        //
    {
        int value = (i % 2) switch
        {
            0 => 1,
            _ => -1
        };

        return value;
    }

    public static int i4_not(int i, int j)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_NOT calculates the NOT of an I4 with respect to a maximum value.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 January 2016
        //
        //  Author:
        //
        //   John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, the value whose NOT is needed.
        //
        //    Input, int J, the maximum value, such as 255.
        //
        //    Output, unsigned int I4_NOT, the NOT of I with respect to J.
        //
    {
        int k = 0;
        int l = 1;

        while (j != 0)
        {
            int i2 = i / 2;
            int j2 = j / 2;

            if (i == 2 * i2)
            {
                k += l;
            }

            i = i2;
            j = j2;
            l = 2 * l;
        }

        return k;
    }

    public static int i4_or(int i, int j)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_OR calculates the inclusive OR of two I4's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 January 2016
        //
        //  Author:
        //
        //   John Burkardt
        //
        //  Parameters:
        //
        //    Input, unsigned int I, J, two values whose inclusive OR is needed.
        //
        //    Output, unsigned int I4_OR, the inclusive OR of I and J.
        //
    {
        int k = 0;
        int l = 1;

        while (i != 0 || j != 0)
        {
            int i2 = i / 2;
            int j2 = j / 2;

            if (i != 2 * i2 || j != 2 * j2)
            {
                k += l;
            }

            i = i2;
            j = j2;
            l = 2 * l;
        }

        return k;
    }

    public static int i4_partition_distinct_count(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_PARTITION_DISTINCT_COUNT returns any value of Q(N).
        //
        //  Discussion:
        //
        //    A partition of an integer N is a representation of the integer
        //    as the sum of nonzero positive integers.  The order of the summands
        //    does not matter.  The number of partitions of N is symbolized
        //    by P(N).  Thus, the number 5 has P(N) = 7, because it has the
        //    following partitions:
        //
        //    5 = 5
        //      = 4 + 1
        //      = 3 + 2
        //      = 3 + 1 + 1
        //      = 2 + 2 + 1
        //      = 2 + 1 + 1 + 1
        //      = 1 + 1 + 1 + 1 + 1
        //
        //    However, if we require that each member of the partition
        //    be distinct, we are computing something symbolized by Q(N).
        //    The number 5 has Q(N) = 3, because it has the following partitions
        //    into distinct parts:
        //
        //    5 = 5
        //      = 4 + 1
        //      = 3 + 2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 February 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //  Parameters:
        //
        //    Input, int N, the integer to be partitioned.
        //
        //    Output, int I4_PARTITION_DISTINCT_COUNT, the number of partitions 
        //    of the integer into distinct parts.
        //
    {
        int i;

        int[] c = new int[n + 1];

        c[0] = 1;

        for (i = 1; i <= n; i++)
        {
            if (i4_is_triangular(i))
            {
                c[i] = 1;
            }
            else
            {
                c[i] = 0;
            }

            int k = 0;
            int k_sign = -1;

            int k2;
            for (;;)
            {
                k += 1;
                k_sign = -k_sign;
                k2 = k * (3 * k + 1);

                if (i < k2)
                {
                    break;
                }

                c[i] += k_sign * c[i - k2];

            }

            k = 0;
            k_sign = -1;

            for (;;)
            {
                k += 1;
                k_sign = -k_sign;
                k2 = k * (3 * k - 1);

                if (i < k2)
                {
                    break;
                }

                c[i] += k_sign * c[i - k2];

            }

        }

        int value = c[n];

        return value;
    }

    public static int i4_rise(int x, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_RISE computes the rising factorial function [X]^N.
        //
        //  Discussion:
        //
        //    [X}^N = X * ( X + 1 ) * ( X + 2 ) * ... * ( X + N - 1 ).
        //
        //    Note that the number of ways of arranging N objects in M ordered
        //    boxes is [M}^N.  (Here, the ordering in each box matters).  Thus,
        //    2 objects in 2 boxes have the following 6 possible arrangements:
        //
        //      -/12, 1/2, 12/-, -/21, 2/1, 21/-.
        //
        //    Moreover, the number of non-decreasing maps from a set of
        //    N to a set of M ordered elements is [M]^N / N!.  Thus the set of
        //    nondecreasing maps from (1,2,3) to (a,b,c,d) is the 20 elements:
        //
        //      aaa, abb, acc, add, aab, abc, acd, aac, abd, aad
        //      bbb, bcc, bdd, bbc, bcd, bbd, ccc, cdd, ccd, ddd.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X, the argument of the rising factorial function.
        //
        //    Input, int N, the order of the rising factorial function.
        //    If N = 0, RISE = 1, if N = 1, RISE = X.  Note that if N is
        //    negative, a "falling" factorial will be computed.
        //
        //    Output, int I4_RISE, the value of the rising factorial function.
        //
    {
        int i;

        int value = 1;

        switch (n)
        {
            case > 0:
            {
                for (i = 1; i <= n; i++)
                {
                    value *= x;
                    x += 1;
                }

                break;
            }
            case < 0:
            {
                for (i = -1; n <= i; i--)
                {
                    value *= x;
                    x -= 1;
                }

                break;
            }
        }

        return value;
    }

    public static void i4_rise_values(ref int n_data, ref int m, ref int n, ref int fmn)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_RISE_VALUES returns values of the integer rising factorial function.
        //
        //  Discussion:
        //
        //    The integer rising factorial function is sometimes symbolized by (m)_n.
        //
        //    The definition is
        //
        //      (m)_n = (m-1+n)! / (m-1)!
        //            = ( m ) * ( m + 1 ) * ( m + 2 ) ... * ( m - 1 + n )
        //            = Gamma ( m + n ) / Gamma ( m )
        //
        //    We assume 0 <= N <= M.
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Pochhammer[m,n]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 December 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //    Stephen Wolfram,
        //    The Mathematica Book,
        //    Fourth Edition,
        //    Cambridge University Press, 1999,
        //    ISBN: 0-521-64314-7,
        //    LC: QA76.95.W65.
        //
        //  Parameters:
        //
        //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, int &M, &N, the arguments of the function.
        //
        //    Output, int &FMN, the value of the function.
        //
    {
        const int N_MAX = 15;

        int[] fmn_vec =
        {
            1, 5, 30, 210, 1680,
            15120, 151200, 1, 10, 4000,
            110, 6840, 840, 970200, 5040
        };

        int[] m_vec =
        {
            5, 5, 5, 5, 5,
            5, 5, 50, 10, 4000,
            10, 18, 4, 98, 1
        };

        int[] n_vec =
        {
            0, 1, 2, 3, 4,
            5, 6, 0, 1, 1,
            2, 3, 4, 3, 7
        };

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

        if (N_MAX < n_data)
        {
            n_data = 0;
            m = 0;
            n = 0;
            fmn = 0;
        }
        else
        {
            m = m_vec[n_data - 1];
            n = n_vec[n_data - 1];
            fmn = fmn_vec[n_data - 1];
        }
    }

    public static int i4_sign3 ( int i )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_SIGN3 returns the three-way sign of an I4.
        //
        //  Discussion:
        //
        //    This is the "three way" sign of an I4.
        //    The value 0 has the sign of 0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, the integer whose sign is desired.
        //
        //    Output, int I4_SIGN3, the sign of I.
    {
        int value = i switch
        {
            < 0 => -1,
            0 => 0,
            _ => 1
        };
        return value;
    }
        
    public static void i4_swap3 ( ref int i, ref int j, ref int k )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_SWAP3 swaps three I4's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, int &I, &J, &K.  On output, the values of I, J, and K
        //    have been interchanged.
        //
    {
        int l = i;
        i = j;
        j = k;
        k = l;
    }
        
}