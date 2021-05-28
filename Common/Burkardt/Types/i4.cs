using System;

namespace Burkardt.Types
{
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
            i4 ival = new i4() {val = 0};
            char c;
            int i;
            int isgn;
            int istate;

            istate = 0;
            isgn = 1;
            i = 0;

            while (i < s.Length)
            {
                c = s[i];
                i = i + 1;
                //
                //  Haven't read anything.
                //
                if (istate == 0)
                {
                    if (c == ' ')
                    {
                    }
                    else if (c == '-')
                    {
                        istate = 1;
                        isgn = -1;
                    }
                    else if (c == '+')
                    {
                        istate = 1;
                        isgn = +1;
                    }
                    else if ('0' <= c && c <= '9')
                    {
                        istate = 2;
                        ival.val = c - '0';
                    }
                    else
                    {
                        ival.error = true;
                        return ival;
                    }
                }
                //
                //  Have read the sign, expecting digits.
                //
                else if (istate == 1)
                {
                    if (c == ' ')
                    {
                    }
                    else if ('0' <= c && c <= '9')
                    {
                        istate = 2;
                        ival.val = c - '0';
                    }
                    else
                    {
                        ival.error = true;
                        return ival;
                    }
                }
                //
                //  Have read at least one digit, expecting more.
                //
                else if (istate == 2)
                {
                    if ('0' <= c && c <= '9')
                    {
                        ival.val = 10 * (ival.val) + c - '0';
                    }
                    else
                    {
                        ival.val = isgn * ival.val;
                        ival.lchar = i - 1;
                        return ival;
                    }
                }
            }

            //
            //  If we read all the characters in the string, see if we're OK.
            //
            if (istate == 2)
            {
                ival.val = isgn * ival.val;
                ival.lchar = typeMethods.s_len_trim(s);
            }
            else
            {
                ival.error = true;
                ival.lchar = 0;
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
            i4vec ret = new i4vec() {ivec = new int[n]};

            int begin = 0;
            int length = s.Length;

            for (int i = 0; i < n; i++)
            {
                i4 res = s_to_i4(s.Substring(begin, length));

                ret.ivec[i] = res.val;
                int lchar = res.lchar;

                if (res.error)
                {
                    return ret;
                }

                begin = begin + lchar;
                length = length - lchar;
            }

            return ret;
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
            int value;

            int jlo = Math.Min(ilo, ihi);
            int jhi = Math.Max(ilo, ihi);

            int wide = jhi + 1 - jlo;

            if (wide == 1)
            {
                value = jlo;
            }
            else
            {
                value = jlo + Math.Abs((ival - jlo) % wide);
            }

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
                for (int jlo = 0; jlo < m; jlo = jlo + 10)
                {
                    int jhi = Math.Min(jlo + 10, m);
                    Console.WriteLine();
                    string cout = "        J:";
                    for (int j = jlo; j < jhi; j++)
                    {
                        cout += "  " + j.ToString().PadLeft(6);
                    }

                    Console.WriteLine(cout);
                    Console.WriteLine("       I:");
                    for (int i = 0; i < l; i++)
                    {
                        cout = "  " + i.ToString().PadLeft(6) + ": ";
                        for (int j = jlo; j < jhi; j++)
                        {
                            cout += "  " + a[i + j * l + k * l * m].ToString().PadLeft(6);
                        }

                        Console.WriteLine(cout);
                    }
                }
            }
        }

        public static void i4vec_print(int n, int[] a, string title)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_PRINT prints an I4VEC.
            //
            //  Discussion:
            //
            //    An I4VEC is a vector of I4's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 November 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of components of the vector.
            //
            //    Input, int A[N], the vector to be printed.
            //
            //    Input, string TITLE, a title.
            //
        {
            Console.WriteLine("");
            Console.WriteLine(title);
            Console.WriteLine("");
            for (int i = 0; i < n; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(8)
                                       + ": " + a[i].ToString().PadLeft(8) + "");
            }
        }

        public static int i4vec_sum(int n, int[] a)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_SUM sums the entries of an I4VEC.
            //
            //  Discussion:
            //
            //    An I4VEC is a vector of I4's.
            //
            //  Example:
            //
            //    Input:
            //
            //      A = ( 1, 2, 3, 4 )
            //
            //    Output:
            //
            //      I4VEC_SUM = 10
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
            //    Input, int N, the number of entries in the vector.
            //
            //    Input, int A[N], the vector to be summed.
            //
            //    Output, int I4VEC_SUM, the sum of the entries of A.
            //
        {
            int i;

            int sum = 0;
            for (i = 0; i < n; i++)
            {
                sum = sum + a[i];
            }

            return sum;
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
            bool value;

            if (n <= 0)
            {
                value = false;
                return value;
            }

            while (1 < n)
            {
                if ((n % 10) != 0)
                {
                    value = false;
                    return value;
                }

                n = n / 10;
            }

            value = true;

            return value;
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

            if (mn < 0)
            {
                value = 0;
            }
            else if (mn == 0)
            {
                value = 1;
            }
            else
            {
                value = mx + 1;

                for (int i = 2; i <= mn; i++)
                {
                    value = (value * (mx + i)) / i;
                }
            }

            return value;
        }
        
        public static double i4_choose_log ( int n, int k )
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
            double value = Helpers.LogGamma ( ( double ) ( n + 1 ) ) 
                           - Helpers.LogGamma ( ( double ) ( k + 1 ) ) 
                           - Helpers.LogGamma ( ( double ) ( n - k + 1 ) );

            return value;
        }

        public static int i4vec_unique_count(int n, int[] a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_UNIQUE_COUNT counts the unique elements in an unsorted I4VEC.
        //
        //  Discussion:
        //
        //    An I4VEC is a vector of I4's.
        //
        //    Because the array is unsorted, this algorithm is O(N^2).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 April 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of elements of A.
        //
        //    Input, int A[N], the array to examine, which does NOT have to
        //    be sorted.
        //
        //    Output, int I4VEC_UNIQUE_COUNT, the number of unique elements of A.
        //
        {
            int unique_num = 0;

            for (int i = 0; i < n; i++)
            {
                unique_num = unique_num + 1;

                for (int j = 0; j < i; j++)
                {
                    if (a[i] == a[j])
                    {
                        unique_num = unique_num - 1;
                        break;
                    }
                }
            }

            return unique_num;
        }

        public static double i4vec_variance(int n, int[] x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_VARIANCE returns the variance of an I4VEC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 May 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vector.
        //
        //    Input, int X[N], the vector whose variance is desired.
        //
        //    Output, double I4VEC_VARIANCE, the variance of the vector entries.
        //
        {
            double mean = i4vec_mean(n, x);

            double variance = 0.0;
            for (int i = 0; i < n; i++)
            {
                variance = variance + ((double) x[i] - mean) * ((double) x[i] - mean);
            }

            if (1 < n)
            {
                variance = variance / (double) (n - 1);
            }
            else
            {
                variance = 0.0;
            }

            return variance;
        }
        
        public static double i4vec_mean ( int n, int[] x )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_MEAN returns the mean of an I4VEC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 May 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vector.
        //
        //    Input, int X[N], the vector whose mean is desired.
        //
        //    Output, double I4VEC_MEAN, the mean, or average, of the vector entries.
        //
        {
            double mean = 0.0;
            for (int i = 0; i < n; i++ )
            {
                mean = mean + ( double ) x[i];
            }

            mean = mean / ( double ) n;

            return mean;
        }
        
        public static int i4_huge ( )
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
        
        public static int i4vec_run_count ( int n, int[] a )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_RUN_COUNT counts runs of equal values in an I4VEC.
        //
        //  Discussion:
        //
        //    An I4VEC is a vector of integer values.
        //
        //    A run is a sequence of equal values.
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
        //    Input, int N, the number of entries in the vector.
        //
        //    Input, int A[N], the vector to be examined.
        //
        //    Output, int I4VEC_RUN_COUNT, the number of runs.
        //
        {
            int run_count = 0;

            if ( n < 1 )
            {
                return run_count;
            }

            int test = 0;

            for (int i = 0; i < n; i++ )
            {
                if ( i == 0 || a[i] != test )
                {
                    run_count = run_count + 1;
                    test = a[i];
                }
            }

            return run_count;
        }
    }
}