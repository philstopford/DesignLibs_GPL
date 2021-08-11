using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void i4vec_concatenate ( int n1, int[] a, int n2, int[] b, ref int[] c )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_CONCATENATE concatenates two I4VEC's.
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
        //    22 November 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, the number of entries in the first vector.
        //
        //    Input, int A[N1], the first vector.
        //
        //    Input, int N2, the number of entries in the second vector.
        //
        //    Input, int B[N2], the second vector.
        //
        //    Output, int C[N1+N2], the concatenated vector.
        //
        {
            int i;

            for ( i = 0; i < n1; i++ )
            {
                c[i] = a[i];
            }
            for ( i = 0; i < n2; i++ )
            {
                c[n1+i] = b[i];
            }
        }
        
        public static int i4vec_product(int n, int[] a, int aIndex = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_PRODUCT multiplies the entries of an I4VEC.
            //
            //  Example:
            //
            //    Input:
            //
            //      A = ( 1, 2, 3, 4 )
            //
            //    Output:
            //
            //      I4VEC_PRODUCT = 24
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vector.
            //
            //    Input, int A[N], the vector
            //
            //    Output, int I4VEC_PRODUCT, the product of the entries of A.
            //
        {
            int i;
            int product;

            product = 1;
            for (i = 0; i < n; i++)
            {
                product = product * a[i + aIndex];
            }

            return product;
        }

        public static void i4vec_copy(int n, int[] a1, ref int[] a2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_COPY copies an I4VEC.
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
            //    25 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vectors.
            //
            //    Input, int A1[N], the vector to be copied.
            //
            //    Output, int A2[N], the copy of A1.
            //
        {
            int i;

            for (i = 0; i < n; i++)
            {
                a2[i] = a1[i];
            }
        }

        public static void i4vec_indicator0 ( int n, ref int[] a )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_INDICATOR0 sets an I4VEC to the indicator vector.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    25 February 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of elements of A.
            //
            //    Output, int A[N], the initialized array.
            //
        {
            int i;

            for ( i = 0; i < n; i++ )
            {
                a[i] = i;
            }
        }
        public static int[] i4vec_indicator0_new(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_INDICATOR0_NEW sets an I4VEC to the indicator vector (0,1,2,...).
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
            //    27 September 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of elements of A.
            //
            //    Output, int I4VEC_INDICATOR0_NEW[N], the array.
            //
        {
            int[] a;
            int i;

            a = new int[n];

            for (i = 0; i < n; i++)
            {
                a[i] = i;
            }

            return a;
        }

        public static int[] i4vec_indicator1_new(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_INDICATOR1_NEW sets an I4VEC to the indicator vector (1,2,3,...).
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
            //    27 September 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of elements of A.
            //
            //    Output, int I4VEC_INDICATOR1_NEW[N], the array.
            //
        {
            int[] a;
            int i;

            a = new int[n];

            for (i = 0; i < n; i++)
            {
                a[i] = i + 1;
            }

            return a;
        }

        public static void i4vec_permute(int n, int[] p, ref int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_PERMUTE permutes an I4VEC in place.
            //
            //  Discussion:
            //
            //    An I4VEC is a vector of I4's.
            //
            //    This routine permutes an array of integer "objects", but the same
            //    logic can be used to permute an array of objects of any arithmetic
            //    type, or an array of objects of any complexity.  The only temporary
            //    storage required is enough to store a single object.  The number
            //    of data movements made is N + the number of cycles of order 2 or more,
            //    which is never more than N + N/2.
            //
            //  Example:
            //
            //    Input:
            //
            //      N = 5
            //      P = (   1,   3,   4,   0,   2 )
            //      A = (   1,   2,   3,   4,   5 )
            //
            //    Output:
            //
            //      A    = (   2,   4,   5,   1,   3 ).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 October 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of objects.
            //
            //    Input, int P[N], the permutation.  P(I) = J means
            //    that the I-th element of the output array should be the J-th
            //    element of the input array.
            //
            //    Input/output, int A[N], the array to be permuted.
            //
        {
            int a_temp;
            int i;
            int iget;
            int iput;
            int istart;

            if (!perm0_check(n, p))
            {
                Console.WriteLine("");
                Console.WriteLine("I4VEC_PERMUTE - Fatal error!");
                Console.WriteLine("  PERM0_CHECK rejects permutation.");
                return;
            }

            //
            //  In order for the sign negation trick to work, we need to assume that the
            //  entries of P are strictly positive.  Presumably, the lowest number is 0.
            //  So temporarily add 1 to each entry to force positivity.
            //
            for (i = 0; i < n; i++)
            {
                p[i] = p[i] + 1;
            }

            //
            //  Search for the next element of the permutation that has not been used.
            //
            for (istart = 1; istart <= n; istart++)
            {
                if (p[istart - 1] < 0)
                {
                    continue;
                }
                else if (p[istart - 1] == istart)
                {
                    p[istart - 1] = -p[istart - 1];
                    continue;
                }
                else
                {
                    a_temp = a[istart - 1];
                    iget = istart;
                    //
                    //  Copy the new value into the vacated entry.
                    //
                    for (;;)
                    {
                        iput = iget;
                        iget = p[iget - 1];

                        p[iput - 1] = -p[iput - 1];

                        if (iget < 1 || n < iget)
                        {
                            Console.WriteLine("");
                            Console.WriteLine("I4VEC_PERMUTE - Fatal error!");
                            Console.WriteLine("  Entry IPUT = " + iput + " of the permutation has");
                            Console.WriteLine("  an illegal value IGET = " + iget + ".");
                            return;
                        }

                        if (iget == istart)
                        {
                            a[iput - 1] = a_temp;
                            break;
                        }

                        a[iput - 1] = a[iget - 1];
                    }
                }
            }

            //
            //  Restore the signs of the entries.
            //
            for (i = 0; i < n; i++)
            {
                p[i] = -p[i];
            }

            //
            //  Restore the entries.
            //
            for (i = 0; i < n; i++)
            {
                p[i] = p[i] - 1;
            }

            return;
        }

        public static int[] i4vec_zero_new(int n)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_ZERO_NEW creates and zeroes an I4VEC.
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
            //    11 July 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vector.
            //
            //    Output, int I4VEC_ZERO_NEW[N], a vector of zeroes.
            //
        {
            int[] a = new int[n];

            for (int i = 0; i < n; i++)
            {
                a[i] = 0;
            }

            return a;
        }

        public static double i4vec_std ( int n, int[] x )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_STD returns the standard deviation of an I4VEC.
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
            //    14 August 2009
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
            //    Output, double I4VEC_STD, the standard deviation of the vector entries.
            //
        {
            int i;
            double mean;
            double std;

            if ( n < 2 )
            {
                std = 0.0;
            }
            else
            {
                mean = 0.0;
                for ( i = 0; i < n; i++ )
                {
                    mean = mean + ( double ) x[i];
                }
                mean = mean / ( double ) n;

                std = 0.0;
                for ( i = 0; i < n; i++ )
                {
                    std = std + Math.Pow ( ( double ) x[i] - mean, 2 );
                }
                std = Math.Sqrt ( std / ( double ) ( n - 1 ) );
            }

            return std;
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

        public static int i4vec_run_count(int n, int[] a)
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

            if (n < 1)
            {
                return run_count;
            }

            int test = 0;

            for (int i = 0; i < n; i++)
            {
                if (i == 0 || a[i] != test)
                {
                    run_count = run_count + 1;
                    test = a[i];
                }
            }

            return run_count;
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
        
        public static int i4vec_maxloc_last ( int n, int[] x )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_MAXLOC_LAST returns the index of the last maximal I4VEC entry.
            //
            //  Example:
            //
            //    X = ( 5, 1, 2, 5, 0, 5, 3 )
            //
            //    I4VEC_MAXLOC_LAST = 5
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
            //    Input, int N, the size of the array.
            //
            //    Input, int X[N], the array to be examined.
            //
            //    Output, int I4VEC_MAXLOC_LAST, the index of the last element of
            //    X of maximal value.
            //
        {
            int i;
            int index;
            int value;

            index = 0;
            value = x[0];

            for ( i = 1; i < n; i++ )
            {
                if ( value <= x[i] )
                {
                    index = i;
                    value = x[i];
                }
            }
            return index;
        }

        public static double i4vec_mean(int n, int[] x)
        {
            if (x.Length <= 0)
            {
                return 0;
            }

            // Limit to the number of items in the array as a maximum
            n = Math.Min(n, x.Length);

            if (n == x.Length)
            {
                return x.Average();
            }

            return x.Take(n).Average();
        }

        public static int i4vec_max(int n, int[] ivec)
        {
            if (ivec.Length <= 0)
            {
                return 0;
            }

            // Limit to the number of items in the array as a maximum
            n = Math.Min(n, ivec.Length);

            if (n == ivec.Length)
            {
                return ivec.Max();
            }

            return ivec.Take(n).Max();
        }

        public static int i4vec_min(int n, int[] ivec)
        {
            if (ivec.Length <= 0)
            {
                return 0;
            }

            // Limit to the number of items in the array as a maximum
            n = Math.Min(n, ivec.Length);

            if (n == ivec.Length)
            {
                return ivec.Min();
            }

            return ivec.Take(n).Min();
        }

        public static void i4vec_negone ( int n, ref int[] a )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_NEGONE sets an I4VEC to -1.
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
            //    12 October 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vector.
            //
            //    Output, int A[N], a vector of -1's.
            //
        {
            int i;

            for ( i = 0; i < n; i++ )
            {
                a[i] = -1;
            }
        }

        public static int[] i4vec_negone_new ( int n )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_NEGONE_NEW creates an I4VEC and sets it to -1.
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
            //    12 October 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vector.
            //
            //    Output, int I4VEC_NEGONE_NEW[N], a vector of -1's.
            //
        {
            int[] a;
            int i;

            a = new int[n];

            for ( i = 0; i < n; i++ )
            {
                a[i] = -1;
            }
            return a;
        }
        
        public static void i4vec_zero(int n, ref int[] a)
        {
            i4vec_zeros(n, ref a);
        }
        
        public static void i4vec_zeros(int n, ref int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_ZEROS zeroes an I4VEC.
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
            //    01 August 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vector.
            //
            //    Output, int A[N], a vector of zeroes.
            //
        {
            int i;

            for (i = 0; i < n; i++)
            {
                a[i] = 0;
            }

            return;
        }

        public static int[] i4vec_zeros_new(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_ZEROS_NEW creates and zeroes an I4VEC.
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
            //    11 July 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vector.
            //
            //    Output, int I4VEC_ZEROS_NEW[N], a vector of zeroes.
            //
        {
            int[] a;
            int i;

            a = new int[n];

            for (i = 0; i < n; i++)
            {
                a[i] = 0;
            }

            return a;
        }

        public static int i4vec_frac(int n, ref int[] a, int k)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_FRAC searches for the K-th smallest entry in an I4VEC.
            //
            //  Discussion:
            //
            //    An I4VEC is a vector of I4's.
            //
            //    Hoare's algorithm is used.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 September 2005
            //
            //  Parameters:
            //
            //    Input, int N, the number of elements of A.
            //
            //    Input/output, int A[N].
            //    On input, A is the array to search.
            //    On output, the elements of A have been somewhat rearranged.
            //
            //    Input, int K, the fractile to be sought.  If K = 1, the minimum
            //    entry is sought.  If K = N, the maximum is sought.  Other values
            //    of K search for the entry which is K-th in size.  K must be at
            //    least 1, and no greater than N.
            //
            //    Output, double I4VEC_FRAC, the value of the K-th fractile of A.
            //
        {
            int frac;
            int i;
            int iryt;
            int j;
            int left;
            int temp;
            int x;

            if (n <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("I4VEC_FRAC - Fatal error!");
                Console.WriteLine("  Illegal nonpositive value of N = " + n + "");
                return (1);
            }

            if (k <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("I4VEC_FRAC - Fatal error!");
                Console.WriteLine("  Illegal nonpositive value of K = " + k + "");
                return (1);
            }

            if (n < k)
            {
                Console.WriteLine("");
                Console.WriteLine("I4VEC_FRAC - Fatal error!");
                Console.WriteLine("  Illegal N < K, K = " + k + "");
                return (1);
            }

            left = 1;
            iryt = n;

            for (;;)
            {
                if (iryt <= left)
                {
                    frac = a[k - 1];
                    break;
                }

                x = a[k - 1];
                i = left;
                j = iryt;

                for (;;)
                {
                    if (j < i)
                    {
                        if (j < k)
                        {
                            left = i;
                        }

                        if (k < i)
                        {
                            iryt = j;
                        }

                        break;
                    }

                    //
                    //  Find I so that X <= A(I).
                    //
                    while (a[i - 1] < x)
                    {
                        i = i + 1;
                    }

                    //
                    //  Find J so that A(J) <= X.
                    //
                    while (x < a[j - 1])
                    {
                        j = j - 1;
                    }

                    if (i <= j)
                    {
                        temp = a[i - 1];
                        a[i - 1] = a[j - 1];
                        a[j - 1] = temp;
                        i = i + 1;
                        j = j - 1;
                    }
                }
            }

            return frac;
        }

        public static int i4vec_median(int n, ref int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_MEDIAN returns the median of an unsorted I4VEC.
            //
            //  Discussion:
            //
            //    An I4VEC is a vector of I4's.
            //
            //    Hoare's algorithm is used.  The values of the vector are
            //    rearranged by this routine.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 September 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of elements of A.
            //
            //    Input/output, int A[N], the array to search.  On output,
            //    the order of the elements of A has been somewhat changed.
            //
            //    Output, int I4VEC_MEDIAN, the value of the median of A.
            //
        {
            int k;
            int median;

            k = (n + 1) / 2;

            median = i4vec_frac(n, ref a, k);

            return median;
        }

        public static void i4vec_heap_d(int n, ref int[] a, int aIndex = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_HEAP_D reorders an I4VEC into a descending heap.
            //
            //  Discussion:
            //
            //    An I4VEC is a vector of I4's.
            //
            //    A heap is an array A with the property that, for every index J,
            //    A[J] >= A[2*J+1] and A[J] >= A[2*J+2], (as long as the indices
            //    2*J+1 and 2*J+2 are legal).
            //
            //  Diagram:
            //
            //                  A(0)
            //
            //            A(1)         A(2)
            //
            //      A(3)       A(4)  A(5) A(6)
            //
            //    A(7) A(8)  A(9) A(10)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 April 1999
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Albert Nijenhuis, Herbert Wilf,
            //    Combinatorial Algorithms,
            //    Academic Press, 1978, second edition,
            //    ISBN 0-12-519260-6.
            //
            //  Parameters:
            //
            //    Input, int N, the size of the input array.
            //
            //    Input/output, int A[N].
            //    On input, an unsorted array.
            //    On output, the array has been reordered into a heap.
            //
        {
            int i;
            int ifree;
            int key;
            int m;
            //
            //  Only nodes (N/2)-1 down to 0 can be "parent" nodes.
            //
            for (i = (n / 2) - 1; 0 <= i; i--)
            {
                //
                //  Copy the value out of the parent node.
                //  Position IFREE is now "open".
                //
                key = a[aIndex + i];
                ifree = i;

                for (;;)
                {
                    //
                    //  Positions 2*IFREE + 1 and 2*IFREE + 2 are the descendants of position
                    //  IFREE.  (One or both may not exist because they equal or exceed N.)
                    //
                    m = 2 * ifree + 1;
                    //
                    //  Does the first position exist?
                    //
                    if (n <= m)
                    {
                        break;
                    }
                    else
                    {
                        //
                        //  Does the second position exist?
                        //
                        if (m + 1 < n)
                        {
                            //
                            //  If both positions exist, take the larger of the two values,
                            //  and update M if necessary.
                            //
                            if (a[aIndex + m] < a[aIndex + m + 1])
                            {
                                m = m + 1;
                            }
                        }

                        //
                        //  If the large descendant is larger than KEY, move it up,
                        //  and update IFREE, the location of the free position, and
                        //  consider the descendants of THIS position.
                        //
                        if (key < a[aIndex + m])
                        {
                            a[aIndex + ifree] = a[aIndex + m];
                            ifree = m;
                        }
                        else
                        {
                            break;
                        }
                    }
                }

                //
                //  When you have stopped shifting items up, return the item you
                //  pulled out back to the heap.
                //
                a[aIndex + ifree] = key;
            }

            return;
        }
        
        public static void i4vec_increment ( int n, ref int[] v )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_INCREMENT increments an I4VEC.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 January 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the size of the array.
            //
            //    Input/output, int V[N], the array to be incremented.
            //
        {
            int i;

            for ( i = 0; i < n; i++ )
            {
                v[i] = v[i] + 1;
            }
        }

        public static int i4vec_index ( int n, int[] a, int aval )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_INDEX returns the location of the first occurrence of a given value.
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
        //    Input, int N, the number of entries in the vector.
        //
        //    Input, int A[N], the vector to be searched.
        //
        //    Input, int AVAL, the value to be indexed.
        //
        //    Output, int I4VEC_INDEX, the first location in A which has the
        //    value AVAL, or -1 if no such index exists.
        //
        {
            int i;

            for ( i = 0; i < n; i++ )
            {
                if ( a[i] == aval )
                {
                    return i;
                }
            }
            return -1;
        }

        public static int[] i4vec_indicator(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_INDICATOR sets an I4VEC to the indicator vector.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    25 February 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of elements of A.
            //
            //    Output, int I4VEC_INDICATOR(N), the initialized array.
            //
        {
            int[] a;
            int i;

            a = new int[n];

            for (i = 0; i < n; i++)
            {
                a[i] = i + 1;
            }

            return a;
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
            int i;

            Console.WriteLine("");
            Console.WriteLine(title + "");
            Console.WriteLine("");
            for (i = 0; i < n; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(8)
                                       + ": " + a[i].ToString().PadLeft(8) + "");
            }
        }

        public static int[] i4vec_sort_heap_index_a(int n, int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an I4VEC.
            //
            //  Discussion:
            //
            //    An I4VEC is a vector of I4's.
            //
            //    The sorting is not actually carried out.  Rather an index array is
            //    created which defines the sorting.  This array may be used to sort
            //    or index the array, or to sort or index related arrays keyed on the
            //    original array.
            //
            //    Once the index array is computed, the sorting can be carried out
            //    "implicitly:
            //
            //      a(indx(*))
            //
            //    or explicitly, by the call
            //
            //      i4vec_permute ( n, indx, a )
            //
            //    after which a(*) is sorted.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the array.
            //
            //    Input, int A[N], an array to be index-sorted.
            //
            //    Output, int I4VEC_SORT_HEAP_INDEX_A[N], contains the sort index.  The
            //    I-th element of the sorted array is A(INDX(I)).
            //
        {
            int aval;
            int i;
            int[] indx;
            int indxt;
            int ir;
            int j;
            int l;

            if (n < 1)
            {
                return null;
            }

            indx = new int[n];

            for (i = 0; i < n; i++)
            {
                indx[i] = i;
            }

            if (n == 1)
            {
                indx[0] = indx[0];
                return indx;
            }

            l = n / 2 + 1;
            ir = n;

            for (;;)
            {

                if (1 < l)
                {
                    l = l - 1;
                    indxt = indx[l - 1];
                    aval = a[indxt];
                }
                else
                {
                    indxt = indx[ir - 1];
                    aval = a[indxt];
                    indx[ir - 1] = indx[0];
                    ir = ir - 1;

                    if (ir == 1)
                    {
                        indx[0] = indxt;
                        break;
                    }
                }

                i = l;
                j = l + l;

                while (j <= ir)
                {
                    if (j < ir)
                    {
                        if (a[indx[j - 1]] < a[indx[j]])
                        {
                            j = j + 1;
                        }
                    }

                    if (aval < a[indx[j - 1]])
                    {
                        indx[i - 1] = indx[j - 1];
                        i = j;
                        j = j + j;
                    }
                    else
                    {
                        j = ir + 1;
                    }
                }

                indx[i - 1] = indxt;
            }

            return indx;
        }

        public static void i4vec_sort_bubble_a ( int n, ref int[] a )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_SORT_BUBBLE_A ascending sorts an integer array using bubble sort.
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
            //    Input, int N, the number of entries in the array.
            //
            //    Input/output, int A[N].
            //    On input, the array to be sorted;
            //    On output, the array has been sorted.
            //
        {
            int i;
            int j;
            int k;

            for ( i = 0; i < n-1; i++ )
            {
                for ( j = i+1; j < n; j++ )
                {
                    if ( a[j] < a[i] )
                    {
                        k = a[i];
                        a[i] = a[j];
                        a[j] = k;
                    }
                }
            }
        }
        
        public static void i4vec_sort_heap_a(int n, ref int[] a, int aIndex = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_SORT_HEAP_A ascending sorts an I4VEC using heap sort.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 April 1999
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    A Nijenhuis and H Wilf,
            //    Combinatorial Algorithms,
            //    Academic Press, 1978, second edition,
            //    ISBN 0-12-519260-6.
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the array.
            //
            //    Input/output, int A[N].
            //    On input, the array to be sorted;
            //    On output, the array has been sorted.
            //
        {
            int n1;
            int temp;

            if (n <= 1)
            {
                return;
            }

            //
            //  1: Put A into descending heap form.
            //
            i4vec_heap_d(n, ref a);
            //
            //  2: Sort A.
            //
            //  The largest object in the heap is in A[0].
            //  Move it to position A[N-1].
            //
            temp = a[aIndex + 0];
            a[aIndex + 0] = a[aIndex + n - 1];
            a[aIndex + n - 1] = temp;
            //
            //  Consider the diminished heap of size N1.
            //
            for (n1 = n - 1; 2 <= n1; n1--)
            {
                //
                //  Restore the heap structure of the initial N1 entries of A.
                //
                i4vec_heap_d(n1, ref a);
                //
                //  Take the largest object from A[0] and move it to A[N1-1].
                //
                temp = a[aIndex + 0];
                a[aIndex + 0] = a[aIndex + n1 - 1];
                a[aIndex + n1 - 1] = temp;

            }

            return;
        }

        public static int i4vec_sorted_unique(int n, ref int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_SORTED_UNIQUE finds the unique elements in a sorted I4VEC.
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
            //    24 August 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of elements in A.
            //
            //    Input/output, int A[N].  On input, the sorted
            //    integer array.  On output, the unique elements in A.
            //
            //    Output, int I4VEC_SORTED_UNIQUE, the number of unique elements in A.
            //
        {
            int i;
            int unique_num;

            unique_num = 0;

            if (n <= 0)
            {
                return unique_num;
            }

            unique_num = 1;

            for (i = 1; i < n; i++)
            {
                if (a[i] != a[unique_num - 1])
                {
                    unique_num = unique_num + 1;
                    a[unique_num - 1] = a[i];
                }
            }

            return unique_num;
        }
        
        public static int i4vec_sorted_unique_count ( int n, int[] a )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_SORTED_UNIQUE_COUNT counts unique elements in a sorted I4VEC.
            //
            //  Discussion:
            //
            //    An I4VEC is a vector of I4's.
            //
            //    Because the array is sorted, this algorithm is O(N).
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
            //    Input, int A[N], the sorted array to examine.
            //
            //    Output, int I4VEC_SORTED_UNIQUE_COUNT, the number of unique elements of A.
            //
        {
            int i;
            int unique_num;

            unique_num = 0;

            if ( n < 1 )
            {
                return unique_num;
            }

            unique_num = 1;

            for ( i = 1; i < n; i++ )
            {
                if ( a[i-1] != a[i] )
                {
                    unique_num = unique_num + 1;
                }
            }

            return unique_num;
        }

        public static int i4vec2_compare(int n, int[] a1, int[] a2, int i, int j)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC2_COMPARE compares pairs of integers stored in two I4VECs.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 September 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of data items.
            //
            //    Input, int A1[N], A2[N], contain the two components of each item.
            //
            //    Input, int I, J, the items to be compared.  These values will be
            //    1-based indices for the arrays A1 and A2.
            //
            //    Output, int I4VEC2_COMPARE, the results of the comparison:
            //    -1, item I < item J,
            //     0, item I = item J,
            //    +1, item J < item I.
            //
        {
            int isgn;

            isgn = 0;

            if (a1[i - 1] < a1[j - 1])
            {
                isgn = -1;
            }
            else if (a1[i - 1] == a1[j - 1])
            {
                if (a2[i - 1] < a2[j - 1])
                {
                    isgn = -1;
                }
                else if (a2[i - 1] < a2[j - 1])
                {
                    isgn = 0;
                }
                else if (a2[j - 1] < a2[i - 1])
                {
                    isgn = +1;
                }
            }
            else if (a1[j - 1] < a1[i - 1])
            {
                isgn = +1;
            }

            return isgn;
        }

        public static void i4vec2_sort_a(int n, ref int[] a1, ref int[] a2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC2_SORT_A ascending sorts a vector of pairs of integers.
            //
            //  Discussion:
            //
            //    Each item to be sorted is a pair of integers (I,J), with the I
            //    and J values stored in separate vectors A1 and A2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 September 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of items of data.
            //
            //    Input/output, int A1[N], A2[N], the data to be sorted..
            //
        {
            int i;
            int indx;
            int isgn;
            int j;
            int temp;
            //
            //  Initialize.
            //
            i = 0;
            indx = 0;
            isgn = 0;
            j = 0;
            //
            //  Call the external heap sorter.
            //
            SortHeapExternalData data = new SortHeapExternalData();
            for (;;)
            {
                Helpers.sort_heap_external(ref data, n, ref indx, ref i, ref j, isgn);
                //
                //  Interchange the I and J objects.
                //
                if (0 < indx)
                {
                    temp = a1[i - 1];
                    a1[i - 1] = a1[j - 1];
                    a1[j - 1] = temp;

                    temp = a2[i - 1];
                    a2[i - 1] = a2[j - 1];
                    a2[j - 1] = temp;
                }
                //
                //  Compare the I and J objects.
                //
                else if (indx < 0)
                {
                    isgn = i4vec2_compare(n, a1, a2, i, j);
                }
                else if (indx == 0)
                {
                    break;
                }

            }
        }

        public static void i4vec2_sorted_unique(int n, ref int[] a1, ref int[] a2, ref int nuniq)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC2_SORTED_UNIQUE finds unique elements in a sorted I4VEC2.
            //
            //  Discussion:
            //
            //    Item I is stored as the pair A1(I), A2(I).
            //
            //    The items must have been sorted, or at least it must be the
            //    case that equal items are stored in adjacent vector locations.
            //
            //    If the items were not sorted, then this routine will only
            //    replace a string of equal values by a single representative.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    09 July 2000
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of items.
            //
            //    Input/output, int A1[N], A2[N].
            //    On input, the array of N items.
            //    On output, an array of NUNIQ unique items.
            //
            //    Output, int *NUNIQ, the number of unique items.
            //
        {
            int itest;

            nuniq = 0;

            if (n <= 0)
            {
                return;
            }

            nuniq = 1;

            for (itest = 1; itest < n; itest++)
            {
                if (a1[itest] != a1[nuniq - 1] ||
                    a2[itest] != a2[nuniq - 1])
                {
                    a1[nuniq] = a1[itest];
                    a2[nuniq] = a2[itest];
                    nuniq = nuniq + 1;
                }
            }

            return;
        }

        public static void i4vec_dec(int n, ref int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_DEC decrements an I4VEC.
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
            //    18 July 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vectors.
            //
            //    Input/output, int A[N], the vector to be decremented.
            //
        {
            for (int i = 0; i < n; i++)
            {
                a[i] = a[i] - 1;
            }
        }
        
        public static int[] i4vec_histogram ( int n, int[] a, int histo_num )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_HISTOGRAM computes a histogram of the elements of an I4VEC.
        //
        //  Discussion:
        //
        //    It is assumed that the entries in the vector A are nonnegative.
        //    Only values between 0 and HISTO_NUM will be histogrammed.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of elements of A.
        //
        //    Input, int A[N], the array to examine.
        //
        //    Input, int HISTO_NUM, the maximum value for which a
        //    histogram entry will be computed.
        //
        //    Output, int I4VEC_HISTOGRAM[HISTO_NUM+1], contains the number of
        //    entries of A with the values of 0 through HISTO_NUM.
        //
        {
            int[] histo_gram;
            int i;

            histo_gram = new int[histo_num+1];

            for ( i = 0; i <= histo_num; i++ )
            {
                histo_gram[i] = 0;
            }

            for ( i = 0; i < n; i++ )
            {
                if ( 0 <= a[i] && a[i] <= histo_num )
                {
                    histo_gram[a[i]] = histo_gram[a[i]] + 1;
                }
            }

            return histo_gram;
        }

        public static void i4vec_inc(int n, ref int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_INC increments an I4VEC.
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
            //    18 July 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vectors.
            //
            //    Input/output, int A[N], the vector to be incremented.
            //
        {
            for (int i = 0; i < n; i++)
            {
                a[i] = a[i] + 1;
            }
        }
        
        public static int[] i4vec_cum ( int n, int[] a )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_CUM computes the cumulutive sum of the entries of an I4VEC.
            //
            //  Discussion:
            //
            //    An I4VEC is a vector of I4's.
            //
            //  Example:
            //
            //    Input:
            //
            //      A = (/ 1, 2, 3, 4 /)
            //
            //    Output:
            //
            //      A_CUM = (/ 1, 3, 6, 10 /)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 December 2010
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
            //    Output, int I4VEC_CUM[N], the cumulative sum of the entries of A.
            //
        {
            int[] a_cum;
            int i;

            a_cum = new int[n+1];

            a_cum[0] = a[0];

            for ( i = 1; i <= n - 1; i++ )
            {
                a_cum[i] = a_cum[i-1] + a[i];
            }

            return a_cum;
        }

        public static void i4vec_data_read(string input_filename, int n, ref int[] table)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_DATA_READ reads data from an I4VEC file.
            //
            //  Discussion:
            //
            //    An I4VEC is a vector of I4's.
            //
            //    The file is assumed to contain one record per line.
            //
            //    Records beginning with '#' are comments, and are ignored.
            //    Blank lines are also ignored.
            //
            //    Each line that is not ignored is assumed to contain exactly one value.
            //
            //    There are assumed to be exactly (or at least) N such records.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 July 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string INPUT_FILENAME, the name of the input file.
            //
            //    Input, int N, the number of points.  The program
            //    will stop reading data once N values have been read.
            //
            //    Output, int TABLE[N], the data.
            //
        {

            try
            {
                string[] input = File.ReadAllLines(input_filename);

                int j = 0;

                foreach (string line in input)
                {
                    if (line[0] == '#' || s_len_trim(line) == 0)
                    {
                        continue;
                    }

                    table[j] = Convert.ToInt32(line);
                    j = j + 1;
                }
            }
            catch (Exception e)
            {
                Console.WriteLine("");
                Console.WriteLine("I4VEC_DATA_READ - Fatal error!");
                Console.WriteLine("  Could not open the input file: \"" + input_filename + "\"");
            }

        }

        public static void i4vec_write(string output_filename, int n, int[] table)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_WRITE writes an I4VEC to a file.
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
            //    12 July 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string OUTPUT_FILENAME, the output filename.
            //
            //    Input, int N, the number of points.
            //
            //    Input, int TABLE[N], the data.
            //
        {
            try
            {
                List<string> output = new List<string>();

                for (int j = 0; j < n; j++)
                {
                    output.Add(table[j].ToString());
                }

                File.WriteAllLines(output_filename, output);
            }
            catch (Exception e)
            {
                Console.WriteLine("");
                Console.WriteLine("I4VEC_WRITE - Fatal error!");
                Console.WriteLine("  Could not open the output file.");
            }
        }

        public static int i4vec2_sorted_unique_count(int n, int[] a1, int[] a2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC2_SORTED_UNIQUE_COUNT counts unique elements in an I4VEC2.
            //
            //  Discussion:
            //
            //    Item I is stored as the pair A1(I), A2(I).
            //
            //    The items must have been sorted, or at least it must be the
            //    case that equal items are stored in adjacent vector locations.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 July 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of items.
            //
            //    Input, int A1[N], A2[N], the array of N items.
            //
            //    Output, int I4VEC2_SORTED_UNIQUE_COUNT, the number of unique items.
            //
        {
            int i;
            int iu;
            int unique_num;

            unique_num = 0;

            if (n <= 0)
            {
                return unique_num;
            }

            iu = 0;
            unique_num = 1;

            for (i = 1; i < n; i++)
            {
                if (a1[i] != a1[iu] ||
                    a2[i] != a2[iu])
                {
                    iu = i;
                    unique_num = unique_num + 1;
                }
            }

            return unique_num;
        }

        public static void i4vec2_sorted_uniquely(int n1, int[] a1, int[] b1, int n2, ref int[] a2,
                ref int[] b2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC2_SORTED_UNIQUELY keeps the unique elements in an I4VEC2.
            //
            //  Discussion:
            //
            //    Item I is stored as the pair A1(I), A2(I).
            //
            //    The items must have been sorted, or at least it must be the
            //    case that equal items are stored in adjacent vector locations.
            //
            //    If the items were not sorted, then this routine will only
            //    replace a string of equal values by a single representative.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 July 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N1, the number of items.
            //
            //    Input, int A1[N1], B1[N1], the input array.
            //
            //    Input, int N2, the number of unique items.
            //
            //    Input, int A2[N2], B2[N2], the output array of unique items.
            //
        {
            int i1;
            int i2;

            i1 = 0;
            i2 = 0;

            if (n1 <= 0)
            {
                return;
            }

            a2[i2] = a1[i1];
            b2[i2] = b1[i1];

            for (i1 = 1; i1 < n1; i1++)
            {
                if (a1[i1] != a2[i2] || b1[i1] != b2[i2])
                {
                    i2 = i2 + 1;
                    a2[i2] = a1[i1];
                    b2[i2] = b1[i1];
                }
            }
        }

        public static int[] i4vec_copy_new(int n, int[] a1)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_COPY_NEW copies an I4VEC.
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
            //    04 July 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vectors.
            //
            //    Input, int A1[N], the vector to be copied.
            //
            //    Output, int I4VEC_COPY_NEW[N], the copy of A1.
            //
        {
            int[] a2 = new int[n];

            for (int i = 0; i < n; i++)
            {
                a2[i] = a1[i];
            }

            return a2;
        }

        public static void i4vec_transpose_print(int n, int[] a, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_TRANSPOSE_PRINT prints an I4VEC "transposed".
            //
            //  Discussion:
            //
            //    An I4VEC is a vector of I4's.
            //
            //  Example:
            //
            //    A = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 }
            //    TITLE = "My vector:  "
            //
            //    My vector:      1    2    3    4    5
            //                    6    7    8    9   10
            //                   11
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 July 2004
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
            //    Input, string TITLE, a title to be printed first.
            //    TITLE may be blank or NULL.
            //
        {
            int i;
            int ihi;
            int ilo;
            int title_len;

            int slt = s_len_trim(title);
            
            if (0 < slt)
            {
                title_len = title.Length;

                for (ilo = 1; ilo <= n; ilo = ilo + 5)
                {
                    string cout = "";
                    ihi = Math.Min(ilo + 5 - 1, n);
                    if (ilo == 1)
                    {
                        cout = title;
                    }
                    else
                    {
                        for (i = 1; i <= title_len; i++)
                        {
                            cout += " ";
                        }
                    }

                    for (i = ilo; i <= ihi; i++)
                    {
                        cout += a[i - 1].ToString().PadLeft(12);
                    }

                    Console.WriteLine(cout);
                }
            }
            else
            {
                for (ilo = 1; ilo <= n; ilo = ilo + 5)
                {
                    string cout = "";
                    ihi = Math.Min(ilo + 5 - 1, n);
                    for (i = ilo; i <= ihi; i++)
                    {
                        cout += a[i - 1].ToString().PadLeft(12);
                    }

                    Console.WriteLine(cout);
                }
            }
        }
        
        public static void i4vec0_print ( int n, int[] a, string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC0_PRINT prints an integer vector (0-based indices).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 June 2003
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
            int i;

            if ( 0 < title.Length )
            {
                Console.WriteLine();
                Console.WriteLine(title);
            }

            for ( i = 0; i <= n-1; i++ ) 
            {
                Console.WriteLine(i.ToString().PadLeft(6)    + "  " 
                    + a[i].ToString().PadLeft(8));
            }
        }

        public static void i4vec1_print ( int n, int[] a, string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC1_PRINT prints an integer vector (one-based indices).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 June 2003
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
            int i;

            if ( 0 < title.Length )
            {
                Console.WriteLine();
                Console.WriteLine(title);
            }

            for ( i = 0; i <= n-1; i++ ) 
            {
                Console.WriteLine("  "
                                  + (i+1).ToString().PadLeft(6)  + "  " 
                                  + a[i].ToString().PadLeft(8));
            }

        }
        
        public static bool i4vec_ascends ( int n, ref int[] x )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_ASCENDS determines if an I4VEC is (weakly) ascending.
            //
            //  Example:
            //
            //    X = ( -8, 1, 2, 3, 7, 7, 9 )
            //
            //    I4VEC_ASCENDS = TRUE
            //
            //    The sequence is not required to be strictly ascending.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 April 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the size of the array.
            //
            //    Input, int X[N], the array to be examined.
            //
            //    Output, bool I4VEC_ASCENDS, is TRUE if the entries of X ascend.
            //
        {
            int i;

            for ( i = 1; i <= n - 1; i++ )
            {
                if ( x[i] < x[i-1] )
                {
                    return false;
                }
            }
            return true;
        }

        public static void i4vec_backtrack(int n, int maxstack, int[] stack, ref int[] x, ref int indx,
                ref int k, ref int nstack, ref int[] ncan)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_BACKTRACK supervises a backtrack search for an I4VEC.
            //
            //  Discussion:
            //
            //    The routine tries to construct an integer vector one index at a time,
            //    using possible candidates as supplied by the user.
            //
            //    At any time, the partially constructed vector may be discovered to be
            //    unsatisfactory, but the routine records information about where the
            //    last arbitrary choice was made, so that the search can be
            //    carried out efficiently, rather than starting out all over again.
            //
            //    First, call the routine with INDX = 0 so it can initialize itself.
            //
            //    Now, on each return from the routine, if INDX is:
            //      1, you've just been handed a complete candidate vector;
            //         Admire it, analyze it, do what you like.
            //      2, please determine suitable candidates for position X(K).
            //         Return the number of candidates in NCAN(K), adding each
            //         candidate to the end of STACK, and increasing NSTACK.
            //      3, you're done.  Stop calling the routine;
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 July 2004
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
            //    C++ version by John Burkardt.
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
            //    Input, int N, the number of positions to be filled in the vector.
            //
            //    Input, int MAXSTACK, the maximum length of the stack.
            //
            //    Input, int STACK[MAXSTACK], a list of all current candidates for
            //    all positions 1 through K.
            //
            //    Input/output, int X[N], the partial or complete candidate vector.
            //
            //    Input/output, int &INDX, a communication flag.
            //    On input,
            //      0 to start a search.
            //    On output:
            //      1, a complete output vector has been determined and returned in X(1:N);
            //      2, candidates are needed for position X(K);
            //      3, no more possible vectors exist.
            //
            //    Input/output, int *&K, if INDX=2, the current vector index being considered.
            //
            //    Input/output, int &NSTACK, the current length of the stack.
            //
            //    Input/output, int NCAN[N], lists the current number of candidates for
            //    positions 1 through K.
            //
        {
            //
            //  If this is the first call, request a candidate for position 1.
            //
            if (indx == 0)
            {
                k = 1;
                nstack = 0;
                indx = 2;
                return;
            }

            //
            //  Examine the stack.
            //
            for (;;)
            {
                //
                //  If there are candidates for position K, take the first available
                //  one off the stack, and increment K.
                //
                //  This may cause K to reach the desired value of N, in which case
                //  we need to signal the user that a complete set of candidates
                //  is being returned.
                //
                if (0 < ncan[k - 1])
                {
                    x[k - 1] = stack[nstack - 1];
                    nstack = nstack - 1;

                    ncan[k - 1] = ncan[k - 1] - 1;

                    if (k != n)
                    {
                        k = k + 1;
                        indx = 2;
                    }
                    else
                    {
                        indx = 1;
                    }

                    break;
                }
                //
                //  If there are no candidates for position K, then decrement K.
                //  If K is still positive, repeat the examination of the stack.
                //
                else
                {
                    k = k - 1;

                    if (k <= 0)
                    {
                        indx = 3;
                        break;
                    }
                }
            }

            return;
        }
        
        public static bool i4vec_pairwise_prime ( int n, int[] a )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_PAIRWISE_PRIME checks whether an I4VEC is pairwise prime.
            //
            //  Discussion:
            //
            //    Two positive integers I and J are pairwise prime if they have no common
            //    factor greater than 1.
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
            //    Input, int N, the number of values to check.
            //
            //    Input, int A[N], the vector of integers.
            //
            //    Output, bool I4VEC_PAIRWISE_PRIME, is TRUE if the vector of integers
            //    is pairwise prime.
            //
        {
            int i;
            int j;

            for ( i = 0; i < n; i++ )
            {
                for ( j = i+1; j < n; j++ )
                {
                    if ( i4_gcd ( a[i], a[j] ) != 1 ) 
                    {
                        return false;
                    }
                }
            }
            return true;
        }

        public static int[] i4vec_part1_new(int n, int npart)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    I4VEC_PART1_NEW partitions an integer N into NPART parts.
            // 
            //  Example:
            // 
            //    Input:
            // 
            //      N = 17, NPART = 5
            // 
            //    Output:
            // 
            //      X = ( 13, 1, 1, 1, 1 ).
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
            //    Input, int N, the integer to be partitioned.  N
            //    may be positive, zero, or negative.
            // 
            //    Input, int NPART, the number of entries in the array.
            //    1 <= NPART <= N.
            // 
            //    Output, int I4VEC_PART1_NEW[NPART], the partition of N.  The entries of
            //    X add up to N.  X(1) = N + 1 - NPART, and all other entries
            //    are equal to 1.
            // 
        {
            int i;
            int[] x;

            if (npart < 1 || n < npart)
            {
                Console.WriteLine("");
                Console.WriteLine("I4VEC_PART1_NEW - Fatal error!");
                Console.WriteLine("  The input value of NPART is illegal.");
                return null;
            }

            x = new int[npart];

            x[0] = n + 1 - npart;
            for (i = 1; i < npart; i++)
            {
                x[i] = 1;
            }

            return x;
        }

        public static void i4vec_part2(int n, int npart, ref int[] x, int xIndex = 0)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    I4VEC_PART2 partitions an integer N into NPART nearly equal parts.
            // 
            //  Discussion:
            //
            //    Thanks to John Nitao for pointing out a typographical error in 
            //    of the form "x[j=1]" for "x[j-1]", 14 April 2013.
            //
            //  Example:
            // 
            //    Input:
            // 
            //      N = 17, NPART = 5
            // 
            //    Output:
            // 
            //      X = ( 4, 4, 3, 3, 3 ).
            // 
            //  Licensing:
            // 
            //    This code is distributed under the GNU LGPL license.
            // 
            //  Modified:
            // 
            //    14 April 2013
            // 
            //  Author:
            // 
            //    John Burkardt
            // 
            //  Parameters:
            // 
            //    Input, int N, the integer to be partitioned.  N
            //    may be positive, zero, or negative.
            // 
            //    Input, int NPART, the number of entries in the array.
            //    1 <= NPART
            // 
            //    Output, int X[NPART], the partition of N.  The entries of
            //    X add up to N.  The entries of X are either all equal, or
            //    differ by at most 1.  The entries of X all have the same sign
            //    as N, and the "largest" entries occur first.
            // 
        {
            int i;
            int j;

            if (npart < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("I4VEC_PART2 - Fatal error!");
                Console.WriteLine("  The input value of NPART is illegal.");
                return;
            }

            for (i = 0; i < npart; i++)
            {
                x[i + xIndex] = 0;
            }

            if (0 < n)
            {
                j = 1;
                for (i = 1; i <= n; i++)
                {
                    x[xIndex + (j - 1)] = x[xIndex + (j - 1)] + 1;
                    j = j + 1;
                    if (npart < j)
                    {
                        j = 1;
                    }
                }
            }
            else if (n < 0)
            {
                j = 1;
                for (i = n; i <= -1; i++)
                {
                    x[xIndex + (j - 1)] = x[xIndex + (j - 1)] - 1;
                    j = j + 1;
                    if (npart < j)
                    {
                        j = 1;
                    }
                }
            }
        }

        public static int[] i4vec_part2_new(int n, int npart)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    I4VEC_PART2_NEW partitions an integer N into NPART nearly equal parts.
            // 
            //  Example:
            // 
            //    Input:
            // 
            //      N = 17, NPART = 5
            // 
            //    Output:
            // 
            //      X = ( 4, 4, 3, 3, 3 ).
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
            //  Parameters:
            // 
            //    Input, int N, the integer to be partitioned.  N
            //    may be positive, zero, or negative.
            // 
            //    Input, int NPART, the number of entries in the array.
            //    1 <= NPART
            // 
            //    Output, int I4VEC_PART2[NPART], the partition of N.  The entries of
            //    X add up to N.  The entries of X are either all equal, or
            //    differ by at most 1.  The entries of X all have the same sign
            //    as N, and the "largest" entries occur first.
            // 
        {
            int[] x;

            x = new int[npart];

            i4vec_part2(n, npart, ref x);

            return x;
        }

        public static void i4vec_reverse(int n, ref int[] a, int aIndex = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_REVERSE reverses the elements of an I4VEC.
            //
            //  Discussion:
            //
            //    An I4VEC is a vector of I4's.
            //
            //  Example:
            //
            //    Input:
            //
            //      N = 5,
            //      A = ( 11, 12, 13, 14, 15 ).
            //
            //    Output:
            //
            //      A = ( 15, 14, 13, 12, 11 ).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 September 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the array.
            //
            //    Input/output, int A[N], the array to be reversed.
            //
        {
            int i;
            int j;

            for (i = 0; i < n / 2; i++)
            {
                j = a[aIndex + i];
                a[aIndex + i] = a[aIndex + (n - 1 - i)];
                a[aIndex + (n - 1 - i)] = j;
            }
        }

        public static int i4vec_search_binary_a(int n, int[] a, int b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_SEARCH_BINARY_A searches an ascending sorted I4VEC for a value.
        //
        //  Discussion:
        //
        //    An I4VEC is a vector of I4's.
        //
        //    Binary search is used.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Donald Kreher, Douglas Simpson,
        //    Algorithm 1.9,
        //    Combinatorial Algorithms,
        //    CRC Press, 1998, page 26.
        //
        //  Parameters:
        //
        //    Input, int N, the number of elements in the vector.
        //
        //    Input, int A[N], the array to be searched.  A must
        //    be sorted in ascending order.
        //
        //    Input, int B, the value to be searched for.
        //
        //    Output, int I4VEC_SEARCH_BINARY_A, the result of the search.
        //    -1, B does not occur in A.
        //    I, A[I] = B.
        //
        {
            int high;
            int index;
            int low;
            int mid;
            //
            //  Check.
            //
            if (n <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("I4VEC_SEARCH_BINARY_A - Fatal error!");
                Console.WriteLine("  The array dimension N is less than 1.");
                return (1);
            }

            index = -1;

            low = 1;
            high = n;

            while (low <= high)
            {
                mid = (low + high) / 2;

                if (a[mid - 1] == b)
                {
                    index = mid;
                    break;
                }
                else if (a[mid - 1] < b)
                {
                    low = mid + 1;
                }
                else if (b < a[mid - 1])
                {
                    high = mid - 1;
                }
            }

            return index;
        }

        public static int i4vec_search_binary_d(int n, int[] a, int b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_SEARCH_BINARY_D searches a descending sorted I4VEC for a value.
        //
        //  Discussion:
        //
        //    An I4VEC is a vector of I4's.
        //
        //    Binary search is used.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Donald Kreher, Douglas Simpson,
        //    Algorithm 1.9,
        //    Combinatorial Algorithms,
        //    CRC Press, 1998, page 26.
        //
        //  Parameters:
        //
        //    Input, int N, the number of elements in the vector.
        //
        //    Input, int A[N], the array to be searched.  A must
        //    be sorted in descending order.
        //
        //    Input, int B, the value to be searched for.
        //
        //    Output, int I4VEC_SEARCH_BINARY_D, the result of the search.
        //    -1, B does not occur in A.
        //    I, A[I] = B.
        //
        {
            int high;
            int index;
            int low;
            int mid;
            //
            //  Check.
            //
            if (n <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("I4VEC_SEARCH_BINARY_D - Fatal error!");
                Console.WriteLine("  The array dimension N is less than 1.");
                return 1;
            }

            index = -1;

            low = 1;
            high = n;

            while (low <= high)
            {
                mid = (low + high) / 2;

                if (a[mid - 1] == b)
                {
                    index = mid;
                    break;
                }
                else if (b < a[mid - 1])
                {
                    low = mid + 1;
                }
                else if (a[mid - 1] < b)
                {
                    high = mid - 1;
                }
            }

            return index;
        }

        public static int[] i4vec_sort_heap_index_d(int n, int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_SORT_HEAP_INDEX_D does an indexed heap descending sort of an I4VEC.
            //
            //  Discussion:
            //
            //    The sorting is not actually carried out.  Rather an index array is
            //    created which defines the sorting.  This array may be used to sort
            //    or index the array, or to sort or index related arrays keyed on the
            //    original array.
            //
            //    Once the index array is computed, the sorting can be carried out
            //    "implicitly:
            //
            //      A(INDX(I)), I = 1 to N is sorted,
            //
            //    after which A(I), I = 1 to N is sorted.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 April 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the array.
            //
            //    Input, int A[N], an array to be index-sorted.
            //
            //    Output, int I4VEC_SORT_HEAP_INDEX_D[N], contains the sort index.  The
            //    I-th element of the sorted array is A(INDX(I)).
            //
        {
            int aval;
            int i;
            int[] indx;
            int indxt;
            int ir;
            int j;
            int l;

            indx = i4vec_indicator1_new(n);

            l = n / 2 + 1;
            ir = n;

            for (;;)
            {
                if (1 < l)
                {
                    l = l - 1;
                    indxt = indx[l - 1];
                    aval = a[indxt - 1];
                }
                else
                {
                    indxt = indx[ir - 1];
                    aval = a[indxt - 1];
                    indx[ir - 1] = indx[0];
                    ir = ir - 1;

                    if (ir == 1)
                    {
                        indx[0] = indxt;
                        for (i = 0; i < n; i++)
                        {
                            indx[i] = indx[i] - 1;
                        }

                        return indx;
                    }

                }

                i = l;
                j = l + l;

                while (j <= ir)
                {
                    if (j < ir)
                    {
                        if (a[indx[j] - 1] < a[indx[j - 1] - 1])
                        {
                            j = j + 1;
                        }
                    }

                    if (a[indx[j - 1] - 1] < aval)
                    {
                        indx[i - 1] = indx[j - 1];
                        i = j;
                        j = j + j;
                    }
                    else
                    {
                        j = ir + 1;
                    }

                }

                indx[i - 1] = indxt;

            }
        }

        public static void i4vec_sort_insert_a(int n, ref int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_SORT_INSERT_A uses an ascending insertion sort on an I4VEC.
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
            //    13 April 1999
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Donald Kreher, Douglas Simpson,
            //    Algorithm 1.1,
            //    Combinatorial Algorithms,
            //    CRC Press, 1998, page 11.
            //
            //  Parameters:
            //
            //    Input, int N, the number of items in the vector.
            //    N must be positive.
            //
            //    Input/output, int A[N].
            //    On input, A contains data to be sorted.
            //    On output, the entries of A have been sorted in ascending order.
            //
        {
            int i;
            int j;
            int x;

            for (i = 1; i < n; i++)
            {
                x = a[i];

                j = i;

                while (1 <= j && x < a[j - 1])
                {
                    a[j] = a[j - 1];
                    j = j - 1;
                }

                a[j] = x;
            }

            return;
        }

        public static void i4vec_sort_insert_d(int n, ref int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_SORT_INSERT_D uses a descending insertion sort on an I4VEC.
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
            //    13 April 1999
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Donald Kreher, Douglas Simpson,
            //    Algorithm 1.1,
            //    Combinatorial Algorithms,
            //    CRC Press, 1998, page 11.
            //
            //  Parameters:
            //
            //    Input, int N, the number of items in the vector.
            //    N must be positive.
            //
            //    Input/output, int A[N].
            //    On input, A contains data to be sorted.
            //    On output, the entries of A have been sorted in ascending order.
            //
        {
            int i;
            int j;
            int x;

            for (i = 1; i < n; i++)
            {
                x = a[i];
                j = i;

                while (1 <= j && a[j - 1] < x)
                {
                    a[j] = a[j - 1];
                    j = j - 1;
                }

                a[j] = x;
            }

            return;
        }
        
        public static void i4vec_indicator1 ( int n, ref int[] a )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_INDICATOR1 sets an I4VEC to the indicator1 vector.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    25 February 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of elements of A.
            //
            //    Output, int A[N], the initialized array.
            //
        {
            int i;

            for ( i = 0; i < n; i++ )
            {
                a[i] = i + 1;
            }
        }
        
        public static void i4vec_decrement ( int n, ref int[] v )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_DECREMENT decrements an I4VEC.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 January 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the size of the array.
            //
            //    Input/output, int V[N], the array to be decremented.
            //
        {
            int i;

            for ( i = 0; i < n; i++ )
            {
                v[i] = v[i] - 1;
            }
        }
        
        public static bool i4vec_descends ( int n, ref int[] x )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_DESCENDS determines if an I4VEC is decreasing.
            //
            //  Example:
            //
            //    X = ( 9, 7, 7, 3, 2, 1, -8 )
            //
            //    I4VEC_DESCENDS = TRUE
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
            //    Input, int N, the size of the array.
            //
            //    Input, int X[N], the array to be examined.
            //
            //    Output, bool I4VEC_DESCEND, is TRUE if the entries of the array descend.
            //
        {
            int i;

            for ( i = 0; i < n - 1; i++ )
            {
                if ( x[i] < x[i+1] )
                {
                    return false;
                }
            }
            return true;
        }
        
        public static int i4vec_dot_product ( int n, int[] x, int[] y )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_DOT_PRODUCT computes the dot product of two I4VEC's.
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
        //    19 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the size of the array.
        //
        //    Input, int X[N], Y[N], the arrays.
        //
        //    Output, int I4VEC_DOT_PRODUCT, the dot product of X and Y.
        //
        {
            int i;
            int value;

            value = 0;
            for ( i = 0; i < n; i++ )
            {
                value = value + x[i] * y[i];
            }

            return value;
        }

    }
}