using System;
using System.Linq;
using Burkardt.Uniform;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {

        public static double r8vec_sum(int n, double[] a)
//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SUM returns the sum of an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A[N], the vector.
//
//    Output, double R8VEC_SUM, the sum of the vector.
//
        {
            double value = 0.0;
            for (int i = 0; i < n; i++)
            {
                value = value + a[i];
            }

            return value;
        }


        public static int[] r8vec_sort_heap_index_a_new(int n, double[] a)

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SORT_HEAP_INDEX_A_NEW does an indexed heap ascending sort of an R8VEC.
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
//    30 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double A[N], an array to be index-sorted.
//
//    Output, int R8VEC_SORT_HEAP_INDEX_A_NEW[N], contains the sort index.  The
//    I-th element of the sorted array is A(INDX(I)).
//
        {
            int[] indx = new int[n];

            for (int i = 1; i <= n; i++)
            {
                indx[i - 1] = i;
            }

            int l = n / 2 + 1;
            int ir = n;

            for (;;)
            {
                double aval;
                int indxt;
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
                        for (int i = 0; i < n; i++)
                        {
                            indx[i] = indx[i] - 1;
                        }

                        break;
                    }
                }

                int i2 = l;
                int j = l + l;

                while (j <= ir)
                {
                    if (j < ir)
                    {
                        if (a[indx[j - 1] - 1] < a[indx[j] - 1])
                        {
                            j = j + 1;
                        }
                    }

                    if (aval < a[indx[j - 1] - 1])
                    {
                        indx[i2 - 1] = indx[j - 1];
                        i2 = j;
                        j = j + j;
                    }
                    else
                    {
                        j = ir + 1;
                    }
                }

                indx[i2 - 1] = indxt;
            }

            return indx;
        }

        public static double r8vec_dot(int n, double[] a1, double[] a2)
//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DOT computes the dot product of a pair of R8VEC's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], A2[N], the two vectors to be considered.
//
//    Output, double R8VEC_DOT, the dot product of the vectors.
//
        {
            double value = 0.0;
            for (int i = 0; i < n; i++)
            {
                value = value + a1[i] * a2[i];
            }

            return value;
        }

        public static double r8vec_dot_product(int n, double[] a1, double[] a2)
//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DOT_PRODUCT computes the dot product of a pair of R8VEC's.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], A2[N], the two vectors to be considered.
//
//    Output, double R8VEC_DOT_PRODUCT, the dot product of the vectors.
//
        {
            double value = 0.0;
            for (int i = 0; i < n; i++)
            {
                value = value + a1[i] * a2[i];
            }

            return value;
        }

        public static void r8vec_print(int n, double[] a, string title)
//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT prints an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
        {

            Console.WriteLine(title);
            Console.WriteLine();
            for (int i = 0; i < n; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(8)
                                       + ": " + a[i].ToString().PadLeft(14) + "");
            }
        }

        public static void r8vec_unit_sum(int n, ref double[] a)
//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIT_SUM normalizes an R8VEC to have unit sum.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input/output, double A[N], the vector to be normalized.
//    On output, the entries of A should have unit sum.  However, if
//    the input vector has zero sum, the routine halts.
//
        {
            double a_sum = 0.0;
            for (int i = 0; i < n; i++)
            {
                a_sum = a_sum + a[i];
            }

            if (a_sum == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8VEC_UNIT_SUM - Fatal error!");
                Console.WriteLine("  The vector entries sum to 0.");
            }

            for (int i = 0; i < n; i++)
            {
                a[i] = a[i] / a_sum;
            }
        }


        public static double r8vec_length(int dim_num, double[] x)
//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_LENGTH returns the Euclidean length of an R8VEC
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double X[DIM_NUM], the vector.
//
//    Output, double R8VEC_LENGTH, the Euclidean length of the vector.
//
        {
            double value = 0.0;
            for (int i = 0; i < dim_num; i++)
            {
                value = value + Math.Pow(x[i], 2);
            }

            value = Math.Sqrt(value);

            return value;
        }

        public static double r8vec_circular_variance(int n, double[] x)
//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_CIRCULAR_VARIANCE returns the circular variance of an R8VEC
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 December 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double X(N), the vector whose variance is desired.
//
//    Output, double R8VEC_CIRCULAR VARIANCE, the circular variance
//    of the vector entries.
//
        {
            double mean = r8vec_mean(n, x);

            double sum_c = 0.0;
            for (int i = 0; i < n; i++)
            {
                sum_c = sum_c + Math.Cos(x[i] - mean);
            }

            double sum_s = 0.0;
            for (int i = 0; i < n; i++)
            {
                sum_s = sum_s + Math.Sin(x[i] - mean);
            }

            double value = Math.Sqrt(sum_c * sum_c + sum_s * sum_s) / (double) n;

            value = 1.0 - value;

            return value;
        }

        public static double r8vec_diff_norm(int n, double[] a, double[] b)
//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DIFF_NORM returns the L2 norm of the difference of R8VEC's.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    The vector L2 norm is defined as:
//
//      R8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 June 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, double A[N], B[N], the vectors.
//
//    Output, double R8VEC_DIFF_NORM, the L2 norm of A - B.
//
        {
            double value = 0.0;

            for (int i = 0; i < n; i++)
            {
                value = value + (a[i] - b[i]) * (a[i] - b[i]);
            }

            value = Math.Sqrt(value);

            return value;
        }

        public static double r8vec_max(int n, double[] dvec)
//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MAX returns the value of the maximum element in an R8VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double *RVEC, a pointer to the first entry of the array.
//
//    Output, double R8VEC_MAX, the value of the maximum element.  This
//    is set to 0.0 if N <= 0.
//
        {
            if (dvec.Length <= 0)
            {
                return 0;
            }

// Limit to the number of items in the array as a maximum
            n = Math.Min(n, dvec.Length);

            if (n == dvec.Length)
            {
                return dvec.Max();
            }

            return dvec.Take(n).Max();
        }

        public static double r8vec_mean(int n, double[] x)
//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MEAN returns the mean of an R8VEC.
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
//    Input, double X[N], the vector whose mean is desired.
//
//    Output, double R8VEC_MEAN, the mean, or average, of the vector entries.
//
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

        public static double r8vec_min(int n, double[] dvec)
//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MIN returns the value of the minimum element in an R8VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double *RVEC, a pointer to the first entry of the array.
//
//    Output, double R8VEC_MIN, the value of the minimum element.  This
//    is set to 0.0 if N <= 0.
//
        {
            if (dvec.Length <= 0)
            {
                return 0;
            }

// Limit to the number of items in the array as a maximum
            n = Math.Min(n, dvec.Length);

            if (n == dvec.Length)
            {
                return dvec.Min();
            }

            return dvec.Take(n).Min();
        }

        public static double r8vec_variance(int n, double[] x)
//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_VARIANCE returns the variance of an R8VEC.
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
//    Input, double X[N], the vector whose variance is desired.
//
//    Output, double R8VEC_VARIANCE, the variance of the vector entries.
//
        {
            double mean = r8vec_mean(n, x);

            double variance = 0.0;
            for (int i = 0; i < n; i++)
            {
                variance = variance + (x[i] - mean) * (x[i] - mean);
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


        public static double[] r8vec_zero_new(int n)
//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_ZERO_NEW creates and zeroes an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Output, double R8VEC_ZERO_NEW[N], a vector of zeroes.
//
        {
            double[] a = new double[n];

            for (int i = 0; i < n; i++)
            {
                a[i] = 0.0;
            }

            return a;
        }


        public static double[] r8vec_linspace_new(int n, double a_first, double a_last)
//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_LINSPACE_NEW creates a vector of linearly spaced values.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
//
//    In other words, the interval is divided into N-1 even subintervals,
//    and the endpoints of intervals are used as the points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 March 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A_FIRST, A_LAST, the first and last entries.
//
//    Output, double R8VEC_LINSPACE_NEW[N], a vector of linearly spaced data.
//
        {
            double[] a = new double[n];

            if (n == 1)
            {
                a[0] = (a_first + a_last) / 2.0;
            }
            else
            {
                for (int i = 0; i < n; i++)
                {
                    a[i] = ((double) (n - 1 - i) * a_first
                            + (double) (i) * a_last)
                           / (double) (n - 1);
                }
            }

            return a;
        }

        public static double[] r8vec_even_new(int n, double alo, double ahi)
//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_EVEN_NEW returns an R8VEC of values evenly spaced between ALO and AHI.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of values.
//
//    Input, double ALO, AHI, the low and high values.
//
//    Output, double R8VEC_EVEN_NEW[N], N evenly spaced values.
//    Normally, A[0] = ALO and A[N-1] = AHI.
//    However, if N = 1, then A[0] = 0.5*(ALO+AHI).
//
        {
            double[] a = new double[n];

            if (n == 1)
            {
                a[0] = 0.5 * (alo + ahi);
            }
            else
            {
                for (int i = 0; i < n; i++)
                {
                    a[i] = ((double) (n - i - 1) * alo
                            + (double) (i) * ahi)
                           / (double) (n - 1);
                }
            }

            return a;
        }

        public static void r8vec_bracket3(int n, double[] t, double tval, ref int left)
//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_BRACKET3 finds the interval containing or nearest a given value.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    The routine always returns the index LEFT of the sorted array
//    T with the property that either
//    *  T is contained in the interval [ T[LEFT], T[LEFT+1] ], or
//    *  T < T[LEFT] = T[0], or
//    *  T > T[LEFT+1] = T[N-1].
//
//    The routine is useful for interpolation problems, where
//    the abscissa must be located within an interval of data
//    abscissas for interpolation, or the "nearest" interval
//    to the (extreme) abscissa must be found so that extrapolation
//    can be carried out.
//
//    This version of the function has been revised so that the value of
//    LEFT that is returned uses the 0-based indexing natural to C++.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, length of the input array.
//
//    Input, double T[N], an array that has been sorted into ascending order.
//
//    Input, double TVAL, a value to be bracketed by entries of T.
//
//    Input/output, int *LEFT.
//    On input, if 0 <= LEFT <= N-2, LEFT is taken as a suggestion for the
//    interval [ T[LEFT-1] T[LEFT] ] in which TVAL lies.  This interval
//    is searched first, followed by the appropriate interval to the left
//    or right.  After that, a binary search is used.
//    On output, LEFT is set so that the interval [ T[LEFT], T[LEFT+1] ]
//    is the closest to TVAL; it either contains TVAL, or else TVAL
//    lies outside the interval [ T[0], T[N-1] ].
//
        {
            int high;
            int low;
            int mid;
//  
//  Check the input data.
//
            if (n < 2)
            {
                Console.WriteLine("");
                Console.WriteLine("R8VEC_BRACKET3 - Fatal error//");
                Console.WriteLine("  N must be at least 2.");
                return;
            }

//
//  If *LEFT is not between 0 and N-2, set it to the middle value.
//
            if (left < 0 || n - 2 < left)
            {
                left = (n - 1) / 2;
            }

//
//  CASE 1: TVAL < T[*LEFT]:
//  Search for TVAL in (T[I],T[I+1]), for I = 0 to *LEFT-1.
//
            if (tval < t[left])
            {
                if (left == 0)
                {
                    return;
                }
                else if (left == 1)
                {
                    left = 0;
                    return;
                }
                else if (t[left - 1] <= tval)
                {
                    left = left - 1;
                    return;
                }
                else if (tval <= t[1])
                {
                    left = 0;
                    return;
                }

// 
//  ...Binary search for TVAL in (T[I],T[I+1]), for I = 1 to *LEFT-2.
//
                low = 1;
                high = left - 2;

                for (;;)
                {
                    if (low == high)
                    {
                        left = low;
                        return;
                    }

                    mid = (low + high + 1) / 2;

                    if (t[mid] <= tval)
                    {
                        low = mid;
                    }
                    else
                    {
                        high = mid - 1;
                    }
                }
            }
// 
//  CASE 2: T[*LEFT+1] < TVAL:
//  Search for TVAL in (T[I],T[I+1]) for intervals I = *LEFT+1 to N-2.
//
            else if (t[left + 1] < tval)
            {
                if (left == n - 2)
                {
                    return;
                }
                else if (left == n - 3)
                {
                    left = left + 1;
                    return;
                }
                else if (tval <= t[left + 2])
                {
                    left = left + 1;
                    return;
                }
                else if (t[n - 2] <= tval)
                {
                    left = n - 2;
                    return;
                }

// 
//  ...Binary search for TVAL in (T[I],T[I+1]) for intervals I = *LEFT+2 to N-3.
//
                low = left + 2;
                high = n - 3;

                for (;;)
                {

                    if (low == high)
                    {
                        left = low;
                        return;
                    }

                    mid = (low + high + 1) / 2;

                    if (t[mid] <= tval)
                    {
                        low = mid;
                    }
                    else
                    {
                        high = mid - 1;
                    }
                }
            }
//
//  CASE 3: T[*LEFT] <= TVAL <= T[*LEFT+1]:
//  T is just where the user said it might be.
//
            else
            {
            }

            return;
        }


        public static void r8vec_bracket4(int nt, double[] t, int ns, double[] s, ref int[] left )
//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_BRACKET4 finds the interval containing or nearest a given value.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    The routine always returns the index LEFT of the sorted array
//    T with the property that either
//    *  T is contained in the interval [ T[LEFT], T[LEFT+1] ], or
//    *  T < T[LEFT] = T[0], or
//    *  T > T[LEFT+1] = T[NT-1].
//
//    The routine is useful for interpolation problems, where
//    the abscissa must be located within an interval of data
//    abscissas for interpolation, or the "nearest" interval
//    to the (extreme) abscissa must be found so that extrapolation
//    can be carried out.
//
//    This version of the function has been revised so that the value of
//    LEFT that is returned uses the 0-based indexing natural to C++.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NT, length of the input array.
//
//    Input, double T[NT], an array that has been sorted
//    into ascending order.
//
//    Input, int NS, the number of points to be bracketed.
//
//    Input, double S[NS], values to be bracketed by entries of T.
//
//    Output, int LEFT[NS].
//    LEFT[I] is set so that the interval [ T[LEFT[I]], T[LEFT[I]+1] ]
//    is the closest to S[I]; it either contains S[I], or else S[I]
//    lies outside the interval [ T[0], T[NT-1] ].
//
        {
            int high;
            int low;
            int mid;
//  
//  Check the input data.
//
            if (nt < 2)
            {
                Console.WriteLine("");
                Console.WriteLine("R8VEC_BRACKET4 - Fatal error!");
                Console.WriteLine("  NT must be at least 2.");
                return;
            }

            for (int i = 0; i < ns; i++)
            {
                left[i] = (nt - 1) / 2;
//
//  CASE 1: S[I] < T[LEFT]:
//  Search for S[I] in (T[I],T[I+1]), for I = 0 to LEFT-1.
//
                if (s[i] < t[left[i]])
                {
                    if (left[i] == 0)
                    {
                        continue;
                    }
                    else if (left[i] == 1)
                    {
                        left[i] = 0;
                        continue;
                    }
                    else if (t[left[i] - 1] <= s[i])
                    {
                        left[i] = left[i] - 1;
                        continue;
                    }
                    else if (s[i] <= t[1])
                    {
                        left[i] = 0;
                        continue;
                    }

// 
//  ...Binary search for S[I] in (T[I],T[I+1]), for I = 1 to *LEFT-2.
//
                    low = 1;
                    high = left[i] - 2;

                    for (;;)
                    {
                        if (low == high)
                        {
                            left[i] = low;
                            break;
                        }

                        mid = (low + high + 1) / 2;

                        if (t[mid] <= s[i])
                        {
                            low = mid;
                        }
                        else
                        {
                            high = mid - 1;
                        }
                    }
                }
// 
//  CASE 2: T[LEFT+1] < S[I]:
//  Search for S[I] in (T[I],T[I+1]) for intervals I = LEFT+1 to NT-2.
//
                else if (t[left[i] + 1] < s[i])
                {
                    if (left[i] == nt - 2)
                    {
                        continue;
                    }
                    else if (left[i] == nt - 3)
                    {
                        left[i] = left[i] + 1;
                        continue;
                    }
                    else if (s[i] <= t[left[i] + 2])
                    {
                        left[i] = left[i] + 1;
                        continue;
                    }
                    else if (t[nt - 2] <= s[i])
                    {
                        left[i] = nt - 2;
                        continue;
                    }

// 
//  ...Binary search for S[I] in (T[I],T[I+1]) for intervals I = LEFT+2 to NT-3.
//
                    low = left[i] + 2;
                    high = nt - 3;

                    for (;;)
                    {

                        if (low == high)
                        {
                            left[i] = low;
                            break;
                        }

                        mid = (low + high + 1) / 2;

                        if (t[mid] <= s[i])
                        {
                            low = mid;
                        }
                        else
                        {
                            high = mid - 1;
                        }
                    }
                }
//
//  CASE 3: T[LEFT] <= S[I] <= T[LEFT+1]:
//
                else
                {
                }
            }

            return;
        }

        public static double[] r8vec_normal_01_new(int n, ref int seed)
//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_NORMAL_01_NEW returns a unit pseudonormal R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    The standard normal probability distribution function (PDF) has
//    mean 0 and standard deviation 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of values desired.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double R8VEC_NORMAL_01_NEW[N], a sample of the standard normal PDF.
//
//  Local parameters:
//
//    Local, double R[N+1], is used to store some uniform random values.
//    Its dimension is N+1, but really it is only needed to be the
//    smallest even number greater than or equal to N.
//
//    Local, int X_LO, X_HI, records the range of entries of
//    X that we need to compute.
//
        {
            int i;
            int m;
            double[] r;
            double[] x;
            int x_hi;
            int x_lo;

            x = new double[n];
//
//  Record the range of X we need to fill in.
//
            x_lo = 1;
            x_hi = n;
//
//  If we need just one new value, do that here to avoid null arrays.
//
            if (x_hi - x_lo + 1 == 1)
            {
                r = UniformRNG.r8vec_uniform_01_new(2, ref seed);

                x[x_hi - 1] = Math.Sqrt(-2.0 * Math.Log(r[0])) * Math.Cos(2.0 * Math.PI * r[1]);

            }
//
//  If we require an even number of values, that's easy.
//
            else if ((x_hi - x_lo + 1) % 2 == 0)
            {
                m = (x_hi - x_lo + 1) / 2;

                r = UniformRNG.r8vec_uniform_01_new(2 * m, ref seed);

                for (i = 0; i <= 2 * m - 2; i = i + 2)
                {
                    x[x_lo + i - 1] = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Cos(2.0 * Math.PI * r[i + 1]);
                    x[x_lo + i] = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Sin(2.0 * Math.PI * r[i + 1]);
                }

            }
//
//  If we require an odd number of values, we generate an even number,
//  and handle the last pair specially, storing one in X(N), and
//  saving the other for later.
//
            else
            {
                x_hi = x_hi - 1;

                m = (x_hi - x_lo + 1) / 2 + 1;

                r = UniformRNG.r8vec_uniform_01_new(2 * m, ref seed);

                for (i = 0; i <= 2 * m - 4; i = i + 2)
                {
                    x[x_lo + i - 1] = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Cos(2.0 * Math.PI * r[i + 1]);
                    x[x_lo + i] = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Sin(2.0 * Math.PI * r[i + 1]);
                }

                i = 2 * m - 2;

                x[x_lo + i - 1] = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Cos(2.0 * Math.PI * r[i + 1]);

            }

            return x;
        }

        public static void r8vec_mesh_2d ( int nx, int ny, double[] xvec, double[] yvec, 
        double[] xmat, double[] ymat )
//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MESH_2D creates a 2D mesh from X and Y vectors.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    NX = 2
//    XVEC = ( 1, 2, 3 )
//    NY = 3
//    YVEC = ( 4, 5 )
//
//    XMAT = (
//      1, 2, 3
//      1, 2, 3 )
//
//    YMAT = (
//      4, 4, 4
//      5, 5, 5 ) 
// 
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 July 2013
//
//  Parameters:
//
//    Input, int NX, NY, the number of X and Y values.
//
//    Input, double XVEC[NX], YVEC[NY], the X and Y coordinate
//    values.
//
//    Output, double XMAT[NX*NY], YMAT[NX*NY], the coordinate
//    values of points on an NX by NY mesh.
//
        {
            for (int j = 0; j < ny; j++ )
            {
                for (int i = 0; i < nx; i++ )
                {
                    xmat[i+j*nx] = xvec[i];
                }
            }

            for (int j = 0; j < ny; j++ )
            {
                for (int i = 0; i < nx; i++ )
                {
                    ymat[i+j*nx] = yvec[j];
                }
            }
        }

        public static void r8vec2_print(int n, double[] a1, double[] a2, string title)
//****************************************************************************80
//
//  Purpose:
//
//    R8VEC2_PRINT prints an R8VEC2.
//
//  Discussion:
//
//    An R8VEC2 is a dataset consisting of N pairs of real values, stored
//    as two separate vectors A1 and A2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 November 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A1[N], double A2[N], the vectors to be printed.
//
//    Input, string TITLE, a title.
//
        {
            Console.WriteLine();
            Console.WriteLine(title);
            Console.WriteLine();
            for (int i = 0; i <= n - 1; i++)
            {
                Console.WriteLine(i.ToString().PadLeft(6)
                                  + ": " + a1[i].ToString().PadLeft(14)
                                  + "  " + a2[i].ToString().PadLeft(14));
            }
        }
        
        public static double r8vec_norm_affine ( int n, double[] v0, double[] v1 )
//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_NORM_AFFINE returns the affine L2 norm of an R8VEC.
//
//  Discussion:
//
//    The affine vector L2 norm is defined as:
//
//      R8VEC_NORM_AFFINE(V0,V1)
//        = sqrt ( sum ( 1 <= I <= N ) ( V1(I) - V0(I) )^2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the vectors.
//
//    Input, double V0[N], the base vector.
//
//    Input, double V1[N], the vector.
//
//    Output, double R8VEC_NORM_AFFINE, the affine L2 norm.
//
        {
            double value = 0.0;

            for (int i = 0; i < n; i++ )
            {
                value = value + ( v1[i] - v0[i] ) * ( v1[i] - v0[i] );
            }
            value = Math.Sqrt ( value );

            return value;
        }
        
        public static void r8vec_copy ( int n, double[] a1, ref double[] a2 )
//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_COPY copies an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], the vector to be copied.
//
//    Output, double A2[N], the copy of A1.
//
        {
            int i;

            for ( i = 0; i < n; i++ )
            {
                a2[i] = a1[i];
            }
        }
        
        public static void r8vec3_print ( int n, double[] a1, double[] a2, double[] a3, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC3_PRINT prints a triple of real vectors.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 September 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A1[N], double A2[N], double A3[N], the vectors
//    to be printed.
//
//    Input, string TITLE, a title.
//
        {
            Console.WriteLine();
            Console.WriteLine(title);
            Console.WriteLine();
            for (int i = 0; i <= n - 1; i++ )
            {
                Console.WriteLine( i.ToString().PadLeft(4)     + ": "
                                  + a1[i].ToString().PadLeft(10) + "  "
                                  + a2[i].ToString().PadLeft(10) + "  "
                                  + a3[i].ToString().PadLeft(10) );
            }
        }
        
        public static double[]  r8vec_copy_new ( int n, double[] a1 )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_COPY_NEW copies an R8VEC to a new R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], the vector to be copied.
//
//    Output, double R8VEC_COPY_NEW[N], the copy of A1.
//
        {
            int i;

            double[] a2 = new double[n];

            for ( i = 0; i < n; i++ )
            {
                a2[i] = a1[i];
            }
            return a2;
        }
        
        public static double r8vec_i4vec_dot_product ( int n, double[] r8vec, int[] i4vec )
//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_I4VEC_DOT_PRODUCT computes the dot product of an R8VEC and an I4VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double R8VEC[N], the first vector.
//
//    Input, int I4VEC[N], the second vector.
//
//    Output, double R8VEC_I4VEC_DOT_PRODUCT, the dot product of the vectors.
//
        {
            int i;

            double value = 0.0;
            for ( i = 0; i < n; i++ )
            {
                value = value + r8vec[i] * ( double ) ( i4vec[i] );
            }
            return value;
        }
        
        public static double r8vec_product ( int n, double[] a )
//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRODUCT returns the product of the entries of an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A[N], the vector.
//
//    Output, double R8VEC_PRODUCT, the product of the vector.
//
        {
            int i;
            double product;

            product = 1.0;
            for ( i = 0; i < n; i++ )
            {
                product = product * a[i];
            }

            return product;
        }
        
        public static double[] r8vec_even ( int n, double alo, double ahi )
//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_EVEN returns N real values, evenly spaced between ALO and AHI.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of values.
//
//    Input, double ALO, AHI, the low and high values.
//
//    Output, double R8VEC_EVEN[N], N evenly spaced values.
//    Normally, A(1) = ALO and A(N) = AHI.
//    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
//
        {
            double[] a = new double[n];

            if ( n == 1 )
            {
                a[0] = 0.5 * ( alo + ahi );
            }
            else
            {
                for (int i = 1; i <= n; i++ )
                {
                    a[i-1] = ( ( double ) ( n - i     ) * alo
                               + ( double ) (     i - 1 ) * ahi )
                             / ( double ) ( n     - 1 );
                }
            }

            return a;
        }
        
        
        public static double[] r8vec_midspace_new ( int n, double a, double b )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_MIDSPACE_NEW creates a vector of linearly spaced values.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //    This function divides the interval [a,b] into n subintervals, and then
            //    returns the midpoints of those subintervals.
            //
            //  Example:
            //
            //    N = 5, A = 10, B = 20
            //    X = [ 11, 13, 15, 17, 19 ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 June 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vector.
            //
            //    Input, double A, B, the endpoints of the interval.
            //
            //    Output, double R8VEC_MIDSPACE_NEW[N], a vector of linearly spaced data.
            //
        {
            double[] x = new double[n];

            for (int i = 0; i < n; i++ )
            {
                x[i] = ( ( double ) ( 2 * n - 2 * i - 1 ) * a 
                         + ( double ) (         2 * i + 1 ) * b ) 
                       / ( double ) ( 2 * n );
            }

            return x;
        }

        public static double[] r8vec_cheby1space_new(int n, double a, double b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_CHEBY1SPACE_NEW creates Type 1 Chebyshev values in [A,B].
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 September 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vector.
            //
            //    Input, double A, B, the interval.
            //
            //    Output, double R8VEC_CHEBY1SPACE_NEW[N], a vector of Type 1
            //    Chebyshev spaced data.
            //
        {
            double c;
            double theta;

            double[] x = new double[n];

            if (n == 1)
            {
                x[0] = (a + b) / 2.0;
            }
            else
            {
                for (int i = 0; i < n; i++)
                {
                    theta = (double) (n - i - 1) * Math.PI / (double) (n - 1);

                    c = Math.Cos(theta);

                    if ((n % 2) == 1)
                    {
                        if (2 * i + 1 == n)
                        {
                            c = 0.0;
                        }
                    }

                    x[i] = ((1.0 - c) * a
                            + (1.0 + c) * b)
                           / 2.0;
                }
            }

            return x;
        }

        public static double[] r8vec_cheby2space_new ( int n, double a, double b )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_CHEBY2SPACE_NEW creates Type 2 Chebyshev values in [A,B].
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 July 2017
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vector.
            //
            //    Input, double A, B, the interval.
            //
            //    Output, double R8VEC_CHEBY2SPACE_NEW[N], a vector of Type 2
            //    Chebyshev spaced data.
            //
        {
            double c;
            int i;
            double theta;

            double[] x = new double[n];

            for ( i = 0; i < n; i++ )
            {
                theta = ( double ) ( n - i ) * Math.PI / ( double ) ( n + 1 );

                c = Math.Cos ( theta );

                x[i] = ( ( 1.0 - c ) * a  
                         + ( 1.0 + c ) * b ) 
                       /   2.0;
            }

            return x;
        }

    }
}