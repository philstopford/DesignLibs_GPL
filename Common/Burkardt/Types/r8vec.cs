using System;
using System.Linq;

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

                public static double r8vec_dot ( int n, double[] a1, double[] a2 )
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
                    for (int i = 0; i < n; i++ )
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

        public static void r8vec_unit_sum ( int n, ref double[] a )
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
            for (int i = 0; i < n; i++ )
            {
                a_sum = a_sum + a[i];
            }

            if ( a_sum == 0.0 )
            {
                Console.WriteLine("");
                Console.WriteLine("R8VEC_UNIT_SUM - Fatal error!");
                Console.WriteLine("  The vector entries sum to 0.");
            }

            for (int i = 0; i < n; i++ )
            {
                a[i] = a[i] / a_sum;
            }
        }

                
                public static double r8vec_length ( int dim_num, double[] x )
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
                    for (int i = 0; i < dim_num; i++ )
                    {
                        value = value + Math.Pow ( x[i], 2 );
                    }
                    value = Math.Sqrt ( value );
        
                    return value;
                }

                public static double r8vec_circular_variance ( int n, double[] x )
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
            double mean = r8vec_mean ( n, x );

            double sum_c = 0.0;
            for (int i = 0; i < n; i++ )
            {
                sum_c = sum_c + Math.Cos ( x[i] - mean );
            }

            double sum_s = 0.0;
            for (int i = 0; i < n; i++ )
            {
                sum_s = sum_s + Math.Sin ( x[i] - mean );
            }

            double value = Math.Sqrt ( sum_c * sum_c + sum_s * sum_s ) / ( double ) n;

            value = 1.0 - value;

            return value;
        }
        
        public static double r8vec_diff_norm ( int n, double[] a, double[] b )
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

            for (int i = 0; i < n; i++ )
            {
                value = value + ( a[i] - b[i] ) * ( a[i] - b[i] );
            }
            value = Math.Sqrt ( value );

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

        public static double r8vec_mean ( int n, double[] x )
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

        public static double r8vec_min ( int n, double[] dvec )
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

        public static double r8vec_variance ( int n, double[] x )
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
            double mean = r8vec_mean ( n, x );

            double variance = 0.0;
            for (int i = 0; i < n; i++ )
            {
                variance = variance + ( x[i] - mean ) * ( x[i] - mean );
            }

            if ( 1 < n )
            {
                variance = variance / ( double ) ( n - 1 );
            }
            else
            {
                variance = 0.0;
            }

            return variance;
        }
        

        public static double[] r8vec_zero_new ( int n )
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

            for (int i = 0; i < n; i++ )
            {
                a[i] = 0.0;
            }
            return a;
        }
 
        
        public static double[] r8vec_linspace_new ( int n, double a_first, double a_last )
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

            if ( n == 1 )
            {
                a[0] = ( a_first + a_last ) / 2.0;
            }
            else
            {
                for (int i = 0; i < n; i++ )
                {
                    a[i] = ( ( double ) ( n - 1 - i ) * a_first 
                             + ( double ) (         i ) * a_last ) 
                           / ( double ) ( n - 1     );
                }
            }
            return a;
        }        
        
    }
}