using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
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

        public static double r8vec_asum ( int n, double[] a, int aIndex = 0 )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_ASUM sums the absolute values of the entries of an R8VEC.
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
            //    24 January 2015
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
            //    Output, double R8VEC_ASUM, the sum of absolute values of the entries.
            //
        {
            int i;

            double value = 0.0;
            for ( i = 0; i < n; i++ )
            {
                value = value + Math.Abs ( a[aIndex + i] );
            }
            return value;
        }
        
        public static double r8vec_sum(int n, double[] a, int aIndex = 0)
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
                value = value + a[aIndex + i];
            }

            return value;
        }
        
        public static int r8vec2_sum_max_index ( int n, double[] a, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC2_SUM_MAX_INDEX returns the index of the maximum sum of two R8VEC's.
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
        //    17 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the array.
        //
        //    Input, double A[N], B[N], two arrays whose sum
        //    is to be examined.
        //
        //    Output, int R8VEC2_SUM_MAX_INDEX, the index of the largest entry in A+B.
        //
        {
            int i;
            double sum_max;
            int sum_max_index;

            if ( n <= 0 )
            {
                sum_max_index = -1;
            }
            else
            {
                sum_max_index = 1;
                sum_max = a[0] + b[0];

                for ( i = 2; i <= n; i++ )
                {
                    if ( sum_max < a[i-1] + b[i-1] )
                    {
                        sum_max = a[i-1] + b[i-1];
                        sum_max_index = i;
                    }
                }
            }
            return sum_max_index;
        }



    }
}