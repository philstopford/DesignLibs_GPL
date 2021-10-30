namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static int r8vec_compare(int n, double[] a, double[] b, int aIndex = 0, int bIndex = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_COMPARE compares two R8VEC's.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //    The lexicographic ordering is used.
            //
            //  Example:
            //
            //    Input:
            //
            //      A1 = ( 2.0, 6.0, 2.0 )
            //      A2 = ( 2.0, 8.0, 12.0 )
            //
            //    Output:
            //
            //      ISGN = -1
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 September 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vectors.
            //
            //    Input, double A[N], B[N], the vectors to be compared.
            //
            //    Output, int R8VEC_COMPARE, the results of the comparison:
            //    -1, A is lexicographically less than B,
            //     0, A is equal to B,
            //    +1, A is lexicographically greater than B.
            //
        {
            int isgn;
            int k;

            isgn = 0;

            for (k = 0; k < n; k++)
            {
                if (a[(k + aIndex) % a.Length ] < b[(k + bIndex) % b.Length])
                {
                    isgn = -1;
                    return isgn;
                }
                else if (b[(k + bIndex) % b.Length] < a[(k + aIndex) % a.Length])
                {
                    isgn = +1;
                    return isgn;
                }
            }

            return isgn;
        }


        public static bool r8vec_eq(int n, double[] a1, double[] a2, int startIndexA1 = 0, int startIndexA2 = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_EQ is true if every pair of entries in two vectors is equal.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 August 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vectors.
            //
            //    Input, double A1[N], A2[N], two vectors to compare.
            //
            //    Output, bool R8VEC_EQ.
            //    R8VEC_EQ is TRUE if every pair of elements A1(I) and A2(I) are equal,
            //    and FALSE otherwise.
            //
        {
            int i;

            for (i = 0; i < n; i++)
            {
                if (a1[startIndexA1 + i] != a2[startIndexA2 + i])
                {
                    return false;
                }
            }

            return true;
        }

        public static bool r8vec_gt(int n, double[] a1, int startIndexA1, double[] a2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_GT == ( A1 > A2 ) for real vectors.
            //
            //  Discussion:
            //
            //    The comparison is lexicographic.
            //
            //    A1 > A2  <=>                              A1(1) > A2(1) or
            //                 ( A1(1)     == A2(1)     and A1(2) > A2(2) ) or
            //                 ...
            //                 ( A1(1:N-1) == A2(1:N-1) and A1(N) > A2(N)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 August 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the dimension of the vectors.
            //
            //    Input, double A1[N], A2[N], the vectors to be compared.
            //
            //    Output, bool R8VEC_GT, is TRUE if and only if A1 > A2.
            //
        {
            int i;

            for (i = 0; i < n; i++)
            {
                if (a2[i] < a1[startIndexA1 + i])
                {
                    return true;
                }
                else if (a1[startIndexA1 + i] < a2[i])
                {
                    return false;
                }
            }

            return false;
        }

        public static bool r8vec_lt(int n, double[] a1, double[] a2, int startIndexA1 = 0, int startIndexA2 = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_LT == ( A1 < A2 ) for real vectors.
            //
            //  Discussion:
            //
            //    The comparison is lexicographic.
            //
            //    A1 < A2  <=>                              A1(1) < A2(1) or
            //                 ( A1(1)     == A2(1)     and A1(2) < A2(2) ) or
            //                 ...
            //                 ( A1(1:N-1) == A2(1:N-1) and A1(N) < A2(N)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 August 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the dimension of the vectors.
            //
            //    Input, double A1[N], A2[N], the vectors to be compared.
            //
            //    Output, bool R8VEC_LT, is TRUE if and only if A1 < A2.
            //
        {
            int i;

            for (i = 0; i < n; i++)
            {
                if (a1[startIndexA1 + i] < a2[startIndexA2 + i])
                {
                    return true;
                }
                else if (a2[startIndexA2 + i] < a1[startIndexA1 + i])
                {
                    return false;
                }

            }

            return false;
        }

    }
}