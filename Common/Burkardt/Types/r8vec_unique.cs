using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double[] r8vec_uniform_unit_new(int m, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_UNIFORM_UNIT_NEW generates a random unit vector.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 October 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the dimension of the space.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, double R8VEC_UNIFORM_UNIT_NEW[M], a random direction vector, 
            //    with unit norm.
            //
        {
            double[] a;
            int i;
            double norm;
            //
            //  Take M random samples from the normal distribution.
            //
            a = r8vec_normal_01_new(m, ref seed);
            //
            //  Compute the norm.
            //
            norm = 0.0;
            for (i = 0; i < m; i++)
            {
                norm = norm + a[i] * a[i];
            }

            norm = Math.Sqrt(norm);
            //
            //  Normalize.
            //
            for (i = 0; i < m; i++)
            {
                a[i] = a[i] / norm;
            }

            return a;
        }

        public static int r8vec_unique_count(int n, double[] a, double tol)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_UNIQUE_COUNT counts the unique elements in an unsorted R8VEC.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
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
            //    Input, double A[N], the array to examine, which does NOT have to
            //    be sorted.
            //
            //    Input, double TOL, a tolerance for checking equality.
            //
            //    Output, int R8VEC_UNIQUE_COUNT, the number of unique elements of A.
            //
        {
            int i;
            int j;
            int unique_num;

            unique_num = 0;

            for (i = 0; i < n; i++)
            {
                unique_num = unique_num + 1;

                for (j = 0; j < i; j++)
                {
                    if (Math.Abs(a[i] - a[j]) <= tol)
                    {
                        unique_num = unique_num - 1;
                        break;
                    }
                }
            }

            return unique_num;
        }

        public static int[] r8vec_unique_index(int n, double[] a, double tol)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_UNIQUE_INDEX indexes the unique occurrence of values in an R8VEC.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //    For element A(I) of the vector, UNIQUE_INDEX(I) is the uniqueness index
            //    of A(I).  That is, if A_UNIQUE contains the unique elements of A,
            //    gathered in order, then
            //
            //      A_UNIQUE ( UNIQUE_INDEX(I) ) = A(I)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 August 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of elements of A.
            //
            //    Input, double A[N], the unsorted array to examine.
            //
            //    Input, double TOL, a tolerance for equality.
            //
            //    Output, int R8VEC_UNIQUE_INDEX[N], the unique index.
            //
        {
            int i;
            int j;
            int[] unique_index;
            int unique_num;

            unique_index = new int[n];

            for (i = 0; i < n; i++)
            {
                unique_index[i] = -1;
            }

            unique_num = 0;

            for (i = 0; i < n; i++)
            {
                if (unique_index[i] == -1)
                {
                    unique_index[i] = unique_num;
                    for (j = i + 1; j < n; j++)
                    {
                        if (Math.Abs(a[i] - a[j]) <= tol)
                        {
                            unique_index[j] = unique_num;
                        }
                    }

                    unique_num = unique_num + 1;
                }
            }

            return unique_index;
        }


    }
}