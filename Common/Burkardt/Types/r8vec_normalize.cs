using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void r8vec_normalize(int n, ref double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_NORMALIZE normalizes an R8VEC.
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
            //    11 October 2010
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
            //    On output, A should have unit Euclidean norm.
            //
        {
            int i;
            double norm;

            norm = 0.0;
            for (i = 0; i < n; i++)
            {
                norm = norm + a[i] * a[i];
            }

            norm = Math.Sqrt(norm);

            if (norm == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8VEC_NORMALIZE - Fatal error!");
                Console.WriteLine("  The vector norm is 0.");
                return ;
            }

            for (i = 0; i < n; i++)
            {
                a[i] = a[i] / norm;
            }

            return;
        }

        public static void r8vec_normalize_l1(int n, ref double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_NORMALIZE_L1 normalizes an R8VEC to have unit sum.
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
            double a_sum;
            int i;

            a_sum = 0.0;
            for (i = 0; i < n; i++)
            {
                a_sum = a_sum + a[i];
            }

            if (a_sum == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8VEC_NORMALIZE_L1 - Fatal error!");
                Console.WriteLine("  The vector entries sum to 0.");
                return;
            }

            for (i = 0; i < n; i++)
            {
                a[i] = a[i] / a_sum;
            }

            return;
        }
    }
}