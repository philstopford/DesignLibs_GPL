using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void r8vec_permute(int n, int[] p, ref double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_PERMUTE permutes an R8VEC in place.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //    This routine permutes an array of real "objects", but the same
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
            //      A = ( 1.0, 2.0, 3.0, 4.0, 5.0 )
            //
            //    Output:
            //
            //      A    = ( 2.0, 4.0, 5.0, 1.0, 3.0 ).
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
            //    Input, int P[N], the permutation.
            //
            //    Input/output, double A[N], the array to be permuted.
            //
        {
            double a_temp;
            int i;
            int iget;
            int iput;
            int istart;

            perm_check0(n, p);
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
                            Console.WriteLine("R8VEC_PERMUTE - Fatal error!");
                            Console.WriteLine("  A permutation index is out of range.");
                            Console.WriteLine("  P(" + iput + ") = " + iget + "");
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

        }

        public static void r8vec_permute_cyclic(int n, int k, ref double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_PERMUTE_CYCLIC performs a cyclic permutation of an R8VEC.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //    For 0 <= K < N, this function cyclically permutes the input vector
            //    to have the form
            //
            //     ( A[K], A[K+1], ..., A[N-1], A[0], ..., A[K-1] )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 August 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of objects.
            //
            //    Input, int K, the increment used.
            //
            //    Input/output, double A[N], the array to be permuted.
            //
        {
            double[] b;
            int i;
            int ipk;

            b = new double[n];

            for (i = 0; i < n; i++)
            {
                ipk = i4_wrap(i + k, 0, n - 1);
                b[i] = a[ipk];
            }

            for (i = 0; i < n; i++)
            {
                a[i] = b[i];
            }
        }

        public static void r8vec_permute_uniform(int n, ref double[] a, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_PERMUTE_UNIFORM randomly permutes an R8VEC.
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
            //    Input/output, double A[N], the array to be permuted.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
        {
            int[] p;

            p = perm0_uniform_new(n, ref seed);

            r8vec_permute(n, p, ref a);

        }


    }
}