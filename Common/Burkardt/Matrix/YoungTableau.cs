using System;
using Burkardt.Uniform;

namespace Burkardt.MatrixNS
{
    public static class YoungTableau
    {
        public static int ytb_enum(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    YTB_ENUM enumerates the Young tableau of size N.
            //
            //  Discussion:
            //
            //    If A(N) is the number of Young tableau of size N, then A(1) = 1,
            //    A(2) = 2, and
            //
            //    A(N) = A(N-1) + (N-1) * A(N-2).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    01 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the integer which is to be partitioned.
            //
            //    Output, int YTB_ENUM, the number of Young tableau of N.
            //
        {
            int a1;
            int a2;
            int a3;
            int i;
            int num;

            if (n <= 0)
            {
                num = 0;
            }
            else if (n == 1)
            {
                num = 1;
            }
            else if (n == 2)
            {
                num = 2;
            }
            else
            {
                a2 = 1;
                a3 = 2;
                for (i = 3; i <= n; i++)
                {
                    a1 = a2;
                    a2 = a3;
                    a3 = a2 + (i - 1) * a1;
                }

                num = a3;
            }

            return num;
        }

        public static void ytb_next(int n, ref int[] lambda, ref int[] a, ref bool more)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    YTB_NEXT computes the next Young tableau for a given shape.
            //
            //  Discussion:
            //
            //    When the routine is called with MORE = .FALSE. (the first time), and
            //    if LAMBDA on this call has M parts, with M<N, then the user
            //    must also make sure that LAMBDA(M+1) = 0.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 May 2003
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
            //    Input, int N, the integer which is to be partitioned.
            //
            //    Output, int LAMBDA[N], contains a partition of N, that is,
            //    the entries are positive integers that sum to N.
            //
            //    Output, int A[N].  A[I] is the row containing I
            //    in the output tableau.
            //
            //    Input/output, bool &MORE.  Set MORE FALSE before first call.
            //    It is reset to TRUE as the program returns a new tableau
            //    on each call, until the last tableau is computed, when
            //    the program also sets MORE = FALSE.
            //
        {
            int i;
            int ir;
            int it;
            int j;
            int k;
            int isave;

            it = n;

            if (more)
            {
                lambda[0] = 1;
                for (i = 1; i < n; i++)
                {
                    lambda[i] = 0;
                }

                isave = 0;

                for (i = 2; i <= n; i++)
                {
                    lambda[a[i - 1] - 1] = lambda[a[i - 1] - 1] + 1;

                    if (a[i - 1] < a[i - 2])
                    {
                        isave = i;
                        break;
                    }

                }

                if (isave == 0)
                {
                    more = false;
                    return;
                }

                it = lambda[a[isave - 1]];

                for (i = n; 1 <= i; i--)
                {
                    if (lambda[i - 1] == it)
                    {
                        a[isave - 1] = i;
                        lambda[i - 1] = lambda[i - 1] - 1;
                        it = isave - 1;
                        break;
                    }

                }

            }

            k = 1;
            ir = 1;

            for (;;)
            {
                if (n < ir)
                {
                    break;
                }

                if (lambda[ir - 1] != 0)
                {
                    a[k - 1] = ir;
                    lambda[ir - 1] = lambda[ir - 1] - 1;
                    k = k + 1;
                    ir = ir + 1;
                    continue;
                }

                if (it < k)
                {
                    break;
                }

                ir = 1;

            }

            if (n == 1)
            {
                more = false;
                return;
            }

            for (j = 1; j < n; j++)
            {
                if (a[j] < a[j - 1])
                {
                    more = true;
                    return;
                }
            }

            more = false;

            return;
        }

        public static void ytb_print(int n, int[] a, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    YTB_PRINT prints a Young tableau.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    01 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the integer that is partitioned.
            //
            //    Input, int A[N], describes the Young tableau.
            //    A[I] is the row of the tableau on which I occurs.
            //
            //    Input, string TITLE, a title.
            //
        {
            int j;
            int[] jarray;
            int row_i;
            int row_length;

            jarray = new int[n];

            if (0 < title.Length)
            {
                Console.WriteLine("");
                Console.WriteLine(title + "");
            }

            row_i = 0;

            for (;;)
            {
                row_i = row_i + 1;

                row_length = 0;

                for (j = 0; j < n; j++)
                {
                    if (a[j] == row_i)
                    {
                        jarray[row_length] = j;
                        row_length = row_length + 1;
                    }

                }

                if (row_length <= 0)
                {
                    break;
                }

                string cout = "";
                
                for (j = 0; j < row_length; j++)
                {
                    cout += (jarray[j] + 1).ToString().PadLeft(6) + "  ";
                }

                Console.WriteLine(cout);

            }
        }

        public static void ytb_random(int n, int[] lambda, ref int seed, ref int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    YTB_RANDOM selects a random Young tableau of a given shape.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 December 2004
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
            //    Input, int N, the integer which has been partitioned.
            //
            //    Input, int LAMBDA[N], is a partition of N, that is,
            //    N = sum ( 0 <= I < N ) LAMBDA[I].
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, int A[N], the vector describing the Young tableau.
            //
        {
            int i;
            int ih;
            int j;
            int k;
            int m;

            for (i = 0; i < n; i++)
            {
                a[i] = 0;
            }

            i = 0;
            k = 0;

            for (;;)
            {
                i = i + 1;

                for (j = 0; j < lambda[i - 1]; j++)
                {
                    a[j] = a[j] + 1;
                    k = k + 1;
                }

                if (n <= k)
                {
                    break;
                }

            }

            for (m = 1; m <= n; m++)
            {

                for (;;)
                {
                    i = UniformRNG.i4_uniform_ab(1, a[0], ref seed);
                    j = UniformRNG.i4_uniform_ab(1, lambda[0], ref seed);

                    if (i <= a[j - 1] && j <= lambda[i - 1])
                    {
                        break;
                    }
                }

                for (;;)
                {
                    ih = a[j - 1] + lambda[i - 1] - i - j;

                    if (ih == 0)
                    {
                        break;
                    }

                    k = UniformRNG.i4_uniform_ab(1, ih, ref seed);

                    if (k <= lambda[i - 1] - j)
                    {
                        j = j + k;
                    }
                    else
                    {
                        i = k - lambda[i - 1] + i + j;
                    }
                }

                lambda[i - 1] = lambda[i - 1] - 1;
                a[j - 1] = a[j - 1] - 1;
                a[n - m] = i;

            }

            for (i = 0; i < n; i++)
            {
                lambda[a[i] - 1] = lambda[a[i] - 1] + 1;
            }
        }
    }
}