using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static bool r8vec_mirror_next(int n, ref double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_MIRROR_NEXT steps through all sign variations of an R8VEC.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //    In normal use, the user would set every element of A to be positive.
            //    The routine will take the input value of A, and output a copy in
            //    which the signs of one or more entries have been changed.  Repeatedly
            //    calling the routine with the output from the previous call will generate
            //    every distinct "variation" of A; that is, all possible sign variations.
            //
            //    When the output variable DONE is TRUE (or equal to 1), then the
            //    output value of A_NEW is the last in the series.
            //
            //    Note that A may have some zero values.  The routine will essentially
            //    ignore such entries; more exactly, it will not stupidly assume that -0
            //    is a proper "variation" of 0.
            //
            //    Also, it is possible to call this routine with the signs of A set
            //    in any way you like.  The routine will operate properly, but it
            //    will nonethess terminate when it reaches the value of A in which
            //    every nonzero entry has negative sign.
            //
            //
            //    More efficient algorithms using the Gray code seem to require internal
            //    memory in the routine, which is not one of MATLAB's strong points,
            //    or the passing back and forth of a "memory array", or the use of
            //    global variables, or unnatural demands on the user.  This form of
            //    the routine is about as clean as I can make it.
            //
            //  Example:
            //
            //      Input         Output
            //    ---------    --------------
            //    A            A         DONE
            //    ---------    --------  ----
            //     1  2  3     -1  2  3  false
            //    -1  2  3      1 -2  3  false
            //     1 -2  3     -1 -2  3  false
            //    -1 -2  3      1  2 -3  false
            //     1  2 -3     -1  2 -3  false
            //    -1  2 -3      1 -2 -3  false
            //     1 -2 -3     -1 -2 -3  false
            //    -1 -2 -3      1  2  3  true
            //
            //     1  0  3     -1  0  3  false
            //    -1  0  3      1  0 -3  false
            //     1  0 -3     -1  0 -3  false
            //    -1  0 -3      1  0  3  true
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 March 2008
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
            //    Input, int N, the number of entries in the vector.
            //
            //    Input/output, double A[N], a vector of real numbers.  On
            //    output, some signs have been changed.
            //
            //    Output, bool R8VEC_MIRROR_NEXT, is TRUE if the input vector A was
            //    the last element
            //    in the series (every entry was nonpositive); the output vector is reset
            //    so that all entries are nonnegative, but presumably the ride is over.
            //
        {
            bool done;
            int i;
            int positive;
            //
            //  Seek the first strictly positive entry of A.
            //
            positive = -1;
            for (i = 0; i < n; i++)
            {
                if (0.0 < a[i])
                {
                    positive = i;
                    break;
                }
            }

            //
            //  If there is no strictly positive entry of A, there is no successor.
            //
            if (positive == -1)
            {
                for (i = 0; i < n; i++)
                {
                    a[i] = -a[i];
                }

                done = true;
                return done;
            }

            //
            //  Otherwise, negate A up to the positive entry.
            //
            for (i = 0; i <= positive; i++)
            {
                a[i] = -a[i];
            }

            done = false;

            return done;
        }

        public static void r8vec_mirror_ab_next(int m, double[] a, double[] b, ref double[] x,
                ref bool done)

            //*****************************************************************************/
            //
            //  Purpose:
            //
            //    R8VEC_MIRROR_AB_NEXT steps through "mirrored" versions of vector X.
            //
            //  Discussion:
            //
            //    X is an M component vector contained in a rectangle described by A and B,
            //    that is, for each index I, we have
            //      A(I) <= X(I) <= B(I).
            //
            //    "Mirrored" versions of the vector X have one or more components
            //    reflected about the A or B limit.  
            //
            //    As long as each component of X is strictly between the limits A and B,
            //    this means there will be 3^M possible versions of the vector.
            //
            //    If one component of X is equal to one limit or the other, this 
            //    suppresses mirroring across that limit.  If one component of
            //    X, A and B are equal, then no mirroring is done at all in that component.
            //
            //  Example:
            //
            //      A = 0, 0, 0
            //      X = 1, 1, 1
            //      B = 2, 2, 2
            //      results in the following sequence of 3x3x3 values:
            //
            //      0  0  0
            //      0  0  1
            //      0  0  2
            //      0  1  0
            //      0  1  1
            //      .......
            //      2  1  1
            //      2  1  2
            //      2  2  0
            //      2  2  1
            //      2  2  2
            //
            //    A = 0 1 0
            //    X = 1 1 1
            //    B = 2 2 2
            //    results in the following sequence of 3x2x3 values:
            //
            //      0 1 0
            //      0 1 1
            //      0 1 2
            //      0 2 0
            //      0 2 1
            //      0 2 2
            //      1 1 0
            //      1 1 1
            //      1 1 2
            //      1 2 0
            //      1 2 1
            //      1 2 2
            //      2 1 0
            //      2 1 1
            //      2 1 2
            //      2 2 0
            //      2 2 1
            //      2 2 2
            //
            //    A = 0 1 0
            //    X = 1 1 1
            //    B = 2 1 2
            //    results in the following sequence of 3x1x3 values:
            //
            //      0 1 0
            //      0 1 1
            //      0 1 2
            //      1 1 0
            //      1 1 1
            //      1 1 2
            //      2 1 0
            //      2 1 1
            //      2 1 2
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 August 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the number of entries in the vector.
            //
            //    Input, double A[M], B[M], the lower and upper limits.
            //
            //    Input/output, double X[M], a vector being manipulated.
            //
            //    Input/output, bool &DONE.  On first call, DONE should be TRUE, and
            //    A(I) <= X(I) <= B(I) for each index I.  On output, if DONE is TRUE,
            //    then the returned value is the last entry in the sequence.
            //
        {
            int i;
            //
            //  First call:
            //
            if (done)
            {
                //
                //  Ensure all A(I) <= X(I) <= B(I).
                //
                for (i = 0; i < m; i++)
                {
                    if (x[i] < a[i])
                    {
                        Console.WriteLine("");
                        Console.WriteLine("R8VEC_MIRROR_AB_NEXT - Fatal error!");
                        Console.WriteLine("  Not every A(I) <= X(I).");
                        return;
                    }

                    if (b[i] < x[i])
                    {
                        Console.WriteLine("");
                        Console.WriteLine("R8VEC_MIRROR_AB_NEXT - Fatal error!");
                        Console.WriteLine("  Not every X(I) <= B(I).");
                        return;
                    }
                }

                //
                //  Set first element of sequence.
                //
                for (i = 0; i < m; i++)
                {
                    x[i] = 2.0 * a[i] - x[i];
                }

                //
                //  Unless A = B, our sequence is not done.
                //
                done = true;
                for (i = 0; i < m; i++)
                {
                    if (a[i] != b[i])
                    {
                        done = false;
                        break;
                    }
                }
            }
            //
            //  Subsequent calls.
            //
            else
            {
                //
                //  Initialize index to last.
                //  loop
                //    if index < 1, set DONE = true and return.
                //    if the index-th value is below B, increment it and return;
                //    otherwise reset index-th value to A.
                //    decrement the index.
                //
                i = m - 1;

                while (0 <= i)
                {
                    if (x[i] < a[i])
                    {
                        x[i] = 2.0 * a[i] - x[i];
                        return;
                    }
                    else if (x[i] < b[i])
                    {
                        x[i] = 2.0 * b[i] - x[i];
                        return;
                    }
                    else
                    {
                        x[i] = x[i] - 2.0 * (b[i] - a[i]);
                    }

                    i = i - 1;
                }

                done = true;
            }

        }
    }
}