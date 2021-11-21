﻿using System;

namespace Burkardt.AppliedStatistics;

public static partial class Algorithms
{
    public static void simplex_lattice_point_next(int n, int t, ref bool more, ref int[] x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SIMPLEX_LATTICE_POINT_NEXT generates lattice points in a simplex.
        //
        //  Discussion:
        //
        //    The simplex is defined by N-dimensional points X such that:
        //
        //        0 <= X(1:N)
        //
        //    and
        //
        //      sum ( X(1:N) ) <= T
        //
        //    where T is an integer.
        //
        //    Lattice points are points X which satisfy the simplex conditions and
        //    for which all the components are integers.
        //
        //    This routine generates all the lattice points in a given simplex, one at 
        //    a time, in a reverse lexicographic order.
        //
        //    To use the routine, initialize by setting N and T to appropriate values, 
        //    and MORE to FALSE.  The initial value of X is not important.
        //
        //    Call the routine. On return, X will contain the first lattice point in 
        //    the simplex.  If MORE is TRUE, then the routine may be called again to 
        //    get the next point.  In fact, as long as the output value of MORE is 
        //    TRUE, there is at least one more lattice point that can be found by 
        //    making another call.  When MORE is returned as FALSE, then there are no 
        //    more lattice points; the value of X returned at that time is the 
        //    "last" such point.
        //
        //    During the computation of a sequence of lattice points, the user should 
        //    not change the values of N, T, MORE or X.  
        //
        //    The output for N = 3, T = 4 would be:
        //
        //       1    4  0  0
        //       2    3  1  0
        //       3    3  0  1
        //       4    2  2  0
        //       5    2  1  1
        //       6    2  0  2
        //       7    1  3  0
        //       8    1  2  1
        //       9    1  1  2
        //      10    1  0  3
        //      11    0  4  0
        //      12    0  3  1
        //      13    0  2  2
        //      14    0  1  3
        //      15    0  0  4
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Scott Chasalow, Richard Brand,
        //    Algorithm AS 299:
        //    Generation of Simplex Lattice Points,
        //    Applied Statistics,
        //    Volume 44, Number 4, 1995, pages 534-545.
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
        //    Input, int N, the spatial dimension.
        //    N must be positive.
        //
        //    Input, int T, the characteristic of the simplex.
        //    T must be nonnegative.
        //
        //    Input/output, int *MORE, initialized to FALSE by the user to
        //    begin a sequence of calculations, returned by the routine as TRUE,
        //    if there are more values of X that can be calculated, or FALSE
        //    if the accompanying value of X is the last one for this sequence.
        //
        //    Input/output, int X[N], not initialized by the user, but not
        //    changed by the user on subsequent calls.  The routine returns
        //    a new point on each call, and on subsequent calls uses the input
        //    value (old point) to compute the output value (next point).
        //
    {
        switch (more)
        {
            case false when n < 1:
                Console.WriteLine("");
                Console.WriteLine("SIMPLEX_LATTICE_POINT_NEXT - Fatal error!");
                Console.WriteLine("  N < 1.");
                return;
            case false when t < 0:
                Console.WriteLine("");
                Console.WriteLine("SIMPLEX_LATTICE_POINT_NEXT - Fatal error!");
                Console.WriteLine("  T < 0.");
                return;
            case false:
            {
                more = true;

                x[0] = t;
                for (int i = 1; i < n; i++)
                {
                    x[i] = 0;
                }

                more = n switch
                {
                    //
                    //  The first point can actually also be the last!
                    //
                    1 => false,
                    _ => more
                };

                break;
            }
            default:
            {
                //
                //  Search X(N-1 down to 1) for the first nonzero element.
                //  If none, then terminate.  (This should not happen!)
                //  Otherwise, set J to this index.
                //  Decrement X(J) by 1.
                //  Set X(J+1:N) to (T-X(1:J),0,0,...0).
                //
                int j = n - 1;

                for (int i = n - 2; 0 <= i; i--)
                {
                    if (0 >= x[i])
                    {
                        continue;
                    }

                    j = i;
                    break;
                }

                if (j == n - 1)
                {
                    Console.WriteLine("");
                    Console.WriteLine("SIMPLEX_LATTICE_POINT_NEXT - Fatal error!");
                    Console.WriteLine("  The input X vector is nonpositive in all entries");
                    Console.WriteLine("  except possibly the last one.");
                    Console.WriteLine("");
                    Console.WriteLine("  Perhaps the user has miscalled the routine");
                    Console.WriteLine("  or altered data between calls.");
                    Console.WriteLine("");
                    Console.WriteLine("ABNORMAL TERMINATION.");
                    return;
                }

                x[j] -= 1;
                x[j + 1] = t;
                for (int i = 0; i <= j; i++)
                {
                    x[j + 1] -= x[i];
                }

                for (int i = j + 2; i < n; i++)
                {
                    x[i] = 0;
                }

                //
                //  Is this the last point?
                //
                if (x[n - 1] == t)
                {
                    more = false;
                }

                break;
            }
        }
    }
}