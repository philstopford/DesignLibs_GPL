﻿using System;
using Burkardt.MatrixNS;

namespace Burkardt.SubsetNS
{
    public static class Triang
    {
        public static void triang(int n, int[] zeta, ref int[] p )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANG renumbers elements in accordance with a partial ordering.
        //
        //  Discussion:
        //
        //    TRIANG is given a partially ordered set.  The partial ordering
        //    is defined by a matrix ZETA, where element I is partially less than
        //    or equal to element J if and only if ZETA(I,J) = 1.
        //
        //    TRIANG renumbers the elements with a permutation P so that if
        //    element I is partially less than element J in the partial ordering,
        //    then P(I) < P(J) in the usual, numerical ordering.
        //
        //    In other words, the elements are relabeled so that their labels
        //    reflect their ordering.  This is equivalent to relabeling the
        //    matrix so that, on unscrambling it, the matrix would be upper
        //    triangular.
        //
        //    Calling I4MAT_2PERM0 or R8MAT_2PERM0 with P used for both the row
        //    and column permutations applied to matrix ZETA will result in
        //    an upper triangular matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 June 2015
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
        //    Input, int N, the number of elements in the set.
        //
        //    Input, int ZETA[N*N], describes the partial ordering.
        //    ZETA[I,J] =:
        //      0, for diagonal elements (I = J), or
        //         for unrelated elements, or
        //         if J << I.
        //      1, if I << J.
        //
        //    Output, int P[N], a permutation of the elements that reflects
        //    their partial ordering.  P[I] is the new label of element I, with
        //    the property that if ZETA[I,J] = 1, that is, I << J,
        //    then P[I] < P[J] (in the usual ordering).
        //
        {
            int i;
            bool error;
            int iq;
            int ir;
            int it;
            int l;
            int m;
            //
            //  Make sure ZETA represents a partially ordered set.  In other words,
            //  if ZETA(I,J) = 1, then ZETA(J,I) must NOT be 1.
            //
            error = PartialOrdering.pord_check(n, zeta);

            if (error)
            {
                Console.WriteLine("");
                Console.WriteLine("TRIANG - Fatal error!");
                Console.WriteLine("  The matrix ZETA does not represent a");
                Console.WriteLine("  partial ordering.");
                return;
            }

            m = 1;
            l = 0;
            for (i = 0; i < n; i++)
            {
                p[i] = 0;
            }

            it = m + 1;
            ir = m + 1;

            for (;;)
            {
                if (ir <= n)
                {
                    if (p[ir - 1] == 0 && zeta[(ir - 1) + (m - 1) * n] != 0)
                    {
                        p[ir - 1] = m;
                        m = ir;
                        ir = it;
                    }
                    else
                    {
                        ir = ir + 1;
                    }
                }
                else
                {
                    l = l + 1;
                    iq = p[m - 1];
                    p[m - 1] = l;

                    if (iq != 0)
                    {

                        ir = m + 1;
                        m = iq;
                    }
                    else if (m == n)
                    {
                        break;
                    }
                    else
                    {
                        for (;;)
                        {
                            m = m + 1;

                            if (p[m - 1] == 0)
                            {
                                break;
                            }

                            if (m == n)
                            {
                                for (i = 0; i < n; i++)
                                {
                                    p[i] = p[i] - 1;
                                }

                                return;
                            }
                        }

                        it = m + 1;
                        ir = m + 1;
                    }
                }
            }

            //
            //  Decrement the elements of the permutation.
            //
            for (i = 0; i < n; i++)
            {
                p[i] = p[i] - 1;
            }
        }
    }
}