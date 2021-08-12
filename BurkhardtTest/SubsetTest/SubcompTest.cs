using System;
using Burkardt;

namespace SubsetTestNS
{
    public static class SubcompTest
    {
        public static void subcomp_next(int n, int k, ref int[] a, ref bool more,
                ref int h, ref int t, ref int n2, ref bool more2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SUBCOMP_NEXT computes the next subcomposition of N into K parts.
            //
            //  Discussion:
            //
            //    A composition of the integer N into K parts is an ordered sequence
            //    of K nonnegative integers which sum to a value of N.
            //
            //    A subcomposition of the integer N into K parts is a composition
            //    of M into K parts, where 0 <= M <= N.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 June 2015
            //
            //  Author:
            //
            //    John Burkardt.
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
            //    Input, int N, the integer whose subcompositions are desired.
            //
            //    Input, int K, the number of parts in the subcomposition.
            //
            //    Input/output, int A[K], the parts of the subcomposition.
            //
            //    Input/output, bool &MORE, set by the user to start the computation,
            //    and by the routine to terminate it.
            //
            //    Input/output, int &H, &T, &N2, internal parameters needed for the
            //    computation.  The user should allocate space for these in the calling
            //    program, include them in the calling sequence, but never alter them!
            //
            //    Input/output, bool &MORE2, an internal parameter needed for the
            //    computation.  The user should allocate space for it in the calling
            //    program, include it in the calling sequence, but never alter it!
            //
        {
            int i;
            //
            //  The first computation.
            //
            if (!more)
            {
                for (i = 0; i < k; i++)
                {
                    a[i] = 0;
                }

                more = true;

                h = 0;
                t = 0;
                n2 = 0;
                more2 = false;
            }
            //
            //  Do the next element at the current value of N.
            //
            else if (more2)
            {
                Comp.comp_next(n2, k, ref a, ref more2, ref h, ref t);
            }
            else
            {
                more2 = false;
                n2 = n2 + 1;

                Comp.comp_next(n2, k, ref a, ref more2, ref h, ref t);
            }

            //
            //  Termination occurs if MORE2 = FALSE and N2 = N.
            //
            if (!more2 && n2 == n)
            {
                more = false;
            }
        }

        public static void subcompnz_next(int n, int k, ref int[] a, ref bool more,
                ref int h, ref int t, ref int n2, ref bool more2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SUBCOMPNZ_NEXT computes the next subcomposition of N into K nonzero parts.
            //
            //  Discussion:
            //
            //    A composition of the integer N into K nonzero parts is an ordered sequence
            //    of K positive integers which sum to a value of N.
            //
            //    A subcomposition of the integer N into K nonzero parts is a composition
            //    of M into K nonzero parts, where 0 < M <= N.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 June 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the integer whose subcompositions are desired.
            //
            //    Input, int K, the number of parts in the subcomposition.
            //    K must be no greater than N.
            //
            //    Input/output, int A[K], the parts of the subcomposition.
            //
            //    Input/output, bool &MORE, set by the user to start the computation,
            //    and by the routine to terminate it.
            //
            //    Input/output, int &H, &T, &N2, internal parameters needed for the
            //    computation.  The user should allocate space for these in the calling
            //    program, include them in the calling sequence, but never alter them!
            //
            //    Input/output, bool &MORE2, an internal parameter needed for the
            //    computation.  The user should allocate space for it in the calling
            //    program, include it in the calling sequence, but never alter it!
            //
        {
            int i;
            CompNZData data = new CompNZData();

            if (n < k)
            {
                for (i = 0; i < k; i++)
                {
                    a[i] = -1;
                }

                return;
            }

            //
            //  The first computation.
            //
            if (!more)
            {
                for (i = 0; i < k; i++)
                {
                    a[i] = 1;
                }

                more = true;
                h = 0;
                t = 0;
                n2 = k;
                more2 = false;
            }
            //
            //  Do the next element at the current value of N.
            //
            else if (more2)
            {
                Comp.compnz_next(ref data, n2, k, ref a, ref more2);
            }
            else
            {
                more2 = false;
                n2 = n2 + 1;

                Comp.compnz_next(ref data, n2, k, ref a, ref more2);
            }

            //
            //  Termination occurs if MORE2 = FALSE and N2 = N.
            //
            if (!more2 && n2 == n)
            {
                more = false;
            }
        }

        public static void subcompnz2_next(int n_lo, int n_hi, int k, ref int[] a, ref bool more,
                ref int h, ref int t, ref int n2, ref bool more2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SUBCOMPNZ2_NEXT computes the next subcomposition of N into K nonzero parts.
            //
            //  Discussion:
            //
            //    A composition of the integer N into K nonzero parts is an ordered sequence
            //    of K positive integers which sum to a value of N.
            //
            //    A subcomposition of the integer N into K nonzero parts is a composition
            //    of M into K nonzero parts, where 0 < M <= N.
            //
            //    This routine computes all compositions of K into nonzero parts which sum
            //    to values between N_LO and N_HI.
            //
            //    The routine SUBCOMPNZ_NEXT can be regarded as a special case 
            //    where N_LO = K.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 June 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N_LO, N_HI, the range of values of N for which compositions
            //    are desired.
            //    N_LO must be no greater than N_HI.
            //
            //    Input, int K, the number of parts in the subcomposition.
            //    K must be no greater than N_HI.
            //
            //    Input/output, int A[K], the parts of the subcomposition.
            //
            //    Input/output, bool &MORE, set by the user to start the computation,
            //    and by the routine to terminate it.
            //
            //    Input/output, int &H, &T, &N2, internal parameters needed for the
            //    computation.  The user should allocate space for these in the calling
            //    program, include them in the calling sequence, but never alter them!
            //
            //    Input/output, bool &MORE2, an internal parameter needed for the
            //    computation.  The user should allocate space for it in the calling
            //    program, include it in the calling sequence, but never alter it!
            //
        {
            int i;
            CompNZData data = new CompNZData();

            if (n_hi < k)
            {
                for (i = 0; i < k; i++)
                {
                    a[i] = -1;
                }

                return;
            }

            if (n_hi < n_lo)
            {
                for (i = 0; i < k; i++)
                {
                    a[i] = -1;
                }

                return;
            }

            //
            //  The first computation.
            //
            if (!more)
            {
                more = true;
                h = 0;
                t = 0;
                n2 = Math.Max(k, n_lo);
                more2 = false;

                Comp.compnz_next(ref data, n2, k, ref a, ref more2);
            }
            //
            //  Do the next element at the current value of N.
            //
            else if (more2)
            {
                Comp.compnz_next(ref data, n2, k, ref a, ref more2);
            }
            else
            {
                n2 = n2 + 1;

                Comp.compnz_next(ref data, n2, k, ref a, ref more2);
            }

            //
            //  Termination occurs if MORE2 = FALSE and N2 = N_HI.
            //
            if (!more2 && n2 == n_hi)
            {
                more = false;
            }
        }

    }
}