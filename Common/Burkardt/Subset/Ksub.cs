using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.SubsetNS
{
    public static class Ksub
    {
        public static void ksub_next(int n, int k, ref int[] a, ref bool more, ref int m, ref int m2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KSUB_NEXT generates the subsets of size K from a set of size N.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 May 2015
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
            //    Input, int N, the size of the set from which subsets are drawn.
            //
            //    Input, int K, the desired size of the subsets.  K must
            //    be between 0 and N.
            //
            //    Output, int A[K].  A[I] is the I-th element of the
            //    subset.  Thus A[I] will be an integer between 1 and N.
            //    Note that the routine will return the values in A
            //    in sorted order: 1 <= A[0] < A[1] < ... < A[K-1] <= N
            //
            //    Input/output, bool &MORE.  Set MORE = FALSE before first call
            //    for a new sequence of subsets.  It then is set and remains
            //    TRUE as long as the subset computed on this call is not the
            //    final one.  When the final subset is computed, MORE is set to
            //    FALSE as a signal that the computation is done.
            //
            //    Input/output, int &M, &M2, two variables used by this
            //    procedure for bookkeeping.  The user must declare these variables,
            //    and the output values from one call must be used as the input values
            //    on the next.  The user should not change these values.
            //
        {
            int j;

            if (k < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("KSUB_NEXT - Fatal error!");
                Console.WriteLine("  K = " + k + "");
                Console.WriteLine("  but 0 <= K is required!");
                return;
            }

            if (n < k)
            {
                Console.WriteLine("");
                Console.WriteLine("KSUB_NEXT - Fatal error!");
                Console.WriteLine("  N = " + n + "");
                Console.WriteLine("  K = " + k + "");
                Console.WriteLine("  but K <= N is required!");
                return;
            }

            if (!more)
            {
                m2 = 0;
                m = k;
            }
            else
            {
                if (m2 < n - m)
                {
                    m = 0;
                }

                m = m + 1;
                m2 = a[k - m];
            }

            for (j = 1; j <= m; j++)
            {
                a[k + j - m - 1] = m2 + j;
            }

            more = (a[0] != (n - k + 1));

        }

        public static void ksub_next2(int n, int k, ref int[] a, ref int in_, ref int iout)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KSUB_NEXT2 generates the subsets of size K from a set of size N.
            //
            //  Discussion:
            //
            //    This routine uses the revolving door method.  It has no "memory".
            //    It simply calculates the successor of the input set,
            //    and will start from the beginning after the last set.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 May 2003
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
            //    Input, int N, the size of the set from which subsets are drawn.
            //    N must be positive.
            //
            //    Input, int K, the size of the desired subset.  K must be
            //    between 0 and N.
            //
            //    Input/output, int A[K].  On input, the user must
            //    supply a subset of size K in A.  That is, A must
            //    contain K unique numbers, in order, between 1 and N.  On
            //    output, A(I) is the I-th element of the output subset.
            //    The output array is also in sorted order.
            //
            //    Output, int &IN, the element of the output subset which
            //    was not in the input set.  Each new subset differs from the
            //    last one by adding one element and deleting another.
            //
            //    Output, int &IOUT, the element of the input subset which
            //    is not in the output subset.
            //
        {
            int j;
            int m;

            if (n <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("KSUB_NEXT2 - Fatal error!");
                Console.WriteLine("  N = " + n + "");
                Console.WriteLine("  but 0 < N is required!");
                return;
            }

            if (k < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("KSUB_NEXT2 - Fatal error!");
                Console.WriteLine("  K = " + k + "");
                Console.WriteLine("  but 0 <= K is required!");
                return;
            }

            if (n < k)
            {
                Console.WriteLine("");
                Console.WriteLine("KSUB_NEXT2 - Fatal error!");
                Console.WriteLine("  N = " + n + "");
                Console.WriteLine("  K = " + k + "");
                Console.WriteLine("  but K <= N is required!");
                return;
            }

            j = 0;

            for (;;)
            {
                if (0 < j || (k % 2) == 0)
                {
                    j = j + 1;

                    if (k < j)
                    {
                        a[k - 1] = k;
                        in_ = k;
                        iout = n;
                        return;
                    }

                    if (a[j - 1] != j)
                    {
                        iout = a[j - 1];
                        in_ = iout - 1;
                        a[j - 1] = in_;

                        if (j != 1)
                        {
                            in_ = j - 1;
                            a[j - 2] = in_;
                        }

                        return;
                    }
                }

                j = j + 1;
                m = n;

                if (j < k)
                {
                    m = a[j] - 1;
                }

                if (m != a[j - 1])
                {
                    break;
                }

            }

            in_ = a[j - 1] + 1;
            a[j - 1] = in_;
            iout = in_ - 1;

            if (j != 1)
            {
                a[j - 2] = iout;
                iout = j - 1;
            }
        }

        public static void ksub_next3(int n, int k, ref int[] a, ref bool more, ref int in_, ref int iout)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KSUB_NEXT3 generates the subsets of size K from a set of size N.
            //
            //  Discussion:
            //
            //    The routine uses the revolving door method.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 May 2003
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
            //    Input, int N, the size of the set from which subsets are drawn.
            //    N must be positive.
            //
            //    Input, int K, the size of the desired subsets.  K must be
            //    between 0 and N.
            //
            //    Output, int A[K].  A(I) is the I-th element of the
            //    output subset.  The elements of A are sorted.
            //
            //    Input/output, bool &MORE.  On first call, set MORE = FALSE
            //    to signal the beginning.  MORE will be set to TRUE, and on
            //    each call, the routine will return another K-subset.
            //    Finally, when the last subset has been returned,
            //    MORE will be set FALSE and you may stop calling.
            //
            //    Output, int &IN, the element of the output subset which
            //    was not in the input set.  Each new subset differs from the
            //    last one by adding one element and deleting another.  IN is not
            //    defined the first time that the routine returns, and is
            //    set to zero.
            //
            //    Output, int &IOUT, the element of the input subset which is
            //    not in the output subset.  IOUT is not defined the first time
            //    the routine returns, and is set to zero.
            //
        {
            int j;
            int m;

            if (n <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("KSUB_NEXT3 - Fatal error!");
                Console.WriteLine("  N = " + n + "");
                Console.WriteLine("  but 0 < N is required!");
                return;
            }

            if (k < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("KSUB_NEXT3 - Fatal error!");
                Console.WriteLine("  K = " + k + "");
                Console.WriteLine("  but 0 <= K is required!");
                return;
            }

            if (n < k)
            {
                Console.WriteLine("");
                Console.WriteLine("KSUB_NEXT3 - Fatal error!");
                Console.WriteLine("  N = " + n + "");
                Console.WriteLine("  K = " + k + "");
                Console.WriteLine("  but K <= N is required!");
                return;
            }

            if (!more)
            {
                in_ = 0;
                iout = 0;
                typeMethods.i4vec_indicator1(k, ref a);
                more = (k != n);
                return;
            }

            j = 0;

            for (;;)
            {
                if (0 < j || (k % 2) == 0)
                {
                    j = j + 1;

                    if (a[j - 1] != j)
                    {
                        iout = a[j - 1];
                        in_ = iout - 1;
                        a[j - 1] = in_;

                        if (j != 1)
                        {
                            in_ = j - 1;
                            a[j - 2] = in_;
                        }

                        if (k != 1)
                        {
                            more = (a[k - 2] == k - 1);
                        }

                        more = (!more) || (a[k - 1] != n);

                        return;
                    }
                }

                j = j + 1;
                m = n;

                if (j < k)
                {
                    m = a[j] - 1;
                }

                if (m != a[j - 1])
                {
                    break;
                }

            }

            in_ = a[j - 1] + 1;
            a[j - 1] = in_;
            iout = in_ - 1;

            if (j != 1)
            {
                a[j - 2] = iout;
                iout = j - 1;
            }

            if (k != 1)
            {
                more = (a[k - 2] == k - 1);
            }

            more = (!more) || (a[k - 1] != n);

        }

        public static void ksub_next4(int n, int k, ref int[] a, ref bool done)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KSUB_NEXT4 generates the subsets of size K from a set of size N.
            //
            //  Discussion:
            //
            //    The subsets are generated one at a time.
            //
            //    The routine should be used by setting DONE to TRUE, and then calling
            //    repeatedly.  Each call returns with DONE equal to FALSE, the array
            //    A contains information defining a new subset.  When DONE returns
            //    equal to TRUE, there are no more subsets.
            //
            //    There are ( N*(N-1)*...*(N+K-1)) / ( K*(K-1)*...*2*1) such subsets.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 August 2018
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
            //    Input, int N, the size of the entire set.
            //
            //    Input, int K, the size of the desired subset.  K must be
            //    between 0 and N.
            //
            //    Input/output, int A[K], contains information about
            //    the subsets.  On the first call with DONE = TRUE, the input contents
            //    of A don't matter.  Thereafter, the input value of A
            //    should be the same as the output value of the previous call.
            //    In other words, leave the array alone!
            //    On output, as long as DONE is returned FALSE, A contains
            //    information defining a subset of K elements of a set of N elements.
            //    In other words, A will contain K distinct numbers (in order)
            //    between 1 and N.
            //
            //    Input/output, bool &DONE.
            //    On the first call, DONE is an input quantity with a value
            //    of TRUE which tells the program to initialize data and
            //    return the first subset.
            //    On return, DONE is an output quantity that is TRUE as long as
            //    the routine is returning another subset, and FALSE when
            //    there are no more.
            //
        {
            int j;
            int jsave;

            if (k < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("KSUB_NEXT4 - Fatal error!");
                Console.WriteLine("  K = " + k + "");
                Console.WriteLine("  but 0 <= K is required!");
                return;
            }

            if (n < k)
            {
                Console.WriteLine("");
                Console.WriteLine("KSUB_NEXT4 - Fatal error!");
                Console.WriteLine("  N = " + n + "");
                Console.WriteLine("  K = " + k + "");
                Console.WriteLine("  but K <= N is required!");
                return;
            }

            //
            //  First call:
            //
            if (done)
            {
                for (j = 0; j < k; j++)
                {
                    a[j] = j + 1;
                }

                done = false;
            }
            //
            //  Empty set returned on previous call.
            //
            else if (n == 0 || k == 0)
            {
                done = true;
            }
            //
            //  Next call.
            //
            else if (a[0] < n - k + 1)
            {
                jsave = k - 1;

                for (j = 0; j < k - 1; j++)
                {
                    if (a[j] + 1 < a[j + 1])
                    {
                        jsave = j;
                        break;
                    }
                }

                for (j = 0; j < jsave; j++)
                {
                    a[j] = j + 1;
                }

                a[jsave] = a[jsave] + 1;
                done = false;
            }
            else
            {
                done = true;
            }
        }

        public static void ksub_random(int n, int k, ref int seed, ref int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KSUB_RANDOM selects a random subset of size K from a set of size N.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 April 2003
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
            //    Input, int N, the size of the set from which subsets are drawn.
            //
            //    Input, int K, number of elements in desired subsets.  K must
            //    be between 0 and N.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, int A[K].  A(I) is the I-th element of the
            //    output set.  The elements of A are in order.
            //
        {
            int i = 0;
            int ids = 0;
            int ihi = 0;
            int ip = 0;
            int ir = 0;
            int is_ = 0;
            int ix = 0;
            int l = 0;
            int ll = 0;
            int m = 0;
            int m0 = 0;

            if (k < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("KSUB_RANDOM - Fatal error!");
                Console.WriteLine("  K = " + k + "");
                Console.WriteLine("  but 0 <= K is required!");
                return;
            }
            else if (n < k)
            {
                Console.WriteLine("");
                Console.WriteLine("KSUB_RANDOM - Fatal error!");
                Console.WriteLine("  N = " + n + "");
                Console.WriteLine("  K = " + k + "");
                Console.WriteLine("  K <= N is required!");
                return;
            }

            if (k == 0)
            {
                return;
            }

            for (i = 1; i <= k; i++)
            {
                a[i - 1] = ((i - 1) * n) / k;
            }

            for (i = 1; i <= k; i++)
            {
                for (;;)
                {
                    ix = UniformRNG.i4_uniform_ab(1, n, ref seed);

                    l = 1 + (ix * k - 1) / n;

                    if (a[l - 1] < ix)
                    {
                        break;
                    }

                }

                a[l - 1] = a[l - 1] + 1;

            }

            ip = 0;
            is_ = k;

            for (i = 1; i <= k; i++)
            {
                m = a[i - 1];
                a[i - 1] = 0;

                if (m != ((i - 1) * n) / k)
                {
                    ip = ip + 1;
                    a[ip - 1] = m;
                }

            }

            ihi = ip;

            for (i = 1; i <= ihi; i++)
            {
                ip = ihi + 1 - i;
                l = 1 + (a[ip - 1] * k - 1) / n;
                ids = a[ip - 1] - ((l - 1) * n) / k;
                a[ip - 1] = 0;
                a[is_ - 1] = l;
                is_ = is_ - ids;
            }

            for (ll = 1; ll <= k; ll++)
            {
                l = k + 1 - ll;

                if (a[l - 1] != 0)
                {
                    ir = l;
                    m0 = 1 + ((a[l - 1] - 1) * n) / k;
                    m = (a[l - 1] * n) / k - m0 + 1;
                }

                //
                //  There is something wrong with this algorithm!
                //  If A[L-1] is zero, then the values of IR, M0, and M are not defined
                //  on this loop iteration, and hence are either STALE values from the
                //  previous iteration, or UNDEFINED if this is the first pass.
                //  JVB, 21 December 2014.
                //
                ix = UniformRNG.i4_uniform_ab(m0, m0 + m - 1, ref seed);

                i = l + 1;

                while (i <= ir)
                {
                    if (ix < a[i - 1])
                    {
                        break;
                    }

                    ix = ix + 1;
                    a[i - 2] = a[i - 1];
                    i = i + 1;
                }

                a[i - 2] = ix;
                m = m - 1;
            }
        }

        public static void ksub_random2(int n, int k, ref int seed, ref int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KSUB_RANDOM2 selects a random subset of size K from a set of size N.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 May 2003
            //
            //  Author:
            //
            //    FORTRAN77 original version by Albert Nijenhuis, Herbert Wilf.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    A Nijenhuis and H Wilf,
            //    Combinatorial Algorithms,
            //    Academic Press, 1978, second edition,
            //    ISBN 0-12-519260-6.
            //
            //  Parameters:
            //
            //    Input, int N, the size of the set from which subsets are drawn.
            //
            //    Input, int K, number of elements in desired subsets.  K must
            //    be between 0 and N.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, int A[K].  A(I) is the I-th element of the
            //    output set.  The elements of A are in order.
            //
        {
            int available;
            int candidate;
            int have;
            int need;
            double r;

            if (k < 0 || n < k)
            {
                Console.WriteLine("");
                Console.WriteLine("KSUB_RANDOM2 - Fatal error!");
                Console.WriteLine("  N = " + n + "");
                Console.WriteLine("  K = " + k + "");
                Console.WriteLine("  but 0 <= K <= N is required!");
                return;
            }

            if (k == 0)
            {
                return;
            }

            need = k;
            have = 0;
            available = n;
            candidate = 0;

            for (;;)
            {
                candidate = candidate + 1;

                r = UniformRNG.r8_uniform_01(ref seed);

                if (r * (double)available <= (double)need)
                {
                    need = need - 1;
                    a[have] = candidate;
                    have = have + 1;

                    if (need <= 0)
                    {
                        break;
                    }

                }

                available = available - 1;

            }
        }

        public static void ksub_random3(int n, int k, ref int seed, ref int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KSUB_RANDOM3 selects a random subset of size K from a set of size N.
            //
            //  Discussion:
            //
            //    This routine uses Floyd's algorithm.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 May 2003
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
            //    Input, int N, the size of the set from which subsets are drawn.
            //
            //    Input, int K, number of elements in desired subsets.  K must
            //    be between 0 and N.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, int A[N].  I is an element of the subset
            //    if A(I) = 1, and I is not an element if A(I)=0.
            //
        {
            int i;
            int j;

            if (k < 0 || n < k)
            {
                Console.WriteLine("");
                Console.WriteLine("KSUB_RANDOM3 - Fatal error!");
                Console.WriteLine("  N = " + n + "");
                Console.WriteLine("  K = " + k + "");
                Console.WriteLine("  but 0 <= K <= N is required!");
                return;
            }

            for (i = 0; i < n; i++)
            {
                a[i] = 0;
            }

            if (k == 0)
            {
                return;
            }

            for (i = n - k + 1; i <= n; i++)
            {
                j = UniformRNG.i4_uniform_ab(1, i, ref seed);

                if (a[j - 1] == 0)
                {
                    a[j - 1] = 1;
                }
                else
                {
                    a[i - 1] = 1;
                }
            }

            return;
        }

        public static void ksub_random4(int n, int k, ref int seed, ref int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KSUB_RANDOM4 selects a random subset of size K from a set of size N.
            //
            //  Discussion:
            //
            //    This routine is somewhat impractical for the given problem, but
            //    it is included for comparison, because it is an interesting
            //    approach that is superior for certain applications.
            //
            //    The approach is mainly interesting because it is "incremental";
            //    it proceeds by considering every element of the set, and does not
            //    need to know how many elements there are.
            //
            //    This makes this approach ideal for certain cases, such as the
            //    need to pick 5 lines at random from a text file of unknown length,
            //    or to choose 6 people who call a certain telephone number on a
            //    given day.  Using this technique, it is possible to make the
            //    selection so that, whenever the input stops, a valid uniformly
            //    random subset has been chosen.
            //
            //    Obviously, if the number of items is known in advance, and
            //    it is easy to extract K items directly, there is no need for
            //    this approach, and it is less efficient since, among other costs,
            //    it has to generate a random number for each item, and make an
            //    acceptance/rejection test.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 July 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Tom Christiansen, Nathan Torkington,
            //    "8.6: Picking a Random Line from a File",
            //    Perl Cookbook, pages 284-285,
            //    O'Reilly, 1999.
            //
            //  Parameters:
            //
            //    Input, int N, the size of the set from which subsets are drawn.
            //
            //    Input, int K, number of elements in desired subsets.  K must
            //    be between 0 and N.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, int A[K], contains the indices of the selected items.
            //
        {
            int i;
            int j;
            int next;
            double r;

            next = 0;
            //
            //  Here, we use a WHILE to suggest that the algorithm
            //  proceeds to the next item, without knowing how many items
            //  there are in total.
            //
            //  Note that this is really the only place where N occurs,
            //  so other termination criteria could be used, and we really
            //  don't need to know the value of N!
            //
            while (next < n)
            {
                next = next + 1;

                if (next <= k)
                {
                    i = next;
                    a[i - 1] = next;
                }
                else
                {
                    r = UniformRNG.r8_uniform_01(ref seed);

                    if (r * (double)next <= (double)k)
                    {
                        i = UniformRNG.i4_uniform_ab(1, k, ref seed);
                        //
                        //  If we slide the current items down, and insert at the end, we preserve order.
                        //
                        for (j = i; j < k; j++)
                        {
                            a[j - 1] = a[j];
                        }

                        a[k - 1] = next;
                        //      a[i-1] = next;
                    }
                }
            }
        }

        public static int[] ksub_random5(int n, int k, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KSUB_RANDOM5 selects a random subset of size K from a set of size N.
            //
            //  Discussion:
            //
            //    Consider the set A = 1, 2, 3, ... N.  
            //    Choose a random index I1 between 1 and N, and swap items A(1) and A(I1).
            //    Choose a random index I2 between 2 and N, and swap items A(2) and A(I2).
            //    repeat K times.
            //    A(1:K) is your random K-subset.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    09 June 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the size of the set from which subsets
            //    are drawn.
            //
            //    Input, int K, number of elements in desired subsets.
            //    1 <= K <= N.
            //
            //    Input/output, int &SEED, a seed for the random
            //    number generator.
            //
            //    Output, int KSUB_RANDOM5[K], the indices of the randomly
            //    chosen elements.  These are 1-based indices.
            //
        {
            int[] a;
            int[] b;
            int base_ = 1;
            int i;
            int j;
            int t;
            //
            //  Let B index the set.
            //
            b = new int[n];

            for (i = 0; i < n; i++)
            {
                b[i] = i + base_;
            }

            //
            //  Choose item 1 from N things,
            //  choose item 2 from N-1 things,
            //  choose item K from N-K+1 things.
            //
            for (i = 0; i < k; i++)
            {
                j = UniformRNG.i4_uniform_ab(i, n - 1, ref seed);
                t = b[i];
                b[i] = b[j];
                b[j] = t;
            }

            //
            //  Copy the first K elements.
            //
            a = new int[k];

            for (i = 0; i < k; i++)
            {
                a[i] = b[i];
            }

            //
            //  Put the elements in ascending order.
            //
            typeMethods.i4vec_sort_heap_a(k, ref a);

            return a;
        }

        public static void ksub_rank(int k, int[] a, ref int rank)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KSUB_RANK computes the rank of a subset of an N set.
            //
            //  Discussion:
            //
            //    The routine accepts an array representing a subset of size K from a set
            //    of size N, and returns the rank (or order) of that subset.  
            //
            //    This is the same order in which routine KSUB_NEXT2 would produce that subset.
            //
            //    Note the value of N is not input, and is not, in fact,
            //    needed.
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
            //    Input, int K, the number of elements in the subset.
            //
            //    Input, int A[K], contains K distinct numbers between
            //    1 and N, in order.
            //
            //    Output, int &RANK, the rank of this subset.
            //
        {
            int i;
            int iprod;
            int j;

            rank = 0;

            for (i = 1; i <= k; i++)
            {
                iprod = 1;

                for (j = i + 1; j <= a[i - 1] - 1; j++)
                {
                    iprod = iprod * j;
                }

                for (j = 1; j <= a[i - 1] - i - 1; j++)
                {
                    iprod = iprod / j;
                }

                if (a[i - 1] == 1)
                {
                    iprod = 0;
                }

                rank = rank + iprod;
            }

            rank = rank + 1;

        }

        public static void ksub_to_comp(int ns, int ks, int[] as_, ref int nc, ref int kc, ref int[] ac)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KSUB_TO_COMP converts a K-subset to a composition.
            //
            //  Discussion:
            //
            //    There is a bijection between K subsets and compositions.
            //
            //    Because we allow a composition to have entries that are 0, we need
            //    to implicitly add 1 to each entry before establishing the bijection.
            //
            //    Let AS be a KS subset of a set of the integers 1 through NS.
            //
            //    Then let 
            //      NC = NS - KS, 
            //      KC = KS + 1, 
            //    and define
            //      AC(1) = AS(1) - 1;
            //      AC(2:KC-1) = AS(2:KC-1) - AS(1:KC-2) - 1;
            //      AC(KC) = NS - AS(KS).
            //
            //    Then AC is a composition of NC into KC parts.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 December 2013
            //
            //  Parameters:
            //
            //    Input, int NS, the size of the set.
            //
            //    Input, int KS, the size of the subset.
            //
            //    Input, int AS[KS], the entries of the K-subset, 
            //    in increasing order.
            //
            //    Output, int &NC, the composition sum.
            //
            //    Output, int &KC, the number of parts of the composition.
            //
            //    Output, int AC[KC], the parts of the composition.
            //
        {
            int i;

            nc = ns - ks;
            kc = ks + 1;

            ac[0] = as_[0] - 1;
            for (i = 1; i < kc - 1; i++)
            {
                ac[i] = as_[i] - as_[i - 1] - 1;
            }

            ac[kc - 1] = ns - as_[ks - 1];

        }

        public static void ksub_to_compnz(int ns, int ks, int[] as_, ref int nc, ref int kc, ref int[] ac)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KSUB_TO_COMPNZ converts a K-subset to a nonzero composition.
            //
            //  Discussion:
            //
            //    There is a bijection between K subsets and compositions.
            //
            //    Let AS be a KS subset of a set of the integers 1 through NS.
            //
            //    Then let 
            //      NC = NS + 1 
            //      KC = KS + 1, 
            //    and define
            //      AC(1) = AS(1);
            //      AC(2:KC-1) = AS(2:KC-1) - AS(1:KC-2);
            //      AC(KC) = NC - AS(KS).
            //
            //    Then AC is a composition of NC into KC parts.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 May 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NS, the size of the set.
            //
            //    Input, int KS, the size of the subset.
            //
            //    Input, int AS[KS], the entries of the K-subset, 
            //    in increasing order.
            //
            //    Output, int &NC, the composition sum.
            //
            //    Output, int &KC, the number of parts of the composition.
            //
            //    Output, int AC[KC], the parts of the composition.
            //
        {
            int i;

            nc = ns + 1;
            kc = ks + 1;

            ac[0] = as_[0];
            for (i = 1; i < kc - 1; i++)
            {
                ac[i] = as_[i] - as_[i - 1];
            }

            ac[kc - 1] = nc - as_[ks - 1];

        }

        public static void ksub_unrank(int k, int rank, ref int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KSUB_UNRANK returns the subset of a given rank.
            //
            //  Discussion:
            //
            //    The routine is given a rank and returns the corresponding subset of K
            //    elements of a set of N elements.  
            //
            //    It uses the same ranking that KSUB_NEXT2 uses to generate all the subsets 
            //    one at a time.  
            //
            //    Note that the value of N itself is not input, nor is it needed.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 June 2004
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
            //    Input, int K, the number of elements in the subset.
            //
            //    Input, int RANK, the rank of the desired subset.
            //    There are ( N*(N-1)*...*(N+K-1)) / ( K*(K-1)*...*2*1) such
            //    subsets, so RANK must be between 1 and that value.
            //
            //    Output, int A[K], K distinct integers in order between
            //    1 and N, which define the subset.
            //
        {
            int i;
            int ip;
            int iprod;
            int jrank;

            jrank = rank - 1;

            for (i = k; 1 <= i; i--)
            {
                ip = i - 1;
                iprod = 1;

                for (;;)
                {
                    ip = ip + 1;

                    if (ip != i)
                    {
                        iprod = (ip * iprod) / (ip - i);
                    }

                    if (jrank < iprod)
                    {
                        break;
                    }
                }

                if (ip != i)
                {
                    iprod = ((ip - i) * iprod) / ip;
                }

                jrank = jrank - iprod;
                a[i - 1] = ip;
            }

        }


    }
}