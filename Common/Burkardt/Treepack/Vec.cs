using System;
using Burkardt.Uniform;

namespace Burkardt.Treepack
{
    public static class Vec
    {
        public static void vec_next(int n, int ibase, ref int[] iarray, ref bool more)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VEC_NEXT generates all N-vectors of integers modulo a given base.
            //
            //  Discussion:
            //
            //    The items are produced one at a time.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 August 2013
            //
            //  Parameters:
            //
            //    Input, int N, the size of the vectors to be used.
            //
            //    Input, int IBASE, the base to be used.  IBASE = 2 will
            //    give vectors of 0's and 1's, for instance.
            //
            //    Input/output, int IARRAY[N].  On each return,
            //    IARRAY will contain entries in the range 0 to IBASE-1.
            //
            //    Input/output, ref int MORE.  Set this variable 0 before
            //    the first call.  Normally, MORE will be returned 1 but
            //    once all the vectors have been generated, MORE will be
            //    reset 0 and you should stop calling the program.
            //
        {
            int i;
            int kount = 0;
            int last = 0;
            int nn;

            if (!more)
            {
                kount = 1;
                last = (int) Math.Pow(ibase, n);
                more = true;
                for (i = 0; i < n; i++)
                {
                    iarray[i] = 0;
                }
            }
            else
            {
                kount = kount + 1;

                if (kount == last)
                {
                    more = false;
                }

                iarray[n - 1] = iarray[n - 1] + 1;

                for (i = 1; i <= n; i++)
                {
                    nn = n - i;

                    if (iarray[nn] < ibase)
                    {
                        return;
                    }

                    iarray[nn] = 0;

                    if (nn != 0)
                    {
                        iarray[nn - 1] = iarray[nn - 1] + 1;
                    }
                }
            }
        }

        public static void vec_random(int n, int base_, ref int seed, int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VEC_RANDOM selects a random N-vector of integers modulo a given base.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    21 March 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the size of the vector to be generated.
            //
            //    Input, int BASE, the base to be used.
            //
            //    Input/output, ref int SEED, a random number seed.
            //
            //    Output, int A[N], a list of N random values between
            //    0 and BASE-1.
            //
        {
            int i;

            for (i = 0; i < n; i++)
            {
                a[i] = UniformRNG.i4_uniform_ab(0, base_ - 1, ref seed);
            }

        }

        public static int[] vec_random_new(int n, int base_, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VEC_RANDOM_NEW selects a random N-vector of integers modulo a given base.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    21 March 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the size of the vector to be generated.
            //
            //    Input, int BASE, the base to be used.
            //
            //    Input/output, ref int SEED, a random number seed.
            //
            //    Output, int VEC_RANDOM_NEW[N], a list of N random values between
            //    0 and BASE-1.
            //
        {
            int[] a;
            int i;

            a = new int[n];

            for (i = 0; i < n; i++)
            {
                a[i] = UniformRNG.i4_uniform_ab(0, base_ - 1, ref seed);
            }

            return a;
        }
    }
}