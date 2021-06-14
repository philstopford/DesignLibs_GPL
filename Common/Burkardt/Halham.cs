using System;

namespace Burkardt
{
    public static class Halham
    {
        public static bool halham_dim_num_check(int dim_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HALHAM_DIM_NUM_CHECK checks DIM_NUM for a Halton or Hammersley sequence.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 September 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //    DIM_NUM must be positive.
            //
            //    Output, bool HALHAM_DIM_NUM_CHECK, is true if DIM_NUM is legal.
            //
        {
            bool value;

            if (dim_num < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("HALHAM_DIM_NUM_CHECK - Fatal error!");
                Console.WriteLine("  DIM_NUM < 0." + "  DIM_NUM = " + dim_num + "");
                return false;
            }

            value = true;

            return value;
        }

        public static bool halham_leap_check(int dim_num, int[] leap)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HALHAM_LEAP_CHECK checks LEAP for a Halton or Hammersley sequence.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 September 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, int LEAP[DIM_NUM], the successive jumps in the sequence.
            //    Each entry must be greater than 0.
            //
            //    Output, bool HALHAM_LEAP_CHECK, is true if LEAP is legal.
            //
        {
            int i;
            bool value;

            value = true;

            for (i = 0; i < dim_num; i++)
            {
                if (leap[i] < 1)
                {
                    Console.WriteLine("");
                    Console.WriteLine("HALHAM_LEAP_CHECK - Fatal error!");
                    Console.WriteLine("  Leap entries must be greater than 0.");
                    Console.WriteLine("  leap[" + i + "] = " + leap[i] + "");
                    return false;
                }
            }

            return value;
        }

        public static bool halham_n_check(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HALHAM_N_CHECK checks N for a Halton or Hammersley sequence.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 September 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of elements of the subsequence.
            //    N must be positive.
            //
            //    Output, bool HALHAM_N_CHECK, is true if N is legal.
            //
        {
            bool value;

            if (n < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("HALHAM_N_CHECK - Fatal error!");
                Console.WriteLine("  N < 0." + "  N = " + n + "");
                return false;
            }

            value = true;

            return value;
        }

        public static bool halham_seed_check(int dim_num, int[] seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HALHAM_SEED_CHECK checks SEED for a Halton or Hammersley sequence.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 September 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, int SEED[DIM_NUM], the sequence index
            //    corresponding to STEP = 0.  Each entry must be 0 or greater.
            //
            //    Output, bool HALHAM_SEED_CHECK, is true if SEED is legal.
            //
        {
            int i;
            bool value;

            value = true;

            for (i = 0; i < dim_num; i++)
            {
                if (seed[i] < 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("HALHAM_SEED_CHECK - Fatal error!");
                    Console.WriteLine("  SEED entries must be nonnegative.");
                    Console.WriteLine("  seed[" + i + "] = " + seed[i] + "");
                    return false;
                }
            }

            return value;
        }

        public static bool halham_step_check(int step)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HALHAM_STEP_CHECK checks STEP for a Halton or Hammersley sequence.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 September 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int STEP, the index of the subsequence element.
            //    STEP must be 1 or greater.
            //
            //    Output, bool HALHAM_STEP_CHECK, is true if STEP is legal.
            //
        {
            bool value;

            if (step < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("HALHAM_STEP_CHECK - Fatal error!");
                Console.WriteLine("  STEP < 0." + "  STEP = " + step + "");
                return false;
            }

            value = true;

            return value;
        }
    }
}