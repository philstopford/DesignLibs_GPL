using System;
using Burkardt.Types;

namespace Burkardt
{
    public static class Hammersley
    {
        public static double[] hammersley(int i, int m, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HAMMERSLEY computes an element of a Hammersley sequence.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 August 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    John Hammersley,
            //    Monte Carlo methods for solving multivariable problems,
            //    Proceedings of the New York Academy of Science,
            //    Volume 86, 1960, pages 844-874.
            //
            //  Parameters:
            //
            //    Input, int I, the index of the element of the sequence.
            //    0 <= I.
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, int N, the "base" for the first component.
            //    1 <= N.
            //
            //    Output, double HAMMERSLEY[M], the element of the sequence with index I.
            //
        {
            int d;
            int j;
            double[] prime_inv;
            double[] r;
            int[] t;

            t = new int[m];

            t[0] = 0;
            for (j = 1; j < m; j++)
            {
                t[j] = i;
            }

            //
            //  Carry out the computation.
            //
            prime_inv = new double[m];

            prime_inv[0] = 1.0;
            for (j = 1; j < m; j++)
            {
                prime_inv[j] = 1.0 / (double) (Prime.prime(j));
            }

            r = new double[m];

            r[0] = (double) (i % (n + 1)) / (double) (n);
            for (j = 1; j < m; j++)
            {
                r[j] = 0.0;
            }

            while (0 < typeMethods.i4vec_sum(m, t))
            {
                for (j = 1; j < m; j++)
                {
                    d = (t[j] % Prime.prime(j));
                    r[j] = r[j] + (double) (d) * prime_inv[j];
                    prime_inv[j] = prime_inv[j] / (double) (Prime.prime(j));
                    t[j] = (t[j] / Prime.prime(j));
                }
            }

            return r;
        }

        public static int hammersley_inverse(double[] r, int m, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HAMMERSLEY_INVERSE inverts an element of the Hammersley sequence.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 August 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double R[M], the I-th element of the Hammersley sequence.
            //    0 <= R < 1.0
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, int N, the "base" for the first component.
            //    1 <= N.
            //
            //    Output, int HAMMERSLEY_INVERSE, the index of the element of the sequence.
            //
        {
            int d;
            int i;
            int j;
            int p;
            double t;

            for (j = 0; j < m; j++)
            {
                if (r[j] < 0.0 || 1.0 < r[j])
                {
                    Console.WriteLine("");
                    Console.WriteLine("HAMMERSLEY_INVERSE - Fatal error!");
                    Console.WriteLine("  0 <= R <= 1.0 is required.");
                    return (1);
                }
            }

            if (m < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("HAMMERSLEY_INVERSE - Fatal error!");
                Console.WriteLine("  1 <= M is required.");
                return (1);
            }

            if (n < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("HAMMERSLEY_INVERSE - Fatal error!'");
                Console.WriteLine("  1 <= N is required.");
                return (1);
            }

            //
            //  Invert using the second component only, because working with base
            //  2 is accurate.
            //
            if (2 <= m)
            {
                i = 0;
                t = r[1];
                p = 1;

                while (t != 0.0)
                {
                    t = t * 2.0;
                    d = (int) (t);
                    i = i + d * p;
                    p = p * 2;
                    t = t - (double) (d);
                }
            }
            else
            {
                i = (int) (Math.Round((double) (n) * r[0]));
            }

            return i;
        }

        public static double[] hammersley_sequence(int i1, int i2, int m, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HAMMERSLEY_SEQUENCE computes elements I1 through I2 of a Hammersley sequence.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 August 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    John Hammersley,
            //    Monte Carlo methods for solving multivariable problems,
            //    Proceedings of the New York Academy of Science,
            //    Volume 86, 1960, pages 844-874.
            //
            //  Parameters:
            //
            //    Input, int I1, I2, the indices of the first and last 
            //    elements of the sequence.  0 <= I1, I2.
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, int N, the "base" for the first component.
            //    1 <= N.
            //
            //    Output, double HAMMERSLEY_SEQUENCE[M*(abs(I1-I2)+1)], the elements of 
            //    the sequence with indices I1 through I2.
            //
        {
            int d;
            int i;
            int i3;
            int j;
            int k;
            int l;
            double[] prime_inv;
            double[] r;
            int[] t;

            if (i1 <= i2)
            {
                i3 = +1;
            }
            else
            {
                i3 = -1;
            }

            prime_inv = new double[m];
            prime_inv[0] = 1.0;

            l = Math.Abs(i2 - i1) + 1;
            r = new double[m * l];

            for (k = 0; k < l; k++)
            {
                for (j = 0; j < m; j++)
                {
                    r[j + k * m] = 0.0;
                }
            }

            i = i1;

            t = new int[m];
            for (k = 0; k < l; k++)
            {
                t[0] = 0;
                for (j = 1; j < m; j++)
                {
                    t[j] = i;
                }

                //
                //  Carry out the computation.
                //
                for (j = 1; j < m; j++)
                {
                    prime_inv[j] = 1.0 / (double) (Prime.prime(j));
                }

                r[0 + k * m] = (double) (i % (n + 1)) / (double) (n);

                while (0 < typeMethods.i4vec_sum(m, t))
                {
                    for (j = 1; j < m; j++)
                    {
                        d = (t[j] % Prime.prime(j));
                        r[j + k * m] = r[j + k * m] + (double) (d) * prime_inv[j];
                        prime_inv[j] = prime_inv[j] / (double) (Prime.prime(j));
                        t[j] = (t[j] / Prime.prime(j));
                    }
                }

                i = i + i3;
            }

            return r;
        }

        public static bool hammersley_base_check(int dim_num, int[] base_)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HAMMERSLEY_BASE_CHECK is TRUE if BASE is legal.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 July 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, int BASE[DIM_NUM], the bases.
            //
            //    Output, bool HAMMERSLEY_BASE_CHECK.
        {
            int i;
            bool value;

            value = true;

            for (i = 0; i < dim_num; i++)
            {
                if (base_[i] == 0 || base_[i] == 1)
                {
                    Console.WriteLine("");
                    Console.WriteLine("HAMMERSLEY_BASE_CHECK - Fatal error!");
                    Console.WriteLine("  Bases may not be 0 or 1.");
                    typeMethods.i4vec_transpose_print(dim_num, base_, "BASE:  ");
                    return false;
                }
            }

            return value;
        }

        public static double[] hammersley_in_cube01(int dim_num, int n, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HAMMERSLEY_IN_CUBE01 computes Hammersley points in the unit hypercube.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    25 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, int N, the number of elements.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, double HAMMERSLEY_IN_CUBE01[DIM_NUM*N], the points.
            //
        {
            int[] base_;
            int i;
            int[] leap;
            int[] seed_vec;
            int step;
            double[] x;

            base_ = new int[dim_num];
            leap = new int[dim_num];
            seed_vec = new int[dim_num];
            x = new double[dim_num * n];

            step = seed;
            for (i = 0; i < dim_num; i++)
            {
                seed_vec[i] = 0;
            }

            for (i = 0; i < dim_num; i++)
            {
                leap[i] = 1;
            }

            base_[0] = -n;
            for (i = 1; i < dim_num; i++)
            {
                base_[i] = Prime.prime(i);
            }

            i4_to_hammersley_sequence(dim_num, n, step, seed_vec, leap, base_, ref x);

            seed = seed + n;

            return x;
        }

        public static void i4_to_hammersley(int dim_num, int step, int[] seed, int[] leap,
                int[] base_, ref double[] r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4_TO_HAMMERSLEY computes one element of a leaped Hammersley subsequence.
            //
            //  Discussion:
            //
            //    The DIM_NUM-dimensional Hammersley sequence is really DIM_NUM separate
            //    sequences, each generated by a particular base.  If the base is
            //    greater than 1, a standard 1-dimensional
            //    van der Corput sequence is generated.  But if the base is
            //    negative, this is a signal that the much simpler sequence J/(-BASE)
            //    is to be generated.  For the standard Hammersley sequence, the
            //    first spatial coordinate uses a base of (-N), and subsequent
            //    coordinates use bases of successive primes (2, 3, 5, 7, 11, ...).
            //    This program allows the user to specify any combination of bases,
            //    included nonprimes and repeated values.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 September 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    J M Hammersley,
            //    Monte Carlo methods for solving multivariable problems,
            //    Proceedings of the New York Academy of Science,
            //    Volume 86, 1960, pages 844-874.
            //
            //    Ladislav Kocis and William Whiten,
            //    Computational Investigations of Low-Discrepancy Sequences,
            //    ACM Transactions on Mathematical Software,
            //    Volume 23, Number 2, 1997, pages 266-294.
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //    1 <= DIM_NUM is required.
            //
            //    Input, int STEP, the index of the subsequence element.
            //    0 <= STEP is required.
            //
            //    Input, int SEED[DIM_NUM], the Hammersley sequence index corresponding
            //    to STEP = 0.
            //    0 <= SEED(1:DIM_NUM) is required.
            //
            //    Input, int LEAP[DIM_NUM], the successive jumps in the Hammersley sequence.
            //    1 <= LEAP(1:DIM_NUM) is required.
            //
            //    Input, int BASE[DIM_NUM], the Hammersley bases.
            //
            //    Output, double R[DIM_NUM], the STEP-th element of the leaped
            //    Hammersley subsequence.
            //
        {
            double FIDDLE = 0.0;

            double base_inv;
            int digit;
            int i;
            int seed2;
            int temp;
            //
            //  Check the input.
            //
            if (!Halham.halham_dim_num_check(dim_num))
            {
                return;
            }

            if (!Halham.halham_step_check(step))
            {
                return;
            }

            if (!Halham.halham_seed_check(dim_num, seed))
            {
                return;
            }

            if (!Halham.halham_leap_check(dim_num, seed))
            {
                return;
            }

            if (!hammersley_base_check(dim_num, base_))
            {
                return;
            }

            //
            //  Calculate the data.
            //
            for (i = 0; i < dim_num; i++)
            {
                if (1 < base_[i])
                {
                    seed2 = seed[i] + step * leap[i];

                    r[i] = 0.0;

                    base_inv = 1.0 / ((double) base_[i]);

                    while (seed2 != 0)
                    {
                        digit = seed2 % base_[i];
                        r[i] = r[i] + ((double) digit) * base_inv;
                        base_inv = base_inv / ((double) base_[i]);
                        seed2 = seed2 / base_[i];
                    }
                }
                //
                //  In the following computation, the value of FIDDLE can be:
                //
                //    0,   for the sequence 0/N, 1/N, ..., N-1/N
                //    1,   for the sequence 1/N, 2/N, ..., N/N
                //    1/2, for the sequence 1/(2N), 3/(2N), ..., (2*N-1)/(2N)
                //
                else
                {
                    temp = (seed[i] + step * leap[i]) % (-base_[i]);
                    r[i] = ((double) (temp) + FIDDLE) / (double) (-base_[i]);
                }
            }

            return;
        }

        public static void i4_to_hammersley_sequence(int dim_num, int n, int step, int[] seed,
                int[] leap, int[] base_, ref double[] r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4_TO_HAMMERSLEY_SEQUENCE computes N elements of a leaped Hammersley subsequence.
            //
            //  Discussion:
            //
            //    The DIM_NUM-dimensional Hammersley sequence is really DIM_NUM separate
            //    sequences, each generated by a particular base.  If the base is
            //    greater than 1, a standard 1-dimensional
            //    van der Corput sequence is generated.  But if the base is
            //    negative, this is a signal that the much simpler sequence J/(-BASE)
            //    is to be generated.  For the standard Hammersley sequence, the
            //    first spatial coordinate uses a base of (-N), and subsequent
            //    coordinates use bases of successive primes (2, 3, 5, 7, 11, ...).
            //    This program allows the user to specify any combination of bases,
            //    included nonprimes and repeated values.
            //
            //    This routine selects elements of a "leaped" subsequence of the
            //    Hammersley sequence.  The subsequence elements are indexed by a
            //    quantity called STEP, which starts at 0.  The STEP-th subsequence
            //    element is simply element
            //
            //      SEED(1:DIM_NUM) + STEP * LEAP(1:DIM_NUM)
            //
            //    of the original Hammersley sequence.
            //
            //
            //    The data to be computed has two dimensions.
            //
            //    The number of data items is DIM_NUM * N, where DIM_NUM is the spatial dimension
            //    of each element of the sequence, and N is the number of elements of the sequence.
            //
            //    The data is stored in a one dimensional array R.  The first element of the
            //    sequence is stored in the first DIM_NUM entries of R, followed by the DIM_NUM entries
            //    of the second element, and so on.
            //
            //    In particular, the J-th element of the sequence is stored in entries
            //    0+(J-1)*DIM_NUM through (DIM_NUM-1) + (J-1)*DIM_NUM.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 September 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    J M Hammersley,
            //    Monte Carlo methods for solving multivariable problems,
            //    Proceedings of the New York Academy of Science,
            //    Volume 86, 1960, pages 844-874.
            //
            //    Ladislav Kocis and William Whiten,
            //    Computational Investigations of Low-Discrepancy Sequences,
            //    ACM Transactions on Mathematical Software,
            //    Volume 23, Number 2, 1997, pages 266-294.
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, int N, the number of elements of the sequence.
            //
            //    Input, int STEP, the index of the subsequence element.
            //    0 <= STEP is required
            //
            //    Input, int SEED[DIM_NUM], the Hammersley sequence index corresponding
            //    to STEP = 0.
            //
            //    Input, int LEAP[DIM_NUM], the succesive jumps in the Hammersley sequence.
            //
            //    Input, int BASE[DIM_NUM], the Hammersley bases.
            //
            //    Output, double R[DIM_NUM*N], the next N elements of the
            //    leaped Hammersley subsequence, beginning with element STEP.
            //
        {
            double FIDDLE = 0.0;

            double base_inv;
            int digit;
            int i;
            int j;
            int[] seed2;
            int temp;
            //
            //  Check the input.
            //
            if (!Halham.halham_dim_num_check(dim_num))
            {
                return;
            }

            if (!Halham.halham_n_check(n))
            {
                return;
            }

            if (!Halham.halham_step_check(step))
            {
                return;
            }

            if (!Halham.halham_seed_check(dim_num, seed))
            {
                return;
            }

            if (!Halham.halham_leap_check(dim_num, leap))
            {
                return;
            }

            if (!hammersley_base_check(dim_num, base_))
            {
                return;
            }

            //
            //  Calculate the data.
            //
            seed2 = new int[n];

            for (i = 0; i < dim_num; i++)
            {
                if (1 < base_[i])
                {
                    for (j = 0; j < n; j++)
                    {
                        seed2[j] = seed[i] + (step + j) * leap[i];
                    }

                    for (j = 0; j < n; j++)
                    {
                        r[i + j * dim_num] = 0.0;
                    }

                    for (j = 0; j < n; j++)
                    {
                        base_inv = 1.0 / ((double) base_[i]);

                        while (seed2[j] != 0)
                        {
                            digit = seed2[j] % base_[i];
                            r[i + j * dim_num] = r[i + j * dim_num] + ((double) digit) * base_inv;
                            base_inv = base_inv / ((double) base_[i]);
                            seed2[j] = seed2[j] / base_[i];
                        }
                    }
                }
                //
                //  In the following computation, the value of FIDDLE can be:
                //
                //    0,   for the sequence 0/N, 1/N, ..., N-1/N
                //    1,   for the sequence 1/N, 2/N, ..., N/N
                //    1/2, for the sequence 1/(2N), 3/(2N), ..., (2*N-1)/(2N)
                //
                else
                {
                    for (j = 0; j < n; j++)
                    {
                        temp = (seed[i] + (step + j) * leap[i]) % (-base_[i]);

                        r[i + j * dim_num] = ((double) (temp) + FIDDLE)
                                             / (double) (-base_[i]);
                    }
                }
            }
        }

    }
}