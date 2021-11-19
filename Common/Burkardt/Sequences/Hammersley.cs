using System;
using Burkardt.Function;
using Burkardt.Types;

namespace Burkardt.Sequence;

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
        int j;

        int[] t = new int[m];

        t[0] = 0;
        for (j = 1; j < m; j++)
        {
            t[j] = i;
        }

        //
        //  Carry out the computation.
        //
        double[] prime_inv = new double[m];

        prime_inv[0] = 1.0;
        for (j = 1; j < m; j++)
        {
            prime_inv[j] = 1.0 / Prime.prime(j);
        }

        double[] r = new double[m];

        r[0] = i % (n + 1) / (double) n;
        for (j = 1; j < m; j++)
        {
            r[j] = 0.0;
        }

        while (0 < typeMethods.i4vec_sum(m, t))
        {
            for (j = 1; j < m; j++)
            {
                int d = t[j] % Prime.prime(j);
                r[j] += d * prime_inv[j];
                prime_inv[j] /= Prime.prime(j);
                t[j] /= Prime.prime(j);
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
        int i;
        int j;

        for (j = 0; j < m; j++)
        {
            switch (r[j])
            {
                case < 0.0:
                case > 1.0:
                    Console.WriteLine("");
                    Console.WriteLine("HAMMERSLEY_INVERSE - Fatal error!");
                    Console.WriteLine("  0 <= R <= 1.0 is required.");
                    return 1;
            }
        }

        switch (m)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("HAMMERSLEY_INVERSE - Fatal error!");
                Console.WriteLine("  1 <= M is required.");
                return 1;
        }

        switch (n)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("HAMMERSLEY_INVERSE - Fatal error!'");
                Console.WriteLine("  1 <= N is required.");
                return 1;
        }

        switch (m)
        {
            //
            //  Invert using the second component only, because working with base
            //  2 is accurate.
            //
            case >= 2:
            {
                i = 0;
                double t = r[1];
                int p = 1;

                while (t != 0.0)
                {
                    t *= 2.0;
                    int d = (int) t;
                    
                    i += d * p;
                    p *= 2;
                    t -= d;
                }

                break;
            }
            default:
                i = (int) Math.Round(n * r[0]);
                break;
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
        int i3;
        int j;
        int k;

        if (i1 <= i2)
        {
            i3 = +1;
        }
        else
        {
            i3 = -1;
        }

        double[] prime_inv = new double[m];
        prime_inv[0] = 1.0;

        int l = Math.Abs(i2 - i1) + 1;
        double[] r = new double[m * l];

        for (k = 0; k < l; k++)
        {
            for (j = 0; j < m; j++)
            {
                r[j + k * m] = 0.0;
            }
        }

        int i = i1;

        int[] t = new int[m];
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
                prime_inv[j] = 1.0 / Prime.prime(j);
            }

            r[0 + k * m] = i % (n + 1) / (double) n;

            while (0 < typeMethods.i4vec_sum(m, t))
            {
                for (j = 1; j < m; j++)
                {
                    int d = t[j] % Prime.prime(j);
                    r[j + k * m] += d * prime_inv[j];
                    prime_inv[j] /= Prime.prime(j);
                    t[j] /= Prime.prime(j);
                }
            }

            i += i3;
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

        const bool value = true;

        for (i = 0; i < dim_num; i++)
        {
            switch (base_[i])
            {
                case 0:
                case 1:
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
        int i;

        int[] base_ = new int[dim_num];
        int[] leap = new int[dim_num];
        int[] seed_vec = new int[dim_num];
        double[] x = new double[dim_num * n];

        int step = seed;
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

        seed += n;

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
        const double FIDDLE = 0.0;

        int i;
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
            switch (base_[i])
            {
                case > 1:
                {
                    int seed2 = seed[i] + step * leap[i];

                    r[i] = 0.0;

                    double base_inv = 1.0 / base_[i];

                    while (seed2 != 0)
                    {
                        int digit = seed2 % base_[i];
                        r[i] += digit * base_inv;
                        base_inv /= base_[i];
                        seed2 /= base_[i];
                    }

                    break;
                }
                //
                default:
                    int temp = (seed[i] + step * leap[i]) % -base_[i];
                    r[i] = (temp + FIDDLE) / -base_[i];
                    break;
            }
        }
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
        const double FIDDLE = 0.0;

        int i;
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
        int[] seed2 = new int[n];

        for (i = 0; i < dim_num; i++)
        {
            int j;
            switch (base_[i])
            {
                case > 1:
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
                        double base_inv = 1.0 / base_[i];

                        while (seed2[j] != 0)
                        {
                            int digit = seed2[j] % base_[i];
                            r[i + j * dim_num] += digit * base_inv;
                            base_inv /= base_[i];
                            seed2[j] /= base_[i];
                        }
                    }

                    break;
                }
                //
                default:
                {
                    for (j = 0; j < n; j++)
                    {
                        int temp = (seed[i] + (step + j) * leap[i]) % -base_[i];

                        r[i + j * dim_num] = (temp + FIDDLE)
                                             / -base_[i];
                    }

                    break;
                }
            }
        }
    }

}