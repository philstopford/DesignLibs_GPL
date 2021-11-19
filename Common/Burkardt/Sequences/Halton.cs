using System;
using Burkardt.Function;
using Burkardt.Types;

namespace Burkardt.Sequence;

public static class Halton
{
    public static double[] halton(int i, int m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HALTON computes an element of a Halton sequence.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    John Halton,
        //    On the efficiency of certain quasi-random sequences of points
        //    in evaluating multi-dimensional integrals,
        //    Numerische Mathematik,
        //    Volume 2, pages 84-90, 1960.
        //
        //  Parameters:
        //
        //    Input, int I, the index of the element of the sequence.
        //    0 <= I.
        //
        //    Input, int M, the spatial dimension.
        //
        //    Output, double HALTON[M], the element of the sequence with index I.
        //
    {
        int j;

        double[] prime_inv = new double[m];
        double[] r = new double[m];
        int[] t = new int[m];

        for (j = 0; j < m; j++)
        {
            t[j] = i;
        }

        //
        //  Carry out the computation.
        //
        for (j = 0; j < m; j++)
        {
            prime_inv[j] = 1.0 / Prime.prime(j + 1);
        }

        for (j = 0; j < m; j++)
        {
            r[j] = 0.0;
        }

        while (0 < typeMethods.i4vec_sum(m, t))
        {
            for (j = 0; j < m; j++)
            {
                int d = t[j] % Prime.prime(j + 1);
                r[j] += d * prime_inv[j];
                prime_inv[j] /= Prime.prime(j + 1);
                t[j] /= Prime.prime(j + 1);
            }
        }

        return r;
    }

    public static double[] halton_base(int i, int m, int[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HALTON_BASE computes an element of a Halton sequence with user bases.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    John Halton,
        //    On the efficiency of certain quasi-random sequences of points
        //    in evaluating multi-dimensional integrals,
        //    Numerische Mathematik,
        //    Volume 2, pages 84-90, 1960.
        //
        //  Parameters:
        //
        //    Input, int I, the index of the element of the sequence.
        //    0 <= I.
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int B[M], the bases to use for each dimension.
        //
        //    Output, double HALTON_BASE[M], the element of the sequence with index I.
        //
    {
        int j;

        double[] b_inv = new double[m];
        double[] r = new double[m];
        int[] t = new int[m];

        for (j = 0; j < m; j++)
        {
            t[j] = i;
        }

        //
        //  Carry out the computation.
        //
        for (j = 0; j < m; j++)
        {
            b_inv[j] = 1.0 / b[j];
        }

        for (j = 0; j < m; j++)
        {
            r[j] = 0.0;
        }

        while (0 < typeMethods.i4vec_sum(m, t))
        {
            for (j = 0; j < m; j++)
            {
                int d = t[j] % b[j];
                r[j] += d * b_inv[j];
                b_inv[j] /= b[j];
                t[j] /= b[j];
            }
        }

        return r;
    }

    public static int halton_inverse(double[] r, int m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HALTON_INVERSE inverts an element of the Halton sequence.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R[M], the I-th element of the Halton sequence.
        //    0 <= R < 1.0
        //
        //    Input, int M, the spatial dimension.
        //
        //    Output, int HALTON_INVERSE, the index of the element of the sequence.
        //
    {
        int j;

        for (j = 0; j < m; j++)
        {
            switch (r[j])
            {
                case < 0.0:
                case >= 1.0:
                    Console.WriteLine("");
                    Console.WriteLine("HALTON_INVERSE - Fatal error");
                    Console.WriteLine("  0 <= R < 1.0 is required.");
                    return 1;
            }
        }

        //
        //  Invert using the first component only, because working with base
        //  2 is accurate.
        //
        int i = 0;
        double t = r[0];
        int p = 1;

        while (t != 0.0)
        {
            t *= 2.0;
            int d = (int) t;
            i += d * p;
            p *= 2;
            t %= 1.0;
        }

        return i;
    }

    public static double[] halton_sequence(int i1, int i2, int m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HALTON_SEQUENCE computes elements I1 through I2 of a Halton sequence.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    John Halton,
        //    On the efficiency of certain quasi-random sequences of points
        //    in evaluating multi-dimensional integrals,
        //    Numerische Mathematik,
        //    Volume 2, pages 84-90, 1960.
        //
        //  Parameters:
        //
        //    Input, int I1, I2, the indices of the first and last 
        //    elements of the sequence.  0 <= I1, I2.
        //
        //    Input, int M, the spatial dimension.
        //    1 <= M <= 100.
        //
        //    Output, double HALTON_SEQUENCE[M*(abs(I1-I2)+1)], the elements of the 
        //    sequence with indices I1 through I2.
        //
    {
        int i;
        int i3;
        int j;
        int k;

        double[] prime_inv = new double[m];
        double[] r = new double[m * (Math.Abs(i1 - i2) + 1)];
        int[] t = new int[m];

        if (i1 <= i2)
        {
            i3 = +1;
        }
        else
        {
            i3 = -1;
        }

        int n = Math.Abs(i2 - i1) + 1;

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                r[i + j * m] = 0.0;
            }
        }

        i = i1;

        for (k = 0; k < n; k++)
        {
            for (j = 0; j < m; j++)
            {
                t[j] = i;
            }

            //
            //  Carry out the computation.
            //
            for (j = 0; j < m; j++)
            {
                prime_inv[j] = 1.0 / Prime.prime(j + 1);
            }

            while (0 < typeMethods.i4vec_sum(m, t))
            {
                for (j = 0; j < m; j++)
                {
                    int d = t[j] % Prime.prime(j + 1);
                    r[j + k * m] += d * prime_inv[j];
                    prime_inv[j] /= Prime.prime(j + 1);
                    t[j] /= Prime.prime(j + 1);
                }
            }

            i += i3;
        }

        return r;
    }

    public static bool halton_base_check(int dim_num, int[] base_)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HALTON_BASE_CHECK is TRUE if BASE is legal.
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
        //    Input, int BASE[DIM_NUM], the Halton bases.
        //    Each base must be greater than 1.
        //
        //    Output, bool HALTON_BASE_CHECK.
        //
    {
        int i;

        const bool value = true;

        for (i = 0; i < dim_num; i++)
        {
            switch (base_[i])
            {
                case <= 1:
                    Console.WriteLine("");
                    Console.WriteLine("HALTON_BASE_CHECK - Fatal error!");
                    Console.WriteLine("  Bases must be greater than 1.");
                    Console.WriteLine("  base[" + i + "] = " + base_[i] + "");
                    return false;
            }
        }

        return value;
    }

    public static double[] halton_in_circle01_accept(int dim_num, int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HALTON_IN_CIRCLE01_ACCEPT accepts Halton points in a unit circle.
        //
        //  Discussion:
        //
        //    The acceptance/rejection method is used.
        //
        //    The unit circle is centered at the origin and has radius 1.
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
        //    Input, int DIM_NUM, the dimension of the space.
        //
        //    Input, int N, the number of points.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double HALTON_IN_CIRCLE01_ACCEPT[DIM_NUM*N], the points.
        //
    {
        int i;

        int[] base_ = new int[dim_num];
        int[] leap = new int[dim_num];
        int[] seed_vec = new int[dim_num];
        double[] u = new double[dim_num];
        double[] x = new double[dim_num * n];

        int have = 0;

        for (i = 0; i < dim_num; i++)
        {
            seed_vec[i] = 0;
        }

        for (i = 0; i < dim_num; i++)
        {
            leap[i] = 1;
        }

        for (i = 0; i < dim_num; i++)
        {
            base_[i] = Prime.prime(i + 1);
        }

        while (have < n)
        {
            int step = seed;

            i4_to_halton(dim_num, step, seed_vec, leap, base_, ref u);

            seed += 1;

            double total = 0.0;
            for (i = 0; i < dim_num; i++)
            {
                u[i] = 2.0 * u[i] - 1.0;
                total += u[i] * u[i];
            }

            switch (total)
            {
                case <= 1.0:
                {
                    for (i = 0; i < dim_num; i++)
                    {
                        x[i + have * dim_num] = u[i];
                    }

                    have += 1;
                    break;
                }
            }
        }

        return x;
    }

    public static double[] halton_in_circle01_map(int dim_num, int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HALTON_IN_CIRCLE01_MAP maps Halton points into a unit circle.
        //
        //  Discussion:
        //
        //    The unit circle is centered at the origin and has radius 1.
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
        //    Input, int DIM_NUM, the dimension of the space.
        //
        //    Input, int N, the number of points.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double HALTON_IN_CIRCLE01_MAP[DIM_NUM*N], the points.
        //
    {

        int[] base_ = new int[1];
        int j;
        int[] leap = new int[1];
        int[] seed_vec = new int[1];

        double[] r = new double[n];
        double[] t = new double[n];
        double[] x = new double[dim_num * n];

        int step = 0;
        seed_vec[0] = seed;
        leap[0] = 1;
        base_[0] = Prime.prime(1);

        i4_to_halton_sequence(1, n, step, seed_vec, leap, base_, ref r);

        for (j = 0; j < n; j++)
        {
            r[j] = Math.Sqrt(r[j]);
        }

        step = 0;
        seed_vec[0] = seed;
        leap[0] = 1;
        base_[0] = Prime.prime(2);

        i4_to_halton_sequence(1, n, step, seed_vec, leap, base_, ref t);

        for (j = 0; j < n; j++)
        {
            t[j] = 2.0 * Math.PI * t[j];
        }

        for (j = 0; j < n; j++)
        {
            x[0 + j * dim_num] = r[j] * Math.Cos(t[j]);
            x[1 + j * dim_num] = r[j] * Math.Sin(t[j]);
        }

        seed += n;

        return x;
    }

    public static double[] halton_in_cube01(int dim_num, int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HALTON_IN_CUBE01 generates Halton points in the unit hypercube.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 August 2004
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
        //    Output, double HALTON_IN_CUBE01[DIM_NUM*N], the points
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

        for (i = 0; i < dim_num; i++)
        {
            base_[i] = Prime.prime(i + 1);
        }

        i4_to_halton_sequence(dim_num, n, step, seed_vec, leap, base_, ref x);

        seed += n;

        return x;
    }
        
    public static void i4_to_halton ( int seed, int[] base_, int ndim, ref double[] r, int rIndex = 0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_TO_HALTON computes an element of a Halton sequence.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 February 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    John Halton,
        //    On the efficiency of certain quasi-random sequences of points
        //    in evaluating multi-dimensional integrals,
        //    Numerische Mathematik,
        //    Volume 2, pages 84-90, 1960.
        //
        //  Parameters:
        //
        //    Input, int SEED, the index of the desired element.
        //    SEED = 0 is allowed, and returns R = 0.
        //
        //    Input, int BASE[NDIM], the Halton bases, which are usually
        //    distinct prime numbers.  Each base must be greater than 1.
        //
        //    Input, int NDIM, the dimension of the elements of the sequence.
        //
        //    Output, double R[NDIM], the SEED-th element of the Halton sequence
        //    for the given bases.
        //
    {
        int i;

        for ( i = 0; i < ndim; i++ )
        {
            switch (base_[i])
            {
                case <= 1:
                    Console.WriteLine("");
                    Console.WriteLine("I4_TO_HALTON - Fatal error!");
                    Console.WriteLine("  An input base is less than or equal to 1.");
                    Console.WriteLine("  BASE[" + i + "] = " + base_[i] + "");
                    return;
            }
        }

        for ( i = 0; i < ndim; i++ )
        {
            int seed2 = seed;
            double base_inv = 1.0 / base_[i];
            r[rIndex + i] = 0.0;

            while ( seed2 != 0 )
            {
                int digit = seed2 % base_[i];
                r[rIndex + i] += digit * base_inv;
                base_inv /= base_[i];
                seed2 /= base_[i];
            }
        }
    }

    public static void i4_to_halton(int dim_num, int step, int[] seed, int[] leap, int[] base_,
            ref double[] r)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_TO_HALTON computes one element of a leaped Halton subsequence.
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
        //  Reference:
        //
        //    J H Halton,
        //    On the efficiency of certain quasi-random sequences of points
        //    in evaluating multi-dimensional integrals,
        //    Numerische Mathematik,
        //    Volume 2, 1960, pages 84-90.
        //
        //    J H Halton and G B Smith,
        //    Algorithm 247: Radical-Inverse Quasi-Random Point Sequence,
        //    Communications of the ACM,
        //    Volume 7, 1964, pages 701-702.
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
        //    Input, int SEED[DIM_NUM], the Halton sequence index corresponding
        //    to STEP = 0.
        //    0 <= SEED(1:DIM_NUM) is required.
        //
        //    Input, int LEAP[DIM_NUM], the successive jumps in the Halton sequence.
        //    1 <= LEAP(1:DIM_NUM) is required.
        //
        //    Input, int BASE[DIM_NUM], the Halton bases.
        //    1 < BASE(1:DIM_NUM) is required.
        //
        //    Output, double R[DIM_NUM], the STEP-th element of the leaped
        //    Halton subsequence.
        //
    {
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

        if (!Halham.halham_leap_check(dim_num, leap))
        {
            return;
        }

        if (!halton_base_check(dim_num, base_))
        {
            return;
        }

        //
        //  Calculate the data.
        //
        for (i = 0; i < dim_num; i++)
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
        }
    }

    public static void i4_to_halton_sequence(int dim_num, int n, int step, int[] seed,
            int[] leap, int[] base_, ref double[] r)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_TO_HALTON_SEQUENCE computes N elements of a leaped Halton subsequence.
        //
        //  Discussion:
        //
        //    The DIM_NUM-dimensional Halton sequence is really DIM_NUM separate
        //    sequences, each generated by a particular base.
        //
        //    This routine selects elements of a "leaped" subsequence of the
        //    Halton sequence.  The subsequence elements are indexed by a
        //    quantity called STEP, which starts at 0.  The STEP-th subsequence
        //    element is simply element
        //
        //      SEED(1:DIM_NUM) + STEP * LEAP(1:DIM_NUM)
        //
        //    of the original Halton sequence.
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
        //    16 July 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    J H Halton,
        //    On the efficiency of certain quasi-random sequences of points
        //    in evaluating multi-dimensional integrals,
        //    Numerische Mathematik,
        //    Volume 2, 1960, pages 84-90.
        //
        //    J H Halton and G B Smith,
        //    Algorithm 247: Radical-Inverse Quasi-Random Point Sequence,
        //    Communications of the ACM,
        //    Volume 7, 1964, pages 701-702.
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
        //    Input, int SEED[DIM_NUM], the Halton sequence index corresponding
        //    to STEP = 0.
        //
        //    Input, int LEAP[DIM_NUM], the succesive jumps in the Halton sequence.
        //
        //    Input, int BASE[DIM_NUM], the Halton bases.
        //
        //    Output, double R[DIM_NUM*N], the next N elements of the
        //    leaped Halton subsequence, beginning with element STEP.
        //
    {
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

        if (!halton_base_check(dim_num, base_))
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
        }
    }
}