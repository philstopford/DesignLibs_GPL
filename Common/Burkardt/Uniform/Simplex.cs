using System;

namespace Burkardt.Uniform;

public static class Simplex
{
    public static double[] uniform_in_simplex01_map(int dim_num, int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_IN_SIMPLEX01 maps uniform points into the unit simplex.
        //
        //  Discussion:
        //
        //    The interior of the unit DIM_NUM dimensional simplex is the set of points X(1:DIM_NUM)
        //    such that each X(I) is nonnegative, and sum(X(1:DIM_NUM)) <= 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 August 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Reuven Rubinstein,
        //    Monte Carlo Optimization, Simulation, and Sensitivity
        //    of Queueing Networks,
        //    Krieger, 1992,
        //    ISBN: 0894647644,
        //    LC: QA298.R79.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the dimension of the space.
        //
        //    Input, int N, the number of points.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double UNIFORM_IN_SIMPLEX01_MAP[DIM_NUM*N], the points.
        //
    {
        int j;
        //
        //  The construction begins by sampling DIM_NUM+1 points from the
        //  exponential distribution with parameter 1.
        //
        double[] e = new double[dim_num + 1];
        double[] x = new double[dim_num * n];

        for (j = 0; j < n; j++)
        {
            UniformRNG.r8vec_uniform_01(dim_num + 1, ref seed, ref e);

            int i;
            for (i = 0; i <= dim_num; i++)
            {
                e[i] = -Math.Log(e[i]);
            }

            double total = 0.0;
            for (i = 0; i <= dim_num; i++)
            {
                total += e[i];
            }

            for (i = 0; i < dim_num; i++)
            {
                x[i + dim_num * j] = e[i] / total;
            }

        }

        return x;
    }

    public static double[] uniform_on_simplex01_map(int dim_num, int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_ON_SIMPLEX01_MAP maps uniform points onto the unit simplex.
        //
        //  Discussion:
        //
        //    The surface of the unit DIM_NUM-dimensional simplex is the set of points
        //    X(1:DIM_NUM) such that each X(I) is nonnegative,
        //    every X(I) is no greater than 1, and
        //
        //    ( X(I) = 0 for some I, or sum ( X(1:DIM_NUM) ) = 1. )
        //
        //    In DIM_NUM dimensions, there are DIM_NUM sides, and one main face.
        //    This code picks a point uniformly with respect to "area".
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 August 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Reuven Rubinstein,
        //    Monte Carlo Optimization, Simulation, and Sensitivity
        //    of Queueing Networks,
        //    Krieger, 1992,
        //    ISBN: 0894647644,
        //    LC: QA298.R79.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the dimension of the space.
        //
        //    Input, int N, the number of points.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double UNIFORM_ON_SIMPLEX01_MAP[DIM_NUM*N], the points.
        //
    {
        int j;
        //
        //  The construction begins by sampling DIM_NUM points from the
        //  exponential distribution with parameter 1.
        //
        double[] e = new double[dim_num];
        double[] x = new double[dim_num * n];

        for (j = 0; j < n; j++)
        {
            UniformRNG.r8vec_uniform_01(dim_num, ref seed, ref e);

            int i;
            for (i = 0; i < dim_num; i++)
            {
                e[i] = -Math.Log(e[i]);
            }

            double total = 0.0;
            for (i = 0; i < dim_num; i++)
            {
                total += e[i];
            }

            //
            //  Based on their relative areas, choose a side of the simplex,
            //  or the main face.
            //
            for (i = 0; i < dim_num; i++)
            {
                x[i + j * dim_num] = e[i] / total;
            }

            double area1 = Math.Sqrt(dim_num);
            double area2 = dim_num;

            double r = UniformRNG.r8_uniform_01(ref seed);

            if (area1 / (area1 + area2) < r)
            {
                i = UniformRNG.i4_uniform_ab(0, dim_num - 1, ref seed);
                x[i + j * dim_num] = 0.0;
            }

        }

        return x;
    }


}