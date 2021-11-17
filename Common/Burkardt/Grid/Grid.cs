using System;
using Burkardt.RandomNS;
using Burkardt.SubsetNS;
using Burkardt.Types;
using Tuple = System.Tuple;

namespace Burkardt.Grid;

public static class Grid
{
    public static double[] grid_in_cube01(int dim_num, int n, int center, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRID_IN_CUBE01 generates a grid dataset in the unit hypercube.
        //
        //  Discussion:
        //
        //    N points are needed in a DIM_NUM dimensional space.
        //
        //    The points are to lie on a uniform grid of side N_SIDE.
        //
        //    Unless the N = N_SIDE^DIM_NUM for some N_SIDE, we can't use all the
        //    points on a grid.  What we do is find the smallest N_SIDE
        //    that's big enough, and randomly omit some points.
        //
        //    If N_SIDE is 4, then the choices in 1D are:
        //
        //    A: 0,   1/3, 2/3, 1
        //    B: 1/5, 2/5, 3/5, 4/5
        //    C: 0,   1/4, 2/4, 3/4
        //    D: 1/4, 2/4, 3/4, 1
        //    E: 1/8, 3/8, 5/8, 7/8
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
        //    Input, int N, the number of points.
        //
        //    Input, int CENTER, specifies the 1D grid centering:
        //    1: first point is 0.0, last point is 1.0;
        //    2: first point is 1/(N+1), last point is N/(N+1);
        //    3: first point is 0, last point is (N-1)/N;
        //    4: first point is 1/N, last point is 1;
        //    5: first point is 1/(2*N), last point is (2*N-1)/(2*N);
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double GRID_IN_CUBE01[DIM_NUM*N], the points.
        //
    {
        int i;
        int j;
        int n_grid;
        int n_side;
        double[] r;
        int rank;
        int[] rank_list;
        int[] tuple;
        //
        //  Find the dimension of the smallest grid with N points.
        //
        n_side = grid_side(dim_num, n);
        //
        //  We need to select N points out of N_SIDE^DIM_NUM set.
        //
        n_grid = (int) Math.Pow(n_side, dim_num);
        //
        //  Generate a random subset of N items from a set of size N_GRID.
        //
        rank_list = new int[n];

        Ksub.ksub_random2(n_grid, n, ref seed, ref rank_list);
        //
        //  Must make one dummy call to TUPLE_NEXT_FAST with RANK = -1.
        //
        rank = -1;
        tuple = new int[dim_num];
        BTupleData data = new() {base_ = new int[dim_num]};
        BTuple.tuple_next_fast(ref data, n_side, dim_num, rank, ref tuple);
        //
        //  Now generate the appropriate indices, and "center" them.
        //
        r = new double[dim_num * n];

        for (j = 0; j < n; j++)
        {
            rank = rank_list[j] - 1;

            BTuple.tuple_next_fast(ref data, n_side, dim_num, rank, ref tuple);

            switch (center)
            {
                case 1:
                {
                    for (i = 0; i < dim_num; i++)
                    {
                        r[i + j * dim_num] = (tuple[i] - 1)
                                             / (double) (n_side - 1);
                    }

                    break;
                }
                case 2:
                {
                    for (i = 0; i < dim_num; i++)
                    {
                        r[i + j * dim_num] = tuple[i]
                                             / (double) (n_side + 1);
                    }

                    break;
                }
                case 3:
                {
                    for (i = 0; i < dim_num; i++)
                    {
                        r[i + j * dim_num] = (tuple[i] - 1)
                                             / (double) n_side;
                    }

                    break;
                }
                case 4:
                {
                    for (i = 0; i < dim_num; i++)
                    {
                        r[i + j * dim_num] = tuple[i]
                                             / (double) n_side;
                    }

                    break;
                }
                case 5:
                {
                    for (i = 0; i < dim_num; i++)
                    {
                        r[i + j * dim_num] = (2 * tuple[i] - 1)
                                             / (double) (2 * n_side);
                    }

                    break;
                }
            }
        }

        return r;
    }

    public static int grid_side(int dim_num, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRID_SIDE finds the smallest DIM_NUM dimensional grid containing at least N points.
        //
        //  Discussion:
        //
        //    Each coordinate of the grid will have N_SIDE distinct values.
        //    Thus the total number of points in the grid is N_SIDE**DIM_NUM.
        //    This routine seeks the smallest N_SIDE such that N <= N_SIDE**DIM_NUM.
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
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int N, the number of points.
        //
        //    Output, int GRID_SIDE, the length of one side of the smallest
        //    grid in DIM_NUM dimensions that contains at least N points.
        //
    {
        double exponent;
        int n_side;

        switch (n)
        {
            case <= 0:
                n_side = 0;
                return n_side;
        }

        switch (dim_num)
        {
            case <= 0:
                n_side = -1;
                return n_side;
        }

        exponent = 1.0 / dim_num;

        n_side = (int) Math.Pow(n, exponent);

        if (Math.Pow(n_side, dim_num) < n)
        {
            n_side += 1;
        }

        return n_side;
    }

}