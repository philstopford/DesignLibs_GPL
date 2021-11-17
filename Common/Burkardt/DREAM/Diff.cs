using Burkardt.Types;

namespace Burkardt.DREAM;

public static class Diff
{
    public static double[] diff_compute(int chain_num, int gen_index, int gen_num,
            int[] jump_dim, int jump_num, int pair_num, int par_num, int[] r,
            double[] z )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIFF_COMPUTE computes the differential evolution.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 May 2013
        //
        //  Author:
        //
        //    Original FORTRAN90 version by Guannan Zhang.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Jasper Vrugt, CJF ter Braak, CGH Diks, Bruce Robinson, James Hyman, 
        //    Dave Higdon,
        //    Accelerating Markov Chain Monte Carlo Simulation by Differential 
        //    Evolution with Self-Adaptive Randomized Subspace Sampling,
        //    International Journal of Nonlinear Sciences and Numerical Simulation,
        //    Volume 10, Number 3, March 2009, pages 271-288.
        //
        //  Parameters:
        //
        //    Input, int CHAIN_NUM, the total number of chains.
        //    3 <= CHAIN_NUM.
        //
        //    Input, int GEN_INDEX, the index of the current generation.
        //    1 <= GEN_INDEX <= GEN_NUM.
        //
        //    Input, int GEN_NUM, the total number of generations.
        //    2 <= GEN_NUM.
        //
        //    Input, int JUMP_DIM[JUMP_NUM], the dimensions in which
        //    a jump is to be made.
        //
        //    Input, int JUMP_NUM, the number of dimensions in which
        //    a jump will be made.  0 <= JUMP_NUM <= PAR_NUM.
        //
        //    Input, int PAIR_NUM, the number of pairs of 
        //    crossover chains.
        //    0 <= PAIR_NUM.
        //
        //    Input, int PAR_NUM, the total number of parameters.
        //    1 <= PAR_NUM.
        //
        //    Input, int R[2*PAIR_NUM], pairs of chains used
        //    to compute differences.
        //
        //    Input, double Z[PAR_NUM*CHAIN_NUM*GEN_NUM], the Markov chain 
        //    sample data.
        //
        //    Output, double DIFF_COMPUTE[PAR_NUM], the vector of pair differences.
        //
    {
        double[] diff;
        int i1;
        int i2;
        int j;
        int k;
        int pair;
        int r1;
        int r2;
        //
        //  Produce the difference of the pairs used for population evolution.
        //
        diff = typeMethods.r8vec_zero_new(par_num);

        for (pair = 0; pair < pair_num; pair++)
        {
            r1 = r[0 + pair * 2];
            r2 = r[1 + pair * 2];
            for (j = 0; j < jump_num; j++)
            {
                k = jump_dim[j];
                i1 = k + r1 * par_num + (gen_index - 1) * par_num * chain_num;
                i2 = k + r2 * par_num + (gen_index - 1) * par_num * chain_num;
                diff[k] += (z[i1] - z[i2]);
            }
        }

        return diff;
    }
}