using System;
using Burkardt.PDFLib;

namespace Burkardt.DREAM;

public static class Sample
{
    public static double[] sample_candidate(int chain_index, int chain_num, double[] cr,
            int cr_index, int cr_num, int gen_index, int gen_num,
            double[] jumprate_table, int jumpstep, double[] limits, int pair_num,
            int par_num, double[] z )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SAMPLE_CANDIDATE generates candidate parameter samples.
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
        //    Input, int CHAIN_INDEX, the chain index.
        //    0 <= CHAIN_INDEX < CHAIN_NUM.
        //
        //    Input, int CHAIN_NUM, the total number of chains.
        //    3 <= CHAIN_NUM.
        //
        //    Input, double CR[CR_NUM], the CR values.
        //
        //    Input, int CR_INDEX, the index of the chosen CR value.
        //    0 <= CR_INDEX < CR_NUM.
        //
        //    Input, int CR_NUM, the total number of CR values.
        //    1 <= CR_NUM.
        //
        //    Input, int GEN_INDEX, the current generation.
        //    0 <= GEN_INDEX < GEN_NUM.
        //
        //    Input, int GEN_NUM, the total number of generations.
        //    2 <= GEN_NUM.
        //
        //    Input, double JUMPRATE_TABLE[PAR_NUM], the jumprate table.
        //
        //    Input, int JUMPSTEP, forces a "long jump" every
        //    JUMPSTEP generations.
        //
        //    Input, double LIMITS[2*PAR_NUM], limits for the parameters.
        //
        //    Input, int PAIR_NUM, the number of pairs of 
        //    crossover chains.
        //    0 <= PAIR_NUM.
        //
        //    Input, int PAR_NUM, the total number of parameters.
        //    1 <= PAR_NUM.
        //
        //    Input, double Z[PAR_NUM*CHAIN_NUM*GEN_NUM], the Markov chain 
        //    sample data.
        //
        //    Output, double SAMPLE_CANDIDATE[PAR_NUM], a candidate parameter sample.
        //
        //  Local parameters:
        //
        //    Input, int JUMP_DIM[JUMP_NUM], the dimensions in which
        //    a jump is to be made.
        //
        //    Local, int JUMP_NUM, the number of dimensions in which
        //    a jump will be made.  0 <= JUMP_NUM <= PAR_NUM.
        //
        //    Local, double JUMPRATE, the jump rate.
        //
    {
        int i;
        int jump_num = 0;
        double jumprate = 0;
        int[] pair = new int[2];
        //
        //  Used to calculate E following a uniform distribution on (-B,+B).
        //  Because B is currently zero, the noise term is suppressed.
        //
        double b = 0.0;
        //
        //  Pick pairs of other chains for crossover.
        //
        int[] r = new int[2 * pair_num];

        for (i = 0; i < pair_num; i++)
        {
            while (true)
            {
                double r2 = PDF.r8_uniform_01_sample();
                pair[0] = (int) (r2 * chain_num);
                r2 = PDF.r8_uniform_01_sample();
                pair[1] = (int) (r2 * chain_num);

                if (pair[0] != pair[1] &&
                    pair[0] != chain_index &&
                    pair[1] != chain_index)
                {
                    break;
                }
            }

            r[0 + i * 2] = pair[0];
            r[1 + i * 2] = pair[1];
        }

        //
        //  Determine the jump rate.
        //
        int[] jump_dim = new int[par_num];

        Jump.jumprate_choose(cr, cr_index, cr_num, gen_index, jump_dim, ref jump_num,
            ref jumprate, jumprate_table, jumpstep, par_num);
        //
        //  Calculate E in equation 4 of Vrugt.
        //
        double[] noise_e = new double[par_num];

        for (i = 0; i < par_num; i++)
        {
            noise_e[i] = b * (2.0 * PDF.r8_uniform_01_sample() - 1.0);
        }

        //
        //  Get epsilon value from multinormal distribution                      
        //
        double[] eps = new double[par_num];

        double av = 0.0;
        double sd = 1.0E-10;
        for (i = 0; i < par_num; i++)
        {
            eps[i] = PDF.r8_normal_sample(av, sd);
        }

        //
        //  Generate the candidate sample ZP based on equation 4 of Vrugt.
        //
        double[] diff = Diff.diff_compute(chain_num, gen_index, gen_num, jump_dim, jump_num,
            pair_num, par_num, r, z);

        double[] zp = new double[par_num];

        for (i = 0; i < par_num; i++)
        {
            zp[i] = z[i + chain_index * par_num + (gen_index - 1) * par_num * chain_num];
        }

        for (i = 0; i < par_num; i++)
        {
            zp[i] = zp[i] + (1.0 + noise_e[i]) * jumprate * diff[i] + eps[i];
        }

        //
        //  Enforce limits on the sample ZP.
        //
        sample_limits(limits, par_num, ref zp);
        return zp;
    }

    public static void sample_limits(double[] limits, int par_num, ref double[] zp )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SAMPLE_LIMITS enforces limits on a sample variable.
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
        //  Parameters:
        //
        //    Input, double LIMITS[2*PAR_NUM], the parameter limits.
        //
        //    Input, int PAR_NUM, the total number of parameters.
        //    1 <= PAR_NUM.
        //
        //    Input/output, double ZP[PAR_NUM], a variable, whose entries,
        //    if necessary, will be "folded" so that they lie within the limits.
        //
    {
        int i;

        for (i = 0; i < par_num; i++)
        {
            double w = limits[1 + i * 2] - limits[0 + i * 2];

            switch (w)
            {
                case 0.0:
                    zp[i] = limits[0 + i * 2];
                    break;
                case < 0.0:
                    Console.WriteLine("");
                    Console.WriteLine("SAMPLE_LIMITS - Fatal error!");
                    Console.WriteLine("  Upper limit less than lower limit.");
                    return;
                default:
                {
                    while (zp[i] < limits[0 + i * 2])
                    {
                        zp[i] += w;
                    }

                    while (limits[1 + i * 2] < zp[i])
                    {
                        zp[i] -= w;
                    }

                    break;
                }
            }
        }
    }
}