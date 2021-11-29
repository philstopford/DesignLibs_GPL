using System;
using Burkardt.PDFLib;
using Burkardt.Types;

namespace Burkardt.DREAM;

public static class CR
{
    public static void cr_dis_update(int chain_index, int chain_num, ref double[] cr_dis,
            int cr_index, int cr_num, ref int[] cr_ups, int gen_index, int gen_num,
            int par_num, double[] z )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CR_DIS_UPDATE updates the CR distance.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 August 2018
        //
        //  Author:
        //
        //    Original FORTRAN90 version by Guannan Zhang.
        //    C++ version by John Burkardt.
        //
        //  Parameters:
        //
        //    Input, int CHAIN_INDEX, the index of the chain.
        //    0 <= CHAIN_INDEX < CHAIN_NUM.
        //
        //    Input, int CHAIN_NUM, the total number of chains.
        //    3 <= CHAIN_NUM.
        //
        //    Input/output, double CR_DIS[CR_NUM], the CR distances.
        //
        //    Input, int CR_INDEX, the index of the CR.
        //    0 <= CR_INDEX < CR_NUM.
        //
        //    Input, int CR_NUM, the total number of CR values.
        //    1 <= CR_NUM.
        //
        //    Input/output, int CR_UPS[CR_NUM], the number of updates 
        //    for each CR.
        //
        //    Input, int GEN_INDEX, the index of the generation.
        //    0 <= GEN_INDEX < GEN_NUM.
        //
        //    Input, int GEN_NUM, the total number of generations.
        //    2 <= GEN_NUM.
        //
        //    Input, int PAR_NUM, the total number of parameters.
        //    1 <= PAR_NUM.
        //
        //    Input, double Z[PAR_NUM*CHAIN_NUM*GEN_NUM], the Markov chain 
        //    sample data.
        //
    {
        int i;
        //
        //  Compute the standard deviations.
        //
        double[] std = StdCompute.std_compute(chain_num, gen_index, gen_num, par_num, z);
        //
        //  Increment the update count.
        //
        cr_ups[cr_index] += 1;
        //
        //  Update the CR distance.
        //
        for (i = 0; i < par_num; i++)
        {
            int i1 = i + chain_index * par_num + gen_index * par_num * chain_num;
            int i2 = i + chain_index * par_num + (gen_index - 1) * par_num * chain_num;
            cr_dis[cr_index] += Math.Pow((z[i2] - z[i1]) / std[i], 2);
        }
    }

    public static int cr_index_choose(int cr_num, double[] cr_prob)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CR_INDEX_CHOOSE chooses a CR index.
        //
        //  Discussion:
        //
        //    Index I is chosen with probability CR_PROB(I).
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
        //    Input, int CR_NUM, the total number of CR values.
        //    1 <= CR_NUM.
        //
        //    Input, double CR_PROB[CR_NUM], the probability of each CR.
        //
        //    Output, int CR_INDEX_CHOOSE, the index of the CR.
        //    0 <= CR_INDEX_CHOOSE < CR_NUM.
        //
    {
        int cr_index = 0;

        switch (cr_num)
        {
            case 1:
                cr_index = 0;
                break;
            default:
            {
                const int n = 1;
                int[] tmp_index = PDF.i4vec_multinomial_sample(n, cr_prob, cr_num);

                int i;
                for (i = 0; i < cr_num; i++)
                {
                    if (tmp_index[i] != 1)
                    {
                        continue;
                    }

                    cr_index = i;
                    break;
                }

                break;
            }
        }

        return cr_index;
    }

    public static void cr_init(ref double[] cr, ref double[] cr_dis, int cr_num, ref double[] cr_prob,
            ref int[] cr_ups )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CR_INIT initializes the crossover probability values.
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
        //    Output, double CR[CR_NUM], the CR values.
        //
        //    Output, double CR_DIS[CR_NUM], the CR distances.
        //
        //    Input, int CR_NUM, the total number of CR values.
        //    1 <= CR_NUM.
        //
        //    Output, double CR_PROB[CR_NUM], the probability of each CR.
        //
        //    Output, int CR_UPS[CR_NUM], the number of updates
        //    for each CR.
        //
    {
        int i;

        for (i = 0; i < cr_num; i++)
        {
            cr[i] = (i + 1) / (double) cr_num;
            cr_dis[i] = 1.0;
            cr_prob[i] = 1.0 / cr_num;
            cr_ups[i] = 1;
        }
    }

    public static void cr_prob_update(double[] cr_dis, int cr_num, ref double[] cr_prob,
            int[] cr_ups )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CR_PROB_UPDATE updates the CR probabilities.
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
        //    Input, double CR_DIS[CR_NUM], the CR distances.
        //
        //    Input, int CR_NUM, the total number of CR values.
        //    1 <= CR_NUM.
        //
        //    Output, double CR_PROB[CR_NUM], the updated CR probabilities.
        //
        //    Input, int CR_UPS[CR_NUM], the number of updates 
        //    for each CR.
        //
    {
        int i;

        for (i = 0; i < cr_num - 1; i++)
        {
            cr_prob[i] = cr_dis[i] / cr_ups[i];
        }

        double cr_prob_sum = typeMethods.r8vec_sum(cr_num, cr_prob);

        for (i = 0; i < cr_num - 1; i++)
        {
            cr_prob[i] /= cr_prob_sum;
        }
    }
}