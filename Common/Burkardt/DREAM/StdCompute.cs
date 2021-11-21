using System;

namespace Burkardt.DREAM;

public static class StdCompute
{
    public static double[] std_compute ( int chain_num, int gen_index, int gen_num, int par_num, 
            double[] z )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STD_COMPUTE computes the current standard deviations, for each parameter.
        //
        //  Discussion:
        //
        //    The computation encompasses all chains and generations up to the
        //    current ones.
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
        //    Input, int CHAIN_NUM, the total number of chains.
        //    3 <= CHAIN_NUM.
        //
        //    Input, int GEN_INDEX, the current generation.
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
        //    Output, double STD_COMPUTE[PAR_NUM], the standard deviations.
        //
    {
        int i;

        double[] std = new double[par_num];

        for ( i = 0; i < par_num; i++ )
        {
            double mean = 0.0;
            int j;
            int k;
            for ( k = 0; k <= gen_index; k++ )
            {
                for ( j = 0; j < chain_num; j++ )
                {
                    mean += z[i+j*par_num+k*par_num*chain_num];
                }
            }
            mean = mean / chain_num / gen_index;

            std[i] = 0.0;
            for ( k = 0; k <= gen_index; k++ )
            {
                for ( j = 0; j < chain_num; j++ )
                {
                    std[i] += Math.Pow ( z[i+j*par_num+k*par_num*chain_num] - mean, 2 );
                }
            }
            std[i] = Math.Sqrt ( std[i] / (chain_num * gen_index - 1) );
        }

        return std;
    }

}