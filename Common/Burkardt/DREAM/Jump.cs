using System;
using Burkardt.PDFLib;

namespace Burkardt.DREAM;

public static class Jump
{
    public static void jumprate_choose(double[] cr, int cr_index, int cr_num, int gen_index,
            int[] jump_dim, ref int jump_num, ref double jumprate, double[] jumprate_table,
            int jumpstep, int par_num )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    JUMPRATE_CHOOSE chooses a jump rate from the jump rate table.
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
        //    Input, double CR[CR_NUM], the CR values.
        //
        //    Input, int CR_INDEX, the index of the CR.
        //    1 <= CR_INDEX <= CR_NUM.
        //
        //    Input, int CR_NUM, the total number of CR values.
        //    1 <= CR_NUM.
        //
        //    Input, int GEN_INDEX, the current generation.
        //    1 <= GEN_INDEX <= GEN_NUM.
        //
        //    Output, int JUMP_DIM[PAR_NUM], the indexes of the
        //    parameters to be updated.
        //
        //    Output, int &JUMP_NUM, the number of dimensions in which
        //    a jump will be made.  0 <= JUMP_NUM <= PAR_NUM.
        //
        //    Output, double &JUMPRATE, the jump rate.
        //
        //    Input, double JUMPRATE_TABLE[PAR_NUM], the jump rate table.
        //
        //    Input, int JUMPSTEP, forces a "long jump" every
        //    JUMPSTEP generations.
        //
        //    Input, int PAR_NUM, the total number of parameters.
        //    1 <= PAR_NUM.
        //
    {
        int i;
        double r;
        //
        //  Determine the dimensions that will be updated.
        //
        jump_num = 0;
        for (i = 0; i < par_num; i++)
        {
            jump_dim[i] = 0;
        }

        for (i = 0; i < par_num; i++)
        {
            r = PDF.r8_uniform_01_sample();

            if (1.0 - cr[cr_index] < r)
            {
                jump_dim[jump_num] = i;
                jump_num += 1;
            }
        }

        jumprate = (gen_index % jumpstep) switch
        {
            //
            //  Determine if a long jump is forced.
            //
            0 => 0.98,
            _ => par_num switch
            {
                //
                //  If parameter dimension is 1, 2, or 3, fix the jump rate to 0.6.
                //
                <= 3 => 0.6,
                _ => jump_num switch
                {
                    //
                    //  Calculate the general jump rate.
                    //
                    0 => 0.0,
                    _ => jumprate_table[jump_num - 1]
                }
            }
        };
    }

    public static double[] jumprate_table_init(int pair_num, int par_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    JUMPRATE_TABLE_INIT initializes the jump rate table.
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
        //    Input, int PAIR_NUM, the number of pairs of 
        //    crossover chains.
        //    0 <= PAIR_NUM.
        //
        //    Input, int PAR_NUM, the total number of parameters.
        //    1 <= PAR_NUM.
        //
        //    Output, double JUMPRATE_TABLE_INIT[PAR_NUM], the jumprate table.
        //
    {
        double c;
        int i;
        double[] jumprate_table;

        jumprate_table = new double[par_num];

        c = 2.38 / Math.Sqrt(2 * pair_num);

        for (i = 0; i < par_num; i++)
        {
            jumprate_table[i] = c / Math.Sqrt(i + 1);
        }

        return jumprate_table;
    }

    public static void jumprate_table_print(double[] jumprate_table, int pair_num, int par_num )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    JUMPRATE_TABLE_PRINT prints the jump rate table.
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
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double JUMPRATE_TABLE[PAR_NUM], the jumprate table.
        //
        //    Input, int PAIR_NUM, the number of pairs of 
        //    crossover chains.
        //    0 <= PAIR_NUM.
        //
        //    Input, int PAR_NUM, the total number of parameters.
        //    1 <= PAR_NUM.
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("JUMPRATE_TABLE_PRINT");
        Console.WriteLine("");
        Console.WriteLine("   I        Jumprate");
        Console.WriteLine("");
        for (i = 0; i < par_num; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(2)
                                   + "  " + jumprate_table[i].ToString().PadLeft(14) + "");
        }
    }
}