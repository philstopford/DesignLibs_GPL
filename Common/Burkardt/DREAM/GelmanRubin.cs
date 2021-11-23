using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using Burkardt.Types;

namespace Burkardt.DREAM;

public static class GelmanRubin
{
    public static void gr_compute(int chain_num, int gen_index, int gen_num, ref double[] gr,
            ref bool gr_conv, ref int gr_count, int gr_num, double gr_threshold, int par_num,
            double[] z )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GR_COMPUTE computes the Gelman Rubin statistics R used to check
        //    convergence.
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
        //    0 < GEN_INDEX < GEN_NUM.
        //
        //    Input, int GEN_NUM, the total number of generations.
        //    2 <= GEN_NUM.
        //
        //    Output, double GR[PAR_NUM*GR_NUM], the Gelman-Rubin R statistic.
        //
        //    Output, int &GR_CONV, the Gelman-Rubin convergence flag.
        //
        //    Input/output, int &GR_COUNT, counts the number of 
        //    generations at which the Gelman-Rubin statistic has been computed.
        //
        //    Input, int GR_NUM, the number of times the Gelman-Rubin
        //    statistic may be computed.
        //
        //    Input, double GR_THRESHOLD, the convergence tolerance for the
        //    Gelman-Rubin statistic.
        //
        //    Input, int PAR_NUM, the total number of parameters.
        //    1 <= PAR_NUM.
        //
        //    Input, double Z[PAR_NUM*CHAIN_NUM*GEN_NUM], the Markov chain 
        //    sample data.
        //
    {
        int par_index;

        int ind0 = (gen_index + 1) / 2 - 1;
        double rnd0 = ind0 + 1;

        double[] mean_chain = new double[chain_num];

        for (par_index = 0; par_index < par_num; par_index++)
        {
            int chain_index;
            int k;
            for (chain_index = 0; chain_index < chain_num; chain_index++)
            {
                mean_chain[chain_index] = 0.0;
                for (k = ind0; k <= gen_index; k++)
                {
                    mean_chain[chain_index] += z[par_index + chain_index * par_num + k * par_num * chain_num];
                }

                mean_chain[chain_index] /= rnd0;
            }

            double mean_all = typeMethods.r8vec_sum(chain_num, mean_chain) / chain_num;

            double b_var = 0.0;
            for (chain_index = 0; chain_index < chain_num; chain_index++)
            {
                b_var += Math.Pow(mean_chain[chain_index] - mean_all, 2);
            }

            b_var = rnd0 * b_var / (chain_num - 1);

            double s_sum = 0.0;
            for (chain_index = 0; chain_index < chain_num; chain_index++)
            {
                double s = 0.0;
                for (k = ind0; k <= gen_index; k++)
                {
                    s += Math.Pow(z[par_index + chain_index * par_num + k * par_num * chain_num]
                                  - mean_chain[chain_index], 2);
                }

                s_sum += s;
            }

            s_sum /= rnd0 - 1.0;

            double w_var = s_sum / chain_num;

            double var = ((rnd0 - 1.0) * w_var + b_var) / rnd0;

            gr[par_index + gr_count * par_num] = Math.Sqrt(var / w_var);
        }

        //
        //  Set the convergence flag.
        //
        gr_conv = true;

        for (par_index = 0; par_index < par_num; par_index++)
        {
            if (gr_threshold < gr[par_index + gr_count * par_num])
            {
                gr_conv = false;
                break;
            }
        }

        switch (gr_conv)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("GR_COMPUTE:");
                Console.WriteLine("  GR convergence at iteration: " + gen_index + "");
                break;
        }
            
        gr_count += 1;
    }

    public static void gr_init(double[] gr, ref bool gr_conv, ref int gr_count, int gr_num,
            int par_num )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GR_INIT initializes Gelman-Rubin variables.
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
        //    Output, double GR[PAR_NUM*GR_NUM], the Gelman-Rubin statistic.
        //
        //    Output, int &GR_CONV, the convergence flag.
        //
        //    Output, int &GR_COUNT, counts the number of generations
        //    at which the Gelman-Rubin statistic has been computed.
        //
        //    Input, int GR_NUM, the number of times the Gelman-Rubin
        //    statistic may be computed.
        //
        //    Input, int PAR_NUM, the number of parameters.
        //    1 <= PAR_NUM.
        //
    {
        int j;

        for (j = 0; j < gr_num; j++)
        {
            int i;
            for (i = 0; i < par_num; i++)
            {
                gr[i + j * par_num] = 0.0;
            }
        }

        gr_conv = false;
        gr_count = 0;
    }

    public static void gr_write(double[] gr, string gr_filename, int gr_num, int par_num,
            int printstep )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GR_WRITE writes Gelman-Rubin R statistics into a file.
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
        //    Input, double GR[PAR_NUM*GR_NUM], the Gelman-Rubin R statistic.
        //
        //    Input, string GR_FILENAME, the Gelman-Rubin filename.
        //
        //    Input, int GR_NUM, the number of times the Gelman-Rubin
        //    statistic may be computed.
        //
        //    Input, int PAR_NUM, the total number of parameters.
        //    1 <= PAR_NUM.
        //
        //    Input, int PRINTSTEP, the interval between generations on 
        //    which the Gelman-Rubin statistic will be computed and written to a file.
        //
    {
        List<string> gr_unit = new();
        int j;


        gr_unit.Add("DREAM.CPP:Monitored_parameter_interchains_Gelman_Rubin_R_statistic");

        for (j = 0; j < gr_num; j++)
        {
            string tmp = (printstep * (j + 1) - 1).ToString(CultureInfo.InvariantCulture);
            int i;
            for (i = 0; i < par_num; i++)
            {
                tmp += "  " + gr[i + j * par_num];
            }

            gr_unit.Add(tmp);
        }

        try
        {
            File.WriteAllLines(gr_filename, gr_unit);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("GR_WRITE - Fatal error!");
            Console.WriteLine("  Could not open the file \"" + gr_filename + "\"");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("GR_WRITE:");
        Console.WriteLine("  Created the file \"" + gr_filename + "\".");
    }
}