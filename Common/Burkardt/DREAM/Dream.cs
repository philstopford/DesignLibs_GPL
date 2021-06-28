using System;
using Burkardt.PDFLib;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.DREAM
{
    public class ProblemSize
    {
        public int chain_num { get; set; }
        public int cr_num { get; set; }
        public int gen_num { get; set; }
        public int pair_num { get; set; }
        public int par_num { get; set; }
    }

    public class ProblemValue
    {
        public string chain_filename { get; set; }
        public string gr_filename { get; set; }
        public double gr_threshold { get; set; }
        public int jumpstep { get; set; }
        public double[] limits { get; set; }
        public int par_num { get; set; }
        public int printstep { get; set; }
        public string restart_read_filename { get; set; }
        public string restart_write_filename { get; set; }
    }

    public class DreamDelegates
    {
        public delegate double[] prior_sample(int x);

        public delegate double prior_density(int a, double[] b, int c);
        
        public delegate double sample_likelihood(int a, double[] b, double c);
    }
    
    public static class Dream
    {
        public static void dream(ref ProblemSize problem_size,
                ref ProblemValue problem_value_data,
                Func<int, double[]> prior_sample,
                Func<int, double[], int, double> prior_density,
                Func<int, double[], double> sample_likelihood)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for DREAM.
            //
            //  Discussion:
            //
            //    The DREAM program was originally developed by Guannan Zhang, of
            //    Oak Ridge National Laboratory (ORNL); it has been incorporated into 
            //    the DAKOTA package of Sandia National Laboratory, and is
            //    intended to form part of the ORNL package known as TASMANIA.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 June 2013
            //
            //  Author:
            //
            //    John Burkardt
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
            //  Local parameters:
            //
            //    Input, string CHAIN_FILENAME, the "base" filename
            //    to be used for the chain files.  If this is ""
            //    then the chain files will not be written.  This name should 
            //    include a string of 0's which will be replaced by the chain 
            //    indices.  For example, "chain000.txt" would work as long as the
            //    number of chains was 1000 or less.
            //
            //    Local, int CHAIN_NUM, the total number of chains.
            //    3 <= CHAIN_NUM.
            //
            //    Local, int CR_NUM, the total number of CR values.
            //    1 <= CR_NUM.
            //
            //    Local, double FIT[CHAIN_NUM*GEN_NUM], the likelihood of
            //    each sample.
            //
            //    Local, int GEN_NUM, the total number of generations.
            //    2 <= GEN_NUM.
            //
            //    Local, double GR[PAR_NUM*GR_NUM], 
            //    the Gelman-Rubin R statistic.
            //
            //    Local, logical GR_CONV, the Gelman-Rubin convergence flag.
            //
            //    Local, int GR_COUNT, counts the number of generations
            //    at which the Gelman-Rubin statistic has been computed.
            //
            //    Local, string GR_FILENAME, the name of the file
            //    in which values of the Gelman-Rubin statistic will be recorded,
            //    or the empty string "" if this file is not to be written.
            //
            //    Local, int GR_NUM, the number of times the Gelman-Rubin
            //    statistic may be computed.
            //
            //    Local, double GR_THRESHOLD, the convergence tolerance for
            //    the Gelman-Rubin statistic.
            //
            //    Local, double JUMPRATE_TABLE[PAR_NUM], the jumprate table.
            //
            //    Local, int JUMPSTEP, forces a "long jump" every
            //    JUMPSTEP generations.
            //
            //    Local, double LIMITS[2*PAR_NUM], lower and upper bounds
            //    for each parameter.
            //
            //    Local, int PAIR_NUM, the number of pairs of 
            //    crossover chains.
            //    0 <= PAIR_NUM.
            //
            //    Local, int PAR_NUM, the total number of parameters.
            //    1 <= PAR_NUM.
            //
            //    Local, int PRINTSTEP, the interval between generations on 
            //    which the Gelman-Rubin statistic will be computed and written to a file.
            //
            //    Local, string RESTART_READ_FILENAME, the name of the file
            //    containing restart information.  If the calculation is not a restart,
            //    then this should be set to "".
            //
            //    Local, string RESTART_WRITE_FILENAME, the name of the file
            //    to be written, containing restart information.  If a restart file is
            //    not to be written, this should be set to "".
            //
            //    Local, double Z[PAR_NUM*CHAIN_NUM*GEN_NUM], the Markov chain 
            //    sample data.
            //
        {
            double[] fit;
            double[] gr;
            bool gr_conv = false;
            int gr_count = 0;
            int gr_num = 0;
            double[] jumprate_table;
            double[] z;
            
            Console.WriteLine("");
            Console.WriteLine("DREAM");
            Console.WriteLine("  MCMC acceleration by Differential Evolution.");
            //
            //  Get the problem sizes.
            //
            //
            //  Decide if the problem sizes are acceptable.
            //
            if (problem_size.chain_num < 3)
            {
                Console.WriteLine("");
                Console.WriteLine("DREAM - Fatal error!");
                Console.WriteLine("  CHAIN_NUM < 1.");
                return;
            }

            if (problem_size.cr_num < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("DREAM - Fatal error!");
                Console.WriteLine("  CR_NUM < 1.");
                return;
            }

            if (problem_size.gen_num < 2)
            {
                Console.WriteLine("");
                Console.WriteLine("DREAM - Fatal error!");
                Console.WriteLine("  GEN_NUM < 2.");
                return;
            }

            if (problem_size.pair_num < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("DREAM - Fatal error!");
                Console.WriteLine("  PAIR_NUM < 0.");
                return;
            }

            if (problem_size.par_num < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("DREAM - Fatal error!");
                Console.WriteLine("  PAR_NUM < 1.");
                return;
            }

            //
            //  Get the problem data values;
            //
            problem_value_data.limits = typeMethods.r8mat_zero_new(2, problem_size.par_num);

            //
            //  Print the data as a job record.
            //
            Input.input_print(problem_value_data.chain_filename, problem_size.chain_num, problem_size.cr_num, problem_value_data.gr_filename, problem_value_data.gr_threshold,
                problem_value_data.jumpstep, problem_value_data.limits, problem_size.gen_num, problem_size.pair_num, problem_size.par_num,
                problem_value_data.printstep, problem_value_data.restart_read_filename, problem_value_data.restart_write_filename);
            //
            //  Allocate and zero out memory.
            //
            gr_num = problem_size.gen_num / problem_value_data.printstep;

            fit = typeMethods.r8mat_zero_new(problem_size.chain_num, problem_size.gen_num);
            gr = typeMethods.r8mat_zero_new(problem_size.par_num, gr_num);
            z = typeMethods.r8block_zero_new(problem_size.par_num, problem_size.chain_num, problem_size.gen_num);
            //
            //  Set the jump rate table.
            //
            jumprate_table = Jump.jumprate_table_init(problem_size.pair_num, problem_size.par_num);

            Jump.jumprate_table_print(jumprate_table, problem_size.pair_num, problem_size.par_num);
            //
            //  Initialize the Gelman-Rubin data.
            //
            GelmanRubin.gr_init(gr, ref gr_conv, ref gr_count, gr_num, problem_size.par_num);

            Console.WriteLine("");
            Console.WriteLine("GR_PRINT:");
            Console.WriteLine("  GR_CONV  = " + gr_conv + "");
            Console.WriteLine("  GR_COUNT = " + gr_count + "");
            Console.WriteLine("  GR_NUM   = " + gr_num + "");
            //
            //  Set the first generation of the chains from restart data, or by sampling.
            //
            if (problem_value_data.restart_read_filename.Length == 0)
            {
                Chain.chain_init(problem_size.chain_num, fit, problem_size.gen_num, problem_size.par_num, prior_sample, sample_likelihood, ref z);
            }
            else
            {
                Restart.restart_read(problem_size.chain_num, ref fit, problem_size.gen_num, problem_size.par_num, problem_value_data.restart_read_filename, ref z);
            }

            Chain.chain_init_print(problem_size.chain_num, fit, problem_size.gen_num, problem_size.par_num, problem_value_data.restart_read_filename,
                z);
            //
            //  Carry out the DREAM algorithm.
            //
            Dream.dream_algm(problem_size.chain_num, problem_size.cr_num, fit, problem_size.gen_num, gr, ref gr_conv, ref gr_count,
                gr_num, problem_value_data.gr_threshold, jumprate_table, problem_value_data.jumpstep, problem_value_data.limits, problem_size.pair_num,
                problem_size.par_num, problem_value_data.printstep, prior_density, sample_likelihood, ref z);
            //
            //  Save Gelman-Rubin statistic values to a file.
            //
            if (problem_value_data.gr_filename.Length > 0)
            {
                GelmanRubin.gr_write(gr, problem_value_data.gr_filename, gr_num, problem_size.par_num, problem_value_data.printstep);
            }

            //
            //  Save parameter values for all chains at last generation.
            //
            if (problem_value_data.restart_write_filename.Length > 0)
            {
                Restart.restart_write(problem_size.chain_num, fit, problem_size.gen_num, problem_size.par_num, problem_value_data.restart_write_filename,
                    z);
            }

            //
            //  Write each chain to a separate file.
            //
            if (problem_value_data.chain_filename.Length > 0)
            {
                Chain.chain_write(problem_value_data.chain_filename, problem_size.chain_num, fit, problem_size.gen_num, problem_size.par_num, z);
            }

            Console.WriteLine("");
            Console.WriteLine("DREAM");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }
        
        
        public static void dream_algm(int chain_num, int cr_num, double[] fit, int gen_num,
        double[] gr, ref bool gr_conv, ref int gr_count, int gr_num, double gr_threshold,
        double[] jumprate_table, int jumpstep, double[] limits, int pair_num,
        int par_num, int printstep,
        Func<int, double[], int, double> prior_density,
        Func<int, double[], double> sample_likelihood,
        ref double[] z )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DREAM_ALGM gets a candidate parameter sample.
        //
        // Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 May 2013
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
        //    Input, int CR_NUM, the total number of CR values.
        //    1 <= CR_NUM.
        //
        //    Input, double FIT[CHAIN_NUM*GEN_NUM], the likelihood of
        //    each sample.
        //
        //    Input, int GEN_NUM, the total number of generations.
        //    2 <= GEN_NUM.
        //
        //    Input, double GR[PAR_NUM*GR_NUM], 
        //    the Gelman-Rubin R statistic.
        //
        //    Input/output, int &GR_CONV, the Gelman-Rubin convergence flag.
        //
        //    Input/output, int &GR_COUNT, counts the number of generations
        //    at which the Gelman-Rubin statistic has been computed.
        //
        //    Input, int GR_NUM, the number of times the Gelman-Rubin
        //    statistic may be computed.
        //
        //    Input, double GR_THRESHOLD, the convergence tolerance for
        //    the Gelman-Rubin statistic.
        //
        //    Input, double JUMPRATE_TABLE[PAR_NUM], the jumprate table.
        //
        //    Input, int JUMPSTEP, forces a "long jump" every
        //    JUMPSTEP generations.
        //
        //    Input, double LIMITS[2*PAR_NUM], lower and upper bounds
        //    for each parameter.
        //
        //    Input, int PAIR_NUM, the number of pairs of 
        //    crossover chains.
        //    0 <= PAIR_NUM.
        //
        //    Input, int PAR_NUM, the total number of parameters.
        //    1 <= PAR_NUM.
        //
        //    Input, int PRINTSTEP, the interval between generations on 
        //    which the Gelman-Rubin statistic will be computed and written to a file.
        //
        //    Input/output, double Z[PAR_NUM*CHAIN_NUM*GEN_NUM], the Markov chain 
        //    sample data.
        //
        //  Local parameters:
        //
        //    Local, int CHAIN_INDEX, the index of the current chain.
        //    1 <= CHAIN_INDEX <= CHAIN_NUM.
        //
        //    Local, double CR[CR_NUM], the CR values.
        //
        //    Local, double CR_DIS[CR_NUM], the CR distances.
        //
        //    Local, int CR_INDEX, the index of the selected CR value.
        //    1 <= CR_INDEX <= CR_NUM.
        //
        //    Local, double CR_PROB[CR_NUM], the probability of each CR.
        //
        //    Local, double CR_UPS[CR_NUM], the number of updates for each CR.
        //
        //    Local, int GEN_INDEX, the index of the current generation.
        //    1 <= GEN_INDEX <= GEN_NUM.
        //
        //    Local, double ZP[PAR_NUM], a candidate sample.
        //
        //    Local, int ZP_ACCEPT, the number of candidates accepted.
        //
        //    Local, double ZP_ACCEPT_RATE, the rate at which generated
        //    candidates were accepted.
        //
        //    Local, int ZP_COUNT, the number of candidates generated.
        //
        //    Local, double ZP_RATIO, the Metropolis ratio for a candidate.
        //
        {
            int chain_index;
            double[] cr;
            double[] cr_dis;
            int cr_index;
            double[] cr_prob;
            int[] cr_ups;
            int gen_index;
            int i;
            double pd1;
            double pd2;
            double r;
            double[] zp;
            int zp_accept;
            double zp_accept_rate;
            int zp_count;
            double zp_fit;
            double[] zp_old;
            double zp_old_fit;
            double zp_ratio;

            zp_old = new double[par_num];
            zp_count = 0;
            zp_accept = 0;
            //
            //  Initialize the CR values.
            //
            cr = new double[cr_num];
            cr_dis = new double[cr_num];
            cr_prob = new double[cr_num];
            cr_ups = new int[cr_num];

            CR.cr_init(ref cr, ref cr_dis, cr_num, ref cr_prob, ref cr_ups);

            for (gen_index = 1; gen_index < gen_num; gen_index++)
            {
                for (chain_index = 0; chain_index < chain_num; chain_index++)
                {
                    //
                    //  Choose CR_INDEX, the index of a CR.
                    //
                    cr_index = CR.cr_index_choose(cr_num, cr_prob);
                    //
                    //  Generate a sample candidate ZP.
                    //
                    zp = Sample.sample_candidate(chain_index, chain_num, cr, cr_index, cr_num,
                        gen_index, gen_num, jumprate_table, jumpstep, limits, pair_num,
                        par_num, z);

                    zp_count = zp_count + 1;
                    //
                    //  Compute the log likelihood function for ZP.
                    //
                    zp_fit = sample_likelihood(par_num, zp);

                    for (i = 0; i < par_num; i++)
                    {
                        zp_old[i] = z[i + chain_index * par_num + (gen_index - 1) * par_num * chain_num];
                    }

                    zp_old_fit = fit[chain_index + (gen_index - 1) * chain_num];
                    //
                    //  Compute the Metropolis ratio for ZP.
                    //
                    pd1 = prior_density(par_num, zp, 0);

                    pd2 = prior_density(par_num,
                        z, + 0 + chain_index * par_num + (gen_index - 1) * par_num * chain_num);

                    zp_ratio = Math.Exp(
                        (zp_fit + Math.Log(pd1)) -
                        (zp_old_fit + Math.Log(pd2)));

                    zp_ratio = Math.Min(zp_ratio, 1.0);
                    //
                    //  Accept the candidate, or copy the value from the previous generation.
                    //
                    r = PDF.r8_uniform_01_sample();

                    if (r <= zp_ratio)
                    {
                        for (i = 0; i < par_num; i++)
                        {
                            z[i + chain_index * par_num + gen_index * par_num * chain_num] = zp[i];
                        }

                        zp_accept = zp_accept + 1;
                        fit[chain_index + gen_index * chain_num] = zp_fit;
                    }
                    else
                    {
                        for (i = 0; i < par_num; i++)
                        {
                            z[i + chain_index * par_num + gen_index * par_num * chain_num] = zp_old[i];
                        }

                        fit[chain_index + gen_index * chain_num] = zp_old_fit;
                    }

                    //
                    //  Update the CR distance.
                    //
                    if (!gr_conv)
                    {
                        if (1 < cr_num)
                        {
                            CR.cr_dis_update(chain_index, chain_num, ref cr_dis, cr_index,
                                cr_num, ref cr_ups, gen_index, gen_num, par_num, z);
                        }
                    }
                }

                //
                //  Update the multinomial distribution of CR.
                //
                if (!gr_conv)
                {
                    if (1 < cr_num)
                    {
                        if ((gen_index + 1) % 10 == 0)
                        {
                            CR.cr_prob_update(cr_dis, cr_num, ref cr_prob, cr_ups);
                        }
                    }
                }

                //
                //  Every PRINTSTEP interval,
                //  * compute the Gelman Rubin R statistic for this generation,
                //    and determine if convergence has occurred.
                //
                if ((gen_index + 1) % printstep == 0)
                {
                    GelmanRubin.gr_compute(chain_num, gen_index, gen_num, ref gr, ref gr_conv, ref gr_count,
                        gr_num, gr_threshold, par_num, z);
                }

                //
                //  Check for outlier chains.
                //
                if (!gr_conv)
                {
                    if ((gen_index + 1) % 10 == 0)
                    {
                        Chain.chain_outliers(chain_num, gen_index, gen_num, par_num, fit, z);
                    }
                }
            }

            //
            //  Compute the acceptance rate.
            //
            zp_accept_rate = (double) (zp_accept) / (double) (zp_count);

            Console.WriteLine("");
            Console.WriteLine("  The acceptance rate is " + zp_accept_rate + "");
        }
    }
}