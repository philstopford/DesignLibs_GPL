using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using Burkardt.IO;
using Burkardt.Types;

namespace Burkardt.DREAM;

public static class Chain
{
    public static void chain_init(int chain_num, double[] fit, int gen_num, int par_num, Func<int, Dream.SampleResult> prior_sample, Func<int, double[], double> sample_likelihood, 
            ref double[] z )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHAIN_INIT starts Markov chains from a prior distribution.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 May 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int CHAIN_NUM, the total number of chains.
        //    3 <= CHAIN_NUM.
        //
        //    Output, double FIT[CHAIN_NUM*GEN_NUM], the likelihood of
        //    each sample.
        //
        //    Input, int GEN_NUM, the total number of generations.
        //    2 <= GEN_NUM.
        //
        //    Input, int PAR_NUM, the total number of parameters.
        //    1 <= PAR_NUM.
        //
        //    Output, double Z[PAR_NUM*CHAIN_NUM*GEN_NUM], the Markov chain 
        //    sample data.
        //
    {
        int c;

        for (c = 0; c < chain_num; c++)
        {
            double[] zp = prior_sample(par_num).result;

            int p;
            for (p = 0; p < par_num; p++)
            {
                z[p + c * par_num + 0 * par_num * chain_num] = zp[p];
            }

            fit[c + 0 * chain_num] = sample_likelihood(par_num, zp);
        }
    }

    public static void chain_init_print(int chain_num, double[] fit, int gen_num, int par_num,
            string restart_read_filename, double[] z )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHAIN_INIT_PRINT prints the initial values for Markov chains.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 May 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int CHAIN_NUM, the total number of chains.
        //    3 <= CHAIN_NUM.
        //
        //    Input, double FIT[CHAIN_NUM*GEN_NUM], the likelihood of
        //    each sample.
        //
        //    Input, int GEN_NUM, the total number of generations.
        //    2 <= GEN_NUM.
        //
        //    Input, int PAR_NUM, the total number of parameters.
        //    1 <= PAR_NUM.
        //
        //    Input, string RESTART_READ_FILENAME, the name of the 
        //    restart file.
        //
        //    Input, double Z[PAR_NUM*CHAIN_NUM*GEN_NUM], the Markov chain 
        //    sample data.
        //
    {
        int j;

        Console.WriteLine("");
        Console.WriteLine("CHAIN_INIT_PRINT");
        Console.WriteLine("  Display initial values of Markov chains.");

        switch (restart_read_filename.Length)
        {
            case > 0:
                Console.WriteLine("  Initialization from restart file \""
                                  + restart_read_filename + "\"");
                break;
            default:
                Console.WriteLine("  Initialization by sampling prior density.");
                break;
        }

        for (j = 0; j < chain_num; j++)
        {
            Console.WriteLine("");
            Console.WriteLine("  Chain " + j + "");
            Console.WriteLine("  Fitness " + fit[j + 0 * chain_num] + "");
            string cout = "";
            int i;
            for (i = 0; i < par_num; i++)
            {
                cout += "  " + z[i + j * par_num + 0 * par_num * chain_num].ToString(CultureInfo.InvariantCulture).PadLeft(14);
                if ((i + 1) % 5 == 0 || i == par_num - 1)
                {
                    Console.WriteLine(cout);
                }
            }
        }
    }

    public static void chain_outliers(int chain_num, int gen_index, int gen_num, int par_num,
            double[] fit, double[] z )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHAIN_OUTLIERS identifies and modifies outlier chains during burn-in.
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
        //    2 <= GEN_INDEX <= GEN_NUM.
        //
        //    Input, int GEN_NUM, the total number of generations.
        //    2 <= GEN_NUM.
        //
        //    Input, int PAR_NUM, the total number of parameters.
        //    1 <= PAR_NUM.
        //
        //    Input/output, double FIT[CHAIN_NUM*GEN_NUM], the likelihood of
        //    each sample.
        //
        //    Input/output, double Z[PAR_NUM*CHAIN_NUM*GEN_NUM], the Markov
        //    chain sample data.
        //
    {
        int j;
        int k;

        int klo = (gen_index + 1) / 2 - 1;
        int knum = gen_index + 1 - klo;

        double[] avg = new double[chain_num];

        for (j = 0; j < chain_num; j++)
        {
            double t = 0.0;
            for (k = klo; k <= gen_index; k++)
            {
                t += fit[j + k * chain_num];
            }

            avg[j] = t / knum;
        }

        //
        //  Set BEST to be the index of the chain with maximum average.
        //
        int best = 0;
        double avg_max = avg[0];
        for (j = 1; j < chain_num; j++)
        {
            if (!(avg_max < avg[j]))
            {
                continue;
            }

            best = j;
            avg_max = avg[j];
        }

        //
        //  Determine the indices of the chains having averages 1/4 "above" 
        //  and "below" the average.
        //
        double[] avg_sorted = typeMethods.r8vec_copy_new(chain_num, avg);

        typeMethods.r8vec_sort_heap_a(chain_num, ref avg_sorted);

        int ind1 = (int)Math.Round(0.25 * chain_num);
        int ind3 = (int)Math.Round(0.75 * chain_num);

        double q1 = avg_sorted[ind1];
        double q3 = avg_sorted[ind3];
        double qr = q3 - q1;

        //
        //  Identify outlier chains, and replace their later samples
        //  with values from the "best" chain.
        //
        int outlier_num = 0;
        for (j = 0; j < chain_num; j++)
        {
            if (!(avg[j] < q1 - 2.0 * qr))
            {
                continue;
            }

            outlier_num += 1;
            int i;
            for (i = 0; i < par_num; i++)
            {
                z[i + j * par_num + gen_index * par_num * chain_num] =
                    z[i + best * par_num + gen_index * par_num * chain_num];
            }

            for (k = klo; k <= gen_index; k++)
            {
                fit[j + k * chain_num] = fit[best + k * chain_num];
            }
        }

        switch (outlier_num)
        {
            //
            //  List the outlier chains.
            //
            case > 0:
            {
                Console.WriteLine("");
                Console.WriteLine("CHAIN_OUTLIERS:");
                Console.WriteLine("  At iteration " + gen_index
                                                    + " found " + outlier_num + " outlier chains,");
                Console.WriteLine("  whose indices appear below, and for which samples");
                Console.WriteLine("  from the chain with the largest log likelihood function,");
                Console.WriteLine("  index number " + best + " will be substituted.");

                for (j = 0; j < chain_num; j++)
                {
                    if (avg[j] < q1 - 2.0 * qr)
                    {
                        Console.WriteLine("  " + j + "");
                    }
                }

                break;
            }
        }
    }

    public static void chain_write(string chain_filename, int chain_num, double[] fit,
            int gen_num, int par_num, double[] z )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHAIN_WRITE writes samples of each chain to separate files.
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
        //    C++ version John Burkardt.
        //
        //  Parameters:
        //
        //    Input, string CHAIN_FILENAME, the "base" filename
        //    to be used for the chain files.  If this is ""
        //    then the chain files will not be written.  This name should 
        //    include a string of 0's which will be replaced by the chain 
        //    indices.  For example, "chain000.txt" would work as long as the
        //    number of chains was 1000 or less.
        //
        //    Input, int CHAIN_NUM, the total number of chains.
        //    3 <= CHAIN_NUM.
        //
        //    Input, double FIT[CHAIN_NUM*GEN_NUM], the likelihood of
        //    each sample.
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
        List<string> chain = new();
        int j;
        //
        //  Make a temporary copy of the filename template, which we can alter.
        //
        string chain_filename2 = chain_filename;
        //
        //  Write parameter samples of all chains.
        //
        Console.WriteLine("");
        Console.WriteLine("CHAIN_WRITE:");

        for (j = 0; j < chain_num; j++)
        {
            chain.Add("DREAM.CPP:Parameters_and_log_likelihood_for_chain_#" + j + "");

            int k;
            for (k = 0; k < gen_num; k++)
            {
                string tmp = "  " + k
                                  + "  " + fit[j + k * chain_num];
                int i;
                for (i = 0; i < par_num; i++)
                {
                    tmp += "  " + z[i + j * par_num + k * par_num * chain_num];
                }

                chain.Add(tmp);
            }

            try
            {
                File.WriteAllLines(chain_filename2, chain);
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("CHAIN_WRITE - Fatal error!");
                Console.WriteLine("  Could not open file \"" + chain_filename2 + "\".");
                return;
            }

            chain.Clear();

            Console.WriteLine("  Created file \"" + chain_filename2 + "\".");

            Files.filename_inc(ref chain_filename2);
        }
    }
}