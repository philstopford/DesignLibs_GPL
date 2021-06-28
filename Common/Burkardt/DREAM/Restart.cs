using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Burkardt.DREAM
{
    public static class Restart
    {
        public static void restart_read(int chain_num, ref double[] fit, int gen_num, int par_num,
        string restart_read_filename, ref double[] z )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RESTART_READ reads parameter sample data from a restart file.
        //
        //  Discussion:
        //
        //    Only a single generation (presumably the last one) was written to the file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 April 2013
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
        //    Output, double FIT[CHAIN_NUM*GEN_NUM], the likelihood of
        //    each sample.
        //
        //    Input, int GEN_NUM, the total number of generations.
        //    2 <= GEN_NUM.
        //
        //    Input, int PAR_NUM, the total number of parameters.
        //    1 <= PAR_NUM.
        //
        //    Input, string RESTART_READ_FILENAME, the name of 
        //    the restart file.
        //
        //    Output, double Z[PAR_NUM*CHAIN_NUM*GEN_NUM], the Markov chain 
        //    sample data.
        //
        {
            int chain_index;
            int gen_index = 0;
            int index;
            string line;
            int par_index;
            string[] restart;

            try
            {
                //
                //  Read and ignore line 1.
                //
                restart = File.ReadAllLines(restart_read_filename).Skip(1).ToArray();
                
                //
                //  Assume only one generation.
                //
                gen_index = 0;
                //
                //  Read the final fitness and parameter values for each chain.
                //
                int lineIndex = 0;
                for (chain_index = 0; chain_index < chain_num; chain_index++)
                {
                    index = chain_index
                            + chain_num * gen_index;
                    fit[index] = Convert.ToDouble(restart[lineIndex]);
                    lineIndex++;
                    for (par_index = 0; par_index < par_num; par_index++)
                    {
                        index = par_index
                                + par_num * chain_index
                                + par_num * chain_num * gen_index;
                        z[index] = Convert.ToDouble(restart[lineIndex]);
                        lineIndex++;
                    }
                }
            }
            catch (Exception e)
            {
                Console.WriteLine("");
                Console.WriteLine("RESTART_READ - Fatal error!");
                Console.WriteLine("  Could not open the file \""
                                  + restart_read_filename + "\".");
            }
        }

        public static void restart_write(int chain_num, double[] fit, int gen_num, int par_num,
        string restart_write_filename, double[] z )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RESTART_WRITE writes a restart file.
        //
        //  Discussion:
        //
        //    Only data for the final generation is written.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 April 2013
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
        //    Input, double FIT[CHAIN_NUM*GEN_NUM], the likelihood of
        //    each sample.
        //
        //    Input, int GEN_NUM, the total number of generations.
        //    2 <= GEN_NUM.
        //
        //    Input, int PAR_NUM, the total number of parameters.
        //    1 <= PAR_NUM.
        //
        //    Input, string RESTART_WRITE_FILENAME, the name of the 
        //    restart file.
        //
        //    Input, double Z[PAR_NUM*CHAIN_NUM*GEN_NUM], the Markov chain 
        //    sample data.
        //
        {
            int c;
            int p;
            List<string> restart = new List<string>();
            

            restart.Add("DREAM.C:Parameter_values_for_restart.");

            for (c = 0; c < chain_num; c++)
            {
                string tmp =  "  " + c
                    + "  " + fit[c + (gen_num - 1) * chain_num];
                for (p = 0; p < par_num; p++)
                {
                    tmp += "  " + z[p + c * par_num + (gen_num - 1) * par_num * chain_num];
                }

                restart.Add(tmp);
            }

            try
            {
                File.WriteAllLines(restart_write_filename, restart);
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("RESTART_WRITE - Fatal error!");
                Console.WriteLine("  Could not open \"" + restart_write_filename + "\".");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("RESTART_WRITE:");
            Console.WriteLine("  Created restart file \"" + restart_write_filename + "\".");

        }
    }
}