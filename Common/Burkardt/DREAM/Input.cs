using System;
using System.Globalization;

namespace Burkardt.DREAM;

public static class Input
{
    public static void input_print(string chain_filename, int chain_num, int cr_num,
            string gr_filename, double gr_threshold, int jumpstep, double[] limits,
            int gen_num, int pair_num, int par_num, int printstep,
            string restart_read_filename, string restart_write_filename )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INPUT_PRINT prints the data from the input file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 May 2013
        //
        //  Author:
        //
        //    John Burkardt
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
        //    Input, int CR_NUM, the total number of CR values.
        //    1 <= CR_NUM.
        //
        //    Input, string GR_FILENAME, the name of the file
        //    in which values of the Gelman-Rubin statistic will be recorded,
        //    or the empty string "" if this file is not to be written.
        //
        //    Input, double GR_THRESHOLD, the convergence tolerance for the
        //    Gelman-Rubin statistic.
        //
        //    Input, int JUMPSTEP, forces a "long jump" every
        //    JUMPSTEP generations.
        //
        //    Input, double LIMITS[2*PAR_NUM], lower and upper limits
        //    for each parameter.
        //
        //    Input, int GEN_NUM, the total number of generations.
        //    2 <= GEN_NUM.
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
        //    Local, string RESTART_READ_FILENAME, the name of the file
        //    containing restart information.  If the calculation is not a restart,
        //    then this should be set to "".
        //
        //    Local, string RESTART_WRITE_FILENAME, the name of the file
        //    to be written, containing restart information.  If a restart file is
        //    not to be written, this should be set to "".
        //
    {
        int j;

        Console.WriteLine("");
        Console.WriteLine("INPUT_PRINT:");
        Console.WriteLine("");
        Console.WriteLine("  Number of parameters");
        Console.WriteLine("  PAR_NUM = " + par_num + "");
        Console.WriteLine("");
        Console.WriteLine("  LIMITS: Lower and upper limits for each parameter:");
        Console.WriteLine("");
        Console.WriteLine("  Index           Lower           Upper");
        Console.WriteLine("");
        for (j = 0; j < par_num; j++)
        {
            Console.WriteLine("  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(5)
                                   + "  " + limits[0 + j * 2].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + limits[1 + j * 2].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Number of generations:");
        Console.WriteLine("  GEN_NUM = " + gen_num + "");
        Console.WriteLine("");
        Console.WriteLine("  Number of simultaneous chains:");
        Console.WriteLine("  CHAIN_NUM = " + chain_num + "");
        Console.WriteLine("");
        Console.WriteLine("  Chain filename (base):");
        switch (chain_filename.Length)
        {
            case 0:
                Console.WriteLine("  CHAIN_FILENAME = \"(Null)\".");
                break;
            default:
                Console.WriteLine("  CHAIN_FILENAME = \"" + chain_filename + "\".");
                break;
        }

        Console.WriteLine("");
        Console.WriteLine("  Number of pairs of chains for crossover:");
        Console.WriteLine("  PAIR_NUM = " + pair_num + "");
        Console.WriteLine("");
        Console.WriteLine("  Number of crossover values:");
        Console.WriteLine("  CR_NUM = " + cr_num + "");
        Console.WriteLine("");
        Console.WriteLine("  Number of steps til a long jump:");
        Console.WriteLine("  JUMPSTEP = " + jumpstep + "");
        Console.WriteLine("");
        Console.WriteLine("  Interval between Gelman-Rubin computations:");
        Console.WriteLine("  PRINTSTEP = " + printstep + "");
        Console.WriteLine("");
        Console.WriteLine("  Gelman-Rubin data filename:");
        switch (gr_filename.Length)
        {
            case 0:
                Console.WriteLine("  GR_FILENAME = \"(Null)\".");
                break;
            default:
                Console.WriteLine("  GR_FILENAME = \"" + gr_filename + "\".");
                break;
        }

        Console.WriteLine("  GR_THRESHOLD = " + gr_threshold + "");
        Console.WriteLine("");
        Console.WriteLine("  Restart read filename:");
        switch (restart_read_filename.Length)
        {
            case 0:
                Console.WriteLine("  RESTART_READ_FILENAME = \"(Null)\".");
                break;
            default:
                Console.WriteLine("  RESTART_READ_FILENAME = \"" + restart_read_filename + "\".");
                break;
        }

        Console.WriteLine("");
        Console.WriteLine("  Restart write filename:");
        switch (restart_write_filename.Length)
        {
            case 0:
                Console.WriteLine("  RESTART_WRITE_FILENAME = \"(Null)\".");
                break;
            default:
                Console.WriteLine("  RESTART_WRITE_FILENAME = \"" + restart_write_filename + "\".");
                break;
        }
    }
}