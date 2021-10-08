using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.FEM;
using Burkardt.Types;

namespace Burkardt.Function
{
    public static class Lebesgue
    {
        public static double lebesgue_constant(int n, double[] x, int nfun, double[] xfun)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEBESGUE_CONSTANT estimates the Lebesgue constant for a set of points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 March 2014
            //
            //  Author:
            //
            //    John Burkardt.
            //
            //  Parameters:
            //
            //    Jean-Paul Berrut, Lloyd Trefethen,
            //    Barycentric Lagrange Interpolation,
            //    SIAM Review,
            //    Volume 46, Number 3, September 2004, pages 501-517.
            //
            //  Parameters:
            //
            //    Input, int N, the number of interpolation points.
            //
            //    Input, double X[N], the interpolation points.
            //
            //    Input, int NFUN, the number of evaluation points.
            //
            //    Input, double XFUN[NFUN], the evaluation points.
            //
            //    Output, double LEBESGUE_CONSTANT, an estimate of the Lebesgue constant 
            //    for the points.
            //
        {
            double[] lfun;
            double lmax;

            lfun = lebesgue_function(n, x, nfun, xfun);

            lmax = typeMethods.r8vec_max(nfun, lfun);

            return lmax;
        }

        public static double[] lebesgue_function(int n, double[] x, int nfun, double[] xfun)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEBESGUE_FUNCTION evaluates the Lebesgue function for a set of points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 March 2014
            //
            //  Author:
            //
            //    John Burkardt.
            //
            //  Parameters:
            //
            //    Jean-Paul Berrut, Lloyd Trefethen,
            //    Barycentric Lagrange Interpolation,
            //    SIAM Review,
            //    Volume 46, Number 3, September 2004, pages 501-517.
            //
            //  Parameters:
            //
            //    Input, int N, the number of interpolation points.
            //
            //    Input, double X[N], the interpolation points.
            //
            //    Input, int NFUN, the number of evaluation points.
            //
            //    Input, double XFUN[NFUN], the evaluation points.
            //
            //    Output, double LEBESGUE_FUNCTION[NFUN], the Lebesgue function values.
            //
        {
            int i;
            int j;
            double[] lfun;
            double[] llfun;
            double t;

            lfun = new double[nfun];
            //
            //  Handle special case.
            //
            if (n == 1)
            {
                for (j = 0; j < nfun; j++)
                {
                    lfun[j] = 1.0;
                }

                return lfun;
            }

            llfun = FEM_1D_Lagrange.lagrange_value_OLD(n, x, nfun, xfun);

            for (j = 0; j < nfun; j++)
            {
                t = 0.0;
                for (i = 0; i < n; i++)
                {
                    t = t + Math.Abs(llfun[i + j * n]);
                }

                lfun[j] = t;
            }

            return lfun;
        }

        public static void lebesgue_plot(int n, double[] x, int nfun, double[] xfun,
                string label, string filename)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEBESGUE_PLOT plots the Lebesgue function for a set of points.
            //
            //  Discussion:
            //
            //    The interpolation interval is assumed to be [min(XFUN), max(XFUN)].
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 March 2014
            //
            //  Author:
            //
            //    John Burkardt.
            //
            //  Parameters:
            //
            //    Jean-Paul Berrut, Lloyd Trefethen,
            //    Barycentric Lagrange Interpolation,
            //    SIAM Review,
            //    Volume 46, Number 3, September 2004, pages 501-517.
            //
            //  Parameters:
            //
            //    Input, int N, the number of interpolation points.
            //
            //    Input, double X[N], the interpolation points.
            //
            //    Input, int NFUN, the number of evaluation points.
            //
            //    Input, double XFUN[NFUN], the evaluation points.  
            //
            //    Input, string LABEL, a title for the plot.
            //
            //    Input, string FILENAME, a partial filename.
            //    The program will create "filename_commands.txt', 'filename_data.txt',
            //    and 'filename.png'.
            //
        {
            string command_filename;
            List<string> command_unit = new List<string>();
            string data_filename;
            List<string> data_unit = new List<string>();
            int i;
            double[] lfun;
            string png_filename;

            lfun = lebesgue_function(n, x, nfun, xfun);
            //
            //  Create data file.
            //
            data_filename = filename + "_data.txt";
            for (i = 0; i < nfun; i++)
            {
                data_unit.Add(xfun[i] + "  "
                                      + lfun[i] + "");
            }

            File.WriteAllLines(data_filename, data_unit);
            Console.WriteLine("");
            Console.WriteLine("  Created graphics data file '" + data_filename + "'");
            //
            //  Create command file.
            //
            command_filename = filename + "_commands.txt";

            command_unit.Add("# " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("# Usage:");
            command_unit.Add("#  gnuplot < " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("set term png");

            png_filename = filename + ".png";
            command_unit.Add("set output '" + png_filename + "'");
            command_unit.Add("set xlabel '<--- X --->'");
            command_unit.Add("set ylabel '<--- Lebesgue(X) --->'");
            command_unit.Add("set title '" + label + "'");
            command_unit.Add("set grid");
            command_unit.Add("set style data lines");
            command_unit.Add("set timestamp");
            command_unit.Add("plot '" + data_filename
                                      + "' using 1:2 lw 3 linecolor rgb 'red'");

            File.WriteAllLines(command_filename, command_unit);
            Console.WriteLine("  Created graphics command file '" + command_filename + "'");

        }
    }
}