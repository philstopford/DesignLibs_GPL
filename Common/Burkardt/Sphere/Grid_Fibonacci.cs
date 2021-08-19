using System;
using System.Collections.Generic;
using System.IO;

namespace Burkardt.SphereNS
{
    public static class Grid_Fibonacci
    {
        public static void sphere_fibonacci_grid_display(int ng, double[] xg, string prefix )

        //*****************************************************************************/
        //
        //  Purpose:
        //
        //    SPHERE_FIBONACCI_GRID_DISPLAY displays sphere points on a Fibonacci spiral.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NG, the number of points.
        //
        //    Input, double XG[3*NG], the Fibonacci spiral points.
        //
        //    Input, string PREFIX, a prefix for the filenames.
        //
        {
            string command_filename;
            List<string> command_unit = new List<string>();
            string data_filename;
            List<string> data_unit= new List<string>();
            int i;
            int j;
            string plot_filename;
            //
            //  Create graphics data file.
            //
            data_filename = prefix + "_data.txt";

            for (j = 0; j < ng; j++)
            {
                string tmp = "";
                for (i = 0; i < 3; i++)
                {
                    tmp +="  " + xg[i + j * 3];
                }

                data_unit.Add(tmp);
            }

            File.WriteAllLines(data_filename, data_unit);

            Console.WriteLine("");
            Console.WriteLine("  Created data file '" + data_filename + "'");
            //
            //  Create graphics command file.
            //
            command_filename = prefix + "_commands.txt";

            command_unit.Add("# " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("# Usage:");
            command_unit.Add("#  gnuplot < " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("set term png");

            plot_filename = prefix + ".png";

            command_unit.Add("set output '" + plot_filename + "'");
            command_unit.Add("set xlabel '<--- X --->'");
            command_unit.Add("set ylabel '<--- Y --->'");
            command_unit.Add("set zlabel '<--- Z --->'");
            command_unit.Add("set title '" + prefix + "'");
            command_unit.Add("set grid");
            command_unit.Add("set key off");
            command_unit.Add("set style data points");
            command_unit.Add("set timestamp");
            command_unit.Add("set view equal xyz");
            command_unit.Add("splot '" + data_filename + "'");
            command_unit.Add("quit");

            File.WriteAllLines(command_filename, command_unit);
            
            Console.WriteLine("  Created command file '" + command_filename + "%s'");

            return;
        }

        public static double[] sphere_fibonacci_grid_points(int ng)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_FIBONACCI_GRID_POINTS computes sphere points on a Fibonacci spiral.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 May 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Richard Swinbank, James Purser,
            //    Fibonacci grids: A novel approach to global modelling,
            //    Quarterly Journal of the Royal Meteorological Society,
            //    Volume 132, Number 619, July 2006 Part B, pages 1769-1793.
            //
            //  Parameters:
            //
            //    Input, int NG, the number of points.
            //
            //    Output, double SPHERE_FIBONACCI_GRID_POINTS[3*NG], the Fibonacci 
            //    spiral points.
            //
        {
            double cphi;
            double i_r8;
            int j;
            double ng_r8;
            double r8_phi;
            const double r8_pi = 3.141592653589793;
            double sphi;
            double theta;
            double[] xyz;

            xyz = new double[3 * ng];

            r8_phi = (1.0 + Math.Sqrt(5.0)) / 2.0;
            ng_r8 = (double)(ng);

            for (j = 0; j < ng; j++)
            {
                i_r8 = (double)(-ng + 1 + 2 * j);
                theta = 2.0 * r8_pi * i_r8 / r8_phi;
                sphi = i_r8 / ng_r8;
                cphi = Math.Sqrt((ng_r8 + i_r8) * (ng_r8 - i_r8)) / ng_r8;
                xyz[0 + j * 3] = cphi * Math.Sin(theta);
                xyz[1 + j * 3] = cphi * Math.Cos(theta);
                xyz[2 + j * 3] = sphi;
            }

            return xyz;
        }
    }
}