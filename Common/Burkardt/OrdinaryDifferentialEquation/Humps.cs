using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Types;

namespace Burkardt
{
    public static class Humps
    {
        public static double humps_antideriv(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    humps_antideriv evaluates the antiderivative of the humps function.
            //
            //  Discussion:
            //
            //    y = 1.0 / ( r8_square ( x - 0.3 ) + 0.01 ) 
            //      + 1.0 / ( r8_square ( x - 0.9 ) + 0.04 ) 
            //      - 6.0
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 August 2019
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Input:
            //
            //    double x: the argument.
            //
            //  Output:
            //
            //    double humps_antideriv: the value of the antiderivative at x.
            //
        {
            double ya;

            ya = (1.0 / 0.1) * Math.Atan((x - 0.3) / 0.1)
                 + (1.0 / 0.2) * Math.Atan((x - 0.9) / 0.2)
                 - 6.0 * x;

            return ya;
        }

        public static double humps_deriv(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    humps_deriv evaluates the derivative of the humps function.
            //
            //  Discussion:
            //
            //    y = 1.0 / ( ( x - 0.3 )^2 + 0.01 )
            //      + 1.0 / ( ( x - 0.9 )^2 + 0.04 )
            //      - 6.0
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 August 2019
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Input:
            //
            //    double x: the argument.
            //
            //  Output:
            //
            //    double humps_deriv: the value of the derivative at x.
            //
        {
            double yp;

            yp = -2.0 * (x - 0.3) / typeMethods.r8_square(typeMethods.r8_square(x - 0.3) + 0.01)
                 - 2.0 * (x - 0.9) / typeMethods.r8_square(typeMethods.r8_square(x - 0.9) + 0.04);

            return yp;
        }

        public static double humps_deriv2(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    humps_deriv2 evaluates the second derivative of the humps function.
            //
            //  Discussion:
            //
            //    y = 1.0 / ( ( x - 0.3 )^2 + 0.01 )
            //      + 1.0 / ( ( x - 0.9 )^2 + 0.04 )
            //      - 6.0;
            //
            //    yp = - 2.0 * ( x - 0.3 ) / ( ( x - 0.3 )^2 + 0.01 )^2
            //         - 2.0 * ( x - 0.9 ) / ( ( x - 0.9 )^2 + 0.04 )^2;
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 August 2019
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Input:
            //
            //    double x: the argument.
            //
            //  Output:
            //
            //    double humps_deriv2: the value of the second derivative at x.
            //
        {
            double u1;
            double u1p;
            double u2;
            double u2p;
            double v1;
            double v1p;
            double v2;
            double v2p;
            double ypp;

            u1 = -2.0 * (x - 0.3);
            v1 = typeMethods.r8_square(typeMethods.r8_square(x - 0.3) + 0.01);
            u2 = -2.0 * (x - 0.9);
            v2 = typeMethods.r8_square(typeMethods.r8_square(x - 0.9) + 0.04);

            u1p = -2.0;
            v1p = 2.0 * (typeMethods.r8_square(x - 0.3) + 0.01) * 2.0 * (x - 0.3);
            u2p = -2.0;
            v2p = 2.0 * (typeMethods.r8_square(x - 0.9) + 0.04) * 2.0 * (x - 0.9);

            ypp = (u1p * v1 - u1 * v1p) / v1 / v1
                  + (u2p * v2 - u2 * v2p) / v2 / v2;

            return ypp;
        }

        public static double humps_fun(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    humps_fun evaluates the humps function.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 August 2019
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Input:
            //
            //    double x: the evaluation point.
            //
            //  Output:
            //
            //    double y: the function value.
            //
        {
            double y;

            y = 1.0 / (typeMethods.r8_square(x - 0.3) + 0.01)
                + 1.0 / (typeMethods.r8_square(x - 0.9) + 0.04)
                - 6.0;

            return y;
        }

        public static double humps_ode(double x, double y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    humps_ode evaluates the derivative of the humps function for an ODE solver.
            //
            //  Discussion:
            //
            //    This verion of "humps_deriv" appends the input argument "y", as expected 
            //    by most ODE solving software.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 August 2019
            //
            //  Author:
            //
            //    John Burkardt
            //
            //
            //  Input:
            //
            //    double x: the argument.
            //
            //    double y: the value of the dependent variable.
            //
            //  Output:
            //
            //    double humps_ode: the value of the derivative of the humps function.
            //
        {
            double yp;

            yp = -1.0 / typeMethods.r8_square(typeMethods.r8_square(x - 0.3) + 0.01)
                 * 2.0 * (x - 0.3)
                 - 1.0 / typeMethods.r8_square(typeMethods.r8_square(x - 0.9) + 0.04)
                 * 2.0 * (x - 0.9);

            return yp;
        }
        
        public static void plot_xy ( int n, double[] x, double[] y, string prefix )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    plot_xy plots xy data.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 August 2019
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Input:
        //
        //    int n, the number of data points.
        //
        //    double x[n], y[n], the data points.
        //
        //    string prefix, the prefix for the plot names.
        //
        {
            string command_filename;
            List<string> command_unit = new List<string>();
            string data_filename;
            List<string> data_unit = new List<string>();
            int i;
            string output_filename;
            string prefix2;
            //
            //  Create the data file.
            //
            data_filename = prefix + "_data.txt";
            for ( i = 0; i < n; i++ )
            {
                data_unit.Add(x[i] + "  " + y[i] + "");
            }
            File.WriteAllLines(data_filename, data_unit);
            Console.WriteLine("");
            Console.WriteLine("  Created graphics data file '" + data_filename + "'.");
            //
            //  Plot the selected data.
            //
            command_filename = prefix + "_commands.txt";

            command_unit.Add("# " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("# Usage:");
            command_unit.Add("#  gnuplot < " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("set term png");
            command_unit.Add("set nokey");
            output_filename = prefix + ".png";
            command_unit.Add("set output '" + output_filename + "'");
            command_unit.Add("set xlabel '<---X--->'");
            command_unit.Add("set ylabel '<---Y(X)--->'");
            prefix2 = typeMethods.s_escape_tex ( prefix );
            command_unit.Add("set title '" + prefix2 + "'");
            command_unit.Add("set grid");
            command_unit.Add("set style data lines");
            command_unit.Add("plot '" + data_filename + "' using 1:2 lw 3 linecolor rgb 'blue'");

            File.WriteAllLines(command_filename, command_unit);
            Console.WriteLine("  Created graphics command file '" + command_filename+ "'.");
        }
    }
}