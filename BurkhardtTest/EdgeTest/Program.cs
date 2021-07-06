using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Edge;
using Burkardt.Types;

namespace EdgeTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EDGE_TEST tests the EDGE library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    21 September 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("EDGE_TEST");
            Console.WriteLine("  Test the EDGE library.");

            test01();
            test02();
            test03();
            test035();
            test036();
            test037();
            test04();
            Console.WriteLine("");
            Console.WriteLine("EDGE_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 plots functions with jump discontinuities.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 September 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Rick Archibald, Anne Gelb, Jungho Yoon,
            //    Polynomial fitting for edge detection in irregularly sampled signals 
            //    and images,
            //    SIAM Journal on Numerical Analysis,
            //    Volume 43, Number 1, 2006, pages 259-279.
            //
        {
            string command_filename;
            List<string> command_unit = new List<string>();
            string data_filename;
            List<string> data_unit = new List<string>();
            double[] f = new double[1];
            string header = "";
            int i;
            int n = 0;
            int test;
            int test_num;
            string title = "";
            double[] x = new double[1];
            double x_max;
            double x_min;

            Console.WriteLine("");
            Console.WriteLine("TEST01:");
            Console.WriteLine("  Plot 1D test functions.");

            test_num = 7;

            for (test = 1; test <= test_num; test++)
            {
                if (test == 1)
                {
                    n = 101;
                    x_min = -1.0;
                    x_max = +1.0;
                    x = typeMethods.r8vec_linspace_new(n, x_min, x_max);
                    header = "fx1";
                    f = OneD.fx1_vec(n, x);
                    title = "1D Test Function #1";
                }
                else if (test == 2)
                {
                    n = 101;
                    x_min = -1.0;
                    x_max = +1.0;
                    x = typeMethods.r8vec_linspace_new(n, x_min, x_max);
                    header = "fx2";
                    f = OneD.fx2_vec(n, x);
                    title = "1D Test Function #2";
                }
                else if (test == 3)
                {
                    n = 101;
                    x_min = -1.0;
                    x_max = +1.0;
                    x = typeMethods.r8vec_linspace_new(n, x_min, x_max);
                    header = "fx3";
                    f = OneD.fx3_vec(n, x);
                    title = "1D Test Function #3";
                }
                else if (test == 4)
                {
                    n = 101;
                    x_min = 0.0;
                    x_max = +1.0;
                    x = typeMethods.r8vec_linspace_new(n, x_min, x_max);
                    header = "fx4";
                    f = OneD.fx4_vec(n, x);
                    title = "1D Test Function #4";
                }
                else if (test == 5)
                {
                    n = 101;
                    x_min = -1.0;
                    x_max = +1.0;
                    x = typeMethods.r8vec_linspace_new(n, x_min, x_max);
                    header = "fx5";
                    f = OneD.fx5_vec(n, x);
                    title = "1D Test Function #5";
                }
                else if (test == 6)
                {
                    n = 101;
                    x_min = 0.0;
                    x_max = +1.0;
                    x = typeMethods.r8vec_linspace_new(n, x_min, x_max);
                    header = "fx6";
                    f = OneD.fx6_vec(n, x);
                    title = "1D Test Function #6";
                }
                else if (test == 7)
                {
                    n = 101;
                    x_min = 0.0;
                    x_max = +1.0;
                    x = typeMethods.r8vec_linspace_new(n, x_min, x_max);
                    header = "fx7";
                    f = OneD.fx7_vec(n, x);
                    title = "1D Test Function #7";
                }

                data_filename = header + "_data.txt";
                for (i = 0; i < n; i++)
                {
                    data_unit.Add(x[i]
                                  + "  " + f[i] + "");
                }

                File.WriteAllLines(data_filename, data_unit);
                Console.WriteLine("  Created data file '" + data_filename + "'");

                command_filename = header + "_commands.txt";
                command_unit.Add("# " + command_filename + "");
                command_unit.Add("#");
                command_unit.Add("# Usage:");
                command_unit.Add("#  gnuplot < " + command_filename + "");
                command_unit.Add("#");
                command_unit.Add("set term png");
                command_unit.Add("set output '" + header + ".png'");
                command_unit.Add("set xlabel '<--- X --->'");
                command_unit.Add("set ylabel '<--- Y --->'");
                command_unit.Add("set title '" + title + "'");
                command_unit.Add("set grid");
                command_unit.Add("set style data lines");
                command_unit.Add("plot '" + data_filename
                                          + "' using 1:2 with points lt 3 pt 4 linecolor rgb 'blue'");
                command_unit.Add("quit");

                File.WriteAllLines(command_filename, command_unit);
                Console.WriteLine("  Created command file '" + command_filename + "'");
            }
        }

        static void test02()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 plots a function with a jump discontinuity along a circle.
            //
            //  Discussion:
            //
            //    This is example 4.1 in the reference.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 February 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Rick Archibald, Anne Gelb, Jungho Yoon,
            //    Polynomial fitting for edge detection in irregularly sampled signals 
            //    and images,
            //    SIAM Journal on Numerical Analysis,
            //    Volume 43, Number 1, 2006, pages 259-279.
            //
        {
            string command_filename;
            List<string> command_unit = new List<string>();
            string data_filename;
            List<string> data_unit = new List<string>();
            double fxy;
            string header;
            int i;
            int j;
            int n;
            string title;
            double[] x;
            double x_max;
            double x_min;
            double[] y;
            double y_max;
            double y_min;

            Console.WriteLine("");
            Console.WriteLine("TEST02:");
            Console.WriteLine("  Plot 2D test function #1 with jump along circle.");

            header = "fxy1";
            title = "2D test function #1 with discontinuity along circle";

            n = 101;
            x_min = -1.0;
            x_max = +1.0;
            x = typeMethods.r8vec_linspace_new(n, x_min, x_max);
            y_min = -1.0;
            y_max = +1.0;
            y = typeMethods.r8vec_linspace_new(n, y_min, y_max);

            data_filename = header + "_data.txt";
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    fxy = TwoD.fxy1(x[i], y[j]);
                    data_unit.Add(x[i]
                                  + "  " + y[j]
                                  + "  " + fxy + "");
                }

                data_unit.Add("");
            }

            File.WriteAllLines(data_filename, data_unit);
            Console.WriteLine("  Created data file '" + data_filename + "'");

            command_filename = header + "_commands.txt";
            command_unit.Add("# " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("# Usage:");
            command_unit.Add("#  gnuplot < " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("set term png");
            command_unit.Add("set output '" + header + ".png'");
            command_unit.Add("set view 120, 77");
            command_unit.Add("set hidden3d");
            command_unit.Add("set timestamp");
            command_unit.Add("set xlabel '<--- X --->'");
            command_unit.Add("set ylabel '<--- Y --->'");
            command_unit.Add("set zlabel '<--- Z --->'");
            command_unit.Add("set title '" + title + "'");
            command_unit.Add("set grid");
            command_unit.Add("set style data lines");
            command_unit.Add("splot '" + data_filename + "' with lines");
            command_unit.Add("quit");
            File.WriteAllLines(command_filename, command_unit);
            Console.WriteLine("  Created command file '" + command_filename + "'");
        }

        static void test03()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST03 plots a function with a jump discontinuity along a circle.
            //
            //  Discussion:
            //
            //    This is example 4.2 in the reference.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 February 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Rick Archibald, Anne Gelb, Jungho Yoon,
            //    Polynomial fitting for edge detection in irregularly sampled signals 
            //    and images,
            //    SIAM Journal on Numerical Analysis,
            //    Volume 43, Number 1, 2006, pages 259-279.
            //
        {
            string command_filename;
            List<string> command_unit = new List<string>();
            string data_filename;
            List<string> data_unit = new List<string>();
            double fxy;
            string header;
            int i;
            int j;
            int n;
            string title;
            double[] x;
            double x_max;
            double x_min;
            double[] y;
            double y_max;
            double y_min;

            Console.WriteLine("");
            Console.WriteLine("TEST03:");
            Console.WriteLine("  Plot 2D test function #2, the Shepp Logan phantom.");

            header = "fxy2";
            title = "2D test function #2, the Shepp Logan phantom";

            n = 101;
            x_min = -1.0;
            x_max = +1.0;
            x = typeMethods.r8vec_linspace_new(n, x_min, x_max);
            y_min = -1.0;
            y_max = +1.0;
            y = typeMethods.r8vec_linspace_new(n, y_min, y_max);

            data_filename = header + "_data.txt";
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    fxy = TwoD.fxy2(x[i], y[j]);
                    data_unit.Add(x[i]
                                  + "  " + y[j]
                                  + "  " + fxy + "");
                }

                data_unit.Add("");
            }

            File.WriteAllLines(data_filename, data_unit);
            Console.WriteLine("  Created data file '" + data_filename + "'");

            command_filename = header + "_commands.txt";
            command_unit.Add("# " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("# Usage:");
            command_unit.Add("#  gnuplot < " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("set term png");
            command_unit.Add("set output '" + header + ".png'");
            command_unit.Add("set view 30, 75");
            command_unit.Add("set hidden3d");
            command_unit.Add("set timestamp");
            command_unit.Add("set xlabel '<--- X --->'");
            command_unit.Add("set ylabel '<--- Y --->'");
            command_unit.Add("set zlabel '<--- Z --->'");
            command_unit.Add("set title '" + title + "'");
            command_unit.Add("set grid");
            command_unit.Add("set style data lines");
            command_unit.Add("splot '" + data_filename + "' with lines");
            command_unit.Add("quit");
            File.WriteAllLines(command_filename, command_unit);
            Console.WriteLine("  Created command file '" + command_filename + "'");
        }

        static void test035()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST035 plots a function with a jump discontinuity along a circle.
            //
            //  Discussion:
            //
            //    This is example 3.2 in the reference.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 September 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Rick Archibald, Anne Gelb, Jungho Yoon,
            //    Determining the location of discontinuities in the derivatives
            //    of functions,
            //    Applied Numerical Mathematics,
            //    Volume 58, 2008, pages 577-592.
            //
        {
            string command_filename;
            List<string> command_unit = new List<string>();
            string data_filename;
            List<string> data_unit = new List<string>();
            double fxy;
            string header;
            int i;
            int j;
            int n;
            string title;
            double[] x;
            double x_max;
            double x_min;
            double[] y;
            double y_max;
            double y_min;

            Console.WriteLine("");
            Console.WriteLine("TEST035:");
            Console.WriteLine("  Plot 2D test function #3, the modified 2D Harten function.");

            header = "fxy3";
            title = "2D test function #3, the modified 2D Harten function";

            n = 101;
            x_min = -1.0;
            x_max = +1.0;
            x = typeMethods.r8vec_linspace_new(n, x_min, x_max);
            y_min = -1.0;
            y_max = +1.0;
            y = typeMethods.r8vec_linspace_new(n, y_min, y_max);

            data_filename = header + "_data.txt";
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    fxy = TwoD.fxy3(x[i], y[j]);
                    data_unit.Add(x[i]
                                  + "  " + y[j]
                                  + "  " + fxy + "");
                }

                data_unit.Add("");
            }

            File.WriteAllLines(data_filename, data_unit);
            Console.WriteLine("  Created data file '" + data_filename + "'");

            command_filename = header + "_commands.txt";
            command_unit.Add("# " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("# Usage:");
            command_unit.Add("#  gnuplot < " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("set term png");
            command_unit.Add("set output '" + header + ".png'");
            command_unit.Add("set view 30, 75");
            command_unit.Add("set hidden3d");
            command_unit.Add("set timestamp");
            command_unit.Add("set xlabel '<--- X --->'");
            command_unit.Add("set ylabel '<--- Y --->'");
            command_unit.Add("set zlabel '<--- Z --->'");
            command_unit.Add("set title '" + title + "'");
            command_unit.Add("set grid");
            command_unit.Add("set style data lines");
            command_unit.Add("splot '" + data_filename + "' with lines");
            command_unit.Add("quit");
            File.WriteAllLines(command_filename, command_unit);
            Console.WriteLine("  Created command file '" + command_filename + "'");
        }

        static void test036()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST036 plots a function with a derivative discontinuity.
            //
            //  Discussion:
            //
            //    This is example 3.1 in the reference.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    21 September 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Rick Archibald, Anne Gelb, Jungho Yoon,
            //    Determining the location of discontinuities in the derivatives
            //    of functions,
            //    Applied Numerical Mathematics,
            //    Volume 58, 2008, pages 577-592.
            //
        {
            string command_filename;
            List<string> command_unit = new List<string>();
            string data_filename;
            List<string> data_unit = new List<string>();
            double fxy;
            string header;
            int i;
            int j;
            int n;
            string title;
            double[] x;
            double x_max;
            double x_min;
            double[] y;
            double y_max;
            double y_min;

            Console.WriteLine("");
            Console.WriteLine("TEST036:");
            Console.WriteLine("  Plot 2D test function #4, the discontinuous medium wave, P(x,t).");

            header = "fxy4";
            title = "2D test function #4, the discontinuous medium wave, P(x,t)";

            n = 101;
            x_min = -1.0;
            x_max = 0.0;
            x = typeMethods.r8vec_linspace_new(n, x_min, x_max);
            y_min = 0.0;
            y_max = 0.1;
            y = typeMethods.r8vec_linspace_new(n, y_min, y_max);

            data_filename = header + "_data.txt";
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    fxy = TwoD.fxy4(x[i], y[j]);
                    data_unit.Add(x[i]
                                  + "  " + y[j]
                                  + "  " + fxy + "");
                }

                data_unit.Add("");
            }

            File.WriteAllLines(data_filename, data_unit);
            Console.WriteLine("  Created data file '" + data_filename + "'");

            command_filename = header + "_commands.txt";
            command_unit.Add("# " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("# Usage:");
            command_unit.Add("#  gnuplot < " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("set term png");
            command_unit.Add("set output '" + header + ".png'");
            command_unit.Add("set view 30, 45");
            command_unit.Add("set hidden3d");
            command_unit.Add("set timestamp");
            command_unit.Add("set xlabel '<--- X --->'");
            command_unit.Add("set ylabel '<--- Y --->'");
            command_unit.Add("set zlabel '<--- Z --->'");
            command_unit.Add("set title '" + title + "'");
            command_unit.Add("set grid");
            command_unit.Add("set style data lines");
            command_unit.Add("splot '" + data_filename + "' with lines");
            command_unit.Add("quit");
            File.WriteAllLines(command_filename, command_unit);
            Console.WriteLine("  Created command file '" + command_filename + "'");
        }

        static void test037()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST037 plots a function with a derivative discontinuity.
            //
            //  Discussion:
            //
            //    This is example 3.1 in the reference.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    21 September 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Rick Archibald, Anne Gelb, Jungho Yoon,
            //    Determining the location of discontinuities in the derivatives
            //    of functions,
            //    Applied Numerical Mathematics,
            //    Volume 58, 2008, pages 577-592.
            //
        {
            string command_filename;
            List<string> command_unit = new List<string>();
            string data_filename;
            List<string> data_unit = new List<string>();
            double fxy;
            string header;
            int i;
            int j;
            int n;
            string title;
            double[] x;
            double x_max;
            double x_min;
            double[] y;
            double y_max;
            double y_min;

            Console.WriteLine("");
            Console.WriteLine("TEST037:");
            Console.WriteLine("  Plot 2D test function #5, the discontinuous medium wave, U(x,t).");

            header = "fxy5";
            title = "2D test function #5, the discontinuous medium wave, U(x,t)";

            n = 101;
            x_min = -1.0;
            x_max = 0.0;
            x = typeMethods.r8vec_linspace_new(n, x_min, x_max);
            y_min = 0.0;
            y_max = 0.1;
            y = typeMethods.r8vec_linspace_new(n, y_min, y_max);

            data_filename = header + "_data.txt";
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    fxy = TwoD.fxy5(x[i], y[j]);
                    data_unit.Add(x[i]
                                  + "  " + y[j]
                                  + "  " + fxy + "");
                }

                data_unit.Add("");
            }

            File.WriteAllLines(data_filename, data_unit);
            Console.WriteLine("  Created data file '" + data_filename + "'");

            command_filename = header + "_commands.txt";
            command_unit.Add("# " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("# Usage:");
            command_unit.Add("#  gnuplot < " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("set term png");
            command_unit.Add("set output '" + header + ".png'");
            command_unit.Add("set view 30, 45");
            command_unit.Add("set hidden3d");
            command_unit.Add("set timestamp");
            command_unit.Add("set xlabel '<--- X --->'");
            command_unit.Add("set ylabel '<--- Y --->'");
            command_unit.Add("set zlabel '<--- Z --->'");
            command_unit.Add("set title '" + title + "'");
            command_unit.Add("set grid");
            command_unit.Add("set style data lines");
            command_unit.Add("splot '" + data_filename + "' with lines");
            command_unit.Add("quit");
            File.WriteAllLines(command_filename, command_unit);
            Console.WriteLine("  Created command file '" + command_filename + "'");
        }

        static void test04()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST04 plots slices of a 3D function.
            //
            //  Discussion:
            //
            //    Although the slice plots look uninteresting, there is a lot of detail
            //    hidden in the data in variations that are not obvious at first.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 February 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Larry Shepp,
            //    Computerized tomography and nuclear magnetic resonance,
            //    Journal of Computer Assisted Tomography,
            //    Volume 4, Number 1, February 1980, pages 94-107.
            //
        {
            string command_filename;
            List<string> command_unit = new List<string>();
            string data_filename;
            List<string> data_unit = new List<string>();
            double fxyz;
            string header = "";
            int i;
            int j;
            int n;
            int test;
            int test_num;
            string title = "";
            double[] x;
            double x_max;
            double x_min;
            double x_val = 0;
            double[] y;
            double y_max;
            double y_min;
            double y_val = 0;
            double[] z;
            double z_max;
            double z_min;
            double z_val = 0;

            Console.WriteLine("");
            Console.WriteLine("TEST04:");
            Console.WriteLine("  Plot 3D test function #1, the Shepp Logan 3D phantom.");

            test_num = 3;

            n = 101;
            x_min = -1.5;
            x_max = +1.5;
            x = typeMethods.r8vec_linspace_new(n, x_min, x_max);
            y_min = -1.5;
            y_max = +1.5;
            y = typeMethods.r8vec_linspace_new(n, y_min, y_max);
            z_min = -1.5;
            z_max = +1.5;
            z = typeMethods.r8vec_linspace_new(n, z_min, z_max);

            for (test = 1; test <= test_num; test++)
            {
                if (test == 1)
                {
                    x_val = 0.0;
                    title = "Slice X = 0.0";
                    header = "fxyz1_x";
                }
                else if (test == 2)
                {
                    y_val = 0.0;
                    title = "Slice Y = 0.0";
                    header = "fxyz1_y";
                }
                else if (test == 3)
                {
                    z_val = -0.1;
                    title = "Slice Z = - 0.1";
                    header = "fxyz1_z";
                }

                data_filename = header + "_data.txt";
                for (i = 0; i < n; i++)
                {
                    for (j = 0; j < n; j++)
                    {
                        if (test == 1)
                        {
                            fxyz = ThreeD.fxyz1(x_val, y[j], z[i]);
                            data_unit.Add(y[j]
                                          + "  " + z[i]
                                          + "  " + fxyz + "");
                        }
                        else if (test == 2)
                        {
                            fxyz = ThreeD.fxyz1(x[j], y_val, z[i]);
                            data_unit.Add(x[j]
                                          + "  " + z[i]
                                          + "  " + fxyz + "");
                        }
                        else if (test == 3)
                        {
                            fxyz = ThreeD.fxyz1(x[j], y[i], z_val);
                            data_unit.Add(x[j]
                                          + "  " + y[i]
                                          + "  " + fxyz + "");
                        }
                    }

                    data_unit.Add("");
                }

                File.WriteAllLines(data_filename, data_unit);
                Console.WriteLine("  Created data file '" + data_filename + "'");

                command_filename = header + "_commands.txt";
                command_unit.Add("# " + command_filename + "");
                command_unit.Add("#");
                command_unit.Add("# Usage:");
                command_unit.Add("#  gnuplot < " + command_filename + "");
                command_unit.Add("#");
                command_unit.Add("set term png");
                command_unit.Add("set output '" + header + ".png'");
                command_unit.Add("set view 20, 75");
                command_unit.Add("set hidden3d");
                command_unit.Add("set timestamp");
                if (test == 1)
                {
                    command_unit.Add("set xlabel '<--- Y --->'");
                    command_unit.Add("set ylabel '<--- Z --->'");
                    command_unit.Add("set zlabel '<--- X --->'");
                }
                else if (test == 2)
                {
                    command_unit.Add("set xlabel '<--- X --->'");
                    command_unit.Add("set ylabel '<--- Z --->'");
                    command_unit.Add("set zlabel '<--- Y --->'");
                }
                else if (test == 3)
                {
                    command_unit.Add("set xlabel '<--- X --->'");
                    command_unit.Add("set ylabel '<--- Y --->'");
                    command_unit.Add("set zlabel '<--- Z --->'");
                }

                command_unit.Add("set title '" + title + "'");
                command_unit.Add("set grid");
                command_unit.Add("set style data lines");
                command_unit.Add("splot '" + data_filename + "' with lines");
                command_unit.Add("quit");
                File.WriteAllLines(command_filename, command_unit);
                Console.WriteLine("  Created command file '" + command_filename + "'");
            }
        }
    }
}