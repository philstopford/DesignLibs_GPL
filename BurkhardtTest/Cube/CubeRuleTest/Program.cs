using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Cube;
using Burkardt.Types;

namespace CubeRuleTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for CUBE_ARBQ_RULE_TEST.
            //
            //  Discussion:
            //
            //    CUBE_ARBQ_RULE_TEST tests the CUBE_ARBQ_RULE library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU GPL license.
            //
            //  Modified:
            //
            //    09 July 2014
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
            //    This C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Hong Xiao, Zydrunas Gimbutas,
            //    A numerical algorithm for the construction of efficient quadrature
            //    rules in two and higher dimensions,
            //    Computers and Mathematics with Applications,
            //    Volume 59, 2010, pages 663-676.
            //
        {
            int degree;
            string header;
            int n;

            Console.WriteLine("");
            Console.WriteLine("CUBE_ARBQ_RULE_TEST");
            Console.WriteLine("  C version");
            Console.WriteLine("  Test the CUBE_ARBQ_RULE library.");

            degree = 8;
            n = QuadratureRule.cube_arbq_size(degree);
            header = "cube08";

            test01(degree, n);

            test02(degree, n, header);

            test03(degree, n, header);

            test04(degree, n);

            Console.WriteLine("");
            Console.WriteLine("CUBE_ARBQ_RULE_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01(int degree, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 calls CUBE_ARBQ for a quadrature rule of given order.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU GPL license.
            //
            //  Modified:
            //
            //    09 July 2014
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
            //    This C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Hong Xiao, Zydrunas Gimbutas,
            //    A numerical algorithm for the construction of efficient quadrature
            //    rules in two and higher dimensions,
            //    Computers and Mathematics with Applications,
            //    Volume 59, 2010, pages 663-676.
            //
            //  Parameters:
            //
            //    Input, int DEGREE, the desired total polynomial degree exactness
            //    of the quadrature rule.
            //
            //    Input, int N, the number of nodes.
            //
        {
            double d;
            int j;
            double volume;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  Symmetric quadrature rule for a cube.");
            Console.WriteLine("  Polynomial exactness degree DEGREE = " + degree + "");

            volume = 8.0;
            //
            //  Retrieve and print a symmetric quadrature rule.
            //
            x = new double[3 * n];
            w = new double[n];

            QuadratureRule.cube_arbq(degree, n, ref x, ref w);

            Console.WriteLine("");
            Console.WriteLine("  Number of nodes N = " + n + "");

            Console.WriteLine("");
            Console.WriteLine("     J  W       X       Y      Z");
            Console.WriteLine("");
            for (j = 0; j < n; j++)
            {
                Console.WriteLine(j.ToString().PadLeft(4) + "  "
                    + w[j].ToString().PadLeft(14) + "  "
                    + x[0 + j * 3].ToString().PadLeft(14) + "  "
                    + x[1 + j * 3].ToString().PadLeft(14) + "  "
                    + x[2 + j * 3].ToString().PadLeft(14) + "");
            }

            d = typeMethods.r8vec_sum(n, w);

            Console.WriteLine("   Sum    " + d + "");
            Console.WriteLine("  Volume  " + volume + "");
        }

        static void test02(int degree, int n, string header)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 gets a rule and writes it to a file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU GPL license.
            //
            //  Modified:
            //
            //    09 July 2014
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
            //    This C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Hong Xiao, Zydrunas Gimbutas,
            //    A numerical algorithm for the construction of efficient quadrature
            //    rules in two and higher dimensions,
            //    Computers and Mathematics with Applications,
            //    Volume 59, 2010, pages 663-676.
            //
            //  Parameters:
            //
            //    Input, int DEGREE, the desired total polynomial degree exactness
            //    of the quadrature rule.  0 <= DEGREE <= 50.
            //
            //    Input, int N, the number of nodes to be used by the rule.
            //
            //    Input, string HEADER, an identifier for the filenames.
            //
        {
            int i;
            List<string> rule_unit = new List<string>();
            string rule_filename;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  Get a quadrature rule for the symmetric cube.");
            Console.WriteLine("  Then write it to a file.");
            Console.WriteLine("  Polynomial exactness degree DEGREE = " + degree + "");
            //
            //  Retrieve a symmetric quadrature rule.
            //
            x = new double[3 * n];
            w = new double[n];

            QuadratureRule.cube_arbq(degree, n, ref x, ref w);
            //
            //  Write the points and weights to a file.
            //
            rule_filename = header + ".txt";

            for (i = 0; i < n; i++)
            {
                rule_unit.Add( + x[0 + i * 3] + "  "
                    + x[1 + i * 3] + "  "
                    + x[2 + i * 3] + "  "
                    + w[i] + "");
            }

            File.WriteAllLines(rule_filename, rule_unit);
            Console.WriteLine("");
            Console.WriteLine("  Quadrature rule written to file '" + rule_filename + "'");
        }

        static void test03(int degree, int n, string header)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST03 gets a rule and creates GNUPLOT input files.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU GPL license.
            //
            //  Modified:
            //
            //    09 July 2014
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
            //    This C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Hong Xiao, Zydrunas Gimbutas,
            //    A numerical algorithm for the construction of efficient quadrature
            //    rules in two and higher dimensions,
            //    Computers and Mathematics with Applications,
            //    Volume 59, 2010, pages 663-676.
            //
            //  Parameters:
            //
            //    Input, int DEGREE, the desired total polynomial degree exactness
            //    of the quadrature rule.  0 <= DEGREE <= 50.
            //
            //    Input, int N, the number of nodes to be used by the rule.
            //
            //    Input, string HEADER, an identifier for the filenames.
            //
        {
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST03");
            Console.WriteLine("  Get a quadrature rule for the symmetric cube.");
            Console.WriteLine("  Set up GNUPLOT graphics input.");
            Console.WriteLine("  Polynomial exactness degree DEGREE = " + degree + "");
            //
            //  Retrieve a symmetric quadrature rule.
            //
            x = new double[3 * n];
            w = new double[n];

            QuadratureRule.cube_arbq(degree, n, ref x, ref w);
            //
            //  Create files for input to GNUPLOT.
            //
            QuadratureRule.cube_arbq_gnuplot(n, x, header);
        }

        static void test04(int degree, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST04 gets a rule and tests its accuracy.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU GPL license.
            //
            //  Modified:
            //
            //    09 July 2014
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
            //    This C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Hong Xiao, Zydrunas Gimbutas,
            //    A numerical algorithm for the construction of efficient quadrature
            //    rules in two and higher dimensions,
            //    Computers and Mathematics with Applications,
            //    Volume 59, 2010, pages 663-676.
            //
            //  Parameters:
            //
            //    Input, int DEGREE, the desired total polynomial degree exactness
            //    of the quadrature rule.  0 <= DEGREE <= 50.
            //
            //    Input, int N, the number of nodes to be used by the rule.
            //
        {
            double d;
            int i;
            int j;
            int npols;
            double[] pols;
            double[] rints;
            double volume;
            double[] w;
            double[] x;
            double[] z = new double[3];

            Console.WriteLine("");
            Console.WriteLine("TEST04");
            Console.WriteLine("  Get a quadrature rule for the symmetric cube.");
            Console.WriteLine("  Test its accuracy.");
            Console.WriteLine("  Polynomial exactness degree DEGREE = " + degree + "");
            //
            //  Retrieve a symmetric quadrature rule.
            //
            x = new double[3 * n];
            w = new double[n];

            QuadratureRule.cube_arbq(degree, n, ref x, ref w);

            npols = ((degree + 1) * (degree + 2) * (degree + 3)) / 6;
            rints = new double[npols];

            for (j = 0; j < npols; j++)
            {
                rints[j] = 0.0;
            }

            for (i = 0; i < n; i++)
            {
                z[0] = x[0 + i * 3];
                z[1] = x[1 + i * 3];
                z[2] = x[2 + i * 3];

                pols = QuadratureRule.lege3eva(degree, z);
                for (j = 0; j < npols; j++)
                {
                    rints[j] = rints[j] + w[i] * pols[j];
                }
            }

            volume = 8.0;

            d = 0.0;
            d = Math.Pow(rints[0] - Math.Sqrt(volume), 2);
            for (i = 1; i < npols; i++)
            {
                d = d + Math.Pow(rints[i], 2);
            }

            d = Math.Sqrt(d) / (double) (npols);

            Console.WriteLine("");
            Console.WriteLine("  RMS error = " + d + "");
        }
    }
}