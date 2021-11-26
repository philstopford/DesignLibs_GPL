using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using Burkardt.PolynomialNS;
using Burkardt.TetrahedronNS;
using Burkardt.Types;

namespace TetrahedronArbQTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TETRAHEDRON_ARBQ_RULE_TEST.
        //
        //  Discussion:
        //
        //    TETRAHEDRON_ARBQ_RULE_TEST tests the TETRAHEDRON_ARBQ_RULE library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU GPL license.
        //
        //  Modified:
        //
        //    11 July 2014
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
        //    C++ version by John Burkardt.
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
        Console.WriteLine("");
        Console.WriteLine("TETRAHEDRON_ARBQ_RULE_TEST");
        Console.WriteLine("  Test the TETRAHEDRON_ARBQ_RULE library.");

        int degree = 8;
        int n = QuadratureRule.tetrahedron_arbq_size(degree);
        string header = "tetrahedron08";

        test01(degree, n);

        test02(degree, n, header);

        test03(degree, n, header);

        test04(degree, n);
        Console.WriteLine("");
        Console.WriteLine("TETRAHEDRON_ARBQ_RULE_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01(int degree, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 calls TETRAHEDRON_ARBQ for a quadrature rule of given order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU GPL license.
        //
        //  Modified:
        //
        //    11 July 2014
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
        //    C++ version by John Burkardt.
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
        //    of the quadrature rule.  1 <= DEGREE <= 15.
        //
        //    Input, int N, the number of nodes.
        //
    {
        int j;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Quadrature rule for a tetrahedron.");
        Console.WriteLine("  Polynomial exactness degree DEGREE = " + degree + "");

        double volume = Math.Sqrt(8.0) / 3.0;
        //
        //  Retrieve and print a symmetric quadrature rule.
        //
        double[] x = new double[3 * n];
        double[] w = new double[n];

        QuadratureRule.tetrahedron_arbq(degree, n, ref x, ref w);

        Console.WriteLine("");
        Console.WriteLine("  Number of nodes N = " + n + "");

        Console.WriteLine("");
        Console.WriteLine("     J  W       X       Y      Z");
        Console.WriteLine("");
        for (j = 0; j < n; j++)
        {
            Console.WriteLine(j.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                                                      + w[j].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "  "
                                                      + x[0 + j * 3].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "  "
                                                      + x[1 + j * 3].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "  "
                                                      + x[2 + j * 3].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        double d = typeMethods.r8vec_sum(n, w);

        Console.WriteLine("   Sum    " + d + "");
        Console.WriteLine("  Volume  " + volume + "");
    }

    private static void test02(int degree, int n, string header)

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
        //    11 July 2014
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
        //    C++ version by John Burkardt.
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
        //    of the quadrature rule.  0 <= DEGREE <= 15.
        //
        //    Input, int N, the number of nodes to be used by the rule.
        //
        //    Input, string HEADER, an identifier for the filenames.
        //
    {
        int i;
        List<string> rule_unit = new();

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Get a quadrature rule for the tetrahedron.");
        Console.WriteLine("  Then write it to a file.");
        Console.WriteLine("  Polynomial exactness degree DEGREE = " + degree + "");
        //
        //  Retrieve a quadrature rule.
        //
        double[] x = new double[3 * n];
        double[] w = new double[n];

        QuadratureRule.tetrahedron_arbq(degree, n, ref x, ref w);
        //
        //  Write the points and weights to a file.
        //
        string rule_filename = header + ".txt";

        for (i = 0; i < n; i++)
        {
            rule_unit.Add(x[0 + i * 3] + "  "
                                       + x[1 + i * 3] + "  "
                                       + x[2 + i * 3] + "  "
                                       + w[i] + "");
        }

        File.WriteAllLines(rule_filename, rule_unit);
        Console.WriteLine("");
        Console.WriteLine("  Quadrature rule written to file '" + rule_filename + "'");
    }

    private static void test03(int degree, int n, string header)

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
        //    11 July 2014
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
        //    C++ version by John Burkardt.
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
        //    of the quadrature rule.  0 <= DEGREE <= 15.
        //
        //    Input, int N, the number of nodes to be used by the rule.
        //
        //    Input, string HEADER, an identifier for the filenames.
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  Get a quadrature rule for the tetrahedron.");
        Console.WriteLine("  Set up GNUPLOT graphics input.");
        Console.WriteLine("  Polynomial exactness degree DEGREE = " + degree + "");
        //
        //  Retrieve a symmetric quadrature rule.
        //
        double[] x = new double[3 * n];
        double[] w = new double[n];

        QuadratureRule.tetrahedron_arbq(degree, n, ref x, ref w);
        //
        //  Create files for input to GNUPLOT.
        //
        QuadratureRule.tetrahedron_arbq_gnuplot(n, x, header);
    }

    private static void test04(int degree, int n)

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
        //    11 July 2014
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
        //    C++ version by John Burkardt.
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
        //    of the quadrature rule.  0 <= DEGREE <= 15.
        //
        //    Input, int N, the number of nodes to be used by the rule.
        //
    {
        int i;
        int j;
        double[] z = new double[3];

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  Get a quadrature rule for the tetrahedron.");
        Console.WriteLine("  Test its accuracy.");
        Console.WriteLine("  Polynomial exactness degree DEGREE = " + degree + "");
        //
        //  Retrieve a quadrature rule.
        //
        double[] x = new double[3 * n];
        double[] w = new double[n];

        QuadratureRule.tetrahedron_arbq(degree, n, ref x, ref w);

        int npols = (degree + 1) * (degree + 2) * (degree + 3) / 6;
        double[] rints = new double[npols];

        for (j = 0; j < npols; j++)
        {
            rints[j] = 0.0;
        }

        for (i = 0; i < n; i++)
        {
            z[0] = x[0 + i * 3];
            z[1] = x[1 + i * 3];
            z[2] = x[2 + i * 3];

            double[] pols = Polynomial.ortho3eva(degree, z);

            for (j = 0; j < npols; j++)
            {
                rints[j] += w[i] * pols[j];
            }
        }

        double volume = Math.Sqrt(8.0) / 3.0;

        double d = 0.0;
        d = Math.Pow(rints[0] - Math.Sqrt(volume), 2);
        for (i = 1; i < npols; i++)
        {
            d += Math.Pow(rints[i], 2);
        }

        d = Math.Sqrt(d) / npols;

        Console.WriteLine("");
        Console.WriteLine("  RMS error = " + d + "");

    }
}