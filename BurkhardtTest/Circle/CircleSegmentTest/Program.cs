﻿using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using Burkardt.CircleNS;
using Burkardt.IO;
using Burkardt.PolynomialNS;
using Burkardt.Quadrature;
using Burkardt.Types;
using Burkardt.Uniform;

namespace CircleSegmentTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for CIRCLE_SEGMENT_TEST.
        //
        //  Discussion:
        //
        //    CIRCLE_SEGMENT_TEST tests the CIRCLE_SEGMENT library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("CIRCLE_SEGMENT_TEST");

        Console.WriteLine("  Test the CIRCLE_SEGMENT library.");

        test01();
        test05();
        test06();
        test07();
        test08();
        test09();
        test11();
        test13();
        test14();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("CIRCLE_SEGMENT_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests CIRCLE_SEGMENT_AREA_FROM_HEIGHT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  CIRCLE_SEGMENT_AREA_FROM_HEIGHT computes the area of a circle segment.");

        Console.WriteLine("");
        Console.WriteLine("          R               H               Area");
        Console.WriteLine("");
        const double r = 1.0;
        double h = 1.0;
        for (i = 0; i <= 10; i++)
        {
            double area = Segment.circle_segment_area_from_height(r, h);
            Console.WriteLine("  " + r.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + h.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + area.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            h /= 2.0;
        }
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests the AREA and HEIGHT functions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double h;
        double r;
        int test;

        Console.WriteLine("");
        Console.WriteLine("CIRCLE_SEGMENT_TEST05");
        Console.WriteLine("  For circle segment with a given radius R,");
        Console.WriteLine("  CIRCLE_SEGMENT_AREA_FROM_HEIGHT computes the area A, given the height.");
        Console.WriteLine("  CIRCLE_SEGMENT_HEIGHT_FROM_AREA computes height H, given the area.");
        Console.WriteLine("  Check that these functions are inverses of each other");
        Console.WriteLine("  using random values of R, A, and H.");

        Console.WriteLine("");
        Console.WriteLine("        R             H      =>     A    =>       H2");
        Console.WriteLine("");

        int seed = 123456789;

        for (test = 1; test <= 5; test++)
        {
            r = 5.0 * UniformRNG.r8_uniform_01(ref seed);
            h = 2.0 * r * UniformRNG.r8_uniform_01(ref seed);
            a = Segment.circle_segment_area_from_height(r, h);
            double h2 = Segment.circle_segment_height_from_area(r, a);
            Console.WriteLine("  " + r.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + h.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + a.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + h2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("        R             A      =>     H    =>       A2");
        Console.WriteLine("");

        for (test = 1; test <= 5; test++)
        {
            r = 5.0 * UniformRNG.r8_uniform_01(ref seed);
            a = Math.PI * r * r * UniformRNG.r8_uniform_01(ref seed);
            h = Segment.circle_segment_height_from_area(r, a);
            double a2 = Segment.circle_segment_area_from_height(r, h);
            Console.WriteLine("  " + r.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + a.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + h.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + a2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 samples using CIRCLE_SEGMENT_SAMPLE_FROM_HEIGHT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int an_num = 51;
        string boundary_filename = "sample00_boundary.txt";
        List<string> boundary_unit = new();
        string command_filename = "sample00_commands.txt";
        List<string> command_unit = new();
        string data_filename = "sample00_data.txt";
        List<string> data_unit = new();
        const int data_num = 100;
        string graphics_filename = "sample00.png";
        int test;

        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("CIRCLE_SEGMENT_TEST06");
        Console.WriteLine("  CIRCLE_SEGMENT_SAMPLE_FROM_HEIGHT samples a circle segment.");
        Console.WriteLine("");
        Console.WriteLine("  Plot " + data_num + " points from several segments.");
        Console.WriteLine("");

        const double r = 1.0;
        double theta = Math.PI;

        for (test = 1; test <= 4; test++)
        {
            double h = Segment.circle_segment_height_from_angle(r, theta);

            double thetah = theta / 2.0;
            //
            //  Create boundary.
            //
            double[] an = typeMethods.r8vec_linspace_new(an_num, -thetah, +thetah);
            int i;
            for (i = 0; i < an_num; i++)
            {
                an[i] += 0.5 * Math.PI;
            }

            double[] boundary_x = new double[an_num + 1];
            double[] boundary_y = new double[an_num + 1];

            for (i = 0; i < an_num; i++)
            {
                boundary_x[i] = r * Math.Cos(an[i]);
                boundary_y[i] = r * Math.Sin(an[i]);
            }

            boundary_x[an_num] = boundary_x[0];
            boundary_y[an_num] = boundary_y[0];

            Files.filename_inc(ref boundary_filename);
            for (i = 0; i <= an_num; i++)
            {
                boundary_unit.Add("  " + boundary_x[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + boundary_y[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }

            File.WriteAllLines(boundary_filename, boundary_unit);
            Console.WriteLine("");
            Console.WriteLine("  Created boundary file \"" + boundary_filename + "\".");
            //
            //  Create data.
            //
            double[] data_x = new double[data_num + 1];
            double[] data_y = new double[data_num + 1];

            Segment.circle_segment_sample_from_height(r, h, data_num, ref seed, ref data_x, ref data_y);

            Files.filename_inc(ref data_filename);
            for (i = 0; i < data_num; i++)
            {
                data_unit.Add( "  " + data_x[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                    + "  " + data_y[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }

            File.WriteAllLines(data_filename, data_unit);
                
            Console.WriteLine("");
            Console.WriteLine("  Created data file \"" + data_filename + "\".");
            //
            //  Create commands.
            //
            Files.filename_inc(ref command_filename);
            command_unit.Add("# " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("# Usage:");
            command_unit.Add("#  gnuplot < " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("set term png");
            Files.filename_inc(ref graphics_filename);
            command_unit.Add("set output '" + graphics_filename + "'");
            command_unit.Add("set xlabel '<--- X --->'");
            command_unit.Add("set ylabel '<--- Y --->'");
            command_unit.Add("set title 'Circle Segment Sample'");
            command_unit.Add("set grid");
            command_unit.Add("set key off");
            command_unit.Add("set size ratio -1");
            command_unit.Add("set style data lines");
            command_unit.Add("plot '" + data_filename
                                      + "' using 1:2 with points lt 3 pt 3,\\");
            command_unit.Add("    '" + boundary_filename
                                     + "' using 1:2 lw 3 linecolor rgb 'black'");
            command_unit.Add("quit");
            File.WriteAllLines(command_filename, command_unit);

            Console.WriteLine("  Created command file \"" + command_filename + "\".");

            theta /= 2.0;
        }
    }

    private static void test07()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 tests the ANGLE and HEIGHT functions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double h;
        double r;
        double t;
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST07");
        Console.WriteLine("  For circle segment with a given radius R,");
        Console.WriteLine("  CIRCLE_SEGMENT_ANGLE_FROM_HEIGHT computes the angle THETA, given the height.");
        Console.WriteLine("  CIRCLE_SEGMENT_HEIGHT_FROM_ANGLE computes height H, given the angle.");
        Console.WriteLine("  Check that these functions are inverses of each other");
        Console.WriteLine("  using random values of R, T, and H.");
        Console.WriteLine("");
        Console.WriteLine("        R             H      =>     T    =>       H2");
        Console.WriteLine("");

        int seed = 123456789;

        for (test = 1; test <= 5; test++)
        {
            r = 5.0 * UniformRNG.r8_uniform_01(ref seed);
            h = 2.0 * r * UniformRNG.r8_uniform_01(ref seed);
            t = Segment.circle_segment_angle_from_height(r, h);
            double h2 = Segment.circle_segment_height_from_angle(r, t);
            Console.WriteLine("  " + r.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + h.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + t.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + h2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("        R             T      =>     H    =>       T2");
        Console.WriteLine("");
        for (test = 1; test <= 5; test++)
        {
            r = 5.0 * UniformRNG.r8_uniform_01(ref seed);
            t = 2.0 * Math.PI * UniformRNG.r8_uniform_01(ref seed);
            h = Segment.circle_segment_height_from_angle(r, t);
            double t2 = Segment.circle_segment_angle_from_height(r, h);
            Console.WriteLine("  " + r.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + t.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + h.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + t2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void test08()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 tests CIRCLE_SEGMENT_CONTAINS_POINT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] c = new double[2];
        const int n = 1000;
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST08");
        Console.WriteLine("  CIRCLE_SEGMENT_CONTAINS_POINT reports whether");
        Console.WriteLine("  a circle segment contains a point.");
        Console.WriteLine("");
        Console.WriteLine("  Pick a circle segment at random.");
        Console.WriteLine("  Compute " + n + " sample points in the surrounding box.");
        Console.WriteLine("  Compare the area of the segment to the percentage of points");
        Console.WriteLine("  contained in the circle segment.");
        Console.WriteLine("");
        Console.WriteLine("       N       Omega1          Omega2           Area         Estimate");
        Console.WriteLine("");

        const double r = 1.0;
        c[0] = 0.0;
        c[1] = 0.0;
        int seed = 123456789;
        int[] inout = new int[n];

        for (test = 1; test <= 5; test++)
        {
            double omega1 = 2.0 * Math.PI * UniformRNG.r8_uniform_01(ref seed);
            double omega2 = 2.0 * Math.PI * UniformRNG.r8_uniform_01(ref seed);

            if (omega2 < omega1)
            {
                omega2 += 2.0 * Math.PI;
            }

            double[] xy = UniformRNG.r8mat_uniform_01_new(2, n, ref seed);
            int j;
            for (j = 0; j < n; j++)
            {
                int i;
                for (i = 0; i < 2; i++)
                {
                    xy[i + j * 2] = 2.0 * xy[i + j * 2] - 1.0;
                }
            }

            for (j = 0; j < n; j++)
            {
                inout[j] = Segment.circle_segment_contains_point(r, c, omega1, omega2, xy, xyIndex: + j * 2);
            }

            double theta = Segment.circle_segment_angle_from_chord_angles(omega1, omega2);
            double area = Segment.circle_segment_area_from_angle(r, theta);
            double area_est = 4.0 * typeMethods.i4vec_sum(n, inout) / n;

            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + omega1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + omega2.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + area.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + area_est.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void test09()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_SEGMENT_TEST09 looks at the area and centroid calculations.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] c = new double[2];
        double[] p1 = new double[2];
        double[] p2 = new double[2];

        Console.WriteLine("");
        Console.WriteLine("CIRCLE_SEGMENT_TEST09");
        Console.WriteLine("  CIRCLE_SEGMENT_AREA_FROM_CHORD and");
        Console.WriteLine("  CIRCLE_SEGMENT_CENTROID_FROM_CHORD evaluate the area");
        Console.WriteLine("  and centroid of a circle segment, given R, C and P1:P2.");
        Console.WriteLine("");
        Console.WriteLine("  CIRCLE_SEGMENT_AREA_FROM_SAMPLE and");
        Console.WriteLine("  CIRCLE_SEGMENT_CENTROID_FROM_SAMPLE give us Monte Carlo estimates.");
        Console.WriteLine("");
        Console.WriteLine("  GQCIRCSEGM can estimate these values by quadrature.");
        Console.WriteLine("");
        Console.WriteLine("  Start easy, with R = 1, C = (0,0), and Theta centered.");

        int seed = 123456789;
        const double r = 1.0;
        c[0] = 0.0;
        c[1] = 0.0;
        const double theta = Math.PI / 4.0;
        double h = Segment.circle_segment_height_from_angle(r, theta);
        const double omega1 = -theta / 2.0;
        const double omega2 = +theta / 2.0;
        p1[0] = c[0] + r * Math.Cos(omega1);
        p1[1] = c[1] + r * Math.Sin(omega1);
        p2[0] = c[0] + r * Math.Cos(omega2);
        p2[1] = c[1] + r * Math.Sin(omega2);

        double a1 = Segment.circle_segment_area_from_chord(r, c, p1, p2);
        double[] d1 = Segment.circle_segment_centroid_from_chord(r, c, p1, p2);

        Console.WriteLine("");
        Console.WriteLine("         Area          CentroidX    CentroidY");
        Console.WriteLine("");
        Console.WriteLine("  " + a1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + d1[0].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + d1[1].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        //
        //  This only works because the centroid of the height-based circle segment 
        //  is easily transformed to the centroid of the chord based circle segment.
        //
        double a2 = Segment.circle_segment_area_from_height(r, h);
        double[] d2 = Segment.circle_segment_centroid_from_height(r, h);
        Console.WriteLine("  " + a2.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + d2[1].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + (-d2[0]).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

        const int n = 10000;
        double a3 = Segment.circle_segment_area_from_sample(r, c, p1, p2, n, ref seed);
        double[] d3 = Segment.circle_segment_centroid_from_sample(r, c, p1, p2, n, ref seed);
        Console.WriteLine("  " + a3.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + d3[0].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + d3[1].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void test11()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST11 demonstrates CIRCLE_SEGMENT_ROTATION_FROM_CHORD.
        //
        //  Discussion:
        //
        //    We make a table of all pairs of angles that are multiples of pi/12.
        //
        //    For each pair, we compute the rotation, that is, the angle of the
        //    central radius of the circle segment.  We print out the result in
        //    terms of multiples of pi/12.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] c = new double[2];
        int i;
        double[] p1 = new double[2];
        double[] p2 = new double[2];

        Console.WriteLine("");
        Console.WriteLine("TEST11:");
        Console.WriteLine("  CIRCLE_SEGMENT_ROTATION_FROM_CHORD is given the endpoints");
        Console.WriteLine("  of a chord, and is asked to determine the angle of the");
        Console.WriteLine("  central radius vector.");
        Console.WriteLine("");
        Console.WriteLine("  We make a table of all pairs of angles that are multiples");
        Console.WriteLine("  of pi/12, determine the corresponding chord endpoints, and");
        Console.WriteLine("  compute the rotation angle, also printed as a multiple of pi/12.");

        const double r = 2.0;
        c[0] = 3.0;
        c[1] = 5.0;
        Console.WriteLine("");
        Console.WriteLine("     0.0   1.0   2.0   3.0   4.0   5.0   6.0   7.0" + 
                          "   8.0   9.0  10.0  11.0  12.0");
        Console.WriteLine("");
        for (i = 0; i <= 12; i++)
        {
            double rho1 = i * Math.PI / 6.0;
            p1[0] = c[0] + r * Math.Cos(rho1);
            p1[1] = c[1] + r * Math.Sin(rho1);
            string cout = i.ToString(CultureInfo.InvariantCulture).PadLeft(2);
            int j;
            for (j = 0; j <= 12; j++)
            {
                double rho2 = j * Math.PI / 6.0;
                p2[0] = c[0] + r * Math.Cos(rho2);
                p2[1] = c[1] + r * Math.Sin(rho2);
                double alpha = Segment.circle_segment_rotation_from_chord(r, c, p1, p2);
                double t = 6.0 * alpha / Math.PI;
                cout += "  " + t.ToString(CultureInfo.InvariantCulture).PadLeft(4);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test13()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST13 demonstrates GAUSS for quadrature rules.
        //
        //  Discussion:
        //
        //    Some recursion coefficients ALPHA and BETA are listed in Kautsky
        //    and Elhay.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference
        //
        //    Jaroslav Kautsky, Sylvan Elhay,
        //    Calculation of the Weights of Interpolatory Quadratures,
        //    Numerische Mathematik,
        //    Volume 40, Number 3, October 1982, pages 407-422.
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST13");
        Console.WriteLine("  GAUSS computes the points and weights for a");
        Console.WriteLine("  Gauss quadrature rule, given the ALPHA and BETA");
        Console.WriteLine("  recursion coefficients.");
        //
        //  Legendre rule.
        //
        int n = 10;

        double[] alpha = new double[n];
        double[] beta = new double[n];

        for (i = 0; i < n; i++)
        {
            alpha[i] = 0.0;
            beta[i] = i switch
            {
                0 => 2.0,
                _ => 1.0 / (4.0 - 1.0 / (i * i))
            };
        }

        double[] x = new double[n];
        double[] w = new double[n];

        GaussQuadrature.gauss(n, alpha, beta, ref x, ref w);

        Console.WriteLine("");
        Console.WriteLine("  LEGENDRE RULE");
        Console.WriteLine("  Point   Weight");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        //
        //  Hermite rule.
        //
        n = 10;

        alpha = new double[n];
        beta = new double[n];

        for (i = 0; i < n; i++)
        {
            alpha[i] = 0.0;
            beta[i] = i switch
            {
                0 => Math.Sqrt(Math.PI),
                _ => i / 2.0
            };
        }

        x = new double[n];
        w = new double[n];

        GaussQuadrature.gauss(n, alpha, beta, ref x, ref w);

        Console.WriteLine("");
        Console.WriteLine("  HERMITE RULE");
        Console.WriteLine("  Point   Weight");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        //
        //  Laguerre rule.
        //
        n = 10;

        alpha = new double[n];
        beta = new double[n];

        for (i = 0; i < n; i++)
        {
            alpha[i] = 2.0 * (i + 1) - 1.0;
            beta[i] = i switch
            {
                0 => 1.0,
                _ => i * i
            };
        }

        x = new double[n];
        w = new double[n];

        GaussQuadrature.gauss(n, alpha, beta, ref x, ref w);

        Console.WriteLine("");
        Console.WriteLine("  LAGUERRE RULE");
        Console.WriteLine("  Point   Weight");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

    }

    private static void test14()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST14 demonstrates R_JACOBI.
        //
        //  Discussion:
        //
        //    R_JACOBI returns recursion coefficients ALPHA and BETA for rules
        //    using a Jacobi type weight w(x) = (1-x)^A * (1+x)^B.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference
        //
        //    Walter Gautschi,
        //    Orthogonal Polynomials: Computation and Approximation,
        //    Oxford, 2004,
        //    ISBN: 0-19-850672-4,
        //    LC: QA404.5 G3555.
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST14");
        Console.WriteLine("  R_JACOBI computes recursion coefficients ALPHA and BETA");
        Console.WriteLine("  Gauss quadrature rule, given the ALPHA and BETA");
        Console.WriteLine("  recursion coefficients.");
        //
        //  Legendre rule.
        //
        int n = 10;

        double a = 0.0;
        double b = 0.0;
        double[] alpha = new double[n];
        double[] beta = new double[n];

        Jacobi.r_jacobi(n, a, b, ref alpha, ref beta);

        Console.WriteLine("");
        Console.WriteLine("  Legendre weight");
        Console.WriteLine("  A = " + a + ",  B = " + b + "");
        Console.WriteLine("  Alpha          Beta");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + alpha[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + beta[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        //
        //  Chebyshev Type 1 rule.
        //
        n = 10;

        a = -0.5;
        b = -0.5;
        alpha = new double[n];
        beta = new double[n];

        Jacobi.r_jacobi(n, a, b, ref alpha, ref beta);

        Console.WriteLine("");
        Console.WriteLine("  Chebyshev Type 1 weight");
        Console.WriteLine("  A = " + a + ",  B = " + b + "");
        Console.WriteLine("  Alpha          Beta");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + alpha[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + beta[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        //
        //  Chebyshev Type 2 rule.
        //
        n = 10;

        a = +0.5;
        b = +0.5;
        alpha = new double[n];
        beta = new double[n];

        Jacobi.r_jacobi(n, a, b, ref alpha, ref beta);

        Console.WriteLine("");
        Console.WriteLine("  Chebyshev Type 2 weight");
        Console.WriteLine("  A = " + a + ",  B = " + b + "");
        Console.WriteLine("  Alpha          Beta");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + alpha[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + beta[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        //
        //  General Jacobi rule.
        //
        n = 10;

        a = +0.5;
        b = +1.5;
        alpha = new double[n];
        beta = new double[n];

        Jacobi.r_jacobi(n, a, b, ref alpha, ref beta);

        Console.WriteLine("");
        Console.WriteLine("  General Jacobi weight");
        Console.WriteLine("  A = " + a + ",  B = " + b + "");
        Console.WriteLine("  Alpha          Beta");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + alpha[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + beta[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }
}