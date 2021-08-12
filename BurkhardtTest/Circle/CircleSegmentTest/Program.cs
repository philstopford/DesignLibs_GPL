using System;
using System.Collections.Generic;
using System.IO;
using Burkardt;
using Burkardt.CircleNS;
using Burkardt.IO;
using Burkardt.Quadrature;
using Burkardt.Types;
using Burkardt.Uniform;

namespace CircleSegmentTest
{
    class Program
    {
        static void Main(string[] args)
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

        static void test01()

            //****************************************************************************80
            ///
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
            double area;
            double h;
            int i;
            double r;

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  CIRCLE_SEGMENT_AREA_FROM_HEIGHT computes the area of a circle segment.");

            Console.WriteLine("");
            Console.WriteLine("          R               H               Area");
            Console.WriteLine("");
            r = 1.0;
            h = 1.0;
            for (i = 0; i <= 10; i++)
            {
                area = Segment.circle_segment_area_from_height(r, h);
                Console.WriteLine("  " + r.ToString().PadLeft(14)
                    + "  " + h.ToString().PadLeft(14)
                    + "  " + area.ToString().PadLeft(14) + "");
                h = h / 2.0;
            }
        }

        static void test05()

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
            double a2;
            double h;
            double h2;
            double r;
            int seed;
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

            seed = 123456789;

            for (test = 1; test <= 5; test++)
            {
                r = 5.0 * UniformRNG.r8_uniform_01(ref seed);
                h = 2.0 * r * UniformRNG.r8_uniform_01(ref seed);
                a = Segment.circle_segment_area_from_height(r, h);
                h2 = Segment.circle_segment_height_from_area(r, a);
                Console.WriteLine("  " + r.ToString().PadLeft(12)
                    + "  " + h.ToString().PadLeft(12)
                    + "  " + a.ToString().PadLeft(12)
                    + "  " + h2.ToString().PadLeft(12) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("        R             A      =>     H    =>       A2");
            Console.WriteLine("");

            for (test = 1; test <= 5; test++)
            {
                r = 5.0 * UniformRNG.r8_uniform_01(ref seed);
                a = Math.PI * r * r * UniformRNG.r8_uniform_01(ref seed);
                h = Segment.circle_segment_height_from_area(r, a);
                a2 = Segment.circle_segment_area_from_height(r, h);
                Console.WriteLine("  " + r.ToString().PadLeft(12)
                    + "  " + a.ToString().PadLeft(12)
                    + "  " + h.ToString().PadLeft(12)
                    + "  " + a2.ToString().PadLeft(12) + "");
            }
        }

        static void test06()

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
            double[] an;
            int an_num = 51;
            string boundary_filename = "sample00_boundary.txt";
            List<string> boundary_unit = new List<string>();
            double[] boundary_x;
            double[] boundary_y;
            string command_filename = "sample00_commands.txt";
            List<string> command_unit = new List<string>();
            string data_filename = "sample00_data.txt";
            List<string> data_unit = new List<string>();
            int data_num = 100;
            double[] data_x;
            double[] data_y;
            string graphics_filename = "sample00.png";
            double h;
            int i;
            double r;
            int seed;
            int test;
            double theta;
            double thetah;

            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("CIRCLE_SEGMENT_TEST06");
            Console.WriteLine("  CIRCLE_SEGMENT_SAMPLE_FROM_HEIGHT samples a circle segment.");
            Console.WriteLine("");
            Console.WriteLine("  Plot " + data_num + " points from several segments.");
            Console.WriteLine("");

            r = 1.0;
            theta = Math.PI;

            for (test = 1; test <= 4; test++)
            {
                h = Segment.circle_segment_height_from_angle(r, theta);

                thetah = theta / 2.0;
                //
                //  Create boundary.
                //
                an = typeMethods.r8vec_linspace_new(an_num, -thetah, +thetah);
                for (i = 0; i < an_num; i++)
                {
                    an[i] = an[i] + 0.5 * Math.PI;
                }

                boundary_x = new double[an_num + 1];
                boundary_y = new double[an_num + 1];

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
                    boundary_unit.Add("  " + boundary_x[i].ToString().PadLeft(14)
                        + "  " + boundary_y[i].ToString().PadLeft(14) + "");
                }

                File.WriteAllLines(boundary_filename, boundary_unit);
                Console.WriteLine("");
                Console.WriteLine("  Created boundary file \"" + boundary_filename + "\".");
                //
                //  Create data.
                //
                data_x = new double[data_num + 1];
                data_y = new double[data_num + 1];

                Segment.circle_segment_sample_from_height(r, h, data_num, ref seed, ref data_x, ref data_y);

                Files.filename_inc(ref data_filename);
                for (i = 0; i < data_num; i++)
                {
                    data_unit.Add( "  " + data_x[i].ToString().PadLeft(14)
                        + "  " + data_y[i].ToString().PadLeft(14) + "");
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

                theta = theta / 2.0;
            }
        }

        static void test07()

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
            double h2;
            double r;
            int seed;
            double t;
            double t2;
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

            seed = 123456789;

            for (test = 1; test <= 5; test++)
            {
                r = 5.0 * UniformRNG.r8_uniform_01(ref seed);
                h = 2.0 * r * UniformRNG.r8_uniform_01(ref seed);
                t = Segment.circle_segment_angle_from_height(r, h);
                h2 = Segment.circle_segment_height_from_angle(r, t);
                Console.WriteLine("  " + r.ToString().PadLeft(12)
                    + "  " + h.ToString().PadLeft(12)
                    + "  " + t.ToString().PadLeft(12)
                    + "  " + h2.ToString().PadLeft(12) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("        R             T      =>     H    =>       T2");
            Console.WriteLine("");
            for (test = 1; test <= 5; test++)
            {
                r = 5.0 * UniformRNG.r8_uniform_01(ref seed);
                t = 2.0 * Math.PI * UniformRNG.r8_uniform_01(ref seed);
                h = Segment.circle_segment_height_from_angle(r, t);
                t2 = Segment.circle_segment_angle_from_height(r, h);
                Console.WriteLine("  " + r.ToString().PadLeft(12)
                    + "  " + t.ToString().PadLeft(12)
                    + "  " + h.ToString().PadLeft(12)
                    + "  " + t2.ToString().PadLeft(12) + "");
            }
        }

        static void test08()

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
            double area;
            double area_est;
            double[] c = new double[2];
            int i;
            int[] inout;
            int j;
            int n = 1000;
            double omega1;
            double omega2;
            double r;
            int seed;
            int test;
            double theta;
            double[] xy;

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

            r = 1.0;
            c[0] = 0.0;
            c[1] = 0.0;
            seed = 123456789;
            inout = new int[n];

            for (test = 1; test <= 5; test++)
            {
                omega1 = 2.0 * Math.PI * UniformRNG.r8_uniform_01(ref seed);
                omega2 = 2.0 * Math.PI * UniformRNG.r8_uniform_01(ref seed);

                if (omega2 < omega1)
                {
                    omega2 = omega2 + 2.0 * Math.PI;
                }

                xy = UniformRNG.r8mat_uniform_01_new(2, n, ref seed);
                for (j = 0; j < n; j++)
                {
                    for (i = 0; i < 2; i++)
                    {
                        xy[i + j * 2] = 2.0 * xy[i + j * 2] - 1.0;
                    }
                }

                for (j = 0; j < n; j++)
                {
                    inout[j] = Segment.circle_segment_contains_point(r, c, omega1, omega2, xy, xyIndex: + j * 2);
                }

                theta = Segment.circle_segment_angle_from_chord_angles(omega1, omega2);
                area = Segment.circle_segment_area_from_angle(r, theta);
                area_est = 4.0 * (double) (typeMethods.i4vec_sum(n, inout)) / (double) (n);

                Console.WriteLine("  " + n.ToString().PadLeft(6)
                    + "  " + omega1.ToString().PadLeft(14)
                    + "  " + omega2.ToString().PadLeft(14)
                    + "  " + area.ToString().PadLeft(14)
                    + "  " + area_est.ToString().PadLeft(14) + "");
            }
        }

        static void test09()

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
            double a1;
            double a2;
            double a3;
            double[] c = new double[2];
            double[] d1;
            double[] d2;
            double[] d3;
            double h;
            int n;
            double omega1;
            double omega2;
            double[] p1 = new double[2];
            double[] p2 = new double[2];
            const double pi = 3.141592653589793;
            double r;
            int seed;
            double theta;

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

            seed = 123456789;
            r = 1.0;
            c[0] = 0.0;
            c[1] = 0.0;
            theta = pi / 4.0;
            h = Segment.circle_segment_height_from_angle(r, theta);
            omega1 = -theta / 2.0;
            omega2 = +theta / 2.0;
            p1[0] = c[0] + r * Math.Cos(omega1);
            p1[1] = c[1] + r * Math.Sin(omega1);
            p2[0] = c[0] + r * Math.Cos(omega2);
            p2[1] = c[1] + r * Math.Sin(omega2);

            a1 = Segment.circle_segment_area_from_chord(r, c, p1, p2);
            d1 = Segment.circle_segment_centroid_from_chord(r, c, p1, p2);

            Console.WriteLine("");
            Console.WriteLine("         Area          CentroidX    CentroidY");
            Console.WriteLine("");
            Console.WriteLine("  " + a1.ToString().PadLeft(14)
                + "  " + d1[0].ToString().PadLeft(14)
                + "  " + d1[1].ToString().PadLeft(14) + "");
            //
            //  This only works because the centroid of the height-based circle segment 
            //  is easily transformed to the centroid of the chord based circle segment.
            //
            a2 = Segment.circle_segment_area_from_height(r, h);
            d2 = Segment.circle_segment_centroid_from_height(r, h);
            Console.WriteLine("  " + a2.ToString().PadLeft(14)
                + "  " + d2[1].ToString().PadLeft(14)
                + "  " + (-d2[0]).ToString().PadLeft(14) + "");

            n = 10000;
            a3 = Segment.circle_segment_area_from_sample(r, c, p1, p2, n, ref seed);
            d3 = Segment.circle_segment_centroid_from_sample(r, c, p1, p2, n, ref seed);
            Console.WriteLine("  " + a3.ToString().PadLeft(14)
                + "  " + d3[0].ToString().PadLeft(14)
                + "  " + d3[1].ToString().PadLeft(14) + "");

        }

        static void test11()

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
            double alpha;
            double[] c = new double[2];
            int i;
            int j;
            double[] p1 = new double[2];
            double[] p2 = new double[2];
            double r;
            double rho1;
            double rho2;
            double t;

            Console.WriteLine("");
            Console.WriteLine("TEST11:");
            Console.WriteLine("  CIRCLE_SEGMENT_ROTATION_FROM_CHORD is given the endpoints");
            Console.WriteLine("  of a chord, and is asked to determine the angle of the");
            Console.WriteLine("  central radius vector.");
            Console.WriteLine("");
            Console.WriteLine("  We make a table of all pairs of angles that are multiples");
            Console.WriteLine("  of pi/12, determine the corresponding chord endpoints, and");
            Console.WriteLine("  compute the rotation angle, also printed as a multiple of pi/12.");

            r = 2.0;
            c[0] = 3.0;
            c[1] = 5.0;
            Console.WriteLine("");
            Console.WriteLine("     0.0   1.0   2.0   3.0   4.0   5.0   6.0   7.0" + 
                              "   8.0   9.0  10.0  11.0  12.0");
            Console.WriteLine("");
            for (i = 0; i <= 12; i++)
            {
                rho1 = (double) (i) * Math.PI / 6.0;
                p1[0] = c[0] + r * Math.Cos(rho1);
                p1[1] = c[1] + r * Math.Sin(rho1);
                string cout = i.ToString().PadLeft(2);
                for (j = 0; j <= 12; j++)
                {
                    rho2 = (double) (j) * Math.PI / 6.0;
                    p2[0] = c[0] + r * Math.Cos(rho2);
                    p2[1] = c[1] + r * Math.Sin(rho2);
                    alpha = Segment.circle_segment_rotation_from_chord(r, c, p1, p2);
                    t = 6.0 * alpha / Math.PI;
                    cout += "  " + t.ToString().PadLeft(4);
                }

                Console.WriteLine(cout);
            }
        }

        static void test13()

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
            double[] alpha;
            double[] beta;
            int i;
            int n;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST13");
            Console.WriteLine("  GAUSS computes the points and weights for a");
            Console.WriteLine("  Gauss quadrature rule, given the ALPHA and BETA");
            Console.WriteLine("  recursion coefficients.");
            //
            //  Legendre rule.
            //
            n = 10;

            alpha = new double[n];
            beta = new double[n];

            for (i = 0; i < n; i++)
            {
                alpha[i] = 0.0;
                if (i == 0)
                {
                    beta[i] = 2.0;
                }
                else
                {
                    beta[i] = 1.0 / (4.0 - 1.0 / (double) (i * i));
                }
            }

            x = new double[n];
            w = new double[n];

            GaussQuadrature.gauss(n, alpha, beta, ref x, ref w);

            Console.WriteLine("");
            Console.WriteLine("  LEGENDRE RULE");
            Console.WriteLine("  Point   Weight");
            Console.WriteLine("");
            for (i = 0; i < n; i++)
            {
                Console.WriteLine("  " + x[i].ToString().PadLeft(14)
                    + "  " + w[i].ToString().PadLeft(14) + "");
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
                if (i == 0)
                {
                    beta[i] = Math.Sqrt(Math.PI);
                }
                else
                {
                    beta[i] = (double) (i) / 2.0;
                }
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
                Console.WriteLine("  " + x[i].ToString().PadLeft(14)
                    + "  " + w[i].ToString().PadLeft(14) + "");
            }

            //
            //  Laguerre rule.
            //
            n = 10;

            alpha = new double[n];
            beta = new double[n];

            for (i = 0; i < n; i++)
            {
                alpha[i] = 2.0 * (double) (i + 1) - 1.0;
                if (i == 0)
                {
                    beta[i] = 1.0;
                }
                else
                {
                    beta[i] = (double) (i * i);
                }
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
                Console.WriteLine("  " + x[i].ToString().PadLeft(14)
                    + "  " + w[i].ToString().PadLeft(14) + "");
            }

        }

        static void test14()

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
            double a;
            double[] alpha;
            double b;
            double[] beta;
            int i;
            int n;

            Console.WriteLine("");
            Console.WriteLine("TEST14");
            Console.WriteLine("  R_JACOBI computes recursion coefficients ALPHA and BETA");
            Console.WriteLine("  Gauss quadrature rule, given the ALPHA and BETA");
            Console.WriteLine("  recursion coefficients.");
            //
            //  Legendre rule.
            //
            n = 10;

            a = 0.0;
            b = 0.0;
            alpha = new double[n];
            beta = new double[n];

            Jacobi.r_jacobi(n, a, b, ref alpha, ref beta);

            Console.WriteLine("");
            Console.WriteLine("  Legendre weight");
            Console.WriteLine("  A = " + a + ",  B = " + b + "");
            Console.WriteLine("  Alpha          Beta");
            Console.WriteLine("");
            for (i = 0; i < n; i++)
            {
                Console.WriteLine("  " + alpha[i].ToString().PadLeft(14)
                    + "  " + beta[i].ToString().PadLeft(14) + "");
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
                Console.WriteLine("  " + alpha[i].ToString().PadLeft(14)
                    + "  " + beta[i].ToString().PadLeft(14) + "");
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
                Console.WriteLine("  " + alpha[i].ToString().PadLeft(14)
                    + "  " + beta[i].ToString().PadLeft(14) + "");
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
                Console.WriteLine("  " + alpha[i].ToString().PadLeft(14)
                    + "  " + beta[i].ToString().PadLeft(14) + "");
            }
        }
    }
}