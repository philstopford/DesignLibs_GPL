using System;
using System.Collections.Generic;
using System.IO;
using Burkardt;
using Burkardt.Interpolation;
using Burkardt.Types;
using InterpTest;

namespace VandermondeInterp1D
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for VANDERMONDE_INTERP_1D_TEST.
            //
            //  Discussion:
            //
            //    VANDERMONDE_INTERP_1D_TEST tests the VANDERMONDE_INTERP_1D library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 July 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int prob;
            int prob_num;

            Console.WriteLine("");
            Console.WriteLine("VANDERMONDE_INTERP_1D_TEST:");
            Console.WriteLine("  Test the VANDERMONDE_INTERP_1D library.");
            Console.WriteLine("  The QR_SOLVE library is needed.");
            Console.WriteLine("  The R8LIB library is needed.");
            Console.WriteLine("  This test needs the CONDITION library.");
            Console.WriteLine("  This test needs the TEST_INTERP library.");

            vandermonde_coef_1d_test();

            vandermonde_matrix_1d_test();

            vandermonde_value_1d_test();

            prob_num = TestInterp.p00_prob_num();

            for (prob = 1; prob <= prob_num; prob++)
            {
                test01(prob);
            }

            for (prob = 1; prob <= prob_num; prob++)
            {
                test02(prob);
            }

            Console.WriteLine("");
            Console.WriteLine("VANDERMONDE_INTERP_1D_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void vandermonde_coef_1d_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VANDERMONDE_COEF_1D_TEST tests VANDERMONDE_COEF_1D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 July 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] cd;
            int nd = 5;
            double[] xd = {0.0, 1.0, 2.0, 3.0, 4.0};
            double[] yd = {24.0, 0.0, 0.0, 0.0, 0.0};

            Console.WriteLine("");
            Console.WriteLine("VANDERMONDE_COEF_1D_TEST");
            Console.WriteLine("  VANDERMONDE_COEF_1D sets the Vandermonde coefficients for 1D interpolation.");

            typeMethods.r8vec2_print(nd, xd, yd, "  Interpolation data:");

            cd = Vandermonde.vandermonde_coef_1d(nd, xd, yd);

            typeMethods.r8vec_print(nd, cd, "  Vandermonde interpolant coefficients:");

            typeMethods.r8poly_print(nd - 1, cd, "  Vandermonde interpolant polynomial:");

        }

        static void vandermonde_matrix_1d_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VANDERMONDE_MATRIX_1D_TEST tests VANDERMONDE_MATRIX_1D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 July 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] ad;
            int nd = 4;
            double[] xd = {-1.0, 2.0, 3.0, 5.0};

            Console.WriteLine("");
            Console.WriteLine("VANDERMONDE_MATRIX_1D_TEST");
            Console.WriteLine("  VANDERMONDE_MATRIX_1D sets the Vandermonde matrix for 1D interpolation.");

            ad = Vandermonde.vandermonde_matrix_1d(nd, xd);

            typeMethods.r8mat_print(nd, nd, ad, "  Vandermonde matrix:");
        }

        static void vandermonde_value_1d_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VANDERMONDE_VALUE_1D_TEST tests VANDERMONDE_VALUE_1D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 July 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] cd = {24.0, -50.0, +35.0, -10.0, 1.0};
            int nd = 5;
            int ni = 16;
            double x_hi;
            double x_lo;
            double[] xi;
            double[] yi;

            Console.WriteLine("");
            Console.WriteLine("VANDERMONDE_VALUE_1D_TEST");
            Console.WriteLine("  VANDERMONDE_VALUE_1D evaluates a Vandermonde interpolant.");

            typeMethods.r8poly_print(nd - 1, cd, "  The polynomial:");

            x_lo = 0.0;
            x_hi = 5.0;
            xi = typeMethods.r8vec_linspace_new(ni, x_lo, x_hi);

            yi = Vandermonde.vandermonde_value_1d(nd, cd, ni, xi);

            typeMethods.r8vec2_print(ni, xi, yi, "  X, P(X):");

        }

        static void test01(int prob)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 tests VANDERMONDE_INTERP_1D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 October 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] ad;
            bool debug = false;
            double[] cd;
            double condition;
            int i;
            double int_error;
            double ld;
            double li;
            int nd;
            int ni;
            double[] xd;
            double[] xi;
            double xmax;
            double xmin;
            double[] xy;
            double[] yd;
            double[] yi;
            double ymax;
            double ymin;

            Console.WriteLine("");
            Console.WriteLine("TEST01:");
            Console.WriteLine("  Interpolate data from TEST_INTERP problem #" + prob + "");

            nd = TestInterp.p00_data_num(prob);
            Console.WriteLine("  Number of data points = " + nd + "");

            xy = TestInterp.p00_data(prob, 2, nd);

            if (debug)
            {
                typeMethods.r8mat_transpose_print(2, nd, xy, "  Data array:");
            }

            xd = new double[nd];
            yd = new double[nd];

            for (i = 0; i < nd; i++)
            {
                xd[i] = xy[0 + i * 2];
                yd[i] = xy[1 + i * 2];
            }

            //
            //  Compute Vandermonde matrix and get condition number.
            //
            ad = Vandermonde.vandermonde_matrix_1d(nd, xd);

            condition = Matrix.condition_hager(nd, ad);

            Console.WriteLine("");
            Console.WriteLine("  Condition of Vandermonde matrix is " + condition + "");
            //
            //  Solve linear system.
            //
            cd = QRSolve.qr_solve(nd, nd, ad, yd);
            //
            //  #1:  Does interpolant match function at interpolation points?
            //
            ni = nd;
            xi = typeMethods.r8vec_copy_new(ni, xd);
            yi = Vandermonde.vandermonde_value_1d(nd, cd, ni, xi);

            int_error = typeMethods.r8vec_norm_affine(ni, yi, yd) / (double) (ni);

            Console.WriteLine("");
            Console.WriteLine("  L2 interpolation error averaged per interpolant node = " + int_error + "");

            //
            //  #2: Compare estimated curve length to piecewise linear (minimal) curve length.
            //  Assume data is sorted, and normalize X and Y dimensions by (XMAX-XMIN) and
            //  (YMAX-YMIN).
            //
            xmin = typeMethods.r8vec_min(nd, xd);
            xmax = typeMethods.r8vec_max(nd, xd);
            ymin = typeMethods.r8vec_min(nd, yd);
            ymax = typeMethods.r8vec_max(nd, yd);

            ni = 501;
            xi = typeMethods.r8vec_linspace_new(ni, xmin, xmax);
            yi = Vandermonde.vandermonde_value_1d(nd, cd, ni, xi);

            ld = 0.0;
            for (i = 0; i < nd - 1; i++)
            {
                ld = ld + Math.Sqrt(Math.Pow((xd[i + 1] - xd[i]) / (xmax - xmin), 2)
                                    + Math.Pow((yd[i + 1] - yd[i]) / (ymax - ymin), 2));
            }

            li = 0.0;
            for (i = 0; i < ni - 1; i++)
            {
                li = li + Math.Sqrt(Math.Pow((xi[i + 1] - xi[i]) / (xmax - xmin), 2)
                                    + Math.Pow((yi[i + 1] - yi[i]) / (ymax - ymin), 2));
            }

            Console.WriteLine("");
            Console.WriteLine("  Normalized length of piecewise linear interpolant = " + ld + "");
            Console.WriteLine("  Normalized length of polynomial interpolant       = " + li + "");

        }

        static void test02(int prob)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 tests VANDERMONDE_INTERP_1D_MATRIX.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 June 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int PROB, the problem index.
            //
        {
            double[] ad;
            double[] cd;
            string command_filename;
            List<string> command_unit = new List<string>();
            string data_filename;
            List<string> data_unit = new List<string>();
            int i;
            string interp_filename;
            List<string> interp_unit = new List<string>();
            int j;
            int nd;
            int ni;
            string output_filename;
            string title;
            double[] xd;
            double[] xi;
            double xmax;
            double xmin;
            double[] xy;
            double[] yd;
            double[] yi;

            Console.WriteLine("");
            Console.WriteLine("TEST02:");
            Console.WriteLine("  VANDERMONDE_INTERP_1D_MATRIX sets the Vandermonde linear system");
            Console.WriteLine("  for the interpolating polynomial.");
            Console.WriteLine("  Interpolate data from TEST_INTERP problem #" + prob + "");

            nd = TestInterp.p00_data_num(prob);
            Console.WriteLine("  Number of data points = " + nd + "");

            xy = TestInterp.p00_data(prob, 2, nd);

            typeMethods.r8mat_transpose_print(2, nd, xy, "  Data array:");

            xd = new double[nd];
            yd = new double[nd];

            for (i = 0; i < nd; i++)
            {
                xd[i] = xy[0 + 2 * i];
                yd[i] = xy[1 + 2 * i];
            }

            //
            //  Compute Vandermonde matrix and get condition number.
            //
            ad = Vandermonde.vandermonde_matrix_1d(nd, xd);
            //
            //  Solve linear system.
            //
            cd = QRSolve.qr_solve(nd, nd, ad, yd);
            //
            //  Create data file.
            //
            data_filename = "data" + (prob) + ".txt";
            for (j = 0; j < nd; j++)
            {
                data_unit.Add("  " + xd[j]
                                   + "  " + yd[j] + "");
            }

            File.WriteAllLines(data_filename, data_unit);
            Console.WriteLine("");
            Console.WriteLine("  Created graphics data file \"" + data_filename + "\".");
            //
            //  Create interp file.
            //
            ni = 501;
            xmin = typeMethods.r8vec_min(nd, xd);
            xmax = typeMethods.r8vec_max(nd, xd);
            xi = typeMethods.r8vec_linspace_new(ni, xmin, xmax);
            yi = Vandermonde.vandermonde_value_1d(nd, cd, ni, xi);

            interp_filename = "interp" + (prob) + ".txt";
            for (j = 0; j < ni; j++)
            {
                interp_unit.Add("  " + xi[j]
                                     + "  " + yi[j] + "");
            }

            File.WriteAllLines(interp_filename, interp_unit);
            Console.WriteLine("  Created graphics interp file \"" + interp_filename + "\".");
            //
            //  Plot the data and the interpolant.
            //
            command_filename = "commands" + (prob) + ".txt";

            output_filename = "plot" + (prob) + ".png";

            command_unit.Add("# " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("# Usage:");
            command_unit.Add("#  gnuplot < " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("set term png");
            command_unit.Add("set output '" + output_filename + "'");
            command_unit.Add("set xlabel '<---X--->'");
            command_unit.Add("set ylabel '<---Y--->'");
            command_unit.Add("set title 'Data versus Vandermonde polynomial interpolant'");
            command_unit.Add("set grid");
            command_unit.Add("set style data lines");
            command_unit.Add("plot '" + data_filename
                                      + "' using 1:2 with points pt 7 ps 2 lc rgb 'blue',\\");
            command_unit.Add("     '" + interp_filename
                                      + "' using 1:2 lw 3 linecolor rgb 'red'");

            File.WriteAllLines(command_filename, command_unit);
            Console.WriteLine("  Created graphics command file \"" + command_filename + "\".");
        }

    }
}