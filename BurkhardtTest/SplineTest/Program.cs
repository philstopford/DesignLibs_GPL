using System;
using Burkardt.Function;
using Burkardt.MatrixNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace SplineTest;

using Polynomial = Burkardt.PolynomialNS;
using Burkardt.Spline;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SPLINE_TEST.
        //
        //  Discussion:
        //
        //    SPLINE_TEST tests SPLINE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("SPLINE_TEST");
        Console.WriteLine("  Test SPLINE.");

        parabola_val2_test();
        test002();
        test003();
        test004();
        test005();
        test006();

        test01();
        test02();
        test03();
        test04();
        test05();
        test06();
        test07();
        test08();
        test09();

        test10();
        test11();
        test115();
        test116();
        test12();
        test125();
        test126();
        test127();
        test13();
        test14();
        test145();
        test15();
        test16();
        test17();
        test18();
        test19();

        test20();
        test205();
        test21();
        test215();
        test22();
        test225();
        test23();
        test235();
        test24();
        Console.WriteLine("");
        Console.WriteLine("SPLINE_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void parabola_val2_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    parabola_val2_test tests parabola_val2().
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 December 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int NDIM = 1;
        int NDATA = 5;
        int i;
        int left;
        double[] xdata = new double[NDATA];
        double xval;
        double[] ydata = new double[NDIM * NDATA];
        double[] yval = new double[NDIM];
        double[] zdata = new double[NDATA];
        double[] zval = new double[NDIM];

        Console.WriteLine("");
        Console.WriteLine("parabola_val2_test");
        Console.WriteLine("  parabola_val2() evaluates parabolas through");
        Console.WriteLine("  3 points in a table");
        Console.WriteLine("");
        Console.WriteLine("  Our data tables will actually be parabolas:");
        Console.WriteLine("    A: 2*x^2 + 3 * x + 1.");
        Console.WriteLine("    B: 4*x^2 - 2 * x + 5.");
        Console.WriteLine("");

        for (i = 0; i < NDATA; i++)
        {
            xval = 2.0 * (i + 1);
            xdata[i] = xval;
            ydata[0 + i * NDIM] = 2.0 * xval * xval + 3.0 * xval + 1.0;
            zdata[i] = 4.0 * xval * xval - 2.0 * xval + 5.0;
            Console.WriteLine(i.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                                      + xdata[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                                                      + ydata[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                                                      + zdata[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Interpolated data:");
        Console.WriteLine("");
        Console.WriteLine("  LEFT, X, Y1, Y2");
        Console.WriteLine("");

        for (i = 1; i <= 5; i++)
        {
            xval = 2 * i - 1;
            left = i;
            if (NDATA - 2 < left)
            {
                left = NDATA - 2;
            }

            left = left switch
            {
                < 1 => 1,
                _ => left
            };

            Parabola.parabola_val2(NDIM, NDATA, xdata, ydata, left, xval, ref yval);
            Parabola.parabola_val2(NDIM, NDATA, xdata, zdata, left, xval, ref zval);

            Console.WriteLine("  "
                              + left.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + xval.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                              + yval[0].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                              + zval[0].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    private static void test002()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST002 tests R8VEC_BRACKET.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 10;
        int NTEST = 6;
        int i;
        int itest;
        int left = 0;
        int right = 0;
        double[] x = new double[N];
        double[] xtest = { -10.0, 1.0, 4.5, 5.0, 10.0, 12.0 };
        double xval;

        Console.WriteLine("");
        Console.WriteLine("TEST002");
        Console.WriteLine("  R8VEC_BRACKET finds a pair of entries in a");
        Console.WriteLine("    sorted real array which bracket a value.");

        for (i = 1; i <= N; i++)
        {
            x[i - 1] = i;
        }

        x[5] = x[4];

        typeMethods.r8vec_print(N, x, "  Sorted array:");

        for (itest = 0; itest < NTEST; itest++)
        {
            xval = xtest[itest];

            Console.WriteLine("");
            Console.WriteLine("  Search for XVAL = " + xval + "");

            typeMethods.r8vec_bracket(N, x, xval, ref left, ref right);

            Console.WriteLine("  X["
                              + left.ToString(CultureInfo.InvariantCulture).PadLeft(2) + "-1] = "
                              + x[left - 1].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");

            Console.WriteLine("  X["
                              + right.ToString(CultureInfo.InvariantCulture).PadLeft(2) + "-1] = "
                              + x[right - 1].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");

        }
    }

    private static void test003()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST003 tests R8VEC_BRACKET3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 10;
        int NTEST = 6;
        int i;
        int itest;
        int left;
        double[] x = new double[N];
        double[] xtest = { -10.0, 1.0, 4.5, 5.0, 10.0, 12.0 };
        double xval;

        Console.WriteLine("");
        Console.WriteLine("TEST003");
        Console.WriteLine("  R8VEC_BRACKET3 finds a pair of entries in a");
        Console.WriteLine("    sorted real array which bracket a value.");

        for (i = 0; i < N; i++)
        {
            x[i] = i + 1;
        }

        x[5] = x[4];

        typeMethods.r8vec_print(N, x, "  Sorted array:");

        left = (N + 1) / 2;

        for (itest = 0; itest < NTEST; itest++)
        {
            xval = xtest[itest];

            Console.WriteLine("");
            Console.WriteLine("  Search for XVAL = " + xval + "");

            Console.WriteLine("  Starting guess for interval is = " + left + "");

            typeMethods.r8vec_bracket3(N, x, xval, ref left);

            Console.WriteLine("  Nearest interval:");
            int index = left % x.Length;
            switch (index)
            {
                case < 0:
                    index += x.Length;
                    break;
            }
            Console.WriteLine("   X[" + left + "-1]= " + x[index] + "");
            Console.WriteLine("   X[" + left + 1 + "-1]= " + x[left] + "");

        }

    }

    private static void test004()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST004 tests R8VEC_ORDER_TYPE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 4;
        int NTEST = 6;
        int itest;
        int j;
        int order;
        double[] x = new double[N];

        Console.WriteLine("");
        Console.WriteLine("TEST004");
        Console.WriteLine("  R8VEC_ORDER_TYPE classifies a real vector as");
        Console.WriteLine("  -1: no order");
        Console.WriteLine("   0: all equal;");
        Console.WriteLine("   1: ascending;");
        Console.WriteLine("   2: strictly ascending;");
        Console.WriteLine("   3: descending;");
        Console.WriteLine("   4: strictly descending.");
        Console.WriteLine("");

        for (itest = 1; itest <= NTEST; itest++)
        {
            switch (itest)
            {
                case 1:
                    x[0] = 1.0;
                    x[1] = 3.0;
                    x[2] = 2.0;
                    x[3] = 4.0;
                    break;
                case 2:
                    x[0] = 2.0;
                    x[1] = 2.0;
                    x[2] = 2.0;
                    x[3] = 2.0;
                    break;
                case 3:
                    x[0] = 1.0;
                    x[1] = 2.0;
                    x[2] = 2.0;
                    x[3] = 4.0;
                    break;
                case 4:
                    x[0] = 1.0;
                    x[1] = 2.0;
                    x[2] = 3.0;
                    x[3] = 4.0;
                    break;
                case 5:
                    x[0] = 4.0;
                    x[1] = 4.0;
                    x[2] = 3.0;
                    x[3] = 1.0;
                    break;
                case 6:
                    x[0] = 9.0;
                    x[1] = 7.0;
                    x[2] = 3.0;
                    x[3] = 0.0;
                    break;
            }

            order = typeMethods.r8vec_order_type(N, x);

            Console.WriteLine("");
            Console.WriteLine("  Vector of order type " + order + ":");
            Console.WriteLine("");

            for (j = 0; j < N; j++)
            {
                Console.WriteLine("  "
                                  + j.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                  + x[j].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
            }

        }

    }

    private static void test005()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST005 tests D3_NP_FS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 10;
        double[] a;
        double[] b;
        int seed;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST005");
        Console.WriteLine("  D3_NP_FS factors and solves a tridiagonal linear system.");
        //
        //  Set the matrix.
        //
        seed = 123456789;
        a = D3.d3_uniform(N, ref seed);
        //
        //  Set the desired solution.
        //
        x = typeMethods.r8vec_indicator_new(N);
        //
        //  Compute b = A * x.
        //
        b = D3.d3_mxv(N, a, x);
        //
        //  Wipe out the solution.
        //  Solve the system.
        //

        x = D3.d3_np_fs(N, a, b);
        //
        //  Print the solution.
        //
        typeMethods.r8vec_print(N, x, "  Computed solution:");

    }

    private static void test006()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST006 tests DATA_TO_DIF and DIF_VAL.
        //
        //  Discussion:
        //
        //    This test demonstrates how divided difference approximation
        //    improves with N.
        //
        //    Evaluate these polynomials at 2.5.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int MAXTAB = 8;

        double[] diftab = new double[MAXTAB];
        double error;
        int j;
        int ntab;
        double true_value;
        double[] xtab = new double[MAXTAB];
        double xval;
        double[] ytab = new double[MAXTAB];
        double yval;

        Console.WriteLine("");
        Console.WriteLine("TEST006");
        Console.WriteLine("  Approximate Y = EXP(X) using orders 1 to " + MAXTAB + ".");

        Console.WriteLine("");
        Console.WriteLine("  Original data:");
        Console.WriteLine("");
        Console.WriteLine("       X          Y");
        Console.WriteLine("");
        for (j = 0; j < MAXTAB; j++)
        {
            xtab[j] = j;
            ytab[j] = Math.Exp(xtab[j]);

            Console.WriteLine("  "
                              + xtab[j].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                              + ytab[j].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        xval = 2.5;
        true_value = Math.Exp(xval);
        Console.WriteLine("");
        Console.WriteLine("  Evaluate at X = " + xval + " where EXP(X) = "
                          + true_value + "");
        Console.WriteLine("");
        Console.WriteLine("  Order  Approximate Y     Error");
        Console.WriteLine("");

        for (ntab = 1; ntab <= MAXTAB; ntab++)
        {

            for (j = 0; j < ntab; j++)
            {
                xtab[j] = j;
                ytab[j] = Math.Exp(xtab[j]);
            }

            Data.data_to_dif(ntab, xtab, ytab, ref diftab);

            yval = Polynomial.Dif.dif_val(ntab, xtab, diftab, xval);

            error = yval - true_value;

            Console.WriteLine("  "
                              + ntab.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + yval.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                              + error.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");

        }
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests BASIS_FUNCTION_B_VAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int NDATA = 5;

        int i;
        int j;
        int jhi;
        char mark;
        int nsample = 4;
        double[] tdata = { 0.0, 1.0, 4.0, 6.0, 10.0 };
        double thi = 0;
        double tlo = 0;
        double tval;
        double yval;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  BASIS_FUNCTION_B_VAL evaluates the ");
        Console.WriteLine("    B spline basis function.");
        Console.WriteLine("");
        Console.WriteLine("           T            B(T)");
        Console.WriteLine("");

        for (i = 0; i <= NDATA; i++)
        {
            switch (i)
            {
                case 0:
                    tlo = tdata[0] - 0.5 * (tdata[1] - tdata[0]);
                    thi = tdata[0];
                    break;
                default:
                {
                    if (i < NDATA)
                    {
                        tlo = tdata[i - 1];
                        thi = tdata[i];
                    }
                    else if (NDATA <= i)
                    {
                        tlo = tdata[NDATA - 1];
                        thi = tdata[NDATA - 1] + 0.5 * (tdata[NDATA - 1] - tdata[NDATA - 2]);
                    }

                    break;
                }
            }

            if (i < NDATA)
            {
                jhi = nsample - 1;
            }
            else
            {
                jhi = nsample;
            }

            for (j = 0; j <= jhi; j++)
            {

                tval = ((nsample - j) * tlo
                        + j * thi)
                       / nsample;

                yval = Basis.basis_function_b_val(tdata, tval);

                mark = ((0 < i) & (j == 0)) switch
                {
                    true => '*',
                    _ => ' '
                };

                Console.WriteLine("  "
                                  + mark + "  "
                                  + tval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                  + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }
        }

    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests BASIS_FUNCTION_BETA_VAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int NDATA = 5;

        double beta1;
        double beta2;
        int i;
        int j;
        int jhi;
        char mark;
        int nsample = 4;
        double[] tdata = { 0.0, 1.0, 4.0, 6.0, 10.0 };
        double thi = 0;
        double tlo = 0;
        double tval;
        double yval;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  BASIS_FUNCTION_BETA_VAL evaluates the ");
        Console.WriteLine("    Beta spline basis function.");

        beta1 = 1.0;
        beta2 = 0.0;

        Console.WriteLine("");
        Console.WriteLine("  BETA1 = " + beta1 + "");
        Console.WriteLine("  BETA2 = " + beta2 + "");
        Console.WriteLine("");
        Console.WriteLine("              T           B(T)");
        Console.WriteLine("");

        for (i = 0; i <= NDATA; i++)
        {
            switch (i)
            {
                case 0:
                    tlo = tdata[0] - 0.5 * (tdata[1] - tdata[0]);
                    thi = tdata[0];
                    break;
                default:
                {
                    if (i < NDATA)
                    {
                        tlo = tdata[i - 1];
                        thi = tdata[i];
                    }
                    else if (NDATA <= i)
                    {
                        tlo = tdata[NDATA - 1];
                        thi = tdata[NDATA - 1] + 0.5 * (tdata[NDATA - 1] - tdata[NDATA - 2]);
                    }

                    break;
                }
            }

            if (i < NDATA)
            {
                jhi = nsample - 1;
            }
            else
            {
                jhi = nsample;
            }

            for (j = 0; j <= jhi; j++)
            {
                tval = ((nsample - j) * tlo
                        + j * thi)
                       / nsample;

                yval = Basis.basis_function_beta_val(beta1, beta2, tdata, tval);

                mark = ((0 < i) & (j == 0)) switch
                {
                    true => '*',
                    _ => ' '
                };

                Console.WriteLine("  "
                                  + mark + "  "
                                  + tval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                  + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }
        }

        beta1 = 1.0;
        beta2 = 100.0;

        Console.WriteLine("");
        Console.WriteLine("  BETA1 = " + beta1 + "");
        Console.WriteLine("  BETA2 = " + beta2 + "");
        Console.WriteLine("");
        Console.WriteLine("              T           B(T)");
        Console.WriteLine("");

        for (i = 0; i <= NDATA; i++)
        {
            switch (i)
            {
                case 0:
                    tlo = tdata[0] - 0.5 * (tdata[1] - tdata[0]);
                    thi = tdata[0];
                    break;
                default:
                {
                    if (i < NDATA)
                    {
                        tlo = tdata[i - 1];
                        thi = tdata[i];
                    }
                    else if (NDATA <= i)
                    {
                        tlo = tdata[NDATA - 1];
                        thi = tdata[NDATA - 1] + 0.5 * (tdata[NDATA - 1] - tdata[NDATA - 2]);
                    }

                    break;
                }
            }

            if (i < NDATA)
            {
                jhi = nsample - 1;
            }
            else
            {
                jhi = nsample;
            }

            for (j = 0; j <= jhi; j++)
            {
                tval = ((nsample - j) * tlo
                        + j * thi)
                       / nsample;

                yval = Basis.basis_function_beta_val(beta1, beta2, tdata, tval);

                mark = ((0 < i) & (j == 0)) switch
                {
                    true => '*',
                    _ => ' '
                };

                Console.WriteLine("  "
                                  + mark + "  "
                                  + tval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                  + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }
        }

        beta1 = 100.0;
        beta2 = 0.0;

        Console.WriteLine("");
        Console.WriteLine("  BETA1 = " + beta1 + "");
        Console.WriteLine("  BETA2 = " + beta2 + "");
        Console.WriteLine("");
        Console.WriteLine("              T           B(T)");
        Console.WriteLine("");

        for (i = 0; i <= NDATA; i++)
        {
            switch (i)
            {
                case 0:
                    tlo = tdata[0] - 0.5 * (tdata[1] - tdata[0]);
                    thi = tdata[0];
                    break;
                default:
                {
                    if (i < NDATA)
                    {
                        tlo = tdata[i - 1];
                        thi = tdata[i];
                    }
                    else if (NDATA <= i)
                    {
                        tlo = tdata[NDATA - 1];
                        thi = tdata[NDATA - 1] + 0.5 * (tdata[NDATA - 1] - tdata[NDATA - 2]);
                    }

                    break;
                }
            }

            if (i < NDATA)
            {
                jhi = nsample - 1;
            }
            else
            {
                jhi = nsample;
            }

            for (j = 0; j <= jhi; j++)
            {
                tval = ((nsample - j) * tlo
                        + j * thi)
                       / nsample;

                yval = Basis.basis_function_beta_val(beta1, beta2, tdata, tval);

                mark = ((0 < i) & (j == 0)) switch
                {
                    true => '*',
                    _ => ' '
                };

                Console.WriteLine("  "
                                  + mark + "  "
                                  + tval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                  + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }
        }

    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests BASIS_MATRIX_B_UNI and BASIS_MATRIX_TMP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 4;
        int NDATA = 4;

        int i;
        int j;
        int jhi;
        int left;
        char mark;
        double[] mbasis;
        int nsample = 4;
        double[] tdata = { -1.0, 0.0, 1.0, 2.0 };
        double thi = 0;
        double tlo = 0;
        double tval;
        double[] ydata = { 4.0, 7.0, 12.0, 19.0 };
        double yval;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  BASIS_MATRIX_B_UNI sets up the basis matrix");
        Console.WriteLine("    for the uniform B spline.");

        mbasis = Basis.basis_matrix_b_uni();

        Console.WriteLine("");
        Console.WriteLine("    TDATA, YDATA");
        Console.WriteLine("");
        for (i = 0; i < NDATA; i++)
        {
            Console.WriteLine(tdata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                              + ydata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("              T      Spline(T)");
        Console.WriteLine("");

        left = 2;

        for (i = 0; i <= NDATA; i++)
        {
            switch (i)
            {
                case 0:
                    tlo = tdata[0] - 0.5 * (tdata[1] - tdata[0]);
                    thi = tdata[0];
                    break;
                default:
                {
                    if (i < NDATA)
                    {
                        tlo = tdata[i - 1];
                        thi = tdata[i];
                    }
                    else if (NDATA <= i)
                    {
                        tlo = tdata[NDATA - 1];
                        thi = tdata[NDATA - 1] + 0.5 * (tdata[NDATA - 1] - tdata[NDATA - 2]);
                    }

                    break;
                }
            }

            if (i < NDATA)
            {
                jhi = nsample - 1;
            }
            else
            {
                jhi = nsample;
            }

            for (j = 0; j <= jhi; j++)
            {
                tval = ((nsample - j) * tlo
                        + j * thi)
                       / nsample;

                yval = Basis.basis_matrix_tmp(left, N, mbasis, NDATA, tdata, ydata, tval);

                mark = ((0 < i) & (j == 0)) switch
                {
                    true => '*',
                    _ => ' '
                };

                Console.WriteLine("  "
                                  + mark + "  "
                                  + tval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                  + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");

            }
        }


    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests BASIS_MATRIX_BETA_UNI and BASIS_MATRIX_TMP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 4;
        int NDATA = 4;

        double beta1;
        double beta2;
        int i;
        int j;
        int jhi;
        int left;
        char mark;
        double[] mbasis;
        int nsample = 4;
        double[] tdata = { -1.0, 0.0, 1.0, 2.0 };
        double thi = 0;
        double tlo = 0;
        double tval;
        double[] ydata = { 4.0, 7.0, 12.0, 19.0 };
        double yval;

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  BASIS_MATRIX_BETA_UNI sets up the basis matrix");
        Console.WriteLine("    for the uniform beta spline.");
        //
        //  First test
        //
        beta1 = 1.0;
        beta2 = 0.0;

        Console.WriteLine("");
        Console.WriteLine("  BETA1 = " + beta1 + "");
        Console.WriteLine("  BETA2 = " + beta2 + "");

        mbasis = Basis.basis_matrix_beta_uni(beta1, beta2);

        Console.WriteLine("");
        Console.WriteLine("    TDATA, YDATA");
        Console.WriteLine("");
        for (i = 0; i < NDATA; i++)
        {
            Console.WriteLine(tdata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                              + ydata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        left = 2;

        Console.WriteLine("");
        Console.WriteLine("              T      Spline(T)");
        Console.WriteLine("");

        for (i = 0; i <= NDATA; i++)
        {
            switch (i)
            {
                case 0:
                    tlo = tdata[0] - 0.5 * (tdata[1] - tdata[0]);
                    thi = tdata[0];
                    break;
                default:
                {
                    if (i < NDATA)
                    {
                        tlo = tdata[i - 1];
                        thi = tdata[i];
                    }
                    else if (NDATA <= i)
                    {
                        tlo = tdata[NDATA - 1];
                        thi = tdata[NDATA - 1] + 0.5 * (tdata[NDATA - 1] - tdata[NDATA - 2]);
                    }

                    break;
                }
            }

            if (i < NDATA)
            {
                jhi = nsample - 1;
            }
            else
            {
                jhi = nsample;
            }

            for (j = 0; j <= jhi; j++)
            {
                tval = ((nsample - j) * tlo
                        + j * thi)
                       / nsample;

                yval = Basis.basis_matrix_tmp(left, N, mbasis, NDATA, tdata, ydata, tval);

                mark = ((0 < i) & (j == 0)) switch
                {
                    true => '*',
                    _ => ' '
                };

                Console.WriteLine("  "
                                  + mark + "  "
                                  + tval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                  + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }
        }

        //
        //  Second test
        //
        beta1 = 1.0;
        beta2 = 100.0;

        Console.WriteLine("");
        Console.WriteLine("  BETA1 = " + beta1 + "");
        Console.WriteLine("  BETA2 = " + beta2 + "");

        mbasis = Basis.basis_matrix_beta_uni(beta1, beta2);

        Console.WriteLine("");
        Console.WriteLine("    TDATA, YDATA");
        Console.WriteLine("");
        for (i = 0; i < NDATA; i++)
        {
            Console.WriteLine(tdata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                              + ydata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        left = 2;

        Console.WriteLine("");
        Console.WriteLine("              T      Spline(T)");
        Console.WriteLine("");

        for (i = 0; i <= NDATA; i++)
        {
            switch (i)
            {
                case 0:
                    tlo = tdata[0] - 0.5 * (tdata[1] - tdata[0]);
                    thi = tdata[0];
                    break;
                default:
                {
                    if (i < NDATA)
                    {
                        tlo = tdata[i - 1];
                        thi = tdata[i];
                    }
                    else if (NDATA <= i)
                    {
                        tlo = tdata[NDATA - 1];
                        thi = tdata[NDATA - 1] + 0.5 * (tdata[NDATA - 1] - tdata[NDATA - 2]);
                    }

                    break;
                }
            }

            if (i < NDATA)
            {
                jhi = nsample - 1;
            }
            else
            {
                jhi = nsample;
            }

            for (j = 0; j <= jhi; j++)
            {
                tval = ((nsample - j) * tlo
                        + j * thi)
                       / nsample;

                yval = Basis.basis_matrix_tmp(left, N, mbasis, NDATA, tdata, ydata, tval);

                mark = ((0 < i) & (j == 0)) switch
                {
                    true => '*',
                    _ => ' '
                };

                Console.WriteLine("  "
                                  + mark + "  "
                                  + tval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                  + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }
        }

        //
        //  Third test
        //
        beta1 = 100.0;
        beta2 = 0.0;

        Console.WriteLine("");
        Console.WriteLine("  BETA1 = " + beta1 + "");
        Console.WriteLine("  BETA2 = " + beta2 + "");

        mbasis = Basis.basis_matrix_beta_uni(beta1, beta2);

        Console.WriteLine("");
        Console.WriteLine("    TDATA, YDATA");
        Console.WriteLine("");
        for (i = 0; i < NDATA; i++)
        {
            Console.WriteLine(tdata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                              + ydata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        left = 2;

        Console.WriteLine("");
        Console.WriteLine("              T      Spline(T)");
        Console.WriteLine("");

        for (i = 0; i <= NDATA; i++)
        {
            switch (i)
            {
                case 0:
                    tlo = tdata[0] - 0.5 * (tdata[1] - tdata[0]);
                    thi = tdata[0];
                    break;
                default:
                {
                    if (i < NDATA)
                    {
                        tlo = tdata[i - 1];
                        thi = tdata[i];
                    }
                    else if (NDATA <= i)
                    {
                        tlo = tdata[NDATA - 1];
                        thi = tdata[NDATA - 1] + 0.5 * (tdata[NDATA - 1] - tdata[NDATA - 2]);
                    }

                    break;
                }
            }

            if (i < NDATA)
            {
                jhi = nsample - 1;
            }
            else
            {
                jhi = nsample;
            }

            for (j = 0; j <= jhi; j++)
            {
                tval = ((nsample - j) * tlo
                        + j * thi)
                       / nsample;

                yval = Basis.basis_matrix_tmp(left, N, mbasis, NDATA, tdata, ydata, tval);

                mark = ((0 < i) & (j == 0)) switch
                {
                    true => '*',
                    _ => ' '
                };

                Console.WriteLine("  "
                                  + mark + "  "
                                  + tval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                  + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }
        }


    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests BASIS_MATRIX_BEZIER and BASIS_MATRIX_TMP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 4;
        int NDATA = 4;

        int i;
        int j;
        int jhi;
        int left;
        char mark;
        double[] mbasis;
        int nsample = 4;
        double[] tdata = { 0.0, 0.0, 1.0, 1.0 };
        double thi = 0;
        double tlo = 0;
        double tval;
        double[] ydata = { 7.0, 8.3333333, 10.0, 12.0 };
        double yval;

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  BASIS_MATRIX_BEZIER sets up the basis");
        Console.WriteLine("    matrix for the uniform Bezier spline.");

        mbasis = Basis.basis_matrix_bezier();

        Console.WriteLine("");
        Console.WriteLine("    TDATA, YDATA");
        Console.WriteLine("");
        for (i = 0; i < NDATA; i++)
        {
            Console.WriteLine(tdata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                              + ydata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        left = 2;

        Console.WriteLine("");
        Console.WriteLine("              T      Spline(T)");
        Console.WriteLine("");


        for (i = 0; i <= NDATA; i++)
        {
            switch (i)
            {
                case 0:
                    tlo = tdata[0] - 0.5 * (tdata[1] - tdata[0]);
                    thi = tdata[0];
                    break;
                default:
                {
                    if (i < NDATA)
                    {
                        tlo = tdata[i - 1];
                        thi = tdata[i];
                    }
                    else if (NDATA <= i)
                    {
                        tlo = tdata[NDATA - 1];
                        thi = tdata[NDATA - 1] + 0.5 * (tdata[NDATA - 1] - tdata[NDATA - 2]);
                    }

                    break;
                }
            }

            if (i < NDATA)
            {
                jhi = nsample - 1;
            }
            else
            {
                jhi = nsample;
            }

            for (j = 0; j <= jhi; j++)
            {
                tval = ((nsample - j) * tlo
                        + j * thi)
                       / nsample;

                yval = Basis.basis_matrix_tmp(left, N, mbasis, NDATA, tdata, ydata, tval);

                mark = ((0 < i) & (j == 0)) switch
                {
                    true => '*',
                    _ => ' '
                };

                Console.WriteLine("  "
                                  + mark + "  "
                                  + tval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                  + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }
        }


    }

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tests BASIS_MATRIX_HERMITE and BASIS_MATRIX_TMP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 4;
        int NDATA = 4;

        int i;
        int j;
        int jhi;
        int left;
        char mark;
        double[] mbasis;
        int nsample = 4;
        double[] tdata = { 0.0, 0.0, 1.0, 1.0 };
        double thi = 0;
        double tlo = 0;
        double tval;
        double[] ydata = { 7.0, 12.0, 4.0, 6.0 };
        double yval;

        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  BASIS_MATRIX_HERMITE sets up the basis matrix");
        Console.WriteLine("    for the Hermite spline.");

        mbasis = Basis.basis_matrix_hermite();

        Console.WriteLine("");
        Console.WriteLine("    TDATA, YDATA");
        Console.WriteLine("");
        for (i = 0; i < NDATA; i++)
        {
            Console.WriteLine(tdata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                              + ydata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        left = 2;

        Console.WriteLine("");
        Console.WriteLine("              T      Spline(T)");
        Console.WriteLine("");


        for (i = 0; i <= NDATA; i++)
        {
            switch (i)
            {
                case 0:
                    tlo = tdata[0] - 0.5 * (tdata[1] - tdata[0]);
                    thi = tdata[0];
                    break;
                default:
                {
                    if (i < NDATA)
                    {
                        tlo = tdata[i - 1];
                        thi = tdata[i];
                    }
                    else if (NDATA <= i)
                    {
                        tlo = tdata[NDATA - 1];
                        thi = tdata[NDATA - 1] + 0.5 * (tdata[NDATA - 1] - tdata[NDATA - 2]);
                    }

                    break;
                }
            }

            if (i < NDATA)
            {
                jhi = nsample - 1;
            }
            else
            {
                jhi = nsample;
            }

            for (j = 0; j <= jhi; j++)
            {
                tval = ((nsample - j) * tlo
                        + j * thi)
                       / nsample;

                yval = Basis.basis_matrix_tmp(left, N, mbasis, NDATA, tdata, ydata, tval);

                mark = ((0 < i) & (j == 0)) switch
                {
                    true => '*',
                    _ => ' '
                };

                Console.WriteLine("  "
                                  + mark + "  "
                                  + tval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                  + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }
        }


    }

    private static void test07()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 tests BASIS_MATRIX_OVERHAUSER_UNI and BASIS_MATRIX_TMP.
        //
        //  Discussion:
        //
        //    YDATA(1:NDATA) = ( TDATA(1:NDATA) + 2 )**2 + 3
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 4;
        int NDATA = 4;

        int i;
        int j;
        int jhi;
        int left;
        char mark;
        double[] mbasis;
        int nsample = 4;
        double[] tdata = { -1.0, 0.0, 1.0, 2.0 };
        double thi = 0;
        double tlo = 0;
        double tval;
        double[] ydata = { 4.0, 7.0, 12.0, 19.0 };
        double yval;

        Console.WriteLine("");
        Console.WriteLine("TEST07");
        Console.WriteLine("  BASIS_MATRIX_OVERHAUSER_UNI sets up the basis");
        Console.WriteLine("    matrix for the uniform Overhauser spline.");

        mbasis = Basis.basis_matrix_overhauser_uni();

        Console.WriteLine("");
        Console.WriteLine("    TDATA, YDATA");
        Console.WriteLine("");
        for (i = 0; i < NDATA; i++)
        {
            Console.WriteLine("  "
                              + tdata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + ydata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        left = 2;

        Console.WriteLine("");
        Console.WriteLine("              T      Spline(T)");
        Console.WriteLine("");

        for (i = 0; i <= NDATA; i++)
        {
            switch (i)
            {
                case 0:
                    tlo = tdata[0] - 0.5 * (tdata[1] - tdata[0]);
                    thi = tdata[0];
                    break;
                default:
                {
                    if (i < NDATA)
                    {
                        tlo = tdata[i - 1];
                        thi = tdata[i];
                    }
                    else if (NDATA <= i)
                    {
                        tlo = tdata[NDATA - 1];
                        thi = tdata[NDATA - 1] + 0.5 * (tdata[NDATA - 1] - tdata[NDATA - 2]);
                    }

                    break;
                }
            }

            if (i < NDATA)
            {
                jhi = nsample - 1;
            }
            else
            {
                jhi = nsample;
            }

            for (j = 0; j <= jhi; j++)
            {
                tval = ((nsample - j) * tlo
                        + j * thi)
                       / nsample;

                yval = Basis.basis_matrix_tmp(left, N, mbasis, NDATA, tdata, ydata, tval);

                mark = ((0 < i) & (j == 0)) switch
                {
                    true => '*',
                    _ => ' '
                };

                Console.WriteLine("  "
                                  + mark + "  "
                                  + tval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                  + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }
        }


    }

    private static void test08()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 tests BASIS_MATRIX_OVERHAUSER_NONUNI and BASIS_MATRIX_TMP.
        //
        //  Discussion:
        //
        //    YDATA(1:NDATA) = ( TDATA(1:NDATA) - 2 )**2 + 3
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 4;
        int NDATA = 4;

        double alpha;
        double beta;
        int i;
        int j;
        int jhi;
        int left;
        char mark;
        double[] mbasis;
        int nsample = 4;
        double[] tdata = new double[NDATA];
        double thi = 0;
        double tlo = 0;
        double tval;
        double[] ydata = new double[NDATA];
        double yval;

        Console.WriteLine("");
        Console.WriteLine("TEST08");
        Console.WriteLine("  BASIS_MATRIX_OVERHAUSER_NONUNI sets up the");
        Console.WriteLine("    basis matrix for the nonuniform Overhauser");
        Console.WriteLine("    spline.");

        tdata[0] = 0.0;
        tdata[1] = 1.0;
        tdata[2] = 2.0;
        tdata[3] = 3.0;

        alpha = (tdata[2] - tdata[1]) / (tdata[2] - tdata[0]);
        beta = (tdata[2] - tdata[1]) / (tdata[3] - tdata[1]);

        Console.WriteLine("");
        Console.WriteLine("  ALPHA = " + alpha + "");
        Console.WriteLine("  BETA  = " + beta + "");

        mbasis = Basis.basis_matrix_overhauser_nonuni(alpha, beta);

        for (i = 0; i < NDATA; i++)
        {
            ydata[i] = Math.Pow(tdata[i] - 2.0, 2) + 3.0;
        }

        Console.WriteLine("");
        Console.WriteLine("    TDATA, YDATA");
        Console.WriteLine("");
        for (i = 0; i < NDATA; i++)
        {
            Console.WriteLine("  "
                              + tdata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + ydata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        left = 2;

        Console.WriteLine("");
        Console.WriteLine("              T      Spline(T)");
        Console.WriteLine("");

        for (i = 0; i <= NDATA; i++)
        {
            switch (i)
            {
                case 0:
                    tlo = tdata[0] - 0.5 * (tdata[1] - tdata[0]);
                    thi = tdata[0];
                    break;
                default:
                {
                    if (i < NDATA)
                    {
                        tlo = tdata[i - 1];
                        thi = tdata[i];
                    }
                    else if (NDATA <= i)
                    {
                        tlo = tdata[NDATA - 1];
                        thi = tdata[NDATA - 1] + 0.5 * (tdata[NDATA - 1] - tdata[NDATA - 2]);
                    }

                    break;
                }
            }

            if (i < NDATA)
            {
                jhi = nsample - 1;
            }
            else
            {
                jhi = nsample;
            }

            for (j = 0; j <= jhi; j++)
            {
                tval = ((nsample - j) * tlo
                        + j * thi)
                       / nsample;

                yval = Basis.basis_matrix_tmp(left, N, mbasis, NDATA, tdata, ydata, tval);

                mark = ((0 < i) & (j == 0)) switch
                {
                    true => '*',
                    _ => ' '
                };

                Console.WriteLine("  "
                                  + mark + "  "
                                  + tval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                  + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }
        }

        tdata[0] = 0.0;
        tdata[1] = 1.0;
        tdata[2] = 2.0;
        tdata[3] = 5.0;

        alpha = (tdata[2] - tdata[1]) / (tdata[2] - tdata[0]);
        beta = (tdata[2] - tdata[1]) / (tdata[3] - tdata[1]);

        Console.WriteLine("");
        Console.WriteLine("  ALPHA = " + alpha + "");
        Console.WriteLine("  BETA  = " + beta + "");

        mbasis = Basis.basis_matrix_overhauser_nonuni(alpha, beta);

        for (i = 0; i < NDATA; i++)
        {
            ydata[i] = Math.Pow(tdata[i] - 2.0, 2) + 3.0;
        }

        Console.WriteLine("");
        Console.WriteLine("    TDATA, YDATA");
        Console.WriteLine("");
        for (i = 0; i < NDATA; i++)
        {
            Console.WriteLine("  "
                              + tdata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + ydata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        left = 2;

        Console.WriteLine("");
        Console.WriteLine("              T      Spline(T)");
        Console.WriteLine("");

        for (i = 0; i <= NDATA; i++)
        {
            switch (i)
            {
                case 0:
                    tlo = tdata[0] - 0.5 * (tdata[1] - tdata[0]);
                    thi = tdata[0];
                    break;
                default:
                {
                    if (i < NDATA)
                    {
                        tlo = tdata[i - 1];
                        thi = tdata[i];
                    }
                    else if (NDATA <= i)
                    {
                        tlo = tdata[NDATA - 1];
                        thi = tdata[NDATA - 1] + 0.5 * (tdata[NDATA - 1] - tdata[NDATA - 2]);
                    }

                    break;
                }
            }

            if (i < NDATA)
            {
                jhi = nsample - 1;
            }
            else
            {
                jhi = nsample;
            }

            for (j = 0; j <= jhi; j++)
            {
                tval = ((nsample - j) * tlo
                        + j * thi)
                       / nsample;

                yval = Basis.basis_matrix_tmp(left, N, mbasis, NDATA, tdata, ydata, tval);

                mark = ((0 < i) & (j == 0)) switch
                {
                    true => '*',
                    _ => ' '
                };

                Console.WriteLine("  "
                                  + mark + "  "
                                  + tval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                  + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }
        }

        tdata[0] = 0.0;
        tdata[1] = 3.0;
        tdata[2] = 4.0;
        tdata[3] = 5.0;

        alpha = (tdata[2] - tdata[1]) / (tdata[2] - tdata[0]);
        beta = (tdata[2] - tdata[1]) / (tdata[3] - tdata[1]);

        Console.WriteLine("");
        Console.WriteLine("  ALPHA = " + alpha + "");
        Console.WriteLine("  BETA  = " + beta + "");

        mbasis = Basis.basis_matrix_overhauser_nonuni(alpha, beta);

        for (i = 0; i < NDATA; i++)
        {
            ydata[i] = Math.Pow(tdata[i] - 2.0, 2) + 3.0;
        }

        Console.WriteLine("");
        Console.WriteLine("    TDATA, YDATA");
        Console.WriteLine("");
        for (i = 0; i < NDATA; i++)
        {
            Console.WriteLine(tdata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                              + ydata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        left = 2;

        Console.WriteLine("");
        Console.WriteLine("              T      Spline(T)");
        Console.WriteLine("");

        for (i = 0; i <= NDATA; i++)
        {
            switch (i)
            {
                case 0:
                    tlo = tdata[0] - 0.5 * (tdata[1] - tdata[0]);
                    thi = tdata[0];
                    break;
                default:
                {
                    if (i < NDATA)
                    {
                        tlo = tdata[i - 1];
                        thi = tdata[i];
                    }
                    else if (NDATA <= i)
                    {
                        tlo = tdata[NDATA - 1];
                        thi = tdata[NDATA - 1] + 0.5 * (tdata[NDATA - 1] - tdata[NDATA - 2]);
                    }

                    break;
                }
            }

            if (i < NDATA)
            {
                jhi = nsample - 1;
            }
            else
            {
                jhi = nsample;
            }

            for (j = 0; j <= jhi; j++)
            {
                tval = ((nsample - j) * tlo
                        + j * thi)
                       / nsample;

                yval = Basis.basis_matrix_tmp(left, N, mbasis, NDATA, tdata, ydata, tval);

                mark = ((0 < i) & (j == 0)) switch
                {
                    true => '*',
                    _ => ' '
                };

                Console.WriteLine("  "
                                  + mark + "  "
                                  + tval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                  + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }
        }


    }

    private static void test09()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST09 tests BASIS_MATRIX_OVERHAUSER_NONUNI and BASIS_MATRIX_TMP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 4;
        int NDATA = 4;

        double alpha;
        double beta;
        int i;
        int j;
        int jhi;
        int left;
        char mark;
        double[] mbasis;
        int nsample = 4;
        double[] tdata = new double[NDATA];
        double thi = 0;
        double tlo = 0;
        double tval;
        double[] ydata = new double[NDATA];
        double yval;

        Console.WriteLine("");
        Console.WriteLine("TEST09");
        Console.WriteLine("  BASIS_MATRIX_OVERHAUSER_NONUNI sets up the");
        Console.WriteLine("    basis matrix for the nonuniform Overhauser ");
        Console.WriteLine("    spline.");
        Console.WriteLine("");
        Console.WriteLine("  First test that the nonuniform code can match");
        Console.WriteLine("  the uniform code.  Compare these results with");
        Console.WriteLine("  the uniform output.");
        Console.WriteLine("");

        tdata[0] = -1.0;
        tdata[1] = 0.0;
        tdata[2] = 1.0;
        tdata[3] = 2.0;

        alpha = (tdata[2] - tdata[1]) / (tdata[2] - tdata[0]);
        beta = (tdata[2] - tdata[1]) / (tdata[3] - tdata[1]);

        Console.WriteLine("");
        Console.WriteLine("  ALPHA = " + alpha + "");
        Console.WriteLine("  BETA  = " + beta + "");

        mbasis = Basis.basis_matrix_overhauser_nonuni(alpha, beta);

        for (i = 0; i < NDATA; i++)
        {
            ydata[i] = Math.Pow(tdata[i] + 2.0, 2) + 3.0;
        }

        Console.WriteLine("");
        Console.WriteLine("    TDATA, YDATA");
        Console.WriteLine("");
        for (i = 0; i < NDATA; i++)
        {
            Console.WriteLine(tdata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                              + ydata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        left = 2;

        Console.WriteLine("");
        Console.WriteLine("              T      Spline(T)");
        Console.WriteLine("");

        for (i = 0; i <= NDATA; i++)
        {
            switch (i)
            {
                case 0:
                    tlo = tdata[0] - 0.5 * (tdata[1] - tdata[0]);
                    thi = tdata[0];
                    break;
                default:
                {
                    if (i < NDATA)
                    {
                        tlo = tdata[i - 1];
                        thi = tdata[i];
                    }
                    else if (NDATA <= i)
                    {
                        tlo = tdata[NDATA - 1];
                        thi = tdata[NDATA - 1] + 0.5 * (tdata[NDATA - 1] - tdata[NDATA - 2]);
                    }

                    break;
                }
            }

            if (i < NDATA)
            {
                jhi = nsample - 1;
            }
            else
            {
                jhi = nsample;
            }

            for (j = 0; j <= jhi; j++)
            {
                tval = ((nsample - j) * tlo
                        + j * thi)
                       / nsample;

                yval = Basis.basis_matrix_tmp(left, N, mbasis, NDATA, tdata, ydata, tval);

                mark = ((0 < i) & (j == 0)) switch
                {
                    true => '*',
                    _ => ' '
                };

                Console.WriteLine("  "
                                  + mark + "  "
                                  + tval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                  + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  Now test that the nonuniform code on a");
        Console.WriteLine("  nonuniform grid.");
        Console.WriteLine("");

        tdata[0] = -4.0;
        tdata[1] = -3.0;
        tdata[2] = -1.0;
        tdata[3] = 2.0;

        alpha = (tdata[2] - tdata[1]) / (tdata[2] - tdata[0]);
        beta = (tdata[2] - tdata[1]) / (tdata[3] - tdata[1]);

        Console.WriteLine("");
        Console.WriteLine("  ALPHA = " + alpha + "");
        Console.WriteLine("  BETA  = " + beta + "");

        mbasis = Basis.basis_matrix_overhauser_nonuni(alpha, beta);

        for (i = 0; i < NDATA; i++)
        {
            ydata[i] = Math.Pow(tdata[i] + 2.0, 2) + 3.0;
        }

        Console.WriteLine("");
        Console.WriteLine("    TDATA, YDATA");
        Console.WriteLine("");
        for (i = 0; i < NDATA; i++)
        {
            Console.WriteLine(tdata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                              + ydata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        left = 2;

        Console.WriteLine("");
        Console.WriteLine("              T      Spline(T)");
        Console.WriteLine("");

        for (i = 0; i <= NDATA; i++)
        {
            switch (i)
            {
                case 0:
                    tlo = tdata[0] - 0.5 * (tdata[1] - tdata[0]);
                    thi = tdata[0];
                    break;
                default:
                {
                    if (i < NDATA)
                    {
                        tlo = tdata[i - 1];
                        thi = tdata[i];
                    }
                    else if (NDATA <= i)
                    {
                        tlo = tdata[NDATA - 1];
                        thi = tdata[NDATA - 1] + 0.5 * (tdata[NDATA - 1] - tdata[NDATA - 2]);
                    }

                    break;
                }
            }

            if (i < NDATA)
            {
                jhi = nsample - 1;
            }
            else
            {
                jhi = nsample;
            }

            for (j = 0; j <= jhi; j++)
            {
                tval = ((nsample - j) * tlo
                        + j * thi)
                       / nsample;

                yval = Basis.basis_matrix_tmp(left, N, mbasis, NDATA, tdata, ydata, tval);

                mark = ((0 < i) & (j == 0)) switch
                {
                    true => '*',
                    _ => ' '
                };

                Console.WriteLine("  "
                                  + mark + "  "
                                  + tval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                  + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }
        }


    }

    private static void test10()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST10 tests BC_VAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 2;

        int i;
        int nsample = 101;
        double t;
        double[] xcon = { 0.0, 0.75, 1.0 };
        double xval = 0;
        double[] ycon = { 1.0, 0.0, 1.0 };
        double yval = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST10");
        Console.WriteLine("  BC_VAL evaluates a general Bezier function.");
        //
        //  One point on the curve should be about (0.75, 0.536).
        //
        Console.WriteLine("");
        Console.WriteLine("        T             X(T)           Y(T)");
        Console.WriteLine("");

        for (i = 1; i <= nsample; i++)
        {
            t = (i - 1) / (double)(nsample - 1);

            Bezier.bc_val(N, t, xcon, ycon, ref xval, ref yval);

            Console.WriteLine(t.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                       + xval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                       + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  The point ( 0.75, 0.536 ) should be on the curve.");

    }

    private static void test11()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST11 tests BEZ_VAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 2;

        double a = 0.0;
        double b = 1.0;
        double bval;
        int i;
        int nsample = 21;
        double x;
        double[] y = { 1.0, 0.0, 1.0 };

        Console.WriteLine("");
        Console.WriteLine("TEST11");
        Console.WriteLine("  BEZ_VAL evaluates a Bezier function.");
        //
        //  One point on the curve should be (0.75, 20/32).
        //
        Console.WriteLine("");
        Console.WriteLine("        I             X           B");
        Console.WriteLine("");

        for (i = 1; i <= nsample; i++)
        {
            x = ((nsample - i) * a
                 + (i - 1) * b)
                / (nsample - 1);

            bval = Bezier.bez_val(N, x, a, b, y);

            Console.WriteLine(i.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                                      + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                      + bval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  When X = " + 0.75 + "");
        Console.WriteLine("  BEZ_VAL(X) should be " + 0.625 + "");

    }

    private static void test115()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST115 tests BP01.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N_MAX = 3;

        double a = 0.0;
        double b = 1.0;
        double[] bern;
        int i;
        int j;
        int n;
        int nsample = 11;
        double x;

        Console.WriteLine("");
        Console.WriteLine("TEST115");
        Console.WriteLine("  BP01 evaluates the Bernstein basis polynomials");
        Console.WriteLine("  for the interval [0,1].");

        for (n = 0; n <= N_MAX; n++)
        {
            Console.WriteLine("");
            Console.WriteLine("  Degree N = " + n + "");
            Console.WriteLine("");
            Console.WriteLine("   X         BERN(N,0,X)  BERN(N,1,X)  BERN(N,2,X)  BERN(N,3,X)");
            Console.WriteLine("");

            for (i = 1; i <= nsample; i++)
            {
                x = ((nsample - i) * a
                     + (i - 1) * b)
                    / (nsample - 1);

                bern = Polynomial.BernsteinPolynomial.bernstein_poly_01(n, x);

                string cout = x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  ";
                for (j = 0; j <= n; j++)
                {
                    cout += bern[j].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  ";
                }

                Console.WriteLine(cout);

            }
        }
    }

    private static void test116()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST116 tests BPAB.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N_MAX = 3;

        double a = 1.0;
        double b = 3.0;
        double[] bern;
        int i;
        int j;
        int n;
        int nsample = 11;
        double x;

        Console.WriteLine("");
        Console.WriteLine("TEST116");
        Console.WriteLine("  BPAB evaluates the Bernstein basis polynomials");
        Console.WriteLine("  for the interval [A,B].");

        for (n = 0; n <= N_MAX; n++)
        {
            Console.WriteLine("");
            Console.WriteLine("  Degree N = " + n + "");
            Console.WriteLine("");
            Console.WriteLine("   X         BERN(N,0,X)  BERN(N,1,X)  BERN(N,2,X)  BERN(N,3,X)");
            Console.WriteLine("");

            for (i = 1; i <= nsample; i++)
            {
                x = ((nsample - i) * a
                     + (i - 1) * b)
                    / (nsample - 1);

                bern = Polynomial.BernsteinPolynomial.bernstein_poly_ab(n, a, b, x);

                string cout = x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  ";
                for (j = 0; j <= n; j++)
                {
                    cout += bern[j].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  ";
                }

                Console.WriteLine(cout);
            }
        }
    }

    private static void test12()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST12 tests BPAB_APPROX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int MAXDATA = 10;

        double a;
        double b;
        int i;
        int ndata;
        int nsample;
        double[] xdata = new double[MAXDATA + 1];
        double xval;
        double[] ydata = new double[MAXDATA + 1];
        double yval;

        Console.WriteLine("");
        Console.WriteLine("TEST12");
        Console.WriteLine("  BPAB_APPROX evaluates the Bernstein polynomial");
        Console.WriteLine("  approximant to a function F(X).");

        a = 1.0;
        b = 3.0;

        for (ndata = 0; ndata <= 9; ndata += 3)
        {
            for (i = 0; i <= ndata; i++)
            {
                xdata[i] = ndata switch
                {
                    0 => 0.5 * (a + b),
                    _ => ((ndata - i) * a + i * b) / ndata
                };

                ydata[i] = Math.Sin(xdata[i]);
            }

            Console.WriteLine("");
            Console.WriteLine("    XDATA    YDATA");
            Console.WriteLine("");
            for (i = 0; i <= ndata; i++)
            {
                Console.WriteLine(xdata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                                  + ydata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("  Bernstein approximant of degree N = " + ndata + "");
            Console.WriteLine("");
            Console.WriteLine("    X      F(X)     BERN(X)    ERROR");
            Console.WriteLine("");

            nsample = 2 * ndata + 1;

            for (i = 1; i <= nsample; i++)
            {
                xval = nsample switch
                {
                    1 => 0.5 * (a + b),
                    _ => ((nsample - i) * a + (i - 1) * b) / (nsample - 1)
                };

                yval = Polynomial.BernsteinPolynomial.bpab_approx(ndata, a, b, ydata, xval);

                Console.WriteLine(xval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                              + Math.Sin(xval).ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                              + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                              + (yval - Math.Sin(xval)).ToString(CultureInfo.InvariantCulture).PadLeft(12) +
                                                              "");
            }
        }
    }

    private static void test125()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST125 tests LEAST_SET_OLD and LEAST_VAL_OLD.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int MAXDEG = 6;
        int NTAB = 21;

        double[] b = new double[MAXDEG];
        double[] c = new double[MAXDEG + 1];
        double[] d = new double[MAXDEG - 1];
        double eps = 0;
        double error;
        int i;
        int ierror = 0;
        int j;
        int jhi;
        int ndeg;
        double[] ptab = new double[NTAB];
        double[] xtab = new double[NTAB];
        double xval;
        double[] ytab = new double[NTAB];
        double ytrue;
        double yval;

        Console.WriteLine("");
        Console.WriteLine("TEST125");
        Console.WriteLine("  LEAST_SET_OLD sets a least squares polynomial,");
        Console.WriteLine("  LEAST_VAL_OLD evaluates it.");

        for (i = 0; i < NTAB; i++)
        {
            xtab[i] = ((NTAB - i - 1) * -1.0
                       + i * +1.0)
                      / (NTAB - 1);
            ytab[i] = (int)(Math.Exp(xtab[i]) * 100.0 + 0.5)
                      / 100.0;
        }

        Console.WriteLine("");
        Console.WriteLine("  The data to be interpolated:");
        Console.WriteLine("");
        Console.WriteLine("  Number of data values = " + NTAB + "");
        Console.WriteLine("");
        Console.WriteLine("       X             Y");
        Console.WriteLine("");

        Console.WriteLine("");
        Console.WriteLine("    XTAB    YTAB");
        Console.WriteLine("");
        for (i = 0; i < NTAB; i++)
        {
            Console.WriteLine(xtab[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                             + ytab[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        for (ndeg = 1; ndeg <= MAXDEG; ndeg++)
        {
            Console.WriteLine("");
            Console.WriteLine("  Use a polynomial of degree: " + ndeg + "");
            Console.WriteLine("");

            Polynomial.LeastSquares.least_set_old(NTAB, xtab, ytab, ndeg, ptab, b, c, d, ref eps, ref ierror);

            Console.WriteLine("");
            Console.WriteLine("  Total approximation error = " + eps + "");
            Console.WriteLine("");
            Console.WriteLine("       X            F(X)          P(X)          Error");
            Console.WriteLine("");

            for (i = 1; i <= NTAB; i++)
            {
                if (i < NTAB)
                {
                    jhi = 2;
                }
                else
                {
                    jhi = 0;
                }

                for (j = 0; j <= jhi; j++)
                {
                    if (i < NTAB)
                    {
                        xval = ((3 - j) * xtab[i - 1]
                                + j * xtab[i])
                               / 3;
                    }
                    else
                    {
                        xval = xtab[i - 1];
                    }

                    yval = Polynomial.LeastSquares.least_val_old(xval, ndeg, b, c, d);

                    ytrue = (int)(Math.Exp(xval) * 100.0 + 0.5)
                            / 100.0;

                    error = yval - ytrue;

                    Console.WriteLine(xval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                                  + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                                  + ytrue.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                                  + error.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
                }
            }
        }
    }

    private static void test126()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST126 tests LEAST_SET and LEAST_VAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int POINT_NUM = 21;
        int NTERMS = 4;

        double[] b = new double[NTERMS];
        double[] c = new double[NTERMS];
        double[] d = new double[NTERMS];
        double[] f = new double[POINT_NUM];
        int i;
        int nterms2;
        double px;
        double[] w = new double[POINT_NUM];
        double[] x = new double[POINT_NUM];

        for (i = 0; i < POINT_NUM; i++)
        {
            w[i] = 1.0;
            x[i] = -1.0 + i / 10.0;
            f[i] = x[i] * x[i] - x[i] - 6.0;
        }

        Polynomial.LeastSquares.least_set(POINT_NUM, x, f, w, NTERMS, ref b, ref c, ref d);

        Console.WriteLine("");
        Console.WriteLine("TEST126");
        Console.WriteLine("  LEAST_SET sets a least squares polynomial,");
        Console.WriteLine("  LEAST_VAL evaluates it.");
        Console.WriteLine("");
        Console.WriteLine("  X, F(X), P(X), Error");
        Console.WriteLine("");
        for (nterms2 = 1; nterms2 <= NTERMS; nterms2++)
        {
            Console.WriteLine("");
            Console.WriteLine("  Using polynomial order = " + nterms2 + "");
            Console.WriteLine("");
            for (i = 0; i < POINT_NUM; i++)
            {
                px = Polynomial.LeastSquares.least_val(nterms2, b, c, d, x[i]);
                Console.WriteLine("  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                       + "  " + f[i].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                       + "  " + px.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                       + "  " + (px - f[i]).ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }
        }
    }

    private static void test127()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST127 tests LEAST_SET and LEAST_VAL2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int POINT_NUM = 21;
        int NTERMS = 4;

        double[] b = new double[NTERMS];
        double[] c = new double[NTERMS];
        double[] d = new double[NTERMS];
        double[] f = new double[POINT_NUM];
        double[] fp = new double[POINT_NUM];
        int i;
        int nterms2;
        double px = 0;
        double pxp = 0;
        double[] w = new double[POINT_NUM];
        double[] x = new double[POINT_NUM];

        for (i = 0; i < POINT_NUM; i++)
        {
            w[i] = 1.0;
            x[i] = -1.0 + i / 10.0;
            f[i] = x[i] * x[i] - x[i] - 6.0;
            fp[i] = 2.0 * x[i] - 1.0;
        }

        Polynomial.LeastSquares.least_set(POINT_NUM, x, f, w, NTERMS, ref b, ref c, ref d);

        Console.WriteLine("");
        Console.WriteLine("TEST127");
        Console.WriteLine("  LEAST_SET sets a least squares polynomial,");
        Console.WriteLine("  LEAST_VAL2 evaluates it.");
        Console.WriteLine("");
        Console.WriteLine("  X, F(X), P(X), FP(X), PP(X)");
        Console.WriteLine("");

        for (nterms2 = 1; nterms2 <= NTERMS; nterms2++)
        {
            Console.WriteLine("");
            Console.WriteLine("  Using polynomial order = " + nterms2 + "");
            Console.WriteLine("");
            for (i = 0; i < POINT_NUM; i++)
            {
                Polynomial.LeastSquares.least_val2(nterms2, b, c, d, x[i], ref px, ref pxp);
                Console.WriteLine("  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                       + "  " + f[i].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                       + "  " + px.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                       + "  " + fp[i].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                       + "  " + pxp.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }
        }

    }

    private static void test13()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST13 tests SPLINE_B_VAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int NDATA = 11;

        int i;
        int j;
        int jhi;
        char mark;
        int nsample = 4;
        double pi = 3.141592653589793;
        double[] tdata = new double[NDATA];
        double thi = 0;
        double tlo = 0;
        double tval;
        double[] ydata = new double[NDATA];
        double yval;

        Console.WriteLine("");
        Console.WriteLine("TEST13");
        Console.WriteLine("  SPLINE_B_VAL evaluates the");
        Console.WriteLine("    B spline.");
        Console.WriteLine("");
        Console.WriteLine("  TDATA   YDATA");
        Console.WriteLine("");

        for (i = 0; i < NDATA; i++)
        {
            tdata[i] = i;
            ydata[i] = Math.Sin(2.0 * pi * tdata[i] / (NDATA - 1));
            Console.WriteLine(tdata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                              + ydata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("    T, Spline(T)");
        Console.WriteLine("");

        for (i = 0; i <= NDATA; i++)
        {
            switch (i)
            {
                case 0:
                    tlo = tdata[1] - 0.5 * (tdata[1] - tdata[0]);
                    thi = tdata[0];
                    break;
                default:
                {
                    if (i < NDATA)
                    {
                        tlo = tdata[i - 1];
                        thi = tdata[i];
                    }
                    else if (NDATA <= i)
                    {
                        tlo = tdata[NDATA - 1];
                        thi = tdata[NDATA - 1] + 0.5 * (tdata[NDATA - 1] - tdata[NDATA - 2]);
                    }

                    break;
                }
            }

            if (i < NDATA)
            {
                jhi = nsample - 1;
            }
            else
            {
                jhi = nsample;
            }

            for (j = 0; j <= jhi; j++)
            {
                tval = ((nsample - j) * tlo
                        + j * thi)
                       / nsample;

                yval = CubicB.spline_b_val(NDATA, tdata, ydata, tval);

                mark = i switch
                {
                    > 0 when j == 0 => '*',
                    _ => ' '
                };

                Console.WriteLine("  "
                                  + mark + "  "
                                  + tval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                  + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");

            }

        }

    }

    private static void test14()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST14 tests SPLINE_BETA_VAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int NDATA = 11;

        double beta1;
        double beta2;
        int i;
        int j;
        int jhi;
        char mark;
        int nsample = 4;
        double pi = 3.141592653589793;
        double[] tdata = new double[NDATA];
        double thi = 0;
        double tlo = 0;
        double tval;
        double[] ydata = new double[NDATA];
        double yval;

        Console.WriteLine("");
        Console.WriteLine("TEST14");
        Console.WriteLine("  SPLINE_BETA_VAL evaluates the BETA spline.");
        Console.WriteLine("");
        Console.WriteLine("       TDATA         YDATA");
        Console.WriteLine("");
        for (i = 0; i < NDATA; i++)
        {
            tdata[i] = i;
            ydata[i] = Math.Sin(2.0 * pi * tdata[i] / (NDATA - 1));
            Console.WriteLine(tdata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                              + ydata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        beta1 = 1.0;
        beta2 = 0.0;

        Console.WriteLine("");
        Console.WriteLine("  BETA1 = " + beta1 + "");
        Console.WriteLine("  BETA2 = " + beta2 + "");
        Console.WriteLine("");
        Console.WriteLine("    T, Spline(T)");
        Console.WriteLine("");

        for (i = 0; i <= NDATA; i++)
        {

            switch (i)
            {
                case 0:
                    tlo = tdata[0] - 0.5 * (tdata[1] - tdata[0]);
                    thi = tdata[0];
                    break;
                default:
                {
                    if (i < NDATA)
                    {
                        tlo = tdata[i - 1];
                        thi = tdata[i];
                    }
                    else if (NDATA <= i)
                    {
                        tlo = tdata[NDATA - 1];
                        thi = tdata[NDATA - 1] + 0.5 * (tdata[NDATA - 1] - tdata[NDATA - 2]);
                    }

                    break;
                }
            }

            if (i < NDATA)
            {
                jhi = nsample - 1;
            }
            else
            {
                jhi = nsample;
            }

            for (j = 0; j <= jhi; j++)
            {
                tval = ((nsample - j) * tlo
                        + j * thi)
                       / nsample;

                yval = Beta.spline_beta_val(beta1, beta2, NDATA, tdata, ydata, tval);

                mark = i switch
                {
                    > 0 when j == 0 => '*',
                    _ => ' '
                };

                Console.WriteLine("  "
                                  + mark + "  "
                                  + tval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                  + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }
        }

        beta1 = 1.0;
        beta2 = 100.0;

        Console.WriteLine("");
        Console.WriteLine("  BETA1 = " + beta1 + "");
        Console.WriteLine("  BETA2 = " + beta2 + "");
        Console.WriteLine("");
        Console.WriteLine("    T, Spline(T)");
        Console.WriteLine("");

        for (i = 0; i <= NDATA; i++)
        {
            switch (i)
            {
                case 0:
                    tlo = tdata[0] - 0.5 * (tdata[1] - tdata[0]);
                    thi = tdata[0];
                    break;
                default:
                {
                    if (i < NDATA)
                    {
                        tlo = tdata[i - 1];
                        thi = tdata[i];
                    }
                    else if (NDATA <= i)
                    {
                        tlo = tdata[NDATA - 1];
                        thi = tdata[NDATA - 1] + 0.5 * (tdata[NDATA - 1] - tdata[NDATA - 2]);
                    }

                    break;
                }
            }

            if (i < NDATA)
            {
                jhi = nsample - 1;
            }
            else
            {
                jhi = nsample;
            }

            for (j = 0; j <= jhi; j++)
            {
                tval = ((nsample - j) * tlo
                        + j * thi)
                       / nsample;

                yval = Beta.spline_beta_val(beta1, beta2, NDATA, tdata, ydata, tval);

                mark = i switch
                {
                    > 0 when j == 0 => '*',
                    _ => ' '
                };

                Console.WriteLine("  "
                                  + mark + "  "
                                  + tval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                  + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }

        }

        beta1 = 100.0;
        beta2 = 0.0;

        Console.WriteLine("");
        Console.WriteLine("  BETA1 = " + beta1 + "");
        Console.WriteLine("  BETA2 = " + beta2 + "");
        Console.WriteLine("");
        Console.WriteLine("    T, Spline(T)");
        Console.WriteLine("");

        for (i = 0; i <= NDATA; i++)
        {
            switch (i)
            {
                case 0:
                    tlo = tdata[0] - 0.5 * (tdata[1] - tdata[0]);
                    thi = tdata[0];
                    break;
                default:
                {
                    if (i < NDATA)
                    {
                        tlo = tdata[i - 1];
                        thi = tdata[i];
                    }
                    else if (NDATA <= i)
                    {
                        tlo = tdata[NDATA - 1];
                        thi = tdata[NDATA - 1] + 0.5 * (tdata[NDATA - 1] - tdata[NDATA - 2]);
                    }

                    break;
                }
            }

            if (i < NDATA)
            {
                jhi = nsample - 1;
            }
            else
            {
                jhi = nsample;
            }

            for (j = 0; j <= jhi; j++)
            {
                tval = ((nsample - j) * tlo
                        + j * thi)
                       / nsample;

                yval = Beta.spline_beta_val(beta1, beta2, NDATA, tdata, ydata, tval);

                mark = i switch
                {
                    > 0 when j == 0 => '*',
                    _ => ' '
                };

                Console.WriteLine("  "
                                  + mark + "  "
                                  + tval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                  + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }
        }

    }

    private static void test145()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST145 tests SPLINE_CONSTANT_VAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int NDATA = 12;
        int N_TEST = 20;

        double fval;
        int i;
        int j;
        int seed = 123456789;
        double[] tdata;
        double thi = 0;
        double tlo = 0;
        double[] t_test;
        double tval = 0;
        double[] ydata = new double[NDATA];
        double yval;

        Console.WriteLine("");
        Console.WriteLine("TEST145");
        Console.WriteLine("  SPLINE_CONSTANT_VAL evaluates a piecewise");
        Console.WriteLine("  constant spline.");
        Console.WriteLine("");
        Console.WriteLine("  Runge's function, evenly spaced knots.");
        //
        //  Set the data.
        //
        tlo = -1.0;
        thi = +1.0;

        tdata = typeMethods.r8vec_even_new(NDATA - 1, tlo, thi);

        for (i = 0; i < NDATA; i++)
        {
            switch (i)
            {
                case 0:
                    tval = tdata[0];
                    break;
                default:
                {
                    if (i < NDATA - 1)
                    {
                        tval = 0.5 * (tdata[i - 1] + tdata[i]);
                    }
                    else if (i == NDATA - 1)
                    {
                        tval = tdata[i - 1];
                    }

                    break;
                }
            }

            ydata[i] = frunge(tval);
        }

        Console.WriteLine("");
        Console.WriteLine("  The data to be interpolated:");
        Console.WriteLine("");
        Console.WriteLine("  Number of data values = " + NDATA + "");
        Console.WriteLine("");
        Console.WriteLine("       T             Y");
        Console.WriteLine("");

        for (i = 0; i < NDATA; i++)
        {
            Console.WriteLine("  *"
                              + "              "
                              + ydata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            if (i < NDATA - 1)
            {
                Console.WriteLine("  *"
                                  + tdata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }
        }

        //
        //  Sample the spline.
        //
        t_test = UniformRNG.r8vec_uniform_new(N_TEST, tlo - 1.0, thi + 1.0, ref seed);

        typeMethods.r8vec_sort_bubble_a(N_TEST, ref t_test);

        Console.WriteLine("");
        Console.WriteLine("     T     Y(interp)    Y(exact)");
        Console.WriteLine("");

        j = 0;
        Console.WriteLine("  *"
                          + "              "
                          + ydata[j].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        j += 1;

        for (i = 0; i < N_TEST; i++)
        {
            tval = t_test[i];

            yval = Constant.spline_constant_val(NDATA, tdata, ydata, tval);

            if (j <= NDATA - 1)
            {
                while (tdata[j - 1] <= tval)
                {
                    fval = frunge(tdata[j - 1]);
                    Console.WriteLine("  *"
                                      + tdata[j - 1].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                      + "            " + "  "
                                      + fval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
                    Console.WriteLine("  *"
                                      + "            " + "  "
                                      + ydata[j].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
                    j += 1;
                    if (NDATA <= j)
                    {
                        break;
                    }
                }
            }

            fval = frunge(tval);

            Console.WriteLine("   "
                              + tval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + fval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");

        }
    }

    private static void test15()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST15 tests SPLINE_CUBIC_SET and SPLINE_CUBIC_VAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 11;

        int i = 0;
        int ibcbeg = 0;
        int ibcend = 0;
        int j = 0;
        int jhi = 0;
        int k = 0;
        double[] t = new double[N];
        double tval = 0;
        double[] y = new double[N];
        double ybcbeg = 0;
        double ybcend = 0;
        double[] ypp;
        double yppval = 0;
        double ypval = 0;
        double yval = 0;
        //
        //  Set up the data.
        //
        Console.WriteLine("");
        Console.WriteLine("TEST15");
        Console.WriteLine("  SPLINE_CUBIC_SET sets up a cubic spline;");
        Console.WriteLine("  SPLINE_CUBIC_VAL evaluates it.");
        Console.WriteLine("");
        Console.WriteLine("  Runge's function, evenly spaced knots.");
        Console.WriteLine("");
        Console.WriteLine("     I     T         Y");
        Console.WriteLine("");

        for (i = 0; i < N; i++)
        {
            t[i] = ((N - i - 1) * -1.0
                    + i * +1.0)
                   / (N - 1);
            y[i] = frunge(t[i]);
            Console.WriteLine(i.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                                      + t[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                                                      + y[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Try various boundary conditions.
        //
        for (k = 0; k <= 4; k++)
        {
            switch (k)
            {
                case 0:
                    ibcbeg = 0;
                    ybcbeg = 0.0;

                    ibcend = 0;
                    ybcend = 0.0;

                    Console.WriteLine("");
                    Console.WriteLine("  Boundary condition 0 at both ends:");
                    Console.WriteLine("  Spline is quadratic in boundary intervals.");
                    break;
                case 1:
                    ibcbeg = 1;
                    ybcbeg = fprunge(t[0]);

                    ibcend = 1;
                    ybcend = fprunge(t[N - 1]);

                    Console.WriteLine("");
                    Console.WriteLine("  Boundary condition 1 at both ends:");
                    Console.WriteLine("  Y'(left) =  " + ybcbeg + "");
                    Console.WriteLine("  Y'(right) = " + ybcend + "");
                    break;
                case 2:
                    ibcbeg = 2;
                    ybcbeg = fpprunge(t[0]);

                    ibcend = 2;
                    ybcend = fpprunge(t[N - 1]);

                    Console.WriteLine("");
                    Console.WriteLine("  Boundary condition 2 at both ends:");
                    Console.WriteLine("  YP''(left) =  " + ybcbeg + "");
                    Console.WriteLine("  YP''(right) = " + ybcend + "");
                    break;
                case 3:
                    ibcbeg = 2;
                    ybcbeg = 0.0;

                    ibcend = 2;
                    ybcend = 0.0;

                    Console.WriteLine("");
                    Console.WriteLine("  Natural spline:");
                    Console.WriteLine("  YP''(left) =  " + ybcbeg + "");
                    Console.WriteLine("  YP''(right) = " + ybcend + "");
                    break;
                case 4:
                    ibcbeg = 3;
                    ibcend = 3;

                    Console.WriteLine("");
                    Console.WriteLine("  \"Not-a-knot\" spline:");
                    break;
            }

            ypp = Cubic.spline_cubic_set(N, t, y, ibcbeg, ybcbeg, ibcend, ybcend);

            Console.WriteLine("");
            Console.WriteLine("  SPLINE''(T), F''(T):");
            Console.WriteLine("");
            for (i = 0; i < N; i++)
            {
                Console.WriteLine(ypp[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                                                                + fpprunge(t[i]).ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("  T, SPLINE(T), F(T)");
            Console.WriteLine("");

            for (i = 0; i <= N; i++)
            {
                switch (i)
                {
                    case 0:
                        jhi = 1;
                        break;
                    default:
                    {
                        if (i < N)
                        {
                            jhi = 2;
                        }
                        else
                        {
                            jhi = 2;
                        }

                        break;
                    }
                }

                for (j = 1; j <= jhi; j++)
                {
                    switch (i)
                    {
                        case 0:
                            tval = t[0] - 1.0;
                            break;
                        default:
                        {
                            if (i < N)
                            {
                                tval = (
                                           (jhi - j + 1) * t[i - 1]
                                           + (j - 1) * t[i])
                                       / jhi;
                            }
                            else
                            {
                                tval = j switch
                                {
                                    1 => t[N - 1],
                                    _ => t[N - 1] + 1.0
                                };
                            }

                            break;
                        }
                    }

                    yval = Cubic.spline_cubic_val(N, t, y, ypp, tval, ref ypval, ref yppval);

                    Console.WriteLine(tval.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                                                                  + yval.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                                                                  + frunge(tval).ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
                }
            }
        }

    }

    private static void test16()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST16 tests SPLINE_CUBIC_SET and SPLINE_CUBIC_VAL2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 11;

        int i = 0;
        int ibcbeg = 0;
        int ibcend = 0;
        int j = 0;
        int jhi = 0;
        int k = 0;
        int left = 0;
        int left_in = 0;
        double[] t = new double[N];
        double tval = 0;
        double[] y = new double[N];
        double ybcbeg = 0;
        double ybcend = 0;
        double[] ypp;
        double yppval = 0;
        double ypval = 0;
        double yval = 0;
        //
        //  Set up the data.
        //
        Console.WriteLine("");
        Console.WriteLine("TEST16");
        Console.WriteLine("  SPLINE_CUBIC_SET sets up a cubic spline;");
        Console.WriteLine("  SPLINE_CUBIC_VAL2 evaluates it.");
        Console.WriteLine("");
        Console.WriteLine("  Runge's function, evenly spaced knots.");
        Console.WriteLine("");
        Console.WriteLine("     I      T       Y");
        Console.WriteLine("");

        for (i = 0; i < N; i++)
        {
            t[i] = ((N - i - 1) * -1.0
                    + i * +1.0)
                   / (N - 1);
            y[i] = frunge(t[i]);
            Console.WriteLine(i.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                                      + t[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                      + y[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        //
        //  Try all three types of boundary condition.
        //
        for (k = 0; k < 3; k++)
        {
            switch (k)
            {
                case 0:
                    ibcbeg = 0;
                    ybcbeg = 0.0;

                    ibcend = 0;
                    ybcend = 0.0;

                    Console.WriteLine("");
                    Console.WriteLine("  Boundary condition 0 at both ends:");
                    Console.WriteLine("  Spline is quadratic in boundary intervals.");
                    break;
                case 1:
                    ibcbeg = 1;
                    ybcbeg = fprunge(t[0]);

                    ibcend = 1;
                    ybcend = fprunge(t[N - 1]);

                    Console.WriteLine("");
                    Console.WriteLine("  Boundary condition 1 at both ends:");
                    Console.WriteLine("  Y'(left) =  " + ybcbeg + "");
                    Console.WriteLine("  Y'(right) = " + ybcend + "");
                    break;
                case 2:
                    ibcbeg = 2;
                    ybcbeg = fpprunge(t[0]);

                    ibcend = 2;
                    ybcend = fpprunge(t[N - 1]);

                    Console.WriteLine("");
                    Console.WriteLine("  Boundary condition 2 at both ends:");
                    Console.WriteLine("  YP''(left) =  " + ybcbeg + "");
                    Console.WriteLine("  YP''(right) = " + ybcend + "");
                    break;
            }

            ypp = Cubic.spline_cubic_set(N, t, y, ibcbeg, ybcbeg, ibcend, ybcend);

            Console.WriteLine("");
            Console.WriteLine("  SPLINE''(T), F''(T):");
            Console.WriteLine("");
            for (i = 0; i < N; i++)
            {
                Console.WriteLine(ypp[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                                + fpprunge(t[i]).ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }

            left = 0;

            Console.WriteLine("");
            Console.WriteLine("  T, SPLINE(T), F(T), LEFT_IN, LEFT_OUT");
            Console.WriteLine("");

            for (i = 0; i <= N; i++)
            {
                switch (i)
                {
                    case 0:
                        jhi = 1;
                        break;
                    default:
                    {
                        if (i < N)
                        {
                            jhi = 2;
                        }
                        else
                        {
                            jhi = 2;
                        }

                        break;
                    }
                }

                for (j = 1; j <= jhi; j++)
                {
                    switch (i)
                    {
                        case 0:
                            tval = t[0] - 1.0;
                            break;
                        default:
                        {
                            if (i < N)
                            {
                                tval = ((jhi - j + 1) * t[i - 1]
                                        + (j - 1) * t[i])
                                       / jhi;
                            }
                            else
                            {
                                tval = j switch
                                {
                                    1 => t[N - 1],
                                    _ => t[N - 1] + 1.0
                                };
                            }

                            break;
                        }
                    }

                    left_in = left;
                    Cubic.spline_cubic_val2(N, t, tval, ref left, y, ypp, ref yval, ref ypval,
                        ref yppval);

                    Console.WriteLine(tval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                                  + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                                  + frunge(tval).ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                                  + left_in.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                                                  + left.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
                }
            }
        }

    }

    private static void test17()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST17 tests SPLINE_CUBIC_SET and SPLINE_CUBIC_VAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 11;

        int i = 0;
        int ibcbeg = 0;
        int ibcend = 0;
        int j = 0;
        int jhi = 0;
        int k = 0;
        double[] t = new double[N];
        double tval = 0;
        double[] y = new double[N];
        double ybcbeg = 0;
        double ybcend = 0;
        double[] ypp;
        double yppval = 0;
        double ypval = 0;
        double yval = 0;
        //
        //  Set up the data.
        //
        Console.WriteLine("");
        Console.WriteLine("TEST17");
        Console.WriteLine("  SPLINE_CUBIC_SET sets up a cubic spline;");
        Console.WriteLine("  SPLINE_CUBIC_VAL evaluates it.");
        Console.WriteLine("");
        Console.WriteLine("  Cubic function, unevenly spaced knots.");
        Console.WriteLine("");
        Console.WriteLine("  T, Y");
        Console.WriteLine("");

        for (i = 0; i < N; i++)
        {
            t[i] = i / (double)(N - 1);
            t[i] *= t[i];
            y[i] = fcube(t[i]);
            Console.WriteLine(i.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                                      + t[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                                                      + y[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Try all three types of boundary condition.
        //
        for (k = 0; k < 3; k++)
        {
            switch (k)
            {
                case 0:
                    ibcbeg = 0;
                    ybcbeg = 0.0;

                    ibcend = 0;
                    ybcend = 0.0;

                    Console.WriteLine("");
                    Console.WriteLine("  Boundary condition 0 at both ends:");
                    Console.WriteLine("  Spline is quadratic in boundary intervals.");
                    break;
                case 1:
                    ibcbeg = 1;
                    ybcbeg = fpcube(t[0]);

                    ibcend = 1;
                    ybcend = fpcube(t[N - 1]);

                    Console.WriteLine("");
                    Console.WriteLine("  Boundary condition 1 at both ends:");
                    Console.WriteLine("  Y'(left) =  " + ybcbeg + "");
                    Console.WriteLine("  Y'(right) = " + ybcend + "");
                    break;
                case 2:
                    ibcbeg = 2;
                    ybcbeg = fppcube(t[0]);

                    ibcend = 2;
                    ybcend = fppcube(t[N - 1]);

                    Console.WriteLine("");
                    Console.WriteLine("  Boundary condition 2 at both ends:");
                    Console.WriteLine("  YP''(left) =  " + ybcbeg + "");
                    Console.WriteLine("  YP''(right) = " + ybcend + "");
                    break;
            }

            ypp = Cubic.spline_cubic_set(N, t, y, ibcbeg, ybcbeg, ibcend, ybcend);

            Console.WriteLine("");
            Console.WriteLine("  SPLINE''(T), F''(T):");
            Console.WriteLine("");
            for (i = 0; i < N; i++)
            {
                Console.WriteLine(ypp[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                                + fppcube(t[i]).ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("  T, SPLINE(T), F(T)");
            Console.WriteLine("");

            for (i = 0; i <= N; i++)
            {
                switch (i)
                {
                    case 0:
                        jhi = 1;
                        break;
                    default:
                    {
                        if (i < N)
                        {
                            jhi = 2;
                        }
                        else
                        {
                            jhi = 2;
                        }

                        break;
                    }
                }

                for (j = 1; j <= jhi; j++)
                {
                    switch (i)
                    {
                        case 0:
                            tval = t[0] - 1.0;
                            break;
                        default:
                        {
                            if (i < N)
                            {
                                tval = ((jhi - j + 1) * t[i - 1]
                                        + (j - 1) * t[i])
                                       / jhi;
                            }
                            else
                            {
                                tval = j switch
                                {
                                    1 => t[N - 1],
                                    _ => t[N - 1] + 1.0
                                };
                            }

                            break;
                        }
                    }

                    yval = Cubic.spline_cubic_val(N, t, y, ypp, tval, ref ypval, ref yppval);

                    Console.WriteLine(tval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                                  + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                                  + fcube(tval).ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
                }
            }
        }

    }

    private static void test18()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST18 tests SPLINE_CUBIC_SET and SPLINE_CUBIC_VAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 11;

        int i = 0;
        int ibcbeg = 0;
        int ibcend = 0;
        int j = 0;
        int jhi = 0;
        int k = 0;
        double[] t = new double[N];
        double tval = 0;
        double[] y = new double[N];
        double ybcbeg = 0;
        double ybcend = 0;
        double[] ypp;
        double yppval = 0;
        double ypval = 0;
        double yval = 0;
        //
        //  Set up the data.
        //
        Console.WriteLine("");
        Console.WriteLine("TEST18");
        Console.WriteLine("  SPLINE_CUBIC_SET sets up a cubic spline;");
        Console.WriteLine("  SPLINE_CUBIC_VAL evaluates it.");
        Console.WriteLine("");
        Console.WriteLine("  Cubic function, evenly spaced knots.");
        Console.WriteLine("");
        Console.WriteLine("        T            Y");
        Console.WriteLine("");
        for (i = 0; i < N; i++)
        {
            t[i] = i / (double)(N - 1);
            y[i] = fcube(t[i]);
            Console.WriteLine(t[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                          + y[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        //
        //  Try all three types of boundary condition.
        //
        for (k = 0; k < 3; k++)
        {
            switch (k)
            {
                case 0:
                    ibcbeg = 0;
                    ybcbeg = 0.0;

                    ibcend = 0;
                    ybcend = 0.0;

                    Console.WriteLine("");
                    Console.WriteLine("  Boundary condition 0 at both ends:");
                    Console.WriteLine("  Spline is quadratic in boundary intervals.");
                    break;
                case 1:
                    ibcbeg = 1;
                    ybcbeg = fpcube(t[0]);

                    ibcend = 1;
                    ybcend = fpcube(t[N - 1]);

                    Console.WriteLine("");
                    Console.WriteLine("  Boundary condition 1 at both ends:");
                    Console.WriteLine("  Y'(left) =  " + ybcbeg + "");
                    Console.WriteLine("  Y'(right) = " + ybcend + "");
                    break;
                case 2:
                    ibcbeg = 2;
                    ybcbeg = fppcube(t[0]);

                    ibcend = 2;
                    ybcend = fppcube(t[N - 1]);

                    Console.WriteLine("");
                    Console.WriteLine("  Boundary condition 2 at both ends:");
                    Console.WriteLine("  YP''(left) =  " + ybcbeg + "");
                    Console.WriteLine("  YP''(right) = " + ybcend + "");
                    break;
            }

            ypp = Cubic.spline_cubic_set(N, t, y, ibcbeg, ybcbeg, ibcend, ybcend);

            Console.WriteLine("");
            Console.WriteLine("     SPLINE''(T)       F''(T):");
            Console.WriteLine("");
            for (i = 0; i < N; i++)
            {
                Console.WriteLine(ypp[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                                + fppcube(t[i]).ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("        T   SPLINE(T)     F(T)");
            Console.WriteLine("");

            for (i = 0; i <= N; i++)
            {
                switch (i)
                {
                    case 0:
                        jhi = 1;
                        break;
                    default:
                    {
                        if (i < N)
                        {
                            jhi = 2;
                        }
                        else
                        {
                            jhi = 2;
                        }

                        break;
                    }
                }

                for (j = 1; j <= jhi; j++)
                {
                    switch (i)
                    {
                        case 0:
                            tval = t[0] - 1.0;
                            break;
                        default:
                        {
                            if (i < N)
                            {
                                tval = ((jhi - j + 1) * t[i - 1]
                                        + (j - 1) * t[i])
                                       / jhi;
                            }
                            else
                            {
                                tval = j switch
                                {
                                    1 => t[N - 1],
                                    _ => t[N - 1] + 1.0
                                };
                            }

                            break;
                        }
                    }

                    yval = Cubic.spline_cubic_val(N, t, y, ypp, tval, ref ypval, ref yppval);

                    Console.WriteLine("");
                    Console.WriteLine(tval.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                      + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                      + fcube(tval).ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
                    Console.WriteLine("            "
                                      + ypval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                      + fpcube(tval).ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
                    Console.WriteLine("            "
                                      + yppval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                      + fppcube(tval).ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
                }
            }
        }

    }

    private static void test19()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST19 tests SPLINE_CUBIC_SET and SPLINE_CUBIC_VAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 2;

        int i = 0;
        int ibcbeg = 0;
        int ibcend = 0;
        int j = 0;
        int jhi = 0;
        int k1 = 0;
        int k2 = 0;
        double[] t = new double[N];
        double tval = 0;
        double[] y = new double[N];
        double ybcbeg = 0;
        double ybcend = 0;
        double[] ypp;
        double yppval = 0;
        double ypval = 0;
        double yval = 0;
        //
        //  Set up the data.
        //
        Console.WriteLine("");
        Console.WriteLine("TEST19");
        Console.WriteLine("  SPLINE_CUBIC_SET sets up a cubic spline;");
        Console.WriteLine("  SPLINE_CUBIC_VAL evaluates it.");
        Console.WriteLine("");
        Console.WriteLine("  Cubic function, evenly spaced knots.");
        Console.WriteLine("  ONLY TWO KNOTS!");
        Console.WriteLine("");
        Console.WriteLine("  The data to be interpolated:");
        Console.WriteLine("");
        Console.WriteLine("  Number of data values = " + N + "");
        Console.WriteLine("");
        Console.WriteLine("           T             Y");
        Console.WriteLine("");

        for (i = 0; i < N; i++)
        {
            t[i] = i / (double)(N - 1);
            y[i] = fcube(t[i]);
            Console.WriteLine(t[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                          + y[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        //
        //  Try all nine pairs of boundary condition.
        //
        for (k1 = 0; k1 < 3; k1++)
        {
            switch (k1)
            {
                case 0:
                    ibcbeg = 0;
                    ybcbeg = 0.0;

                    Console.WriteLine("");
                    Console.WriteLine("  Boundary condition 0 at left end.");
                    break;
                case 1:
                    ibcbeg = 1;
                    ybcbeg = fpcube(t[0]);

                    Console.WriteLine("");
                    Console.WriteLine("  Boundary condition 1 at left end.");
                    Console.WriteLine("  Y'(left) =  " + ybcbeg + "");
                    break;
                case 2:
                    ibcbeg = 2;
                    ybcbeg = fppcube(t[0]);

                    Console.WriteLine("");
                    Console.WriteLine("  Boundary condition 2 at left end.");
                    Console.WriteLine("  YP''(left) =  " + ybcbeg + "");
                    break;
            }

            for (k2 = 0; k2 < 3; k2++)
            {
                switch (k2)
                {
                    case 0:
                        ibcend = 0;
                        ybcend = 0.0;

                        Console.WriteLine("  Boundary condition 0 at right end.");
                        break;
                    case 1:
                        ibcend = 1;
                        ybcend = fpcube(t[N - 1]);

                        Console.WriteLine("");
                        Console.WriteLine("  Boundary condition 1 at right end.");
                        Console.WriteLine("  Y'(right) = " + ybcend + "");
                        break;
                    case 2:
                        ibcend = 2;
                        ybcend = fppcube(t[N - 1]);

                        Console.WriteLine("");
                        Console.WriteLine("  Boundary condition 2 at right end.");
                        Console.WriteLine("  YP''(right) = " + ybcend + "");
                        break;
                }

                ypp = Cubic.spline_cubic_set(N, t, y, ibcbeg, ybcbeg, ibcend, ybcend);

                Console.WriteLine("");
                Console.WriteLine("  SPLINE''(T)        F''(T)");
                Console.WriteLine("");
                for (i = 0; i < N; i++)
                {
                    Console.WriteLine(ypp[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                                    + fppcube(t[i]).ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
                }

                Console.WriteLine("");
                Console.WriteLine("           T    SPLINE(T)         F(T)");
                Console.WriteLine("");

                for (i = 0; i <= N; i++)
                {
                    switch (i)
                    {
                        case 0:
                            jhi = 1;
                            break;
                        default:
                        {
                            if (i < N)
                            {
                                jhi = 2;
                            }
                            else
                            {
                                jhi = 2;
                            }

                            break;
                        }
                    }

                    for (j = 1; j <= jhi; j++)
                    {

                        switch (i)
                        {
                            case 0:
                                tval = t[0] - 1.0;
                                break;
                            default:
                            {
                                if (i < N)
                                {
                                    tval = ((jhi - j + 1) * t[i - 1]
                                            + (j - 1) * t[i])
                                           / jhi;
                                }
                                else
                                {
                                    tval = j switch
                                    {
                                        1 => t[N - 1],
                                        _ => t[N - 1] + 1.0
                                    };
                                }

                                break;
                            }
                        }

                        yval = Cubic.spline_cubic_val(N, t, y, ypp, tval, ref ypval, ref yppval);

                        Console.WriteLine("");
                        Console.WriteLine(tval.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                          + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                          + fcube(tval).ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
                        Console.WriteLine("            "
                                          + ypval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                          + fpcube(tval).ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
                        Console.WriteLine("            "
                                          + yppval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                          + fppcube(tval).ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
                    }
                }
            }
        }

    }

    private static void test20()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST20 tests SPLINE_HERMITE_SET and SPLINE_HERMITE_VAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int NDATA = 4;

        double[] c;
        double fpval;
        double fval;
        int i;
        int j;
        int jhi;
        char mark;
        double pi = 3.141592653589793;
        double[] tdata = new double[NDATA];
        double tval;
        double[] ydata = new double[NDATA];
        double[] ypdata = new double[NDATA];
        double ypval = 0;
        double yval = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST20");
        Console.WriteLine("  SPLINE_HERMITE_SET sets up a Hermite spline;");
        Console.WriteLine("  SPLINE_HERMITE_VAL evaluates it.");
        //
        //  Set the data.
        //
        for (i = 0; i < NDATA; i++)
        {
            tdata[i] = 0.5 * i * pi / (NDATA - 1);
            ydata[i] = Math.Sin(tdata[i]);
            ypdata[i] = Math.Cos(tdata[i]);
        }

        Console.WriteLine("");
        Console.WriteLine("  Data");
        Console.WriteLine("");
        Console.WriteLine("     TDATA(I)     YDATA[I]     Y'DATA[I]");
        Console.WriteLine("");

        for (i = 0; i < NDATA; i++)
        {
            Console.WriteLine(tdata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                              + ydata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                              + ypdata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        //
        //  Set up the spline.
        //
        c = Hermite.spline_hermite_set(NDATA, tdata, ydata, ypdata);
        //
        //  Now evaluate the spline all over the place.
        //
        Console.WriteLine("");
        Console.WriteLine("              T     Y(hermite)     Y(exact)   Y'(hermite)     Y'(exact)");
        Console.WriteLine("");

        for (i = 0; i < NDATA; i++)
        {
            if (i == NDATA - 1)
            {
                jhi = 0;
            }
            else
            {
                jhi = 2;
            }

            for (j = 0; j <= jhi; j++)
            {
                tval = 0.5 * (3 * i + j) * pi
                       / (3 * (NDATA - 1));

                fval = Math.Sin(tval);
                fpval = Math.Cos(tval);

                Hermite.spline_hermite_val(NDATA, tdata, c, tval, ref yval, ref ypval);

                mark = j switch
                {
                    0 => '*',
                    _ => ' '
                };

                Console.WriteLine("  "
                                  + mark + "  "
                                  + tval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                  + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                  + fval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                  + ypval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                  + fpval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }

        }

    }

    private static void test205()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST205 tests SPLINE_LINEAR_INT and SPLINE_LINEAR_INTSET.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 4;

        double a;
        double b;
        double[] data_x = new double[N];
        double[] data_y = new double[N];
        int i;
        double[] int_x = { 0.0, 1.0, 4.0, 5.0, 10.0 };
        double[] int_v = { 10.0, 2.0, 8.0, 27.5 };
        double value = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST205");
        Console.WriteLine("  SPLINE_LINEAR_INTSET is given some interval endpoints,");
        Console.WriteLine("  and a value associated with each interval.");
        Console.WriteLine("");
        Console.WriteLine("  It determines a linear spline, with breakpoints");
        Console.WriteLine("  at the centers of each interval, whose integral");
        Console.WriteLine("  over each interval is equal to the given value.");

        typeMethods.r8vec_print(N + 1, int_x, "  The interval end points:");
        typeMethods.r8vec_print(N, int_v, "  The desired interval integral values:");

        Linear.spline_linear_intset(N, int_x, int_v, ref data_x, ref data_y);

        typeMethods.r8vec_print(N, data_x, "  The spline break points:");
        typeMethods.r8vec_print(N, data_y, "  The spline data values: ");

        Console.WriteLine("");
        Console.WriteLine("  As a check, call SPLINE_LINEAR_INT to compute");
        Console.WriteLine("  the integral of the spline over each interval,");
        Console.WriteLine("  and compare to the desired value.");
        Console.WriteLine("");
        Console.WriteLine("       A         B       Desired      Computed");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            a = int_x[i - 1];
            b = int_x[i];
            value = Linear.spline_linear_int(N, data_x, data_y, a, b);
            Console.WriteLine(a.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                                                      + b.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                                                      + int_v[i - 1].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                      + value.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

    }

    private static void test21()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST21 tests SPLINE_LINEAR_VAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 11;

        double fval = 0;
        int i = 0;
        int j = 0;
        int jhi = 0;
        double[] t = new double[N];
        double tval = 0;
        double[] y = new double[N];
        double ypval = 0;
        double yval = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST21");
        Console.WriteLine("  SPLINE_LINEAR_VAL evaluates a linear spline.");
        Console.WriteLine("");
        Console.WriteLine("  Runge's function, evenly spaced knots.");

        for (i = 0; i < N; i++)
        {
            t[i] = ((N - i - 1) * -1.0
                    + i * +1.0)
                   / (N - 1);
            y[i] = frunge(t[i]);
        }

        Console.WriteLine("");
        Console.WriteLine("  The data to be interpolated:");
        Console.WriteLine("");
        Console.WriteLine("  Number of data values = " + N + "");
        Console.WriteLine("");
        Console.WriteLine("        T             Y");
        Console.WriteLine("");

        for (i = 0; i < N; i++)
        {
            Console.WriteLine(t[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                          + y[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Interpolation:");
        Console.WriteLine("");
        Console.WriteLine("       T             Y            Yexact");
        Console.WriteLine("");
        for (i = 0; i <= N; i++)
        {
            switch (i)
            {
                case 0:
                    jhi = 1;
                    break;
                default:
                {
                    if (i < N)
                    {
                        jhi = 2;
                    }
                    else
                    {
                        jhi = 2;
                    }

                    break;
                }
            }

            for (j = 1; j <= jhi; j++)
            {
                switch (i)
                {
                    case 0:
                        tval = t[0] - 1.0;
                        break;
                    default:
                    {
                        if (i < N)
                        {
                            tval = ((jhi - j + 1) * t[i - 1]
                                    + (j - 1) * t[i])
                                   / jhi;
                        }
                        else
                        {
                            tval = j switch
                            {
                                1 => t[N - 1],
                                _ => t[N - 1] + 1.0
                            };
                        }

                        break;
                    }
                }

                Linear.spline_linear_val(N, t, y, tval, ref yval, ref ypval);

                fval = frunge(tval);

                Console.WriteLine(tval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                              + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                              + fval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");

            }
        }

    }

    private static void test215()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST215 tests SPLINE_LINEAR_INT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 3;

        double a;
        double b;
        int i;
        double int_val;
        double[] t = { 2.0, 4.5, 7.5 };
        double[] y = { 3.0, 3.75, 5.5 };

        Console.WriteLine("");
        Console.WriteLine("TEST215");
        Console.WriteLine("  SPLINE_LINEAR_INT computes the integral of a linear spline.");
        Console.WriteLine("");

        Console.WriteLine("");
        Console.WriteLine("  The data to be interpolated:");
        Console.WriteLine("");
        Console.WriteLine("  Number of data values = " + N + "");
        Console.WriteLine("");
        Console.WriteLine("	  T		Y");
        Console.WriteLine("");

        for (i = 0; i < N; i++)
        {
            Console.WriteLine(t[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                          + y[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("    A             B           Integral");
        Console.WriteLine("");

        for (i = 1; i <= 5; i++)
        {
            switch (i)
            {
                case 1:
                    a = 0.0;
                    b = 4.0;
                    break;
                case 2:
                    a = 4.0;
                    b = 5.0;
                    break;
                case 3:
                    a = 5.0;
                    b = 10.0;
                    break;
                case 4:
                    a = 0.0;
                    b = 10.0;
                    break;
                default:
                    a = 10.0;
                    b = 0.0;
                    break;
            }

            int_val = Linear.spline_linear_int(N, t, y, a, b);

            Console.WriteLine("  "
                              + a.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + b.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + int_val.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

    }

    private static void test22()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST22 tests SPLINE_OVERHAUSER_UNI_VAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int NDATA = 11;

        int i;
        int j;
        int jhi;
        char mark;
        int nsample = 4;
        double pi = 3.141592653589793;
        double[] tdata = new double[NDATA];
        double thi = 0;
        double tlo = 0;
        double tval;
        double[] ydata = new double[NDATA];
        double yval;

        Console.WriteLine("");
        Console.WriteLine("TEST22");
        Console.WriteLine("  SPLINE_OVERHAUSER_UNI_VAL evaluates the");
        Console.WriteLine("    uniform Overhauser spline.");

        for (i = 0; i < NDATA; i++)
        {
            tdata[i] = i;
            ydata[i] = Math.Sin(2.0 * pi * tdata[i] / (NDATA - 1));
        }

        Console.WriteLine("");
        Console.WriteLine("  The data to be interpolated:");
        Console.WriteLine("");
        Console.WriteLine("  Number of data values = " + NDATA + "");
        Console.WriteLine("");
        Console.WriteLine("       T             Y");
        Console.WriteLine("");
        for (i = 0; i < NDATA; i++)
        {
            Console.WriteLine("  "
                              + tdata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + ydata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("    T, Spline(T)");
        Console.WriteLine("");

        for (i = 0; i <= NDATA; i++)
        {
            switch (i)
            {
                case 0:
                    tlo = tdata[0] - 0.5 * (tdata[1] - tdata[0]);
                    thi = tdata[0];
                    break;
                default:
                {
                    if (i < NDATA)
                    {
                        tlo = tdata[i - 1];
                        thi = tdata[i];
                    }
                    else if (NDATA <= i)
                    {
                        tlo = tdata[NDATA - 1];
                        thi = tdata[NDATA - 1] + 0.5 * (tdata[NDATA - 1] - tdata[NDATA - 2]);
                    }

                    break;
                }
            }

            if (i < NDATA)
            {
                jhi = nsample - 1;
            }
            else
            {
                jhi = nsample;
            }

            for (j = 0; j <= jhi; j++)
            {
                tval = ((nsample - j) * tlo
                        + j * thi)
                       / nsample;

                yval = Overhauser.spline_overhauser_uni_val(NDATA, tdata, ydata, tval);

                mark = ((0 < i) & (j == 0)) switch
                {
                    true => '*',
                    _ => ' '
                };

                Console.WriteLine("  "
                                  + mark + "  "
                                  + tval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                  + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");

            }

        }

    }

    private static void test225()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST225 tests SPLINE_OVERHAUSER_NONUNI_VAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int NDATA = 11;

        int i;
        int j;
        int jhi;
        char mark;
        int nsample = 4;
        double pi = 3.141592653589793;
        double[] tdata = new double[NDATA];
        double thi = 0;
        double tlo = 0;
        double tval;
        double[] ydata = new double[NDATA];
        double yval;

        Console.WriteLine("");
        Console.WriteLine("TEST225");
        Console.WriteLine("  SPLINE_OVERHAUSER_NONUNI_VAL evaluates the");
        Console.WriteLine("    nonuniform Overhauser spline.");
        Console.WriteLine("");
        Console.WriteLine("  In this draft of a test, we simply repeat the work");
        Console.WriteLine("  for the uniform test.");

        for (i = 0; i < NDATA; i++)
        {
            tdata[i] = i;
            ydata[i] = Math.Sin(2.0 * pi * tdata[i] / (NDATA - 1));
        }

        Console.WriteLine("");
        Console.WriteLine("  The data to be interpolated:");
        Console.WriteLine("");
        Console.WriteLine("  Number of data values = " + NDATA + "");
        Console.WriteLine("");
        Console.WriteLine("       T             Y");
        Console.WriteLine("");
        for (i = 0; i < NDATA; i++)
        {
            Console.WriteLine("  "
                              + tdata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + ydata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("    T, Spline(T)");
        Console.WriteLine("");

        for (i = 0; i <= NDATA; i++)
        {
            switch (i)
            {
                case 0:
                    tlo = tdata[0] - 0.5 * (tdata[1] - tdata[0]);
                    thi = tdata[0];
                    break;
                default:
                {
                    if (i < NDATA)
                    {
                        tlo = tdata[i - 1];
                        thi = tdata[i];
                    }
                    else if (NDATA <= i)
                    {
                        tlo = tdata[NDATA - 1];
                        thi = tdata[NDATA - 1] + 0.5 * (tdata[NDATA - 1] - tdata[NDATA - 2]);
                    }

                    break;
                }
            }

            if (i < NDATA)
            {
                jhi = nsample - 1;
            }
            else
            {
                jhi = nsample;
            }

            for (j = 0; j <= jhi; j++)
            {
                tval = ((nsample - j) * tlo
                        + j * thi)
                       / nsample;

                yval = Overhauser.spline_overhauser_nonuni_val(NDATA, tdata, ydata, tval);

                mark = ((0 < i) & (j == 0)) switch
                {
                    true => '*',
                    _ => ' '
                };

                Console.WriteLine("  "
                                  + mark + "  "
                                  + tval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                  + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");

            }

        }

    }

    private static void test23()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST23 tests SPLINE_OVERHAUSER_VAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int NDATA = 4;
        int NDIM = 1;

        int i;
        double[] tdata = new double[NDATA];
        double tval;
        double[] ydata = new double[NDIM * NDATA];
        double[] yval = new double[NDIM];
        double[] zdata = new double[NDATA];
        double[] zval = new double[NDIM];

        Console.WriteLine("");
        Console.WriteLine("TEST23");
        Console.WriteLine("  SPLINE_OVERHAUSER_VAL evaluates the");
        Console.WriteLine("    Overhauser spline.");
        //
        //  Set the data.
        //
        tdata[0] = 1.0;
        ydata[0 + 0 * NDIM] = 0.0;
        zdata[0 + 0 * NDIM] = 0.0;

        tdata[1] = 2.0;
        ydata[0 + 1 * NDIM] = 1.0;
        zdata[0 + 1 * NDIM] = 1.0;

        tdata[2] = 3.0;
        ydata[0 + 2 * NDIM] = 2.0;
        zdata[0 + 2 * NDIM] = -1.0;

        tdata[3] = 4.0;
        ydata[0 + 3 * NDIM] = 3.0;
        zdata[0 + 3 * NDIM] = 0.0;

        Console.WriteLine("");
        Console.WriteLine("  Data");
        Console.WriteLine("  TDATA[I], YDATA[I], ZDATA[I]");
        Console.WriteLine("");
        for (i = 0; i < NDATA; i++)
        {
            Console.WriteLine("  "
                              + tdata[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + ydata[0 + i * NDIM].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + zdata[0 + i * NDIM].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        //
        //  Now evaluate the spline all over the place.
        //
        Console.WriteLine("");
        Console.WriteLine("  T, Spline value");
        Console.WriteLine("");

        for (i = 0; i <= 6 * NDATA + 3; i++)
        {
            tval = i / 6.0;
            Overhauser.spline_overhauser_val(NDIM, NDATA, tdata, ydata, tval, yval);
            Overhauser.spline_overhauser_val(NDIM, NDATA, tdata, zdata, tval, zval);

            Console.WriteLine("  "
                              + tval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + yval[0].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + zval[0].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void test235()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST235 tests SPLINE_PCHIP_SET and SPLINE_PCHIP_VAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 21;
        int NE = 101;

        double[] d = new double[N];
        double diff;
        double[] f = new double[N];
        double[] fe = new double[NE];
        int i;
        double[] x = new double[N];
        double[] xe = new double[NE];

        Console.WriteLine("");
        Console.WriteLine("TEST235");
        Console.WriteLine("  SPLINE_PCHIP_SET sets up a piecewise cubic");
        Console.WriteLine("    Hermite interpolant.");
        Console.WriteLine("  SPLINE_PCHIP_VAL evaluates the interpolant.");
        Console.WriteLine("");
        //
        //  Compute Runge's function at N points in [-1,1].
        //
        for (i = 0; i < N; i++)
        {
            x[i] = -1.0 + i / 10.0;
            f[i] = frunge(x[i]);
        }

        //
        //  SPLINE_PCHIP_SET takes the data in X and F, and constructs a table in D
        //  that defines the interpolant.
        //
        CubicHermite.spline_pchip_set(N, x, f, ref d);
        //
        //  Evaluate the interpolant and derivative at NE points from -1 to 0.
        //
        for (i = 0; i < NE; i++)
        {
            xe[i] = -1.0 + i / (double)(NE - 1);
        }

        CubicHermite.spline_pchip_val(N, x, f, d, NE, xe, ref fe);
        //
        //  Print the table of X, F(exact) and F(interpolated)
        //
        for (i = 0; i < NE; i++)
        {
            diff = fe[i] - frunge(xe[i]);

            Console.WriteLine("  " + xe[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + frunge(xe[i]).ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + fe[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + diff.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void test24()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST24 tests SPLINE_QUADRATIC_VAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 11;

        double fval;
        int i;
        int j;
        int jhi;
        double[] t = new double[N];
        double tval;
        double[] y = new double[N];
        double ypval = 0;
        double yval = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST24");
        Console.WriteLine("  SPLINE_QUADRATIC_VAL evaluates a");
        Console.WriteLine("    quadratic spline.");
        Console.WriteLine("");
        Console.WriteLine("  Runge''s function, evenly spaced knots.");

        for (i = 0; i < N; i++)
        {
            t[i] = ((N - i - 1) * -1.0
                    + i * +1.0)
                   / (N - 1);
            y[i] = frunge(t[i]);
        }

        Console.WriteLine("");
        Console.WriteLine("  The data to be interpolated:");
        Console.WriteLine("");
        Console.WriteLine("  Number of data values = " + N + "");
        Console.WriteLine("");
        Console.WriteLine("       T             Y");
        Console.WriteLine("");
        for (i = 0; i < N; i++)
        {
            Console.WriteLine("  "
                              + t[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + y[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Interpolated values");
        Console.WriteLine("");
        Console.WriteLine("       T             Y           Y(exact)");
        Console.WriteLine("");

        for (i = 0; i <= N; i++)
        {
            switch (i)
            {
                case 0:
                    jhi = 1;
                    break;
                default:
                {
                    if (i < N)
                    {
                        jhi = 2;
                    }
                    else
                    {
                        jhi = 2;
                    }

                    break;
                }
            }

            for (j = 1; j <= jhi; j++)
            {
                switch (i)
                {
                    case 0:
                        tval = t[0] - 1.0;
                        break;
                    default:
                    {
                        if (i < N)
                        {
                            tval = ((jhi - j + 1) * t[i - 1]
                                    + (j - 1) * t[i])
                                   / jhi;
                        }
                        else
                        {
                            tval = j switch
                            {
                                1 => t[N - 1],
                                _ => t[N - 1] + 1.0
                            };
                        }

                        break;
                    }
                }

                Quadratic.spline_quadratic_val(N, t, y, tval, ref yval, ref ypval);

                fval = frunge(tval);

                Console.WriteLine(tval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                              + yval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                                              + fval.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }
        }

    }

    private static void parabola_formula(double x, ref double y, ref double yp, ref double ypp)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARABOLA_FORMULA evaluates a parabola for us.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        y = 2.0 * x * x + 3.0 * x + 1.0;
        yp = 4.0 * x + 3.0;
        ypp = 4.0;

    }

    private static double frunge(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FRUNGE sets the Runge data values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx;

        fx = 1.0 / (1.0 + 25.0 * x * x);

        return fx;
    }

    private static double fprunge(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FPRUNGE sets the Runge derivative values at the endpoints.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double bot;
        double fx;

        bot = 1.0 + 25.0 * x * x;
        fx = -50.0 * x / (bot * bot);

        return fx;
    }

    private static double fpprunge(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FPPRUNGE sets the Runge second derivative values at the endpoints.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double bot;
        double fx;

        bot = 1.0 + 25.0 * x * x;
        fx = (-50.0 + 3750.0 * x * x) / (bot * bot * bot);

        return fx;
    }

    private static double fcube(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FCUBE evaluates a cubic function.
        //
        //  Discussion:
        //
        //    Y(X) = ( ( X + 2 ) * X + 3 ) * X + 4
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 February 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the point at which the function is evaluated.
        //
        //    Output, double FCUBE, the value of the function.
        //
    {
        double fx;

        fx = ((x + 2.0) * x + 3.0) * x + 4.0;

        return fx;
    }

    private static double fpcube(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FPCUBE evaluates the derivative of a cubic function.
        //
        //  Discussion:
        //
        //    Y(X) = ( ( X + 2 ) * X + 3 ) * X + 4
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 February 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the point at which the function is evaluated.
        //
        //    Output, double FPCUBE, the value of the derivative of the function.
        //
    {
        double fx;

        fx = (3.0 * x + 4.0) * x + 3.0;

        return fx;
    }

    private static double fppcube(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FPPCUBE evaluates the second derivative of a cubic function.
        //
        //  Discussion:
        //
        //    Y(X) = ( ( X + 2 ) * X + 3 ) * X + 4
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 February 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the point at which the function is evaluated.
        //
        //    Output, double FPPCUBE, the value of the second derivative of the function.
        //
    {
        double fx;

        fx = 6.0 * x + 4.0;

        return fx;
    }
}