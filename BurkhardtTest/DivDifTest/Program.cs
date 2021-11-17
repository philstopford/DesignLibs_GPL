using System;
using Burkardt;
using Burkardt.MatrixNS;
using Burkardt.PolynomialNS;
using Burkardt.Quadrature;
using Burkardt.RootsNS;
using Burkardt.Types;

namespace DivDifTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for divDif.dif_test.
        //
        //  Discussion:
        //
        //    DIVDIF_TEST tests DIVDIF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("DIVDIF_TEST");
        Console.WriteLine("  Test DIVDIF.");

        test01();
        test02();

        dif_basis_test();

        test05();
        test06();

        r8poly_basis_test();
        r8poly_shift_test();
        ncc_rule_test();
        nco_rule_test();

        test18();

        roots_to_r8poly_test();
        dif_derivk_table_test();
        dif_basis_deriv_test();
        dif_basis_derivk_test();

        Console.WriteLine("");
        Console.WriteLine("divDif.dif_test");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests DIF*;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int MAXTAB = 10;

        double[] diftab = new double[MAXTAB];
        double[] diftab2 = new double[MAXTAB];
        double[] diftab3 = new double[MAXTAB];
        int i;
        int ntab;
        int ntab2;
        int ntab3 = 0;
        double[] xtab = new double[MAXTAB];
        double[] xtab2 = new double[MAXTAB];
        double[] xtab3 = new double[MAXTAB];
        double xval;
        double[] ytab = new double[MAXTAB];
        double yval;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  DATA_TO_DIF_DISPLAY sets up a difference table");
        Console.WriteLine("  and displays intermediate calculations;");
        Console.WriteLine("  DIF_APPEND appends a new data point;");
        Console.WriteLine("  DIF_ANTIDERIV computes the antiderivative;");
        Console.WriteLine("  DIF_DERIV_TABLE computes the derivative;");
        Console.WriteLine("  DIF_SHIFT_ZERO shifts all the abscissas to 0;");
        Console.WriteLine("  DIF_VAL evaluates at a point;");
        Console.WriteLine("");
        //
        //  Set XTAB, YTAB to X, X^2.
        //
        ntab = 4;

        for (i = 0; i < ntab; i++)
        {
            xtab[i] = i + 1;
            ytab[i] = xtab[i] * xtab[i];
        }

        Data.data_to_dif_display(ntab, xtab, ytab, ref diftab);

        Dif.dif_print(ntab, xtab, diftab, "  The divided difference polynomial:");
        //
        //  Add (5,25) to the table.
        //
        Console.WriteLine("");
        Console.WriteLine("  DIF_APPEND can add the data (5,25) to the table.");
        Console.WriteLine("");

        xval = 5.0;
        yval = 25.0;

        Dif.dif_append(ntab, xtab, diftab, xval, yval, ref ntab, ref xtab, ref diftab);

        Dif.dif_print(ntab, xtab, diftab,
            "  The updated divided difference polynomial:");
        //
        //  Evaluate the polynomial at 2.5.
        //
        Console.WriteLine("");
        Console.WriteLine("  DIF_VAL can evaluate the table at a point.");
        Console.WriteLine("");

        xval = 2.5;

        yval = Dif.dif_val(ntab, xtab, diftab, xval);

        Console.WriteLine("");
        Console.WriteLine("  DIF_VAL reports P(" + xval + ") = " + yval + "");
        //
        //  Shift the base to zero.
        //
        Dif.dif_shift_zero(ntab, ref xtab, ref diftab);

        Dif.dif_print(ntab, xtab, diftab,
            "  The divided difference table after DIF_SHIFT_ZERO:");
        //
        //  Compute a table for the derivative.
        //
        ntab2 = ntab - 1;
        Dif.dif_deriv_table(ntab, xtab, diftab, ref xtab2, diftab2);

        Dif.dif_print(ntab2, xtab2, diftab2,
            "  The divided difference table for the derivative:");

        yval = Dif.dif_val(ntab2, xtab2, diftab2, xval);

        Console.WriteLine("");
        Console.WriteLine("  DIF_VAL reports P'(" + xval + ") = " + yval + "");
        //
        //  Compute the antiderivative.
        //
        Dif.dif_antideriv(ntab, xtab, diftab, ref ntab3, xtab3, ref diftab3);

        Dif.dif_print(ntab3, xtab3, diftab3,
            "  The divided difference table for the antiderivative:");

        yval = Dif.dif_val(ntab3, xtab3, diftab3, xval);

        Console.WriteLine("");
        Console.WriteLine("  DIF_VAL reports (Anti)P(" + xval + ") = " + yval + "");
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests DATA_TO_DIF and DIF_VAL.
        //
        //  Discussion:
        //
        //    This routine demonstrates how divided difference approximation
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
        //    18 January 2007
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
        Console.WriteLine("TEST02");
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
                              + xtab[j].ToString().PadLeft(10) + "  "
                              + ytab[j].ToString().PadLeft(10) + "");
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

            yval = Dif.dif_val(ntab, xtab, diftab, xval);

            error = yval - true_value;

            Console.WriteLine("  "
                              + ntab.ToString().PadLeft(6) + "  "
                              + yval.ToString().PadLeft(10) + "  "
                              + error.ToString().PadLeft(10) + "");

        }
    }

    private static void dif_basis_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    Dif.dif_basis_test tests Dif.dif_basis().
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int NTAB = 5;

        double[] diftab;
        int i;
        int j;
        int nstep = 9;
        double[] pointer;
        int pointerIndex = 0;
        int set = 3;
        double xhi;
        double xlo;
        double[] xtab = new double[NTAB];
        double xval;
        double yval;

        Console.WriteLine("");
        Console.WriteLine("Dif.dif_basis_test");
        Console.WriteLine("  Dif.dif_basis() computes Lagrange basis polynomials");
        Console.WriteLine("  in difference form.");
        Console.WriteLine("");
        //
        //  Set the base points.
        //
        typeMethods.r8vec_indicator(NTAB, ref xtab);

        typeMethods.r8vec_print(NTAB, xtab, "  The base points:");
        //
        //  Get the difference tables for the basis polynomials and print them.
        //
        diftab = new double[NTAB * NTAB];

        Dif.dif_basis(NTAB, xtab, ref diftab);

        Console.WriteLine("");
        Console.WriteLine("  The table of difference vectors defining the basis polynomials.");
        Console.WriteLine("  Each ROW represents a polynomial.");
        Console.WriteLine("");

        pointer = diftab;
            
        string cout = "";
        for (i = 0; i < NTAB; i++)
        {
            cout = "  ";
            for (j = 0; j < NTAB; j++)
            {
                cout += pointerIndex.ToString().PadLeft(10) + "  ";
                pointerIndex++;
            }

            Console.WriteLine(cout);
        }

        //
        //  Evaluate basis polynomial 3 at a set of points.
        //
        Console.WriteLine("");
        Console.WriteLine("  Evaluate basis polynomial #" + set + " at a set of points.");
        Console.WriteLine("");
        Console.WriteLine("      X        Y");
        Console.WriteLine("");

        xhi = NTAB;
        xlo = 1.0;
        //
        //  Advance pointer to beginning of data for basis polynomial SET.
        //
        pointer = diftab;
        pointerIndex = 0;
        for (i = 1; i < set; i++)
        {
            for (j = 1; j <= NTAB; j++)
            {
                pointerIndex++;
            }
        }

        for (i = 1; i <= nstep; i++)
        {

            xval = ((nstep - i) * xlo
                    + (i - 1) * xhi)
                   / (nstep - 1);

            yval = Dif.dif_val(NTAB, xtab, pointer, xval, diftabIndex:pointerIndex);

            Console.WriteLine("  "
                              + xval.ToString().PadLeft(10) + "  "
                              + yval.ToString().PadLeft(10) + "");

        }
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests DATA_TO_DIF_DISPLAY, DIF_PRINT, DIF_SHIFT_ZERO, DIF_TO_R8POLY;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int MAXTAB = 10;

        double[] c = new double[MAXTAB];
        double[] diftab1 = new double[MAXTAB];
        double[] diftab2 = new double[MAXTAB];
        int i;
        int ntab;
        double x;
        double[] xtab1 = new double[MAXTAB];
        double[] xtab2 = new double[MAXTAB];
        double[] ytab1 = new double[MAXTAB];
        double[] ytab2 = new double[MAXTAB];

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  DIF_TO_R8POLY converts a difference table to a polynomial;");
        Console.WriteLine("  DIF_SHIFT_ZERO shifts a divided difference table to 0;");
        Console.WriteLine("");
        Console.WriteLine("  These are equivalent operations");
        Console.WriteLine("");
        //
        //  Set XTAB, YTAB to X, F(X)
        //
        ntab = 4;
        for (i = 0; i < ntab; i++)
        {
            x = i + 1;
            xtab1[i] = x;
            xtab2[i] = x;
            ytab1[i] = -4.0 + x * (3.0 + x * (-2.0 + x));
            ytab2[i] = ytab1[i];
        }

        //
        //  Compute and display the finite difference table.
        //
        Data.data_to_dif_display(ntab, xtab1, ytab1, ref diftab1);

        Data.data_to_dif_display(ntab, xtab2, ytab2, ref diftab2);
        //
        //  Examine the corresponding polynomial.
        //
        Dif.dif_print(ntab, xtab1, diftab1, "  The divided difference table:");
        //
        //  Shift to zero using DIF_SHIFT_ZERO.
        //
        Dif.dif_shift_zero(ntab, ref xtab1, ref diftab1);

        typeMethods.r8poly_print(ntab, diftab1, "  The polynomial using DIF_SHIFT_ZERO:");
        //
        //  Shift to zero using DIF_TO_R8POLY.
        //
        Dif.dif_to_r8poly(ntab, xtab2, diftab2, ref c);

        typeMethods.r8poly_print(ntab, c, "  The polynomial using DIF_TO_R8POLY:");
    }

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tests R8POLY*.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 5;

        int i;
        double[] poly_cof = new double[N];
        double[] poly_cof2 = new double[N + 1];
        double[] poly_cof3 = new double[N - 1];
        double xval;
        double yval;
        double yval2;
        double yval3;

        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  R8POLY_ANT_COF computes the coefficients of the");
        Console.WriteLine("  antiderivative of a polynomial;");
        Console.WriteLine("  R8POLY_ANT_VAL evaluates the antiderivative of");
        Console.WriteLine("  a polynomial;");
        Console.WriteLine("  R8POLY_DER_COF computes the coefficients of the");
        Console.WriteLine("  derivative of a polynomial;");
        Console.WriteLine("  R8POLY_DER_VAL evaluates the derivative of");
        Console.WriteLine("  a polynomial;");
        Console.WriteLine("  R8POLY_PRINT prints a polynomial;");
        Console.WriteLine("  R8POLY_VAL evaluates a polynomial.");

        for (i = 0; i < N; i++)
        {
            poly_cof[i] = i + 1;
        }

        typeMethods.r8poly_print(N, poly_cof, "  The initial polynomial:");

        typeMethods.r8poly_ant_cof(N, poly_cof, ref poly_cof2);

        typeMethods.r8poly_print(N + 1, poly_cof2, "  The antiderivative polynomial:");

        typeMethods.r8poly_der_cof(N, poly_cof, ref poly_cof3);

        typeMethods.r8poly_print(N - 1, poly_cof3, "  The derivative polynomial:");

        Console.WriteLine("");
        Console.WriteLine("  Evaluate the polynomial, antiderivative and");
        Console.WriteLine("  derivative, using only the original polynomial");
        Console.WriteLine("  coefficients:");
        Console.WriteLine("");
        Console.WriteLine("  X   P(X)   Anti_P(X)     P'(X)");
        Console.WriteLine("");

        for (i = 0; i <= 2; i++)
        {

            xval = i;


            yval = typeMethods.r8poly_val_horner(N, poly_cof, xval);

            yval2 = typeMethods.r8poly_ant_val(N, poly_cof, xval);

            yval3 = typeMethods.r8poly_der_val(N, poly_cof, xval);

            Console.WriteLine("  "
                              + xval.ToString().PadLeft(10) + "  "
                              + yval.ToString().PadLeft(10) + "  "
                              + yval2.ToString().PadLeft(10) + "  "
                              + yval3.ToString().PadLeft(10) + "");

        }
    }

    private static void r8poly_basis_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    typeMethods.r8poly_basis_test tests typeMethods.r8poly_basis().
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int NTAB = 5;

        int i;
        int j;
        int nstep = 9;
        double[] pointer;
        int pointerIndex = 0;
        double[] polcof;
        int set = 3;
        double xhi;
        double xlo;
        double[] xtab = new double[NTAB];
        double xval;
        double yval;

        Console.WriteLine("");
        Console.WriteLine("typeMethods.r8poly_basis_test");
        Console.WriteLine("  typeMethods.r8poly_basis() computes Lagrange basis polynomials");
        Console.WriteLine("  in standard form.");
        Console.WriteLine("");

        polcof = new double[NTAB * NTAB];
        //
        //  Set the base points.
        //
        typeMethods.r8vec_indicator(NTAB, ref xtab);
        //
        //  Get the difference tables for the basis polynomials and print them.
        //
        typeMethods.r8poly_basis(NTAB, xtab, ref polcof);

        Console.WriteLine("");
        Console.WriteLine("  The table of difference vectors defining the basis polynomials.");
        Console.WriteLine("  Each ROW represents a polynomial.");
        Console.WriteLine("");

        pointer = polcof;
            
        string cout = "";

        for (i = 0; i < NTAB; i++)
        {
            for (j = 0; j < NTAB; j++)
            {
                cout += pointerIndex.ToString().PadLeft(10) + "  ";
                pointerIndex++;
            }

            Console.WriteLine(cout);
        }

        //
        //  Advance the pointer to the beginning of the data for basis polynomial SET.
        //
        pointer = polcof;
        pointerIndex = 0;
        for (i = 1; i < set; i++)
        {
            for (j = 1; j <= NTAB; j++)
            {
                pointerIndex++;
            }
        }

        //
        //  Print basis polynomial SET in polynomial form.
        //
        typeMethods.r8poly_print(NTAB, pointer, "  One basis polynomial in standard form:", aIndex: pointerIndex);
        //
        //  Evaluate basis polynoimial SET at a set of points.
        //
        Console.WriteLine("");
        Console.WriteLine("  Evaluate basis polynomial #" + set + " at a set of points.");
        Console.WriteLine("");
        Console.WriteLine("      X        Y");
        Console.WriteLine("");

        xhi = NTAB;
        xlo = 1.0;

        for (i = 1; i <= nstep; i++)
        {

            xval = ((nstep - i) * xlo
                    + (i - 1) * xhi)
                   / (nstep - 1);

            yval = typeMethods.r8poly_val_horner(NTAB, pointer, xval, polyCofIndex:pointerIndex);

            Console.WriteLine("  "
                              + xval.ToString().PadLeft(10) + "  "
                              + yval.ToString().PadLeft(10) + "");

        }
    }

    private static void r8poly_shift_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    typeMethods.r8poly_shift_test tests typeMethods.r8poly_shift().
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 3;

        int i;
        double[] poly_cof = new double[N];
        double scale;
        double shift;

        scale = 2.0;
        shift = +3.0;
        poly_cof[0] = +6.0;
        poly_cof[1] = -1.0;
        poly_cof[2] = 2.0;

        Console.WriteLine("");
        Console.WriteLine("typeMethods.r8poly_shift_test");
        Console.WriteLine("  typeMethods.r8poly_shift shifts polynomial coefficients.");
        Console.WriteLine("");
        Console.WriteLine("  Polynomial coefficients for argument X");
        Console.WriteLine("");

        for (i = 0; i < N; i++)
        {
            Console.WriteLine("  "
                              + i.ToString().PadLeft(6) + "  "
                              + poly_cof[i].ToString().PadLeft(10) + "");
        }

        typeMethods.r8poly_shift(scale, shift, N, ref poly_cof);

        Console.WriteLine("");
        Console.WriteLine("  SCALE = " + scale + "");
        Console.WriteLine("  SHIFT = " + shift + "");
        Console.WriteLine("");
        Console.WriteLine("  Polynomial coefficients for argument Z = SCALE * X + SHIFT:");
        Console.WriteLine("");

        for (i = 0; i < N; i++)
        {
            Console.WriteLine("  "
                              + i.ToString().PadLeft(6) + "  "
                              + poly_cof[i].ToString().PadLeft(10) + "");
        }
    }

    private static void ncc_rule_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ncc_rule_test tests ncc_rule();
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int NORDER = 8;

        int i;
        double[] weight = new double[NORDER];
        double[] xtab = new double[NORDER];

        Console.WriteLine("");
        Console.WriteLine("ncc_rule_test");
        Console.WriteLine("  ncc_rule() computes closed Newton Cotes formulas;");
        Console.WriteLine("");

        NewtonCotesQuadrature.ncc_rule(NORDER, ref xtab, ref weight);

        Console.WriteLine("");
        Console.WriteLine("  Newton-Cotes Closed Quadrature Rule:");
        Console.WriteLine("");
        Console.WriteLine("    Abscissa       Weight");
        Console.WriteLine("");

        for (i = 0; i < NORDER; i++)
        {
            Console.WriteLine("  "
                              + (i + 1).ToString().PadLeft(6) + "  "
                              + xtab[i].ToString().PadLeft(10) + "  "
                              + weight[i].ToString().PadLeft(10) + "");
        }
    }

    private static void nco_rule_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    nco_rule_test tests nco_rule().
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 October 2019
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int NORDER = 8;

        int i;
        double[] weight = new double[NORDER];
        double[] xtab = new double[NORDER];

        Console.WriteLine("");
        Console.WriteLine("nco_rule_test");
        Console.WriteLine("  nco_rule() computes open Newton Cotes formulas.");
        Console.WriteLine("");

        NewtonCotesQuadrature.nco_rule(NORDER, ref xtab, ref weight);

        Console.WriteLine("");
        Console.WriteLine("  Newton-Cotes Open Quadrature Rule:");
        Console.WriteLine("");
        Console.WriteLine("    Abscissa       Weight");
        Console.WriteLine("");

        for (i = 0; i < NORDER; i++)
        {
            Console.WriteLine("  "
                              + (i + 1).ToString().PadLeft(6) + "  "
                              + xtab[i].ToString().PadLeft(10) + "  "
                              + weight[i].ToString().PadLeft(10) + "");
        }
    }

    private static void test18()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST18 tests ROOTS_TO_DIF and DIF_TO_R8POLY.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int MAXROOTS = 4;

        double[] c = new double[MAXROOTS + 1];
        double[] diftab = new double[MAXROOTS + 1];
        int nroots;
        int ntab = 0;
        double[] roots = new double[MAXROOTS];
        double[] xtab = new double[MAXROOTS + 1];

        Console.WriteLine("");
        Console.WriteLine("TEST18");
        Console.WriteLine("  ROOTS_TO_DIF computes a divided difference");
        Console.WriteLine("  polynomial with the given roots;");
        Console.WriteLine("  DIF_TO_R8POLY converts it to a standard form polynomial.");
        Console.WriteLine("");

        nroots = 1;
        roots[0] = 3.0;
        typeMethods.r8vec_print(nroots, roots, "  The roots:");

        Roots.roots_to_dif(nroots, roots, ref ntab, ref xtab, ref diftab);
        Dif.dif_to_r8poly(ntab, xtab, diftab, ref c);
        typeMethods.r8poly_print(ntab, c, "  The polynomial:");

        nroots = 2;
        roots[0] = 3.0;
        roots[1] = 1.0;
        typeMethods.r8vec_print(nroots, roots, "  The roots:");

        Roots.roots_to_dif(nroots, roots, ref ntab, ref xtab, ref diftab);
        Dif.dif_to_r8poly(ntab, xtab, diftab, ref c);
        typeMethods.r8poly_print(ntab, c, "  The polynomial:");

        nroots = 3;
        roots[0] = 3.0;
        roots[1] = 1.0;
        roots[2] = 2.0;
        typeMethods.r8vec_print(nroots, roots, "  The roots:");

        Roots.roots_to_dif(nroots, roots, ref ntab, ref xtab, ref diftab);
        Dif.dif_to_r8poly(ntab, xtab, diftab, ref c);
        typeMethods.r8poly_print(ntab, c, "  The polynomial:");

        nroots = 4;
        roots[0] = 3.0;
        roots[1] = 1.0;
        roots[2] = 2.0;
        roots[3] = 4.0;
        typeMethods.r8vec_print(nroots, roots, "  The roots:");

        Roots.roots_to_dif(nroots, roots, ref ntab, ref xtab, ref diftab);
        Dif.dif_to_r8poly(ntab, xtab, diftab, ref c);
        typeMethods.r8poly_print(ntab, c, "  The polynomial:");
    }

    private static void roots_to_r8poly_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    roots_to_typeMethods.r8poly_test tests roots_to_r8poly();
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int MAXROOT = 5;

        double[] c = new double[MAXROOT + 1];
        int nc = 0;
        int nroot;
        double[] roots = new double[MAXROOT];

        Console.WriteLine("");
        Console.WriteLine("roots_to_typeMethods.r8poly_test");
        Console.WriteLine("  roots_to_r8poly() computes polynomial coefficients from roots.");
        Console.WriteLine("");

        nroot = 1;
        roots[0] = 3.0;
        typeMethods.r8vec_print(nroot, roots, "  The roots:");

        Roots.roots_to_r8poly(nroot, roots, ref nc, ref c);

        typeMethods.r8poly_print(nc, c, "  The polynomial:");

        nroot = 2;
        roots[0] = 3.0;
        roots[1] = 1.0;
        typeMethods.r8vec_print(nroot, roots, "  The roots:");

        Roots.roots_to_r8poly(nroot, roots, ref nc, ref c);

        typeMethods.r8poly_print(nc, c, "  The polynomial:");

        nroot = 3;
        roots[0] = 3.0;
        roots[1] = 1.0;
        roots[2] = 2.0;
        typeMethods.r8vec_print(nroot, roots, "  The roots:");

        Roots.roots_to_r8poly(nroot, roots, ref nc, ref c);

        typeMethods.r8poly_print(nc, c, "  The polynomial:");

        nroot = 4;
        roots[0] = 3.0;
        roots[1] = 1.0;
        roots[2] = 2.0;
        roots[3] = 4.0;
        typeMethods.r8vec_print(nroot, roots, "  The roots:");

        Roots.roots_to_r8poly(nroot, roots, ref nc, ref c);

        typeMethods.r8poly_print(nc, c, "  The polynomial:");
    }

    private static void dif_derivk_table_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    Dif.dif_derivk_table_test tests Dif.dif_derivk_table;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] c0;
        double[] d0;
        double[] d1;
        double[] d2;
        double[] d3;
        double[] d4;
        double[] f0;
        int i;
        int j;
        int k;
        int n0;
        int n1;
        int n2;
        int n3;
        int n4;
        double x;
        double[] x0;
        double[] x1;
        double[] x2;
        double[] x3;
        double[] x4;
        double y0;
        double y1;
        double y2;
        double y3;
        double y4;

        Console.WriteLine("");
        Console.WriteLine("Dif.dif_derivk_table_test");
        Console.WriteLine("  Dif.dif_derivk_table() computes the K-th derivative");
        Console.WriteLine("  for a divided difference table.");
        //
        //  Set the 0 data points.
        //
        n0 = 5;
        x0 = new double[n0];
        f0 = new double[n0];
        for (i = 0; i < 5; i++)
        {
            x0[i] = i - 2;
        }

        //
        //  Set data for x^4/24+x^3/3+x^2/2+x+1
        //
        for (j = 0; j < n0; j++)
        {
            f0[j] = 1.0;
            for (i = 4; 1 <= i; i--)
            {
                f0[j] = f0[j] * x0[j] / i + 1.0;
            }
        }

        //
        //  Compute the difference table.
        //
        d0 = new double[n0];
        Data.data_to_dif(n0, x0, f0, ref d0);
        Dif.dif_print(n0, x0, d0, "  The divided difference polynomial P0:");

        c0 = new double[n0];
        Dif.dif_to_r8poly(n0, x0, d0, ref c0);

        typeMethods.r8poly_print(n0, c0, "  Using DIF_TO_R8POLY");
        //
        //  Compute the difference table for the K=1 derivative.
        //
        k = 1;
        n1 = n0 - k;
        x1 = new double[n1];
        d1 = new double[n1];
        Dif.dif_derivk_table(n0, x0, d0, k, x1, ref d1);
        //
        //  Compute the difference table for the K=2 derivative.
        //
        k = 2;
        n2 = n0 - k;
        x2 = new double[n2];
        d2 = new double[n2];
        Dif.dif_derivk_table(n0, x0, d0, k, x2, ref d2);
        //
        //  Compute the difference table for the K=3 derivative.
        //
        k = 3;
        n3 = n0 - k;
        x3 = new double[n3];
        d3 = new double[n3];
        Dif.dif_derivk_table(n0, x0, d0, k, x3, ref d3);
        //
        //  Compute the difference table for the K=4 derivative.
        //
        k = 4;
        n4 = n0 - k;
        x4 = new double[n4];
        d4 = new double[n4];
        Dif.dif_derivk_table(n0, x0, d0, k, x4, ref d4);
        //
        //  Evaluate all 5 polynomials.
        //
        Console.WriteLine("");
        Console.WriteLine("  Evaluate difference tables for the function P0");
        Console.WriteLine("  and its first four derivatives, P1...P4.");
        Console.WriteLine("");
        Console.WriteLine("      X         P0        P1        P2        P3        P4");
        Console.WriteLine("");

        for (i = 0; i <= 10; i++)
        {
            x = i / 5.0;
            y0 = Dif.dif_val(n0, x0, d0, x);
            y1 = Dif.dif_val(n1, x1, d1, x);
            y2 = Dif.dif_val(n2, x2, d2, x);
            y3 = Dif.dif_val(n3, x3, d3, x);
            y4 = Dif.dif_val(n4, x4, d4, x);
            Console.WriteLine("  " + x.ToString().PadLeft(8)
                                   + "  " + y0.ToString().PadLeft(8)
                                   + "  " + y1.ToString().PadLeft(8)
                                   + "  " + y2.ToString().PadLeft(8)
                                   + "  " + y3.ToString().PadLeft(8)
                                   + "  " + y4.ToString().PadLeft(8) + "");
        }
    }

    private static void dif_basis_deriv_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    Dif.dif_basis_deriv_test tests Dif.dif_basis_deriv;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] c;
        double[] ddp;
        int nd = 3;
        double[] xd;
        double[] xdp;

        Console.WriteLine("");
        Console.WriteLine("Dif.dif_basis_deriv_test");
        Console.WriteLine("  Dif.dif_basis_deriv() computes difference tables for");
        Console.WriteLine("  the first derivative of each Lagrange basis.");

        xd = new double[nd];
        xdp = new double[nd - 1];
        ddp = new double[(nd - 1) * nd];

        xd[0] = -2.0;
        xd[1] = 1.0;
        xd[2] = 5.0;

        Dif.dif_basis_deriv(nd, xd, ref xdp, ref ddp);
        //
        //  Because the difference tables were shifted to all 0 abscissas,
        //  they contain the polynomial coefficients.
        //
        typeMethods.r8mat_transpose_print(nd - 1, nd, ddp,
            "  Lagrange basis derivative polynomial coefficients:");

        c = new double[nd - 1];
        Dif.dif_to_r8poly(nd - 1, xdp, ddp, ref c, diftabIndex: + 0 * (nd - 1));
        typeMethods.r8poly_print(nd - 1, c, "  P1'=-(2x-6)/21");

        Dif.dif_to_r8poly(nd - 1, xdp, ddp, ref c, diftabIndex:  + 1 * (nd - 1));
        typeMethods.r8poly_print(nd - 1, c, "  P2'=-(2x-3)/12");

        Dif.dif_to_r8poly(nd - 1, xdp, ddp, ref c, diftabIndex:  + 2 * (nd - 1));
        typeMethods.r8poly_print(nd - 1, c, "  P3'=(2x+1)/28");
    }

    private static void dif_basis_derivk_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    Dif.dif_basis_derivk_test tests DIF_BASIS_DERIVK;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] c;
        double[] ddp;
        int k = 2;
        int nd = 5;
        double[] xd;
        double[] xdp;

        Console.WriteLine("");
        Console.WriteLine("Dif.dif_basis_derivk_test");
        Console.WriteLine("  DIF_BASIS_DERIVK computes difference tables for");
        Console.WriteLine("  the K-th derivative of each Lagrange basis.");

        xd = new double[nd];
        xdp = new double[nd - k];
        ddp = new double[(nd - k) * nd];

        xd[0] = 1.0;
        xd[1] = 2.0;
        xd[2] = 3.0;
        xd[3] = 4.0;
        xd[4] = 5.0;

        Dif.dif_basis_derivk(nd, xd, k, ref xdp, ref ddp);
        //
        //  Because the difference tables were shifted to all 0 abscissas,
        //  they contain the polynomial coefficients.
        //
        typeMethods.r8mat_transpose_print(nd - k, nd, ddp,
            "  Lagrange basis K-th derivative polynomial coefficients:");

        c = new double[nd - k];

        Dif.dif_to_r8poly(nd - k, xdp, ddp, ref c, diftabIndex: + 0 * (nd - k));
        typeMethods.r8poly_print(nd - k, c, "  P1''=(12x^2-84x+142)/24");

        Dif.dif_to_r8poly(nd - k, xdp, ddp, ref c, diftabIndex:  + 1 * (nd - k));
        typeMethods.r8poly_print(nd - k, c, "  P2''=-2x^2+13x-59/3");

        Dif.dif_to_r8poly(nd - k, xdp, ddp, ref c, diftabIndex:  + 2 * (nd - k));
        typeMethods.r8poly_print(nd - k, c, "  P3''=3x^2-18x+49/2");

        Dif.dif_to_r8poly(nd - k, xdp, ddp, ref c, diftabIndex:  + 3 * (nd - k));
        typeMethods.r8poly_print(nd - k, c, "  P4''=-2x^2+11x-41/3");

        Dif.dif_to_r8poly(nd - k, xdp, ddp, ref c, diftabIndex:  + 4 * (nd - k));
        typeMethods.r8poly_print(nd - k, c, "  P5''=(6x^2-30x+35)/12");
    }
}