﻿using System;
using System.Globalization;
using Burkardt.PDFLib;
using Burkardt.Types;

namespace PDFLibTest;

internal static class Program
{
    private static void Main()
//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PDFLIB_TEST.
//
//  Discussion:
//
//    PDFLIB_TEST tests the PDFLIB library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 January 2018
//
//  Author:
//
//    John Burkardt
//
    {
        Console.WriteLine("");
        Console.WriteLine("PDFLIB_TEST");
        Console.WriteLine("  Test the PDFLIB library.");

        i4_binomial_pdf_test();
        i4_binomial_sample_test();
        i4_uniform_sample_test();
        r8_chi_sample_test();
        r8po_fa_test();
        r8vec_multinormal_pdf_test();

        Console.WriteLine("");
        Console.WriteLine("PDFLIB_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void i4_binomial_pdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    I4_BINOMIAL_PDF_TEST calls I4_BINOMIAL_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 January 2018
//
//  Author:
//
//    John Burkardt
//
    {
        Console.WriteLine("");
        Console.WriteLine("I4_BINOMIAL_PDF_TEST");
        Console.WriteLine("  I4_BINOMIAL_PDF reports");
        Console.WriteLine("  PROB, the probability that");
        Console.WriteLine("  N trials, with");
        Console.WriteLine("  P probability of success result in");
        Console.WriteLine("  K successes.");
        Console.WriteLine("");
        Console.WriteLine("   N         P   K        PROB");
        Console.WriteLine("");

        int n = 5;
        double p = 0.25;

        for (int k = 0; k <= n; k++)
        {
            double prob = PDF.i4_binomial_pdf(n, p, k);
            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + p.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + k.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + prob.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void i4_binomial_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    I4_BINOMIAL_SAMPLE_TEST calls I4_BINOMIAL_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 January 2018
//
//  Author:
//
//    John Burkardt
//
    {
        Console.WriteLine("");
        Console.WriteLine("I4_BINOMIAL_SAMPLE_TEST");
        Console.WriteLine("  I4_BINOMIAL_SAMPLE samples the binomial distribution.");
        Console.WriteLine("");
        Console.WriteLine("   N         P   K        PDF");
        Console.WriteLine("");

        for (int i = 1; i <= 10; i++)
        {
            int n = PDF.i4_uniform_sample(1, 20);
            double p = PDF.r8_uniform_sample(0.0, 1.0);
            int k = PDF.i4_binomial_sample(n, p);
            double pdf = PDF.i4_binomial_pdf(n, p, k);
            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + p.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + k.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void i4_uniform_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM_SAMPLE_TEST calls I4_UNIFORM_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 January 2018
//
//  Author:
//
//    John Burkardt
//
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("I4_UNIFORM_SAMPLE_TEST");
        Console.WriteLine("  I4_UNIFORM_SAMPLE samples the uniform distribution on integers.");
        Console.WriteLine("  Generate integer C between limits A and B.");
        Console.WriteLine("");
        Console.WriteLine("    A    B   C");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            int a = PDF.i4_uniform_sample(-10, +10);
            int b = PDF.i4_uniform_sample(a, 20);
            int c = PDF.i4_uniform_sample(a, b);
            Console.WriteLine("  " + a.ToString(CultureInfo.InvariantCulture).PadLeft(3)
                                   + "  " + b.ToString(CultureInfo.InvariantCulture).PadLeft(3)
                                   + "  " + c.ToString(CultureInfo.InvariantCulture).PadLeft(3) + "");
        }
    }

    private static void r8_chi_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    R8_CHI_SAMPLE_TEST calls R8_CHI_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
    {
        Console.WriteLine("");
        Console.WriteLine("R8_CHI_SAMPLE_TEST");
        Console.WriteLine("  R8_CHI_SAMPLE ( DF ) samples the Chi distribution with");
        Console.WriteLine("  DF degrees of freedom.");
//
//  Set the current generator index to #2, which (this being C#) has index 1!.
//
        int g = 1;
        // cgn_set(g);
        Console.WriteLine("");
        Console.WriteLine("  Current generator index = " + g + "");
//
//  Repeatedly call R8_CHI_SAMPLE ( DF ).
//
        Console.WriteLine("");
        Console.WriteLine("   I       DF       R8_CHI_SAMPLE ( DF )");
        Console.WriteLine("");

        for (int i = 0; i <= 10; i++)
        {
            double df = 5.0 * PDF.r8_uniform_01_sample() + 1.0;
            double u = PDF.r8_chi_sample(df);
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + df.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + u.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void r8po_fa_test()

//****************************************************************************80
//
//  Purpose:
//
//    R8PO_FA_TEST tests R8PO_FA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 June 2013
//
//  Author:
//
//    John Burkardt
//
    {
        int j;
        const int n = 5;

        Console.WriteLine("");
        Console.WriteLine("R8PO_FA_TEST");
        Console.WriteLine("  R8PO_FA computes the Cholesky factor R of a");
        Console.WriteLine("  positive definite matrix A, so that A = R' * R.");
        Console.WriteLine("");
        Console.WriteLine("  Start with random R1;");
        Console.WriteLine("  Compute A = R1' * R1.");
        Console.WriteLine("  Call R8MAT_POFAC and see if you recover R2 = R1.");
//
//  Generate a random upper triangular matrix with positive diagonal.
//
        double[] r1 = new double[n * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i <= j; i++)
            {
                r1[i + j * n] = PDF.r8_uniform_01_sample();
            }

            for (i = j + 1; i < n; i++)
            {
                r1[i + j * n] = 0.0;
            }
        }

        typeMethods.r8ge_print(n, n, r1, "  R1:");
//
//  Compute a positive definite symmetric matrix A.
//
        double[] a = typeMethods.r8ge_mtm(n, r1, r1);

        typeMethods.r8ge_print(n, n, a, "  A:");

        double[] r2 = typeMethods.r8po_fa(n, a);

        double diff = typeMethods.r8mat_norm_fro_affine(n, n, r1, r2);

        Console.WriteLine("");
        Console.WriteLine("  Frobenius difference between R1 and R2 = " + diff + "");
    }

    public static void r8vec_multinormal_pdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MULTINORMAL_PDF_TEST tests R8VEC_MULTINORMAL_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 August 2015
//
//  Author:
//
//    John Burkardt
//
    {
        int i;
        int j;
        const int n = 5;

        Console.WriteLine("");
        Console.WriteLine("R8VEC_MULTINORMAL_PDF_TEST");
        Console.WriteLine("  R8VEC_MULTINORMAL_PDF evaluates the PDF for the");
        Console.WriteLine("  multinormal distribution.");
        Console.WriteLine("");
        Console.WriteLine("  The covariance matrix is C.");
        Console.WriteLine("  The definition uses the inverse of C;");
        Console.WriteLine("  R8VEC_MULTINORMAL_PDF uses the Cholesky factor.");
        Console.WriteLine("  Verify that the algorithms are equivalent.");
//
//  Generate a random upper triangular matrix with positive diagonal.
//
        double[] r1 = new double[n * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i <= j; i++)
            {
                if (i == j)
                {
                    r1[i + j * n] = Math.Abs(PDF.r8_uniform_01_sample());
                }
                else
                {
                    r1[i + j * n] = PDF.r8_uniform_01_sample();
                }
            }

            for (i = j + 1; i < n; i++)
            {
                r1[i + j * n] = 0.0;
            }
        }

        typeMethods.r8ge_print(n, n, r1, "  R1:");
//
//  Compute a positive definite symmetric matrix C.
//
        double[] c = typeMethods.r8ge_mtm(n, r1, r1);
        typeMethods.r8ge_print(n, n, c, "  C:");
//
//  Compute the Cholesky factor.
//
        double[] r2 = typeMethods.r8mat_pofac(n, c);
        typeMethods.r8ge_print(n, n, r2, "  R2:");
//
//  Compute the determinant of C.
//
        double c_det = typeMethods.r8mat_podet(n, r2);
        Console.WriteLine("");
        Console.WriteLine("  Determinant of C = " + c_det + "");
//
//  Compute the inverse of C.
//
        double[] c_inv = typeMethods.r8mat_poinv(n, r2);
//
//  Compute a random set of means.
//
        double[] mu = new double[n];
        for (i = 0; i < n; i++)
        {
            mu[i] = PDF.r8_normal_01_sample();
        }

//
//  Compute X as small variations from MU.
//
        double[] x = new double[n];
        for (i = 0; i < n; i++)
        {
            double eps = 0.01 * PDF.r8_normal_01_sample();
            x[i] = (1.0 + eps) * mu[i];
        }

//
//  Compute PDF1 from the function.
//
        double pdf1 = PDF.r8vec_multinormal_pdf(n, mu, r2, c_det, x);
//
//  Compute PDF2 from the definition.
//
        double[] y = new double[n];
        for (i = 0; i < n; i++)
        {
            y[i] = x[i] - mu[i];
        }

        double xcx = 0.0;
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                if (i <= j)
                {
                    xcx += y[i] * c_inv[i + j * n] * y[j];
                }
                else
                {
                    xcx += y[i] * c_inv[j + i * n] * y[j];
                }
            }
        }

        double pdf2 = 1.0 / Math.Sqrt(Math.Pow(2.0 * Math.PI, n))
                      * 1.0 / Math.Sqrt(c_det)
                      * Math.Exp(-0.5 * xcx);

        Console.WriteLine("");
        Console.WriteLine("  PDF1 = " + pdf1 + "");
        Console.WriteLine("  PDF2 = " + pdf2 + "");
            
    }
}