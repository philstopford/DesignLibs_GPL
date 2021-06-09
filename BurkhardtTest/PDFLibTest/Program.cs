using System;
using Burkardt.PDFLib;
using Burkardt.Types;

namespace PDFLibTest
{
    class Program
    {
        static void Main(string[] args)
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

        static void i4_binomial_pdf_test()

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
                Console.WriteLine("  " + n.ToString().PadLeft(2)
                                  + "  " + p.ToString().PadLeft(8)
                                  + "  " + k.ToString().PadLeft(2)
                                  + "  " + prob.ToString().PadLeft(14) + "");
            }
        }

        static void i4_binomial_sample_test()

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
                Console.WriteLine("  " + n.ToString().PadLeft(2)
                    + "  " + p.ToString().PadLeft(8)
                    + "  " + k.ToString().PadLeft(2)
                    + "  " + pdf.ToString().PadLeft(14) + "");
            }
        }

        static void i4_uniform_sample_test()

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
            int a;
            int b;
            int c;
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
                a = PDF.i4_uniform_sample(-10, +10);
                b = PDF.i4_uniform_sample(a, 20);
                c = PDF.i4_uniform_sample(a, b);
                Console.WriteLine("  " + a.ToString().PadLeft(3)
                    + "  " + b.ToString().PadLeft(3)
                    + "  " + c.ToString().PadLeft(3) + "");
            }

            return;
        }

        static void r8_chi_sample_test()

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
                Console.WriteLine("  " + i.ToString().PadLeft(2)
                    + "  " + df.ToString().PadLeft(14)
                    + "  " + u.ToString().PadLeft(14) + "");
            }
        }

        static void r8po_fa_test()

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
            double[] a;
            double diff;
            int i;
            int j;
            int n = 5;
            double[] r1;
            double[] r2;

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
            r1 = new double[n * n];

            for (j = 0; j < n; j++)
            {
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
            a = typeMethods.r8ge_mtm(n, r1, r1);

            typeMethods.r8ge_print(n, n, a, "  A:");

            r2 = typeMethods.r8po_fa(n, a);

            diff = typeMethods.r8mat_norm_fro_affine(n, n, r1, r2);

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
            double[] c;
            double c_det;
            double[] c_inv;
            double eps;
            int i;
            int j;
            double[] mu;
            int n = 5;
            double pdf1;
            double pdf2;
            const double r8_pi = 3.141592653589793;
            double[] r1;
            double[] r2;
            double[] x;
            double xcx;
            double[] y;

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
            r1 = new double[n * n];

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
            c = typeMethods.r8ge_mtm(n, r1, r1);
            typeMethods.r8ge_print(n, n, c, "  C:");
//
//  Compute the Cholesky factor.
//
            r2 = typeMethods.r8mat_pofac(n, c);
            typeMethods.r8ge_print(n, n, r2, "  R2:");
//
//  Compute the determinant of C.
//
            c_det = typeMethods.r8mat_podet(n, r2);
            Console.WriteLine("");
            Console.WriteLine("  Determinant of C = " + c_det + "");
//
//  Compute the inverse of C.
//
            c_inv = typeMethods.r8mat_poinv(n, r2);
//
//  Compute a random set of means.
//
            mu = new double[n];
            for (i = 0; i < n; i++)
            {
                mu[i] = PDF.r8_normal_01_sample();
            }

//
//  Compute X as small variations from MU.
//
            x = new double[n];
            for (i = 0; i < n; i++)
            {
                eps = 0.01 * PDF.r8_normal_01_sample();
                x[i] = (1.0 + eps) * mu[i];
            }

//
//  Compute PDF1 from the function.
//
            pdf1 = PDF.r8vec_multinormal_pdf(n, mu, r2, c_det, x);
//
//  Compute PDF2 from the definition.
//
            y = new double[n];
            for (i = 0; i < n; i++)
            {
                y[i] = x[i] - mu[i];
            }

            xcx = 0.0;
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    if (i <= j)
                    {
                        xcx = xcx + y[i] * c_inv[i + j * n] * y[j];
                    }
                    else
                    {
                        xcx = xcx + y[i] * c_inv[j + i * n] * y[j];
                    }
                }
            }

            pdf2 = 1.0 / Math.Sqrt(Math.Pow(2.0 * r8_pi, n))
                   * 1.0 / Math.Sqrt(c_det)
                   * Math.Exp(-0.5 * xcx);

            Console.WriteLine("");
            Console.WriteLine("  PDF1 = " + pdf1 + "");
            Console.WriteLine("  PDF2 = " + pdf2 + "");
            
        }
    }
}