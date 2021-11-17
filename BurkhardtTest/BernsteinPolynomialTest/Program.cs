using System;
using Burkardt;
using Burkardt.PolynomialNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace BernsteinPolynomialTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for BERNSTEIN_POLYNOMIAL_TEST.
        //
        //  Discussion:
        //
        //    BERNSTEIN_POLYNOMIAL_TEST tests the BERNSTEIN_POLYNOMIAL library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 March 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("BERNSTEIN_POLYNOMIAL_TEST");
        Console.WriteLine("  Test the BERNSTEIN_POLYNOMIAL library.");

        bernstein_matrix_test();
        bernstein_matrix_determinant_test();
        bernstein_matrix_inverse_test();
        bernstein_poly_01_test();
        bernstein_poly_01_test2();
        bernstein_poly_01_matrix_test();
        bernstein_poly_ab_test();
        bernstein_poly_ab_approx_test();
        bernstein_to_legendre_test();
        bernstein_to_power_test();
        bernstein_vandermonde_test();

        Console.WriteLine("");
        Console.WriteLine("BERNSTEIN_POLYNOMIAL_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine(" ");
    }

    private static void bernstein_matrix_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNSTEIN_MATRIX_TEST tests BERNSTEIN_MATRIX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("BERNSTEIN_MATRIX_TEST");
        Console.WriteLine("  BERNSTEIN_MATRIX returns a matrix A which transforms a");
        Console.WriteLine("  polynomial coefficient vector from the power basis to");
        Console.WriteLine("  the Bernstein basis.");

        int n = 5;
        double[] a = BernsteinPolynomial.bernstein_matrix(n);
        typeMethods.r8mat_print(n, n, a, "  Bernstein matrix A of order 5:");
    }

    private static void bernstein_matrix_test2()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNSTEIN_MATRIX_TEST2 tests BERNSTEIN_MATRIX.
        //
        //  Discussion:
        //
        //    Here we use the Bernstein matrix to describe a Bernstein polynomial
        //    in terms of the standard monomials.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  BERNSTEIN_MATRIX returns a matrix A which");
        Console.WriteLine("  transforms a polynomial coefficient vector");
        Console.WriteLine("  from the the Bernstein basis to the power basis.");
        Console.WriteLine("  We can use this to get explicit values of the");
        Console.WriteLine("  4-th degree Bernstein polynomial coefficients as");
        Console.WriteLine("");
        Console.WriteLine("    b(4,K)(X) = C4 * x^4");
        Console.WriteLine("              + C3 * x^3");
        Console.WriteLine("              + C2 * x^2");
        Console.WriteLine("              + C1 * x");
        Console.WriteLine("              + C0 * 1");

        int n = 5;
        Console.WriteLine("");
        Console.WriteLine("     K       C4           C3            C2" +
                          "            C1             C0");
        Console.WriteLine("");

        double[] a = BernsteinPolynomial.bernstein_matrix(n);
        double[] x = new double[n];

        for (int k = 0; k < n; k++)
        {
            for (int i = 0; i < n; i++)
            {
                x[i] = 0.0;
            }

            x[k] = 1.0;

            double[] ax = typeMethods.r8mat_mv_new(n, n, a, x);

            string cout = "  " + k.ToString().PadLeft(4) + "  ";
            for (int i = 0; i < n; i++)
            {
                cout += "  " + ax[i].ToString().PadLeft(14);
            }

            Console.WriteLine(cout);
        }
    }

    private static void bernstein_matrix_determinant_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNSTEIN_MATRIX_DETERMINANT_TEST tests BERNSTEIN_MATRIX_DETERMINANT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("BERNSTEIN_MATRIX_DETERMINANT_TEST");
        Console.WriteLine("  BERNSTEIN_MATRIX_DETERMINANT computes the determinant of");
        Console.WriteLine("  the Bernstein matrix.");
        Console.WriteLine("");
        Console.WriteLine("     N         ||A||          det(A)");
        Console.WriteLine("                              computed");
        Console.WriteLine("");

        for (int n = 5; n <= 15; n++)
        {
            double[] a = BernsteinPolynomial.bernstein_matrix(n);
            double a_norm_frobenius = typeMethods.r8mat_norm_fro(n, n, a);

            double d1 = BernsteinPolynomial.bernstein_matrix_determinant(n);

            Console.WriteLine("  " + n.ToString().PadLeft(4)
                                   + "  " + a_norm_frobenius.ToString().PadLeft(14)
                                   + "  " + d1.ToString().PadLeft(14) + "");

        }

    }

    private static void bernstein_matrix_inverse_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNSTEIN_MATRIX_INVERSE_TEST tests BERNSTEIN_MATRIX_INVERSE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("BERNSTEIN_MATRIX_INVERSE_TEST");
        Console.WriteLine("  BERNSTEIN_MATRIX_INVERSE computes the inverse of the");
        Console.WriteLine("  Bernstein matrix A.");
        Console.WriteLine("");
        Console.WriteLine("     N     ||A||       ||inv(A)|| ||I-A*inv(A)||");
        Console.WriteLine("");

        for (int n = 5; n <= 15; n++)
        {
            double[] a = BernsteinPolynomial.bernstein_matrix(n);
            double a_norm_frobenius = typeMethods.r8mat_norm_fro(n, n, a);

            double[] b = BernsteinPolynomial.bernstein_matrix_inverse(n);
            double b_norm_frobenius = typeMethods.r8mat_norm_fro(n, n, b);

            double[] c = typeMethods.r8mat_mm_new(n, n, n, a, b);
            double error_norm_frobenius = typeMethods.r8mat_is_identity(n, c);

            Console.WriteLine("  " + n.ToString().PadLeft(4)
                                   + "  " + a_norm_frobenius.ToString().PadLeft(14)
                                   + "  " + b_norm_frobenius.ToString().PadLeft(14)
                                   + "  " + error_norm_frobenius.ToString().PadLeft(14) + "");

        }
    }

    private static void bernstein_poly_01_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNSTEIN_POLY_01_TEST tests BERNSTEIN_POLY_01.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 July 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double b = 0;
        int k = 0;
        int n = 0;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("BERNSTEIN_POLY_01_TEST:");
        Console.WriteLine("  BERNSTEIN_POLY_01 evaluates the Bernstein polynomials");
        Console.WriteLine("  based on the interval [0,1].");
        Console.WriteLine("");
        Console.WriteLine("     N     K     X       Exact         BP01(N,K)(X)");
        Console.WriteLine("");

        int n_data = 0;

        while (true)
        {
            Burkardt.Values.Bernstein.bernstein_poly_01_values(ref n_data, ref n, ref k, ref x, ref b);

            if (n_data == 0)
            {
                break;
            }

            double[] bvec = BernsteinPolynomial.bernstein_poly_01(n, x);

            Console.WriteLine("  " + n.ToString().PadLeft(4)
                                   + "  " + k.ToString().PadLeft(4)
                                   + "  " + x.ToString().PadLeft(7)
                                   + "  " + b.ToString().PadLeft(14)
                                   + "  " + bvec[k].ToString().PadLeft(14) + "");
        }
    }

    private static void bernstein_poly_01_test2()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNSTEIN_POLY_01_TEST2 tests BERNSTEIN_POLY_01.
        //
        //  Discussion:
        //
        //    Here we test the Partition-of-Unity property.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("BERNSTEIN_POLY_01_TEST2:");
        Console.WriteLine("  BERNSTEIN_POLY_01 evaluates the Bernstein polynomials");
        Console.WriteLine("  based on the interval [0,1].");
        Console.WriteLine("");
        Console.WriteLine("  Here we test the partition of unity property.");
        Console.WriteLine("");
        Console.WriteLine("     N     X          Sum ( 0 <= K <= N ) BP01(N,K)(X)");
        Console.WriteLine("");

        int seed = 123456789;

        for (int n = 0; n <= 10; n++)
        {
            double x = UniformRNG.r8_uniform_01(ref seed);

            double[] bvec = BernsteinPolynomial.bernstein_poly_01(n, x);

            Console.WriteLine("  " + n.ToString().PadLeft(4)
                                   + "  " + x.ToString().PadLeft(7)
                                   + "  " + typeMethods.r8vec_sum(n + 1, bvec).ToString().PadLeft(14) + "");

        }
    }

    private static void bernstein_poly_01_matrix_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNSTEIN_POLY_01_MATRIX_TEST tests BERNSTEIN_POLY_01_MATRIX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("BERNSTEIN_POLY_01_MATRIX_TEST");
        Console.WriteLine("  BERNSTEIN_POLY_01_MATRIX is given M data values X,");
        Console.WriteLine("  and a degree N, and returns an Mx(N+1) matrix B such that");
        Console.WriteLine("  B(i,j) is the j-th Bernstein polynomial evaluated at the.");
        Console.WriteLine("  i-th data value.");

        int m = 5;
        double[] x = typeMethods.r8vec_linspace_new(m, 0.0, 1.0);
        int n = 1;
        double[] b = BernsteinPolynomial.bernstein_poly_01_matrix(m, n, x);
        typeMethods.r8mat_print(m, n + 1, b, "  B(5,1+1):");

        m = 5;
        x = typeMethods.r8vec_linspace_new(m, 0.0, 1.0);
        n = 4;
        b = BernsteinPolynomial.bernstein_poly_01_matrix(m, n, x);
        typeMethods.r8mat_print(m, n + 1, b, "  B(5,4+1):");

        m = 10;
        x = typeMethods.r8vec_linspace_new(m, 0.0, 1.0);
        n = 4;
        b = BernsteinPolynomial.bernstein_poly_01_matrix(m, n, x);
        typeMethods.r8mat_print(m, n + 1, b, "  B(10,4+1):");

        m = 3;
        x = typeMethods.r8vec_linspace_new(m, 0.0, 1.0);
        n = 5;
        b = BernsteinPolynomial.bernstein_poly_01_matrix(m, n, x);
        typeMethods.r8mat_print(m, n + 1, b, "  B(3,5+1):");
    }

    private static void bernstein_poly_ab_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNSTEIN_POLY_AB_TEST tests BERNSTEIN_POLY_AB.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 July 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n = 10;

        Console.WriteLine("");
        Console.WriteLine("BERNSTEIN_POLY_AB_TEST");
        Console.WriteLine("  BERNSTEIN_POLY_AB evaluates Bernstein polynomials over an");
        Console.WriteLine("  arbitrary interval [A,B].");
        Console.WriteLine("");
        Console.WriteLine("  Here, we demonstrate that ");
        Console.WriteLine("    BPAB(N,K,A1,B1)(X1) = BPAB(N,K,A2,B2)(X2)");
        Console.WriteLine("  provided only that");
        Console.WriteLine("    (X1-A1)/(B1-A1) = (X2-A2)/(B2-A2).");

        double x = 0.3;
        double a = 0.0;
        double b = 1.0;
        double[] bern = BernsteinPolynomial.bernstein_poly_ab(n, a, b, x);

        Console.WriteLine("");
        Console.WriteLine("     N     K     A        B        X       BPAB(N,K,A,B)(X)");
        Console.WriteLine("");
        for (int k = 0; k <= n; k++)
        {
            Console.WriteLine("  " + n.ToString().PadLeft(4)
                                   + "  " + k.ToString().PadLeft(4)
                                   + "  " + a.ToString().PadLeft(7)
                                   + "  " + b.ToString().PadLeft(7)
                                   + "  " + x.ToString().PadLeft(7)
                                   + "  " + bern[k].ToString().PadLeft(14) + "");
        }
            
        x = 1.3;
        a = 1.0;
        b = 2.0;
        bern = BernsteinPolynomial.bernstein_poly_ab(n, a, b, x);

        Console.WriteLine("");
        Console.WriteLine("     N     K     A        B        X       BPAB(N,K,A,B)(X)");
        Console.WriteLine("");
        for (int k = 0; k <= n; k++)
        {
            Console.WriteLine("  " + n.ToString().PadLeft(4)
                                   + "  " + k.ToString().PadLeft(4)
                                   + "  " + a.ToString().PadLeft(7)
                                   + "  " + b.ToString().PadLeft(7)
                                   + "  " + x.ToString().PadLeft(7)
                                   + "  " + bern[k].ToString().PadLeft(14) + "");
        }
            
        x = 2.6;
        a = 2.0;
        b = 4.0;
        bern = BernsteinPolynomial.bernstein_poly_ab(n, a, b, x);

        Console.WriteLine("");
        Console.WriteLine("     N     K     A        B        X       BPAB(N,K,A,B)(X)");
        Console.WriteLine("");

        for (int k = 0; k <= n; k++)
        {
            Console.WriteLine("  " + n.ToString().PadLeft(4)
                                   + "  " + k.ToString().PadLeft(4)
                                   + "  " + a.ToString().PadLeft(7)
                                   + "  " + b.ToString().PadLeft(7)
                                   + "  " + x.ToString().PadLeft(7)
                                   + "  " + bern[k].ToString().PadLeft(14) + "");
        }
    }

    private static void bernstein_poly_ab_approx_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNSTEIN_POLY_AB_APPROX_TEST tests BERNSTEIN_POLY_AB_APPROX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int maxdata = 20;
        int nval = 501;

        Console.WriteLine("");
        Console.WriteLine("BERNSTEIN_POLY_AB_APPROX_TEST");
        Console.WriteLine("  BERNSTEIN_POLY_AB_APPROX evaluates the Bernstein polynomial");
        Console.WriteLine("  approximant to a function F(X).");

        double a = 1.0;
        double b = 3.0;

        Console.WriteLine("");
        Console.WriteLine("     N      Max Error");
        Console.WriteLine("");

        for (int ndata = 0; ndata <= maxdata; ndata++)
        {
            //
            //  Generate data values.
            //
            double[] xdata = new double[ndata + 1];
            double[] ydata = new double[ndata + 1];
            for (int i = 0; i <= ndata; i++)
            {
                xdata[i] = ndata switch
                {
                    0 => 0.5 * (a + b),
                    _ => ((ndata - i) * a + i * b) / ndata
                };

                ydata[i] = Math.Sin(xdata[i]);
            }

            //
            //  Compare the true function and the approximant.
            //
            double[] xval = typeMethods.r8vec_linspace_new(nval, a, b);

            double error_max = 0.0;

            double[] yval = BernsteinPolynomial.bernstein_poly_ab_approx(ndata, a, b, ydata, nval, xval);

            error_max = 0.0;
            for (int i = 0; i < nval; i++)
            {
                error_max = Math.Max(error_max, Math.Abs(yval[i] - Math.Sin(xval[i])));
            }

            Console.WriteLine("  " + ndata.ToString().PadLeft(4)
                                   + "  " + error_max.ToString().PadLeft(14) + "");

        }
    }

    private static void bernstein_to_legendre_test()

        //*****************************************************************************/
        //
        //  Purpose:
        //
        //    BERNSTEIN_TO_LEGENDRE_TEST tests BERNSTEIN_TO_LEGENDRE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 March 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n = 5;

        Console.WriteLine("");
        Console.WriteLine("BERNSTEIN_TO_LEGENDRE_TEST:");
        Console.WriteLine("  BERNSTEIN_TO_LEGENDRE returns the matrix A which maps");
        Console.WriteLine("  polynomial coefficients from Bernstein to Legendre form.");

        double[] a = BernsteinPolynomial.bernstein_to_legendre(n);
        typeMethods.r8mat_print(n + 1, n + 1, a, "  A = bernstein_to_legendre(5):");

        double[] b = BernsteinPolynomial.legendre_to_bernstein(n);
        typeMethods.r8mat_print(n + 1, n + 1, b, "  B = legendre_to_bernstein(5):");

        double[] c = typeMethods.r8mat_mm_new(n + 1, n + 1, n + 1, a, b);
        double e = typeMethods.r8mat_is_identity(n + 1, c);
        Console.WriteLine("");
        Console.WriteLine("  ||A*B-I|| = " + e + "");

    }

    private static void bernstein_to_power_test()

        //*****************************************************************************/
        //
        //  Purpose:
        //
        //    BERNSTEIN_TO_POWER_TEST tests BERNSTEIN_TO_POWER.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 March 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n = 5;

        Console.WriteLine("");
        Console.WriteLine("BERNSTEIN_TO_POWER_TEST:");
        Console.WriteLine("  BERNSTEIN_TO_POWER returns the matrix A which maps");
        Console.WriteLine("  polynomial coefficients from Bernstein to Power form.");

        double[] a = BernsteinPolynomial.bernstein_to_power(n);
        typeMethods.r8mat_print(n + 1, n + 1, a, "  A = bernstein_to_power(5):");

        double[] b = BernsteinPolynomial.power_to_bernstein(n);
        typeMethods.r8mat_print(n + 1, n + 1, b, "  B = power_to_bernstein(5):");

        double[] c = typeMethods.r8mat_mm_new(n + 1, n + 1, n + 1, a, b);
        double e = typeMethods.r8mat_is_identity(n + 1, c);
        Console.WriteLine("");
        Console.WriteLine("  ||A*B-I|| = " + e + "");
    }

    private static void bernstein_vandermonde_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNSTEIN_VANDERMONDE_TEST tests BERNSTEIN_VANDERMONDE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a;
        int n;

        Console.WriteLine("");
        Console.WriteLine("BERNSTEIN_VANDERMONDE_TEST");
        Console.WriteLine("  BERNSTEIN_VANDERMONDE returns an NxN matrix whose (I,J) entry");
        Console.WriteLine("  is the value of the J-th Bernstein polynomial of degree N-1");
        Console.WriteLine("  evaluated at the I-th equally spaced point in [0,1].");

        n = 8;
        a = BernsteinPolynomial.bernstein_vandermonde(n);
        typeMethods.r8mat_print(n, n, a, "  Bernstein Vandermonde ( 8 ):");
    }
}