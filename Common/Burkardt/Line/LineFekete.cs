using System;

namespace Burkardt.LineNS
{
    public static class LineFekete
    {
        public static void line_fekete_chebyshev(int m, double a, double b, int n, double[] x,
                ref int nf, ref double[] xf, ref double[] wf)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LINE_FEKETE_CHEBYSHEV: approximate Fekete points in an interval [A,B].
            //
            //  Discussion:
            //
            //    We use the Chebyshev basis.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 April 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Len Bos, Norm Levenberg,
            //    On the calculation of approximate Fekete points: the univariate case,
            //    Electronic Transactions on Numerical Analysis, 
            //    Volume 30, pages 377-397, 2008.
            //   
            //  Parameters:
            //
            //    Input, int M, the number of basis polynomials.
            //
            //    Input, double A, B, the endpoints of the interval.
            //
            //    Input, int N, the number of sample points.
            //    M <= N.
            //
            //    Input, double X(N), the coordinates of the sample points.
            //
            //    Output, int &NF, the number of Fekete points.
            //    If the computation is successful, NF = M.
            //
            //    Output, double XF(NF), the coordinates of the Fekete points.
            //
            //    Output, double WF(NF), the weights of the Fekete points.
            //
        {
            int i;
            int j;
            double[] mom;
            double[] v;
            double[] w;

            if (n < m)
            {
                Console.WriteLine("");
                Console.WriteLine("LINE_FEKETE_CHEBYSHEV - Fatal error!");
                Console.WriteLine("  N < M.");
                return;
            }

            //
            //  Compute the Chebyshev-Vandermonde matrix.
            //
            v = VandermondeMatrix.cheby_van1(m, a, b, n, x);
            //
            //  MOM(I) = Integral ( A <= x <= B ) Tab(A,B,I;x) dx
            //
            mom = new double[m];

            mom[0] = Math.PI * (b - a) / 2.0;
            for (i = 1; i < m; i++)
            {
                mom[i] = 0.0;
            }

            //
            //  Solve the system for the weights W.
            //
            w = QRSolve.qr_solve(m, n, v, mom);
            //
            //  Extract the data associated with the nonzero weights.
            //
            nf = 0;
            for (j = 0; j < n; j++)
            {
                if (w[j] != 0.0)
                {
                    if (nf < m)
                    {
                        xf[nf] = x[j];
                        wf[nf] = w[j];
                        nf = nf + 1;
                    }
                }
            }
        }

        public static void line_fekete_legendre(int m, double a, double b, int n, double[] x,
                ref int nf, ref double[] xf, ref double[] wf)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LINE_FEKETE_LEGENDRE: approximate Fekete points in an interval [A,B].
            //
            //  Discussion:
            //
            //    We use the uniform weight and the Legendre basis:
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 April 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Len Bos, Norm Levenberg,
            //    On the calculation of approximate Fekete points: the univariate case,
            //    Electronic Transactions on Numerical Analysis, 
            //    Volume 30, pages 377-397, 2008.
            //    
            //  Parameters:
            //
            //    Input, int M, the number of basis polynomials.
            //
            //    Input, double A, B, the endpoints of the interval.
            //
            //    Input, int N, the number of sample points.
            //    M <= N.
            //
            //    Input, double X(N), the coordinates of the sample points.
            //
            //    Output, int &NF, the number of Fekete points.
            //    If the computation is successful, NF = M.
            //
            //    Output, double XF(NF), the coordinates of the Fekete points.
            //
            //    Output, double WF(NF), the weights of the Fekete points.
            //
        {
            int i;
            int j;
            double[] mom;
            double[] v;
            double[] w;

            if (n < m)
            {
                Console.WriteLine("");
                Console.WriteLine("LINE_FEKETE_LEGENDRE - Fatal error!");
                Console.WriteLine("  N < M.");
                return;
            }

            //
            //  Compute the Legendre-Vandermonde matrix.
            //
            v = VandermondeMatrix.legendre_van(m, a, b, n, x);
            //
            //  MOM(i) = integral ( A <= X <= B ) Lab(A,B,I;X) dx
            //
            mom = new double[m];
            mom[0] = b - a;
            for (i = 1; i < m; i++)
            {
                mom[i] = 0.0;
            }

            //
            //  Solve the system for the weights W.
            //
            w = QRSolve.qr_solve(m, n, v, mom);
            //
            //  Extract the data associated with the nonzero weights.
            //
            nf = 0;
            for (j = 0; j < n; j++)
            {
                if (w[j] != 0.0)
                {
                    if (nf < m)
                    {
                        xf[nf] = x[j];
                        wf[nf] = w[j];
                        nf = nf + 1;
                    }
                }
            }

        }

        public static void line_fekete_monomial(int m, double a, double b, int n, double[] x,
                ref int nf, ref double[] xf, ref double[] wf)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LINE_FEKETE_MONOMIAL: approximate Fekete points in an interval [A,B].
            //
            //  Discussion:
            //
            //    We use the uniform weight and the monomial basis:
            //
            //      P(j) = x^(j-1)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 April 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Alvise Sommariva, Marco Vianello,
            //    Computing approximate Fekete points by QR factorizations of Vandermonde 
            //    matrices,
            //    Computers and Mathematics with Applications,
            //    Volume 57, 2009, pages 1324-1336.
            //    
            //  Parameters:
            //
            //    Input, int M, the number of basis polynomials.
            //
            //    Input, double A, B, the endpoints of the interval.
            //
            //    Input, int N, the number of sample points.
            //    M <= N.
            //
            //    Input, double X(N), the coordinates of the sample points.
            //
            //    Output, int &NF, the number of Fekete points.
            //    If the computation is successful, NF = M.
            //
            //    Output, double XF(NF), the coordinates of the Fekete points.
            //
            //    Output, double WF(NF), the weights of the Fekete points.
            //
        {
            int i;
            int j;
            double[] mom;
            double[] v;
            double[] w;

            if (n < m)
            {
                Console.WriteLine("");
                Console.WriteLine("LINE_FEKETE_MONOMIAL - Fatal error!");
                Console.WriteLine("  N < M.");
                return;
            }

            //
            //  Form the moments.
            //
            mom = Monomial.line_monomial_moments(a, b, m);
            //
            //  Form the rectangular Vandermonde matrix V for the polynomial basis.
            //
            v = new double[m * n];
            for (j = 0; j < n; j++)
            {
                v[0 + j * m] = 1.0;
                for (i = 1; i < m; i++)
                {
                    v[i + j * m] = v[i - 1 + j * m] * x[j];
                }
            }

            //
            //  Solve the system for the weights W.
            //
            w = QRSolve.qr_solve(m, n, v, mom);
            //
            //  Extract the data associated with the nonzero weights.
            //
            nf = 0;
            for (j = 0; j < n; j++)
            {
                if (w[j] != 0.0)
                {
                    if (nf < m)
                    {
                        xf[nf] = x[j];
                        wf[nf] = w[j];
                        nf = nf + 1;
                    }
                }
            }
        }
    }
}