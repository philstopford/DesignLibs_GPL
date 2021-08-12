using System;
using Burkardt.IntegralNS;
using Burkardt.MatrixNS;
using Burkardt.Types;

namespace Burkardt
{
    public static class Jacobi
    {
        public static double[] jacobi1(int n, double[] a, double[] b, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    JACOBI1 carries out one step of the Jacobi iteration.
            //
            //  Discussion:
            //
            //    The linear system A*x=b is to be solved.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 January 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrix.
            //
            //    Input, double A[N,N], the matrix.
            //
            //    Input, double B[N], the right hand side.
            //
            //    Input, double X[N], the current solution estimate.
            //
            //    Output, double JACOBI1[N], the solution estimate updated by
            //    one step of the Jacobi iteration.
            //
        {
            int i;
            int j;
            double[] x_new;

            x_new = new double[n];

            for (i = 0; i < n; i++)
            {
                x_new[i] = b[i];
                for (j = 0; j < n; j++)
                {
                    if (j != i)
                    {
                        x_new[i] = x_new[i] - a[i + j * n] * x[j];
                    }
                }

                x_new[i] = x_new[i] / a[i + i * n];
            }

            return x_new;
        }

        public static void jacobi_eigenvalue(int n, double[] a, int it_max, ref double[] v,
                ref double[] d, ref int it_num, ref int rot_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    JACOBI_EIGENVALUE carries out the Jacobi eigenvalue iteration.
            //
            //  Discussion:
            //
            //    This function computes the eigenvalues and eigenvectors of a
            //    real symmetric matrix, using Rutishauser's modfications of the classical
            //    Jacobi rotation method with threshold pivoting. 
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 July 2013
            //
            //  Author:
            //
            //    C++ version by John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrix.
            //
            //    Input, double A[N*N], the matrix, which must be square, real,
            //    and symmetric.
            //
            //    Input, int IT_MAX, the maximum number of iterations.
            //
            //    Output, double V[N*N], the matrix of eigenvectors.
            //
            //    Output, double D[N], the eigenvalues, in descending order.
            //
            //    Output, int &IT_NUM, the total number of iterations.
            //
            //    Output, int &ROT_NUM, the total number of rotations.
            //
        {
            double[] bw;
            double c;
            double g;
            double gapq;
            double h;
            int i;
            int j;
            int k;
            int l;
            int m;
            int p;
            int q;
            double s;
            double t;
            double tau;
            double term;
            double termp;
            double termq;
            double theta;
            double thresh;
            double w;
            double[] zw;

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    if (i == j)
                    {
                        v[i + j * n] = 1.0;
                    }
                    else
                    {
                        v[i + j * n] = 0.0;
                    }
                }
            }

            for (i = 0; i < n; i++)
            {
                d[i] = a[i + i * n];
            }

            bw = new double[n];
            zw = new double[n];

            for (i = 0; i < n; i++)
            {
                bw[i] = d[i];
                zw[i] = 0.0;
            }

            it_num = 0;
            rot_num = 0;

            while (it_num < it_max)
            {
                it_num = it_num + 1;
                //
                //  The convergence threshold is based on the size of the elements in
                //  the strict upper triangle of the matrix.
                //
                thresh = 0.0;
                for (j = 0; j < n; j++)
                {
                    for (i = 0; i < j; i++)
                    {
                        thresh = thresh + a[i + j * n] * a[i + j * n];
                    }
                }

                thresh = Math.Sqrt(thresh) / (double) (4 * n);

                if (thresh == 0.0)
                {
                    break;
                }

                for (p = 0; p < n; p++)
                {
                    for (q = p + 1; q < n; q++)
                    {
                        gapq = 10.0 * Math.Abs(a[p + q * n]);
                        termp = gapq + Math.Abs(d[p]);
                        termq = gapq + Math.Abs(d[q]);
                        //
                        //  Annihilate tiny offdiagonal elements.
                        //
                        if (4 < it_num &&
                            termp == Math.Abs(d[p]) &&
                            termq == Math.Abs(d[q]))
                        {
                            a[p + q * n] = 0.0;
                        }
                        //
                        //  Otherwise, apply a rotation.
                        //
                        else if (thresh <= Math.Abs(a[p + q * n]))
                        {
                            h = d[q] - d[p];
                            term = Math.Abs(h) + gapq;

                            if (term == Math.Abs(h))
                            {
                                t = a[p + q * n] / h;
                            }
                            else
                            {
                                theta = 0.5 * h / a[p + q * n];
                                t = 1.0 / (Math.Abs(theta) + Math.Sqrt(1.0 + theta * theta));
                                if (theta < 0.0)
                                {
                                    t = -t;
                                }
                            }

                            c = 1.0 / Math.Sqrt(1.0 + t * t);
                            s = t * c;
                            tau = s / (1.0 + c);
                            h = t * a[p + q * n];
                            //
                            //  Accumulate corrections to diagonal elements.
                            //
                            zw[p] = zw[p] - h;
                            zw[q] = zw[q] + h;
                            d[p] = d[p] - h;
                            d[q] = d[q] + h;

                            a[p + q * n] = 0.0;
                            //
                            //  Rotate, using information from the upper triangle of A only.
                            //
                            for (j = 0; j < p; j++)
                            {
                                g = a[j + p * n];
                                h = a[j + q * n];
                                a[j + p * n] = g - s * (h + g * tau);
                                a[j + q * n] = h + s * (g - h * tau);
                            }

                            for (j = p + 1; j < q; j++)
                            {
                                g = a[p + j * n];
                                h = a[j + q * n];
                                a[p + j * n] = g - s * (h + g * tau);
                                a[j + q * n] = h + s * (g - h * tau);
                            }

                            for (j = q + 1; j < n; j++)
                            {
                                g = a[p + j * n];
                                h = a[q + j * n];
                                a[p + j * n] = g - s * (h + g * tau);
                                a[q + j * n] = h + s * (g - h * tau);
                            }

                            //
                            //  Accumulate information in the eigenvector matrix.
                            //
                            for (j = 0; j < n; j++)
                            {
                                g = v[j + p * n];
                                h = v[j + q * n];
                                v[j + p * n] = g - s * (h + g * tau);
                                v[j + q * n] = h + s * (g - h * tau);
                            }

                            rot_num = rot_num + 1;
                        }
                    }
                }

                for (i = 0; i < n; i++)
                {
                    bw[i] = bw[i] + zw[i];
                    d[i] = bw[i];
                    zw[i] = 0.0;
                }
            }

            //
            //  Restore upper triangle of input matrix.
            //
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < j; i++)
                {
                    a[i + j * n] = a[j + i * n];
                }
            }

            //
            //  Ascending sort the eigenvalues and eigenvectors.
            //
            for (k = 0; k < n - 1; k++)
            {
                m = k;
                for (l = k + 1; l < n; l++)
                {
                    if (d[l] < d[m])
                    {
                        m = l;
                    }
                }

                if (m != k)
                {
                    t = d[m];
                    d[m] = d[k];
                    d[k] = t;
                    for (i = 0; i < n; i++)
                    {
                        w = v[i + m * n];
                        v[i + m * n] = v[i + k * n];
                        v[i + k * n] = w;
                    }
                }
            }
        }

        public static void r_jacobi(int n, double a, double b, ref double[] alpha, ref double[] beta)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R_JACOBI computes recurrence coefficients for monic Jacobi polynomials.
            //
            //  Discussion:
            //
            //    This function generates the first N recurrence coefficients for monic 
            //    Jacobi polynomials with parameters A and B. 
            //
            //    These polynomials are orthogonal on [-1,1] relative to the weight
            //
            //      w(x) = (1.0-x)^A * (1.0+x)^B. 
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
            //    Original MATLAB version by Dirk Laurie, Walter Gautschi.
            //    C version by John Burkardt.
            //
            //  Reference:
            //
            //    Walter Gautschi,
            //    Orthogonal Polynomials: Computation and Approximation,
            //    Oxford, 2004,
            //    ISBN: 0-19-850672-4,
            //    LC: QA404.5 G3555.
            //
            //  Parameters:
            //
            //    Input, int N, the number of coefficients desired.
            //
            //    Input, double A, B, the parameters for the Jacobi polynomial.
            //    -1.0 < A, -1.0 < B.
            //
            //    Output, double ALPHA[N], BETA[N], the first N recurrence
            //    coefficients.
            //
        {
            int i;
            double i_r8;
            double mu;
            double nab;
            double nu;

            if (a <= -1.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R_JACOBI - Fatal error!");
                Console.WriteLine("  Illegal value of A.");
                return;
            }

            if (b <= -1.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R_JACOBI - Fatal error!");
                Console.WriteLine("  Illegal value of B.");
                return;
            }

            nu = (b - a) / (a + b + 2.0);

            mu = Math.Pow(2.0, a + b + 1.0)
                 * typeMethods.r8_gamma(a + 1.0)
                 * typeMethods.r8_gamma(b + 1.0)
                 / typeMethods.r8_gamma(a + b + 2.0);

            alpha[0] = nu;
            beta[0] = mu;

            if (n == 1)
            {
                return;
            }

            for (i = 1; i < n; i++)
            {
                i_r8 = (double) (i + 1);
                alpha[i] = (b - a) * (b + a)
                           / (2.0 * (i_r8 - 1.0) + a + b)
                           / (2.0 * i_r8 + a + b);
            }

            beta[1] = 4.0 * (a + 1.0) * (b + 1.0)
                      / (a + b + 2.0) / (a + b + 2.0)
                      / (a + b + 3.0);

            for (i = 2; i < n; i++)
            {
                i_r8 = (double) (i + 1);
                nab = 2.0 * (i_r8 - 1.0) + a + b;
                beta[i] = 4.0 * (i_r8 - 1.0 + a) * (i_r8 - 1.0 + b)
                          * (i_r8 - 1.0) * (i_r8 - 1.0 + a + b)
                          / nab / nab
                          / (nab + 1.0)
                          / (nab - 1.0);
            }
        }

        public static double monomial_quadrature_jacobi(int expon, double alpha, double beta,
                int order, double[] w, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MONOMIAL_QUADRATURE_JACOBI applies a quadrature rule to a monomial.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    22 January 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int EXPON, the exponent.
            //
            //    Input, double ALPHA, the exponent of (1-X) in the weight factor.
            //
            //    Input, double BETA, the exponent of (1+X) in the weight factor.
            //
            //    Input, int ORDER, the number of points in the rule.
            //
            //    Input, double W[ORDER], the quadrature weights.
            //
            //    Input, double X[ORDER], the quadrature points.
            //
            //    Output, double MONOMIAL_QUADRATURE_JACOBI, the quadrature error.
            //
        {
            double exact;
            int i;
            double quad;
            double quad_error;
            //
            //  Get the exact value of the integral of the unscaled monomial.
            //
            exact = Integral.jacobi_integral(expon, alpha, beta);
            //
            //  Evaluate the unweighted monomial at the quadrature points.
            //
            quad = 0.0;
            for (i = 0; i < order; i++)
            {
                quad = quad + w[i] * Math.Pow(x[i], expon);
            }

            //
            //  Absolute error for cases where exact integral is zero,
            //  Relative error otherwise.
            //
            if (exact == 0.0)
            {
                quad_error = Math.Abs(quad);
            }
            else
            {
                quad_error = Math.Abs(quad - exact) / Math.Abs(exact);
            }

            return quad_error;
        }

        public static double j_double_product_integral(int i, int j, double a, double b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    J_DOUBLE_PRODUCT_INTEGRAL: integral of J(i,x)*J(j,x)*(1-x)^a*(1+x)^b.
            //
            //  Discussion:
            //
            //    VALUE = integral ( -1 <= x <= +1 ) J(i,x)*J(j,x)*(1-x)^a*(1+x)^b dx
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 April 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int I, J, the polynomial indices.
            //
            //    Input, double A, B, the parameters.
            //    -1 < A, B.
            //
            //    Output, double VALUE, the value of the integral.
            //
        {
            double i_r8;
            double value;

            if (i != j)
            {
                value = 0.0;
            }
            else
            {
                i_r8 = (double) (i);

                value = Math.Pow(2, a + b + 1.0)
                        / (2.0 * i_r8 + a + b + 1.0)
                        * Helpers.Gamma(i_r8 + a + 1.0)
                        * Helpers.Gamma(i_r8 + b + 1.0)
                        / typeMethods.r8_factorial(i)
                        / Helpers.Gamma(i_r8 + a + b + 1.0);
            }

            return value;
        }

        public static double j_integral(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    J_INTEGRAL evaluates a monomial integral associated with J(n,a,b,x).
            //
            //  Discussion:
            //
            //    The integral:
            //
            //      integral ( -1 <= x < +1 ) x^n dx
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 April 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the exponent.
            //    0 <= N.
            //
            //    Output, double J_INTEGRAL, the value of the integral.
            //
        {
            double value;

            if ((n % 2) == 1)
            {
                value = 0.0;
            }
            else
            {
                value = 2.0 / (double) (n + 1);
            }

            return value;
        }

        public static double[] j_polynomial(int m, int n, double alpha, double beta, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    JACOBI_POLY evaluates the Jacobi polynomial J(n,a,b,x).
            //
            //  Differential equation:
            //
            //    (1-X*X) Y'' + (BETA-ALPHA-(ALPHA+BETA+2) X) Y' + N (N+ALPHA+BETA+1) Y = 0
            //
            //  Recursion:
            //
            //    P(0,ALPHA,BETA,X) = 1,
            //
            //    P(1,ALPHA,BETA,X) = ( (2+ALPHA+BETA)*X + (ALPHA-BETA) ) / 2
            //
            //    P(N,ALPHA,BETA,X)  =
            //      (
            //        (2*N+ALPHA+BETA-1)
            //        * ((ALPHA^2-BETA**2)+(2*N+ALPHA+BETA)*(2*N+ALPHA+BETA-2)*X)
            //        * P(N-1,ALPHA,BETA,X)
            //        -2*(N-1+ALPHA)*(N-1+BETA)*(2*N+ALPHA+BETA) * P(N-2,ALPHA,BETA,X)
            //      ) / 2*N*(N+ALPHA+BETA)*(2*N-2+ALPHA+BETA)
            //
            //  Restrictions:
            //
            //    -1 < ALPHA
            //    -1 < BETA
            //
            //  Norm:
            //
            //    Integral ( -1 <= X <= 1 ) ( 1 - X )^ALPHA * ( 1 + X )^BETA
            //      * P(N,ALPHA,BETA,X)^2 dX
            //    = 2^(ALPHA+BETA+1) * Gamma ( N + ALPHA + 1 ) * Gamma ( N + BETA + 1 ) /
            //      ( 2 * N + ALPHA + BETA ) * N! * Gamma ( N + ALPHA + BETA + 1 )
            //
            //  Special values:
            //
            //    P(N,ALPHA,BETA,1) = (N+ALPHA)!/(N!*ALPHA!) for integer ALPHA.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Milton Abramowitz, Irene Stegun,
            //    Handbook of Mathematical Functions,
            //    National Bureau of Standards, 1964,
            //    ISBN: 0-486-61272-4,
            //    LC: QA47.A34.
            //
            //  Parameters:
            //
            //    Input, int M, the number of evaluation points.
            //
            //    Input, int N, the highest order polynomial to compute.  Note
            //    that polynomials 0 through N will be computed.
            //
            //    Input, double ALPHA, one of the parameters defining the Jacobi
            //    polynomials, ALPHA must be greater than -1.
            //
            //    Input, double BETA, the second parameter defining the Jacobi
            //    polynomials, BETA must be greater than -1.
            //
            //    Input, double X[M], the evaluation points.
            //
            //    Output, double J_POLYNOMIAL[M*(N+1)], the values.
            //
        {
            double c1;
            double c2;
            double c3;
            double c4;
            int i;
            int j;
            double[] v;

            if (alpha <= -1.0)
            {
                Console.WriteLine("");
                Console.WriteLine("J_POLYNOMIAL - Fatal error!");
                Console.WriteLine("  Illegal input value of ALPHA = " + alpha + "");
                Console.WriteLine("  But ALPHA must be greater than -1.");
                return null;
            }

            if (beta <= -1.0)
            {
                Console.WriteLine("");
                Console.WriteLine("J_POLYNOMIAL - Fatal error!");
                Console.WriteLine("  Illegal input value of BETA = " + beta + "");
                Console.WriteLine("  But BETA must be greater than -1.");
                return null;
            }

            if (n < 0)
            {
                return null;
            }

            v = new double[m * (n + 1)];

            for (i = 0; i < m; i++)
            {
                v[i + 0 * m] = 1.0;
            }

            if (n == 0)
            {
                return v;
            }

            for (i = 0; i < m; i++)
            {
                v[i + 1 * m] = (1.0 + 0.5 * (alpha + beta)) * x[i]
                               + 0.5 * (alpha - beta);
            }

            for (i = 0; i < m; i++)
            {
                for (j = 2; j <= n; j++)
                {
                    c1 = 2.0 * (double) (j) * ((double) (j) + alpha + beta)
                         * ((double) (2 * j - 2) + alpha + beta);

                    c2 = ((double) (2 * j - 1) + alpha + beta)
                         * ((double) (2 * j) + alpha + beta)
                         * ((double) (2 * j - 2) + alpha + beta);

                    c3 = ((double) (2 * j - 1) + alpha + beta)
                         * (alpha + beta) * (alpha - beta);

                    c4 = -(double) (2) * ((double) (j - 1) + alpha)
                                       * ((double) (j - 1) + beta)
                                       * ((double) (2 * j) + alpha + beta);

                    v[i + j * m] = ((c3 + c2 * x[i]) * v[i + (j - 1) * m] + c4 * v[i + (j - 2) * m]) / c1;
                }
            }

            return v;
        }

        public static void j_polynomial_values(ref int n_data, ref int n, ref double a, ref double b, ref double x,
        ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    J_POLYNOMIAL_VALUES returns some values of the Jacobi polynomial.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      JacobiP[ n, a, b, x ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //    Stephen Wolfram,
        //    The Mathematica Book,
        //    Fourth Edition,
        //    Cambridge University Press, 1999,
        //    ISBN: 0-521-64314-7,
        //    LC: QA76.95.W65.
        //
        //  Parameters:
        //
        //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, int &N, the degree of the polynomial.
        //
        //    Output, double &A, &B, parameters of the function.
        //
        //    Output, double &X, the argument of the function.
        //
        //    Output, double &FX, the value of the function.
        //
        {
            int N_MAX = 26;

            double[] a_vec =
            {
                0.0, 0.0, 0.0, 0,
                0.0, 0.0, 1.0, 2,
                3.0, 4.0, 5.0, 0,
                0.0, 0.0, 0.0, 0,
                0.0, 0.0, 0.0, 0,
                0.0, 0.0, 0.0, 0,
                0.0, 0.0
            }
            ;

            double[] b_vec =
            {
                1.0, 1.0, 1.0, 1.0,
                1.0, 1.0, 1.0, 1.0,
                1.0, 1.0, 1.0, 2.0,
                3.0, 4.0, 5.0, 1.0,
                1.0, 1.0, 1.0, 1.0,
                1.0, 1.0, 1.0, 1.0,
                1.0, 1.0
            }
            ;

            double[] fx_vec =
            {
                0.1000000000000000E+01,
                0.2500000000000000E+00,
                -0.3750000000000000E+00,
                -0.4843750000000000E+00,
                -0.1328125000000000E+00,
                0.2753906250000000E+00,
                -0.1640625000000000E+00,
                -0.1174804687500000E+01,
                -0.2361328125000000E+01,
                -0.2616210937500000E+01,
                0.1171875000000000E+00,
                0.4218750000000000E+00,
                0.5048828125000000E+00,
                0.5097656250000000E+00,
                0.4306640625000000E+00,
                -0.6000000000000000E+01,
                0.3862000000000000E-01,
                0.8118400000000000E+00,
                0.3666000000000000E-01,
                -0.4851200000000000E+00,
                -0.3125000000000000E+00,
                0.1891200000000000E+00,
                0.4023400000000000E+00,
                0.1216000000000000E-01,
                -0.4396200000000000E+00,
                0.1000000000000000E+01
            }
            ;

            int[] n_vec =
            {
                0, 1, 2, 3,
                4, 5, 5, 5,
                5, 5, 5, 5,
                5, 5, 5, 5,
                5, 5, 5, 5,
                5, 5, 5, 5,
                5, 5
            }
            ;

            double[] x_vec =
            {
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                -1.0E+00,
                -0.8E+00,
                -0.6E+00,
                -0.4E+00,
                -0.2E+00,
                0.0E+00,
                0.2E+00,
                0.4E+00,
                0.6E+00,
                0.8E+00,
                1.0E+00
            }
            ;

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                n = 0;
                a = 0.0;
                b = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                n = n_vec[n_data - 1];
                a = a_vec[n_data - 1];
                b = b_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static double[] j_polynomial_zeros(int n, double alpha, double beta)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    J_POLYNOMIAL_ZEROS: zeros of Jacobi polynomial J(n,a,b,x).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 April 2012
            //
            //  Author:
            //
            //    John Burkardt.
            //
            //  Reference:
            //
            //    Sylvan Elhay, Jaroslav Kautsky,
            //    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
            //    Interpolatory Quadrature,
            //    ACM Transactions on Mathematical Software,
            //    Volume 13, Number 4, December 1987, pages 399-415.
            //
            //  Parameters:
            //
            //    Input, int, N, the order.
            //
            //    Input, double, ALPHA, BETA, the parameters.
            //    -1 < ALPHA, BETA.
            //
            //    Output, double J_POLYNOMIAL_ZEROS[N], the zeros.
            //
        {
            double a2b2;
            double ab;
            double abi;
            double[] bj;
            int i;
            double i_r8;
            double[] w;
            double[] x;
            double zemu;

            ab = alpha + beta;
            abi = 2.0 + ab;
            //
            //  Define the zero-th moment.
            //
            zemu = Math.Pow(2.0, ab + 1.0) * Helpers.Gamma(alpha + 1.0)
                                      * Helpers.Gamma(beta + 1.0) / Helpers.Gamma(abi);
            //
            //  Define the Jacobi matrix.
            //
            x = new double[n];
            x[0] = (beta - alpha) / abi;
            for (i = 1; i < n; i++)
            {
                x[i] = 0.0;
            }

            bj = new double[n];

            bj[0] = 4.0 * (1.0 + alpha) * (1.0 + beta)
                    / ((abi + 1.0) * abi * abi);
            for (i = 1; i < n; i++)
            {
                bj[i] = 0.0;
            }

            a2b2 = beta * beta - alpha * alpha;

            for (i = 1; i < n; i++)
            {
                i_r8 = (double) (i + 1);
                abi = 2.0 * i_r8 + ab;
                x[i] = a2b2 / ((abi - 2.0) * abi);
                abi = abi * abi;
                bj[i] = 4.0 * i_r8 * (i_r8 + alpha) * (i_r8 + beta)
                    * (i_r8 + ab) / ((abi - 1.0) * abi);
            }

            for (i = 0; i < n; i++)
            {
                bj[i] = Math.Sqrt(bj[i]);
            }

            w = new double[n];

            w[0] = Math.Sqrt(zemu);
            for (i = 1; i < n; i++)
            {
                w[i] = 0.0;
            }

            //
            //  Diagonalize the Jacobi matrix.
            //
            IMTQLX.imtqlx(n, ref x, ref bj, ref w);

            return x;
        }

        public static void j_quadrature_rule(int n, double alpha, double beta, ref double[] x,
        ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    J_QUADRATURE_RULE: Gauss-Jacobi quadrature based on J(n,a,b,x).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 April 2012
        //
        //  Author:
        //
        //    John Burkardt.
        //
        //  Reference:
        //
        //    Sylvan Elhay, Jaroslav Kautsky,
        //    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
        //    Interpolatory Quadrature,
        //    ACM Transactions on Mathematical Software,
        //    Volume 13, Number 4, December 1987, pages 399-415.
        //
        //  Parameters:
        //
        //    Input, int, N, the order.
        //
        //    Input, double, ALPHA, BETA, the parameters.
        //    -1 < ALPHA, BETA.
        //
        //    Output, double X[N], the abscissas.
        //
        //    Output, double W[N], the weights.
        //
        {
            double a2b2;
            double ab;
            double abi;
            double[] bj;
            int i;
            double i_r8;
            double zemu;

            ab = alpha + beta;
            abi = 2.0 + ab;
            //
            //  Define the zero-th moment.
            //
            zemu = Math.Pow(2.0, ab + 1.0) * Helpers.Gamma(alpha + 1.0)
                                      * Helpers.Gamma(beta + 1.0) / Helpers.Gamma(abi);
            //
            //  Define the Jacobi matrix.
            //
            x[0] = (beta - alpha) / abi;
            for (i = 1; i < n; i++)
            {
                x[i] = 0.0;
            }

            bj = new double[n];

            bj[0] = 4.0 * (1.0 + alpha) * (1.0 + beta)
                    / ((abi + 1.0) * abi * abi);
            for (i = 1; i < n; i++)
            {
                bj[i] = 0.0;
            }

            a2b2 = beta * beta - alpha * alpha;

            for (i = 1; i < n; i++)
            {
                i_r8 = (double) (i + 1);
                abi = 2.0 * i_r8 + ab;
                x[i] = a2b2 / ((abi - 2.0) * abi);
                abi = abi * abi;
                bj[i] = 4.0 * i_r8 * (i_r8 + alpha) * (i_r8 + beta)
                    * (i_r8 + ab) / ((abi - 1.0) * abi);
            }

            for (i = 0; i < n; i++)
            {
                bj[i] = Math.Sqrt(bj[i]);
            }

            w[0] = Math.Sqrt(zemu);
            for (i = 1; i < n; i++)
            {
                w[i] = 0.0;
            }

            //
            //  Diagonalize the Jacobi matrix.
            //
            IMTQLX.imtqlx(n, ref x, ref bj, ref w);

            for (i = 0; i < n; i++)
            {
                w[i] = w[i] * w[i];
            }
        }

        public static void kjacopols(double x, double a, double b, int n, ref double[] pols, int polsIndex = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KJACOPOLS evaluates Jacobi polynomials.
            //
            //  Discussion:
            //
            //    This routine evaluates Jacobi polynomials.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU GPL license.
            //
            //  Modified:
            //
            //    30 June 2014
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
            //    This C++ version by John Burkardt.
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
            //    Input, double X, the evaluation point.
            //
            //    Input, double A, B, the parameter values.
            //
            //    Input, int N, the highest degree to be evaluated.
            //
            //    Output, double POLS[N+1], the polynomial values.
            //
        {
            double alpha;
            double beta;
            int k;
            double pk;
            double pkm1;
            double pkp1;

            pkp1 = 1.0;
            pols[polsIndex + 0] = pkp1;

            if (n == 0)
            {
                return;
            }

            pk = pkp1;
            pkp1 = (a / 2.0 - b / 2.0)
                   + (1.0 + a / 2.0 + b / 2.0) * x;
            pols[polsIndex + 1] = pkp1;

            if (n == 1)
            {
                return;
            }

            for (k = 2; k <= n; k++)
            {
                pkm1 = pk;
                pk = pkp1;

                alpha = (2.0 * k + a + b - 1.0)
                        * (a * a - b * b + (2.0 * k + a + b - 2.0)
                            * (2.0 * k + a + b) * x);

                beta = 2.0 * (k + a - 1.0) * (k + b - 1.0)
                       * (2.0 * k + a + b);

                pkp1 = (alpha * pk - beta * pkm1)
                       / (2.0 * k * (k + a + b)
                          * (2.0 * k + a + b - 2.0));

                pols[polsIndex + k] = pkp1;
            }
        }

        public static void kjacopols2(double x, double a, double b, int n, ref double[] pols,
                ref double[] ders)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KJACOPOLS2 evaluates Jacobi polynomials and derivatives.
            //
            //  Discussion:
            //
            //    This routine evaluates Jacobi polynomials and  derivatives.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU GPL license.
            //
            //  Modified:
            //
            //    30 June 2014
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
            //    This C++ version by John Burkardt.
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
            //    Input, double X, the evaluation point.
            //
            //    Input, double A, B, the parameter values.
            //
            //    Input, int N, the highest degree to be evaluated.
            //
            //    Output, double POLS[N+1], the polynomial values.
            //
            //    Output, double DERS[N+1], the polynomial derivative values.
            //
        {
            double alpha1;
            double alpha2;
            double beta;
            double dk;
            double dkm1;
            double dkp1;
            double gamma;
            int k;
            double pk;
            double pkm1;
            double pkp1;

            pkp1 = 1.0;
            pols[0] = pkp1;

            dkp1 = 0.0;
            ders[0] = dkp1;

            if (n == 0)
            {
                return;
            }

            pk = pkp1;
            pkp1 = (a / 2.0 - b / 2.0)
                   + (1.0 + a / 2.0 + b / 2.0) * x;
            pols[1] = pkp1;

            dk = dkp1;
            dkp1 = (1.0 + a / 2.0 + b / 2.0);
            ders[1] = dkp1;

            if (n == 1)
            {
                return;
            }

            for (k = 2; k <= n; k++)
            {
                pkm1 = pk;
                pk = pkp1;
                dkm1 = dk;
                dk = dkp1;

                alpha1 = (2.0 * k + a + b - 1.0) * (a * a - b * b);
                alpha2 = (2.0 * k + a + b - 1.0)
                         * ((2.0 * k + a + b - 2.0)
                            * (2.0 * k + a + b));
                beta = 2.0 * (k + a - 1.0) * (k + b - 1.0)
                       * (2.0 * k + a + b);
                gamma = (2.0 * k * (k + a + b)
                         * (2.0 * k + a + b - 2.0));
                pkp1 = ((alpha1 + alpha2 * x) * pk - beta * pkm1) / gamma;
                dkp1 = ((alpha1 + alpha2 * x) * dk
                    - beta * dkm1 + alpha2 * pk) / gamma;

                pols[k] = pkp1;
                ders[k] = dkp1;
            }
        }

        public static void jacobi_recur(ref double p2, ref double dp2, ref double p1, double x, int order,
                double alpha, double beta, double[] b, double[] c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    JACOBI_RECUR evaluates a Jacobi polynomial.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    18 February 2008
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Arthur Stroud, Don Secrest,
            //    Gaussian Quadrature Formulas,
            //    Prentice Hall, 1966,
            //    LC: QA299.4G3S7.
            //
            //  Parameters:
            //
            //    Output, double *P2, the value of J(ORDER)(X).
            //
            //    Output, double *DP2, the value of J'(ORDER)(X).
            //
            //    Output, double *P1, the value of J(ORDER-1)(X).
            //
            //    Input, double X, the point at which polynomials are evaluated.
            //
            //    Input, int ORDER, the order of the polynomial to be computed.
            //
            //    Input, double ALPHA, BETA, the exponents of (1+X) and
            //    (1-X) in the quadrature rule.
            //
            //    Input, double B[ORDER], C[ORDER], the recursion coefficients.
            //
        {
            double dp0;
            double dp1;
            int i;
            double p0;

            p1 = 1.0;
            dp1 = 0.0;

            p2 = x + (alpha - beta) / (alpha + beta + 2.0);
            dp2 = 1.0;

            for (i = 2; i <= order; i++)
            {
                p0 = p1;
                dp0 = dp1;

                p1 = p2;
                dp1 = dp2;

                p2 = (x - b[i - 1]) * (p1) - c[i - 1] * p0;
                dp2 = (x - b[i - 1]) * dp1 + (p1) - c[i - 1] * dp0;
            }
        }

        public static void jacobi_root(ref double x, int order, double alpha, double beta,
                ref double dp2, ref double p1, double[] b, double[] c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    JACOBI_ROOT improves an approximate root of a Jacobi polynomial.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    18 February 2008
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Arthur Stroud, Don Secrest,
            //    Gaussian Quadrature Formulas,
            //    Prentice Hall, 1966,
            //    LC: QA299.4G3S7.
            //
            //  Parameters:
            //
            //    Input/output, double *X, the approximate root, which
            //    should be improved on output.
            //
            //    Input, int ORDER, the order of the polynomial to be computed.
            //
            //    Input, double ALPHA, BETA, the exponents of (1+X) and
            //    (1-X) in the quadrature rule.
            //
            //    Output, double *DP2, the value of J'(ORDER)(X).
            //
            //    Output, double *P1, the value of J(ORDER-1)(X).
            //
            //    Input, double B[ORDER], C[ORDER], the recursion coefficients.
            //
        {
            double d;
            double eps;
            double p2 = 0;
            int step;
            int step_max = 10;

            eps = double.Epsilon;

            for (step = 1; step <= step_max; step++)
            {
                jacobi_recur(ref p2, ref dp2, ref p1, x, order, alpha, beta, b, c);

                d = p2 / (dp2);
                x = x - d;

                if (Math.Abs(d) <= eps * (Math.Abs(x) + 1.0))
                {
                    return;
                }
            }
        }
    }
}