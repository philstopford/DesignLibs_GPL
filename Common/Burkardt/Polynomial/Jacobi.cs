using System;
using Burkardt.MatrixNS;
using Burkardt.Types;

namespace Burkardt.PolynomialNS;

public static class Jacobi
{
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

        switch (a)
        {
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("R_JACOBI - Fatal error!");
                Console.WriteLine("  Illegal value of A.");
                return;
        }

        switch (b)
        {
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("R_JACOBI - Fatal error!");
                Console.WriteLine("  Illegal value of B.");
                return;
        }

        double nu = (b - a) / (a + b + 2.0);

        double mu = Math.Pow(2.0, a + b + 1.0)
                    * typeMethods.r8_gamma(a + 1.0)
                    * typeMethods.r8_gamma(b + 1.0)
                    / typeMethods.r8_gamma(a + b + 2.0);

        alpha[0] = nu;
        beta[0] = mu;

        switch (n)
        {
            case 1:
                return;
        }

        for (i = 1; i < n; i++)
        {
            i_r8 = i + 1;
            alpha[i] = (b - a) * (b + a)
                       / (2.0 * (i_r8 - 1.0) + a + b)
                       / (2.0 * i_r8 + a + b);
        }

        beta[1] = 4.0 * (a + 1.0) * (b + 1.0)
                  / (a + b + 2.0) / (a + b + 2.0)
                  / (a + b + 3.0);

        for (i = 2; i < n; i++)
        {
            i_r8 = i + 1;
            double nab = 2.0 * (i_r8 - 1.0) + a + b;
            beta[i] = 4.0 * (i_r8 - 1.0 + a) * (i_r8 - 1.0 + b)
                      * (i_r8 - 1.0) * (i_r8 - 1.0 + a + b)
                      / nab / nab
                      / (nab + 1.0)
                      / (nab - 1.0);
        }
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
        int i;

        switch (alpha)
        {
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("J_POLYNOMIAL - Fatal error!");
                Console.WriteLine("  Illegal input value of ALPHA = " + alpha + "");
                Console.WriteLine("  But ALPHA must be greater than -1.");
                return null;
        }

        switch (beta)
        {
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("J_POLYNOMIAL - Fatal error!");
                Console.WriteLine("  Illegal input value of BETA = " + beta + "");
                Console.WriteLine("  But BETA must be greater than -1.");
                return null;
        }

        switch (n)
        {
            case < 0:
                return null;
        }

        double[] v = new double[m * (n + 1)];

        for (i = 0; i < m; i++)
        {
            v[i + 0 * m] = 1.0;
        }

        switch (n)
        {
            case 0:
                return v;
        }

        for (i = 0; i < m; i++)
        {
            v[i + 1 * m] = (1.0 + 0.5 * (alpha + beta)) * x[i]
                           + 0.5 * (alpha - beta);
        }

        for (i = 0; i < m; i++)
        {
            int j;
            for (j = 2; j <= n; j++)
            {
                double c1 = 2.0 * j * (j + alpha + beta)
                            * (2 * j - 2 + alpha + beta);

                double c2 = (2 * j - 1 + alpha + beta)
                            * (2 * j + alpha + beta)
                            * (2 * j - 2 + alpha + beta);

                double c3 = (2 * j - 1 + alpha + beta)
                            * (alpha + beta) * (alpha - beta);

                double c4 = -(double) 2 * (j - 1 + alpha)
                                        * (j - 1 + beta)
                                        * (2 * j + alpha + beta);

                v[i + j * m] = ((c3 + c2 * x[i]) * v[i + (j - 1) * m] + c4 * v[i + (j - 2) * m]) / c1;
            }
        }

        return v;
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
        //    C++ version by John Burkardt.
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
        int k;

        double pkp1 = 1.0;
        pols[polsIndex + 0] = pkp1;

        switch (n)
        {
            case 0:
                return;
        }

        double pk = pkp1;
        pkp1 = a / 2.0 - b / 2.0
               + (1.0 + a / 2.0 + b / 2.0) * x;
        pols[polsIndex + 1] = pkp1;

        switch (n)
        {
            case 1:
                return;
        }

        for (k = 2; k <= n; k++)
        {
            double pkm1 = pk;
            pk = pkp1;

            double alpha = (2.0 * k + a + b - 1.0)
                           * (a * a - b * b + (2.0 * k + a + b - 2.0)
                               * (2.0 * k + a + b) * x);

            double beta = 2.0 * (k + a - 1.0) * (k + b - 1.0)
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
        //    C++ version by John Burkardt.
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
        int k;

        double pkp1 = 1.0;
        pols[0] = pkp1;

        double dkp1 = 0.0;
        ders[0] = dkp1;

        switch (n)
        {
            case 0:
                return;
        }

        double pk = pkp1;
        pkp1 = a / 2.0 - b / 2.0
               + (1.0 + a / 2.0 + b / 2.0) * x;
        pols[1] = pkp1;

        double dk = dkp1;
        dkp1 = 1.0 + a / 2.0 + b / 2.0;
        ders[1] = dkp1;

        switch (n)
        {
            case 1:
                return;
        }

        for (k = 2; k <= n; k++)
        {
            double pkm1 = pk;
            pk = pkp1;
            double dkm1 = dk;
            dk = dkp1;

            double alpha1 = (2.0 * k + a + b - 1.0) * (a * a - b * b);
            double alpha2 = (2.0 * k + a + b - 1.0)
                            * ((2.0 * k + a + b - 2.0)
                               * (2.0 * k + a + b));
            double beta = 2.0 * (k + a - 1.0) * (k + b - 1.0)
                          * (2.0 * k + a + b);
            double gamma = 2.0 * k * (k + a + b)
                           * (2.0 * k + a + b - 2.0);
            pkp1 = ((alpha1 + alpha2 * x) * pk - beta * pkm1) / gamma;
            dkp1 = ((alpha1 + alpha2 * x) * dk
                - beta * dkm1 + alpha2 * pk) / gamma;

            pols[k] = pkp1;
            ders[k] = dkp1;
        }
    }

    public static double[] jacobi_poly(int n, double alpha, double beta, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    JACOBI_POLY evaluates the Jacobi polynomials at X.
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
        //        * ((ALPHA^2-BETA^2)+(2*N+ALPHA+BETA)*(2*N+ALPHA+BETA-2)*X)
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
        //    20 April 2012
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
        //    Input, int N, the highest order polynomial to compute.  Note
        //    that polynomials 0 through N will be computed.
        //
        //    Input, double ALPHA, one of the parameters defining the Jacobi
        //    polynomials, ALPHA must be greater than -1.
        //
        //    Input, double BETA, the second parameter defining the Jacobi
        //    polynomials, BETA must be greater than -1.
        //
        //    Input, double X, the point at which the polynomials are to be evaluated.
        //
        //    Output, double JACOBI_POLY[N+1], the values of the first N+1 Jacobi
        //    polynomials at the point X.
        //
    {
        int i;

        switch (alpha)
        {
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("JACOBI_POLY - Fatal error!");
                Console.WriteLine("  Illegal input value of ALPHA = " + alpha + "");
                Console.WriteLine("  But ALPHA must be greater than -1.");
                return null;
        }

        switch (beta)
        {
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("JACOBI_POLY - Fatal error!");
                Console.WriteLine("  Illegal input value of BETA = " + beta + "");
                Console.WriteLine("  But BETA must be greater than -1.");
                return null;
        }

        switch (n)
        {
            case < 0:
                return null;
        }

        double[] cx = new double[n + 1];

        cx[0] = 1.0;

        switch (n)
        {
            case 0:
                return cx;
        }

        cx[1] = (1.0 + 0.5 * (alpha + beta)) * x
                + 0.5 * (alpha - beta);

        for (i = 2; i <= n; i++)
        {
            double c1 = 2.0 * i * (i + alpha + beta)
                        * (2 * i - 2 + alpha + beta);

            double c2 = (2 * i - 1 + alpha + beta)
                        * (2 * i + alpha + beta)
                        * (2 * i - 2 + alpha + beta);

            double c3 = (2 * i - 1 + alpha + beta)
                        * (alpha + beta) * (alpha - beta);

            double c4 = -(double)2 * (i - 1 + alpha)
                                   * (i - 1 + beta)
                                   * (2 * i + alpha + beta);

            cx[i] = ((c3 + c2 * x) * cx[i - 1] + c4 * cx[i - 2]) / c1;
        }

        return cx;
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
        int i;

        p1 = 1.0;
        double dp1 = 0.0;

        p2 = x + (alpha - beta) / (alpha + beta + 2.0);
        dp2 = 1.0;

        for (i = 2; i <= order; i++)
        {
            double p0 = p1;
            double dp0 = dp1;

            p1 = p2;
            dp1 = dp2;

            p2 = (x - b[i - 1]) * p1 - c[i - 1] * p0;
            dp2 = (x - b[i - 1]) * dp1 + p1 - c[i - 1] * dp0;
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
        double p2 = 0;
        int step;
        const int step_max = 10;

        double eps = typeMethods.r8_epsilon();

        for (step = 1; step <= step_max; step++)
        {
            jacobi_recur(ref p2, ref dp2, ref p1, x, order, alpha, beta, b, c);

            double d = p2 / dp2;
            x -= d;

            if (Math.Abs(d) <= eps * (Math.Abs(x) + 1.0))
            {
                return;
            }
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
        int i;

        double ab = alpha + beta;
        double abi = 2.0 + ab;
        //
        //  Define the zero-th moment.
        //
        double zemu = Math.Pow(2.0, ab + 1.0) * Helpers.Gamma(alpha + 1.0)
                                              * Helpers.Gamma(beta + 1.0) / Helpers.Gamma(abi);
        //
        //  Define the Jacobi matrix.
        //
        double[] x = new double[n];
        x[0] = (beta - alpha) / abi;
        for (i = 1; i < n; i++)
        {
            x[i] = 0.0;
        }

        double[] bj = new double[n];

        bj[0] = 4.0 * (1.0 + alpha) * (1.0 + beta)
                / ((abi + 1.0) * abi * abi);
        for (i = 1; i < n; i++)
        {
            bj[i] = 0.0;
        }

        double a2b2 = beta * beta - alpha * alpha;

        for (i = 1; i < n; i++)
        {
            double i_r8 = i + 1;
            abi = 2.0 * i_r8 + ab;
            x[i] = a2b2 / ((abi - 2.0) * abi);
            abi *= abi;
            bj[i] = 4.0 * i_r8 * (i_r8 + alpha) * (i_r8 + beta)
                * (i_r8 + ab) / ((abi - 1.0) * abi);
        }

        for (i = 0; i < n; i++)
        {
            bj[i] = Math.Sqrt(bj[i]);
        }

        double[] w = new double[n];

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
}