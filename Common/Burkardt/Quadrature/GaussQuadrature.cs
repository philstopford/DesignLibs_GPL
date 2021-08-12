using System;

namespace Burkardt.Quadrature
{
    public static class GaussQuadrature
    {
        public static void gauss(int n, double[] alpha, double[] beta, ref double[] x, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GAUSS computes a Gauss quadrature rule.
        //
        //  Discussion:
        //
        //    Given a weight function W encoded by the first N recurrence coefficients 
        //    ALPHA and BETA for the associated orthogonal polynomials, the call 
        //      call gauss ( n, alpha, beta, x, w ) 
        //    generates the nodes and weights of the N-point Gauss quadrature rule 
        //    for the weight function W.
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
        //    Original MATLAB version by Walter Gautschi.
        //    C++ version by John Burkardt.
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
        //    Input, int N, the order of the desired quadrature rule.
        //
        //    Input, double ALPHA[N], BETA[N], the alpha and beta recurrence 
        //    coefficients for the othogonal polynomials associated with the
        //    weight function.
        //
        //    Output, double X[N], W[N], the nodes and  weights of the desired 
        //    quadrature rule.  The nodes are listed in increasing order.
        //
        {
            double[] a;
            int i;
            int it_max;
            int it_num = 0;
            int j;
            int rot_num = 0;
            double[] v;
            //
            //  Define the tridiagonal Jacobi matrix.
            //
            a = new double[n * n];

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    if (i == j)
                    {
                        a[i + j * n] = alpha[i];
                    }
                    else if (i == j - 1)
                    {
                        a[i + j * n] = Math.Sqrt(beta[j]);
                    }
                    else if (i - 1 == j)
                    {
                        a[i + j * n] = Math.Sqrt(beta[i]);
                    }
                    else
                    {
                        a[i + j * n] = 0.0;
                    }
                }
            }

            //
            //  Get the eigenvectors and eigenvalues.
            //
            it_max = 100;

            v = new double[n * n];

            Jacobi.jacobi_eigenvalue(n, a, it_max, ref v, ref x, ref it_num, ref rot_num);

            for (j = 0; j < n; j++)
            {
                w[j] = beta[0] * v[0 + j * n] * v[0 + j * n];
            }
        }
    }
}