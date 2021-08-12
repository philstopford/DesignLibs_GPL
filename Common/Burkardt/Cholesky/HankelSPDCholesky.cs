using Burkardt.Types;

namespace Burkardt.CholeskyNS
{
    public static class HankelSPDCholesky
    {
        public static double[] hankel_spd_cholesky_lower(int n, double[] lii, double[] liim1 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    hankel_spd_cholesky_lower returns L such that L*L' is Hankel PDS.
        //
        //  Discussion:
        //
        //    In other words, H = L * L' is a symmetric positive definite matrix
        //    with the property that H is constant along antidiagonals, so that
        //
        //      H(I+J) = h(k-1), for 1 <= I, J <= N, 1 <= K <= 2*N-1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 January 2017
        //
        //  Author:
        //
        //    S Al-Homidan, M Alshahrani.
        //    C++ implementation by John Burkardt.
        //
        //  Reference:
        //
        //    S Al-Homidan, M Alshahrani,
        //    Positive Definite Hankel Matrices Using Cholesky Factorization,
        //    Computational Methods in Applied Mathematics,
        //    Volume 9, Number 3, pages 221-225, 2009.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double LII[N], values to be used in L(I,I), 
        //    for 1 <= I <= N.
        //
        //    Input, double LIIM1[N-1], values to be used in L(I+1,I) 
        //    for 1 <= I <= N-1.
        //
        //    Output, double HANKEL_spd_CHOLESKY_LOWER[N,N], the lower 
        //    Cholesky factor.
        //
        {
            double alpha;
            double beta;
            int i;
            int j;
            double[] l;
            int q;
            int r;
            int s;
            int t;

            l = typeMethods.r8mat_zeros_new(n, n);

            for (i = 0; i < n; i++)
            {
                l[i + i * n] = lii[i];
            }

            for (i = 0; i < n - 1; i++)
            {
                l[i + 1 + i * n] = liim1[i];
            }

            for (i = 2; i < n; i++)
            {
                for (j = 0; j < i - 1; j++)
                {
                    if (((i + j) % 2) == 0)
                    {
                        q = (i + j) / 2;
                        r = q;
                    }
                    else
                    {
                        q = (i + j - 1) / 2;
                        r = q + 1;
                    }

                    alpha = 0.0;
                    for (s = 0; s <= q; s++)
                    {
                        alpha = alpha + l[q + s * n] * l[r + s * n];
                    }

                    beta = 0.0;
                    for (t = 0; t < j; t++)
                    {
                        beta = beta + l[i + t * n] * l[j + t * n];
                    }

                    l[i + j * n] = (alpha - beta) / l[j + j * n];
                }
            }

            return l;
        }
    }
}