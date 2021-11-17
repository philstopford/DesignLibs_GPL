namespace Burkardt.MatrixNS;

public static class Fibonacci
{
    public static double[] fibonacci2(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FIBONACCI2 returns the FIBONACCI2 matrix.
        //
        //  Example:
        //
        //    N = 5
        //
        //    0 1 0 0 0
        //    1 1 0 0 0
        //    0 1 1 0 0
        //    0 0 1 1 0
        //    0 0 0 1 1
        //
        //  Properties:
        //
        //    A is generally not symmetric: A' /= A.
        //
        //    A is tridiagonal.
        //
        //    Because A is tridiagonal, it has property A (bipartite).
        //
        //    A is banded, with bandwidth 3.
        //
        //    A is integral, therefore det ( A ) is integral, and 
        //    det ( A ) * inverse ( A ) is integral.
        //
        //    A is a zero/one matrix.
        //
        //    If N = 1 then
        //      det ( A ) = 0
        //    else
        //      det ( A ) = -1
        //
        //    If 1 < N, then A is unimodular.
        //
        //    When applied to a Fibonacci1 matrix B, the Fibonacci2 matrix
        //    A produces the "next" Fibonacci1 matrix C = A*B.
        //
        //    Let PHI be the golden ratio (1+sqrt(5))/2.
        //
        //    For 2 <= N, the eigenvalues and eigenvectors are:
        //
        //    LAMBDA(1)     = PHI,     vector = (1,PHI,PHI^2,...PHI^(N-1));
        //    LAMBDA(2:N-1) = 1        vector = (0,0,0,...,0,1);
        //    LAMBDA(N)     = 1 - PHI. vector = ((-PHI)^(N-1),(-PHI)^(N-2),...,1)
        //
        //    Note that there is only one eigenvector corresponding to 1.
        //    Hence, for 3 < N, the matrix is defective.  This fact means, 
        //    for instance, that the convergence of the eigenvector in the power 
        //    method will be very slow.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 May 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, double FIBONACCI2[N*N], the matrix.
        //
    {
        double[] a;
        int i;
        int j;

        a = new double[n * n];

        for (i = 1; i <= n; i++)
        {
            for (j = 1; j <= n; j++)
            {
                switch (i)
                {
                    case 1 when j == 2:
                        a[i - 1 + (j - 1) * n] = 1.0;
                        break;
                    case 1:
                        a[i - 1 + (j - 1) * n] = 0.0;
                        break;
                    default:
                    {
                        if (j == i - 1 || j == i)
                        {
                            a[i - 1 + (j - 1) * n] = 1.0;
                        }
                        else
                        {
                            a[i - 1 + (j - 1) * n] = 0.0;
                        }

                        break;
                    }
                }
            }
        }

        return a;
    }
}