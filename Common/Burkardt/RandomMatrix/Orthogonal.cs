﻿using Burkardt.Types;

namespace Burkardt.RandomMatrix;

public static class Orthogonal
{
    public static double[] orth_random(int n, ref typeMethods.r8NormalData data, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ORTH_RANDOM returns the ORTH_RANDOM matrix.
        //
        //  Discussion:
        //
        //    The matrix is a random orthogonal matrix.
        //
        //  Properties:
        //
        //    The inverse of A is equal to A'.
        //    A is orthogonal: A * A' = A' * A = I.
        //    Because A is orthogonal, it is normal: A' * A = A * A'.
        //    Columns and rows of A have unit Euclidean norm.
        //    Distinct pairs of columns of A are orthogonal.
        //    Distinct pairs of rows of A are orthogonal.
        //    The L2 vector norm of A*x = the L2 vector norm of x for any vector x.
        //    The L2 matrix norm of A*B = the L2 matrix norm of B for any matrix B.
        //    det ( A ) = +1 or -1.
        //    A is unimodular.
        //    All the eigenvalues of A have modulus 1.
        //    All singular values of A are 1.
        //    All entries of A are between -1 and 1.
        //
        //  Discussion:
        //
        //    Thanks to Eugene Petrov, B I Stepanov Institute of Physics,
        //    National Academy of Sciences of Belarus, for convincingly
        //    pointing out the severe deficiencies of an earlier version of
        //    this routine.
        //
        //    Essentially, the computation involves saving the Q factor of the
        //    QR factorization of a matrix whose entries are normally distributed.
        //    However, it is only necessary to generate this matrix a column at
        //    a time, since it can be shown that when it comes time to annihilate
        //    the subdiagonal elements of column K, these (transformed) elements of
        //    column K are still normally distributed random values.  Hence, there
        //    is no need to generate them at the beginning of the process and
        //    transform them K-1 times.
        //
        //    For computational efficiency, the individual Householder transformations
        //    could be saved, as recommended in the reference, instead of being
        //    accumulated into an explicit matrix format.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 July 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Pete Stewart,
        //    Efficient Generation of Random Orthogonal Matrices With an Application
        //    to Condition Estimators,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 17, Number 3, June 1980, pages 403-409.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input/output, int &SEED, a seed for the random number
        //    generator.
        //
        //    Output, double ORTH_RANDOM[N*N] the matrix.
        //
    {
        int i;
        int j;
        //
        //  Start with A = the identity matrix.
        //
        double[] a = typeMethods.r8mat_zero_new(n, n);

        for (i = 0; i < n; i++)
        {
            a[i + i * n] = 1.0;
        }

        //
        //  Now behave as though we were computing the QR factorization of
        //  some other random matrix.  Generate the N elements of the first column,
        //  compute the Householder matrix H1 that annihilates the subdiagonal elements,
        //  and set A := A * H1' = A * H.
        //
        //  On the second step, generate the lower N-1 elements of the second column,
        //  compute the Householder matrix H2 that annihilates them,
        //  and set A := A * H2' = A * H2 = H1 * H2.
        //
        //  On the N-1 step, generate the lower 2 elements of column N-1,
        //  compute the Householder matrix HN-1 that annihilates them, and
        //  and set A := A * H(N-1)' = A * H(N-1) = H1 * H2 * ... * H(N-1).
        //  This is our random orthogonal matrix.
        //
        double[] x = new double[n];

        for (j = 0; j < n - 1; j++)
        {
            //
            //  Set the vector that represents the J-th column to be annihilated.
            //
            for (i = 0; i < j; i++)
            {
                x[i] = 0.0;
            }

            for (i = j; i < n; i++)
            {
                x[i] = typeMethods.r8_normal_01(ref data, ref seed);
            }

            //
            //  Compute the vector V that defines a Householder transformation matrix
            //  H(V) that annihilates the subdiagonal elements of X.
            //
            //  The COLUMN argument here is 1-based.
            //
            double[] v = typeMethods.r8vec_house_column(n, x, j + 1);
            //
            //  Postmultiply the matrix A by H'(V) = H(V).
            //
            typeMethods.r8mat_house_axh(n, ref a, v);

        }
            
        return a;
    }
}