using System;
using Burkardt.MatrixNS;
using Burkardt.Types;

namespace Burkardt
{
    public static class SingleValueDecomposition
    {
        public static void svd_truncated_u(int m, int n, double[] a, ref double[] un, ref double[] sn,
                ref double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SVD_TRUNCATED_U gets the truncated SVD when N <= M
            //
            //  Discussion:
            //
            //    A(mxn) = U(mxm)  * S(mxn)  * V(nxn)'
            //           = Un(mxn) * Sn(nxn) * V(nxn)'
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 September 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns in the matrix A.
            //
            //    Input, double A[M*N], the matrix whose singular value
            //    decomposition we are investigating.
            //
            //    Output, double UN[M*N], SN[N*N], V[N*N], the factors
            //    that form the singular value decomposition of A.
            //
        {
            double[] a_copy;
            double[] e;
            int i;
            int info;
            int j;
            int lda;
            int ldu;
            int ldv;
            int job;
            double[] sdiag;
            double[] work;
            //
            //  The correct size of E and SDIAG is min ( m+1, n).
            //
            a_copy = new double[m * n];
            e = new double[m + n];
            sdiag = new double[m + n];
            work = new double[m];
            //
            //  Compute the eigenvalues and eigenvectors.
            //
            job = 21;
            lda = m;
            ldu = m;
            ldv = n;
            //
            //  The input matrix is destroyed by the routine.  Since we need to keep
            //  it around, we only pass a copy to the routine.
            //
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    a_copy[i + j * m] = a[i + j * m];
                }
            }

            info = DSVDC.dsvdc(ref a_copy, lda, m, n, ref sdiag, ref e, ref un, ldu, ref v, ldv, work, job);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("SVD_TRUNCATED_U - Failure!");
                Console.WriteLine("  The SVD could not be calculated.");
                Console.WriteLine("  LINPACK routine DSVDC returned a nonzero");
                Console.WriteLine("  value of the error flag, INFO = " + info + "");
                return;
            }

            //
            //  Make the NxN matrix S from the diagonal values in SDIAG.
            //
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    if (i == j)
                    {
                        sn[i + j * n] = sdiag[i];
                    }
                    else
                    {
                        sn[i + j * n] = 0.0;
                    }
                }
            }
        }

        public static void svd_truncated_v ( int m, int n, double[] a, ref double[] u, ref double[] sm, 
                ref double[] vm )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SVD_TRUNCATED_V gets the truncated SVD when M <= N.
            //
            //  Discussion:
            //
            //    A(mxn) = U(mxm) * S(mxn)  * V(nxn)'
            //           = U(mxm) * Sm(mxm) * Vm(mxn)'
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 March 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns in the matrix A.
            //
            //    Input, double A[M*N], the matrix whose singular value
            //    decomposition we are investigating.
            //
            //    Output, double U[M*M], SM[M*M], VM[M*N], the factors
            //    that form the singular value decomposition of A.
            //
        {
            double[] a2;
            int m2;
            int n2;
            //
            //  Transpose the matrix!
            //
            a2 = typeMethods.r8mat_transpose_new ( m, n, a );
            m2 = n;
            n2 = m;

            svd_truncated_u ( m2, n2, a2, ref vm, ref sm, ref u );

        }

        public static void singular_vectors(int m, int n, int basis_num, ref double[] a, ref double[] sval)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SINGULAR_VECTORS computes the desired singular values.
            //
            //  Discussion:
            //
            //    The LINPACK SVD routine DSVDC is used to compute the singular
            //    value decomposition:
            //
            //      A = U * S * V'
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 May 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the number of spatial dimensions.
            //
            //    Input, int N, the number of data points.
            //
            //    Input, int BASIS_NUM, the number of basis vectors to be extracted.
            //
            //    Input/output, double A[M*N]; on input, the matrix whose 
            //    singular values are to be computed.  On output, A(M,1:BASIS_NUM)
            //    contains the first BASIS_NUM left singular vectors.
            //
            //    Output, double SVAL[BASIS_NUM], the first BASIS_NUM
            //    singular values.
            //
        {
            double[] e;
            int i;
            int info;
            int lda;
            int ldu;
            int ldv;
            int job;
            double[] s;
            double[] u;
            double[] v;
            double[] work;

            Console.WriteLine("");
            Console.WriteLine("SINGULAR_VECTORS");
            Console.WriteLine("  For an MxN matrix A in general storage,");
            Console.WriteLine("  The LINPACK routine DSVDC computes the");
            Console.WriteLine("  singular value decomposition:");
            Console.WriteLine("");
            Console.WriteLine("    A = U * S * V'");
            Console.WriteLine("");
            //
            //  Compute the eigenvalues and eigenvectors.
            //
            s = new double[Math.Min(m + 1, n)];
            e = new double[n];
            lda = m;
            u = a;
            ldu = m;
            v = null;
            ldv = n;
            work = new double[m];
            job = 20;

            info = DSVDC.dsvdc(ref a, lda, m, n, ref s, ref e, ref u, ldu, ref v, ldv, work, job);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("SINGULAR_VECTORS - Warning:");
                Console.WriteLine("  DSVDC returned nonzero INFO = " + info + "");
                ;
                return;
            }

            for (i = 0; i < basis_num; i++)
            {
                sval[i] = s[i];
            }

            Console.WriteLine("");
            Console.WriteLine("  The leading singular values:");
            Console.WriteLine("");

            for (i = 0; i < basis_num; i++)
            {
                Console.WriteLine("  "
                                  + (i + 1).ToString().PadLeft(4) + "  "
                                  + sval[i].ToString().PadLeft(16) + "");
            }

        }

        public static void svd_product_test(int m, int n, double[] a, double[] u,
                double[] s, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SVD_PRODUCT_TEST tests that A = U * S * V'.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 September 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns in the matrix A.
            //
            //    Input, double A[M*N], the matrix whose singular value
            //    decomposition we are investigating.
            //
            //    Input, double U[M*M], S[M*N], V[N*N], the factors
            //    that form the singular value decomposition of A.
            //
        {
            double a_norm;
            double dif_norm;
            int i;
            int j;
            int k;
            double[] svt;
            double[] usvt;

            a_norm = typeMethods.r8mat_norm_fro(m, n, a);

            svt = new double[m * n];
            for (i = 0; i < m; i++)
            {
                for (j = 0; j < n; j++)
                {
                    svt[i + j * m] = 0.0;
                    for (k = 0; k < n; k++)
                    {
                        svt[i + j * m] = svt[i + j * m] + s[i + k * m] * v[j + k * n];
                    }
                }
            }

            usvt = new double[m * n];

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < n; j++)
                {
                    usvt[i + j * m] = 0.0;
                    for (k = 0; k < m; k++)
                    {
                        usvt[i + j * m] = usvt[i + j * m] + u[i + k * m] * svt[k + j * m];
                    }
                }
            }

            typeMethods.r8mat_print(m, n, usvt, "  The product U * S * V':");

            dif_norm = typeMethods.r8mat_dif_fro(m, n, a, usvt);

            Console.WriteLine("");
            Console.WriteLine("  Frobenius Norm of A, A_NORM = " + a_norm + "");
            Console.WriteLine("");
            Console.WriteLine("  ABSOLUTE ERROR for A = U*S*V'");
            Console.WriteLine("  Frobenius norm of difference A-U*S*V' = " + dif_norm + "");
            Console.WriteLine("");
            Console.WriteLine("  RELATIVE ERROR for A = U*S*V':");
            Console.WriteLine("  Ratio of DIF_NORM / A_NORM = " + dif_norm / a_norm + "");

        }

        public static void r8mat_svd_linpack(int m, int n, double[] a, ref double[] u, ref double[] s,
                ref double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_SVD_LINPACK gets the SVD of a matrix using a call to LINPACK.
            //
            //  Discussion:
            //
            //    The singular value decomposition of a real MxN matrix A has the form:
            //
            //      A = U * S * V'
            //
            //    where 
            //
            //      U is MxM orthogonal,
            //      S is MxN, and entirely zero except for the diagonal;
            //      V is NxN orthogonal.
            //
            //    Moreover, the nonzero entries of S are positive, and appear
            //    in order, from largest magnitude to smallest.
            //
            //    This routine calls the LINPACK routine DSVDC to compute the
            //    factorization.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 September 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns in the matrix A.
            //
            //    Input, double A[M*N], the matrix whose singular value
            //    decomposition we are investigating.
            //
            //    Output, double U[M*M], S[M*N], V[N*N], the factors
            //    that form the singular value decomposition of A.
            //
        {
            double[] a_copy;
            double[] e;
            int i;
            int info;
            int j;
            int lda;
            int ldu;
            int ldv;
            int job;
            double[] sdiag;
            double[] work;
            //
            //  The correct size of E and SDIAG is min ( m+1, n).
            //
            a_copy = new double[m * n];
            e = new double[m + n];
            sdiag = new double[m + n];
            work = new double[m];
            //
            //  Compute the eigenvalues and eigenvectors.
            //
            job = 11;
            lda = m;
            ldu = m;
            ldv = n;
            //
            //  The input matrix is destroyed by the routine.  Since we need to keep
            //  it around, we only pass a copy to the routine.
            //
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    a_copy[i + j * m] = a[i + j * m];
                }
            }

            info = DSVDC.dsvdc(ref a_copy, lda, m, n, ref sdiag, ref e, ref u, ldu, ref v, ldv, work, job);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8MAT_SVD_LINPACK - Failure!");
                Console.WriteLine("  The SVD could not be calculated.");
                Console.WriteLine("  LINPACK routine DSVDC returned a nonzero");
                Console.WriteLine("  value of the error flag, INFO = " + info + "");
                return;
            }

            //
            //  Make the MxN matrix S from the diagonal values in SDIAG.
            //
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    if (i == j)
                    {
                        s[i + j * m] = sdiag[i];
                    }
                    else
                    {
                        s[i + j * m] = 0.0;
                    }
                }
            }
            //
            //  Note that we do NOT need to transpose the V that comes out of LINPACK!
            //

            return;
        }
    }
}