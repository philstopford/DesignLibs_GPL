namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void r8mat_house_axh(int n, ref double[] a, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_HOUSE_AXH computes A*H where H is a compact Householder matrix.
            //
            //  Discussion: 							    
            //
            //    An R8MAT is a doubly dimensioned array of double precision values, which
            //    may be stored as a vector in column-major order.
            //
            //    The Householder matrix H(V) is defined by
            //
            //      H(V) = I - 2 * v * v' / ( v' * v )
            //
            //    This routine is not particularly efficient.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 July 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of A.
            //
            //    Input/output, double A[N*N], on input, the matrix to be postmultiplied.
            //    On output, A has been replaced by A*H.
            //
            //    Input, double V[N], a vector defining a Householder matrix.
            //
        {
            double[] ah;
            int i;
            int j;
            int k;
            double v_normsq;

            v_normsq = 0.0;
            for (i = 0; i < n; i++)
            {
                v_normsq = v_normsq + v[i] * v[i];
            }

            //
            //  Compute A*H' = A*H
            //
            ah = new double[n * n];

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    ah[i + j * n] = a[i + j * n];
                    for (k = 0; k < n; k++)
                    {
                        ah[i + j * n] = ah[i + j * n] - 2.0 * a[i + k * n] * v[k] * v[j] / v_normsq;
                    }
                }
            }

            //
            //  Copy A = AH;
            //
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    a[i + j * n] = ah[i + j * n];
                }
            }
        }

        public static double[] r8mat_house_axh_new(int n, double[] a, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_HOUSE_AXH_NEW computes A*H where H is a compact Householder matrix.
            //
            //  Discussion:
            //
            //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
            //    in column-major order.
            //
            //    The Householder matrix H(V) is defined by
            //
            //      H(V) = I - 2 * v * v' / ( v' * v )
            //
            //    This routine is not particularly efficient.
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
            //  Parameters:
            //
            //    Input, int N, the order of A.
            //
            //    Input, double A[N*N], the matrix to be postmultiplied.
            //
            //    Input, double V[N], a vector defining a Householder matrix.
            //
            //    Output, double R8MAT_HOUSE_AXH[N*N], the product A*H.
            //
        {
            double[] ah;
            int i;
            int j;
            int k;
            double v_normsq;

            v_normsq = 0.0;
            for (i = 0; i < n; i++)
            {
                v_normsq = v_normsq + v[i] * v[i];
            }

            //
            //  Compute A*H' = A*H
            //
            ah = new double[n * n];

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    ah[i + j * n] = a[i + j * n];
                    for (k = 0; k < n; k++)
                    {
                        ah[i + j * n] = ah[i + j * n] - 2.0 * a[i + k * n] * v[k] * v[j] / v_normsq;
                    }
                }
            }

            return ah;
        }

        public static double[] r8mat_house_form(int n, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_HOUSE_FORM constructs a Householder matrix from its compact form.
            //
            //  Discussion:
            //
            //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
            //    in column-major order.
            //
            //    H(v) = I - 2 * v * v' / ( v' * v )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 October 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrix.
            //
            //    Input, double V[N], the vector defining the Householder matrix.
            //
            //    Output, double R8MAT_HOUSE_FORM[N*N], the Householder matrix.
            //
        {
            double beta;
            double[] h;
            int i;
            int j;
            //
            //  Compute the L2 norm of V.
            //
            beta = 0.0;
            for (i = 0; i < n; i++)
            {
                beta = beta + v[i] * v[i];
            }

            //
            //  Form the matrix H.
            //
            h = r8mat_identity_new(n);

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    h[i + j * n] = h[i + j * n] - 2.0 * v[i] * v[j] / beta;
                }
            }

            return h;
        }

        public static double[] r8mat_house_hxa(int n, double[] a, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_HOUSE_HXA computes H*A where H is a compact Householder matrix.
            //
            //  Discussion:
            //
            //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
            //    in column-major order.
            //
            //    The Householder matrix H(V) is defined by
            //
            //      H(V) = I - 2 * v * v' / ( v' * v )
            //
            //    This routine is not particularly efficient.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 October 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of A.
            //
            //    Input, double A[N*N], the matrix to be premultiplied.
            //
            //    Input, double V[N], a vector defining a Householder matrix.
            //
            //    Output, double R8MAT_HOUSE_HXA[N*N], the product H*A.
            //
        {
            double[] ha;
            int i;
            int j;
            int k;
            double v_normsq;

            v_normsq = 0.0;
            for (i = 0; i < n; i++)
            {
                v_normsq = v_normsq + v[i] * v[i];
            }

            //
            //  Compute A*H' = A*H
            //
            ha = new double[n * n];

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    ha[i + j * n] = a[i + j * n];
                    for (k = 0; k < n; k++)
                    {
                        ha[i + j * n] = ha[i + j * n] - 2.0 * v[i] * v[k] * a[k + j * n] / v_normsq;
                    }
                }
            }

            return ha;
        }

        public static double[] r8mat_house_post(int n, double[] a, int row, int col)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_HOUSE_POST computes a Householder post-multiplier matrix.
            //
            //  Discussion:
            //
            //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
            //    in column-major order.
            //
            //    H(ROW,COL) has the property that the ROW-th column of
            //    A*H(ROW,COL) is zero from entry COL+1 to the end.
            //
            //    In the most common case, where a QR factorization is being computed,
            //    ROW = COL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 October 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrices.
            //
            //    Input, double A[N*N], the matrix whose Householder matrix
            //    is to be computed.
            //
            //    Input, int ROW, COL, specify the location of the
            //    entry of the matrix A which is to be preserved.  The entries in
            //    the same row, but higher column, will be zeroed out if
            //    A is postmultiplied by H.
            //
            //    Output, double R8MAT_HOUSE_POST[N*N], the Householder matrix.
            //
        {
            double[] a_row;
            double[] h;
            int j;
            double[] v;
            //
            //  Extract the ROW-th row of A.
            //
            a_row = new double[n];

            for (j = 0; j < col - 1; j++)
            {
                a_row[j] = 0.0;
            }

            for (j = col - 1; j < n; j++)
            {
                a_row[j] = a[row + j * n];
            }

            //
            //  Set up the vector V.
            //
            v = r8vec_house_column(n, a_row, col);
            //
            //  Form the matrix H(V).
            //
            h = r8mat_house_form(n, v);

            return h;
        }

        public static double[] r8mat_house_pre(int n, double[] a, int row, int col)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_HOUSE_PRE computes a Householder pre-multiplier matrix.
            //
            //  Discussion:
            //
            //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
            //    in column-major order.
            //
            //    H(ROW,COL) has the property that the COL-th column of
            //    H(ROW,COL)*A is zero from entry ROW+1 to the end.
            //
            //    In the most common case, where a QR factorization is being computed,
            //    ROW = COL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 October 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrices.
            //
            //    Input, double A[N*N], the matrix whose Householder matrix
            //    is to be computed.
            //
            //    Input, int ROW, COL, specify the location of the
            //    entry of the matrix A which is to be preserved.  The entries in
            //    the same column, but higher rows, will be zeroed out if A is
            //    premultiplied by H.
            //
            //    Output, double R8MAT_HOUSE_PRE[N*N], the Householder matrix.
            //
        {
            double[] a_col;
            double[] h;
            int i;
            double[] v;
            //
            //  Extract the COL-th column of A.
            //
            a_col = new double[n];

            for (i = 0; i < row - 1; i++)
            {
                a_col[i] = 0.0;
            }

            for (i = row - 1; i < n; i++)
            {
                a_col[i] = a[i + col * n];
            }

            //
            //  Set up the vector V.
            //
            v = r8vec_house_column(n, a_col, row);
            //
            //  Form the matrix H(V).
            //
            h = r8mat_house_form(n, v);

            return h;
        }
        
    }
}