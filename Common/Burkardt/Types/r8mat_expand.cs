namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double[] r8mat_expand_linear(int m, int n, double[] x, int mfat, int nfat)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_EXPAND_LINEAR linearly interpolates new data into an R8MAT.
            //
            //  Discussion:
            //
            //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
            //    in column-major order.
            //
            //    In this routine, the expansion is specified by giving the number
            //    of intermediate values to generate between each pair of original
            //    data rows and columns.
            //
            //    The interpolation is not actually linear.  It uses the functions
            //
            //      1, x, y, and xy.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 October 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns of input data.
            //
            //    Input, double X[M*N], the original data.
            //
            //    Input, int MFAT, NFAT, the number of data values to interpolate
            //    between each row, and each column, of original data values.
            //
            //    Output, double XFAT[M2*N2], the fattened data, where
            //    M2 = (M-1)*(MFAT+1)+1,
            //    N2 = (N-1)*(NFAT+1)+1.
            //
        {
            int i;
            int ihi;
            int ii;
            int iii;
            int ip1;
            int j;
            int jhi;
            int jj;
            int jjj;
            int jp1;
            int m2;
            int n2;
            double s;
            double t;
            double x00;
            double x01;
            double x10;
            double x11;
            double[] xfat;

            m2 = (m - 1) * (mfat + 1) + 1;
            n2 = (n - 1) * (nfat + 1) + 1;

            xfat = new double[m2 * n2];

            for (i = 1; i <= m; i++)
            {
                if (i < m)
                {
                    ihi = mfat;
                }
                else
                {
                    ihi = 0;
                }

                for (j = 1; j <= n; j++)
                {
                    if (j < n)
                    {
                        jhi = nfat;
                    }
                    else
                    {
                        jhi = 0;
                    }

                    if (i < m)
                    {
                        ip1 = i + 1;
                    }
                    else
                    {
                        ip1 = i;
                    }

                    if (j < n)
                    {
                        jp1 = j + 1;
                    }
                    else
                    {
                        jp1 = j;
                    }

                    x00 = x[i - 1 + (j - 1) * m];
                    x10 = x[ip1 - 1 + (j - 1) * m];
                    x01 = x[i - 1 + (jp1 - 1) * m];
                    x11 = x[ip1 - 1 + (jp1 - 1) * m];

                    for (ii = 0; ii <= ihi; ii++)
                    {
                        s = (double) (ii) / (double) (ihi + 1);

                        for (jj = 0; jj <= jhi; jj++)
                        {
                            t = (double) (jj) / (double) (jhi + 1);

                            iii = 1 + (i - 1) * (mfat + 1) + ii;
                            jjj = 1 + (j - 1) * (nfat + 1) + jj;

                            xfat[iii - 1 + (jjj - 1) * m2] =
                                x00
                                + s * (x10 - x00)
                                + t * (x01 - x00)
                                + s * t * (x11 - x10 - x01 + x00);
                        }
                    }
                }
            }

            return xfat;
        }

        public static double[] r8mat_expand_linear2(int m, int n, double[] a, int m2, int n2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_EXPAND_LINEAR2 expands an R8MAT by linear interpolation.
            //
            //  Discussion:
            //
            //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
            //    in column-major order.
            //
            //    In this version of the routine, the expansion is indicated
            //    by specifying the dimensions of the expanded array.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 October 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns in A.
            //
            //    Input, double A(M,N), a "small" M by N array.
            //
            //    Input, int M2, N2, the number of rows and columns in A2.
            //
            //    Output, double R8MAT_EXPAND_LINEAR2[M2*N2], the expanded array,
            //    which contains an interpolated version of the data in A.
            //
        {
            double[] a2;
            int i;
            int i1;
            int i2;
            int j;
            int j1;
            int j2;
            double r;
            double r1;
            double r2;
            double s;
            double s1;
            double s2;

            a2 = new double[m2 * n2];

            for (i = 1; i <= m2; i++)
            {
                if (m2 == 1)
                {
                    r = 0.5;
                }
                else
                {
                    r = (double) (i - 1) / (double) (m2 - 1);
                }

                i1 = 1 + (int) (r * (double) (m - 1));
                i2 = i1 + 1;

                if (m < i2)
                {
                    i1 = m - 1;
                    i2 = m;
                }

                r1 = (double) (i1 - 1) / (double) (m - 1);
                r2 = (double) (i2 - 1) / (double) (m - 1);

                for (j = 1; j <= n2; j++)
                {
                    if (n2 == 1)
                    {
                        s = 0.5;
                    }
                    else
                    {
                        s = (double) (j - 1) / (double) (n2 - 1);
                    }

                    j1 = 1 + (int) (s * (double) (n - 1));
                    j2 = j1 + 1;

                    if (n < j2)
                    {
                        j1 = n - 1;
                        j2 = n;
                    }

                    s1 = (double) (j1 - 1) / (double) (n - 1);
                    s2 = (double) (j2 - 1) / (double) (n - 1);

                    a2[i - 1 + (j - 1) * m2] =
                        ((r2 - r) * (s2 - s) * a[i1 - 1 + (j1 - 1) * m]
                         + (r - r1) * (s2 - s) * a[i2 - 1 + (j1 - 1) * m]
                         + (r2 - r) * (s - s1) * a[i1 - 1 + (j2 - 1) * m]
                         + (r - r1) * (s - s1) * a[i2 - 1 + (j2 - 1) * m])
                        / ((r2 - r1) * (s2 - s1));
                }
            }

            return a2;
        }
        
    }
}