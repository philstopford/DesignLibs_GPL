namespace Burkardt.Types;

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

        int m2 = (m - 1) * (mfat + 1) + 1;
        int n2 = (n - 1) * (nfat + 1) + 1;

        double[] xfat = new double[m2 * n2];

        for (i = 1; i <= m; i++)
        {
            int ihi = i < m ? mfat : 0;

            int j;
            for (j = 1; j <= n; j++)
            {
                int jhi = j < n ? nfat : 0;

                int ip1;
                if (i < m)
                {
                    ip1 = i + 1;
                }
                else
                {
                    ip1 = i;
                }

                int jp1;
                if (j < n)
                {
                    jp1 = j + 1;
                }
                else
                {
                    jp1 = j;
                }

                double x00 = x[i - 1 + (j - 1) * m];
                double x10 = x[ip1 - 1 + (j - 1) * m];
                double x01 = x[i - 1 + (jp1 - 1) * m];
                double x11 = x[ip1 - 1 + (jp1 - 1) * m];

                int ii;
                for (ii = 0; ii <= ihi; ii++)
                {
                    double s = ii / (double) (ihi + 1);

                    int jj;
                    for (jj = 0; jj <= jhi; jj++)
                    {
                        double t = jj / (double) (jhi + 1);

                        int iii = 1 + (i - 1) * (mfat + 1) + ii;
                        int jjj = 1 + (j - 1) * (nfat + 1) + jj;

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
        int i;

        double[] a2 = new double[m2 * n2];

        for (i = 1; i <= m2; i++)
        {
            double r = m2 switch
            {
                1 => 0.5,
                _ => (i - 1) / (double) (m2 - 1)
            };

            int i1 = 1 + (int) (r * (m - 1));
            int i2 = i1 + 1;

            if (m < i2)
            {
                i1 = m - 1;
                i2 = m;
            }

            double r1 = (i1 - 1) / (double) (m - 1);
            double r2 = (i2 - 1) / (double) (m - 1);

            int j;
            for (j = 1; j <= n2; j++)
            {
                double s = n2 switch
                {
                    1 => 0.5,
                    _ => (j - 1) / (double) (n2 - 1)
                };

                int j1 = 1 + (int) (s * (n - 1));
                int j2 = j1 + 1;

                if (n < j2)
                {
                    j1 = n - 1;
                    j2 = n;
                }

                double s1 = (j1 - 1) / (double) (n - 1);
                double s2 = (j2 - 1) / (double) (n - 1);

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