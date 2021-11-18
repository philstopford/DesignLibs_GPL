using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8block_delete(int l, int m, int n, ref double[][][] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BLOCK_DELETE frees memory associated with an R8BLOCK.
        //
        //  Discussion:
        //
        //    This function releases the memory associated with an array that was 
        //    created by a command like
        //      double ***a;
        //      a = r8block_new ( l, m, n );
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 March 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int L, M, N, the number of rows, columns, and layers in the array.
        //
        //    Input, double ***A, the pointer to the data.
        //
    {
        a = null;
    }

    public static double[] r8block_expand_linear(int l, int m, int n, double[] x, int lfat,
            int mfat, int nfat )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BLOCK_EXPAND_LINEAR linearly interpolates new data into a 3D block.
        //
        //  Discussion:
        //
        //    In this routine, the expansion is specified by giving the number
        //    of intermediate values to generate between each pair of original
        //    data rows and columns.
        //
        //    The interpolation is not actually linear.  It uses the functions
        //
        //      1, x, y, z, xy, xz, yz, xyz.
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
        //    Input, int L, M, N, the dimensions of the input data.
        //
        //    Input, double X[L*M*N], the original data.
        //
        //    Input, int LFAT, MFAT, NFAT, the number of data values to interpolate
        //    original data values in the first, second and third dimensions.
        //
        //    Output, double XFAT[L2*M2*N2], the fattened data, where
        //    L2 = (L-1)*(LFAT+1)+1,
        //    M2 = (M-1)*(MFAT+1)+1,
        //    N2 = (N-1)*(NFAT+1)+1.
        //
    {
        int i;

        int l2 = (l - 1) * (lfat + 1) + 1;
        int m2 = (m - 1) * (mfat + 1) + 1;
        int n2 = (n - 1) * (nfat + 1) + 1;

        double[] xfat = new double[l2 * m2 * n2];

        for (i = 1; i <= l; i++)
        {
            int ihi = i < l ? lfat : 0;

            int j;
            for (j = 1; j <= m; j++)
            {
                int jhi = j < m ? mfat : 0;

                int k;
                for (k = 1; k <= n; k++)
                {
                    int khi = k < n ? nfat : 0;

                    int ip1;
                    if (i < l)
                    {
                        ip1 = i + 1;
                    }
                    else
                    {
                        ip1 = i;
                    }

                    int jp1;
                    if (j < m)
                    {
                        jp1 = j + 1;
                    }
                    else
                    {
                        jp1 = j;
                    }

                    int kp1;
                    if (k < n)
                    {
                        kp1 = k + 1;
                    }
                    else
                    {
                        kp1 = k;
                    }

                    double x000 = x[i - 1 + (j - 1) * l + (k - 1) * l * m];
                    double x001 = x[i - 1 + (j - 1) * l + (kp1 - 1) * l * m];
                    double x100 = x[ip1 - 1 + (j - 1) * l + (k - 1) * l * m];
                    double x101 = x[ip1 - 1 + (j - 1) * l + (kp1 - 1) * l * m];
                    double x010 = x[i - 1 + (jp1 - 1) * l + (k - 1) * l * m];
                    double x011 = x[i - 1 + (jp1 - 1) * l + (kp1 - 1) * l * m];
                    double x110 = x[ip1 - 1 + (jp1 - 1) * l + (k - 1) * l * m];
                    double x111 = x[ip1 - 1 + (jp1 - 1) * l + (kp1 - 1) * l * m];

                    int ii;
                    for (ii = 0; ii <= ihi; ii++)
                    {
                        double r = ii / (double) (ihi + 1);

                        int jj;
                        for (jj = 0; jj <= jhi; jj++)
                        {
                            double s = jj / (double) (jhi + 1);

                            int kk;
                            for (kk = 0; kk <= khi; kk++)
                            {
                                double t = kk / (double) (khi + 1);

                                int iii = 1 + (i - 1) * (lfat + 1) + ii;
                                int jjj = 1 + (j - 1) * (mfat + 1) + jj;
                                int kkk = 1 + (k - 1) * (nfat + 1) + kk;

                                xfat[iii - 1 + (jjj - 1) * l2 + (kkk - 1) * l2 * m2] =
                                    x000 * (1.0 - r) * (1.0 - s) * (1.0 - t)
                                    + x001 * (1.0 - r) * (1.0 - s) * t
                                    + x010 * (1.0 - r) * s * (1.0 - t)
                                    + x011 * (1.0 - r) * s * t
                                    + x100 * r * (1.0 - s) * (1.0 - t)
                                    + x101 * r * (1.0 - s) * t
                                    + x110 * r * s * (1.0 - t)
                                    + x111 * r * s * t;
                            }
                        }
                    }
                }
            }
        }

        return xfat;
    }

    public static double[][][] r8block_new(int l, int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BLOCK_NEW allocates a new R8BLOCK.
        //
        //  Discussion:
        //
        //    A declaration of the form
        //      double ***a;
        //    is necesary.  Then an assignment of the form:
        //      a = r8block_new ( l, m, n );
        //    allows the user to assign entries to the matrix using typical
        //    3D array notation:
        //      a[2][3][4] = 17.0;
        //      y = a[1][0][3];
        //    and so on.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 March 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int L, M, N, the number of rows, columns and layers.
        //
        //    Output, double R8BLOCK_NEW[L][M][N], a new block.
        //
    {
        int i;

        double[][][] a = new double[l][][];
            
        for (i = 0; i < l; i++)
        {
            a[i] = new double[m][];
        }

        for (i = 0; i < l; i++)
        {
            int j;
            for (j = 0; j < m; j++)
            {
                a[i][j] = new double[n];
            }
        }

        return a;
    }

    public static void r8block_print(int l, int m, int n, double[] a, string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BLOCK_PRINT prints an R8BLOCK block (a 3D matrix).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int L, M, N, the dimensions of the block.
        //
        //    Input, double A[L*M*N], the matrix to be printed.
        //
        //    Input, string TITLE, a title.
        //
    {
        int k;

        Console.WriteLine("");
        Console.WriteLine(title + "");

        for (k = 1; k <= n; k++)
        {
            Console.WriteLine("");
            Console.WriteLine("  K = " + k + "");
            Console.WriteLine("");
            int jlo;
            for (jlo = 1; jlo <= m; jlo += 5)
            {
                int jhi = Math.Min(jlo + 4, m);
                Console.WriteLine("");
                string cout = "      ";
                int j;
                for (j = jlo; j <= jhi; j++)
                {
                    cout += j.ToString().PadLeft(7) + "       ";
                }

                Console.WriteLine("");
                Console.WriteLine("");
                int i;
                for (i = 1; i <= l; i++)
                {
                    cout += i.ToString().PadLeft(5) + ":";
                    for (j = jlo; j <= jhi; j++)
                    {
                        cout += "  " + a[i - 1 + (j - 1) * l + (k - 1) * l * m].ToString().PadLeft(12);
                    }

                    Console.WriteLine(cout);
                }
            }
        }
    }

    public static double[] r8block_zero_new(int l, int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BLOCK_ZEROS_NEW returns a new zeroed R8BLOCK.
        //
        //  Discussion:
        //
        //    An R8BLOCK is a triple dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 April 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int L, M, N, the number of rows and columns.
        //
        //    Output, double R8BLOCK_ZEROS_NEW[L*M*N], the new zeroed matrix.
        //
    {
        int k;

        double[] a = new double[l * m * n];

        for (k = 0; k < n; k++)
        {
            int j;
            for (j = 0; j < m; j++)
            {
                int i;
                for (i = 0; i < l; i++)
                {
                    a[i + j * l + k * l * m] = 0.0;
                }
            }
        }

        return a;
    }

}