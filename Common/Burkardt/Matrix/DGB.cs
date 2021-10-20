using System;

namespace Burkardt.MatrixNS
{
    public static partial class Matrix
    {
        public static double[] dgb_mxv(int m, int n, int ml, int mu, double[] a, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DGB_MXV multiplies a DGB matrix times a vector.
            //
            //  Discussion:
            //
            //    The DGB storage format is for an M by N banded matrix, with lower
            //    bandwidth ML and upper bandwidth MU.  Storage includes room for ML
            //    extra superdiagonals, which may be required to store nonzero entries
            //    generated during Gaussian elimination.
            //
            //    The original M by N matrix is "collapsed" downward, so that diagonals
            //    become rows of the storage array, while columns are preserved.  The
            //    collapsed array is logically 2*ML+MU+1 by N.
            //
            //    LINPACK and LAPACK storage of general band matrices requires
            //    an extra ML upper diagonals for possible fill in entries during
            //    Gauss elimination.  This routine does not access any entries
            //    in the fill in diagonals, because it assumes that the matrix
            //    has NOT had Gauss elimination applied to it.  If the matrix
            //    has been Gauss eliminated, then the routine DGB_MU must be
            //    used instead.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 January 1999
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Dongarra, Bunch, Cleve Moler, Stewart,
            //    LINPACK User's Guide,
            //    SIAM, 1979.
            //
            //  Parameters:
            //
            //    Input, int M, the number of rows of the matrix.
            //    M must be positive.
            //
            //    Input, int N, the number of columns of the matrix.
            //    N must be positive.
            //
            //    Input, int ML, MU, the lower and upper bandwidths.
            //    ML and MU must be nonnegative, and no greater than min(M,N)-1.
            //
            //    Input, double A[(2*ML+MU+1)*N], the DGB matrix.
            //
            //    Input, double X[N], the vector to be multiplied by A.
            //
            //    Output, double DGB_MXV[M], the product A * x.
            //
        {
            double[] b;
            int i;
            int j;
            int jhi;
            int jlo;

            b = new double[m];

            for (i = 1; i <= m; i++)
            {
                b[i - 1] = 0.0;
                jlo = Math.Max(1, i - ml);
                jhi = Math.Min(n, i + mu);
                for (j = jlo; j <= jhi; j++)
                {
                    b[i - 1] = b[i - 1] + a[i - j + ml + mu + (j - 1) * (2 * ml + mu + 1)] * x[j - 1];
                }
            }

            return b;
        }

        public static void dgb_print_some(int m, int n, int ml, int mu, double[] a, int ilo,
                int jlo, int ihi, int jhi, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DGB_PRINT_SOME prints some of a DGB matrix.
            //
            //  Discussion:
            //
            //    The DGB storage format is used for an M by N banded matrix, with lower
            //    bandwidth ML and upper bandwidth MU.  Storage includes room for ML extra
            //    superdiagonals, which may be required to store nonzero entries generated
            //    during Gaussian elimination.
            //
            //    The original M by N matrix is "collapsed" downward, so that diagonals
            //    become rows of the storage array, while columns are preserved.  The
            //    collapsed array is logically 2*ML+MU+1 by N.
            //
            //    The two dimensional array can be further reduced to a one dimensional
            //    array, stored by columns.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 April 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the number of rows of the matrix.
            //    M must be positive.
            //
            //    Input, int N, the number of columns of the matrix.
            //    N must be positive.
            //
            //    Input, int ML, MU, the lower and upper bandwidths.
            //    ML and MU must be nonnegative, and no greater than min(M,N)-1..
            //
            //    Input, double A[(2*ML+MU+1)*N], the DGB matrix.
            //
            //    Input, int ILO, JLO, IHI, JHI, designate the first row and
            //    column, and the last row and column to be printed.
            //
            //    Input, string TITLE, a title.
            //
        {
            int INCX = 5;

            int col = 2 * ml + mu + 1;
            int i;
            int i2hi;
            int i2lo;
            int j;
            int j2hi;
            int j2lo;

            Console.WriteLine("");
            Console.WriteLine(title + "");
            //
            //  Print the columns of the matrix, in strips of 5.
            //
            for (j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX)
            {
                j2hi = j2lo + INCX - 1;
                j2hi = Math.Min(j2hi, n);
                j2hi = Math.Min(j2hi, jhi);

                Console.WriteLine("");
                string cout = "  Col: ";
                for (j = j2lo; j <= j2hi; j++)
                {
                    cout += j.ToString().PadLeft(7) + "       ";
                }

                Console.WriteLine(cout);
                Console.WriteLine("  Row");
                Console.WriteLine("  ---");
                //
                //  Determine the range of the rows in this strip.
                //
                i2lo = Math.Max(ilo, 1);
                i2lo = Math.Max(i2lo, j2lo - mu - ml);

                i2hi = Math.Min(ihi, m);
                i2hi = Math.Min(i2hi, j2hi + ml);

                for (i = i2lo; i <= i2hi; i++)
                {
                    //
                    //  Print out (up to) 5 entries in row I, that lie in the current strip.
                    //
                    cout = i.ToString().PadLeft(6) + "  ";
                    for (j = j2lo; j <= j2hi; j++)
                    {
                        if (i < j - mu - ml || j + ml < i)
                        {
                            cout += "            ";
                        }
                        else
                        {
                            cout += a[i - j + ml + mu + (j - 1) * col].ToString().PadLeft(10) + "  ";
                        }
                    }

                    Console.WriteLine(cout);
                }
            }

            Console.WriteLine("");
        }
    }
}