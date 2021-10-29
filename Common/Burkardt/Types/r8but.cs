using System;
using Burkardt.Uniform;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double r8but_det(int n, int mu, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8BUT_DET computes the determinant of an R8BUT matrix.
            //
            //  Discussion:
            //
            //    The R8BUT storage format is used for a banded upper triangular matrix.
            //    The matrix is assumed to be zero above the MU-th superdiagonal.
            //    The matrix is stored in an MU+1 by N array.
            //    Columns are preserved.
            //
            //    The diagonal is stored in row MU+1 of the array.
            //    The first superdiagonal in row MU, columns 2 through N.
            //    The second superdiagonal in row MU-1, columns 3 through N.
            //    The MU-th superdiagonal in row 1, columns MU+1 through N.
            //
            //  Example:
            //
            //    N = 5, MU = 2
            //
            //    A11 A12 A13   0   0
            //      0 A22 A23 A24   0
            //      0   0 A33 A34 A35
            //      0   0   0 A44 A45
            //      0   0   0   0 A55
            //                --- ---
            //                    ---
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 October 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrix.
            //
            //    Input, int MU, the upper bandwidth.
            //
            //    Input, double A[(MU+1)*N], the R8BUT matrix.
            //
            //    Output, double R8BUT_DET, the determinant of A.
            //
        {
            double det;
            int j;

            det = 1.0;
            for (j = 1; j <= n; j++)
            {
                det = det * a[(mu + 1 - 1) + (j - 1) * (mu + 1)];
            }

            return det;
        }

        public static double[] r8but_indicator(int n, int mu)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8BUT_INDICATOR sets up an R8BUT indicator matrix.
            //
            //  Discussion:
            //
            //    The R8BUT storage format is used for a banded upper triangular matrix.
            //    The matrix is assumed to be zero above the MU-th superdiagonal.
            //    The matrix is stored in an MU+1 by N array.
            //    Columns are preserved.
            //
            //    The diagonal is stored in row MU+1 of the array.
            //    The first superdiagonal in row MU, columns 2 through N.
            //    The second superdiagonal in row MU-1, columns 3 through N.
            //    The MU-th superdiagonal in row 1, columns MU+1 through N.
            //
            //  Example:
            //
            //    N = 5, MU = 2
            //
            //    A11 A12 A13   0   0
            //      0 A22 A23 A24   0
            //      0   0 A33 A34 A35
            //      0   0   0 A44 A45
            //      0   0   0   0 A55
            //                --- ---
            //                    ---
            //
            //    The indicator matrix is stored as:
            //
            //       0   0  13  24  35
            //       0  12  23  34  45
            //      11  22  33  44  55
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 October 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of columns of the matrix.
            //
            //    Input, int MU, the upper bandwidth.
            //
            //    Output, double A[(MU+1)*N], the R8BUT matrix.
            //
        {
            double[] a;
            int fac;
            int i;
            int j;

            a = new double[(mu + 1) * n];

            fac = (int) Math.Pow(10, (int) Math.Log10(n) + 1);

            for (i = 1; i <= n; i++)
            {
                for (j = i; j <= Math.Min(n, i + mu); j++)
                {
                    a[i - j + mu + 1 - 1 + (j - 1) * (mu + 1)] = (double) (fac * i + j);
                }
            }

            for (i = 1; i <= mu; i++)
            {
                for (j = 1; j <= mu + 1 - i; j++)
                {
                    a[i - 1 + (j - 1) * (mu + 1)] = 0.0;
                }
            }

            return a;
        }

        public static double[] r8but_mtv(int n, int mu, double[] a, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8BUT_MTV multiplies a vector by an R8BUT matrix.
            //
            //  Discussion:
            //
            //    The R8BUT storage format is used for a banded upper triangular matrix.
            //    The matrix is assumed to be zero above the MU-th superdiagonal.
            //    The matrix is stored in an MU+1 by N array.
            //    Columns are preserved.
            //
            //    The diagonal is stored in row MU+1 of the array.
            //    The first superdiagonal in row MU, columns 2 through N.
            //    The second superdiagonal in row MU-1, columns 3 through N.
            //    The MU-th superdiagonal in row 1, columns MU+1 through N.
            //
            //  Example:
            //
            //    N = 5, MU = 2
            //
            //    A11 A12 A13   0   0
            //      0 A22 A23 A24   0
            //      0   0 A33 A34 A35
            //      0   0   0 A44 A45
            //      0   0   0   0 A55
            //                --- ---
            //                    ---
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 October 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrix.
            //
            //    Input, int MU, the upper bandwidth.
            //
            //    Input, double A[(MU+1)*N], the R8BUT matrix.
            //
            //    Input, double X[N], the vector to be multiplied by A.
            //
            //    Output, double R8BUT_MTV(N), the product X*A.
            //
        {
            double[] b;
            int i;
            int ilo;
            int j;

            b = r8vec_zeros_new(n);

            for (i = 1; i <= n; i++)
            {
                ilo = Math.Max(1, i - mu);
                for (j = ilo; j <= i; j++)
                {
                    b[i - 1] = b[i - 1] + x[j - 1] * a[j - i + mu + (i - 1) * (mu + 1)];
                }
            }

            return b;
        }

        public static double[] r8but_mv(int n, int mu, double[] a, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8BUT_MV multiplies an R8BUT matrix times a vector.
            //
            //  Discussion:
            //
            //    The R8BUT storage format is used for a banded upper triangular matrix.
            //    The matrix is assumed to be zero above the MU-th superdiagonal.
            //    The matrix is stored in an MU+1 by N array.
            //    Columns are preserved.
            //
            //    The diagonal is stored in row MU+1 of the array.
            //    The first superdiagonal in row MU, columns 2 through N.
            //    The second superdiagonal in row MU-1, columns 3 through N.
            //    The MU-th superdiagonal in row 1, columns MU+1 through N.
            //
            //  Example:
            //
            //    N = 5, MU = 2
            //
            //    A11 A12 A13   0   0
            //      0 A22 A23 A24   0
            //      0   0 A33 A34 A35
            //      0   0   0 A44 A45
            //      0   0   0   0 A55
            //                --- ---
            //                    ---
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 October 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrix.
            //
            //    Input, int MU, the upper bandwidth.
            //
            //    Input, double A[(MU+1)*N], the R8BUT matrix.
            //
            //    Input, double X[N], the vector to be multiplied by A.
            //
            //    Output, double R8BUT_MV[N], the product A * x.
            //
        {
            double[] b;
            int i;
            int j;

            b = r8vec_zeros_new(n);

            for (i = 1; i <= n; i++)
            {
                for (j = i; j <= Math.Min(n, i + mu); j++)
                {
                    b[i - 1] = b[i - 1] + a[i - j + mu + 1 - 1 + (j - 1) * (mu + 1)] * x[j - 1];
                }
            }

            return b;
        }

        public static void r8but_print(int n, int mu, double[] a, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8BUT_PRINT prints an R8BUT matrix.
            //
            //  Discussion:
            //
            //    The R8BUT storage format is used for a banded upper triangular matrix.
            //    The matrix is assumed to be zero above the MU-th superdiagonal.
            //    The matrix is stored in an MU+1 by N array.
            //    Columns are preserved.
            //
            //    The diagonal is stored in row MU+1 of the array.
            //    The first superdiagonal in row MU, columns 2 through N.
            //    The second superdiagonal in row MU-1, columns 3 through N.
            //    The MU-th superdiagonal in row 1, columns MU+1 through N.
            //
            //  Example:
            //
            //    N = 5, MU = 2
            //
            //    A11 A12 A13   0   0
            //      0 A22 A23 A24   0
            //      0   0 A33 A34 A35
            //      0   0   0 A44 A45
            //      0   0   0   0 A55
            //                --- ---
            //                    ---
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
            //    Input, int N, the order of the matrix.
            //
            //    Input, int MU, the upper bandwidth.
            //
            //    Input, double A[(MU+1)*N], the R8BUT matrix.
            //
            //    Input, string TITLE, a title.
            //
        {
            r8but_print_some(n, mu, a, 1, 1, n, n, title);

            return;
        }

        public static void r8but_print_some(int n, int mu, double[] a, int ilo, int jlo,
                int ihi, int jhi, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8BUT_PRINT_SOME prints some of an R8BUT matrix.
            //
            //  Discussion:
            //
            //    The R8BUT storage format is used for a banded upper triangular matrix.
            //    The matrix is assumed to be zero above the MU-th superdiagonal.
            //    The matrix is stored in an MU+1 by N array.
            //    Columns are preserved.
            //
            //    The diagonal is stored in row MU+1 of the array.
            //    The first superdiagonal in row MU, columns 2 through N.
            //    The second superdiagonal in row MU-1, columns 3 through N.
            //    The MU-th superdiagonal in row 1, columns MU+1 through N.
            //
            //  Example:
            //
            //    N = 5, MU = 2
            //
            //    A11 A12 A13   0   0
            //      0 A22 A23 A24   0
            //      0   0 A33 A34 A35
            //      0   0   0 A44 A45
            //      0   0   0   0 A55
            //                --- ---
            //                    ---
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
            //    Input, int N, the order of the matrix.
            //
            //    Input, int MU, the upper bandwidth.
            //
            //    Input, double A[(MU+1)*N], the R8BUT matrix.
            //
            //    Input, int ILO, JLO, IHI, JHI, the first row and
            //    column, and the last row and column to be printed.
            //
            //    Input, string TITLE, a title.
            //
        {
            int INCX = 5;

            int i;
            int i2hi;
            int i2lo;
            int inc;
            int j;
            int j2;
            int j2hi;
            int j2lo;
            string cout = "";

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

                inc = j2hi + 1 - j2lo;

                Console.WriteLine("");
                cout = "  Col: ";

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
                i2lo = Math.Max(i2lo, j2lo);
                i2hi = Math.Min(ihi, n);
                i2hi = Math.Min(i2hi, j2hi + mu);

                for (i = i2lo; i <= i2hi; i++)
                {

                    cout = i.ToString().PadLeft(4) + "  ";
                    //
                    //  Print out (up to) 5 entries in row I, that lie in the current strip.
                    //
                    for (j2 = 1; j2 <= inc; j2++)
                    {
                        j = j2lo - 1 + j2;

                        if (i <= j && j <= i + mu)
                        {
                            cout += a[i - j + mu + 1 - 1 + (j - 1) * (mu + 1)].ToString().PadLeft(12) + "  ";
                        }
                        else
                        {
                            cout += "              ";
                        }
                    }

                    Console.WriteLine(cout);
                }
            }

            Console.WriteLine("");

        }

        public static double[] r8but_random(int n, int mu, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8BUT_RANDOM randomizes an R8BUT matrix.
            //
            //  Discussion:
            //
            //    The R8BUT storage format is used for a banded upper triangular matrix.
            //    The matrix is assumed to be zero above the MU-th superdiagonal.
            //    The matrix is stored in an MU+1 by N array.
            //    Columns are preserved.
            //
            //    The diagonal is stored in row MU+1 of the array.
            //    The first superdiagonal in row MU, columns 2 through N.
            //    The second superdiagonal in row MU-1, columns 3 through N.
            //    The MU-th superdiagonal in row 1, columns MU+1 through N.
            //
            //  Example:
            //
            //    N = 5, MU = 2
            //
            //    A11 A12 A13   0   0
            //      0 A22 A23 A24   0
            //      0   0 A33 A34 A35
            //      0   0   0 A44 A45
            //      0   0   0   0 A55
            //                --- ---
            //                    ---
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 October 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of columns of the matrix.
            //
            //    Input, int MU, the upper bandwidth.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, double R8BUT_RANDOM[(MU+1)*N], the R8BUT matrix.
            //
        {
            double[] a;
            int i;
            int j;

            a = new double[(mu + 1) * n];

            for (i = 1; i <= mu + 1; i++)
            {
                for (j = 1; j <= mu + 1 - i; j++)
                {
                    a[i - 1 + (j - 1) * (mu + 1)] = 0.0;
                }

                for (j = Math.Max(1, mu + 2 - i); j <= n; j++)
                {
                    a[i - 1 + (j - 1) * (mu + 1)] = UniformRNG.r8_uniform_01(ref seed);
                }

            }

            return a;
        }

        public static double[] r8but_sl(int n, int mu, double[] a, double[] b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8BUT_SL solves A*x=b, where A is an R8BUT matrix.
            //
            //  Discussion:
            //
            //    The R8BUT storage format is used for a banded upper triangular matrix.
            //    The matrix is assumed to be zero above the MU-th superdiagonal.
            //    The matrix is stored in an MU+1 by N array.
            //    Columns are preserved.
            //
            //    The diagonal is stored in row MU+1 of the array.
            //    The first superdiagonal in row MU, columns 2 through N.
            //    The second superdiagonal in row MU-1, columns 3 through N.
            //    The MU-th superdiagonal in row 1, columns MU+1 through N.
            //
            //  Example:
            //
            //    N = 5, MU = 2
            //
            //    A11 A12 A13   0   0
            //      0 A22 A23 A24   0
            //      0   0 A33 A34 A35
            //      0   0   0 A44 A45
            //      0   0   0   0 A55
            //                --- ---
            //                    ---
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 October 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrix.
            //
            //    Input, int MU, the upper bandwidth.
            //
            //    Input, double A[(MU+1)*N], the R8BUT matrix.
            //
            //    Input, double B[N], the right hand side.
            //
            //    Output, double R8BUT_SL[N], the solution vector.
            //
        {
            int i;
            int j;
            int jlo;
            double[] x;

            x = new double[n];

            for (i = 0; i < n; i++)
            {
                x[i] = b[i];
            }

            for (j = n; 1 <= j; j--)
            {
                x[j - 1] = x[j - 1] / a[j - j + mu + (j - 1) * (mu + 1)];
                jlo = Math.Max(1, j - mu);
                for (i = jlo; i <= j - 1; i++)
                {
                    x[i - 1] = x[i - 1] - a[i - j + mu + (j - 1) * (mu + 1)] * x[j - 1];
                }
            }

            return x;
        }

        public static double[] r8but_slt(int n, int mu, double[] a, double[] b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8BUT_SLT solves A'*x=b, where A is an R8BUT matrix.
            //
            //  Discussion:
            //
            //    The R8BUT storage format is used for a banded upper triangular matrix.
            //    The matrix is assumed to be zero above the MU-th superdiagonal.
            //    The matrix is stored in an MU+1 by N array.
            //    Columns are preserved.
            //
            //    The diagonal is stored in row MU+1 of the array.
            //    The first superdiagonal in row MU, columns 2 through N.
            //    The second superdiagonal in row MU-1, columns 3 through N.
            //    The MU-th superdiagonal in row 1, columns MU+1 through N.
            //
            //  Example:
            //
            //    N = 5, MU = 2
            //
            //    A11 A12 A13   0   0
            //      0 A22 A23 A24   0
            //      0   0 A33 A34 A35
            //      0   0   0 A44 A45
            //      0   0   0   0 A55
            //                --- ---
            //                    ---
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 October 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrix.
            //
            //    Input, int MU, the upper bandwidth.
            //
            //    Input, double A[(MU+1)*N], the R8BUT matrix.
            //
            //    Input, double B[N], the right hand side.
            //
            //    Output, double R8BUT_SLT[N], the solution vector.
            //
        {
            int i;
            int ihi;
            int j;
            double[] x;

            x = new double[n];

            for (i = 0; i < n; i++)
            {
                x[i] = b[i];
            }

            for (j = 1; j <= n; j++)
            {
                x[j - 1] = x[j - 1] / a[j - j + mu + (j - 1) * (mu + 1)];
                ihi = Math.Min(n, j + mu);
                for (i = j + 1; i <= ihi; i++)
                {
                    x[i - 1] = x[i - 1] - a[j - i + mu + (i - 1) * (mu + 1)] * x[j - 1];
                }
            }

            return x;
        }

        public static double[] r8but_to_r8ge(int n, int mu, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8BUT_TO_R8GE copies an R8BUT matrix to an R8GE matrix.
            //
            //  Discussion:
            //
            //    The R8BUT storage format is for a banded upper triangular matrix.
            //
            //    To save storage, only the diagonal and upper triangle of A is stored,
            //    in a compact diagonal format that preserves columns.
            //
            //    The diagonal is stored in row MU+1 of the array.
            //    The first superdiagonal in row MU, columns 2 through N.
            //    The second superdiagonal in row MU-1, columns 3 through N.
            //    The MU-th superdiagonal in row 1, columns MU+1 through N.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 October 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrices.
            //    N must be positive.
            //
            //    Input, int MU, the upper bandwidth of A.
            //    MU must be nonnegative, and no greater than N-1.
            //
            //    Input, double A[(MU+1)*N], the R8BUT matrix.
            //
            //    Output, double R8BUT_TO_R8GE[N*N], the R8GE matrix.
            //
        {
            double[] b;
            int i;
            int j;

            b = new double[n * n];

            for (i = 1; i <= n; i++)
            {
                for (j = 1; j <= n; j++)
                {
                    if (i <= j && j <= i + mu)
                    {
                        b[i - 1 + (j - 1) * n] = a[mu + i - j + (j - 1) * (mu + 1)];
                    }
                    else
                    {
                        b[i - 1 + (j - 1) * n] = 0.0;
                    }
                }
            }

            return b;
        }

        public static double[] r8but_zeros(int n, int mu)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8BUT_ZEROS zeros an R8BUT matrix.
            //
            //  Discussion:
            //
            //    The R8BUT storage format is used for a banded upper triangular matrix.
            //    The matrix is assumed to be zero above the MU-th superdiagonal.
            //    The matrix is stored in an MU+1 by N array.
            //    Columns are preserved.
            //
            //    The diagonal is stored in row MU+1 of the array.
            //    The first superdiagonal in row MU, columns 2 through N.
            //    The second superdiagonal in row MU-1, columns 3 through N.
            //    The MU-th superdiagonal in row 1, columns MU+1 through N.
            //
            //  Example:
            //
            //    N = 5, MU = 2
            //
            //    A11 A12 A13   0   0
            //      0 A22 A23 A24   0
            //      0   0 A33 A34 A35
            //      0   0   0 A44 A45
            //      0   0   0   0 A55
            //                --- ---
            //                    ---
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    18 August 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of columns of the matrix.
            //
            //    Input, int MU, the upper bandwidth.
            //
            //    Output, double R8BUT_ZEROS[(MU+1)*N], the R8BUT matrix.
            //
        {
            double[] a;
            int i;
            int j;

            a = new double[(mu + 1) * n];

            for (i = 0; i < mu + 1; i++)
            {
                for (j = 0; j < n; j++)
                {
                    a[i + j * (mu + 1)] = 0.0;
                }
            }

            return a;
        }

    }
}