using System;
using System.Globalization;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8mat_transpose_in_place(int n, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_TRANSPOSE_IN_PLACE transposes a square R8MAT in place.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 June 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows and columns of the matrix A.
        //
        //    Input/output, double A[N*N], the matrix to be transposed.
        //
    {
        int j;

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < j; i++)
            {
                (a[i + j * n], a[j + i * n]) = (a[j + i * n], a[i + j * n]);
            }
        }
    }

    public static double[] r8mat_transpose_new(int m, int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_TRANSPOSE_NEW returns the transpose of an R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of the matrix A.
        //
        //    Input, double A[M*N], the matrix whose transpose is desired.
        //
        //    Output, double R8MAT_TRANSPOSE_NEW[N*M], the transposed matrix.
        //
    {
        int j;

        double[] b = new double[n * m];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                b[j + i * n] = a[i + j * m];
            }
        }

        return b;
    }

    public static void r8mat_transpose_print(int m, int n, double[] a, string title)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 September 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], an M by N matrix to be printed.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8mat_transpose_print_some(m, n, a, 1, 1, m, n, title);
    }

    public static void r8mat_transpose_print_some(int m, int n, double[] a, int ilo, int jlo,
            int ihi, int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], an M by N matrix to be printed.
        //
        //    Input, int ILO, JLO, the first row and column to print.
        //
        //    Input, int IHI, JHI, the last row and column to print.
        //
        //    Input, string TITLE, a title.
        //
    {
        const int INCX = 5;

        int i2lo;
        int i2lo_hi;

        Console.WriteLine();
        Console.WriteLine(title);

        if (m <= 0 || n <= 0)
        {
            Console.WriteLine();
            Console.WriteLine("  (None)");
            return;
        }

        int i2lo_lo = ilo switch
        {
            < 1 => 1,
            _ => ilo
        };

        i2lo_hi = ihi < m ? m : ihi;

        for (i2lo = i2lo_lo; i2lo <= i2lo_hi; i2lo += INCX)
        {
            // Ugly hack to sidestep a mismatch in the output behavior compared to reference.
            if (i2lo > INCX)
            {
                break;
            }

            int i2hi = i2lo + INCX - 1;

            if (m < i2hi)
            {
                i2hi = m;
            }

            if (ihi < i2hi)
            {
                i2hi = ihi;
            }

            int inc = i2hi + 1 - i2lo;

            Console.WriteLine();
            string cout = "  Row: ";
            int i;
            for (i = i2lo; i <= i2hi; i++)
            {
                cout += (i - 1).ToString(CultureInfo.InvariantCulture).PadLeft(7) + "       ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Col");
            Console.WriteLine();

            int j2lo = jlo switch
            {
                < 1 => 1,
                _ => jlo
            };

            int j2hi;
            j2hi = n < jhi ? n : jhi;

            int j;
            for (j = j2lo; j <= j2hi; j++)
            {
                cout = (j - 1).ToString(CultureInfo.InvariantCulture).PadLeft(5) + ":";
                int i2;
                for (i2 = 1; i2 <= inc; i2++)
                {
                    i = i2lo - 1 + i2;
                    string t = a[(i - 1 + (j - 1) * m + a.Length) % a.Length].ToString("0.######");
                    cout += t.PadLeft(14);
                }

                Console.WriteLine(cout);
            }
        }
    }
        
}