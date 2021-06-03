using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double[] r8ge_mtm(int n, double[] a, double[] b )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_MTM computes C=A'*B for R8GE matrices.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 August 2015
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
//    Input, double A[N*N], B[N*N], the factors.
//
//    Output, double C[N*N], the product.
//
        {
            double[] c = new double[n * n];

            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < n; i++)
                {
                    c[i + j * n] = 0.0;
                    for (int k = 0; k < n; k++)
                    {
                        c[i + j * n] = c[i + j * n] + a[k + i * n] * b[k + j * n];
                    }
                }
            }

            return c;
        }

        public static void r8ge_print(int m, int n, double[] a, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_PRINT prints an R8GE matrix.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
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
//    Input, double A[M*N], the R8GE matrix.
//
//    Input, string TITLE, a title.
//
        {
            r8ge_print_some(m, n, a, 1, 1, m, n, title);
        }

        public static void r8ge_print_some(int m, int n, double[] a, int ilo, int jlo, int ihi,
        int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_PRINT_SOME prints some of an R8GE matrix.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
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
//    Input, double A[M*N], the R8GE matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
        {
            int INCX = 5;

            int i2hi;
            int i2lo;
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
//
//  For each column J in the current range...
//
//  Write the header.
//
                string cout = "  Col:    ";
                for (int j = j2lo; j <= j2hi; j++)
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
                i2hi = Math.Min(ihi, m);

                for (int i = i2lo; i <= i2hi; i++)
                {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
                    cout = i.ToString().PadLeft(5) + "  ";
                    for (int j = j2lo; j <= j2hi; j++)
                    {
                        cout += a[i - 1 + (j - 1) * m].ToString().PadLeft(12) + "  ";
                    }

                    Console.WriteLine(cout);
                }
            }
        }
    }
}