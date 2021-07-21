using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void r8pp_delete ( int m, int n, ref double[][] a )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8PP_DELETE frees the memory set aside by R8PP_NEW.
            //
            //  Discussion:
            //
            //    An R8PP is a pointer to pointers to R8's, and is a sort of
            //    variably-dimensioned matrix.
            //
            //    This function releases the memory associated with an array that was 
            //    created by a command like:
            //
            //      double **a;
            //      a = r8pp_new ( m, n );
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 November 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns in the array.
            //
            //    Input, double **A, the pointer to the pointers.
            //
        {
            a = null;
        }

        public static void r8pp_print(int n, double[] a, string title)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8PP_PRINT prints a R8PP matrix.
            //
            //  Discussion:
            //
            //    The R8PP storage format is appropriate for a symmetric positive
            //    definite matrix.  Only the upper triangle of the matrix is stored,
            //    by successive partial columns, in an array of length (N*(N+1))/2,
            //    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
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
            //    N must be positive.
            //
            //    Input, double A[(N*(N+1))/2], the R8PP matrix.
            //
            //    Input, string TITLE, a title.
            //
        {
            r8pp_print_some(n, a, 1, 1, n, n, title);
        }

        public static void r8pp_print_some(int n, double[] a, int ilo, int jlo, int ihi,
                int jhi, string title)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8PP_PRINT_SOME prints some of a R8PP matrix.
            //
            //  Discussion:
            //
            //    The R8PP storage format is appropriate for a symmetric positive
            //    definite matrix.  Only the upper triangle of the matrix is stored,
            //    by successive partial columns, in an array of length (N*(N+1))/2,
            //    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
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
            //    N must be positive.
            //
            //    Input, double A[(N*(N+1))/2], the R8PP matrix.
            //
            //    Input, int ILO, JLO, IHI, JHI, designate the first row and
            //    column, and the last row and column to be printed.
            //
            //    Input, string TITLE, a title.
            //
        {
            int INCX = 5;

            Console.WriteLine("");
            Console.WriteLine(title + "");
            //
            //  Print the columns of the matrix, in strips of 5.
            //
            for (int j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX)
            {
                int j2hi = j2lo + INCX - 1;
                j2hi = Math.Min(j2hi, n);
                j2hi = Math.Min(j2hi, jhi);

                Console.WriteLine("");
                string cout = "  Col: ";
                int j;
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
                int i2lo = Math.Max(ilo, 1);
                int i2hi = Math.Min(ihi, n);

                for (int i = i2lo; i <= i2hi; i++)
                {
                    cout = i.ToString().PadLeft(6) + "  ";
                    //
                    //  Print out (up to) 5 entries in row I, that lie in the current strip.
                    //
                    for (j = j2lo; j <= j2hi; j++)
                    {
                        double aij;
                        if (i <= j)
                        {
                            aij = a[i - 1 + (j * (j - 1)) / 2];
                        }
                        else
                        {
                            aij = a[j - 1 + (i * (i - 1)) / 2];
                        }

                        cout += aij.ToString().PadLeft(12) + "  ";
                    }

                    Console.WriteLine(cout);
                }
            }
        }

        
        public static double[][] r8pp_new(int m, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8PP_NEW allocates a new R8PP.
            //
            //  Discussion:
            //
            //    An R8PP is a pointer to pointers to R8's, and is a sort of
            //    variably-dimensioned matrix.
            //
            //    A declaration of the form
            //      double **a;
            //    is necesary.  Then an assignment of the form:
            //      a = r8pp_new ( m, n );
            //    allows the user to assign entries to the matrix using typical
            //    2D array notation:
            //      a[2][3] = 17;
            //      y = a[1][0];
            //    and so on.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 November 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns in the matrix.
            //
            //    Output, double **R8PP_NEW, a pointer to the pointers to the M by N array.
            //
        {
            double[][] a;
            int i;

            a = new double [m][];

            for (i = 0; i < m; i++)
            {
                a[i] = new double[n];
            }

            return a;
        }
    }
}