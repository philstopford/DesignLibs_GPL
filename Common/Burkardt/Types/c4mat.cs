using System;
using System.Numerics;

namespace Burkardt.Types
{

    public static partial class typeMethods
    {
        public static void c4mat_print(int m, int n, Complex[] a, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C4MAT_PRINT prints a C4MAT.
            //
            //  Discussion:
            //
            //    A C4MAT is an array of complex <float> values.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    11 June 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns in the matrix.
            //
            //    Input, complex <float> A[M*N], the matrix.
            //
            //    Input, string TITLE, a title.
            //
        {
            c4mat_print_some(m, n, a, 1, 1, m, n, title);
        }

        public static void c4mat_print_some(int m, int n, Complex[] a, int ilo, int jlo,
                int ihi, int jhi, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C4MAT_PRINT_SOME prints some of a C4MAT.
            //
            //  Discussion:
            //
            //    A C4MAT is an array of complex <float> values.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 April 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns in the matrix.
            //
            //    Input, complex <float> A[M*N], the matrix.
            //
            //    Input, int ILO, JLO, IHI, JHI, the first row and
            //    column, and the last row and column to be printed.
            //
            //    Input, string TITLE, a title.
            //
        {
            Complex c;
            int i;
            int i2hi;
            int i2lo;
            int inc;
            int incx = 4;
            int j;
            int j2;
            int j2hi;
            int j2lo;

            Console.WriteLine("");
            Console.WriteLine(title + "");

            if (m <= 0 || n <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("  (None)");
                return;
            }

            //
            //  Print the columns of the matrix, in strips of INCX.
            //
            for (j2lo = jlo; j2lo <= Math.Min(jhi, n); j2lo = j2lo + incx)
            {
                j2hi = j2lo + incx - 1;
                j2hi = Math.Min(j2hi, n);
                j2hi = Math.Min(j2hi, jhi);

                inc = j2hi + 1 - j2lo;

                Console.WriteLine("");
                string cout = "  Col: ";
                for (j = j2lo; j <= j2hi; j++)
                {
                    j2 = j + 1 - j2lo;
                    cout += "     " + j.ToString().PadLeft(10) + "     ";
                }

                Console.WriteLine(cout);
                Console.WriteLine("  Row");
                Console.WriteLine("  ---");
                //
                //  Determine the range of the rows in this strip.
                //
                i2lo = Math.Max(ilo, 1);
                i2hi = Math.Min(ihi, m);

                for (i = i2lo; i <= i2hi; i++)
                {
                    cout = i.ToString().PadLeft(5) + ":";
                    //
                    //  Print out (up to) INCX entries in row I, that lie in the current strip.
                    //
                    for (j2 = 1; j2 <= inc; j2++)
                    {
                        j = j2lo - 1 + j2;
                        c = a[i - 1 + (j - 1) * m];
                        cout += "  " + c.Real.ToString().PadLeft(8)
                                     + "  " + c.Imaginary.ToString().PadLeft(8);
                    }

                    Console.WriteLine(cout);
                }
            }
        }

        public static Complex[] c4mat_test(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C4MAT_TEST returns a test matrix.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 April 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrix.
            //
            //    Output, complex <float> C4MAT_TEST[N*N], the matrix.
            //
        {
            Complex I = new Complex(0.0, 1.0);

            Complex[] a = new Complex[n * n];

            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < n; i++)
                {
                    float angle = (float) (2.0 * Math.PI * (float) (i * j) / (float) (n));

                    a[i + j * n] = Complex.Exp(I * angle) / Math.Sqrt((float) (n));
                }
            }

            return a;
        }

        public static Complex[] c4mat_test_inverse(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C4MAT_TEST_INVERSE returns the inverse of a test matrix.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 April 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrix.
            //
            //    Output, complex <float> C4MAT_TEST_INVERSE[N*N], the matrix.
            //
        {
            Complex[] a = c4mat_test(n);

            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < j; i++)
                {
                    Complex t = Complex.Conjugate(a[i + j * n]);
                    a[i + j * n] = Complex.Conjugate(a[j + i * n]);
                    a[j + i * n] = t;
                }
            }

            return a;
        }
    }
}