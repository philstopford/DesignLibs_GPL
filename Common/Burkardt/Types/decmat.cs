using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void decmat_det(int n, int[] atop, int[] abot, int dec_digit,
                ref int dtop, ref int dbot)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DECMAT_DET finds the determinant of an N by N matrix of decimal entries.
            //
            //  Discussion:
            //
            //    The brute force method is used.  The routine should only be used for
            //    small matrices, since this calculation requires the summation of N!
            //    products of N numbers.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 October 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of rows and columns of A.
            //
            //    Input, int ATOP[N*N], ABOT[N*N], the decimal
            //    representation of the matrix.
            //
            //    Input, int DEC_DIGIT, the number of decimal digits.
            //
            //    Output, int &DTOP, &DBOT, the decimal determinant of the matrix.
            //
        {
            bool even = false;
            int i;
            int[] iarray;
            int ibot;
            int ibot1;
            int ibot2;
            int itop;
            int itop1;
            int itop2;
            bool more = false;

            dtop = 0;
            dbot = 1;

            iarray = new int[n];
            //
            //  Compute the next permutation.
            //
            more = false;

            for (;;)
            {
                Permutation.perm0_next(n, iarray, ref more, ref even);
                //
                //  The sign of this term depends on the sign of the permutation.
                //
                if (even)
                {
                    itop = 1;
                }
                else
                {
                    itop = -1;
                }

                ibot = 0;
                //
                //  Choose one item from each row, as specified by the permutation,
                //  and multiply them.
                //
                for (i = 0; i < n; i++)
                {
                    itop1 = itop;
                    ibot1 = ibot;
                    itop2 = atop[i + iarray[i] * n];
                    ibot2 = abot[i + iarray[i] * n];

                    dec_mul(itop1, ibot1, itop2, ibot2, dec_digit, ref itop, ref ibot);

                }

                //
                //  Add this term to the total.
                //
                dec_add(itop, ibot, dtop, dbot, dec_digit, ref dtop, ref dbot);

                if (!more)
                {
                    break;
                }
            }
        }

        public static void decmat_print(int m, int n, int[] a, int[] b, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DECMAT_PRINT prints out decimal vectors and matrices.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns in the matrix.
            //
            //    Input, int A[M*N], B[M*N], the decimal matrix.
            //
            //    Input, string TITLE, a title.
            //
        {
            int exponent;
            int i;
            int j;
            int jmax;
            int jmin;
            int kmax;
            int mantissa;
            int ncolum = 80;
            int npline;
            string s;
            //
            //  Figure out how wide we must make each column.
            //
            kmax = 0;

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < n; j++)
                {
                    kmax = Math.Max(kmax, dec_width(a[i + j * m], b[i + j * m]));
                }
            }

            kmax = kmax + 2;
            npline = ncolum / kmax;
            //
            //  Now do the printing.
            //
            for (jmin = 0; jmin < n; jmin = jmin + npline)
            {
                jmax = Math.Min(jmin + npline - 1, n - 1);

                if (0 < title.Length)
                {
                    Console.WriteLine("");
                    Console.WriteLine(title + "");
                }

                if (0 < jmin || jmax < n - 1)
                {
                    Console.WriteLine("Columns " + jmin + " to " + jmax + "");
                    Console.WriteLine("");
                }

                for (i = 0; i < m; i++)
                {
                    string cout = "";
                    for (j = jmin; j <= jmax; j++)
                    {
                        mantissa = a[i + j * m];
                        exponent = b[i + j * m];
                        s = dec_to_s(mantissa, exponent);
                        cout += s.PadLeft(kmax) + "  ";
                    }

                    Console.WriteLine(cout);
                }
            }
        }
    }
}