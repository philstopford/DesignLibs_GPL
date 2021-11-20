using System;
using Burkardt.Types;

namespace Burkardt.MatrixNS;

public static class DecimalMatrix
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

        dtop = 0;
        dbot = 1;

        int[] iarray = new int[n];
        //
        //  Compute the next permutation.
        //
        bool more = false;

        for (;;)
        {
            Permutation.perm0_next(n, iarray, ref more, ref even);
            int itop = even switch
            {
                //
                //  The sign of this term depends on the sign of the permutation.
                //
                true => 1,
                _ => -1
            };

            int ibot = 0;
            //
            //  Choose one item from each row, as specified by the permutation,
            //  and multiply them.
            //
            int i;
            for (i = 0; i < n; i++)
            {
                int itop1 = itop;
                int ibot1 = ibot;
                int itop2 = atop[i + iarray[i] * n];
                int ibot2 = abot[i + iarray[i] * n];

                typeMethods.dec_mul(itop1, ibot1, itop2, ibot2, dec_digit, ref itop, ref ibot);

            }

            //
            //  Add this term to the total.
            //
            typeMethods.dec_add(itop, ibot, dtop, dbot, dec_digit, ref dtop, ref dbot);

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
        int i;
        int j;
        int jmin;
        const int ncolum = 80;
        //
        //  Figure out how wide we must make each column.
        //
        int kmax = 0;

        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
            {
                kmax = Math.Max(kmax, typeMethods.dec_width(a[i + j * m], b[i + j * m]));
            }
        }

        kmax += 2;
        int npline = ncolum / kmax;
        //
        //  Now do the printing.
        //
        for (jmin = 0; jmin < n; jmin += npline)
        {
            int jmax = Math.Min(jmin + npline - 1, n - 1);

            switch (title.Length)
            {
                case > 0:
                    Console.WriteLine("");
                    Console.WriteLine(title + "");
                    break;
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
                    int mantissa = a[i + j * m];
                    int exponent = b[i + j * m];
                    string s = typeMethods.dec_to_s(mantissa, exponent);
                    cout += s.PadLeft(kmax) + "  ";
                }

                Console.WriteLine(cout);
            }
        }
    }
}