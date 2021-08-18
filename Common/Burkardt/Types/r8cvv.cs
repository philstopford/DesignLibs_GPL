using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static int[] r8cvv_offset(int m, int[] nr)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8CVV_OFFSET determines the row offsets of an R8CVV.
            //
            //  Discussion:
            //
            //    An R8CVV is a "vector of vectors" of R8's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 November 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the number of rows in the array.
            //
            //    Input, int NR[M], the row sizes.
            //
            //    Output, int R8CVV_OFFSET[M+1], the row offsets.
            //
        {
            int i;
            int[] roff;

            roff = new int[m + 1];

            roff[0] = 0;
            for (i = 0; i < m; i++)
            {
                roff[i + 1] = roff[i] + nr[i];
            }

            return roff;
        }

        public static void r8cvv_print(int mn, double[] a, int m, int[] roff, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8CVV_PRINT prints an R8CVV.
            //
            //  Discussion:
            //
            //    An R8CVV is a "vector of vectors" of R8's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 November 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int MN, the size of the cell array.
            //
            //    Input, double A[MN], the cell array.
            //
            //    Input, int M, the number of rows in the array.
            //
            //    Input, int ROFF[M+1], the row offsets.
            //
            //    Input, string TITLE, a title.
            //
        {
            int i;
            int k;
            int k1;
            int k2;
            int khi;
            int klo;

            Console.WriteLine("");
            Console.WriteLine(title + "");
            Console.WriteLine("");

            for (i = 0; i < m; i++)
            {
                k1 = roff[i];
                k2 = roff[i + 1];

                for (klo = k1; klo < k2; klo = klo + 5)
                {
                    string cout = "";
                    khi = Math.Min(klo + 5, k2);
                    if (klo == k1)
                    {
                        cout += i.ToString().PadLeft(5);
                    }
                    else
                    {
                        cout += "     ";
                    }

                    cout += "  ";
                    for (k = klo; k < khi; k++)
                    {
                        cout += a[k].ToString().PadLeft(14);
                    }

                    Console.WriteLine(cout);
                }
            }
        }

        public static double[] r8cvv_rget_new(int mn, double[] a, int m, int[] roff, int i)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8CVV_RGET_NEW gets row I from an R8CVV.
            //
            //  Discussion:
            //
            //    An R8CVV is a "vector of vectors" of R8's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 November 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int MN, the size of the cell array.
            //
            //    Input, double A[MN], the cell array.
            //
            //    Input, int M, the number of rows in the array.
            //
            //    Input, int ROFF[M+1], the row offsets.
            //
            //    Input, int I, the row.
            //    0 <= I < M.
            //
            //    Output, double R8CVV_RGET_NEW[NR[I]], the value of A(I,*).
            //
        {
            double[] ai;
            int j;
            int k1;
            int k2;
            int nv;

            k1 = roff[i];
            k2 = roff[i + 1];
            nv = k2 - k1;
            ai = new double[nv];
            for (j = 0; j < nv; j++)
            {
                ai[j] = a[k1 + j];
            }

            return ai;
        }

        public static void r8cvv_rset(int mn, double[] a, int m, int[] roff, int i, double[] ai)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8CVV_RSET sets row I from an R8CVV.
            //
            //  Discussion:
            //
            //    An R8CVV is a "vector of vectors" of R8's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 November 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int MN, the size of the cell array.
            //
            //    Input/output, double A[MN], the cell array.
            //
            //    Input, int M, the number of rows in the array.
            //
            //    Input, int ROFF[M+1], the row offsets.
            //
            //    Input, int I, the row.
            //    0 <= I < M.
            //
            //    Input, double AI[NR[I]], the new value of A(I,*).
            //
        {
            int j;
            int k1;
            int k2;
            int nv;

            k1 = roff[i];
            k2 = roff[i + 1];
            nv = k2 - k1;
            for (j = 0; j < nv; j++)
            {
                a[k1 + j] = ai[j];
            }

        }
    }
}