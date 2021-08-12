namespace Burkardt.MatrixNS
{
    public static class AlternatingSignMatrix
    {
        public static int asm_enum(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ASM_ENUM returns the number of alternating sign matrices of a given order.
            //
            //  Discussion:
            //
            //    N     ASM_NUM
            //
            //    0       1
            //    1       1
            //    2       2
            //    3       7
            //    4      42
            //    5     429
            //    6    7436
            //    7  218348
            //
            //    A direct formula is
            //
            //      ASM_NUM ( N ) = product ( 0 <= I <= N-1 ) ( 3 * I + 1 )! / ( N + I )!
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrices.
            //
            //    Output, int ASM_ENUM, the number of alternating sign matrices of
            //    order N.
            //
        {
            int[] a;
            int asm_num;
            int[] b;
            int[] c;
            int i;
            int nn;

            if (n + 1 <= 0)
            {
                return 0;
            }

            //
            //  Row 1
            //
            if (n + 1 == 1)
            {
                return 1;
            }

            //
            //  Row 2
            //
            if (n + 1 == 2)
            {
                return 1;
            }

            a = new int[n + 1];
            b = new int[n + 1];
            c = new int[n + 1];

            b[0] = 2;
            c[0] = 2;
            a[0] = 1;
            a[1] = 1;
            //
            //  Row 3 and on.
            //
            for (nn = 3; nn <= n; nn++)
            {
                b[nn - 2] = nn;
                for (i = nn - 2; 2 <= i; i--)
                {
                    b[i - 1] = b[i - 1] + b[i - 2];
                }

                b[0] = 2;

                c[nn - 2] = 2;
                for (i = nn - 2; 2 <= i; i--)
                {
                    c[i - 1] = c[i - 1] + c[i - 2];
                }

                c[0] = nn;

                for (i = 2; i <= nn - 1; i++)
                {
                    a[0] = a[0] + a[i - 1];
                }

                for (i = 2; i <= nn; i++)
                {
                    a[i - 1] = a[i - 2] * c[i - 2] / b[i - 2];
                }
            }

            asm_num = 0;
            for (i = 0; i < n; i++)
            {
                asm_num = asm_num + a[i];
            }

            return asm_num;
        }

        public static void asm_triangle(int n, ref int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ASM_TRIANGLE returns a row of the alternating sign matrix triangle.
            //
            //  Discussion:
            //
            //    The first seven rows of the triangle are as follows:
            //
            //          1      2      3      4      5      6     7
            //
            //    0     1
            //    1     1      1
            //    2     2      3      2
            //    3     7     14     14      7
            //    4    42    105    135    105     42
            //    5   429   1287   2002   2002   1287    429
            //    6  7436  26026  47320  56784  47320  26026  7436
            //
            //    For a given N, the value of A(J) represents entry A(I,J) of
            //    the triangular matrix, and gives the number of alternating sign matrices
            //    of order N in which the (unique) 1 in row 1 occurs in column J.
            //
            //    Thus, of alternating sign matrices of order 3, there are
            //    2 with a leading 1 in column 1:
            //
            //      1 0 0  1 0 0
            //      0 1 0  0 0 1
            //      0 0 1  0 1 0
            //
            //    3 with a leading 1 in column 2, and
            //
            //      0 1 0  0 1 0  0 1 0
            //      1 0 0  0 0 1  1-1 1
            //      0 0 1  1 0 0  0 1 0
            //
            //    2 with a leading 1 in column 3:
            //
            //      0 0 1  0 0 1
            //      1 0 0  0 1 0
            //      0 1 0  1 0 0
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the desired row.
            //
            //    Output, int A[N+1], the entries of the row.
            //
        {
            int[] b;
            int[] c;
            int i;
            int nn;
            //
            if (n + 1 <= 0)
            {
                return;
            }

            //
            //  Row 1
            //
            a[0] = 1;

            if (n + 1 == 1)
            {
                return;
            }

            //
            //  Row 2
            //
            a[0] = 1;
            a[1] = 1;

            if (n + 1 == 2)
            {
                return;
            }

            //
            //  Row 3 and on.
            //
            b = new int[n + 1];
            c = new int[n + 1];

            b[0] = 2;
            c[0] = 2;

            for (nn = 3; nn <= n + 1; nn++)
            {

                b[nn - 2] = nn;
                for (i = nn - 2; 2 <= i; i--)
                {
                    b[i - 1] = b[i - 1] + b[i - 2];
                }

                b[0] = 2;

                c[nn - 2] = 2;
                for (i = nn - 2; 2 <= i; i--)
                {
                    c[i - 1] = c[i - 1] + c[i - 2];
                }

                c[0] = nn;

                for (i = 2; i <= nn - 1; i++)
                {
                    a[0] = a[0] + a[i - 1];
                }

                for (i = 2; i <= nn; i++)
                {
                    a[i - 1] = a[i - 2] * c[i - 2] / b[i - 2];
                }

            }
        }
    }
}