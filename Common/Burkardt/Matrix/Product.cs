using Burkardt.Types;

namespace Burkardt.MatrixNS;

public static class MatrixProduct
{
    public static void matrix_product_opt(int n, int[] rank, ref int cost, ref int[] order)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MATRIX_PRODUCT_OPT determines the optimal cost of a matrix product.
        //
        //  Discussion:
        //
        //    The cost of multiplying an LxM matrix by an M by N matrix is
        //    assessed as L*M*N.
        //
        //    Any particular order of multiplying a set of N matrices is equivalent
        //    to parenthesizing an expression of N objects.
        //
        //    The actual number of ways of parenthesizing an expression
        //    of N objects is C(N), the N-th Catalan number.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 June 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Robert Sedgewick,
        //    Algorithms,
        //    Addison-Wesley, 1984, pages 486-489.
        //
        //  Parameters:
        //
        //    Input, int N, the number of matrices to be multiplied.
        //
        //    Input, int RANK[N+1], the rank information for the matrices.
        //    Matrix I has RANK[I] rows and RANK[I+1] columns.
        //
        //    Output, int &COST, the cost of the multiplication if the optimal
        //    order is used.
        //
        //    Output, int ORDER[N-1], indicates the order in which the N-1
        //    multiplications are to be carried out.  ORDER[0] is the first
        //    multiplication to do, and so on.
        //
    {
        const int STACK_MAX = 100;

        int i;
        int j;
        int[] stack = new int[STACK_MAX];
        //
        //  Initialize the cost matrix.
        //
        int[] best = new int[n * n];
        int[] cost2 = new int[n * n];

        for (i = 0; i < n; i++)
        {
            for (j = 0; j <= i; j++)
            {
                cost2[i + j * n] = 0;
            }

            for (j = i + 1; j < n; j++)
            {
                cost2[i + j * n] = typeMethods.i4_huge();
            }
        }

        //
        //  Initialize the BEST matrix.
        //
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                best[i + j * n] = 0;
            }
        }

        //
        //  Compute the cost and best matrices.
        //
        for (j = 1; j <= n - 1; j++)
        {
            for (i = 1; i <= n - j; i++)
            {
                int k;
                for (k = i + 1; k <= i + j; k++)
                {
                    int cost3 = cost2[i - 1 + (k - 2) * n] + cost2[k - 1 + (i + j - 1) * n]
                                                           + rank[i - 1] * rank[k - 1] * rank[i + j];

                    if (cost3 >= cost2[i - 1 + (i + j - 1) * n])
                    {
                        continue;
                    }

                    cost2[i - 1 + (i + j - 1) * n] = cost3;
                    best[i - 1 + (i + j - 1) * n] = k;
                }
            }
        }

        //
        //  Pick off the optimal cost.
        //
        cost = cost2[0 + (n - 1) * n];
        //
        //  Backtrack to determine the optimal order.
        //
        int stack_num = 0;

        int i1 = 1;
        int i2 = n;

        if (i1 + 1 < i2)
        {
            stack[stack_num] = i1;
            stack_num += 1;
            stack[stack_num] = i2;
            stack_num += 1;
        }

        int step = n - 1;
        //
        //  Take an item off the stack.
        //
        while (0 < stack_num)
        {
            stack_num -= 1;
            int i3 = stack[stack_num];
            stack_num -= 1;
            i1 = stack[stack_num];

            i2 = best[i1 - 1 + (i3 - 1) * n];

            step -= 1;
            order[step] = i2 - 1;
            //
            //  The left chunk is matrices (I1...I2-1)
            //
            if (i1 == i2 - 1)
            {
            }
            else if (i1 + 1 == i2 - 1)
            {
                step -= 1;
                order[step] = i2 - 2;
            }
            else
            {
                stack[stack_num] = i1;
                stack_num += 1;
                stack[stack_num] = i2 - 1;
                stack_num += 1;
            }

            //
            //  The right chunk is matrices (I2...I3)
            //
            if (i2 == i3)
            {
            }
            else if (i2 + 1 == i3)
            {
                step -= 1;
                order[step] = i2;
            }
            else
            {
                stack[stack_num] = i2;
                stack_num += 1;
                stack[stack_num] = i3;
                stack_num += 1;
            }

        }
    }
}