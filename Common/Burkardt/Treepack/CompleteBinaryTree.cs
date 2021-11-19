using System;
using System.Globalization;

namespace Burkardt.Treepack;

public static class CompleteBinaryTree
{
    public static void cbt_traverse(int depth)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CBT_TRAVERSE traverses a complete binary tree of given depth.
        //
        //  Discussion:
        //
        //    There will be 2^DEPTH terminal nodes of the complete binary tree.
        //
        //    This function traverses the tree, and prints out a binary code of 0's
        //    and 1's each time it encounters a terminal node.  This results in a 
        //    printout of the binary digits from 0 to 2^DEPTH - 1.
        //
        //    The function is intended as a framework to be used to traverse a binary
        //    tree.  Thus, in practice, a user would insert some action when a terminal
        //    node is encountered.
        //
        //    Another use would occur when a combinatorial search is being made, for
        //    example in a knapsack problem.  Each binary string then represents which
        //    objects are to be included in the knapsack.  In that case, the traversal
        //    could be speeded up by noticing cases where a nonterminal node has been
        //    reached, but the knapsack is already full, in which case the only solution
        //    uses none of the succeeding items, or overfull, in which case no solutions
        //    exist that include this initial path segment.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 December 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DEPTH, the depth of the tree.
        //
    {
        const int DOWNLEFT = 1;
        int i;
        const int UP = 3;
        const int UPDOWNRIGHT = 2;

        switch (depth)
        {
            case < 1:
                return;
        }

        int[] b = new int[depth + 1];

        for (i = 0; i <= depth; i++)
        {
            b[i] = 0;
        }

        int p = 0;
        int direction = DOWNLEFT;
        int k = 0;

        for (;;)
        {
            //
            //  Try going in direction DOWNLEFT.
            //
            string cout;
            if (direction == DOWNLEFT)
            {
                p += 1;
                b[p - 1] = 0;
                if (p < depth)
                {
                    cout = "           ";
                    for (i = 0; i < p; i++)
                    {
                        cout += b[i];
                    }

                    Console.WriteLine(cout);
                }
                else
                {
                    cout = "  (  " + k.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  ";
                    for (i = 0; i < p; i++)
                    {
                        cout += b[i];
                    }

                    Console.WriteLine(cout);
                    k += 1;
                    direction = UPDOWNRIGHT;
                }
            }

            //
            //  Try going in direction UPDOWNRIGHT.
            //
            if (direction == UPDOWNRIGHT)
            {
                b[p - 1] = +1;
                if (p < depth)
                {
                    cout = "           ";
                    for (i = 0; i < p; i++)
                    {
                        cout += b[i];
                    }

                    Console.WriteLine(cout);
                    direction = DOWNLEFT;
                }
                else
                {
                    cout = "  )  " + k.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  ";
                    for (i = 0; i < p; i++)
                    {
                        cout += b[i];
                    }

                    Console.WriteLine(cout);
                    k += 1;
                    direction = UP;
                }
            }

            //
            //  Try going in direction UP.
            //
            if (direction != UP)
            {
                continue;
            }

            p -= 1;
            if (1 <= p)
            {
                cout = "           ";
                for (i = 0; i < p; i++)
                {
                    cout += b[i];
                }

                Console.WriteLine(cout);
                direction = b[p - 1] switch
                {
                    0 => UPDOWNRIGHT,
                    _ => direction
                };
            }
            else
            {
                break;
            }
        }
    }
}