using System;

namespace Burkardt.RankingNS
{
    public static partial class Ranking
    {
        public static bool tree_check(int n, int[] t)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    TREE_CHECK checks a tree.
            // 
            //  Licensing:
            // 
            //    This code is distributed under the GNU LGPL license.
            // 
            //  Modified:
            // 
            //    29 July 2011
            // 
            //  Author:
            // 
            //    John Burkardt
            // 
            //  Reference:
            // 
            //    Donald Kreher, Douglas Simpson,
            //    Combinatorial Algorithms,
            //    CRC Press, 1998,
            //    ISBN: 0-8493-3988-X,
            //    LC: QA164.K73.
            // 
            //  Parameters:
            // 
            //    Input, int N, the number of nodes in the tree.
            //    N must be positive.
            // 
            //    Input, int T[2*(N-1)], describes the edges of the tree
            //    as pairs of nodes.
            //
            //    Output, bool TREE_CHECK.
            //    TRUE, the data is legal.
            //    FALSE, the data is not legal.
        {
            bool check;
            int[] d;
            int i;
            int j;
            int k;
            int x;
            int y;

            check = true;

            if (n < 1)
            {
                check = false;
                return check;
            }

            for (i = 0; i < 2; i++)
            {
                for (j = 0; j < n - 1; j++)
                {
                    if (t[i + j * 2] < 1 || n < t[i + j * 2])
                    {
                        check = false;
                        return check;
                    }
                }
            }

            // 
            //  Compute the degree of each node.
            //
            d = edge_degree(n, n - 1, t);
            // 
            //  Delete a node of degree 1, N-1 times.
            // 
            for (k = 1; k <= n - 1; k++)
            {
                x = 1;

                while (d[x - 1] != 1)
                {
                    x = x + 1;
                    if (n < x)
                    {
                        check = false;
                        return check;
                    }
                }

                // 
                //  Find its neighbor.
                // 
                j = 1;

                for (;;)
                {
                    if (t[0 + (j - 1) * 2] == x)
                    {
                        y = t[1 + (j - 1) * 2];
                        break;
                    }

                    if (t[1 + (j - 1) * 2] == x)
                    {
                        y = t[0 + (j - 1) * 2];
                        break;
                    }

                    j = j + 1;

                    if (n - 1 < j)
                    {
                        check = false;
                        return check;
                    }
                }

                // 
                //  Delete the edge.
                // 
                t[0 + (j - 1) * 2] = -t[0 + (j - 1) * 2];
                t[1 + (j - 1) * 2] = -t[1 + (j - 1) * 2];

                d[x - 1] = d[x - 1] - 1;
                d[y - 1] = d[y - 1] - 1;
            }

            for (j = 0; j < n - 1; j++)
            {
                for (i = 0; i < 2; i++)
                {
                    t[i + j * 2] = -t[i + j * 2];
                }
            }

            return check;
        }

        public static int tree_enum(int n)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    TREE_ENUM enumerates the trees on N nodes.
            // 
            //  Licensing:
            // 
            //    This code is distributed under the GNU LGPL license.
            // 
            //  Modified:
            // 
            //    24 July 2011
            // 
            //  Author:
            // 
            //    John Burkardt
            // 
            //  Parameters:
            // 
            //    Input, int N, the number of nodes in each tree.
            //    N must normally be at least 3, but for this routine,
            //    any value of N is allowed.
            // 
            //    Output, int TREE_ENUM, the number of distinct elements.
            // 
        {
            int value;

            if (n < 1)
            {
                value = 0;
            }
            else if (n == 1)
            {
                value = 1;
            }
            else if (n == 2)
            {
                value = 1;
            }
            else
            {
                value = (int)Math.Pow(n, n - 2);
            }

            return value;
        }

        public static int tree_rank(int n, int[] t)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    TREE_RANK ranks a tree.
            // 
            //  Licensing:
            // 
            //    This code is distributed under the GNU LGPL license.
            // 
            //  Modified:
            // 
            //    25 July 2011
            // 
            //  Author:
            // 
            //    John Burkardt
            // 
            //  Reference:
            // 
            //    Donald Kreher, Douglas Simpson,
            //    Combinatorial Algorithms,
            //    CRC Press, 1998,
            //    ISBN: 0-8493-3988-X,
            //    LC: QA164.K73.
            // 
            //  Parameters:
            // 
            //    Input, int N, the number of nodes in the tree.
            //    N must be at least 3.
            // 
            //    Input, int T[2*(N-1)], describes the edges of the tree
            //    as pairs of nodes.
            // 
            //    Output, int RANK, the rank of the tree.
            // 
        {
            bool check;
            int[] p;
            int rank;
            // 
            //  Check the tree.
            // 
            check = tree_check(n, t);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("TREE_RANK - Fatal error!");
                Console.WriteLine("  The tree is illegal.");
                return (1);
            }

            // 
            //  Convert the tree to a Pruefer code.
            // 
            p = tree_to_pruefer(n, t);
            // 
            //  Find the rank of the Pruefer code.
            // 
            rank = pruefer_rank(n, p);
            
            return rank;
        }

        public static void tree_successor(int n, ref int[] t, ref int rank )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    TREE_SUCCESSOR returns the successor of a tree.
        // 
        //  Licensing:
        // 
        //    This code is distributed under the GNU LGPL license.
        // 
        //  Modified:
        // 
        //    26 July 2011
        // 
        //  Author:
        // 
        //    John Burkardt
        // 
        //  Reference:
        // 
        //    Donald Kreher, Douglas Simpson,
        //    Combinatorial Algorithms,
        //    CRC Press, 1998,
        //    ISBN: 0-8493-3988-X,
        //    LC: QA164.K73.
        // 
        //  Parameters:
        // 
        //    Input, int N, the number of nodes in the tree.
        //    N must be at least 3.
        // 
        //    Input/output, int T[2*(N-1)], describes the edges of the
        //    tree as pairs of nodes.  On output, the input tree has been replaced
        //    by its successor.
        // 
        //    Input/output, int &RANK, the rank of the tree.
        // 
        {
            bool check;
            int i;
            int[] p;
            // 
            //  Return the first element.
            // 
            if (rank == -1)
            {
                p = new int[n - 2];

                for (i = 0; i < n - 2; i++)
                {
                    p[i] = 1;
                }

                pruefer_to_tree(n, p, ref t);
                rank = 0;
                return;
            }

            // 
            //  Check the tree.
            // 
            check = tree_check(n, t);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("TREE_SUCCESSOR - Fatal error!");
                Console.WriteLine("  The tree is illegal.");
                return;
            }

            // 
            //  Convert the tree to a Pruefer code.
            // 
            p = tree_to_pruefer(n, t);
            // 
            //  Find the successor of the Pruefer code.
            // 
            pruefer_successor(n, ref p, ref rank);
            // 
            //  Convert the Pruefer code to the tree.
            // 
            pruefer_to_tree(n, p, ref t);
            
            return;
        }

        public static int[] tree_to_pruefer(int n, int[] t)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    TREE_TO_PRUEFER converts a tree to a Pruefer code.
            // 
            //  Licensing:
            // 
            //    This code is distributed under the GNU LGPL license.
            // 
            //  Modified:
            // 
            //    25 July 2011
            // 
            //  Author:
            // 
            //    John Burkardt
            // 
            //  Reference:
            // 
            //    Donald Kreher, Douglas Simpson,
            //    Combinatorial Algorithms,
            //    CRC Press, 1998,
            //    ISBN: 0-8493-3988-X,
            //    LC: QA164.K73.
            // 
            //  Parameters:
            // 
            //    Input, int N, the number of nodes in the tree.
            //    N must be positive.
            // 
            //    Input, int T[2*(N-1)], describes the edges of the tree
            //    as pairs of nodes.
            // 
            //    Output, int TREE_TO_PRUEFER[N-2], the Pruefer code for the tree.
            // 
        {
            bool check;
            int[] d;
            int i;
            int j;
            int k;
            int[] p;
            int x;
            int y;
            // 
            //  Check.
            // 
            check = tree_check(n, t);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("TREE_TO_PRUEFER - Fatal error!");
                Console.WriteLine("  The tree is illegal.");
                return null;
            }

            // 
            //  Compute the degree of each node.
            // 
            d = edge_degree(n, n - 1, t);

            p = new int[n - 2];

            for (j = 1; j <= n - 2; j++)
            {
                // 
                //  Find a node of degree 1.
                // 
                x = n;
                while (d[x - 1] != 1)
                {
                    x = x - 1;
                }

                // 
                //  Find its neighbor.
                // 
                k = 1;

                for (;;)
                {
                    if (t[0 + (k - 1) * 2] == x)
                    {
                        y = t[1 + (k - 1) * 2];
                        break;
                    }

                    if (t[1 + (k - 1) * 2] == x)
                    {
                        y = t[0 + (k - 1) * 2];
                        break;
                    }

                    k = k + 1;
                }

                // 
                //  Store the neighbor.
                // 
                p[j - 1] = y;
                // 
                //  Delete the edge from the tree.
                // 
                d[x - 1] = d[x - 1] - 1;
                d[y - 1] = d[y - 1] - 1;

                t[0 + (k - 1) * 2] = -t[0 + (k - 1) * 2];
                t[1 + (k - 1) * 2] = -t[1 + (k - 1) * 2];
            }

            // 
            //  Remove the negative signs from the first N-2 columns of the tree.
            // 
            for (j = 0; j < n - 2; j++)
            {
                for (i = 0; i < 2; i++)
                {
                    t[i + j * 2] = -t[i + j * 2];
                }
            }

            return p;
        }

        public static int[] tree_unrank(int rank, int n)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    TREE_UNRANK unranks a tree.
            // 
            //  Licensing:
            // 
            //    This code is distributed under the GNU LGPL license.
            // 
            //  Modified:
            // 
            //    27 July 2011
            // 
            //  Author:
            // 
            //    John Burkardt
            // 
            //  Reference:
            // 
            //    Donald Kreher, Douglas Simpson,
            //    Combinatorial Algorithms,
            //    CRC Press, 1998,
            //    ISBN: 0-8493-3988-X,
            //    LC: QA164.K73.
            // 
            //  Parameters:
            // 
            //    Input, int RANK, the rank of the tree.
            // 
            //    Input, int N, the number of nodes in the tree.
            //    N must be at least 3.
            // 
            //    Output, int T[2*(N-1)], describes the edges of the tree
            //    as pairs of nodes.
            // 
        {
            int[] p;
            int[] t;
            int tree_num;
            // 
            //  Check.
            // 
            if (n < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("TREE_UNRANK - Fatal error!");
                Console.WriteLine("  Input N is illegal.");
                return null;
            }

            tree_num = tree_enum(n);

            if (rank < 0 || tree_num < rank)
            {
                Console.WriteLine("");
                Console.WriteLine("TREE_UNRANK - Fatal error!");
                Console.WriteLine("  The input rank is illegal.");
                return null;
            }

            // 
            //  Unrank the Pruefer code.
            // 
            p = pruefer_unrank(rank, n);
            // 
            //  Convert the Pruefer code to a tree.
            // 
            t = pruefer_to_tree_new(n, p);

            return t;
        }
    }
}