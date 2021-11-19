using System;
using Burkardt.Sequence;
using Burkardt.Treepack;
using Burkardt.Types;

namespace TreepackTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TREEPACK_TEST.
        //
        //  Discussion:
        //
        //    TREEPACK_TEST tests the TREEPACK library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TREEPACK_TEST");
        Console.WriteLine("  Test the TREEPACK library.");

        test005();
        test006();
        test01();
        test02();
        test025();
        test03();
        test04();
        test05();
        test06();
        test07();
        test08();
        test09();
        test10();
        test11();
        test12();
        test13();
        test14();
        test15();
        test16();
        test17();
        test18();

        Console.WriteLine("");
        Console.WriteLine("TREEPACK_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test005()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST005 tests CATALAN and CATALAN_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int c = 0;
        int[] c2;
        int n = 0;
        int n_data;

        Console.WriteLine("");
        Console.WriteLine("TEST005");
        Console.WriteLine("  CATALAN computes Catalan numbers.");
        Console.WriteLine("  CATALAN_VALUES returns some exact values.");
        Console.WriteLine("");
        Console.WriteLine("  N  exact C(I)  computed C(I)");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Burkardt.Values.Catalan.catalan_values(ref n_data, ref n, ref c);

            if (n_data == 0)
            {
                break;
            }

            c2 = Catalan.catalan(n);

            Console.WriteLine("  " + n.ToString().PadLeft(4)
                                   + "  " + c.ToString().PadLeft(6)
                                   + "  " + c2[n].ToString().PadLeft(6) + "");

        }
    }

    private static void test006()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST006 tests CBT_TRAVERSE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int depth = 4;

        Console.WriteLine("");
        Console.WriteLine("TEST006");
        Console.WriteLine("  CBT_TRAVERSE traverses a complete binary tree.");
        Console.WriteLine("");
        Console.WriteLine("  For this demonstration, we simply print our path.");
        Console.WriteLine("  The tree depth is " + depth + "");
        Console.WriteLine("");

        CompleteBinaryTree.cbt_traverse(depth);
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests PRUEFER_TO_TREE_ARC.
        //
        //  Discussion:
        //
        //    The tree is
        //
        //          5
        //          |
        //    2-3-6-8-1-9
        //      |   |
        //      7   4
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] code = {1, 3, 8, 8, 3, 6, 8};
        int[] inode;
        int[] jnode;
        int nnode = 9;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  PRUEFER_TO_TREE_ARC computes a tree from its Pruefer code.");
        Console.WriteLine("");
        Console.WriteLine("          5");
        Console.WriteLine("          |");
        Console.WriteLine("    2-3-6-8-1-9");
        Console.WriteLine("      |   |");
        Console.WriteLine("      7   4");

        typeMethods.i4vec_print(nnode - 2, code, "  The Pruefer code:");

        inode = new int[nnode - 1];
        jnode = new int[nnode - 1];

        Pruefer.pruefer_to_tree_arc(nnode, code, ref inode, ref jnode);

        Graph.graph_arc_print(nnode - 1, inode, jnode, "  The graph:");
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests PRUEFER_TO_TREE_2_NEW.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int NNODE = 9;

        int[] code = {1, 3, 8, 8, 3, 6, 8, 0, 0};
        int[] itree;
        int nnode = NNODE;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  PRUEFER_TO_TREE_2_NEW produces a tree from its Pruefer code");

        typeMethods.i4vec_print(nnode - 2, code, "  The Pruefer code:");

        itree = Pruefer.pruefer_to_tree_2_new(nnode, code);

        typeMethods.i4vec_print(nnode - 1, itree, "  The edge list of the tree:");
    }

    private static void test025()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST025 tests PRUEFER_TO_TREE_2.
        //
        //  Discussion:
        //
        //    This example is used to illustrate the Nijenhuis and Wilf algorithm
        //    LBLTRE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int NNODE = 4;

        int[] code = new int[NNODE];
        int i;
        int[] itree = new int[NNODE];
        int j;
        int nnode = NNODE;

        Console.WriteLine("");
        Console.WriteLine("TEST025");
        Console.WriteLine("  PRUEFER_TO_TREE_2 produces a tree from its Pruefer code");
        Console.WriteLine("");
        Console.WriteLine("   Code      Tree");
        Console.WriteLine("");
        for (j = 1; j <= nnode; j++)
        {
            code[1] = j;
            for (i = 1; i <= nnode; i++)
            {
                code[0] = i;
                Pruefer.pruefer_to_tree_2(nnode, code, ref itree);
                Console.WriteLine("  " + code[0].ToString().PadLeft(2)
                                       + "  " + code[1].ToString().PadLeft(2)
                                       + "  "
                                       + "  " + itree[0].ToString().PadLeft(2)
                                       + "  " + itree[1].ToString().PadLeft(2)
                                       + "  " + itree[2].ToString().PadLeft(2) + "");
            }
        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests TREE_ARC_TO_PRUEFER.
        //
        //  Discussion:
        //
        //    The tree is
        //
        //          5
        //          |
        //    2-3-6-8-1-9
        //      |   |
        //      7   4
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] code;
        int[] inode = {2, 3, 3, 6, 8, 8, 8, 1};
        int[] jnode = {3, 7, 6, 8, 4, 5, 1, 9};
        int nnode = 9;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  TREE_ARC_TO_PRUEFER: Tree => Pruefer code");
        Console.WriteLine("");
        Console.WriteLine("          5");
        Console.WriteLine("          |");
        Console.WriteLine("    2-3-6-8-1-9");
        Console.WriteLine("      |   |");
        Console.WriteLine("      7   4");

        Graph.graph_arc_print(nnode - 1, inode, jnode, "  The graph:");

        code = Tree.tree_arc_to_pruefer(nnode, inode, jnode);

        typeMethods.i4vec_print(nnode - 2, code, "  The Pruefer code:");
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests TREE_ARC_CENTER.
        //
        //  Discussion:
        //
        //    The tree is
        //
        //                4
        //                |
        //    2---3---6---8---1---9
        //        |       |
        //        7       5
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int NNODE = 9;

        int[] center = new int[2];
        int eccent = 0;
        int[] inode = {2, 3, 3, 6, 8, 8, 8, 1};
        int[] jnode = {3, 7, 6, 8, 4, 5, 1, 9};
        int nnode = NNODE;
        int parity = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  TREE_ARC_CENTER computes the center of a tree.");

        Graph.graph_arc_print(nnode - 1, inode, jnode, "  The edge list of the tree:");

        Tree.tree_arc_center(nnode, inode, jnode, ref center, ref eccent, ref parity);

        Console.WriteLine("");
        Console.WriteLine("  Parity = " + parity + "");
        Console.WriteLine("  Eccentricity is " + eccent + "");

        switch (parity)
        {
            case 0:
                Console.WriteLine("  No center node (degenerate case).");
                break;
            case 1:
                Console.WriteLine("  Center node: " + center[0] + "");
                break;
            default:
                Console.WriteLine("  Center nodes: " + center[0] + "  " + center[1] + "");
                break;
        }
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests TREE_ARC_CENTER.
        //
        //  Discussion:
        //
        //    Compare:
        //
        //    2--1--3
        //
        //    1--2--3
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int NNODE = 3;

        int[] center = new int[2];
        int eccent = 0;
        int[] inode = new int[NNODE - 1];
        int[] jnode = new int[NNODE - 1];
        int nnode = NNODE;
        int parity = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  TREE_ARC_CENTER computes the center of a tree.");

        inode[0] = 1;
        inode[1] = 1;
        jnode[0] = 2;
        jnode[1] = 3;

        Graph.graph_arc_print(nnode - 1, inode, jnode, "  The edge list of the tree:");

        Tree.tree_arc_center(nnode, inode, jnode, ref center, ref eccent, ref parity);

        Console.WriteLine("");
        Console.WriteLine("  Parity = " + parity + "");
        Console.WriteLine("  Eccentricity is " + eccent + "");

        switch (parity)
        {
            case 0:
                Console.WriteLine("  No center node (degenerate case).");
                break;
            case 1:
                Console.WriteLine("  Center node: " + center[0] + "");
                break;
            default:
                Console.WriteLine("  Center nodes: " + center[0] + "  " + center[1] + "");
                break;
        }

        inode[0] = 2;
        inode[1] = 2;
        jnode[0] = 1;
        jnode[1] = 3;

        Graph.graph_arc_print(nnode - 1, inode, jnode, "  The edge list of the tree:");

        Tree.tree_arc_center(nnode, inode, jnode, ref center, ref eccent, ref parity);

        Console.WriteLine("");
        Console.WriteLine("  Parity = " + parity + "");
        Console.WriteLine("  Eccentricity is " + eccent + "");

        switch (parity)
        {
            case 0:
                Console.WriteLine("  No center node (degenerate case).");
                break;
            case 1:
                Console.WriteLine("  Center node: " + center[0] + "");
                break;
            default:
                Console.WriteLine("  Center nodes: " + center[0] + "  " + center[1] + "");
                break;
        }
    }

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tests TREE_ARC_CENTER.
        //
        //  Discussion:
        //
        //    The tree is
        //
        //     4     8     7
        //     |     |     |
        //  5--1-----2-----3--9
        //     |     |     |
        //     6     10   11
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int NNODE = 11;

        int[] center = new int[2];
        int eccent = 0;
        int[] inode = {1, 1, 1, 2, 2, 3, 3, 3, 1, 2};
        int[] jnode = {4, 5, 6, 8, 10, 7, 9, 11, 2, 3};
        int nnode = NNODE;
        int parity = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  TREE_ARC_CENTER computes the center of a tree.");

        Graph.graph_arc_print(nnode - 1, inode, jnode, "  The edge list of the tree:");

        Tree.tree_arc_center(nnode, inode, jnode, ref center, ref eccent, ref parity);

        Console.WriteLine("");
        Console.WriteLine("  Parity = " + parity + "");
        Console.WriteLine("  Eccentricity is " + eccent + "");

        switch (parity)
        {
            case 0:
                Console.WriteLine("  No center node (degenerate case).");
                break;
            case 1:
                Console.WriteLine("  Center node: " + center[0] + "");
                break;
            default:
                Console.WriteLine("  Center nodes: " + center[0] + "  " + center[1] + "");
                break;
        }
    }

    private static void test07()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 tests TREE_ARC_DIAM.
        //
        //  Discussion:
        //
        //    The tree is:
        //
        //                5
        //                |
        //    2---3---6---8---1---9
        //        |       |
        //        7       4
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int diam = 0;
        int[] inode = {2, 3, 3, 6, 8, 8, 8, 1};
        int[] jnode = {3, 7, 6, 8, 4, 5, 1, 9};
        int[] label = new int[9];
        int nnode = 9;
        int nnode1 = 0;
        int nnode2 = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST07");
        Console.WriteLine("  TREE_ARC_DIAM computes the diameter of a tree.");

        Graph.graph_arc_print(nnode - 1, inode, jnode, "  The edge list of the tree:");

        Tree.tree_arc_diam(nnode, inode, jnode, ref diam, ref label, ref nnode1, ref nnode2);

        Console.WriteLine("");
        Console.WriteLine("  This tree has a diameter of " + diam + "");
        Console.WriteLine("  between nodes " + nnode1 + " and " + nnode2 + "");

        typeMethods.i4vec_print(nnode, label, "  Nodes and labels:");
    }

    private static void test08()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 tests TREE_ARC_RANDOM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int[] icode = new int[2];
        int[] inode = new int[3];
        int[] jnode = new int[3];
        int nnode = 4;
        int seed;

        seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("TEST08");
        Console.WriteLine("  TREE_ARC_RANDOM produces a random labeled");
        Console.WriteLine("  tree and its Pruefer code.");
        Console.WriteLine("");

        for (i = 1; i <= 5; i++)
        {
            Tree.tree_arc_random(nnode, ref seed, ref icode, ref inode, ref jnode);

            Graph.graph_arc_print(nnode - 1, inode, jnode, "  The random tree:");

            typeMethods.i4vec_print(nnode - 2, icode, "  The Pruefer code:");
        }
    }

    private static void test09()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST09 tests TREE_ENUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int nnode;
        int num;

        Console.WriteLine("");
        Console.WriteLine("TEST09");
        Console.WriteLine("  TREE_ENUM enumerates the labeled trees on a given");
        Console.WriteLine("  number of nodes.");
        Console.WriteLine("");

        for (nnode = 0; nnode <= 10; nnode++)
        {
            num = Tree.tree_enum(nnode);
            Console.WriteLine("  " + nnode.ToString().PadLeft(8)
                                   + "  " + num.ToString().PadLeft(10) + "");
        }
    }

    private static void test10()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST10 tests TREE_PARENT_NEXT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int NNODE = 4;

        int[] icode = new int[NNODE];
        int[] itree = new int[NNODE];
        bool more;
        int nnode = NNODE;

        Console.WriteLine("");
        Console.WriteLine("TEST10");
        Console.WriteLine("  TREE_PARENT_NEXT finds all labeled trees of a given ");
        Console.WriteLine("  order, and their Pruefer codes.");
        Console.WriteLine("");
        Console.WriteLine("  Pruefer code     Tree");
        Console.WriteLine("");

        more = false;

        TreeNextData data = new();
            
        for (;;)
        {
            Tree.tree_parent_next(ref data, nnode, ref icode, ref itree, ref more);

            Console.WriteLine("  " + icode[0]
                                   + "  " + icode[1]
                                   + "            "
                                   + "  " + itree[0].ToString().PadLeft(2)
                                   + "  " + itree[1].ToString().PadLeft(2)
                                   + "  " + itree[2].ToString().PadLeft(2) + "");

            if (!more)
            {
                break;
            }
        }
    }

    private static void test11()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST11 tests TREE_RB_ENUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int nnode;
        int num;

        Console.WriteLine("");
        Console.WriteLine("TEST11");
        Console.WriteLine("  TREE_RB_ENUM enumerates the rooted binary trees on a ");
        Console.WriteLine("  given number of nodes.");
        Console.WriteLine("");

        for (nnode = 0; nnode <= 11; nnode++)
        {
            num = Tree.tree_rb_enum(nnode);

            Console.WriteLine("  " + nnode.ToString().PadLeft(8)
                                   + "  " + num.ToString().PadLeft(8) + "");
        }
    }

    private static void test12()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST12 tests TREE_RB_LEX_NEXT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] a = new int[11];
        int i;
        int j;
        bool more;
        int n = 11;

        Console.WriteLine("");
        Console.WriteLine("TEST12");
        Console.WriteLine("  TREE_RB_LEX_NEXT produces all rooted binary trees with");
        Console.WriteLine("  a given number of nodes, in lexicographic order, using");
        Console.WriteLine("  the preorder traversal representation.");
        Console.WriteLine("");
        Console.WriteLine("  The number of nodes N = " + n + "");
        Console.WriteLine("");

        more = false;
        i = 0;

        for (;;)
        {
            Tree.tree_rb_lex_next(n, a, ref more);

            if (!more)
            {
                break;
            }

            i += 1;
            string cout = "  " + i.ToString().PadLeft(2) + "  ";
            for (j = 0; j < n; j++)
            {
                cout += a[j];
            }

            Console.WriteLine(cout);
        }
    }

    private static void test13()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST13 tests TREE_RB_LEX_NEXT, TREE_RB_TO_PARENT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] a = new int[11];
        int i;
        int j;
        bool more;
        int n = 11;
        int[] parent;

        Console.WriteLine("");
        Console.WriteLine("TEST13");
        Console.WriteLine("  TREE_RB_LEX_NEXT produces all rooted binary trees with");
        Console.WriteLine("  a given number of nodes, in lexicographic order,");
        Console.WriteLine("  using the preorder traversal representation.");
        Console.WriteLine("  TREE_RB_TO_PARENT converts the preorder traversal form");
        Console.WriteLine("  to the more comprehensible parent node representation.");
        Console.WriteLine("");
        Console.WriteLine("  The number of nodes N = " + n + "");
        Console.WriteLine("");

        more = false;
        i = 0;

        for (;;)
        {
            Tree.tree_rb_lex_next(n, a, ref more);

            if (!more)
            {
                break;
            }

            parent = Tree.tree_rb_to_parent(n, a);

            i += 1;
            string cout = "  " + i.ToString().PadLeft(2) + "  ";
            for (j = 0; j < n; j++)
            {
                cout += parent[j].ToString().PadLeft(3);
            }

            Console.WriteLine(cout);

        }
    }

    private static void test14()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST14 tests TREE_RB_YULE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] a = new int[11];
        int i;
        int j;
        int n;
        const int N_MAX = 11;
        int seed;

        seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("TEST14");
        Console.WriteLine("  TREE_RB_YULE carries out one step of the Yule model");
        Console.WriteLine("  on a rooted binary tree stored in preorder traversal form.");
        Console.WriteLine("");
        Console.WriteLine("  Each call adds two children to an arbitary leaf.");

        for (i = 1; i <= 5; i++)
        {
            Console.WriteLine("");
            Console.WriteLine("  Simulation " + i + "");
            Console.WriteLine("");
            Console.WriteLine("  Nodes  Preorder code");
            Console.WriteLine("");

            n = 0;

            for (;;)
            {
                Tree.tree_rb_yule(ref n, ref seed, a);

                string cout = "  " + n.ToString().PadLeft(2) + "  ";
                for (j = 0; j < n; j++)
                {
                    cout += a[j];
                }

                Console.WriteLine(cout);

                if (N_MAX < n + 2)
                {
                    break;
                }
            }
        }
    }

    private static void test15()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST15 tests TREE_ROOTED_CODE.
        //
        //  Discussion:
        //
        //      1
        //      |:
        //      | :
        //      |  :
        //      2   3
        //     :|:  |:
        //    4 5 6 7 8
        //     :|  :
        //    9 10  11
        //      |
        //      12
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] code;
        int nnode = 12;
        int[] parent = {0, 1, 1, 2, 2, 2, 3, 3, 5, 5, 6, 10};

        Console.WriteLine("");
        Console.WriteLine("TEST15");
        Console.WriteLine("  TREE_ROOTED_CODE: code of a rooted tree.");

        typeMethods.i4vec_print(nnode, parent, "  Parent vector for tree:");

        code = Tree.tree_rooted_code(nnode, parent);

        typeMethods.i4vec_print(2 * nnode, code, "  The tree code:");
    }

    private static void test16()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST16 tests TREE_ROOTED_DEPTH.
        //
        //  Discussion:
        //
        //      1
        //      |:
        //      | :
        //      |  :
        //   4--2-+ 3--8
        //      | |  |
        //   9--5 |  7 
        //      | 6-----11
        //      10  
        //      |
        //      12
        //
        //    Depths
        //
        //    1  2  3  4  5  6  7  8  9 10 11 12
        //    0  1  1  2  2  2  2  2  3  3  3  4
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int depth = 0;
        int[] depth_node;
        int nnode = 12;
        int[] parent = {0, 1, 1, 2, 2, 2, 3, 3, 5, 5, 6, 10};

        Console.WriteLine("");
        Console.WriteLine("TEST16");
        Console.WriteLine("  TREE_ROOTED_DEPTH: depth of a rooted tree.");

        typeMethods.i4vec_print(nnode, parent, "  Parent vector for tree:");

        depth_node = new int[nnode];

        Tree.tree_rooted_depth(nnode, parent, ref depth, ref depth_node);

        typeMethods.i4vec_print(nnode, depth_node, "  Individual node depths:");

        Console.WriteLine("");
        Console.WriteLine("  Overall rooted tree depth: " + depth + "");
    }

    private static void test17()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST17 tests TREE_ROOTED_ENUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int nnode = 10;
        int[] ntree;

        Console.WriteLine("");
        Console.WriteLine("TEST17");
        Console.WriteLine("  TREE_ROOTED_ENUM counts unlabeled rooted trees.");

        ntree = Tree.tree_rooted_enum(nnode);

        typeMethods.i4vec_print(nnode, ntree,
            "  Number of trees with given number of nodes:");
    }

    private static void test18()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST18 tests TREE_ROOTED_RANDOM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int[] itree;
        int j;
        int nnode = 5;
        int seed;

        seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("TEST18");
        Console.WriteLine("  TREE_ROOTED_RANDOM: random unlabeled rooted trees.");
        Console.WriteLine("");
        Console.WriteLine("  Selecting random trees, rooted at 1");
        Console.WriteLine("  Number of nodes is NNODE = " + nnode + "");
        Console.WriteLine("");
        Console.WriteLine("  Each tree is described by the nodes that");
        Console.WriteLine("  connect nodes 2 through NNODE.");
        Console.WriteLine("");
        for (i = 1; i <= 5; i++)
        {
            itree = Tree.tree_rooted_random(nnode, ref seed);

            string cout = "  ";
            for (j = 1; j < nnode; j++)
            {
                cout += itree[j].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
        }
    }
}