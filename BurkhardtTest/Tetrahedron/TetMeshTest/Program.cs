using System;
using System.Globalization;
using Burkardt.TetrahedronNS;
using Burkardt.Types;
using Burkardt.Uniform;
using Tetrahedron = Burkardt.TetrahedronNS.Tetrahedron;

namespace TetMeshTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TET_MESH_TEST
        //
        //  Discussion:
        //
        //    TET_MESH_TEST tests the TET_MESH library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TET_MESH_TEST");
        Console.WriteLine("  Test the TET_MESH library.");

        test001();
        test002();
        test003();
        test004();
        test005();
        test006();
        test007();

        Console.WriteLine("");
        Console.WriteLine("TET_MESH_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test001()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST001 tests R8MAT_SOLVE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 3;
        const int RHS_NUM = 2;

        double[] a =
        {
            1.0, 4.0, 7.0,
            2.0, 5.0, 8.0,
            3.0, 6.0, 0.0,
            14.0, 32.0, 23.0,
            7.0, 16.0, 7.0
        };
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST001");
        Console.WriteLine("  R8MAT_SOLVE solves linear systems.");
        //
        //  Print out the matrix to be inverted.
        //
        typeMethods.r8mat_print(N, N + RHS_NUM, a, "  The linear system:");
        //
        //  Solve the systems.
        //
        int info = typeMethods.r8mat_solve(N, RHS_NUM, ref a);

        if (info != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("  The input matrix was singular.");
            Console.WriteLine("  The solutions could not be computed.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("  The computed solutions:");
        Console.WriteLine("");
        for (i = 0; i < N; i++)
        {
            string cout = "";
            int j;
            for (j = N; j < N + RHS_NUM; j++)
            {
                cout += a[i + j * N].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  ";
            }

            Console.WriteLine(cout);
        }
    }

    private static void test002()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST002 tests TETRAHEDRON_ORDER4_PHYSICAL_TO_REFERENCE,
        //    TETRAHEDRON_ORDER4_REFERENCE_TO_PHYSICAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 10;

        int j;
        double[] phy = new double[3 * N];
        double[] ref_ = new double[3 * N];
        double[] ref2 = new double[3 * N];
        double[] t =
        {
            5.0, 0.0, 0.0,
            8.0, 0.0, 0.0,
            5.0, 2.0, 0.0,
            6.0, 1.0, 2.0
        };

        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("TEST002");
        Console.WriteLine("  For an order 4 tetrahedron,");
        Console.WriteLine("  TETRAHEDRON_ORDER4_PHYSICAL_TO_REFERENCE");
        Console.WriteLine("  maps a physical point to a reference point.");
        Console.WriteLine("  TETRAHEDRON_ORDER4_REFERENCE_TO_PHYSICAL ");
        Console.WriteLine("  maps a reference point to a physical point.");
        Console.WriteLine("");
        Console.WriteLine("     ( R, S, T )          ==>  ( X, Y, Z )           ==> ( R2, S2, T2 )");
        Console.WriteLine("");

        Tetrahedron.tetrahedron_reference_sample(N, ref seed, ref ref_);

        Tetrahedron.tetrahedron_order4_reference_to_physical(t, N, ref_, ref phy);
        Tetrahedron.tetrahedron_order4_physical_to_reference(t, N, phy, ref ref2);

        for (j = 0; j < N; j++)
        {
            Console.WriteLine("  " + ref_[0 + j * 3].ToString("0.####").PadLeft(8)
                                   + "  " + ref_[1 + j * 3].ToString("0.####").PadLeft(8)
                                   + "  " + ref_[2 + j * 3].ToString("0.####").PadLeft(8)
                                   + "  "
                                   + "  " + phy[0 + j * 3].ToString("0.####").PadLeft(8)
                                   + "  " + phy[1 + j * 3].ToString("0.####").PadLeft(8)
                                   + "  " + phy[2 + j * 3].ToString("0.####").PadLeft(8)
                                   + "  "
                                   + "  " + ref2[0 + j * 3].ToString("0.####").PadLeft(8)
                                   + "  " + ref2[1 + j * 3].ToString("0.####").PadLeft(8)
                                   + "  " + ref2[2 + j * 3].ToString("0.####").PadLeft(8) + "");
        }
    }

    private static void test003()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST003 tests TETRAHEDRON_ORDER10_TO_ORDER4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int node_num1 = 0;
        int node_num2 = 0;
        int tet_num1 = 0;
        int tet_num2 = 0;
        const int tet_order1 = 10;
        const int tet_order2 = 4;

        Console.WriteLine("");
        Console.WriteLine("TEST003");
        Console.WriteLine("  For an order 10 tet mesh,");
        Console.WriteLine("  TETRAHEDRON_ORDER10_TO_ORDER4");
        Console.WriteLine("  makes a linear (order 4) tet mesh by using");
        Console.WriteLine("  the existing nodes, and replacing each");
        Console.WriteLine("  quadratic tetrahedron by 8 linear tetrahedrons.");

        TetMesh.tet_mesh_order10_example_size(ref node_num1, ref tet_num1);

        double[] node_xyz = new double[3 * node_num1];
        int[] tet_node1 = new int[tet_order1 * tet_num1];

        TetMesh.tet_mesh_order10_example_set(node_num1, tet_num1,
            ref node_xyz, ref tet_node1);

        typeMethods.i4mat_transpose_print_some(tet_order1, tet_num1, tet_node1,
            1, 1, tet_order1, 5, "  First 5 quadratic tetrahedrons:");

        TetMesh.tet_mesh_order10_to_order4_size(node_num1, tet_num1,
            ref node_num2, ref tet_num2);

        Console.WriteLine("");
        Console.WriteLine("  Quadratic mesh size is       " + tet_num1 + "");
        Console.WriteLine("  Linearized mesh size will be " + tet_num2 + "");

        int[] tet_node2 = new int[tet_order2 * tet_num2];

        TetMesh.tet_mesh_order10_to_order4_compute(tet_num1, tet_node1,
            tet_num2, ref tet_node2);

        typeMethods.i4mat_transpose_print_some(tet_order2, tet_num2, tet_node2,
            1, 1, tet_order2, 5, "  First 5 linear tetrahedrons:");
    }

    private static void test004()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST004 tests TETRAHEDRON_ORDER10_TO_ORDER4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int node_num = 0;
        int tet_num = 0;
        const int tet_order = 10;

        Console.WriteLine("");
        Console.WriteLine("TEST004");
        Console.WriteLine("  TET_MESH_NODE_ORDER determines the order of ");
        Console.WriteLine("  each node in a tet mesh.");
        Console.WriteLine("");
        Console.WriteLine("  The order of a node is the number of tetrahedrons");
        Console.WriteLine("  that use the node as part of their definition.");

        TetMesh.tet_mesh_order10_example_size(ref node_num, ref tet_num);

        Console.WriteLine("");
        Console.WriteLine("  This mesh has tetrahedron order " + tet_order + "");
        Console.WriteLine("  The number of tetrahedrons is   " + tet_num + "");

        double[] node_xyz = new double[3 * node_num];
        int[] tet_node = new int[tet_order * tet_num];

        TetMesh.tet_mesh_order10_example_set(node_num, tet_num,
            ref node_xyz, ref tet_node);

        typeMethods.i4mat_transpose_print(tet_order, tet_num, tet_node,
            "  The tet mesh:");

        int[] node_order = TetMesh.tet_mesh_node_order(tet_order, tet_num, tet_node, node_num);

        typeMethods.i4vec_print(node_num, node_order, "  Node orders:");

        Console.WriteLine("");
        Console.WriteLine("  Check that the following are equal:");
        Console.WriteLine("");
        Console.WriteLine("  Number of tetrahedrons * order = " + tet_num * tet_order + "");
        Console.WriteLine("  Sum of node orders             = " + typeMethods.i4vec_sum(node_num, node_order) + "");
    }

    private static void test005()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST005 tests TETRAHEDRON_BARYCENTRIC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int test1;
        const int test1_num = 3;
        const int test2_num = 5;

        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("TEST005");
        Console.WriteLine("  TETRAHEDRON_BARYCENTRIC computes the barycentric");
        Console.WriteLine("  coordinates of a point.");
        //
        //  Choose a random tetrahedron.
        //
        for (test1 = 1; test1 <= test1_num; test1++)
        {
            double[] tet_xyz = UniformRNG.r8mat_uniform_01_new(3, 4, ref seed);

            typeMethods.r8mat_transpose_print(3, 4, tet_xyz, "  Random tetrahedron:");
            //
            //  Choose barycentric coordinates C1 at random.
            //
            //  Define a point P.
            //
            //  Have TETRAHEDRON_BARYCENTRIC compute C2, the barycentric coordinates of P.
            //
            int test2;
            for (test2 = 1; test2 <= test2_num; test2++)
            {
                double[] c1 = UniformRNG.r8vec_uniform_01_new(4, ref seed);
                double c1_sum = typeMethods.r8vec_sum(4, c1);
                int i;
                for (i = 0; i < 4; i++)
                {
                    c1[i] /= c1_sum;
                }

                double[] p = typeMethods.r8mat_mv_new(3, 4, tet_xyz, c1);

                double[] c2 = Tetrahedron.tetrahedron_barycentric(tet_xyz, p);

                Console.WriteLine("");
                string cout = "  C1 = ";
                for (i = 0; i < 4; i++)
                {
                    cout += "  " + c1[i].ToString("0.######").PadLeft(14);
                }

                Console.WriteLine(cout);
                cout = "  C2 = ";
                for (i = 0; i < 4; i++)
                {
                    cout += "  " + c2[i].ToString("0.######").PadLeft(14);
                }

                Console.WriteLine(cout);
            }
        }
    }

    private static void test006()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST006 tests TET_MESH_TET_NEIGHBORS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int node_num = 0;
        int tet_num = 0;
        const int tet_order = 4;

        Console.WriteLine("");
        Console.WriteLine("TEST006");
        Console.WriteLine("  TET_MESH_TET_NEIGHBORS computes the 4 neighboring");
        Console.WriteLine("  tetrahedrons of each tetrahedron in a tet mesh.");
        Console.WriteLine("  containing a point.");
        //
        //  Set up the example tetrahedron mesh.
        //
        TetMesh.tet_mesh_order4_example_size(ref node_num, ref tet_num);

        Console.WriteLine("");
        Console.WriteLine("  This mesh has tetrahedron order " + tet_order + "");
        Console.WriteLine("  The number of tetrahedrons is   " + tet_num + "");

        double[] node_xyz = new double[3 * node_num];
        int[] tet_node = new int[tet_order * tet_num];

        TetMesh.tet_mesh_order4_example_set(node_num, tet_num, ref node_xyz, ref tet_node);
        //
        //  Print the tets.
        //
        typeMethods.i4mat_transpose_print_some(tet_order, tet_num, tet_node,
            1, 1, tet_order, 10, "  First 10 Tets:");
        //
        //  The TET_NEIGHBOR array is needed by TET_MESH_DELAUNAY_SEARCH.
        //
        int[] tet_neighbor = TetMesh.tet_mesh_neighbor_tets(tet_order, tet_num, tet_node);

        typeMethods.i4mat_transpose_print_some(4, tet_num, tet_neighbor,
            1, 1, 4, 10, "  First 10 Tet Neighbors:");

    }

    private static void test007()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST007 tests TET_MESH_SEARCH_NAIVE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int face = 0;
        int node_num = 0;
        double[] p = new double[3];
        int step_num = 0;
        int test;
        const int test_num = 5;
        int tet_num = 0;
        const int tet_order = 4;
        double[] tet_xyz = new double[3 * 4];

        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("TEST007");
        Console.WriteLine("  TET_MESH_SEARCH_NAIVE uses a naive algorithm");
        Console.WriteLine("  to search a tetrahedral mesh for the tetrahedron");
        Console.WriteLine("  containing a point.");
        //
        //  Set up the example tetrahedron mesh.
        //
        TetMesh.tet_mesh_order4_example_size(ref node_num, ref tet_num);

        Console.WriteLine("");
        Console.WriteLine("  This mesh has tetrahedron order " + tet_order + "");
        Console.WriteLine("  The number of tetrahedrons is   " + tet_num + "");

        double[] node_xyz = new double[3 * node_num];
        int[] tet_node = new int[tet_order * tet_num];

        TetMesh.tet_mesh_order4_example_set(node_num, tet_num, ref node_xyz, ref tet_node);
        //
        //  The TET_NEIGHBOR array is needed for the Delaunay search.
        //
        int[] tet_neighbor = TetMesh.tet_mesh_neighbor_tets(tet_order, tet_num, tet_node);

        for (test = 1; test <= test_num; test++)
        {
            //
            //  Choose a tetrahedron at random.
            //
            int tet1 = UniformRNG.i4_uniform_ab(0, tet_num - 1, ref seed);

            Console.WriteLine("");
            Console.WriteLine("  Point was chosen from tetrahedron    " + tet1.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");

            int j;
            for (j = 0; j < 4; j++)
            {
                int k = tet_node[j + tet1 * 4];
                int i;
                for (i = 0; i < 3; i++)
                {
                    tet_xyz[i + j * 3] = node_xyz[i + k * 3];
                }
            }

            //
            //  Choose a point in the tetrahedron at random.
            //
            Tetrahedron.tetrahedron_sample(tet_xyz, 1, ref seed, ref p);
            //
            //  Naive search.
            //
            int tet2 = TetMesh.tet_mesh_search_naive(node_num, node_xyz, tet_order, tet_num,
                tet_node, p, ref step_num);

            Console.WriteLine("  Naive search ended in tetrahedron    " + tet2.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                                                        + ", number of steps = " + step_num + "");
            //
            //  Delaunay search.
            //
            int tet3 = TetMesh.tet_mesh_search_delaunay(node_num, node_xyz, tet_order,
                tet_num, tet_node, tet_neighbor, p, ref face, ref step_num);

            Console.WriteLine("  Delaunay search ended in tetrahedron " + tet3.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                                                        + ", number of steps = " + step_num + "");
        }
    }
}