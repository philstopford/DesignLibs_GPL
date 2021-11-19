using System;
using Burkardt;
using Burkardt.FEM;
using Burkardt.Graph;
using Burkardt.MatrixNS;
using Burkardt.NavierStokesNS;
using Burkardt.Types;
using Burkardt.Uniform;
using Triangle = Burkardt.TriangleNS.Triangle;

namespace FEM2DPackTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for FEM2D_PACK_TEST.
        //
        //  Discussion:
        //
        //    FEM2D_PACK_TEST tests the FEM2D_PACK library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 January 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("FEM2D_PACK_TEST:");
        Console.WriteLine("  Test the FEM2D_PACK library.");

        test01();
        test02();
        test03();
        test04();
        test05();
        test07();
        test08();
        test09();

        test10();
        test105();
        test11();
        test12();
        test13();
        test135();
        test14();
        test15();
        test16();
        test18();
        test19();

        test20();
        test21();
        test22();
        test23();
        test24();

        Console.WriteLine("");
        Console.WriteLine("FEM2D_PACK_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests BANDWIDTH_MESH.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] element_node;
        int element_num;
        int element_order;
        int m = 0;
        int ml = 0;
        int mu = 0;
        int nelemx;
        int nelemy;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  BANDWIDTH_MESH computes the geometric bandwidth:");
        Console.WriteLine("  of a finite element mesh.");

        nelemx = 2;
        nelemy = 6;

        Console.WriteLine("");
        Console.WriteLine("  NELEMX = " + nelemx + "");
        Console.WriteLine("  NELEMY = " + nelemy + "");

        element_order = 6;
        element_num = Grid.grid_element_num("T6", nelemx, nelemy);

        Console.WriteLine("");
        Console.WriteLine("  ELEMENT_ORDER = " + element_order + "");
        Console.WriteLine("  ELEMENT_NUM   = " + element_num + "");

        element_node = Grid.grid_t6_element(nelemx, nelemy);

        Grid.grid_print(element_order, element_num, element_node);

        Bandwidth.bandwidth_mesh(element_order, element_num, element_node, ref ml, ref mu, ref m);

        Console.WriteLine("");
        Console.WriteLine("  Lower bandwidth ML = " + ml + "");
        Console.WriteLine("  Upper bandwidth MU = " + mu + "");
        Console.WriteLine("  Total bandwidth M  = " + m + "");
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests BANDWIDTH_VAR, NS_T6_VAR_COUNT, NS_T6_VAR_SET.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] element_node;
        int element_num = 0;
        int element_order = 0;
        int i;
        int ihi;
        int ilo;
        int m = 0;
        int ml = 0;
        int mu = 0;
        int nelemx;
        int nelemy;
        int node;
        int node_num;
        int[] var;
        int[] var_node;
        int var_num;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  For the Navier Stokes variables associated with");
        Console.WriteLine("  a mesh of T6 elements,");
        Console.WriteLine("  NS_T6_VAR_COUNT counts variables,");
        Console.WriteLine("  NS_T6_VAR_SET sets them,");
        Console.WriteLine("  BANDWIDTH_VAR computes the variable bandwidth.");

        nelemx = 2;
        nelemy = 6;

        Console.WriteLine("");
        Console.WriteLine("  NELEMX = " + nelemx + "");
        Console.WriteLine("  NELEMY = " + nelemy + "");

        element_order = 6;
        element_num = Grid.grid_element_num("T6", nelemx, nelemy);
        node_num = Grid.grid_node_num("T6", nelemx, nelemy);

        Console.WriteLine("");
        Console.WriteLine("  ELEMENT_ORDER = " + element_order + "");
        Console.WriteLine("  ELEMENT_NUM   = " + element_num + "");
        Console.WriteLine("  NODE_NUM      = " + node_num + "");

        element_node = Grid.grid_t6_element(nelemx, nelemy);

        Grid.grid_print(element_order, element_num, element_node);

        var_node = new int[node_num + 1];

        var_num = NavierStokes.ns_t6_var_count(element_num, element_node, node_num, ref var_node);

        Console.WriteLine("");
        Console.WriteLine("  Number of variables VAR_NUM = " + var_num + "");

        typeMethods.i4vec_print(node_num + 1, var_node, "  VAR_NODE pointer vector:");

        var = NavierStokes.ns_t6_var_set(element_num, element_node, node_num, var_node,
            var_num);

        Console.WriteLine("");
        Console.WriteLine("  Node    U_Var    V_Var    P_Var");
        Console.WriteLine("");

        for (node = 0; node < node_num; node++)
        {
            ilo = var_node[node];
            ihi = var_node[node + 1] - 1;

            string cout = "  " + node.ToString().PadLeft(8);
            cout += "  ";

            for (i = ilo; i <= ihi; i++)
            {
                cout += "  " + var[i - 1].ToString().PadLeft(8);
            }

            Console.WriteLine(cout);
        }

        Bandwidth.bandwidth_var(element_order, element_num, element_node,
            node_num, var_node, var_num, var, ref ml, ref mu, ref m);

        Console.WriteLine("");
        Console.WriteLine("  Lower bandwidth ML = " + ml + "");
        Console.WriteLine("  Upper bandwidth MU = " + mu + "");
        Console.WriteLine("  Total bandwidth M  = " + m + "");
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests BASIS_11_**_TEST.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 March 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  BASIS_11_T3_TEST - Test the T3 basis functions.");
        Console.WriteLine("  BASIS_11_T4_TEST - Test the T4 basis functions.");
        Console.WriteLine("  BASIS_11_T6_TEST - Test the T6 basis functions.");

        Basis11.basis_11_t3_test();

        Basis11.basis_11_t4_test();

        Basis11.basis_11_t6_test();
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests BASIS_MN_**_TEST.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 February 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  Test the computation of basis functions by evaluating them");
        Console.WriteLine("  at the nodes that define the basis functions.");
        Console.WriteLine("");
        Console.WriteLine("  BASIS_MN_Q4_TEST - for the Q4 element.");
        Console.WriteLine("  BASIS_MN_T3_TEST - for the T3 element.");
        Console.WriteLine("  BASIS_MN_T4_TEST - for the T4 element.");
        Console.WriteLine("  BASIS_MN_T6_TEST - for the T6 element.");

        Basis_mn.basis_mn_q4_test();

        Basis_mn.basis_mn_t3_test();

        Basis_mn.basis_mn_t4_test();

        Basis_mn.basis_mn_t6_test();
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 demonstrates DERIVATIVE_AVERAGE_T3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 June 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int NELEMX = 7;
        int NELEMY = 5;

        double angle;
        double[] c;
        int col;
        double[] dcdx;
        double dcdx_exact;
        double[] dcdy;
        double dcdy_exact;
        int[] element_node;
        int element_num;
        int node;
        int node_num;
        double[] node_xy;
        double r;
        int row;
        double x;
        double y;

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  DERIVATIVE_AVERAGE_T3 averages the spatial derivatives");
        Console.WriteLine("  of a finite element function at the nodes.");
        //
        //  How many elements are there?
        //
        element_num = Grid.grid_t3_element_num(NELEMX, NELEMY);
        //
        //  How many nodes are there?
        //
        node_num = Grid.grid_t3_node_num(NELEMX, NELEMY);

        c = new double[node_num];
        dcdx = new double[node_num];
        dcdy = new double[node_num];
        node_xy = new double[2 * node_num];
        //
        //  Get the nodes that make up each element.
        //
        element_node = Grid.grid_t3_element(NELEMX, NELEMY);
        //
        //  Generate the coordinates of the nodes.
        //
        node = 0;

        for (row = 0; row <= NELEMY; row++)
        {
            r = ((NELEMY - row) * 1.0
                 + +row * 3.0)
                / NELEMY;

            for (col = 0; col <= NELEMX; col++)
            {
                angle = ((NELEMX - col) * 135.0
                         + +col * 45.0)
                        / NELEMX;

                angle = Helpers.degrees_to_radians(angle);

                node_xy[0 + node * 2] = r * Math.Cos(angle);
                node_xy[1 + node * 2] = r * Math.Sin(angle);

                node += 1;
            }
        }

        //
        //  Set the finite element function.
        //
        for (node = 0; node < node_num; node++)
        {
            x = node_xy[0 + node * 2];
            y = node_xy[1 + node * 2];
            c[node] = Math.Sin(x) * (1.0 + y * y);
        }

        Derivative.derivative_average_t3(node_num, node_xy, element_num,
            element_node, c, dcdx, dcdy);

        Console.WriteLine("");
        Console.WriteLine("  C         X               Y");
        Console.WriteLine("         dCdX(computed)  dCdY(computed)");
        Console.WriteLine("         dCdX(exact)     dCdY(exact)");
        Console.WriteLine("");

        for (node = 0; node < node_num; node++)
        {
            x = node_xy[0 + node * 2];
            y = node_xy[1 + node * 2];

            dcdx_exact = Math.Cos(x) * (1.0 * y * y);
            dcdy_exact = Math.Sin(x) * 2.0 * y;

            Console.WriteLine("");
            Console.WriteLine("  " + c[node].ToString().PadLeft(14)
                                   + "  " + node_xy[0 + node * 2].ToString().PadLeft(14)
                                   + "  " + node_xy[1 + node * 2].ToString().PadLeft(14) + node_xy[1 + node * 2] + "");
            Console.WriteLine("  " + "              "
                                   + "  " + dcdx[node].ToString().PadLeft(14)
                                   + "  " + dcdy[node].ToString().PadLeft(14) + "");
            Console.WriteLine("  " + "              "
                                   + "  " + dcdx_exact.ToString().PadLeft(14)
                                   + "  " + dcdy_exact.ToString().PadLeft(14) + dcdy_exact + "");
        }
    }

    private static void test07()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 tests ELEMENT_EPS using Q4 elements.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 February 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int NELEMX = 7;
        int NELEMY = 5;

        int ELEMENT_NUM = NELEMX* NELEMY;
        int NODE_NUM = ( NELEMX + 1 ) *(NELEMY + 1);

        double angle;
        string code = "Q4";
        int col;
        int element;
        bool[] element_mask = new bool[ELEMENT_NUM];
        int[] element_node;
        int element_show = 2;
        string file_name = "fem2d_pack_prb_q4.eps";
        int node;
        int node_show = 2;
        double[] node_xy = new double[2 * NODE_NUM];
        double r;
        int row;

        Console.WriteLine("");
        Console.WriteLine("TEST07");
        Console.WriteLine("  ELEMENTS_EPS creates an Encapsulated PostScript");
        Console.WriteLine("  file containing an image of a mesh.");

        element_node = Grid.grid_q4_element(NELEMX, NELEMY);

        node = 0;

        for (row = 0; row <= NELEMY; row++)
        {
            r = ((NELEMY - row) * 1.0
                 + +row * 3.0)
                / NELEMY;

            for (col = 0; col <= NELEMX; col++)
            {
                angle = ((NELEMX - col) * 135.0
                         + +col * 45.0)
                        / NELEMX;

                angle = Helpers.degrees_to_radians(angle);

                node_xy[0 + node * 2] = r * Math.Cos(angle);
                node_xy[1 + node * 2] = r * Math.Sin(angle);

                node += 1;
            }

        }

        for (element = 0; element < ELEMENT_NUM; element++)
        {
            element_mask[element] = true;
        }

        Element.elements_eps(file_name, NODE_NUM, node_xy, code,
            ELEMENT_NUM, element_mask, element_node, node_show, element_show);
    }

    private static void test08()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 tests ELEMENT_EPS, using T3 elements.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int NELEMX = 7;
        int NELEMY = 5;

        int ELEMENT_NUM = 2 * NELEMX * NELEMY;
        int NODE_NUM = ( NELEMX + 1 ) *(NELEMY + 1);

        double angle;
        string code = "T3";
        int col;
        int element;
        bool[] element_mask = new bool[ELEMENT_NUM];
        int[] element_node;
        int element_show = 2;
        string file_name = "fem2d_pack_prb_t3.eps";
        int node;
        int node_show = 2;
        double[] node_xy = new double[2 * NODE_NUM];
        double r;
        int row;

        Console.WriteLine("");
        Console.WriteLine("TEST08");
        Console.WriteLine("  ELEMENTS_EPS creates an Encapsulated PostScript");
        Console.WriteLine("  file containing an image of a mesh.");

        element_node = Grid.grid_t3_element(NELEMX, NELEMY);

        node = 0;

        for (row = 0; row <= NELEMY; row++)
        {
            r = ((NELEMY - row) * 1.0
                 + +row * 3.0)
                / NELEMY;

            for (col = 0; col <= NELEMX; col++)
            {
                angle = ((NELEMX - col) * 135.0
                         + +col * 45.0)
                        / NELEMX;

                angle = Helpers.degrees_to_radians(angle);

                node_xy[0 + node * 2] = r * Math.Cos(angle);
                node_xy[1 + node * 2] = r * Math.Sin(angle);

                node += 1;
            }
        }

        for (element = 0; element < ELEMENT_NUM; element++)
        {
            element_mask[element] = true;
        }

        Element.elements_eps(file_name, NODE_NUM, node_xy, code,
            ELEMENT_NUM, element_mask, element_node, node_show, element_show);
    }

    private static void test09()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST09 tests ELEMENT_EPS, using T4 elements.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 March 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int NELEMX = 7;
        int NELEMY = 5;

        string code = "T4";
        int col;
        int element;
        bool[] element_mask;
        int[] element_node;
        int element_num;
        int element_show = 1;
        string file_name = "fem2d_pack_prb_t4.eps";
        int i;
        int node;
        int node_num;
        int node_show = 2;
        double[] node_xy;
        int row;
        double x;
        double y;

        Console.WriteLine("");
        Console.WriteLine("TEST09");
        Console.WriteLine("  ELEMENTS_EPS creates an Encapsulated PostScript");
        Console.WriteLine("  file containing an image of a T4 mesh.");
        //
        //  How many elements are there?
        //
        element_num = Grid.grid_t4_element_num(NELEMX, NELEMY);

        element_mask = new bool[element_num];
        //
        //  How many nodes are there?
        //
        node_num = Grid.grid_t4_node_num(NELEMX, NELEMY);

        node_xy = new double[2 * node_num];
        //
        //  Get the nodes that make up each element.
        //
        element_node = Grid.grid_t4_element(NELEMX, NELEMY);

        node = 0;

        for (row = 0; row <= NELEMY; row++)
        {
            y = ((3 * NELEMY - row) * 0.0
                 + +row * 6.0)
                / (3 * NELEMY);

            for (col = 0; col <= NELEMX; col++)
            {
                x = ((2 * NELEMX - col) * 0.0
                     + +col * 6.0)
                    / (2 * NELEMX);

                node_xy[0 + node * 2] = x;
                node_xy[1 + node * 2] = y;

                node += 1;
            }

            //
            //  Skip over the two rows of interior nodes.
            //
            node += NELEMX;
            node += NELEMX;
        }

        //
        //  The coordinates of interior nodes are the average of the vertices.
        //
        for (element = 0; element < element_num; element++)
        {
            x = 0.0;
            y = 0.0;
            for (i = 0; i < 3; i++)
            {
                node = element_node[i + element * 4];
                x += node_xy[0 + (node - 1) * 2];
                y += node_xy[1 + (node - 1) * 2];
            }

            node = element_node[3 + element * 4];
            node_xy[0 + (node - 1) * 2] = x / 3.0;
            node_xy[1 + (node - 1) * 2] = y / 3.0;
        }

        for (element = 0; element < element_num; element++)
        {
            element_mask[element] = true;
        }

        Element.elements_eps(file_name, node_num, node_xy, code,
            element_num, element_mask, element_node, node_show, element_show);
    }

    private static void test10()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST10 tests ELEMENT_EPS, using T6 elements.
        //
        //  Discussion:
        //
        //    We generate a big grid of T6 elements, but we only want to
        //    look at the six elements shared by node 85.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 February 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int NELEMX = 6;
        int NELEMY = 6;

        int ELEMENT_NUM = 2 * NELEMX * NELEMY;
        int NODE_NUM = (2 * NELEMX + 1) * (2 * NELEMY + 1);

        string code = "T6";
        int col;
        int element;
        bool[] element_mask = new bool[ELEMENT_NUM];
        int[] element_node;
        int element_show = 2;
        string file_name = "fem2d_pack_prb_t6.eps";
        int node;
        int node_show = 2;
        double[] node_xy = new double[2 * NODE_NUM];
        double x;
        double y;
        int row;

        Console.WriteLine("");
        Console.WriteLine("TEST10");
        Console.WriteLine("  ELEMENTS_EPS creates an Encapsulated PostScript");
        Console.WriteLine("  file containing an image of a mesh.");

        element_node = Grid.grid_t6_element(NELEMX, NELEMY);

        node = 0;

        for (row = 0; row <= 2 * NELEMY; row++)
        {
            y = ((2 * NELEMY - row) * 0.0
                 + +row * 6.0)
                / (2 * NELEMY);

            for (col = 0; col <= 2 * NELEMX; col++)
            {
                x = ((2 * NELEMX - col) * 0.0
                     + +col * 6.0)
                    / (2 * NELEMX);

                node_xy[0 + node * 2] = x;
                node_xy[1 + node * 2] = y;

                node += 1;
            }

        }

        for (element = 0; element < ELEMENT_NUM; element++)
        {
            element_mask[element] = false;
        }

        element_mask[30 - 1] = true;
        element_mask[31 - 1] = true;
        element_mask[32 - 1] = true;
        element_mask[41 - 1] = true;
        element_mask[42 - 1] = true;
        element_mask[43 - 1] = true;

        Element.elements_eps(file_name, NODE_NUM, node_xy, code,
            ELEMENT_NUM, element_mask, element_node, node_show, element_show);
    }

    private static void test105()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST105 tests GRID_NODES_01.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 May 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int node;
        int node_num;
        double[] node_xy;
        int num_x = 5;
        int num_y = 3;

        Console.WriteLine("");
        Console.WriteLine("TEST105");
        Console.WriteLine("  GRID_NODES_01 computes a regular grid in the unit square.");
        Console.WriteLine("");
        Console.WriteLine("  NUM_X =    " + num_x + "");
        Console.WriteLine("  NUM_Y =    " + num_y + "");
        node_num = num_x * num_y;
        Console.WriteLine("  NODE_NUM = " + node_num + "");
        Console.WriteLine("");

        node_xy = Grid.grid_nodes_01(num_x, num_y);

        for (node = 0; node < node_num; node++)
        {
            Console.WriteLine("  " + node.ToString().PadLeft(8)
                                   + "  " + node_xy[0 + node * 2].ToString().PadLeft(14)
                                   + "  " + node_xy[1 + node * 2].ToString().PadLeft(14) + "");
        }
    }

    private static void test11()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST11 tests GRID_TEST.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 February 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TEST11");
        Console.WriteLine("  GRID_TEST tests the grid routines.");

        Grid.grid_test("Q4");

        Grid.grid_test("Q8");

        Grid.grid_test("Q9");

        Grid.grid_test("Q12");

        Grid.grid_test("Q16");

        Grid.grid_test("QL");

        Grid.grid_test("T3");

        Grid.grid_test("T6");

        Grid.grid_test("T10");
    }

    private static void test12()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST12 tests INTERP_TEST.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TEST12");
        Console.WriteLine("  INTERP_TEST tests the interpolating Math.Power");
        Console.WriteLine("  of the element.");

        Interp.interp_test("Q4");

        Interp.interp_test("Q8");

        Interp.interp_test("Q9");

        Interp.interp_test("Q12");

        Interp.interp_test("Q16");

        Interp.interp_test("QL");

        Interp.interp_test("T3");

        Interp.interp_test("T4");

        Interp.interp_test("T6");

        Interp.interp_test("T10");
    }

    private static void test13()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST13 tests MAP_TEST.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 February 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TEST13");
        Console.WriteLine("  MAP_TEST tests the map routines.");

        Map.map_test("Q4");

        Map.map_test("Q8");

        Map.map_test("Q9");

        Map.map_test("Q12");

        Map.map_test("Q16");

        Map.map_test("QL");

        Map.map_test("T3");

        Map.map_test("T6");

        Map.map_test("T10");
    }

    private static void test135()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST135 tests MASS_MATRIX_T3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 January 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int ELEMENT_NUM = 8;
        int NODE_NUM = 9;

        double[] a;
        int[] element_node =  {
                1, 4, 2,
                5, 2, 4,
                4, 7, 5,
                8, 5, 7,
                2, 5, 3,
                6, 3, 5,
                5, 8, 6,
                9, 6, 8
            }
            ;
        double[] node_xy =  {
                0.0, 0.0,
                0.0, 0.5,
                0.0, 1.0,
                0.5, 0.0,
                0.5, 0.5,
                0.5, 1.0,
                1.0, 0.0,
                1.0, 0.5,
                1.0, 1.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TEST135");
        Console.WriteLine("  MASS_MATRIX_T3 computes the mass matrix for");
        Console.WriteLine("  a finite element system using T3 elements");
        Console.WriteLine("  (linear triangles).");

        a = Mass.mass_matrix_t3(NODE_NUM, ELEMENT_NUM, element_node, node_xy);

        typeMethods.r8mat_print(NODE_NUM, NODE_NUM, a, "  The T3 mass matrix:");
    }

    private static void test14()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST14 tests MASS_MATRIX_T6.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 February 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int ELEMENT_NUM = 2;
        int NODE_NUM = 9;

        double[] a;
        int[] element_node =  {
                1, 3, 7, 2, 5, 4,
                9, 7, 3, 8, 5, 6
            }
            ;
        double[] node_xy =  {
                0.0, 0.0,
                0.0, 0.5,
                0.0, 1.0,
                0.5, 0.0,
                0.5, 0.5,
                0.5, 1.0,
                1.0, 0.0,
                1.0, 0.5,
                1.0, 1.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TEST14");
        Console.WriteLine("  MASS_MATRIX_T6 computes the mass matrix for");
        Console.WriteLine("  a finite element system using T6 elements");
        Console.WriteLine("  (quadratic triangles).");

        a = Mass.mass_matrix_t6(NODE_NUM, ELEMENT_NUM, element_node, node_xy);

        typeMethods.r8mat_print(NODE_NUM, NODE_NUM, a, "  The T6 mass matrix:");

    }

    private static void test15()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST15 tests PHYSICAL_TO_REFERENCE_T3 and REFERENCE_TO_PHYSICAL_T3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 February 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 10;

        int i;
        int j;
        double[] phy = new double[2 * N];
        double[] ref_ = new double[ 2 * N];
        double[] ref2 = new double[2 * N];
        int seed;
        double[] t =  {
                1.0, 1.0,
                3.0, 1.0,
                2.0, 5.0
            }
            ;

        seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("TEST15");
        Console.WriteLine("  For an order 3 triangle,");
        Console.WriteLine("  PHYSICAL_TO_REFERENCE_T3 maps a physical point to");
        Console.WriteLine("    a reference point.");
        Console.WriteLine("  REFERENCE_TO_PHYSICAL_T3 maps a reference point to");
        Console.WriteLine("    a physical point.");
        Console.WriteLine("");
        Console.WriteLine("      XSI     ETA  ==>  X       Y    ==>  XSI2    ETA2");
        Console.WriteLine("");

        for (j = 0; j < N; j++)
        {
            for (i = 0; i < 2; i++)
            {
                ref_[i +j * 2] = UniformRNG.r8_uniform_01(ref seed);
            }

            switch (ref_[
                        0 + j * 2] + ref_[
                        1 + j * 2])
            {
                case > 1.0:
                    ref_[
                        0 + j * 2] = 1.0 - ref_[
                        0 + j * 2];
                    ref_[
                        1 + j * 2] = 1.0 - ref_[
                        1 + j * 2];
                    break;
            }
        }

        Reference.reference_to_physical_t3(t, N, ref_, ref phy);

        PhysicalToRef.physical_to_reference_t3(t, N, phy, ref ref2);

        for (j = 0; j < N; j++)
        {
            Console.WriteLine("  " + ref_[0 + j * 2].ToString().PadLeft(10)
                                   + "  " + ref_[1 + j * 2].ToString().PadLeft(10)
                                   + "  " + phy[0 + j * 2].ToString().PadLeft(10)  
                                   + "  " + phy[1 + j * 2].ToString().PadLeft(10) 
                                   + "  " + ref2[0 + j * 2].ToString().PadLeft(10)
                                   + "  " + ref2[1 + j * 2].ToString().PadLeft(10) + "");
        }
    }

    private static void test16()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST16 tests REFERENCE_TO_PHYSICAL_T6.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 February 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 16;

        int j;
        double[] phy = new double[2 * N];
        double[] ref_ = {
                0.00, 0.00,
                1.00, 0.00,
                0.00, 1.00,
                0.50, 0.00,
                0.50, 0.50,
                0.00, 0.50,
                0.25, 0.75,
                0.75, 0.25,
                0.40, 0.10,
                0.30, 0.20,
                0.20, 0.30,
                0.10, 0.40,
                0.10, 0.10,
                0.20, 0.20,
                0.30, 0.30,
                0.40, 0.40
            }
            ;
        double[] t =  {
                0.0, 0.0,
                2.0, 0.0,
                0.0, 4.0,
                1.0, 0.0,
                1.0, 1.0,
                0.0, 2.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TEST16");
        Console.WriteLine("  For an order 6 triangle,");
        Console.WriteLine("  REFERENCE_TO_PHYSICAL_T6 maps a reference point to");
        Console.WriteLine("    a physical point.");
        Console.WriteLine("");
        Console.WriteLine("      XSI     ETA  ==>  X       Y");
        Console.WriteLine("");

        Reference.reference_to_physical_t6(t, N, ref_, phy);

        for (j = 0; j < N; j++)
        {
            Console.WriteLine("  " + ref_[0 + j * 2].ToString().PadLeft(8) 
                                   + "  " + ref_[1 + j * 2].ToString().PadLeft(8)
                                   + "  " + phy[0 + j * 2].ToString().PadLeft(8)
                                   + "  " + phy[1 + j * 2].ToString().PadLeft(8) + "");
        }
    }

    private static void test18()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST18 tests the shape routines.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 February 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TEST18");
        Console.WriteLine("  SHAPE_TEST tests the shape routines.");

        Shape.shape_test("Q4");

        Shape.shape_test("Q8");

        Shape.shape_test("Q9");

        Shape.shape_test("Q12");

        Shape.shape_test("Q16");

        Shape.shape_test("QL");

        Shape.shape_test("T3");

        Shape.shape_test("T6");

        Shape.shape_test("T10");
    }

    private static void test19()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST19 tests SPHERE_GRID_Q4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int element;
        int[] element_node;
        int element_num;
        int element_order = 4;
        int nelemx = 8;
        int nelemy = 8;
        int node;
        int node_num;
        double[] node_xyz;
        int order;

        Console.WriteLine("");
        Console.WriteLine("TEST19");
        Console.WriteLine("  SPHERE_GRID_Q4_ELEMENT sets up a grid of");
        Console.WriteLine("    Q4 quadrilaterals on a sphere.");
        Console.WriteLine("  SPHERE_GRID_Q4_ELEMENT_NUM returns the number");
        Console.WriteLine("    of elements in the grid");
        Console.WriteLine("  SPHERE_GRID_Q4_NODE_NUM returns the number");
        Console.WriteLine("    of nodes in the grid.");
        Console.WriteLine("  SPHERE_GRID_Q4_NODE_XYZ returns the coordinates");
        Console.WriteLine("    of nodes in the grid.");

        element_num = Burkardt.SphereNS.Grid.sphere_grid_q4_element_num(nelemx, nelemy);
        node_num = Burkardt.SphereNS.Grid.sphere_grid_q4_node_num(nelemx, nelemy);

        Console.WriteLine("");
        Console.WriteLine("  Expected number of nodes =    " + node_num + "");
        Console.WriteLine("  Expected number of elements = " + element_num + "");

        element_node = Burkardt.SphereNS.Grid.sphere_grid_q4_element(nelemx, nelemy);

        Console.WriteLine("");
        Console.WriteLine("  The elements and their nodes:");
        Console.WriteLine("");

        for (element = 0; element < element_num; element++)
        {
            string cout = "  " + (element + 1).ToString().PadLeft(4) + "  ";
            for (order = 0; order < element_order; order++)
            {
                cout += "  " + element_node[order + element * element_order].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
        }

        node_xyz = Burkardt.SphereNS.Grid.sphere_grid_q4_node_xyz(nelemx, nelemy);

        Console.WriteLine("");
        Console.WriteLine("  The node coordinates:");
        Console.WriteLine("");

        for (node = 0; node < node_num; node++)
        {
            Console.WriteLine("  " + (node + 1).ToString().PadLeft(4)
                                   + "  " + node_xyz[0 + node * 3].ToString().PadLeft(12)
                                   + "  " + node_xyz[1 + node * 3].ToString().PadLeft(12)
                                   + "  " + node_xyz[2 + node * 3].ToString().PadLeft(12) + "");
        }

        //
        //  Write the elements and nodes to files.
        //
        typeMethods.r8mat_write("sphere_q4_nodes.txt", 3, node_num, node_xyz);

        typeMethods.i4mat_write("sphere_q4_elements.txt", element_order, element_num,
            element_node);
    }

    private static void test20()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST20 tests SPHERE_GRID_Q9.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int element;
        int[] element_node;
        int element_num;
        int element_order = 9;
        int nelemx = 3;
        int nelemy = 4;
        int node;
        int node_num;
        double[] node_xyz;
        int order;

        Console.WriteLine("");
        Console.WriteLine("TEST20");
        Console.WriteLine("  SPHERE_GRID_Q9_ELEMENT sets up a grid of");
        Console.WriteLine("    Q9 quadrilaterals on a sphere.");
        Console.WriteLine("  SPHERE_GRID_Q9_ELEMENT_NUM returns the number");
        Console.WriteLine("    of elements in the grid");
        Console.WriteLine("  SPHERE_GRID_Q9_NODE_NUM returns the number");
        Console.WriteLine("    of nodes in the grid.");
        Console.WriteLine("  SPHERE_GRID_Q9_NODE_XYZ returns the coordinates");
        Console.WriteLine("    of nodes in the grid.");

        element_num = Burkardt.SphereNS.Grid.sphere_grid_q9_element_num(nelemx, nelemy);
        node_num = Burkardt.SphereNS.Grid.sphere_grid_q9_node_num(nelemx, nelemy);

        Console.WriteLine("");
        Console.WriteLine("  Expected number of nodes =    " + node_num + "");
        Console.WriteLine("  Expected number of elements = " + element_num + "");

        element_node = Burkardt.SphereNS.Grid.sphere_grid_q9_element(nelemx, nelemy);

        Console.WriteLine("");
        Console.WriteLine("  The elements and their nodes:");
        Console.WriteLine("");

        for (element = 0; element < element_num; element++)
        {
            string cout = "  " + (element + 1).ToString().PadLeft(4) + "  ";
            for (order = 0; order < element_order; order++)
            {
                cout += "  " + element_node[order + element * element_order].ToString().PadLeft(4) ;
            }

            Console.WriteLine(cout);
        }

        node_xyz = Burkardt.SphereNS.Grid.sphere_grid_q9_node_xyz(nelemx, nelemy);

        Console.WriteLine("");
        Console.WriteLine("  The node coordinates:");
        Console.WriteLine("");

        for (node = 0; node < node_num; node++)
        {
            Console.WriteLine("  " + (node + 1).ToString().PadLeft(4)
                                   + "  " + node_xyz[0 + node * 3].ToString().PadLeft(12)
                                   + "  " + node_xyz[1 + node * 3].ToString().PadLeft(12)
                                   + "  " + node_xyz[2 + node * 3].ToString().PadLeft(12) + "");
        }

        //
        //  Write the elements and nodes to files.
        //
        typeMethods.r8mat_write("sphere_q9_nodes.txt", 3, node_num, node_xyz);

        typeMethods.i4mat_write("sphere_q9_elements.txt", element_order, element_num,
            element_node);
    }

    private static void test21()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST21 tests SPHERE_GRID_Q16.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int element;
        int[] element_node;
        int element_num;
        int element_order = 16;
        int i;
        int ilo;
        int j;
        int k;
        int nelemx = 2;
        int nelemy = 2;
        int node;
        int node_num;
        double[] node_xyz;

        Console.WriteLine("");
        Console.WriteLine("TEST21");
        Console.WriteLine("  SPHERE_GRID_Q16_ELEMENT sets up a grid of");
        Console.WriteLine("    Q16 quadrilaterals on a sphere.");
        Console.WriteLine("  SPHERE_GRID_Q16_ELEMENT_NUM returns the number");
        Console.WriteLine("    of elements in the grid");
        Console.WriteLine("  SPHERE_GRID_Q16_NODE_NUM returns the number");
        Console.WriteLine("    of nodes in the grid.");
        Console.WriteLine("  SPHERE_GRID_Q16_NODE_XYZ returns the coordinates");
        Console.WriteLine("    of nodes in the grid.");

        element_num = Burkardt.SphereNS.Grid.sphere_grid_q16_element_num(nelemx, nelemy);
        node_num = Burkardt.SphereNS.Grid.sphere_grid_q16_node_num(nelemx, nelemy);

        Console.WriteLine("");
        Console.WriteLine("  Expected number of nodes =    " + node_num + "");
        Console.WriteLine("  Expected number of elements = " + element_num + "");

        element_node = Burkardt.SphereNS.Grid.sphere_grid_q16_element(nelemx, nelemy);

        Console.WriteLine("");
        Console.WriteLine("  The elements and their nodes, listed in a way");
        Console.WriteLine("  that suggests their geometry:");
        Console.WriteLine("");

        element = element_num;

        for (j = 1; j <= nelemy; j++)
        {
            for (i = 1; i <= nelemx; i++)
            {
                element -= 1;
                Console.WriteLine("");

                    
                for (ilo = 12; 0 <= ilo; ilo -= 4)
                {
                    string cout = ilo switch
                    {
                        12 => "  " + (element + 1).ToString().PadLeft(4),
                        _ => "      "
                    };

                    for (k = ilo; k <= ilo + 3; k++)
                    {
                        cout += "  " + element_node[k + element * element_order].ToString().PadLeft(4);
                    }

                    Console.WriteLine(cout);
                }
            }
        }

        node_xyz = Burkardt.SphereNS.Grid.sphere_grid_q16_node_xyz(nelemx, nelemy);

        Console.WriteLine("");
        Console.WriteLine("  The node coordinates:");
        Console.WriteLine("");

        for (node = 0; node < node_num; node++)
        {
            Console.WriteLine("  " + (node + 1).ToString().PadLeft(4)
                                   + "  " + node_xyz[0 + node * 3].ToString().PadLeft(12) 
                                   + "  " + node_xyz[1 + node * 3].ToString().PadLeft(12)
                                   + "  " + node_xyz[2 + node * 3].ToString().PadLeft(12) + "");
        }

        //
        //  Write the elements and nodes to files.
        //
        typeMethods.r8mat_write("sphere_q16_nodes.txt", 3, node_num, node_xyz);

        typeMethods.i4mat_write("sphere_q16_elements.txt", element_order, element_num,
            element_node);

    }

    private static void test22()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST22 tests SPHERE_GRID_T3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int dim;
        int dim_num = 3;
        int element;
        string element_file_name = "sphere_t3_elements.txt";
        int[] element_node;
        int element_num;
        int element_order = 3;
        int nelemx = 8;
        int nelemy = 8;
        int node;
        int node_num;
        string node_file_name = "sphere_t3_nodes.txt";
        double[] node_xyz;
        int order;

        Console.WriteLine("");
        Console.WriteLine("TEST22");
        Console.WriteLine("  SPHERE_GRID_T3_ELEMENT sets up a grid of T3 triangles");
        Console.WriteLine("    on a sphere.");
        Console.WriteLine("  SPHERE_GRID_T3_ELEMENT_NUM returns the number");
        Console.WriteLine("    of elements in the grid");
        Console.WriteLine("  SPHERE_GRID_T3_NODE_NUM returns the number");
        Console.WriteLine("    of nodes in the grid.");
        Console.WriteLine("  SPHERE_GRID_T3_NODE_XYZ returns the coordinates");
        Console.WriteLine("    of nodes in the grid.");

        element_num = Burkardt.SphereNS.Grid.sphere_grid_t3_element_num(nelemx, nelemy);
        node_num = Burkardt.SphereNS.Grid.sphere_grid_t3_node_num(nelemx, nelemy);

        Console.WriteLine("");
        Console.WriteLine("  Expected number of nodes =    " + node_num + "");
        Console.WriteLine("  Expected number of elements = " + element_num + "");
        //
        //  Generate the ELEMENT_NODE array, print it, and write it to a file.
        //
        element_node = Burkardt.SphereNS.Grid.sphere_grid_t3_element(nelemx, nelemy);

        Console.WriteLine("");
        Console.WriteLine("  The elements and their nodes:");
        Console.WriteLine("");

        for (element = 0; element < element_num; element++)
        {
            string cout = "  " + (element + 1).ToString().PadLeft(4) + "  ";
            for (order = 0; order < element_order; order++)
            {
                cout += "  " + element_node[order + element * element_order].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
        }

        typeMethods.i4mat_write(element_file_name, element_order, element_num,
            element_node);
        //
        //  Generate the NODE_XYZ array, print it, and write it to a file.
        //
        node_xyz = Burkardt.SphereNS.Grid.sphere_grid_t3_node_xyz(nelemx, nelemy);

        Console.WriteLine("");
        Console.WriteLine("  The node coordinates:");
        Console.WriteLine("");

        for (node = 0; node < node_num; node++)
        {
            string cout1 = "  " + (node + 1).ToString().PadLeft(4);
            for (dim = 0; dim < dim_num; dim++)
            {
                cout1 += "  " + node_xyz[dim + node * 3].ToString().PadLeft(12);
            }

            Console.WriteLine(cout1);
        }

        //
        //  Write the elements and nodes to files.
        //
        typeMethods.r8mat_write(node_file_name, dim_num, node_num, node_xyz);
    }

    private static void test23()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST23 tests SPHERE_GRID_T6.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int element;
        int[] element_node;
        int element_num;
        int element_order = 6;
        int nelemx = 3;
        int nelemy = 4;
        int node;
        int node_num;
        double[] node_xyz;
        int order;

        Console.WriteLine("");
        Console.WriteLine("TEST23");
        Console.WriteLine("  SPHERE_GRID_T6_ELEMENT sets up a grid of T6 triangles");
        Console.WriteLine("    on a sphere.");
        Console.WriteLine("  SPHERE_GRID_T6_ELEMENT_NUM returns the number");
        Console.WriteLine("    of elements in the grid");
        Console.WriteLine("  SPHERE_GRID_T6_NODE_NUM returns the number");
        Console.WriteLine("    of nodes in the grid");
        Console.WriteLine("  SPHERE_GRID_T6_NODE_XYZ returns the coordinates");
        Console.WriteLine("    of nodes in the grid.");

        element_num = Burkardt.SphereNS.Grid.sphere_grid_t6_element_num(nelemx, nelemy);
        node_num = Burkardt.SphereNS.Grid.sphere_grid_t6_node_num(nelemx, nelemy);

        Console.WriteLine("");
        Console.WriteLine("  Expected number of nodes =    " + node_num + "");
        Console.WriteLine("  Expected number of elements = " + element_num + "");

        element_node = Burkardt.SphereNS.Grid.sphere_grid_t6_element(nelemx, nelemy);

        Console.WriteLine("");
        Console.WriteLine("  The elements and their nodes:");
        Console.WriteLine("");

        for (element = 0; element < element_num; element++)
        {
            string cout = "  " + (element + 1).ToString().PadLeft(4) + "  ";
            for (order = 0; order < element_order; order++)
            {
                cout += "  " + element_node[order + element * element_order].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
        }

        node_xyz = Burkardt.SphereNS.Grid.sphere_grid_t6_node_xyz(nelemx, nelemy);

        Console.WriteLine("");
        Console.WriteLine("  The node coordinates:");
        Console.WriteLine("");

        for (node = 0; node < node_num; node++)
        {
            Console.WriteLine("  " + (node + 1).ToString().PadLeft(4)
                                   + "  " + node_xyz[0 + node * 3].ToString().PadLeft(12)
                                   + "  " + node_xyz[1 + node * 3].ToString().PadLeft(12)
                                   + "  " + node_xyz[2 + node * 3].ToString().PadLeft(12) + "");
        }

        //
        //  Write the elements and nodes to files.
        //
        typeMethods.r8mat_write("sphere_t6_nodes.txt", 3, node_num, node_xyz);

        typeMethods.i4mat_write("sphere_t6_elements.txt", element_order, element_num,
            element_node);

    }

    private static void test24()

        //*****************************************************************************
        //
        //  Purpose:
        //
        //    TEST24 tests TRIANGLE_UNIT_SET.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 February 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int ORDER_MAX = 64;

        int a;
        int b;
        double coef;
        double err;
        double exact;
        int i;
        int order;
        double quad;
        int rule;
        int rule_max = 20;
        double value = 0;
        double[] weight = new double[ORDER_MAX];
        double[] xtab = new double[ORDER_MAX];
        double[] ytab = new double[ORDER_MAX];

        Console.WriteLine("");
        Console.WriteLine("TEST24");
        Console.WriteLine("  TRIANGLE_UNIT_SET sets up a quadrature");
        Console.WriteLine("  in the unit triangle,");
        Console.WriteLine("");

        for (a = 0; a <= 10; a++)
        {
            for (b = 0; b <= 10 - a; b++)
            {
                coef = (a + b + 2) * (double) (a + b + 1);
                for (i = 1; i <= b; i++)
                {
                    coef = coef * (a + i) / i;
                }

                Console.WriteLine("");
                Console.WriteLine("  A = " + a + "  B = " + b + "");
                Console.WriteLine("");
                Console.WriteLine("  Rule       QUAD           ERROR");
                Console.WriteLine("");

                for (rule = 1; rule <= rule_max; rule++)
                {
                    order = Triangle.triangle_unit_size(rule);

                    Triangle.triangle_unit_set(rule, xtab, ytab, weight);

                    quad = 0.0;

                    for (i = 0; i < order; i++)
                    {
                        switch (a)
                        {
                            case 0 when b == 0:
                                value = coef;
                                break;
                            case 0 when b != 0:
                                value = coef * Math.Pow(ytab[i], b);
                                break;
                            default:
                            {
                                if (a != 0 && b == 0)
                                {
                                    value = coef * Math.Pow(xtab[i], a);
                                }
                                else if (a != 0 && b != 0)
                                {
                                    value = coef * Math.Pow(xtab[i], a) * Math.Pow(ytab[i], b);
                                }

                                break;
                            }
                        }

                        quad += 0.5 * weight[i] * value;
                    }

                    exact = 1.0;
                    err = Math.Abs(exact - quad);

                    Console.WriteLine("  " + rule.ToString().PadLeft(4)
                                           + "  " + quad.ToString().PadLeft(14)
                                           + "  " + err.ToString().PadLeft(14) + "");
                }
            }
        }
    }
}