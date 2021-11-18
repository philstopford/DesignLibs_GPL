using System;
using System.Globalization;
using System.Linq;
using Burkardt.TetrahedronNS;
using Burkardt.Types;

namespace GeometryTest;

public static class TetrahedronTest
{
    public static void test203()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST203 tests TETRAHEDRON_CENTROID_3D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;

        double[] centroid;
        double[] tetra =
        {
            0.000000, 0.942809, -0.333333,
            -0.816496, -0.816496, -0.333333,
            0.816496, -0.816496, -0.333333,
            0.000000, 0.000000, 1.000000
        };

        Console.WriteLine("");
        Console.WriteLine("TEST203");
        Console.WriteLine("  For a tetrahedron in 3D,");
        Console.WriteLine("  TETRAHEDRON_CENTROID_3D computes the centroid;");

        typeMethods.r8mat_transpose_print(DIM_NUM, 4, tetra, "  Tetrahedron vertices:");

        centroid = Geometry.tetrahedron_centroid_3d(tetra);

        typeMethods.r8vec_print(DIM_NUM, centroid, "  Centroid:");

    }

    public static void test2031()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST2031 tests TETRAHEDRON_CONTAINS_POINT_3D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int TEST_NUM = 3;

        double[] c;
        double[] c_test =
        {
            0.0, 0.1, 0.2, 0.7,
            -1.3, 2.0, 0.2, 0.1,
            0.8, 0.6, -0.5, 0.1
        };
        int i;
        bool inside;
        int j;
        double[] p = new double[DIM_NUM];
        int test;
        double[] tetra =
        {
            0.000000, 0.942809, -0.333333,
            -0.816496, -0.816496, -0.333333,
            0.816496, -0.816496, -0.333333,
            0.000000, 0.000000, 1.000000
        };

        Console.WriteLine("");
        Console.WriteLine("TEST2031");
        Console.WriteLine("  For a tetrahedron in 3D,");
        Console.WriteLine("  TETRAHEDRON_CONTAINS_POINT_3D finds if a point ");
        Console.WriteLine("  is inside;");

        typeMethods.r8mat_transpose_print(DIM_NUM, 4, tetra, "  Tetrahedron vertices:");

        Console.WriteLine("");
        Console.WriteLine("  P     Inside_Tetra?");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            c = c_test.Skip(+test * 4).ToArray();

            typeMethods.r8vec_zero(DIM_NUM, ref p);
            for (i = 0; i < DIM_NUM; i++)
            {
                for (j = 0; j < 4; j++)
                {
                    p[i] += tetra[i + j * DIM_NUM] * c[j];
                }
            }

            inside = Geometry.tetrahedron_contains_point_3d(tetra, p);

            Console.WriteLine("  " + p[0].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + p[1].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + p[2].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + inside.ToString(CultureInfo.InvariantCulture).PadLeft(1) + "");
        }

    }

    public static void test2032()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST2032 tests TETRAHEDRON_CIRCUMSPHERE_3D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;

        double[] pc = new double[DIM_NUM];
        double r = 0;
        double[] tetra =
        {
            0.577350269189626, 0.0, 0.0,
            -0.288675134594813, 0.5, 0.0,
            -0.288675134594813, -0.5, 0.0,
            0.0, 0.0, 0.816496580927726
        };

        Console.WriteLine("");
        Console.WriteLine("TEST2032");
        Console.WriteLine("  For a tetrahedron in 3D,");
        Console.WriteLine("  TETRAHEDRON_CIRCUMSPHERE_3D computes the circumsphere;");

        typeMethods.r8mat_transpose_print(DIM_NUM, 4, tetra, "  Tetrahedron vertices:");

        Geometry.tetrahedron_circumsphere_3d(tetra, ref r, ref pc);

        typeMethods.r8vec_print(DIM_NUM, pc, "  Circumsphere center:");

        Console.WriteLine("");
        Console.WriteLine("  Circumsphere radius is " + r + "");

    }

    public static void test20321()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST20321 tests TETRAHEDRON_EDGE_LENGTH_3D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;

        double[] edge_length;
        double[] tetra =
        {
            0.577350269189626, 0.0, 0.0,
            -0.288675134594813, 0.5, 0.0,
            -0.288675134594813, -0.5, 0.0,
            0.0, 0.0, 0.816496580927726
        };

        Console.WriteLine("");
        Console.WriteLine("TEST20321");
        Console.WriteLine("  For a tetrahedron in 3D,");
        Console.WriteLine("  TETRAHEDRON_EDGE_LENGTH_3D computes the edge lengths;");

        typeMethods.r8mat_transpose_print(DIM_NUM, 4, tetra, "  Tetrahedron vertices:");

        edge_length = Geometry.tetrahedron_edge_length_3d(tetra);

        typeMethods.r8vec_print(6, edge_length, "  Edge lengths:");


    }

    public static void test20322()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST20322 tests TETRAHEDRON_INSPHERE_3D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;

        double[] pc = new double[DIM_NUM];
        double r = 0;
        double[] tetra =
        {
            0.577350269189626, 0.0, 0.0,
            -0.288675134594813, 0.5, 0.0,
            -0.288675134594813, -0.5, 0.0,
            0.0, 0.0, 0.816496580927726
        };

        Console.WriteLine("");
        Console.WriteLine("TEST20322");
        Console.WriteLine("  For a tetrahedron in 3D,");
        Console.WriteLine("  TETRAHEDRON_INSPHERE_3D computes the insphere;");

        typeMethods.r8mat_transpose_print(DIM_NUM, 4, tetra, "  Tetrahedron vertices:");

        Geometry.tetrahedron_insphere_3d(tetra, ref r, ref pc);

        typeMethods.r8vec_print(DIM_NUM, pc, "  Insphere center:");

        Console.WriteLine("");
        Console.WriteLine("  Insphere radius is " + r + "");

    }

    public static void tetrahedron_lattice_layer_point_next_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_LATTICE_LAYER_POINT_NEXT_TEST tests TETRAHEDRON_LATTICE_LAYER_POINT_NEXT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] c = new int[4];
        int i;
        int j;
        int layer;
        bool more;
        int n = 3;
        int[] v = new int[3];
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("TETRAHEDRON_LATTICE_LAYER_POINT_NEXT_TEST");
        Console.WriteLine("  TETRAHEDRON_LATTICE_LAYER_POINT_NEXT returns the next");
        Console.WriteLine("  point in a tetrahedron lattice layer defined by:");
        Console.WriteLine("");
        Console.WriteLine("    C[3] - 1 < X[0]/C[0] + X[1]/C[1] +X[2]/C[2] <= C[3].");

        c[0] = 2;
        c[1] = 3;
        c[2] = 4;
        v[0] = 0;
        v[1] = 0;
        v[2] = 0;

        Console.WriteLine("");
        Console.WriteLine("  N = " + n + "");
        cout = "  C =       ";
        for (i = 0; i < n; i++)
        {
            cout += "  " + c[i].ToString(CultureInfo.InvariantCulture).PadLeft(4);
        }

        Console.WriteLine(cout);

        for (layer = 0; layer <= 2; layer++)
        {
            Console.WriteLine("");
            Console.WriteLine("  Layer " + layer + "");
            Console.WriteLine("");

            c[3] = layer;
            more = false;
            i = 0;

            for (;;)
            {
                Geometry.tetrahedron_lattice_layer_point_next(c, ref v, ref more);
                if (!more)
                {
                    Console.WriteLine("  No more.");
                    break;
                }

                i += 1;
                cout = "  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4);
                for (j = 0; j < n; j++)
                {
                    cout += "  " + v[j].ToString(CultureInfo.InvariantCulture).PadLeft(4);
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static void test203225()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST203225 tests TETRAHEDRON_LATTICE_POINT_NEXT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 3;

        int[] c = new int[N + 1];
        int i;
        int j;
        bool more;
        int n = N;
        int[] v = new int[N];
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("TEST203225");
        Console.WriteLine("  TETRAHEDRON_LATTICE_POINT_NEXT returns the next lattice");
        Console.WriteLine("  point in a tetrahedron defined by:");
        Console.WriteLine("");
        Console.WriteLine("    0 <= X(1)/C(1) + X(2)/C(2) + X(3)/C(3) <= C(4).");

        for (i = 0; i < n + 1; i++)
        {
            c[i] = n + 1 - i;
        }

        for (i = 0; i < n; i++)
        {
            v[i] = 0;
        }

        more = false;

        Console.WriteLine("");
        Console.WriteLine("  N = " + n + "");
        cout = "  C = ";
        for (i = 0; i < n + 1; i++)
        {
            cout += "  " + c[i].ToString(CultureInfo.InvariantCulture).PadLeft(4);
        }

        Console.WriteLine(cout);
        Console.WriteLine("");

        i = 0;

        for (;;)
        {
            Geometry.tetrahedron_lattice_point_next(c, ref v, ref more);
            if (!more)
            {
                Console.WriteLine("  No more.");
                break;
            }

            i += 1;
            cout = "  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4);
            for (j = 0; j < n; j++)
            {
                cout += "  " + v[j].ToString(CultureInfo.InvariantCulture).PadLeft(4);
            }

            Console.WriteLine(cout);
        }
    }

    public static void test20323()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST20323 tests TETRAHEDRON_QUALITY1_3D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int TEST_NUM = 2;

        double quality;
        int test;
        double[] tetra;
        double[] tetra_test =
        {
            0.577350269189626, 0.0, 0.0,
            -0.288675134594813, 0.5, 0.0,
            -0.288675134594813, -0.5, 0.0,
            0.0, 0.0, 0.816496580927726,
            0.577350269189626, 0.0, 0.0,
            -0.288675134594813, 0.5, 0.0,
            -0.288675134594813, -0.5, 0.0,
            0.0, 0.0, 0.408248290463863
        };

        Console.WriteLine("");
        Console.WriteLine("TEST20323");
        Console.WriteLine("  For a tetrahedron in 3D,");
        Console.WriteLine("  TETRAHEDRON_QUALITY1_3D computes quality measure #1;");

        for (test = 0; test < TEST_NUM; test++)
        {
            tetra = tetra_test.Skip(+test * DIM_NUM * 4).ToArray();

            typeMethods.r8mat_transpose_print(DIM_NUM, 4, tetra, "  Tetrahedron vertices:");

            quality = Geometry.tetrahedron_quality1_3d(tetra);

            Console.WriteLine("");
            Console.WriteLine("  Tetrahedron quality is " + quality + "");
        }

    }

    public static void test203232()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST203232 tests TETRAHEDRON_QUALITY2_3D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int TEST_NUM = 2;

        double quality2;
        int test;
        double[] tetra;
        double[] tetra_test =
        {
            0.577350269189626, 0.0, 0.0,
            -0.288675134594813, 0.5, 0.0,
            -0.288675134594813, -0.5, 0.0,
            0.0, 0.0, 0.816496580927726,
            0.577350269189626, 0.0, 0.0,
            -0.288675134594813, 0.5, 0.0,
            -0.288675134594813, -0.5, 0.0,
            0.0, 0.0, 0.408248290463863
        };

        Console.WriteLine("");
        Console.WriteLine("TEST203232");
        Console.WriteLine("  For a tetrahedron in 3D,");
        Console.WriteLine("  TETRAHEDRON_QUALITY2_3D computes quality measure #2;");

        for (test = 0; test < TEST_NUM; test++)
        {
            tetra = tetra_test.Skip(+test * DIM_NUM * 4).ToArray();

            typeMethods.r8mat_transpose_print(DIM_NUM, 4, tetra, "  Tetrahedron vertices:");

            quality2 = Geometry.tetrahedron_quality2_3d(tetra);

            Console.WriteLine("");
            Console.WriteLine("  Tetrahedron quality is " + quality2 + "");
        }

    }

    public static void test203233()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST203233 tests TETRAHEDRON_QUALITY3_3D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int TEST_NUM = 2;

        double quality3;
        int test;
        double[] tetra;
        double[] tetra_test =
        {
            0.577350269189626, 0.0, 0.0,
            -0.288675134594813, 0.5, 0.0,
            -0.288675134594813, -0.5, 0.0,
            0.0, 0.0, 0.816496580927726,
            0.577350269189626, 0.0, 0.0,
            -0.288675134594813, 0.5, 0.0,
            -0.288675134594813, -0.5, 0.0,
            0.0, 0.0, 0.408248290463863
        };

        Console.WriteLine("");
        Console.WriteLine("TEST203233");
        Console.WriteLine("  For a tetrahedron in 3D,");
        Console.WriteLine("  TETRAHEDRON_QUALITY3_3D computes quality measure #3;");

        for (test = 0; test < TEST_NUM; test++)
        {
            tetra = tetra_test.Skip(+test * DIM_NUM * 4).ToArray();

            typeMethods.r8mat_transpose_print(DIM_NUM, 4, tetra, "  Tetrahedron vertices:");

            quality3 = Geometry.tetrahedron_quality3_3d(tetra);

            Console.WriteLine("");
            Console.WriteLine("  Tetrahedron quality is " + quality3 + "");
        }

    }

    public static void test203234()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST203234 tests TETRAHEDRON_QUALITY4_3D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int TEST_NUM = 2;

        double quality4;
        int test;
        double[] tetra;
        double[] tetra_test =
        {
            0.577350269189626, 0.0, 0.0,
            -0.288675134594813, 0.5, 0.0,
            -0.288675134594813, -0.5, 0.0,
            0.0, 0.0, 0.816496580927726,
            0.577350269189626, 0.0, 0.0,
            -0.288675134594813, 0.5, 0.0,
            -0.288675134594813, -0.5, 0.0,
            0.0, 0.0, 0.408248290463863
        };

        Console.WriteLine("");
        Console.WriteLine("TEST203234");
        Console.WriteLine("  For a tetrahedron in 3D,");
        Console.WriteLine("  TETRAHEDRON_QUALITY4_3D computes quality measure #4;");

        for (test = 0; test < TEST_NUM; test++)
        {
            tetra = tetra_test.Skip(+test * DIM_NUM * 4).ToArray();

            typeMethods.r8mat_transpose_print(DIM_NUM, 4, tetra, "  Tetrahedron vertices:");

            quality4 = Geometry.tetrahedron_quality4_3d(tetra);

            Console.WriteLine("");
            Console.WriteLine("  Tetrahedron quality is " + quality4 + "");
        }

    }

    public static void test203235()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST203235 tests TETRAHEDRON_RHOMBIC_SIZE_3D, TETRAHEDRON_RHOMBIC_SHAPE_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;

        int edge_num = 0;
        int face_num = 0;
        int[] face_order;
        int face_order_max = 0;
        int[] face_point;
        int point_num = 0;
        double[] point_coord;

        Console.WriteLine("");
        Console.WriteLine("TEST203235");
        Console.WriteLine("  For the rhombic tetrahedron,");
        Console.WriteLine("  TETRAHEDRON_RHOMBIC_SIZE_3D returns dimension information;");
        Console.WriteLine("  TETRAHEDRON_RHOMBIC_SHAPE_3D returns face and order information.");
        Console.WriteLine("  SHAPE_PRINT_3D prints this information.");
        //
        //  Get the data sizes.
        //
        Geometry.tetrahedron_rhombic_size_3d(ref point_num, ref edge_num, ref face_num,
            ref face_order_max);

        Console.WriteLine("");
        Console.WriteLine("    Number of vertices: " + point_num + "");
        Console.WriteLine("    Number of edges   : " + edge_num + "");
        Console.WriteLine("    Number of faces   : " + face_num + "");
        Console.WriteLine("    Maximum face order: " + face_order_max + "");
        //
        //  Make room for the data.
        //
        face_order = new int[face_num];
        face_point = new int[face_order_max * face_num];
        point_coord = new double[DIM_NUM * point_num];
        //
        //  Get the data.
        //
        Geometry.tetrahedron_rhombic_shape_3d(point_num, face_num, face_order_max,
            ref point_coord, ref face_order, ref face_point);
        //
        //  Print the data.
        //
        Burkardt.Geometry.Shape.shape_print_3d(point_num, face_num, face_order_max,
            point_coord, face_order, face_point);

    }

    public static void test20324()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST20324 tests TETRAHEDRON_SAMPLE_3D, TETRAHEDRON_BARYCENTRIC_3D.
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
        int DIM_NUM = 3;
        int TEST_NUM = 10;

        double[] p = new double[DIM_NUM];
        int seed = 123456789;
        double[] t =
        {
            1.0, 4.0, 3.0,
            2.0, 4.0, 3.0,
            1.0, 6.0, 3.0,
            1.0, 4.0, 4.0
        };
        int test;
        double[] xsi;

        Console.WriteLine("");
        Console.WriteLine("TEST20324");
        Console.WriteLine("  TETRAHEDRON_SAMPLE_3D samples a tetrahedron.");
        Console.WriteLine("  TETRAHEDRON_BARYCENTRIC_3D converts Cartesian to");
        Console.WriteLine("  barycentric coordinates.");
        Console.WriteLine("");
        Console.WriteLine("  We are computing the barycentric coordinates just to");
        Console.WriteLine("  verify that the points are inside the tetrahedron.");

        typeMethods.r8mat_transpose_print(DIM_NUM, DIM_NUM + 1, t, "  Tetrahedron vertices");

        Console.WriteLine("");
        Console.WriteLine("     P      Barycentric:");
        Console.WriteLine("");

        for (test = 1; test <= TEST_NUM; test++)
        {
            Geometry.tetrahedron_sample_3d(t, 1, ref seed, ref p);

            string cout = "  " + p[0].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                               + "  " + p[1].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                               + "  " + p[2].ToString(CultureInfo.InvariantCulture).PadLeft(8);

            xsi = Geometry.tetrahedron_barycentric_3d(t, p);

            Console.WriteLine(cout + "    "
                                   + "  " + xsi[0].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + xsi[1].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + xsi[2].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + xsi[3].ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");

        }

    }

    public static void test20325()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST20325 tests TETRAHEDRON_SIZE_3D, TETRAHEDRON_SHAPE_3D, SHAPE_PRINT_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;

        int edge_num = 0;
        int face_num = 0;
        int[] face_order;
        int face_order_max = 0;
        int[] face_point;
        int point_num = 0;
        double[] point_coord;

        Console.WriteLine("");
        Console.WriteLine("TEST20325");
        Console.WriteLine("  For the tetrahedron,");
        Console.WriteLine("  TETRAHEDRON_SIZE_3D returns dimension information;");
        Console.WriteLine("  TETRAHEDRON_SHAPE_3D returns face and order information.");
        Console.WriteLine("  SHAPE_PRINT_3D prints this information.");
        //
        //  Get the data sizes.
        //
        Geometry.tetrahedron_size_3d(ref point_num, ref edge_num, ref face_num, ref face_order_max);

        Console.WriteLine("");
        Console.WriteLine("    Number of vertices: " + point_num + "");
        Console.WriteLine("    Number of edges   : " + edge_num + "");
        Console.WriteLine("    Number of faces   : " + face_num + "");
        Console.WriteLine("    Maximum face order: " + face_order_max + "");
        //
        //  Make room for the data.
        //
        face_order = new int[face_num];
        face_point = new int[face_order_max * face_num];
        point_coord = new double[DIM_NUM * point_num];
        //
        //  Get the data.
        //
        Geometry.tetrahedron_shape_3d(point_num, face_num, face_order_max, ref point_coord,
            ref face_order, ref face_point);
        //
        //  Print the data.
        //
        Burkardt.Geometry.Shape.shape_print_3d(point_num, face_num, face_order_max,
            point_coord, face_order, face_point);


    }

    public static void tetrahedron_solid_angles_3d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    Geometry.tetrahedron_solid_angles_3d tests TETRAHEDRON_VOLUME_3D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] angle;
        double[] t1 =
        {
            0.000000, 0.942809, -0.333333,
            -0.816496, -0.816496, -0.333333,
            0.816496, -0.816496, -0.333333,
            0.000000, 0.000000, 1.000000
        };
        double[] t2 =
        {
            0.000000, 0.000000, 0.000000,
            1.000000, 0.000000, 0.000000,
            0.000000, 1.000000, 0.000000,
            0.000000, 0.000000, 1.000000
        };
        double[] t3 =
        {
            0.000000, 0.000000, 0.000000,
            1.000000, 0.000000, 0.000000,
            0.000000, 2.000000, 0.000000,
            0.000000, 0.000000, 4.000000
        };
        double[] t4 =
        {
            0.000000, 0.000000, 0.000000,
            1.000000, 0.000000, 0.000000,
            0.000000, 1.000000, 0.000000,
            1.000000, 1.000000, 1.000000
        };

        Console.WriteLine("");
        Console.WriteLine("Geometry.tetrahedron_solid_angles_3d_TEST");
        Console.WriteLine("  Geometry.tetrahedron_solid_angles_3d computes the solid angles");
        Console.WriteLine("  associated with the vertices of a tetrahedron in 3D.");

        typeMethods.r8mat_transpose_print(3, 4, t1, "  Tetrahedron #1");
        angle = Geometry.tetrahedron_solid_angles_3d(t1);
        typeMethods.r8vec_print(4, angle, "  Solid angles for tetrahedron #1");

        typeMethods.r8mat_transpose_print(3, 4, t2, "  Tetrahedron #2");
        angle = Geometry.tetrahedron_solid_angles_3d(t2);
        typeMethods.r8vec_print(4, angle, "  Solid angles for tetrahedron #2");

        typeMethods.r8mat_transpose_print(3, 4, t3, "  Tetrahedron #3");
        angle = Geometry.tetrahedron_solid_angles_3d(t3);
        typeMethods.r8vec_print(4, angle, "  Solid angles for tetrahedron #3");

        typeMethods.r8mat_transpose_print(3, 4, t4, "  Tetrahedron #4");
        angle = Geometry.tetrahedron_solid_angles_3d(t4);
        typeMethods.r8vec_print(4, angle, "  Solid angles for tetrahedron #4");

    }

    public static void test2033()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST2033 tests TETRAHEDRON_VOLUME_3D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;

        double[] tetra =
        {
            0.000000, 0.942809, -0.333333,
            -0.816496, -0.816496, -0.333333,
            0.816496, -0.816496, -0.333333,
            0.000000, 0.000000, 1.000000
        };
        double volume;

        Console.WriteLine("");
        Console.WriteLine("TEST2033");
        Console.WriteLine("  For a tetrahedron in 3D,");
        Console.WriteLine("  TETRAHEDRON_VOLUME_3D computes the volume;");

        typeMethods.r8mat_transpose_print(DIM_NUM, DIM_NUM + 1, tetra, "  Tetrahedron vertices");

        volume = Geometry.tetrahedron_volume_3d(tetra);

        Console.WriteLine("");
        Console.WriteLine("  Volume = " + volume + "");

    }

}