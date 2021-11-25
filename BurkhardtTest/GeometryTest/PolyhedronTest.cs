using System;
using System.Globalization;
using Burkardt.Polyhedron;
using Burkardt.Types;
using Burkardt.Uniform;

namespace GeometryTest;

public static class PolyhedronTest
{
    public static void polyhedron_area_3d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYHEDRON_AREA_3D_TEST tests POLYHEDRON_AREA_3D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int FACE_NUM = 4;
        const int ORDER_MAX = 3;
        const int NODE_NUM = 4;

        const double area_exact = 2.366025;
        double[] coord =
        {
            0.0, 0.0, 0.0,
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0
        };
        int i;
        int j;
        int[] node =
        {
            2, 1, 0,
            0, 1, 3,
            0, 3, 2,
            1, 2, 3
        };
        int[] order = {3, 3, 3, 3};
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("POLYHEDRON_AREA_3D_TEST");
        Console.WriteLine("  For a polyhedron in 3D:");
        Console.WriteLine("  POLYHEDRON_AREA_3D computes surface area;");
        Console.WriteLine("");
        Console.WriteLine("  Number of faces is " + FACE_NUM + "");

        typeMethods.i4vec_print(FACE_NUM, order, "  Order of each face:");

        Console.WriteLine("");
        Console.WriteLine("  Nodes per face:");
        Console.WriteLine("");
        for (i = 0; i < FACE_NUM; i++)
        {
            cout = "  " + i.ToString().PadLeft(4);
            for (j = 0; j < order[i]; j++)
            {
                cout += "  " + node[i * ORDER_MAX + j].ToString().PadLeft(10);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  Nodal coordinates:");
        Console.WriteLine("");
        for (j = 0; j < NODE_NUM; j++)
        {
            cout = "  " + j.ToString().PadLeft(4);
            for (i = 0; i < 3; i++)
            {
                cout += "  " + coord[i + j * 3].ToString(CultureInfo.InvariantCulture).PadLeft(10);
            }

            Console.WriteLine(cout);
        }

        double area = Geometry.polyhedron_area_3d(coord, ORDER_MAX, FACE_NUM, node, NODE_NUM,
            order);

        Console.WriteLine("");
        Console.WriteLine("  Surface area = " + area + "");
        Console.WriteLine("  Exact area =   " + area_exact + "");

    }

    public static void polyhedron_centroid_3d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYHEDRON_CENTROID_3D_TEST tests POLYHEDRON_CENTROID_3D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int DIM_NUM = 3;
        const int FACE_NUM = 4;
        const int ORDER_MAX = 3;
        const int NODE_NUM = 4;

        double[] centroid_exact = {0.25, 0.25, 0.25};
        double[] coord =
        {
            0.0, 0.0, 0.0,
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0
        };
        int i;
        int j;
        int[] node =
        {
            3, 2, 1,
            1, 2, 4,
            1, 4, 3,
            2, 3, 4
        };
        int[] order = {3, 3, 3, 3};
        string cout;

        Console.WriteLine("");
        Console.WriteLine("POLYHEDRON_CENTROID_3D_TEST");
        Console.WriteLine("  For a polyhedron in 3D:");
        Console.WriteLine("  POLYHEDRON_CENTROID_3D computes the centroid;");
        Console.WriteLine("");
        Console.WriteLine("  Number of faces is " + FACE_NUM + "");

        typeMethods.i4vec_print(FACE_NUM, order, "  Order of each face:");

        Console.WriteLine("");
        Console.WriteLine("  Nodes per face:");
        Console.WriteLine("");
        for (j = 0; j < FACE_NUM; j++)
        {
            cout = "  " + j.ToString().PadLeft(4);
            for (i = 0; i < order[j]; i++)
            {
                cout += "  " + node[i + j * ORDER_MAX].ToString().PadLeft(10);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  Nodal coordinates:");
        Console.WriteLine("");
        for (j = 0; j < NODE_NUM; j++)
        {
            cout = "  " + j.ToString().PadLeft(4);
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += "  " + coord[i + j * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(10);
            }

            Console.WriteLine(cout);
        }

        double[] centroid = Geometry.polyhedron_centroid_3d(coord, ORDER_MAX, FACE_NUM, node,
            NODE_NUM, order);

        Console.WriteLine("");
        Console.WriteLine("  The computed centroid is "
                          + centroid[0].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                          + centroid[1].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                          + centroid[2].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        Console.WriteLine("  The correct centroid is  "
                          + centroid_exact[0].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                          + centroid_exact[1].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                          + centroid_exact[2].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");

    }

    public static void test0825()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0825 tests POLYHEDRON_CONTAINS_POINT_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int DIM_NUM = 3;
        const int FACE_NUM = 4;
        const int FACE_ORDER_MAX = 3;
        const int NODE_NUM = 4;
        const int TEST_NUM = 100;

        int[] face_order = {3, 3, 3, 3};
        int[] face_point =
        {
            1, 2, 4,
            1, 3, 2,
            1, 4, 3,
            2, 3, 4
        };
        int test;
        double[] v =
        {
            0.0, 0.0, 0.0,
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0
        };

        Console.WriteLine("");
        Console.WriteLine("TEST0825");
        Console.WriteLine("  POLYHEDRON_CONTAINS_POINT_3D determines if a point is ");
        Console.WriteLine("  inside a polyhedron.");
        Console.WriteLine("");
        Console.WriteLine("  We test this routine by using a tetrahedron as the polyhedron.");
        Console.WriteLine("  For this shape, an independent check can be made,");
        Console.WriteLine("  using barycentric coordinates.");
        Console.WriteLine("");
        Console.WriteLine("  We label these checks IN1 and IN2, and we expect them to agree.");

        typeMethods.r8mat_transpose_print(DIM_NUM, NODE_NUM, v, "  The vertices:");

        typeMethods.i4vec_print(FACE_NUM, face_order, "  The face orders:");

        typeMethods.i4mat_transpose_print(FACE_ORDER_MAX, FACE_NUM, face_point,
            "  The nodes making each face:");

        Console.WriteLine("");
        Console.WriteLine("      X           Y           Z      IN1 IN2");
        Console.WriteLine("");

        int seed = 123456789;

        for (test = 1; test <= TEST_NUM; test++)
        {
            double[] p = UniformRNG.r8vec_uniform_01_new(DIM_NUM, ref seed);

            bool inside1 = Geometry.polyhedron_contains_point_3d(NODE_NUM, FACE_NUM,
                FACE_ORDER_MAX, v, face_order, face_point, p);

            double[] c = Burkardt.TetrahedronNS.Geometry.tetrahedron_barycentric_3d(v, p);

            bool inside2 = 0.0 <= c[0] && c[0] <= 1.0 &&
                           0.0 <= c[1] && c[1] <= 1.0 &&
                           0.0 <= c[2] && c[2] <= 1.0 &&
                           0.0 <= c[3] && c[3] <= 1.0 &&
                           c[0] + c[1] + c[2] + c[3] <= 1.0;

            Console.WriteLine("  " + p[0].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + p[1].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + p[2].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + inside1.ToString().PadLeft(1)
                                   + "  " + inside2.ToString().PadLeft(1) + "");

            if (inside1 != inside2)
            {
                Console.WriteLine("??? Disagreement!  Barycentric coordinates:");
                Console.WriteLine("  " + c[0].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + c[1].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + c[2].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + c[3].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
            }

        }
    }

    public static void test083()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST083 tests POLYHEDRON_VOLUME_3D and POLYHEDRON_VOLUME_3D_2;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int FACE_NUM = 4;
        int NODE_NUM = 4;
        int ORDER_MAX = 3;

        double[] coord =
        {
            0.0, 0.0, 0.0,
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0
        };
        int i;
        int j;
        int[] node =
        {
            2, 1, 0,
            0, 1, 3,
            0, 3, 2,
            1, 2, 3
        };
        int[] order = {3, 3, 3, 3};
        const double volume_exact = 1.0 / 6.0;
        string cout;

        Console.WriteLine("");
        Console.WriteLine("TEST083");
        Console.WriteLine("  For a polyhedron in 3D:");
        Console.WriteLine("  POLYHEDRON_VOLUME_3D computes the volume;");
        Console.WriteLine("  POLYHEDRON_VOLUME_3D_2 computes the volume;");
        Console.WriteLine("");
        Console.WriteLine("  Number of faces is " + FACE_NUM + "");

        typeMethods.i4vec_print(FACE_NUM, order, "  Order of each face:");

        Console.WriteLine("");
        Console.WriteLine("  Nodes per face:");
        Console.WriteLine("");
        for (j = 0; j < FACE_NUM; j++)
        {
            cout = "  " + j.ToString().PadLeft(4);
            for (i = 0; i < order[j]; i++)
            {
                cout += "  " + node[i + j * ORDER_MAX].ToString().PadLeft(10);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  Nodal coordinates:");
        Console.WriteLine("");
        for (j = 0; j < NODE_NUM; j++)
        {
            cout = "  " + j.ToString().PadLeft(4);
            for (i = 0; i < 3; i++)
            {
                cout += "  " + coord[i + j * 3].ToString(CultureInfo.InvariantCulture).PadLeft(10);
            }

            Console.WriteLine(cout);
        }

        double volume1 = Geometry.polyhedron_volume_3d(coord, ORDER_MAX, FACE_NUM, node, NODE_NUM,
            order);

        double volume2 = Geometry.polyhedron_volume_3d_2(coord, ORDER_MAX, FACE_NUM, node, NODE_NUM,
            order);

        Console.WriteLine("");
        Console.WriteLine("  Volume (method 1) = " + volume1 + "");
        Console.WriteLine("  Volume (method 2) = " + volume2 + "");
        Console.WriteLine("  Volume (exact)    = " + volume_exact + "");

    }

}