using System;
using Burkardt.IO;
using Burkardt.Types;

namespace STLASCIITest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for STLA_IO_TEST.
        //
        //  Discussion:
        //
        //    STLA_IO_TEST tests the STLA_IO library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("STLA_IO_TEST:");
        Console.WriteLine("  Test the STLA_IO library.");

        test01("cube.stla");
        test02("cube.stla");
        test03("cube.stla");
        test04("cube_new.stla");
        test05();
        Console.WriteLine("");
        Console.WriteLine("STLA_IO_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");

    }

    private static void test01(string input_file_name)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests STLA_CHECK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  STLA_CHECK makes some simple checks on a file.");

        Console.WriteLine("");
        if (STL_ASCII.stla_check(input_file_name))
        {
            Console.WriteLine("  The file \"" + input_file_name
                                              + "\" seems to be a legal ASCII STL file.");
        }
        else
        {
            Console.WriteLine("  The file \"" + input_file_name
                                              + "\" does NOT seem to be a legal ASCII STL file.");
        }

    }

    private static void test02(string input_file_name)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests STLA_SIZE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int face_num = 0;
        int node_num = 0;
        int solid_num = 0;
        int text_num = 0;
        STL_ASCII.STLData data = new();

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  STLA_SIZE determines the size of various objects");
        Console.WriteLine("  in an ASCII STL file.");

        STL_ASCII.stla_size(input_file_name, ref solid_num, ref node_num, ref face_num, ref text_num);

        STL_ASCII.stla_size_print(ref data, input_file_name, solid_num, node_num, face_num, text_num);

    }

    private static void test03(string input_file_name)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests STLA_READ.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int face_num = 0;
        int node_num = 0;
        int solid_num = 0;
        int text_num = 0;
        STL_ASCII.STLData data = new();

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  STLA_READ reads an object in an ASCII STL file.");

        STL_ASCII.stla_size(input_file_name, ref solid_num, ref node_num, ref face_num, ref text_num);

        int[] face_node = new int[3 * face_num];
        double[] face_normal = new double[3 * face_num];
        double[] node_xyz = new double[3 * node_num];

        bool error = STL_ASCII.stla_read(ref data, input_file_name, node_num, face_num, ref node_xyz,
            ref face_node, ref face_normal);

        switch (error)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("  STLA_READ returned error flag.");
                return;
        }

        STL_ASCII.stla_size_print(ref data, input_file_name, solid_num, node_num, face_num, text_num);

        STL_ASCII.stla_face_node_print(ref data, face_num, face_node);
        STL_ASCII.stla_face_normal_print(face_num, face_normal);
        STL_ASCII.stla_node_xyz_print(node_num, node_xyz);

    }

    private static void test04(string output_file_name)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests STLA_WRITE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int FACE_NUM = 12;
        const int NODE_NUM = 8;

        int[] face_node =
        {
            1, 3, 2,
            2, 3, 4,
            1, 6, 5,
            1, 2, 6,
            3, 7, 4,
            4, 7, 8,
            5, 6, 8,
            5, 8, 7,
            1, 5, 7,
            1, 7, 3,
            2, 4, 6,
            6, 4, 8
        };
        double[] face_normal =
        {
            0.0, 0.0, -1.0,
            0.0, 0.0, -1.0,
            0.0, -1.0, 0.0,
            0.0, -1.0, 0.0,
            0.0, +1.0, 0.0,
            0.0, +1.0, 0.0,
            0.0, 0.0, +1.0,
            0.0, 0.0, +1.0,
            -1.0, 0.0, 0.0,
            -1.0, 0.0, 0.0,
            +1.0, 0.0, 0.0,
            +1.0, 0.0, 0.0
        };
        double[] node_xyz =
        {
            0.0, 0.0, 0.0,
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            1.0, 1.0, 0.0,
            0.0, 0.0, 1.0,
            1.0, 0.0, 1.0,
            0.0, 1.0, 1.0,
            1.0, 1.0, 1.0
        };
        STL_ASCII.STLData data = new();

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  STLA_WRITE writes an ASCII STL file.");

        const int offset = 1;
        STL_ASCII.stla_offset_set(ref data, offset);

        STL_ASCII.stla_write(ref data, output_file_name, NODE_NUM, FACE_NUM, node_xyz,
            face_node, face_normal);

        Console.WriteLine("");
        Console.WriteLine("  Graphics data was written to the STLA file \""
                          + output_file_name + "\".");

    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests STLA_FACE_NORMAL_COMPUTE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int FACE_NUM = 12;
        const int NODE_NUM = 8;

        int face;
        int[] face_node =
        {
            1, 3, 2,
            2, 3, 4,
            1, 6, 5,
            1, 2, 6,
            3, 7, 4,
            4, 7, 8,
            5, 6, 8,
            5, 8, 7,
            1, 5, 7,
            1, 7, 3,
            2, 4, 6,
            6, 4, 8
        };
        double[] face_normal =
        {
            0.0, 0.0, -1.0,
            0.0, 0.0, -1.0,
            0.0, -1.0, 0.0,
            0.0, -1.0, 0.0,
            0.0, +1.0, 0.0,
            0.0, +1.0, 0.0,
            0.0, 0.0, +1.0,
            0.0, 0.0, +1.0,
            -1.0, 0.0, 0.0,
            -1.0, 0.0, 0.0,
            +1.0, 0.0, 0.0,
            +1.0, 0.0, 0.0
        };
        double[] node_xyz =
        {
            0.0, 0.0, 0.0,
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            1.0, 1.0, 0.0,
            0.0, 0.0, 1.0,
            1.0, 0.0, 1.0,
            0.0, 1.0, 1.0,
            1.0, 1.0, 1.0
        };
        STL_ASCII.STLData data = new();

        const int offset = 1;
        STL_ASCII.stla_offset_set(ref data, offset);

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  STLA_FACE_NORMAL_COMPUTE computes the face normal");
        Console.WriteLine("  vectors for an STLA file.");
        Console.WriteLine("");
        Console.WriteLine("  We have an STLA solid, and its exact normals.");
        Console.WriteLine("  We now call STLA_FACE_NORMAL_COMPUTE to");
        Console.WriteLine("  recompute the normals.");

        double[] face_normal2 = STL_ASCII.stla_face_normal_compute(ref data, NODE_NUM, FACE_NUM, node_xyz,
            face_node);

        Console.WriteLine("");
        Console.WriteLine("  We print out the maximum error, defined as");
        Console.WriteLine("    |1 - dot ( n1, n2 )|");
        Console.WriteLine("  where n1 and n2 are the exact and computed normals.");

        double dot_max = 0.0;

        for (face = 0; face < FACE_NUM; face++)
        {
            dot_max = Math.Max(dot_max, Math.Abs(1.0 -
                                                 typeMethods.r8vec_dot(3, face_normal, face_normal2, +face * 3,
                                                     +face * 3)));
        }

        Console.WriteLine("");
        Console.WriteLine("  Maximum error = " + dot_max + "");

    }
}