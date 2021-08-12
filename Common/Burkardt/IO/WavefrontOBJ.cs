using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Types;

namespace Burkardt.IO
{
    public static class WavefrontOBJ
    {
        public static void obj_face_node_print(int face_num, int order_max, int[] face_order,
                int[] face_node)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    OBJ_FACE_NODE_PRINT prints the node indices for each face.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    25 January 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int FACE_NUM, the number of faces.
            //
            //    Input, int ORDER_MAX, the maximum number of vertices
            //    per face.
            //
            //    Input, int FACE_ORDER[FACE_NUM], the number of vertices
            //    per face.
            //
            //    Input, int FACE_NODE[ORDER_MAX*FACE_NUM], the nodes that
            //    make up each face.
            //
        {
            int face;
            int i;
            int order;

            Console.WriteLine("");
            Console.WriteLine("    Face   Order      Nodes");
            Console.WriteLine("");

            for (face = 0; face < face_num; face++)
            {
                string cout = "";
                order = face_order[face];
                cout += "  " + face.ToString().PadLeft(6)
                             + "  " + order.ToString().PadLeft(6);
                for (i = 0; i < order; i++)
                {
                    cout += "  " + face_node[i + face * order_max].ToString().PadLeft(6);
                }

                Console.WriteLine(cout);
            }
        }

        public static void obj_normal_vector_print(int normal_num, double[] normal_vector)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    OBJ_NORMAL_VECTOR_PRINT prints the normal vectors.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 January 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NORMAL_NUM, the number of normal vectors.
            //
            //    Input, double NORMAL_VECTOR[3*NORMAL_NUM], the normal vectors.
            //
        {
            int i;
            int normal;

            Console.WriteLine("");
            Console.WriteLine("  Normal Vectors:");
            Console.WriteLine("");

            for (normal = 0; normal < normal_num; normal++)
            {
                string cout = "";
                cout += "  " + normal.ToString().PadLeft(6);
                for (i = 0; i < 3; i++)
                {
                    cout += "  " + normal_vector[i + normal * 3].ToString().PadLeft(14);
                }

                Console.WriteLine(cout);
            }
        }

        public static void obj_node_xyz_print(int node_num, double[] node_xyz)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    OBJ_NODE_XYZ_PRINT prints the node coordinates.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 January 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Input, double NODE_XYZ[3*NODE_NUM], the coordinates
            //    of the nodes.
            //
        {
            int i;
            int node;

            Console.WriteLine("");
            Console.WriteLine("    Node         Coordinates");
            Console.WriteLine("");

            for (node = 0; node < node_num; node++)
            {
                string cout = "";
                cout += "  " + node.ToString().PadLeft(6);
                for (i = 0; i < 3; i++)
                {
                    cout += "  " + node_xyz[i + node * 3].ToString().PadLeft(14);
                }

                Console.WriteLine(cout);
            }
        }

        public static void obj_size(string filename, ref int node_num, ref int face_num, ref int normal_num,
                ref int order_max)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    OBJ_SIZE determines sizes of graphics objects in an Alias OBJ file.
            //
            //  Example:
            //
            //    #  magnolia.obj
            //
            //    v -3.269770 -39.572201 0.876128
            //    v -3.263720 -39.507999 2.160890
            //    ...
            //    v 0.000000 -9.988540 0.000000
            //
            //    vn 1.0 0.0 0.0
            //    ...
            //    vn 0.0 1.0 0.0
            //
            //    f 8 9 11 10
            //    f 12 13 15 14
            //    ...
            //    f 788 806 774
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    31 March 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string FILENAME, the input file name.
            //
            //    Output, int *NODE_NUM, the number of points.
            //
            //    Output, int *FACE_NUM, the number of faces.
            //
            //    Output, int *NORMAL_NUM, the number of normal vectors.
            //
            //    Output, int *ORDER_MAX, the maximum face order.
            //
        {
            string[] input;
            int n;
            int text_num;
            //
            //  Initialize.
            //
            node_num = 0;
            face_num = 0;
            normal_num = 0;
            order_max = 0;
            text_num = 0;

            try
            {
                input = File.ReadAllLines(filename);
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("OBJ_SIZE - Fatal error!");
                Console.WriteLine("  Could not open file.");
                return;
            }

            text_num = 0;

            foreach (string text in input)
            {

                text_num = text_num + 1;

                if (typeMethods.s_len_trim(text) == 0)
                {
                    continue;
                }

                if (text[0] == '#')
                {
                    continue;
                }

                if (text[0] == '$')
                {
                    continue;
                }

                //
                //  F V1 V2 ... VN
                //  Face.
                //
                if (text[0] == 'f' || text[0] == 'F')
                {
                    string[] tokens = text.Split(new char[] {' '});
                    n = tokens.Length;
                    order_max = Math.Max(order_max, n - 1);
                    face_num = face_num + 1;
                }
                //
                //  VN X Y Z
                //  Vertex normals.
                //
                else if ((text[0] == 'v' || text[0] == 'V') &&
                         (text[1] == 'n' || text[1] == 'N'))
                {
                    normal_num = normal_num + 1;
                }
                //
                //  V X Y Z W
                //  Geometric vertex.
                //
                else if (text[0] == 'v' || text[0] == 'V')
                {
                    node_num = node_num + 1;
                }

            }

            Console.WriteLine("");
            Console.WriteLine("  Read " + text_num + " lines from \"" + filename + "\".");
        }

        public static void obj_size_print(string filename, int node_num, int face_num,
                int normal_num, int order_max)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    OBJ_SIZE_PRINT prints sizes associated with an OBJ file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    25 January 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string FILENAME, the name of the file.
            //
            //    Input, int NODE_NUM, the number of vertices defined.
            //
            //    Input, int FACE_NUM, the number of faces defined.
            //
            //    Input, int NORMAL_NUM, the number of normal
            //    vectors defined.
            //
            //    Input, int ORDER_MAX, the maximum number of vertices
            //    per face.
            //
        {
            Console.WriteLine("");
            Console.WriteLine("  Object sizes for OBJ file \"" + filename + "\".");
            Console.WriteLine("");
            Console.WriteLine("  Nodes              = " + node_num + "");
            Console.WriteLine("  Faces              = " + face_num + "");
            Console.WriteLine("  Maximum face order = " + order_max + "");
            Console.WriteLine("  Normal vectors     = " + normal_num + "");

            return;
        }

        public static void obj_vertex_normal_print(int order_max, int face_num, int[] face_order,
                int[] vertex_normal)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    OBJ_VERTEX_NORMAL_PRINT prints the normal vectors indices per vertex.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    27 January 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int ORDER_MAX, the maximum number of vertices
            //    per face.
            //
            //    Input, int FACE_NUM, the number of faces.
            //
            //    Input, int FACE_ORDER[FACE_NUM], the number of vertices
            //    per face.
            //
            //    Input, int VERTEX_NORMAL[ORDER_MAX*FACE_NUM], the
            //    indices of normal vectors per vertex.
            //
        {
            int face;
            int i;
            int order;

            Console.WriteLine("");
            Console.WriteLine("  Normal Vector Indices:");
            Console.WriteLine("");
            Console.WriteLine("    Face   Order");
            Console.WriteLine("");

            for (face = 0; face < face_num; face++)
            {
                string cout = "";
                order = face_order[face];
                cout += "  " + face.ToString().PadLeft(6)
                             + "  " + order.ToString().PadLeft(6)
                             + "  ";
                for (i = 0; i < order; i++)
                {
                    cout += "  " + vertex_normal[i + face * order_max].ToString().PadLeft(6);
                }

                Console.WriteLine(cout);
            }
        }

        public static void obj_write(string output_filename, int node_num, int face_num,
                int normal_num, int order_max, double[] node_xyz, int[] face_order,
                int[] face_node, double[] normal_vector, int[] vertex_normal)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    OBJ_WRITE writes graphics information to an Alias OBJ file.
            //
            //  Discussion:
            //
            //    If no normal vectors are supplied (NORMAL_NUM <= 0) then
            //    a simple format is used for the "F" records.  Otherwise,
            //    the "v//vn" format is used.
            //
            //  Example:
            //
            //    #  no_normals.obj
            //
            //    g Group002
            //
            //    v -3.269770 -39.572201 0.876128
            //    v -3.263720 -39.507999 2.160890
            //    ...
            //    v 0.000000 -9.988540 0.000000
            //
            //    f 8 9 11 10
            //    f 12 13 15 14
            //    ...
            //    f 788 806 774
            //
            //    #  normals_supplied.obj
            //
            //    g Group001
            //
            //    v -3.269770 -39.572201 0.876128
            //    v -3.263720 -39.507999 2.160890
            //    ...
            //    v 0.000000 -9.988540 0.000000
            //
            //    vn 0.0 1.0 0.0
            //    vn 1.0 0.0 0.0
            //    ...
            //    vn 0.0 0.0 1.0
            //
            //    f 8//1 9//2 11//3 10//4
            //    f 12//5 13//6 15//7 14//8
            //    ...
            //    f 788//800 806//803 774//807
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    25 January 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string OUTPUT_FILENAME, the name of the output file.
            //
            //    Input, int NODE_NUM, the number of points.
            //
            //    Input, int FACE_NUM, the number of faces.
            //
            //    Input, int NORMAL_NUM, the number of normal vectors.
            //
            //    Input, int ORDER_MAX, the maximum number of vertices
            //    per face.
            //
            //    Input, double NODE_XYZ[3*NODE_NUM], the coordinates of points.
            //
            //    Input, int FACE_ORDER[FACE_NUM], the number of vertices
            //    per face.
            //
            //    Input, int FACE_NODE[ORDER_MAX*FACE_NUM], the nodes
            //    making faces.
            //
            //    Input, double NORMAL_VECTOR[3*NORMAL_NUM], normal vectors.
            //
            //    Input, int VERTEX_NORMAL[ORDER_MAX*FACE_NUM], the
            //    indices of normal vectors per vertex.
            //
        {
            int face;
            int i;
            int node;
            int normal;
            List<string> output = new List<string>();
            int text_num;
            int vertex;
            double w;

            text_num = 0;

            output.Add("# " + output_filename + "");
            output.Add("# created by obj_io::obj_write.C");
            output.Add("");
            output.Add("g Group001");

            text_num = text_num + 4;
            //
            //  V: vertex coordinates.
            //  For some reason, a fourth "coordinate" may be recommended.
            //  What is its meaning?
            //
            if (0 < node_num)
            {
                output.Add("");
                text_num = text_num + 1;
            }

            w = 1.0;
            for (node = 0; node < node_num; node++)
            {
                string tmp = "v";
                for (i = 0; i < 3; i++)
                {
                    tmp += "  " + node_xyz[i + 3 * node];
                }

                output.Add(tmp + "  " + w + "");
                text_num = text_num + 1;
            }

            //
            //  VN: normal vectors.
            //
            if (0 < normal_num)
            {
                output.Add("");
                text_num = text_num + 1;

                for (normal = 0; normal < normal_num; normal++)
                {
                    string tmp = "vn";
                    for (i = 0; i < 3; i++)
                    {
                        tmp += "  " + normal_vector[i + normal * 3];
                    }

                    output.Add(tmp);
                    text_num = text_num + 1;
                }
            }

            //
            //  F: Faces, specified as a list of triples, one triple for each vertex:
            //     vertex index/vertex texture index/vertex normal index
            //
            if (0 < face_num)
            {
                output.Add("");
                text_num = text_num + 1;
            }

            for (face = 0; face < face_num; face++)
            {
                string tmp = "f";
                for (vertex = 0; vertex < face_order[face]; vertex++)
                {
                    tmp += "  " + face_node[vertex + face * order_max];
                    if (0 < normal_num)
                    {
                        tmp += "//" + vertex_normal[vertex + face * order_max];
                    }
                }

                output.Add(tmp);
                text_num = text_num + 1;
            }

            try
            {
                File.WriteAllLines(output_filename, output);
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("OBJ_WRITE - Fatal error!");
                Console.WriteLine("  Could not open the output file \"" + output_filename + "\".");
            }

            //
            //  Report.
            //
            if (false)
            {
                Console.WriteLine("");
                Console.WriteLine("OBJ_WRITE:");
                Console.WriteLine("  Wrote " + text_num + " text lines to \""
                                  + output_filename + "\"");
            }
        }
    }
}