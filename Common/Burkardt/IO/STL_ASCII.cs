using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Types;

namespace Burkardt.IO
{
    public static class STL_ASCII
    {
        public class STLData
        {
            public int stla_offset_value;
        }
        
        public static bool stla_check(string input_file_name)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    STLA_CHECK checks an ASCII StereoLithography file.
            //
            //  Example:
            //
            //    solid MYSOLID
            //      facet normal 0.4 0.4 0.2
            //        outerloop
            //          vertex  1.0 2.1 3.2
            //          vertex  2.1 3.7 4.5
            //          vertex  3.1 4.5 6.7
            //        end loop
            //      end facet
            //      ...
            //      facet normal 0.2 0.2 0.4
            //        outerloop
            //          vertex  2.0 2.3 3.4
            //          vertex  3.1 3.2 6.5
            //          vertex  4.1 5.5 9.0
            //        end loop
            //      end facet
            //    end solid MYSOLID
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    22 September 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    3D Systems, Inc,
            //    Stereolithography Interface Specification,
            //    October 1989.
            //
            //  Parameters:
            //
            //    Input, string INPUT_FILE_NAME, the name of the input file.
            //
            //    Output, bool STLA_CHECK, is TRUE if the file is legal.
            //
        {
            bool check;
            bool done = false;
            bool error = false;
            int i;
            string[] input;
            int lchar = 0;
            int state;
            string text;
            int text_num;
            int vertex = 0;
            string word1;
            string word2;
            typeMethods.WordData data = new typeMethods.WordData();

            state = 0;
            text_num = 0;
            //
            //  Open the file.
            //
            try
            {
                input = File.ReadAllLines(input_file_name);
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("STLA_CHECK - Fatal error!");
                Console.WriteLine("  Could not open the file \"" + input_file_name + "\".");
                check = false;
                return check;
            }

            //
            //  Read the next line of text.
            //
            int index = 0;
            for (;;)
            {
                try
                {
                    text = input[index];
                    index++;
                }
                catch (Exception e)
                {
                    if (state != 0 &&
                        state != 1)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("STLA_CHECK - Fatal error!");
                        Console.WriteLine("  File line number = " + text_num + "");
                        Console.WriteLine("  End-of-file, but model not finished.");
                        check = false;
                        return check;
                    }

                    break;
                }

                text_num = text_num + 1;

                done = true;
                //
                //  Read the first word in the line.
                //
                word1 = typeMethods.word_next_read(ref data, text, ref done);

                if (done)
                {
                    Console.WriteLine("");
                    Console.WriteLine("STLA_CHECK - Fatal error!");
                    Console.WriteLine("  File line number = " + text_num + "");
                    Console.WriteLine("  No information on line.");
                    check = false;
                    return check;
                }

                //
                //  "Doctor" the text, changing a beginning occurrence of:
                //
                //      END FACET to ENDFACET
                //      END LOOP to ENDLOOP
                //      END SOLID to ENDSOLID
                //      FACET NORMAL to FACETNORMAL
                //      OUTER LOOP to OUTERLOOP
                //
                if (typeMethods.s_eqi(word1, "END"))
                {
                    word2 = typeMethods.word_next_read(ref data, text, ref done);

                    if (!typeMethods.s_eqi(word2, "FACET") &&
                        !typeMethods.s_eqi(word2, "LOOP") &&
                        !typeMethods.s_eqi(word2, "SOLID"))
                    {
                        Console.WriteLine("");
                        Console.WriteLine("STLA_CHECK - Fatal error!");
                        Console.WriteLine("  File line number = " + text_num + "");
                        Console.WriteLine("  The tag END was followed by an illegal word:");
                        Console.WriteLine("  \"" + word2 + "\"");
                        Console.WriteLine("  when expecting \"FACET\", \"LOOP\", or \"SOLID\".");
                        check = false;
                        return check;
                    }

                    word1 = word1 + word2;

                }
                else if (typeMethods.s_eqi(word1, "FACET"))
                {
                    word2 = typeMethods.word_next_read(ref data, text, ref done);

                    if (!typeMethods.s_eqi(word2, "NORMAL"))
                    {
                        Console.WriteLine("");
                        Console.WriteLine("STLA_CHECK - Fatal error!");
                        Console.WriteLine("  File line number = " + text_num + "");
                        Console.WriteLine("  The tag FACET was followed by an illegal word:");
                        Console.WriteLine("  \"" + word2 + "\"");
                        Console.WriteLine("  when expecting \"NORMAL\".");
                        check = false;
                        return check;
                    }

                    word1 = word1 + word2;
                }
                else if (typeMethods.s_eqi(word1, "OUTER"))
                {
                    word2 = typeMethods.word_next_read(ref data, text, ref done);

                    if (!typeMethods.s_eqi(word2, "LOOP"))
                    {
                        Console.WriteLine("");
                        Console.WriteLine("STLA_CHECK - Fatal error!");
                        Console.WriteLine("  File line number = " + text_num + "");
                        Console.WriteLine("  The tag OUTER was followed by an illegal word:");
                        Console.WriteLine("  \"" + word2 + "\"");
                        Console.WriteLine("  when expecting \"LOOP\".");
                        check = false;
                        return check;
                    }

                    word1 = word1 + word2;
                }

                //
                //  This first word tells us what to do.
                //
                //  SOLID - begin a new solid.
                //    Valid in state 0, moves to state 1.
                //  ENDSOLID - end current solid.
                //    Valid in state 1, moves to state 0.
                //
                //  FACETNORMAL - begin a new facet.
                //    Valid in state 0 or 1, moves to state 2.
                //  ENDFACET - end current facet.
                //    Valid in state 2, moves to state 1.
                //
                //  OUTERLOOP - begin a list of vertices.
                //    Valid in state 2, moves to state 3.
                //  ENDLOOP - end vertex list.
                //    Valid in state 3, moves to state 2.
                //
                //  VERTEX - give coordinates of next vertex.
                //    Valid in state 3 if current vertex count is 0, 1 or 2.
                //
                //  End of file -
                //    Valid in state 0 or 1.
                //
                if (typeMethods.s_eqi(word1, "SOLID"))
                {
                    if (state != 0)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("STLA_CHECK - Fatal error!");
                        Console.WriteLine("  File line number = " + text_num + "");
                        Console.WriteLine("  A new SOLID statement was encountered, but we");
                        Console.WriteLine("  have not finished processing the current solid.");
                        check = false;
                        return check;
                    }

                    state = 1;
                }
                else if (typeMethods.s_eqi(word1, "ENDSOLID"))
                {
                    if (state != 1)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("STLA_CHECK - Fatal error!");
                        Console.WriteLine("  File line number = " + text_num + "");
                        Console.WriteLine("  An END SOLID statement was encountered, but");
                        Console.WriteLine("  either we have not begun a solid at all, or we");
                        Console.WriteLine("  are not at an appropriate point to finish the");
                        Console.WriteLine("  current solid.");
                        check = false;
                        return check;
                    }

                    state = 0;
                }
                else if (typeMethods.s_eqi(word1, "FACETNORMAL"))
                {
                    if (state != 0 && state != 1)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("STLA_CHECK - Fatal error!");
                        Console.WriteLine("  File line number = " + text_num + "");
                        Console.WriteLine("  Model not in right state for FACET.");
                        check = false;
                        return check;
                    }

                    state = 2;

                    for (i = 1; i <= 3; i++)
                    {
                        word2 = typeMethods.word_next_read(ref data, text, ref done);

                        if (done)
                        {
                            Console.WriteLine("");
                            Console.WriteLine("STLA_CHECK - Fatal error!");
                            Console.WriteLine("  File line number = " + text_num + "");
                            Console.WriteLine("  End of information while reading a component");
                            Console.WriteLine("  of the normal vector.");
                            check = false;
                            return check;
                        }

                        typeMethods.s_to_r8(word2, ref lchar, ref error);

                        if (error)
                        {
                            Console.WriteLine("");
                            Console.WriteLine("STLA_CHECK - Fatal error!");
                            Console.WriteLine("  File line number = " + text_num + "");
                            Console.WriteLine("  Error while reading a component of the normal vector.");
                            check = false;
                            return check;
                        }
                    }
                }
                else if (typeMethods.s_eqi(word1, "ENDFACET"))
                {
                    if (state != 2)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("STLA_CHECK - Fatal error!");
                        Console.WriteLine("  File line number = " + text_num + "");
                        Console.WriteLine("  Model not in right state for ENDFACET.");
                        check = false;
                        return check;
                    }

                    state = 1;
                }
                else if (typeMethods.s_eqi(word1, "OUTERLOOP"))
                {
                    if (state != 2)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("STLA_CHECK - Fatal error!");
                        Console.WriteLine("  File line number = " + text_num + "");
                        Console.WriteLine("  Model not in right state for OUTERLOOP.");
                        check = false;
                        return check;
                    }

                    state = 3;
                    vertex = 0;
                }
                else if (typeMethods.s_eqi(word1, "ENDLOOP"))
                {
                    if (state != 3)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("STLA_CHECK - Fatal error!");
                        Console.WriteLine("  File line number = " + text_num + "");
                        Console.WriteLine("  Model not in right state for ENDLOOP.");
                        check = false;
                        return check;
                    }

                    state = 2;
                }
                else if (typeMethods.s_eqi(word1, "VERTEX"))
                {
                    if (state != 3)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("STLA_CHECK - Fatal error!");
                        Console.WriteLine("  File line number = " + text_num + "");
                        Console.WriteLine("  Model not in right state for VERTEX.");
                        check = false;
                        return check;
                    }

                    if (3 <= vertex)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("STLA_CHECK - Fatal error!");
                        Console.WriteLine("  File line number = " + text_num + "");
                        Console.WriteLine("  More than 3 vertices specified for a face.");
                        check = false;
                        return check;
                    }

                    for (i = 1; i <= 3; i++)
                    {
                        word2 = typeMethods.word_next_read(ref data, text, ref done);

                        if (done)
                        {
                            Console.WriteLine("");
                            Console.WriteLine("STLA_CHECK - Fatal error!");
                            Console.WriteLine("  File line number = " + text_num + "");
                            Console.WriteLine("  The value of a vertex coordinate is missing.");
                            check = false;
                            return check;
                        }

                        typeMethods.s_to_r8(word2, ref lchar, ref error);

                        if (error)
                        {
                            Console.WriteLine("");
                            Console.WriteLine("STLA_CHECK - Fatal error!");
                            Console.WriteLine("  File line number = " + text_num + "");
                            Console.WriteLine("  The value of a vertex coordinate makes no sense.");
                            check = false;
                            return check;
                        }
                    }

                    vertex = vertex + 1;
                }
                else
                {
                    Console.WriteLine("");
                    Console.WriteLine("STLA_CHECK - Fatal error!");
                    Console.WriteLine("  File line number = " + text_num + "");
                    Console.WriteLine("  Unrecognized line in file.");
                    check = false;
                    return check;
                }
            }

            check = true;

            return check;
        }

        public static void stla_face_node_print(ref STLData data, int face_num, int[] face_node)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    STLA_FACE_NODE_PRINT prints the node indices for each face.
            //
            //  Discussion:
            //
            //    If the global variable OFFSET is set to 1, then it is assumed that 
            //    all indices are 1-based.  In that case, this routine will print
            //    face numbers from 1 to FACE_NUM.
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
            //  Reference:
            //
            //    3D Systems, Inc,
            //    Stereolithography Interface Specification,
            //    October 1989.
            //
            //  Parameters:
            //
            //    Input, int FACE_NUM, the number of faces.
            //
            //    Input, int FACE_NODE[3*FACE_NUM], the nodes that make up each face.
            //
        {
            int face;
            int offset;
            int vertex;

            offset = data.stla_offset_value;

            Console.WriteLine("");
            Console.WriteLine("    Face         Nodes");
            Console.WriteLine("");

            for (face = 0; face < face_num; face++)
            {
                string cout = "  " + (face + offset).ToString().PadLeft(6);
                for (vertex = 0; vertex < 3; vertex++)
                {
                    cout += "  " + face_node[vertex + face * 3].ToString().PadLeft(6);
                }

                Console.WriteLine(cout);
            }
        }

        public static double[] stla_face_normal_compute(ref STLData data, int node_num, int face_num, double[] node_xyz,
                int[] face_node)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    STLA_FACE_NORMAL_COMPUTE computes normal vectors for an ASCII StereoLithography file.
            //
            //  Discussion:
            //
            //    This routine computes the normal vector to each triangular face
            //    in the STLA solid.  If the nodes of each triangular face are
            //    listed in counterclockwise order (as seen from outside the solid),
            //    then the normal vectors will be properly outward facing.
            //
            //    The normal vectors will have unit Euclidean norm.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    22 September 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    3D Systems, Inc,
            //    Stereolithography Interface Specification,
            //    October 1989.
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Input, int FACE_NUM, the number of faces.
            //
            //    Input, double NODE_XYZ[3*NODE_NUM], the node coordinates.
            //
            //    Input, int FACE_NODE[3*FACE_NUM], the nodes making faces.
            //
            //    Input, int FACE_MAX, the maximum number of faces.
            //
            //    Output, double STLA_FACE_NORMAL_COMPUTE[3*FACE_NUM], the normal 
            //    vector at each face.
            //
        {
            int face;
            double[] face_normal;
            int i;
            int n1;
            int n2;
            int n3;
            double norm;
            int offset;
            double[] v1 = new double[3];
            double[] v2 = new double[3];
            double[] v3;

            offset = data.stla_offset_value;

            face_normal = new double[3 * face_num];

            for (face = 0; face < face_num; face++)
            {
                n1 = face_node[0 + face * 3] - offset;
                n2 = face_node[1 + face * 3] - offset;
                n3 = face_node[2 + face * 3] - offset;

                for (i = 0; i < 3; i++)
                {
                    v1[i] = node_xyz[i + n2 * 3] - node_xyz[i + n1 * 3];
                }

                for (i = 0; i < 3; i++)
                {
                    v2[i] = node_xyz[i + n3 * 3] - node_xyz[i + n1 * 3];
                }

                v3 = typeMethods.r8vec_cross_3d(v1, v2);

                norm = typeMethods.r8vec_length(3, v3);

                if (norm != 0.0)
                {
                    for (i = 0; i < 3; i++)
                    {
                        face_normal[i + face * 3] = v3[i] / norm;
                    }
                }
                else
                {
                    for (i = 0; i < 3; i++)
                    {
                        face_normal[i + face * 3] = v3[i];
                    }
                }
            }

            return face_normal;
        }

        public static void stla_face_normal_print(int face_num, double[] face_normal)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    STLA_FACE_NORMAL_PRINT prints the normal vectors.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    18 September 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    3D Systems, Inc,
            //    Stereolithography Interface Specification,
            //    October 1989.
            //
            //  Parameters:
            //
            //    Input, int FACE_NUM, the number of faces.
            //
            //    Input, double FACE_NORMAL[3*FACE_NUM], the normal vector at each face.
            //
        {
            int face;
            int i;

            Console.WriteLine("");
            Console.WriteLine("    Face         Normal Vectors");
            Console.WriteLine("");

            for (face = 0; face < face_num; face++)
            {
                string cout = "  " + face.ToString().PadLeft(6);
                for (i = 0; i < 3; i++)
                {
                    cout += "  " + face_normal[i + face * 3].ToString().PadLeft(14);
                }

                Console.WriteLine(cout);
            }
        }

        public static void stla_node_xyz_print(int node_num, double[] node_xyz)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    STLA_NODE_XYZ_PRINT prints the node coordinates.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    22 September 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    3D Systems, Inc,
            //    Stereolithography Interface Specification,
            //    October 1989.
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Input, double NODE_XYZ[3*FACE_NUM], the normal vector at each face.
            //
        {
            int i;
            int node;

            Console.WriteLine("");
            Console.WriteLine("    Node         Coordinates");
            Console.WriteLine("");

            for (node = 0; node < node_num; node++)
            {
                string cout = "  " + node.ToString().PadLeft(6);
                for (i = 0; i < 3; i++)
                {
                    cout += "  " + node_xyz[i + node * 3].ToString().PadLeft(14);
                }

                Console.WriteLine(cout);
            }
        }
        
        public static void stla_offset_set(ref STLData data, int offset)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    STLA_OFFSET_SET sets the STLA offset.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    22 September 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    3D Systems, Inc,
            //    Stereolithography Interface Specification,
            //    October 1989.
            //
            //  Parameters:
            //
            //    Input, int OFFSET, the new value for the STLA offset.
            //    This should only be 0 or 1.
            //
        {
            if (offset != 0 && offset != 1)
            {
                Console.WriteLine("");
                Console.WriteLine("STLA_OFFSET_SET - Fatal error!");
                Console.WriteLine("  Input values of OFFSET must be 0 or 1.");
                Console.WriteLine("  Illegal input value was " + offset + "");
                return;
            }

            data.stla_offset_value = offset;

        }

        public static bool stla_read(ref STLData stldata, string input_file_name, int node_num, int face_num,
                ref double[] node_xy, ref int[] face_node, ref double[] face_normal)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    STLA_READ reads graphics information from an ASCII StereoLithography file.
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
            //  Reference:
            //
            //    3D Systems, Inc,
            //    Stereolithography Interface Specification,
            //    October 1989.
            //
            //  Parameters:
            //
            //    Input, string INPUT_FILE_NAME, the name of the input file.
            //
            //    Input, int NODE_NUM, the number of vertices defined.
            //
            //    Input, int FACE_NUM, the number of faces defined.
            //
            //    Output, double NODE_XY[3*NODE_NUM], the coordinates of points.
            //
            //    Output, int FACE_NODE[3*FACE_NUM], the nodes that make up each face.
            //
            //    Output, double FACE_NORMAL[3*FACE_NUM], the normal vector 
            //    at each face.
            //
            //    Output, bool STLA_READ, is TRUE if an error occurred.
            //
        {
            bool done;
            double dval;
            bool error;
            int face;
            int i;
            string[] input;
            int lchar = 0;
            int node;
            int offset;
            int state;
            double[] temp = new double[3];
            string text;
            int text_num;
            int vertex = 0;
            string word1;
            string word2;
            typeMethods.WordData data = new typeMethods.WordData();

            error = false;
            state = 0;
            offset = stldata.stla_offset_value;
            text_num = 0;

            face = 0;
            node = 0;
            //
            //  Open the file.
            //
            try
            {
                input = File.ReadAllLines(input_file_name);
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("STLA_READ - Fatal error!");
                Console.WriteLine("  Could not open the file \"" + input_file_name + "\".");
                error = true;
                return error;
            }

            //
            //  Read the next line of text.
            //
            int index = 0;
            for (;;)
            {
                try
                {
                    text = input[index];
                    index++;
                }
                catch (Exception e)
                {
                    if (state != 0 && state != 1)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("STLA_READ - Fatal error!");
                        Console.WriteLine("  File line number = " + text_num + "");
                        Console.WriteLine("  End-of-file, but model not finished.");
                        error = true;
                        return error;
                    }

                    break;
                }

                text_num = text_num + 1;

                done = true;
                //
                //  Read the first word in the line.
                //
                word1 = typeMethods.word_next_read(ref data, text, ref done);

                if (done)
                {
                    Console.WriteLine("");
                    Console.WriteLine("STLA_READ - Fatal error!");
                    Console.WriteLine("  File line number = " + text_num + "");
                    Console.WriteLine("  No information on line.");
                    error = true;
                    return error;
                }

                //
                //  "Doctor" the text, changing a beginning occurrence of:
                //
                //      END FACET to ENDFACET
                //      END LOOP to ENDLOOP
                //      END SOLID to ENDSOLID
                //      FACET NORMAL to FACETNORMAL
                //      OUTER LOOP to OUTERLOOP
                //
                if (typeMethods.s_eqi(word1, "END"))
                {
                    word2 = typeMethods.word_next_read(ref data, text, ref done);

                    if (!typeMethods.s_eqi(word2, "FACET") &&
                        !typeMethods.s_eqi(word2, "LOOP") &&
                        !typeMethods.s_eqi(word2, "SOLID"))
                    {
                        Console.WriteLine("");
                        Console.WriteLine("STLA_READ - Fatal error!");
                        Console.WriteLine("  File line number = " + text_num + "");
                        Console.WriteLine("  The tag END was followed by an illegal word:");
                        Console.WriteLine("  \"" + word2 + "\"");
                        Console.WriteLine("  when expecting \"FACET\", \"LOOP\", or \"SOLID\".");
                        error = true;
                        return error;
                    }

                    word1 = word1 + word2;
                }
                else if (typeMethods.s_eqi(word1, "FACET"))
                {
                    word2 = typeMethods.word_next_read(ref data, text, ref done);

                    if (!typeMethods.s_eqi(word2, "NORMAL"))
                    {
                        Console.WriteLine("");
                        Console.WriteLine("STLA_READ - Fatal error!");
                        Console.WriteLine("  File line number = " + text_num + "");
                        Console.WriteLine("  The tag FACET was followed by an illegal word:");
                        Console.WriteLine("  \"" + word2 + "\"");
                        Console.WriteLine("  when expecting \"NORMAL\".");
                        error = true;
                        return error;
                    }

                    word1 = word1 + word2;
                }
                else if (typeMethods.s_eqi(word1, "OUTER"))
                {
                    word2 = typeMethods.word_next_read(ref data, text, ref done);

                    if (!typeMethods.s_eqi(word2, "LOOP"))
                    {
                        Console.WriteLine("");
                        Console.WriteLine("STLA_READ - Fatal error!");
                        Console.WriteLine("  File line number = " + text_num + "");
                        Console.WriteLine("  The tag OUTER was followed by an illegal word:");
                        Console.WriteLine("  \"" + word2 + "\"");
                        Console.WriteLine("  when expecting \"LOOP\".");
                        error = true;
                        return error;
                    }

                    word1 = word1 + word2;
                }

                //
                //  This first word tells us what to do.
                //
                //  SOLID - begin a new solid.
                //    Valid in state 0, moves to state 1.
                //  ENDSOLID - end current solid.
                //    Valid in state 1, moves to state 0.
                //
                //  FACETNORMAL - begin a new facet.
                //    Valid in state 0 or 1, moves to state 2.
                //  ENDFACET - end current facet.
                //    Valid in state 2, moves to state 1.
                //
                //  OUTERLOOP - begin a list of vertices.
                //    Valid in state 2, moves to state 3.
                //  ENDLOOP - end vertex list.
                //    Valid in state 3, moves to state 2.
                //
                //  VERTEX - give coordinates of next vertex.
                //    Valid in state 3 if current vertex count is 0, 1 or 2.
                //
                //  End of file -
                //    Valid in state 0 or 1.
                //
                if (typeMethods.s_eqi(word1, "SOLID"))
                {
                    if (state != 0)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("STLA_READ - Warning!");
                        Console.WriteLine("  Model not in right state for SOLID.");
                        Console.WriteLine("  File line number = " + text_num + "");
                        error = true;
                        return error;
                    }

                    state = 1;
                }
                else if (typeMethods.s_eqi(word1, "ENDSOLID"))
                {
                    if (state != 1)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("STLA_READ - Warning!");
                        Console.WriteLine("  Model not in right state for ENDSOLID.");
                        Console.WriteLine("  File line number = " + text_num + "");
                        error = true;
                        return error;
                    }

                    state = 0;
                }
                else if (typeMethods.s_eqi(word1, "FACETNORMAL"))
                {
                    if (state != 0 && state != 1)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("STLA_READ - Warning!");
                        Console.WriteLine("  Model not in right state for FACET.");
                        Console.WriteLine("  File line number = " + text_num + "");
                        error = true;
                        return error;
                    }

                    state = 2;

                    if (face_num <= face)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("STLA_READ - Warning!");
                        Console.WriteLine("  More faces being read than expected.");
                        Console.WriteLine("  File line number = " + text_num + "");
                        error = true;
                        return error;
                    }

                    for (i = 0; i < 3; i++)
                    {
                        face_normal[i + face * 3] = 0.0;
                        word2 = typeMethods.word_next_read(ref data, text, ref done);
                        if (!done)
                        {
                            dval = typeMethods.s_to_r8(word2, ref lchar, ref error);
                            if (error)
                            {
                                return error;
                            }

                            face_normal[i + face * 3] = dval;
                        }
                    }
                }
                else if (typeMethods.s_eqi(word1, "ENDFACET"))
                {
                    if (state != 2)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("STLA_READ - Warning!");
                        Console.WriteLine("  Model not in right state for ENDFACET.");
                        Console.WriteLine("  File line number = " + text_num + "");
                        error = true;
                        return error;
                    }

                    face = face + 1;
                    state = 1;
                }
                else if (typeMethods.s_eqi(word1, "OUTERLOOP"))
                {
                    if (state != 2)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("STLA_READ - Warning!");
                        Console.WriteLine("  Model not in right state for OUTERLOOP.");
                        Console.WriteLine("  File line number = " + text_num + "");
                        error = true;
                        return error;
                    }

                    state = 3;
                    vertex = 0;
                }
                else if (typeMethods.s_eqi(word1, "ENDLOOP"))
                {
                    if (state != 3)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("STLA_READ - Warning!");
                        Console.WriteLine("  Model not in right state for ENDLOOP.");
                        Console.WriteLine("  File line number = " + text_num + "");
                        error = true;
                        return error;
                    }

                    state = 2;
                }
                else if (typeMethods.s_eqi(word1, "VERTEX"))
                {
                    if (state != 3)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("STLA_READ - Warning!");
                        Console.WriteLine("  Model not in right state for VERTEX.");
                        Console.WriteLine("  File line number = " + text_num + "");
                        error = true;
                        return error;
                    }

                    if (3 <= vertex)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("STLA_READ - Warning!");
                        Console.WriteLine("  Too many vertices for face.");
                        Console.WriteLine("  File line number = " + text_num + "");
                        error = true;
                        return error;
                    }

                    for (i = 0; i < 3; i++)
                    {
                        word2 = typeMethods.word_next_read(ref data, text, ref done);
                        if (done)
                        {
                            error = true;
                            return error;
                        }

                        dval = typeMethods.s_to_r8(word2, ref lchar, ref error);
                        if (error)
                        {
                            return error;
                        }

                        temp[i] = dval;
                    }

                    if (node_num <= node)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("STLA_READ - Warning!");
                        Console.WriteLine("  More nodes being read than expected.");
                        Console.WriteLine("  File line number = " + text_num + "");
                        error = true;
                        return error;
                    }

                    for (i = 0; i < 3; i++)
                    {
                        node_xy[i + node * 3] = temp[i];
                    }

                    face_node[vertex + face * 3] = node + offset;

                    node = node + 1;
                    vertex = vertex + 1;
                }
                else
                {
                    Console.WriteLine("");
                    Console.WriteLine("STLA_READ - Warning!");
                    Console.WriteLine("  Unrecognized line in file.");
                    Console.WriteLine("  File line number = " + text_num + "");
                    error = true;
                    return error;
                }
            }

            return error;
        }


        public static void stla_size(string input_file_name, ref int solid_num, ref int node_num,
                ref int face_num, ref int text_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    STLA_SIZE determines sizes associated with an STLA file.
            //
            //  Discussion:
            //
            //    This routine assumes that the file is a legal STLA file.
            //
            //    To perform checks on the file, call STLA_CHECK first.
            //
            //    Note that the counts for the number of nodes and edges are
            //    overestimates, since presumably, most nodes will be defined several
            //    times, once for each face they are part of, and most edges will
            //    be defined twice.
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
            //  Reference:
            //
            //    3D Systems, Inc,
            //    Stereolithography Interface Specification,
            //    October 1989.
            //
            //  Parameters:
            //
            //    Input, string INPUT_FILE_NAME, the name of the input file.
            //
            //    Output, int *SOLID_NUM, the number of solids defined.
            //    Presumably, this is 1.
            //
            //    Output, int *NODE_NUM, the number of vertices defined.
            //
            //    Output, int *FACE_NUM, the number of faces defined.
            //
            //    Output, int *TEXT_NUM, the number of lines of text.
            //
        {
            bool done = false;
            bool error = false;
            int i;
            string[] input;
            int lchar = 0;
            int state = 0;
            string text = "";
            int vertex = 0;
            string word1;
            string word2;
            typeMethods.WordData data = new typeMethods.WordData();

            state = 0;

            text_num = 0;
            solid_num = 0;
            node_num = 0;
            face_num = 0;
            //
            //  Open the file.
            //
            try
            {
                input = File.ReadAllLines(input_file_name);
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("STLA_SIZE - Fatal error!");
                Console.WriteLine("  Could not open the file \"" + input_file_name + "\".");
                return;
            }

            //
            //  Read the next line of text.
            //
            int index = 0;
            for (;;)
            {
                try
                {
                    text = input[index];
                    index++;
                }
                catch (Exception e)
                {
                    if (state != 0 &&
                        state != 1)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("STLA_SIZE - Fatal error!");
                        Console.WriteLine("  File line number = " + text_num + "");
                        Console.WriteLine("  End-of-file, but model not finished.");
                        return;
                    }

                    break;
                }

                text_num = text_num + 1;

                done = true;
                //
                //  Read the first word in the line.
                //
                word1 = typeMethods.word_next_read(ref data, text, ref done);

                if (done)
                {
                    Console.WriteLine("");
                    Console.WriteLine("STLA_CHECK - Fatal error!");
                    Console.WriteLine("  File line number = " + text_num + "");
                    Console.WriteLine("  No information on line.");
                    return;
                }

                //
                //  "Doctor" the text, changing a beginning occurrence of:
                //
                //      END FACET to ENDFACET
                //      END LOOP to ENDLOOP
                //      END SOLID to ENDSOLID
                //      FACET NORMAL to FACETNORMAL
                //      OUTER LOOP to OUTERLOOP
                //
                if (typeMethods.s_eqi(word1, "END"))
                {
                    word2 = typeMethods.word_next_read(ref data, text, ref done);

                    if (!typeMethods.s_eqi(word2, "FACET") &&
                        !typeMethods.s_eqi(word2, "LOOP") &&
                        !typeMethods.s_eqi(word2, "SOLID"))
                    {
                        Console.WriteLine("");
                        Console.WriteLine("STLA_CHECK - Fatal error!");
                        Console.WriteLine("  File line number = " + text_num + "");
                        Console.WriteLine("  The tag END was followed by an illegal word:");
                        Console.WriteLine("  \"" + word2 + "\"");
                        Console.WriteLine("  when expecting \"FACET\", \"LOOP\", or \"SOLID\".");
                        return;
                    }

                    word1 = word1 + word2;
                }
                else if (typeMethods.s_eqi(word1, "FACET"))
                {
                    word2 = typeMethods.word_next_read(ref data, text, ref done);

                    if (!typeMethods.s_eqi(word2, "NORMAL"))
                    {
                        Console.WriteLine("");
                        Console.WriteLine("STLA_CHECK - Fatal error!");
                        Console.WriteLine("  File line number = " + text_num + "");
                        Console.WriteLine("  The tag FACET was followed by an illegal word:");
                        Console.WriteLine("  \"" + word2 + "\"");
                        Console.WriteLine("  when expecting \"NORMAL\".");
                        return;
                    }

                    word1 = word1 + word2;
                }
                else if (typeMethods.s_eqi(word1, "OUTER"))
                {
                    word2 = typeMethods.word_next_read(ref data, text, ref done);

                    if (!typeMethods.s_eqi(word2, "LOOP"))
                    {
                        Console.WriteLine("");
                        Console.WriteLine("STLA_CHECK - Fatal error!");
                        Console.WriteLine("  File line number = " + text_num + "");
                        Console.WriteLine("  The tag OUTER was followed by an illegal word:");
                        Console.WriteLine("  \"" + word2 + "\"");
                        Console.WriteLine("  when expecting \"LOOP\".");
                        return;
                    }

                    word1 = word1 + word2;
                }

                //
                //  This first word tells us what to do.
                //
                //  SOLID - begin a new solid.
                //    Valid in state 0, moves to state 1.
                //  ENDSOLID - end current solid.
                //    Valid in state 1, moves to state 0.
                //
                //  FACETNORMAL - begin a new facet.
                //    Valid in state 0 or 1, moves to state 2.
                //  ENDFACET - end current facet.
                //    Valid in state 2, moves to state 1.
                //
                //  OUTERLOOP - begin a list of vertices.
                //    Valid in state 2, moves to state 3.
                //  ENDLOOP - end vertex list.
                //    Valid in state 3, moves to state 2.
                //
                //  VERTEX - give coordinates of next vertex.
                //    Valid in state 3 if current vertex count is 0, 1 or 2.
                //
                //  End of file -
                //    Valid in state 0 or 1.
                //
                if (typeMethods.s_eqi(word1, "SOLID"))
                {
                    if (state != 0)
                    {
                        return;
                    }

                    state = 1;
                }
                else if (typeMethods.s_eqi(word1, "ENDSOLID"))
                {
                    if (state != 1)
                    {
                        return;
                    }

                    state = 0;

                    solid_num = solid_num + 1;
                }
                else if (typeMethods.s_eqi(word1, "FACETNORMAL"))
                {
                    if (state != 0 && state != 1)
                    {
                        return;
                    }

                    state = 2;

                    for (i = 1; i <= 3; i++)
                    {
                        word2 = typeMethods.word_next_read(ref data, text, ref done);

                        if (done)
                        {
                            return;
                        }

                        typeMethods.s_to_r8(word2, ref lchar, ref error);

                        if (error)
                        {
                            return;
                        }
                    }
                }
                else if (typeMethods.s_eqi(word1, "ENDFACET"))
                {
                    if (state != 2)
                    {
                        return;
                    }

                    state = 1;
                    face_num = face_num + 1;
                }
                else if (typeMethods.s_eqi(word1, "OUTERLOOP"))
                {
                    if (state != 2)
                    {
                        return;
                    }

                    state = 3;
                    vertex = 0;
                }
                else if (typeMethods.s_eqi(word1, "ENDLOOP"))
                {
                    if (state != 3)
                    {
                        return;
                    }

                    state = 2;
                }
                else if (typeMethods.s_eqi(word1, "VERTEX"))
                {
                    if (state != 3)
                    {
                        return;
                    }

                    if (3 <= vertex)
                    {
                        return;
                    }

                    for (i = 1; i <= 3; i++)
                    {
                        word2 = typeMethods.word_next_read(ref data, text, ref done);

                        if (done)
                        {
                            return;
                        }

                        typeMethods.s_to_r8(word2, ref lchar, ref error);

                        if (error)
                        {
                            return;
                        }

                    }

                    vertex = vertex + 1;
                    node_num = node_num + 1;
                }
                else
                {
                    return;
                }
            }
        }


        public static void stla_size_print(ref STLData data, string input_file_name, int solid_num, int node_num,
                int face_num, int text_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    STLA_SIZE_PRINT prints sizes associated with an STLA file.
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
            //  Reference:
            //
            //    3D Systems, Inc,
            //    Stereolithography Interface Specification,
            //    October 1989.
            //
            //  Parameters:
            //
            //    Input, string INPUT_FILE_NAME, the name of the input file, or the
            //    name of the object.
            //
            //    Input, int SOLID_NUM, the number of solids defined.
            //
            //    Input, int NODE_NUM, the number of vertices defined.
            //
            //    Input, int FACE_NUM, the number of faces defined.
            //
            //    Input, int TEXT_NUM, the number of lines of text in the file.
            //
        {
            Console.WriteLine("");
            Console.WriteLine("  Sizes for STLA object \"" + input_file_name + "\".");
            Console.WriteLine("");
            Console.WriteLine("  Solids =                   " + solid_num + "");
            Console.WriteLine("  Nodes (may be repeated) =  " + node_num + "");
            Console.WriteLine("  Faces (triangular only) =  " + face_num + "");
            Console.WriteLine("");
            Console.WriteLine("  The index offset value =   " + data.stla_offset_value + "");
            Console.WriteLine("  Number of lines of text =  " + text_num + "");

        }

        public static void stla_write(ref STLData data, string output_file_name, int node_num, int face_num,
                double[] node_xyz, int[] face_node, double[] face_normal)

            //****************************************************************************
            //
            //  Purpose:
            //   
            //    STLA_WRITE writes an ASCII STL (stereolithography) file.
            //
            //  Example:
            //
            //    solid MYSOLID
            //      facet normal 0.4 0.4 0.2
            //        outerloop
            //          vertex  1.0 2.1 3.2
            //          vertex  2.1 3.7 4.5
            //          vertex  3.1 4.5 6.7
            //        end loop
            //      end facet
            //      ...
            //      facet normal 0.2 0.2 0.4
            //        outerloop
            //          vertex  2.0 2.3 3.4
            //          vertex  3.1 3.2 6.5
            //          vertex  4.1 5.5 9.0
            //        end loop
            //      end facet
            //    end solid MYSOLID
            //
            //  Discussion:
            //
            //    The polygons in an STL file should only be triangular.  This routine 
            //    will try to automatically decompose higher-order polygonal faces into 
            //    suitable triangles, without actually modifying the internal graphics 
            //    data.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    22 September 1998
            //
            //  Author:
            // 
            //    John Burkardt
            //
            //  Reference:
            //
            //    3D Systems, Inc,
            //    Stereolithography Interface Specification,
            //    October 1989.
            //
            //  Parameters:
            //
            //    Input, string OUTPUT_FILE_NAME, the name of the output file.
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Input, int FACE_NUM, the number of faces.
            //
            //    Input, double NODE_XYZ[3*NODE_NUM], the node coordinates.
            //
            //    Input, int FACE_NODE[3*FACE_NUM], the nodes making faces.
            //
            //    Input, int FACE_MAX, the maximum number of faces.
            //
            //    Input, double FACE_NORMAL[3*FACE_NUM], the normal vector at each face.
            //
        {
            int face;
            int i;
            int node;
            int offset;
            List<string> output_unit = new List<string>();
            int vertex;

            offset = data.stla_offset_value;

            //
            //  Initialize.
            //
            output_unit.Add("solid MYSOLID");

            for (face = 0; face < face_num; face++)
            {
                string cout = "  facet normal";
                for (i = 0; i < 3; i++)
                {
                    cout += "  " + face_normal[i + face * 3].ToString().PadLeft(10);
                }

                output_unit.Add(cout);
                output_unit.Add("    outer loop");
                for (vertex = 0; vertex < 3; vertex++)
                {
                    node = face_node[vertex + face * 3] - offset;
                    cout = "      vertex  ";
                    for (i = 0; i < 3; i++)
                    {
                        cout += "  " + node_xyz[i + node * 3].ToString().PadLeft(10);
                    }

                    output_unit.Add(cout);
                }

                output_unit.Add("    end loop");
                output_unit.Add("  end facet");
            }

            output_unit.Add("end solid MYSOLID");


            try
            {
                File.WriteAllLines(output_file_name, output_unit);
            }
            catch (Exception e)
            {
                Console.WriteLine("");
                Console.WriteLine("STLA_WRITE - Fatal error!");
                Console.WriteLine("  Could not open the file \"" + output_file_name + "\".");
            }
        }
    }
}