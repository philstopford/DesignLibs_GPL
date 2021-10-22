using System;
using Burkardt.Types;

namespace Burkardt.Geometry
{
    public static class Dual
    {
        public static void dual_shape_3d(int point_num, int face_num, int face_order_max,
                double[] point_coord, int[] face_order, int[] face_point, int point_num2,
                int face_num2, int face_order_max2, ref double[] point_coord2, ref int[] face_order2,
                ref int[] face_point2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DUAL_SHAPE_3D constructs the dual of a shape in 3D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 August 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int POINT_NUM, the number of points.
            //
            //    Input, int FACE_NUM, the number of faces.
            //
            //    Input, int FACE_ORDER_MAX, the maximum number of vertices per face.
            //
            //    Input, double POINT_COORD[3*POINT_NUM], the point coordinates.
            //
            //    Input, int FACE_ORDER[FACE_NUM], the number of vertices per face.
            //
            //    Input, int FACE_POINT[FACE_ORDER_MAX*FACE_NUM]; FACE_POINT(I,J)
            //    contains the index of the I-th point in the J-th face.  The
            //    points are listed in the counter clockwise direction defined
            //    by the outward normal at the face.
            //
            //    Input, int POINT_NUM2, the number of points in the dual.
            //
            //    Input, int FACE_NUM2, the number of faces in the dual.
            //
            //    Input, int FACE_ORDER_MAX2, the maximum number of vertices per face
            //    in the dual.
            //
            //    Output, double POINT_COORD2[3*POINT_NUM2], the point coordinates
            //    of the dual.
            //
            //    Input, int FACE_ORDER2[FACE_NUM2], the number of vertices
            //    per face.
            //
            //    Output, int FACE_POINT2[FACE_ORDER_MAX2*FACE_NUM2], the vertices
            //    of each face in the dual.
            //
        {
            int col = 0;
            int face;
            int i;
            int inext;
            int iprev;
            int istop;
            int j;
            int k;
            double norm;
            int row = 0;
            double x;
            double y;
            double z;
            //
            //  This computation should really compute the center of gravity
            //  of the face, in the general case.
            //
            //  We'll also assume the vertices of the original and the dual
            //  are to lie on the unit sphere, so we can normalize the
            //  position vector of the vertex.
            //
            for (face = 0; face < face_num; face++)
            {
                x = 0.0;
                y = 0.0;
                z = 0.0;
                for (j = 0; j < face_order[face]; j++)
                {
                    k = face_point[j + face * face_order_max];
                    x = x + point_coord[0 + (k - 1) * 3];
                    y = y + point_coord[1 + (k - 1) * 3];
                    z = z + point_coord[2 + (k - 1) * 3];
                }

                norm = Math.Sqrt(x * x + y * y + z * z);

                point_coord2[0 + face * face_order_max2] = x / norm;
                point_coord2[1 + face * face_order_max2] = y / norm;
                point_coord2[2 + face * face_order_max2] = z / norm;
            }

            //
            //  Now build the face in the dual associated with each node FACE.
            //
            for (face = 1; face <= face_num2; face++)
            {
                //
                //  Initialize the order.
                //
                face_order2[face - 1] = 0;
                //
                //  Find the first occurrence of FACE in an edge of polyhedron.
                //  ROW and COL are 1-based indices.
                //
                typeMethods.i4col_find_item(face_order_max, face_num, face_point, face,
                    ref row, ref col);

                if (row <= 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("DUAL_SHAPE_3D - Fatal error!");
                    Console.WriteLine("  Could not find an edge using node " + face + "");
                    return;
                }

                //
                //  Save the following node as ISTOP.
                //  When we encounter ISTOP again, this will mark the end of our search.
                //
                i = row + 1;
                if (face_order[col - 1] < i)
                {
                    i = 1;
                }

                istop = face_point[i - 1 + (col - 1) * face_order_max];
                //
                //  Save the previous node as INEXT.
                //
                for (;;)
                {
                    i = row - 1;
                    if (i < 1)
                    {
                        i = i + face_order[col - 1];
                    }

                    inext = face_point[i - 1 + (col - 1) * face_order_max];

                    face_order2[face - 1] = face_order2[face - 1] + 1;
                    face_point2[face_order2[face - 1] - 1 + (face - 1) * face_order_max2] = col;
                    //
                    //  If INEXT != ISTOP, continue.
                    //
                    if (inext == istop)
                    {
                        break;
                    }

                    //
                    //  Set IPREV:= INEXT.
                    //
                    iprev = inext;
                    //
                    //  Search for the occurrence of the edge FACE-IPREV.
                    //  ROW and COL are 1-based indices.
                    //
                    typeMethods.i4col_find_pair_wrap(face_order_max, face_num, face_point, face,
                        iprev, ref row, ref col);

                    if (row <= 0)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("DUAL_SHAPE_3D - Fatal error!");
                        Console.WriteLine("  No edge from node " + iprev + "");
                        Console.WriteLine("  to node " + face + "");
                        return;
                    }
                }
            }

            return;
        }

        public static void dual_size_3d(int point_num, int edge_num, int face_num,
                int face_order_max, double[] point_coord, int[] face_order, int[] face_point,
                ref int point_num2, ref int edge_num2, ref int face_num2, ref int face_order_max2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DUAL_SIZE_3D determines sizes for a dual of a shape in 3D.
            //
            //  Discussion:
            //
            //    We don't actually need FACE_POINT as input here.  But since the
            //    three arrays occur together everywhere else, it seems unnecessarily
            //    user-confusing to vary the usage here!
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
            //  Parameters:
            //
            //    Input, int POINT_NUM, the number of points.
            //
            //    Input, int EDGE_NUM, the number of edges.
            //
            //    Input, int FACE_NUM, the number of faces.
            //
            //    Input, int FACE_ORDER_MAX, the maximum number of vertices per face.
            //
            //    Input, double POINT_COORD[3*POINT_NUM], the point coordinates.
            //
            //    Input, int FACE_ORDER[FACE_NUM], the number of vertices per face.
            //
            //    Input, int FACE_POINT[FACE_ORDER_MAX*FACE_NUM]; FACE_POINT(I,J)
            //    contains the index of the I-th point in the J-th face.  The
            //    points are listed in the counter clockwise direction defined
            //    by the outward normal at the face.
            //
            //    Output, int *POINT_NUM2, the number of points in the dual.
            //
            //    Output, int *EDGE_NUM2, the number of edges in the dual.
            //
            //    Output, int *FACE_NUM2, the number of faces in the dual.
            //
            //    Output, int *FACE_ORDER_MAX2, the maximum number of vertices per face
            //    in the dual.
            //
        {
            int i;
            int face;
            int[] face_order2;
            int face2;
            //
            //  These values are easy to compute:
            //
            point_num2 = face_num;
            edge_num2 = edge_num;
            face_num2 = point_num;
            //
            //  To determine FACE_ORDER_MAX2 is not so easy.
            //  You have to construct the FACE_ORDER array for the dual shape.
            //  The order of a dual face is the number of edges that the vertex occurs in.
            //  But then all we have to do is count how many times each item shows up
            //  in the FACE_POINT array.
            //
            face_order2 = new int[(face_num2)];

            for (i = 0; i < face_num2; i++)
            {
                face_order2[i] = 0;
            }

            for (face = 0; face < face_num; face++)
            {
                for (i = 0; i < face_order[face]; i++)
                {
                    face2 = face_point[i + face * face_order_max];
                    if (face2 == 0)
                    {
                        Console.WriteLine("WHOA!");
                    }

                    ;
                    face_order2[face2 - 1] = face_order2[face2 - 1] + 1;
                }
            }

            face_order_max2 = 0;
            for (i = 0; i < face_num2; i++)
            {
                face_order_max2 = Math.Max(face_order_max2, face_order2[i]);
            }
        }
    }
}