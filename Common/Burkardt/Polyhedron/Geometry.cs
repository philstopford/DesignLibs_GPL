using System;
using Burkardt.Types;

namespace Burkardt.Polyhedron
{
    public static class Geometry
    {
        public static double polyhedron_area_3d(double[] coord, int order_max, int face_num,
                int[] node, int node_num, int[] order)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYHEDRON_AREA_3D computes the surface area of a polyhedron in 3D.
            //
            //  Restriction:
            //
            //    The computation is not valid unless the faces of the polyhedron
            //    are planar polygons.
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
            //  Reference:
            //
            //    Allen Van Gelder,
            //    Efficient Computation of Polygon Area and Polyhedron Volume,
            //    Graphics Gems V, edited by Alan Paeth,
            //    AP Professional, 1995, T385.G6975.
            //
            //  Parameters:
            //
            //    Input, double COORD[NODE_NUM*3], the 3D coordinates of the vertices.
            //
            //    Input, int ORDER_MAX, the maximum number of vertices that make
            //    up a face of the polyhedron.
            //
            //    Input, int FACE_NUM, the number of faces of the polyhedron.
            //
            //    Input, int NODE[FACE_NUM*ORDER_MAX].  Face I is defined by
            //    the vertices NODE(I,0) through NODE(I,ORDER(I)-1).  These vertices
            //    are listed in neighboring order.
            //
            //    Input, int NODE_NUM, the number of points stored in COORD.
            //
            //    Input, int ORDER[FACE_NUM], the number of vertices making up each face.
            //
            //    Output, double POLYHEDRON_AREA_3D, the surface area of the polyhedron.
            //
        {
            int DIM_NUM = 3;

            double ainc;
            double area;
            int face;
            int j;
            int k;
            double[] p1 = new double[DIM_NUM];
            double[] p2 = new double[DIM_NUM];
            double[] p3;
            double[] p4 = new double[DIM_NUM];

            area = 0.0;
            //
            //  For each face
            //
            for (face = 0; face < face_num; face++)
            {
                typeMethods.r8vec_zero(DIM_NUM, ref p4);
                //
                //  For each triangle in the face, compute the normal vector.
                //
                for (j = 0; j < order[face]; j++)
                {
                    k = node[j + face * order_max];
                    p1[0] = coord[0 + k * 3];
                    p1[1] = coord[1 + k * 3];
                    p1[2] = coord[2 + k * 3];

                    if (j + 1 < order[face])
                    {
                        k = node[j + 1 + face * order_max];
                    }
                    else
                    {
                        k = node[0 + face * order_max];
                    }

                    p2[0] = coord[0 + k * 3];
                    p2[1] = coord[1 + k * 3];
                    p2[2] = coord[2 + k * 3];

                    p3 = typeMethods.r8vec_cross_product_3d(p1, p2);

                    p4[0] = p4[0] + p3[0];
                    p4[1] = p4[1] + p3[1];
                    p4[2] = p4[2] + p3[2];

                }

                //
                //  Add the magnitude of the normal vector to the sum.
                //
                ainc = typeMethods.r8vec_norm(DIM_NUM, p4);

                area = area + ainc;
            }

            area = 0.5 * area;

            return area;
        }

        public static double[] polyhedron_centroid_3d(double[] coord, int order_max, int face_num,
                int[] node, int node_num, int[] order)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYHEDRON_CENTROID_3D computes the centroid of a polyhedron in 3D.
            //
            //  Discussion:
            //
            //    The centroid can be computed as the volume-weighted average of
            //    the centroids of the tetrahedra defined by choosing a point in
            //    the interior of the polyhedron, and using as a base every triangle
            //    created by triangulating the faces of the polyhedron.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 August 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double COORD[NODE_NUM*3], the 3D coordinates of the vertices.
            //    The vertices may be listed in any order.
            //
            //    Input, int ORDER_MAX, the maximum number of vertices that make
            //    up a face of the polyhedron.
            //
            //    Input, int FACE_NUM, the number of faces of the polyhedron.
            //
            //    Input, int NODE[FACE_NUM*ORDER_MAX].  Face I is defined by
            //    the vertices NODE(I,1) through NODE(I,ORDER(I)).  These vertices
            //    are listed in neighboring order.
            //
            //    Input, int NODE_NUM, the number of points stored in COORD.
            //
            //    Input, int ORDER[FACE_NUM], the number of vertices making up
            //    each face.
            //
            //    Output, double POLYHEDRON_CENTROID_3D[3], the centroid of the polyhedron.
            //
        {
            int DIM_NUM = 3;

            double area;
            double[] centroid;
            int face;
            int i;
            int j;
            int n;
            int n1;
            int n2;
            int n3;
            double[] normal = new double[DIM_NUM];
            double[] point = new double[DIM_NUM];
            double polygon_area;
            double[] polygon_centroid;
            double[] tet = new double[DIM_NUM * 4];
            double[] tetrahedron_centroid;
            double tetrahedron_volume;
            int v;
            double[] vert;
            int vert_num;
            double volume;
            //
            //  Compute a point in the interior.
            //  We take the area-weighted centroid of each face.
            //
            typeMethods.r8vec_zero(DIM_NUM, ref point);

            vert = new double[DIM_NUM * order_max];

            area = 0.0;

            for (face = 0; face < face_num; face++)
            {
                vert_num = order[face];

                for (j = 0; j < vert_num; j++)
                {
                    n = node[j + face * order_max];

                    vert[0 + j * DIM_NUM] = coord[0 + (n - 1) * DIM_NUM];
                    vert[1 + j * DIM_NUM] = coord[1 + (n - 1) * DIM_NUM];
                    vert[2 + j * DIM_NUM] = coord[2 + (n - 1) * DIM_NUM];
                }

                polygon_area = Burkardt.Polygon.Geometry.polygon_area_3d(vert_num, vert, ref normal);

                polygon_centroid = Burkardt.Polygon.Geometry.polygon_centroid_3d(vert_num, vert);

                for (i = 0; i < DIM_NUM; i++)
                {
                    point[i] = point[i] + polygon_area * polygon_centroid[i];
                }

                area = area + polygon_area;

            }


            point[0] = point[0] / area;
            point[1] = point[1] / area;
            point[2] = point[2] / area;
            //
            //  Now triangulate each face.
            //  For each triangle, consider the tetrahedron created by including POINT.
            //
            centroid = new double[DIM_NUM];

            typeMethods.r8vec_zero(DIM_NUM, ref centroid);

            volume = 0.0;

            for (face = 0; face < face_num; face++)
            {
                n3 = node[order[face] - 1 + face * order_max];

                typeMethods.r8vec_copy(DIM_NUM, coord, ref tet, a1index: +(n3 - 1) * 3, a2index: +2 * 3);

                for (v = 0; v < order[face] - 2; v++)
                {
                    n1 = node[v + face * order_max];
                    n2 = node[v + 1 + face * order_max];

                    typeMethods.r8vec_copy(DIM_NUM, coord, ref tet, a1index: +(n1 - 1) * 3, a2index: +0 * 3);
                    typeMethods.r8vec_copy(DIM_NUM, coord, ref tet, a1index: +(n2 - 1) * 3, a2index: +1 * 3);
                    typeMethods.r8vec_copy(DIM_NUM, point, ref tet, a2index: +3 * 3);

                    tetrahedron_volume = Burkardt.TetrahedronNS.Geometry.tetrahedron_volume_3d(tet);

                    tetrahedron_centroid = Burkardt.TetrahedronNS.Geometry.tetrahedron_centroid_3d(tet);

                    for (i = 0; i < DIM_NUM; i++)
                    {
                        centroid[i] = centroid[i]
                                      + tetrahedron_volume * tetrahedron_centroid[i];
                    }

                    volume = volume + tetrahedron_volume;

                }
            }

            for (i = 0; i < 3; i++)
            {
                centroid[i] = centroid[i] / volume;
            }

            return centroid;
        }

        public static bool polyhedron_contains_point_3d(int node_num, int face_num,
                int face_order_max, double[] v, int[] face_order, int[] face_point,
                double[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYHEDRON_CONTAINS_POINT_3D determines if a point is inside a polyhedron.
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
            //  Reference:
            //
            //    Paulo Cezar Pinto Carvalho, Paulo Roma Cavalcanti,
            //    Point in Polyhedron Testing Using Spherical Polygons,
            //    in Graphics Gems V,
            //    edited by Alan Paeth,
            //    Academic Press, 1995, T385.G6975.
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of vertices.
            //
            //    Input, int FACE_NUM, the number of faces.
            //
            //    Input, int FACE_ORDER_MAX, the maximum order of any face.
            //
            //    Input, double V[3*NODE_NUM], the coordinates of the vertices.
            //
            //    Input, int FACE_ORDER[FACE_NUM], the order of each face.
            //
            //    Input, int FACE_POINT[FACE_ORDER_MAX*FACE_NUM], the indices of the
            //    nodes that make up each face.
            //
            //    Input, double P[3], the point to be tested.
            //
            //    Output, bool POLYHEDRON_CONTAINS_POINT_3D, is true if the point
            //    is inside the polyhedron.
            //
        {
            int DIM_NUM = 3;

            double area;
            int face;
            int i;
            bool inside;
            int k;
            int node;
            int node_num_face;
            double[] v_face;

            v_face = new double[DIM_NUM * face_order_max];

            area = 0.0;

            for (face = 0; face < face_num; face++)
            {
                node_num_face = face_order[face];

                for (k = 0; k < node_num_face; k++)
                {
                    node = face_point[k + face * face_order_max];

                    for (i = 0; i < DIM_NUM; i++)
                    {
                        v_face[i + k * DIM_NUM] = v[i + (node - 1) * DIM_NUM];
                    }
                }

                area = area + Burkardt.Polygon.Geometry.polygon_solid_angle_3d(node_num_face, v_face, p);
            }

            //
            //  AREA should be -4*PI, 0, or 4*PI.
            //  So this test should be quite safe!
            //
            if (area < -2.0 * Math.PI || 2.0 * Math.PI < area)
            {
                inside = true;
            }
            else
            {
                inside = false;
            }


            return inside;
        }

        public static double polyhedron_volume_3d(double[] coord, int order_max, int face_num,
                int[] node, int node_num, int[] order)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYHEDRON_VOLUME_3D computes the volume of a polyhedron in 3D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    21 August 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double COORD[NODE_NUM*3], the 3D coordinates of the vertices.
            //    The vertices may be listed in any order.
            //
            //    Input, int ORDER_MAX, the maximum number of vertices that make
            //    up a face of the polyhedron.
            //
            //    Input, int FACE_NUM, the number of faces of the polyhedron.
            //
            //    Input, int NODE[FACE_NUM*ORDER_MAX].  Face I is defined by
            //    the vertices NODE(I,1) through NODE(I,ORDER(I)).  These vertices
            //    are listed in neighboring order.
            //
            //    Input, int NODE_NUM, the number of points stored in COORD.
            //
            //    Input, int ORDER[FACE_NUM], the number of vertices making up
            //    each face.
            //
            //    Output, double POLYHEDRON_VOLUME_3D, the volume of the polyhedron.
            //
        {
            int DIM_NUM = 3;

            int face;
            int n1;
            int n2;
            int n3;
            double term;
            int v;
            double volume;
            double x1;
            double x2;
            double x3;
            double y1;
            double y2;
            double y3;
            double z1;
            double z2;
            double z3;
            //
            volume = 0.0;
            //
            //  Triangulate each face.
            //
            for (face = 0; face < face_num; face++)
            {
                n3 = node[order[face] - 1 + face * order_max];
                x3 = coord[0 + n3 * 3];
                y3 = coord[1 + n3 * 3];
                z3 = coord[2 + n3 * 3];

                for (v = 0; v < order[face] - 2; v++)
                {
                    n1 = node[v + face * order_max];
                    x1 = coord[0 + n1 * 3];
                    y1 = coord[1 + n1 * 3];
                    z1 = coord[2 + n1 * 3];

                    n2 = node[v + 1 + face * order_max];
                    x2 = coord[0 + n2 * 3];
                    y2 = coord[1 + n2 * 3];
                    z2 = coord[2 + n2 * 3];

                    term = x1 * y2 * z3 - x1 * y3 * z2
                        + x2 * y3 * z1 - x2 * y1 * z3
                        + x3 * y1 * z2 - x3 * y2 * z1;

                    volume = volume + term;
                }

            }

            volume = volume / 6.0;

            return volume;
        }

        public static double polyhedron_volume_3d_2(double[] coord, int order_max, int face_num,
                int[] node, int node_num, int[] order)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYHEDRON_VOLUME_3D_2 computes the volume of a polyhedron in 3D.
            //
            //  Discussion:
            //
            //    The computation is not valid unless the faces of the polyhedron
            //    are planar polygons.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 July 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Allen Van Gelder,
            //    Efficient Computation of Polygon Area and Polyhedron Volume,
            //    Graphics Gems V, edited by Alan Paeth,
            //    AP Professional, 1995, T385.G6975.
            //
            //  Parameters:
            //
            //    Input, double COORD[3*NODE_NUM], the 3D coordinates of the vertices.
            //    The vertices may be listed in any order.
            //
            //    Input, int ORDER_MAX, the maximum number of vertices that make
            //    up a face of the polyhedron.
            //
            //    Input, int FACE_NUM, the number of faces of the polyhedron.
            //
            //    Input, int NODE[FACE_NUM*ORDER_MAX].  Face I is defined by
            //    the vertices NODE(I,1) through NODE(I,ORDER(I)).  These vertices
            //    are listed in neighboring order.
            //
            //    Input, int NODE_NUM, the number of points stored in COORD.
            //
            //    Input, int ORDER[FACE_NUM], the number of vertices making up
            //    each face.
            //
            //    Output, double POLYHEDRON_VOLUME_3D_2, the volume of the polyhedron.
            //
        {
            int DIM_NUM = 3;

            int face;
            int j;
            int k;
            double volume;
            double[] v1 = new double[DIM_NUM];
            double[] v2 = new double[DIM_NUM];
            double[] v3;
            double[] v4 = new double[DIM_NUM];

            volume = 0.0;

            for (face = 0; face < face_num; face++)
            {
                typeMethods.r8vec_zero(DIM_NUM, ref v4);
                //
                //  Compute the area vector for this face.
                //
                for (j = 0; j < order[face]; j++)
                {
                    k = node[j + face * order_max];
                    v1[0] = coord[0 + k * 3];
                    v1[1] = coord[1 + k * 3];
                    v1[2] = coord[2 + k * 3];

                    if (j + 1 < order[face])
                    {
                        k = node[j + 1 + face * order_max];
                    }
                    else
                    {
                        k = node[0 + face * order_max];
                    }

                    v2[0] = coord[0 + k * 3];
                    v2[1] = coord[1 + k * 3];
                    v2[2] = coord[2 + k * 3];

                    v3 = typeMethods.r8vec_cross_product_3d(v1, v2);

                    v4[0] = v4[0] + v3[0];
                    v4[1] = v4[1] + v3[1];
                    v4[2] = v4[2] + v3[2];

                }

                //
                //  Area vector dot any vertex.
                //
                k = node[0 + face * order_max];

                volume = volume + v4[0] * coord[0 + k * 3]
                                + v4[1] * coord[1 + k * 3]
                                + v4[2] * coord[2 + k * 3];

            }

            volume = volume / 6.0;

            return volume;
        }

    }
}