﻿using System;
using Burkardt.Types;

namespace Burkardt.SphereNS;

public static class Icosphere
{
    public static int sphere_icos_edge_num(int factor)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_ICOS_EDGE_NUM sizes an icosahedral grid on a sphere.
        //
        //  Discussion:
        //
        //    With FACTOR = 1, the grid has 20 triangular faces, 30 edges, and 12 nodes.
        //
        //    With FACTOR = 2, each triangle of the icosahedron is subdivided into
        //    2x2 subtriangles, resulting in 80 faces, 120 edges, and 
        //    42 = 12 + 20 * 3 * (1)/2 + 20 * 0 ) nodes.
        //
        //    With FACTOR = 3, each triangle of the icosahedron is subdivided into
        //    3x3 subtriangles, resulting in 180 faces, 270 edges and 
        //    92 ( = 12 + 20 * 3 * (2)/2 + 20 * 1 ) nodes.
        //
        //    In general, each triangle is subdivided into FACTOR*FACTOR subtriangles,
        //    resulting in 20 * FACTOR * FACTOR faces, 30 * FACTOR * FACTOR edges, and
        //      12 
        //    + 20 * 3          * (FACTOR-1) / 2 
        //    + 20 * (FACTOR-2) * (FACTOR-1) / 2 nodes.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int FACTOR, the subdivision factor, which must
        //    be at least 1.
        //
        //    Output, int SPHERE_ICOS_EDGE_NUM, the number of edges.
        //
    {
        int edge_num = 30 * factor * factor;

        return edge_num;
    }

    public static int sphere_icos_face_num(int factor)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_ICOS_FACE_NUM sizes an icosahedral grid on a sphere.
        //
        //  Discussion:
        //
        //    With FACTOR = 1, the grid has 20 triangular faces, 30 edges, and 12 nodes.
        //
        //    With FACTOR = 2, each triangle of the icosahedron is subdivided into
        //    2x2 subtriangles, resulting in 80 faces, 120 edges, and 
        //    42 = 12 + 20 * 3 * (1)/2 + 20 * 0 ) nodes.
        //
        //    With FACTOR = 3, each triangle of the icosahedron is subdivided into
        //    3x3 subtriangles, resulting in 180 faces, 270 edges and 
        //    92 ( = 12 + 20 * 3 * (2)/2 + 20 * 1 ) nodes.
        //
        //    In general, each triangle is subdivided into FACTOR*FACTOR subtriangles,
        //    resulting in 20 * FACTOR * FACTOR faces, 30 * FACTOR * FACTOR edges, and
        //      12 
        //    + 20 * 3          * (FACTOR-1) / 2 
        //    + 20 * (FACTOR-2) * (FACTOR-1) / 2 nodes.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int FACTOR, the subdivision factor, which must
        //    be at least 1.
        //
        //    Output, int SPHERE_ICOS_FACE_NUM, the number of faces.
        //
    {
        int face_num = 20 * factor * factor;

        return face_num;
    }

    public static int sphere_icos_point_num(int factor)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_ICOS_POINT_NUM sizes an icosahedral grid on a sphere.
        //
        //  Discussion:
        //
        //    With FACTOR = 1, the grid has 20 triangular faces, 30 edges, and 12 nodes.
        //
        //    With FACTOR = 2, each triangle of the icosahedron is subdivided into
        //    2x2 subtriangles, resulting in 80 faces, 120 edges, and 
        //    42 = 12 + 20 * 3 * (1)/2 + 20 * 0 ) nodes.
        //
        //    With FACTOR = 3, each triangle of the icosahedron is subdivided into
        //    3x3 subtriangles, resulting in 180 faces, 270 edges and 
        //    92 ( = 12 + 20 * 3 * (2)/2 + 20 * 1 ) nodes.
        //
        //    In general, each triangle is subdivided into FACTOR*FACTOR subtriangles,
        //    resulting in 20 * FACTOR * FACTOR faces, 30 * FACTOR * FACTOR edges, and
        //      12 
        //    + 20 * 3          * (FACTOR-1) / 2 
        //    + 20 * (FACTOR-2) * (FACTOR-1) / 2 nodes.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int FACTOR, the subdivision factor, which must
        //    be at least 1.
        //
        //    Output, int SPHERE_ICOS_POINT_NUM, the number of points.
        //
    {
        int point_num = 12
                        + 10 * 3 * (factor - 1)
                        + 10 * (factor - 2) * (factor - 1);

        return point_num;
    }

    public static double[] sphere_icos1_points(int factor, int node_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_ICOS1_POINTS returns icosahedral grid points on a sphere.
        //
        //  Discussion:
        //
        //    With FACTOR = 1, the grid has 20 triangular faces and 12 nodes.
        //
        //    With FACTOR = 2, each triangle of the icosahedron is subdivided into
        //    2x2 subtriangles, resulting in 80 faces and 
        //    42 = 12 + 20 * 3 * (1)/2 + 20 * 0 ) nodes.
        //
        //    With FACTOR = 3, each triangle of the icosahedron is subdivided into
        //    3x3 subtriangles, resulting in 180 faces and 
        //    92 ( = 12 + 20 * 3 * (2)/2 + 20 * 1 ) nodes.
        //
        //    In general, each triangle is subdivided into FACTOR*FACTOR subtriangles,
        //    resulting in 20 * FACTOR * FACTOR faces and
        //      12 
        //    + 20 * 3          * (FACTOR-1) / 2 
        //    + 20 * (FACTOR-2) * (FACTOR-1) / 2 nodes.
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
        //    Input, int FACTOR, the subdivision factor, which must
        //    be at least 1.
        //
        //    Input, int NODE_NUM, the number of nodes, as reported
        //    by SPHERE_GRID_ICOS_NUM.
        //
        //    Output, double SPHERE_ICOS1_POINTS[3*NODE_NUM], the node coordinates.
        //
        //  Local Parameters:
        //
        //    POINT_NUM, EDGE_NUM, FACE_NUM and FACE_ORDER_MAX are counters 
        //    associated with the icosahedron, and POINT_COORD, EDGE_POINT, 
        //    FACE_ORDER and FACE_POINT are data associated with the icosahedron.
        //    We need to refer to this data to generate the grid.
        //
        //    NODE counts the number of nodes we have generated so far.  At the
        //    end of the routine, it should be equal to NODE_NUM.
        //
    {
        int a;
        int b;
        int dim;
        int edge;
        int edge_num = 0;
        int face;
        int face_num = 0;
        int face_order_max = 0;
        double node_norm;
        int point;
        int point_num = 0;

        double[] node_xyz = new double[3 * node_num];
        //
        //  Size the icosahedron.
        //
        Icosahedron.Geometry.icos_num(ref point_num, ref edge_num, ref face_num, ref face_order_max);
        //
        //  Set the icosahedron.
        //
        double[] point_coord = new double[3 * point_num];
        int[] edge_point = new int[2 * edge_num];
        int[] face_order = new int[face_num];
        int[] face_point = new int[face_order_max * face_num];

        Icosahedron.Geometry.icos_shape(point_num, edge_num, face_num, face_order_max,
            ref point_coord, ref edge_point, ref face_order, ref face_point);
        //
        //  Generate the point coordinates.
        //
        //  A.  Points that are the icosahedral vertices.
        //
        int node = 0;
        for (point = 0; point < point_num; point++)
        {
            for (dim = 0; dim < 3; dim++)
            {
                node_xyz[dim + node * 3] = point_coord[dim + point * 3];
            }

            node += 1;
        }

        //
        //  B. Points in the icosahedral edges, at 
        //  1/FACTOR, 2/FACTOR, ..., (FACTOR-1)/FACTOR.
        //
        for (edge = 0; edge < edge_num; edge++)
        {
            a = edge_point[0 + edge * 2];
            b = edge_point[1 + edge * 2];

            int f;
            for (f = 1; f < factor; f++)
            {
                for (dim = 0; dim < 3; dim++)
                {
                    node_xyz[dim + node * 3] =
                        ((factor - f) * point_coord[dim + a * 3]
                         + f * point_coord[dim + b * 3])
                        / factor;
                }

                node_norm = typeMethods.r8vec_norm(3, node_xyz, +node * 3);

                for (dim = 0; dim < 3; dim++)
                {
                    node_xyz[dim + node * 3] /= node_norm;
                }

                node += 1;
            }
        }

        //
        //  C.  Points in the icosahedral faces.
        //
        for (face = 0; face < face_num; face++)
        {
            a = face_point[0 + face * 3];
            b = face_point[1 + face * 3];
            int c = face_point[2 + face * 3];

            int f1;
            for (f1 = 1; f1 < factor; f1++)
            {
                int f2;
                for (f2 = 1; f2 < factor - f1; f2++)
                {
                    for (dim = 0; dim < 3; dim++)
                    {
                        node_xyz[dim + node * 3] =
                            ((factor - f1 - f2) * point_coord[dim + a * 3]
                             + f1 * point_coord[dim + b * 3]
                             + f2 * point_coord[dim + c * 3])
                            / factor;
                    }

                    node_norm = typeMethods.r8vec_norm(3, node_xyz, +node * 3);

                    for (dim = 0; dim < 3; dim++)
                    {
                        node_xyz[dim + node * 3] /= node_norm;
                    }

                    node += 1;
                }
            }
        }

        return node_xyz;
    }

    public static double[] sphere_icos2_points(int factor, int node_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_ICOS2_POINTS returns icosahedral grid points on a sphere.
        //
        //  Discussion:
        //
        //    With FACTOR = 1, the grid has 20 triangular faces and 12 nodes.
        //
        //    With FACTOR = 2, each triangle of the icosahedron is subdivided into
        //    2x2 subtriangles, resulting in 80 faces and 
        //    42 = 12 + 20 * 3 * (1)/2 + 20 * 0 ) nodes.
        //
        //    With FACTOR = 3, each triangle of the icosahedron is subdivided into
        //    3x3 subtriangles, resulting in 180 faces and 
        //    92 ( = 12 + 20 * 3 * (2)/2 + 20 * 1 ) nodes.
        //
        //    In general, each triangle is subdivided into FACTOR*FACTOR subtriangles,
        //    resulting in 20 * FACTOR * FACTOR faces and
        //      12 
        //    + 20 * 3          * (FACTOR-1) / 2 
        //    + 20 * (FACTOR-2) * (FACTOR-1) / 2 nodes.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int FACTOR, the subdivision factor, which must
        //    be at least 1.
        //
        //    Input, int NODE_NUM, the number of nodes, as reported
        //    by SPHERE_IMP_GRID_ICOS_NUM.
        //
        //    Output, double SPHERE_ICOS2_POINTS[3*NODE_NUM], the node coordinates.
        //
        //  Local Parameters:
        //
        //    POINT_NUM, EDGE_NUM, FACE_NUM and FACE_ORDER_MAX are counters 
        //    associated with the icosahedron, and POINT_COORD, EDGE_POINT, 
        //    FACE_ORDER and FACE_POINT are data associated with the icosahedron.
        //    We need to refer to this data to generate the grid.
        //
        //    NODE counts the number of nodes we have generated so far.  At the
        //    end of the routine, it should be equal to NODE_NUM.
        //
    {
        int a;
        double angle;
        double[] ab = new double[3];
        double[] ac = new double[3];
        double[] acn = new double[3];
        double[] acp = new double[3];
        int b;
        double[] bn = new double[3];
        double bn_norm;
        double[] bp = new double[3];
        double[] cn = new double[3];
        double[] cp = new double[3];
        int edge;
        int edge_num = 0;
        int face;
        int face_num = 0;
        int face_order_max = 0;
        int i;
        int j;
        int point_num = 0;

        double[] node_xyz = new double[3 * node_num];
        //
        //  Size the icosahedron.
        //
        Icosahedron.Geometry.icos_num(ref point_num, ref edge_num, ref face_num, ref face_order_max);
        //
        //  Set the icosahedron.
        //
        double[] point_coord = new double[3 * point_num];
        int[] edge_point = new int[2 * edge_num];
        int[] face_order = new int[face_num];
        int[] face_point = new int[face_order_max * face_num];

        Icosahedron.Geometry.icos_shape(point_num, edge_num, face_num, face_order_max,
            ref point_coord, ref edge_point, ref face_order, ref face_point);
        //
        //  Generate the point coordinates.
        //
        //  A.  Points that are the icosahedral vertices.
        //
        int node = 0;
        for (j = 0; j < point_num; j++)
        {
            for (i = 0; i < 3; i++)
            {
                node_xyz[i + j * 3] = point_coord[i + j * 3];
            }

            node += 1;
        }

        //
        //  B. Points in the icosahedral edges, at 
        //  1/FACTOR, 2/FACTOR, ..., (FACTOR-1)/FACTOR.
        //
        for (edge = 0; edge < edge_num; edge++)
        {
            a = edge_point[0 + edge * 2];

            b = edge_point[1 + edge * 2];
            //
            //  Determine the "distance" = angle between points A and B.
            //
            double theta = Distance.sphere_distance_xyz(point_coord, point_coord, index1: +a * 3, index2: +b * 3);
            //
            //  Polarize B into BP + BN and normalize BN.
            //
            typeMethods.r8vec_polarize(3, point_coord, point_coord, ref bn, ref bp, aIndex: +b * 3, pIndex: +a * 3);
            bn_norm = typeMethods.r8vec_norm(3, bn);
            for (i = 0; i < 3; i++)
            {
                bn[i] /= bn_norm;
            }

            //
            //  March from A to B, by taking equally spaced angles from 0 to THETA.
            //  F = 0      => ANGLE = 0     => A
            //  F = FACTOR => ANGLE = THETA => B
            //
            int f;
            for (f = 1; f < factor; f++)
            {
                angle = f * theta / factor;

                for (i = 0; i < 3; i++)
                {
                    node_xyz[i + node * 3] = Math.Cos(angle) * point_coord[i + a * 3]
                                             + Math.Sin(angle) * bn[i];
                }

                node += 1;
            }
        }

        //
        //  C.  Points in the icosahedral faces.
        //
        for (face = 0; face < face_num; face++)
        {
            a = face_point[0 + face * 3];
            b = face_point[1 + face * 3];
            int c = face_point[2 + face * 3];
            //
            //  Determine the "distance" = angle between points A and B, A and C.
            //
            double theta_ab = Distance.sphere_distance_xyz(point_coord, point_coord, +a * 3, +b * 3);
            double theta_ac = Distance.sphere_distance_xyz(point_coord, point_coord, +a * 3, +c * 3);
            //
            //  Polarize B = BP + BN and normalize BN, C = CP + CN, and normalize CN.
            //
            typeMethods.r8vec_polarize(3, point_coord, point_coord, ref bn, ref bp, +b * 3, +a * 3);
            bn_norm = typeMethods.r8vec_norm(3, bn);
            for (i = 0; i < 3; i++)
            {
                bn[i] /= bn_norm;
            }

            typeMethods.r8vec_polarize(3, point_coord, point_coord, ref cn, ref cp, +c * 3, +a * 3);
            double cn_norm = typeMethods.r8vec_norm(3, cn);
            for (i = 0; i < 3; i++)
            {
                cn[i] /= cn_norm;
            }

            //
            //  March AB from A to B:
            //    FA = 0      => ANGLE = 0        => AB = A
            //    FA = FACTOR => ANGLE = THETA_AB => AB = B
            //
            //  March AC from A to C:
            //    FA = 0      => ANGLE = 0        => AC = A
            //    FA = FACTOR => ANGLE = THETA_AC => AC = C
            //
            int fa;
            for (fa = 2; fa < factor; fa++)
            {
                //
                //  Determine points AB and AC that use Math.Cos ( FA / FACTOR ) of A 
                //  and Math.Cos ( ( FACTOR - FA ) / FACTOR ) of B or C.
                //
                angle = fa * theta_ab / factor;
                for (i = 0; i < 3; i++)
                {
                    ab[i] = Math.Cos(angle) * point_coord[i + a * 3] + Math.Sin(angle) * bn[i];
                }

                angle = fa * theta_ac / factor;
                for (i = 0; i < 3; i++)
                {
                    ac[i] = Math.Cos(angle) * point_coord[i + a * 3] + Math.Sin(angle) * cn[i];
                }

                //
                //  Determine the "distance" = angle between points AB and AC.
                //
                double theta_bc = Distance.sphere_distance_xyz(ab, ac);
                //
                //  Polarize AC into ACP + ACN and normalize ACN.
                //
                typeMethods.r8vec_polarize(3, ac, ab, ref acn, ref acp);
                double acn_norm = typeMethods.r8vec_norm(3, acn);
                for (i = 0; i < 3; i++)
                {
                    acn[i] /= acn_norm;
                }

                //
                //  The interval between AB and AC is broken into FA intervals.
                //  Go from 1 to FA - 1.
                //
                int fbc;
                for (fbc = 1; fbc < fa; fbc++)
                {
                    angle = fbc * theta_bc / fa;
                    for (i = 0; i < 3; i++)
                    {
                        node_xyz[i + node * 3] = Math.Cos(angle) * ab[i] + Math.Sin(angle) * acn[i];
                    }

                    node += 1;
                }
            }
        }

        return node_xyz;
    }

}