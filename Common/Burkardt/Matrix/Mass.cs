using System;
using Burkardt.FEM;
using Triangle = Burkardt.TriangleNS.Triangle;

namespace Burkardt.MatrixNS
{
    public class Mass
    {
        public static double[] mass_matrix_t3(int node_num, int element_num, int[] element_node,
        double[] node_xy )

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    MASS_MATRIX_T3 computes the mass matrix, using 3-node triangles.
        //
        //  Discussion:
        //
        //    The mass matrix to be estimated has the form:
        //
        //      A(I,J) = integral ( PHI(I)(X,Y) * PHI(J)(X,Y) ) d Region
        //
        //    where PHI(I) and PHI(J) are the shape functions associated with
        //    the I-th and J-th variables.
        //
        //  Element T3:
        //
        //    |
        //    1  3
        //    |  ..
        //    |  . .
        //    S  .  .
        //    |  .   .
        //    |  .    .
        //    0  1-----2
        //    |
        //    +--0--R--1-->
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
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        //
        //    Input, int ELEMENT_NODE[3*ELEMENT_NUM], the nodes that make up each element.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the nodes.
        //
        //    Output, double MASS_MATRIX_T3= new double[NODE_NUM*NODE_NUM], the mass matrix.
        //
        {
            int ELEMENT_ORDER = 3;

            double[] a;
            double area;
            double[] dwdr = new double[ELEMENT_ORDER];
            double[] dwds = new double[ELEMENT_ORDER];
            int element;
            int global1;
            int global2;
            int local1;
            int local2;
            int p1;
            int p2;
            int p3;
            int quad;
            int quad_num;
            double r;
            double[] rtab;
            int rule;
            double s;
            double[] stab;
            double[] w = new double[ELEMENT_ORDER];
            double[] weight;

            a = new double[node_num * node_num];
            //
            //  Zero out the matrix.
            //
            for (global2 = 0; global2 < node_num; global2++)
            {
                for (global1 = 0; global1 < node_num; global1++)
                {
                    a[global1 + global2 * node_num] = 0.0;
                }
            }

            //
            //  Get the weights and abscissas for a unit triangle.
            //
            rule = 4;
            quad_num = Triangle.triangle_unit_size(rule);

            rtab = new double[quad_num];
            stab = new double[quad_num];
            weight = new double[quad_num];

            Triangle.triangle_unit_set(rule, rtab, stab, weight);
            //
            //  For each element.
            //
            for (element = 0; element < element_num; element++)
            {
                p1 = element_node[0 + element * 3] - 1;
                p2 = element_node[1 + element * 3] - 1;
                p3 = element_node[2 + element * 3] - 1;

                area = 0.5 * Math.Abs(
                    node_xy[0 + p1 * 2] * (node_xy[1 + p2 * 2] - node_xy[1 + p3 * 2])
                    + node_xy[0 + p2 * 2] * (node_xy[1 + p3 * 2] - node_xy[1 + p1 * 2])
                    + node_xy[0 + p3 * 2] * (node_xy[1 + p1 * 2] - node_xy[1 + p2 * 2]));

                if (area == 0.0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("MASS_MATRIX_T3 - Fatal error!");
                    Console.WriteLine("  Zero area for element " + element + "");
                    Console.WriteLine("  Node 1 = " + p1 + "");
                    Console.WriteLine("  X = " + node_xy[0 + p1 * 2] + "");
                    Console.WriteLine("  Y = " + node_xy[1 + p1 * 2] + "");
                    Console.WriteLine("  Node 2 = " + p2 + "");
                    Console.WriteLine("  X = " + node_xy[0 + p2 * 2] + "");
                    Console.WriteLine("  Y = " + node_xy[1 + p2 * 2] + "");
                    Console.WriteLine("  Node 3 = " + p3 + "");
                    Console.WriteLine("  X = " + node_xy[0 + p3 * 2] + "");
                    Console.WriteLine("  Y = " + node_xy[1 + p3 * 2] + "");
                    return null;
                }

                //
                //  For each quadrature point in the element...
                //
                for (quad = 0; quad < quad_num; quad++)
                {
                    r = rtab[quad];
                    s = stab[quad];

                    Shape.shape_t3(r, s, w, ref dwdr, ref dwds);
                    //
                    //  For each basis function PHI(I) associated with a node in the element,
                    //
                    for (local1 = 0; local1 < 3; local1++)
                    {
                        global1 = element_node[local1 + element * 3] - 1;
                        //
                        //  For each "neighbor" basis function PHI(J) associated with a node in
                        //  the element.
                        //
                        for (local2 = 0; local2 < 3; local2++)
                        {
                            global2 = element_node[local2 + element * 3] - 1;

                            a[global1 + global2 * node_num] = a[global1 + global2 * node_num]
                                                              + area * weight[quad] * w[local1] * w[local2];
                        }
                    }
                }
            }

            return a;
        }

        public static double[] mass_matrix_t6(int node_num, int element_num, int[] element_node,
        double[] node_xy )

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    MASS_MATRIX_T6 computes the mass matrix, using 6-node triangles.
        //
        //  Discussion:
        //
        //    The mass matrix to be estimated has the form:
        //
        //      A(I,J) = integral ( PHI(I)(X,Y) * PHI(J)(X,Y) ) d Region
        //
        //    where PHI(I) and PHI(J) are the shape functions associated with
        //    the I-th and J-th variables.
        //
        //  Element T6:
        //
        //    |
        //    1  3
        //    |  ..
        //    |  . .
        //    S  6  5
        //    |  .   .
        //    |  .    .
        //    0  1--4--2
        //    |
        //    +--0--R--1-->
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        //
        //    Input, int ELEMENT_NODE[6*ELEMENT_NUM], the nodes that make up each element.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the nodes.
        //
        //    Output, double MASS_MATRIX_T6= new double[NODE_NUM*NODE_NUM], the mass matrix.
        //
        {
            int ELEMENT_ORDER = 6;

            double[] a;
            double area;
            double[] dwdr = new double[ELEMENT_ORDER];
            double[] dwds = new double[ELEMENT_ORDER];
            int element;
            int global1;
            int global2;
            int local1;
            int local2;
            int p1;
            int p2;
            int p3;
            int quad;
            int quad_num;
            double r;
            double[] rtab;
            int rule;
            double s;
            double[] stab;
            double[] w = new double[ELEMENT_ORDER];
            double[] weight;

            a = new double[node_num * node_num];
            //
            //  Zero out the matrix.
            //
            for (global2 = 0; global2 < node_num; global2++)
            {
                for (global1 = 0; global1 < node_num; global1++)
                {
                    a[global1 + global2 * node_num] = 0.0;
                }
            }

            //
            //  Get the weights and abscissas for a unit triangle.
            //
            rule = 12;
            quad_num = Triangle.triangle_unit_size(rule);

            rtab = new double[quad_num];
            stab = new double[quad_num];
            weight = new double[quad_num];

            Triangle.triangle_unit_set(rule, rtab, stab, weight);
            //
            //  For each element.
            //
            for (element = 0; element < element_num; element++)
            {
                p1 = element_node[0 + element * 6] - 1;
                p2 = element_node[1 + element * 6] - 1;
                p3 = element_node[2 + element * 6] - 1;

                area = 0.5 * Math.Abs(
                    node_xy[0 + p1 * 2] * (node_xy[1 + p2 * 2] - node_xy[1 + p3 * 2])
                    + node_xy[0 + p2 * 2] * (node_xy[1 + p3 * 2] - node_xy[1 + p1 * 2])
                    + node_xy[0 + p3 * 2] * (node_xy[1 + p1 * 2] - node_xy[1 + p2 * 2]));

                if (area == 0.0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("MASS_MATRIX_T6 - Fatal error!");
                    Console.WriteLine("  Zero area for element " + element + "");
                    Console.WriteLine("  Node 1 = " + p1 + "");
                    Console.WriteLine("  X = " + node_xy[0 + p1 * 2] + "");
                    Console.WriteLine("  Y = " + node_xy[1 + p1 * 2] + "");
                    Console.WriteLine("  Node 2 = " + p2 + "");
                    Console.WriteLine("  X = " + node_xy[0 + p2 * 2] + "");
                    Console.WriteLine("  Y = " + node_xy[1 + p2 * 2] + "");
                    Console.WriteLine("  Node 3 = " + p3 + "");
                    Console.WriteLine("  X = " + node_xy[0 + p3 * 2] + "");
                    Console.WriteLine("  Y = " + node_xy[1 + p3 * 2] + "");
                    return null;
                }

                //
                //  For each quadrature point in the element...
                //
                for (quad = 0; quad < quad_num; quad++)
                {
                    r = rtab[quad];
                    s = stab[quad];

                    Shape.shape_t6(r, s, w, ref dwdr, ref dwds);
                    //
                    //  For each basis function PHI(I) associated with a node in the element,
                    //
                    for (local1 = 0; local1 < 6; local1++)
                    {
                        global1 = element_node[local1 + element * 6] - 1;
                        //
                        //  For each "neighbor" basis function PHI(J) associated with a node in
                        //  the element.
                        //
                        for (local2 = 0; local2 < 6; local2++)
                        {
                            global2 = element_node[local2 + element * 6] - 1;

                            a[global1 + global2 * node_num] = a[global1 + global2 * node_num]
                                                              + area * weight[quad] * w[local1] * w[local2];
                        }
                    }
                }
            }

            return a;
        }
    }
}