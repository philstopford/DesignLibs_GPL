using System;

namespace Burkardt.FEM
{
    public class Basis_mn
    {
        public static void basis_mn_q4(double[] q, int n, double[] p, ref double[] phi,
                ref double[] dphidx, ref double[] dphidy)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BASIS_MN_Q4: all bases at N points for a Q4 element.
            //
            //  Discussion:
            //
            //    The routine is given the coordinates of the vertices of a quadrilateral.
            //    It works directly with these coordinates, and does not refer to a 
            //    reference element.
            //
            //    The sides of the element are presumed to lie along coordinate axes.
            //
            //    The routine evaluates the basis functions associated with each corner,
            //    and their derivatives with respect to X and Y.
            //
            //  Physical Element Q4:
            //
            //    |
            //    |  4-----3
            //    |  |     |
            //    |  |     |
            //    Y  |     |
            //    |  |     |
            //    |  |     |
            //    |  1-----2
            //    |
            //    +-----X------>
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
            //  Parameters:
            //
            //    Input, double Q[2*4], the coordinates of the vertices.
            //    It is common to list these points in counter clockwise order.
            //
            //    Input, int N, the number of evaluation points.
            //
            //    Input, double P[2*N], the evaluation points.
            //
            //    Output, double PHI[4*N], the bases at the evaluation points.
            //
            //    Output, double DPHIDX[4*N], DPHIDY[4*N], the derivatives of the
            //    bases at the evaluation points.
            //
        {
            double area;
            int i;
            int j;

            area = (q[0 + 2 * 2] - q[0 + 0 * 2])
                   * (q[1 + 2 * 2] - q[1 + 0 * 2]);

            for (j = 0; j < n; j++)
            {
                phi[0 + j * 4] = (q[0 + 2 * 2] - p[0 + j * 2])
                                 * (q[1 + 2 * 2] - p[1 + j * 2]);
                phi[1 + j * 4] = (p[0 + j * 2] - q[0 + 0 * 2])
                                 * (q[1 + 2 * 2] - p[1 + j * 2]);
                phi[2 + j * 4] = (p[0 + j * 2] - q[0 + 0 * 2])
                                 * (p[1 + j * 2] - q[1 + 0 * 2]);
                phi[3 + j * 4] = (q[0 + 2 * 2] - p[0 + j * 2])
                                 * (p[1 + j * 2] - q[1 + 0 * 2]);

                dphidx[0 + j * 4] = -(q[1 + 2 * 2] - p[1 + j * 2]);
                dphidx[1 + j * 4] = (q[1 + 2 * 2] - p[1 + j * 2]);
                dphidx[2 + j * 4] = (p[1 + j * 2] - q[1 + 0 * 2]);
                dphidx[3 + j * 4] = -(p[1 + j * 2] - q[1 + 0 * 2]);

                dphidy[0 + j * 4] = -(q[0 + 2 * 2] - p[0 + j * 2]);
                dphidy[1 + j * 4] = -(p[0 + j * 2] - q[0 + 0 * 2]);
                dphidy[2 + j * 4] = (p[0 + j * 2] - q[0 + 0 * 2]);
                dphidy[3 + j * 4] = (q[0 + 2 * 2] - p[0 + j * 2]);
            }

            //
            //  Normalize.
            //
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < 4; i++)
                {
                    phi[i + j * 4] = phi[i + j * 4] / area;
                    dphidx[i + j * 4] = dphidx[i + j * 4] / area;
                    dphidy[i + j * 4] = dphidy[i + j * 4] / area;
                }
            }

            return;
        }

        public static void basis_mn_q4_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BASIS_MN_Q4_TEST verifies BASIS_MN_Q4.
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
            //  Parameters:
            //
            //    None
            //
        {
            int NODE_NUM = 4;

            double[] dphidx = new double[NODE_NUM * NODE_NUM];
            double[] dphidy = new double[NODE_NUM * NODE_NUM];
            int i;
            int j;
            double[] phi = new double[NODE_NUM * NODE_NUM];
            double[] q =
                {
                    3.0, 1.0,
                    5.0, 1.0,
                    5.0, 4.0,
                    3.0, 4.0
                }
                ;
            double sum_x;
            double sum_y;

            Console.WriteLine("");
            Console.WriteLine("BASIS_MN_Q4_TEST:");
            Console.WriteLine("  Verify basis functions for element Q4.");
            Console.WriteLine("");
            Console.WriteLine("  Number of nodes = " + NODE_NUM + "");

            Console.WriteLine("");
            Console.WriteLine("  Physical Nodes:");
            Console.WriteLine("");
            for (i = 0; i < NODE_NUM; i++)
            {
                Console.WriteLine("  "
                                  + q[0 + i * 2].ToString().PadLeft(10) + "  "
                                  + q[1 + i * 2].ToString().PadLeft(10) + "");
            }

            basis_mn_q4(q, NODE_NUM, q, ref phi, ref dphidx, ref dphidy);

            Console.WriteLine("");
            Console.WriteLine("  The basis function values at basis nodes");
            Console.WriteLine("  should form the identity matrix.");
            Console.WriteLine("");

            for (i = 0; i < NODE_NUM; i++)
            {
                string cout = "";
                for (j = 0; j < NODE_NUM; j++)
                {
                    cout += "  " + phi[i + j * NODE_NUM].ToString().PadLeft(10);
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("  The X and Y derivatives should sum to 0.");
            Console.WriteLine("");
            Console.WriteLine("  dPhidX sum, dPhidY sum:");
            Console.WriteLine("");

            for (j = 0; j < NODE_NUM; j++)
            {
                sum_x = 0.0;
                for (i = 0; i < NODE_NUM; i++)
                {
                    sum_x = sum_x + dphidx[i + j * NODE_NUM];
                }

                sum_y = 0.0;
                for (i = 0; i < NODE_NUM; i++)
                {
                    sum_y = sum_y + dphidy[i + j * NODE_NUM];
                }

                Console.WriteLine("  "
                                  + sum_x.ToString().PadLeft(10) + "  "
                                  + sum_y.ToString().PadLeft(10) + "");
            }

        }

        public static void basis_mn_t3(double[] t, int n, double[] p, ref double[] phi,
                ref double[] dphidx, ref double[] dphidy)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BASIS_MN_T3: all bases at N points for a T3 element.
            //
            //  Discussion:
            //
            //    The routine is given the coordinates of the vertices of a triangle.
            //    It works directly with these coordinates, and does not refer to a 
            //    reference element.
            //
            //    The sides of the triangle DO NOT have to lie along a coordinate
            //    axis.
            //
            //    The routine evaluates the basis functions associated with each vertex,
            //    and their derivatives with respect to X and Y.
            //
            //  Physical Element T3: 
            //       
            //            3
            //           . .
            //          .   .
            //         .     .
            //        .       .
            //       1---------2
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 February 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double[] t, the coordinates of the vertices
            //    of the triangle.  It is common to list these points in counter clockwise
            //    order.
            //
            //    Input, int N, the number of evaluation points.
            //
            //    Input, double P[2*N], the points where the basis functions 
            //    are to be evaluated.
            //
            //    Output, double PHI[3*N], the value of the basis functions 
            //    at the evaluation points.
            //
            //    Output, double DPHIDX[3*N], DPHIDY[3*N], the value of the 
            //    derivatives at the evaluation points.
            //
            //  Local parameters:
            //
            //    Local, double AREA, is (twice) the area of the triangle.
            //
        {
            double area;
            int i;
            int j;

            area = t[0 + 0 * 2] * (t[1 + 1 * 2] - t[1 + 2 * 2])
                   + t[0 + 1 * 2] * (t[1 + 2 * 2] - t[1 + 0 * 2])
                   + t[0 + 2 * 2] * (t[1 + 0 * 2] - t[1 + 1 * 2]);

            if (area == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("BASIS_MN_T3 - Fatal error!");
                Console.WriteLine("  Element has zero area.");
                return;
            }

            for (j = 0; j < n; j++)
            {
                phi[0 + j * 3] = ((t[0 + 2 * 2] - t[0 + 1 * 2]) * (p[1 + j * 2] - t[1 + 1 * 2])
                                  - (t[1 + 2 * 2] - t[1 + 1 * 2]) * (p[0 + j * 2] - t[0 + 1 * 2]));
                dphidx[0 + j * 3] = -(t[1 + 2 * 2] - t[1 + 1 * 2]);
                dphidy[0 + j * 3] = (t[0 + 2 * 2] - t[0 + 1 * 2]);

                phi[1 + j * 3] = ((t[0 + 0 * 2] - t[0 + 2 * 2]) * (p[1 + j * 2] - t[1 + 2 * 2])
                                  - (t[1 + 0 * 2] - t[1 + 2 * 2]) * (p[0 + j * 2] - t[0 + 2 * 2]));
                dphidx[1 + j * 3] = -(t[1 + 0 * 2] - t[1 + 2 * 2]);
                dphidy[1 + j * 3] = (t[0 + 0 * 2] - t[0 + 2 * 2]);

                phi[2 + j * 3] = ((t[0 + 1 * 2] - t[0 + 0 * 2]) * (p[1 + j * 2] - t[1 + 0 * 2])
                                  - (t[1 + 1 * 2] - t[1 + 0 * 2]) * (p[0 + j * 2] - t[0 + 0 * 2]));
                dphidx[2 + j * 3] = -(t[1 + 1 * 2] - t[1 + 0 * 2]);
                dphidy[2 + j * 3] = (t[0 + 1 * 2] - t[0 + 0 * 2]);
            }

            //
            //  Normalize.
            //
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < 3; i++)
                {
                    phi[i + j * 3] = phi[i + j * 3] / area;
                    dphidx[i + j * 3] = dphidx[i + j * 3] / area;
                    dphidy[i + j * 3] = dphidy[i + j * 3] / area;
                }
            }
        }

        public static void basis_mn_t3_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BASIS_MN_T3_TEST verifies BASIS_MN_T3.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    10 February 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    None.
            //
        {
            int NODE_NUM = 3;

            double[] dphidx = new double[NODE_NUM * NODE_NUM];
            double[] dphidy = new double[NODE_NUM * NODE_NUM];
            int i;
            int j;
            double[] phi = new double[NODE_NUM * NODE_NUM];
            double sum_x;
            double sum_y;
            double[] t =
            {
                2.0, 0.0,
                4.0, 3.0,
                0.0, 4.0
            };

            Console.WriteLine("");
            Console.WriteLine("BASIS_MN_T3_TEST:");
            Console.WriteLine("  Verify basis functions for element T3.");
            Console.WriteLine("");
            Console.WriteLine("  Number of nodes = " + NODE_NUM + "");

            Console.WriteLine("");
            Console.WriteLine("  Physical Nodes:");
            Console.WriteLine("");
            for (j = 0; j < NODE_NUM; j++)
            {
                Console.WriteLine("  "
                                  + t[0 + j * 2].ToString().PadLeft(10) + "  "
                                  + t[1 + j * 2].ToString().PadLeft(10) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("  The basis function values at basis nodes");
            Console.WriteLine("  should form the identity matrix.");
            Console.WriteLine("");

            basis_mn_t3(t, NODE_NUM, t, ref phi, ref dphidx, ref dphidy);

            for (j = 0; j < NODE_NUM; j++)
            {
                string cout = "";
                for (i = 0; i < NODE_NUM; i++)
                {
                    cout += "  " + phi[i + j * NODE_NUM].ToString().PadLeft(10);
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("  The X and Y derivatives should sum to 0.");
            Console.WriteLine("");
            Console.WriteLine("  dPhidX sum, dPhidY sum:");
            Console.WriteLine("");

            for (j = 0; j < NODE_NUM; j++)
            {
                sum_x = 0.0;
                sum_y = 0.0;
                for (i = 0; i < NODE_NUM; i++)
                {
                    sum_x = sum_x + dphidx[i + j * NODE_NUM];
                    sum_y = sum_y + dphidy[i + j * NODE_NUM];
                }

                Console.WriteLine("  "
                                  + sum_x.ToString().PadLeft(10) + "  "
                                  + sum_y.ToString().PadLeft(10) + "");
            }

        }

        public static void basis_mn_t4(double[] t, int n, double[] p, ref double[] phi,
                ref double[] dphidx, ref double[] dphidy)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BASIS_MN_T4: all bases at N points for a T4 element.
            //
            //  Discussion:
            //
            //    The T4 element is the cubic bubble triangle.
            //
            //    The routine is given the coordinates of the vertices of a triangle.
            //    It works directly with these coordinates, and does not refer to a 
            //    reference element.
            //
            //    The sides of the triangle DO NOT have to lie along a coordinate
            //    axis.
            //
            //    The routine evaluates the basis functions associated with each vertex,
            //    and their derivatives with respect to X and Y.
            //
            //  Physical Element T4: 
            //       
            //            3
            //           . .
            //          .   .
            //         .  4  .
            //        .       .
            //       1---------2
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
            //  Parameters:
            //
            //    Input, double T[2*4], the coordinates of the vertices
            //    of the triangle, and the coordinates of the centroid.  
            //    It is common to list the first three points in counter clockwise
            //    order.
            //
            //    Input, int N, the number of evaluation points.
            //
            //    Input, double P[2*N], the points where the basis functions 
            //    are to be evaluated.
            //
            //    Output, double PHI[4*N], the value of the basis functions 
            //    at the evaluation points.
            //
            //    Output, double DPHIDX[4*N], DPHIDY[4*N], the value of the 
            //    derivatives at the evaluation points.
            //
            //  Local parameters:
            //
            //    Local, double AREA, is (twice) the area of the triangle.
            //
        {
            double area;
            int i;
            int j;

            area = t[0 + 0 * 2] * (t[1 + 1 * 2] - t[1 + 2 * 2])
                   + t[0 + 1 * 2] * (t[1 + 2 * 2] - t[1 + 0 * 2])
                   + t[0 + 2 * 2] * (t[1 + 0 * 2] - t[1 + 1 * 2]);

            for (j = 0; j < n; j++)
            {
                phi[0 + j * 4] = ((t[0 + 2 * 2] - t[0 + 1 * 2]) * (p[1 + j * 2] - t[1 + 1 * 2])
                                  - (t[1 + 2 * 2] - t[1 + 1 * 2]) * (p[0 + j * 2] - t[0 + 1 * 2]));
                dphidx[0 + j * 4] = -(t[1 + 2 * 2] - t[1 + 1 * 2]);
                dphidy[0 + j * 4] = (t[0 + 2 * 2] - t[0 + 1 * 2]);

                phi[1 + j * 4] = ((t[0 + 0 * 2] - t[0 + 2 * 2]) * (p[1 + j * 2] - t[1 + 2 * 2])
                                  - (t[1 + 0 * 2] - t[1 + 2 * 2]) * (p[0 + j * 2] - t[0 + 2 * 2]));
                dphidx[1 + j * 4] = -(t[1 + 0 * 2] - t[1 + 2 * 2]);
                dphidy[1 + j * 4] = (t[0 + 0 * 2] - t[0 + 2 * 2]);

                phi[2 + j * 4] = ((t[0 + 1 * 2] - t[0 + 0 * 2]) * (p[1 + j * 2] - t[1 + 0 * 2])
                                  - (t[1 + 1 * 2] - t[1 + 0 * 2]) * (p[0 + j * 2] - t[0 + 0 * 2]));
                dphidx[2 + j * 4] = -(t[1 + 1 * 2] - t[1 + 0 * 2]);
                dphidy[2 + j * 4] = (t[0 + 1 * 2] - t[0 + 0 * 2]);
            }

            //
            //  Normalize the first three functions.
            //
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < 3; i++)
                {
                    phi[i + j * 4] = phi[i + j * 4] / area;
                    dphidx[i + j * 4] = dphidx[i + j * 4] / area;
                    dphidy[i + j * 4] = dphidy[i + j * 4] / area;
                }
            }

            //
            //  Compute the cubic bubble function.
            //
            for (j = 0; j < n; j++)
            {
                phi[3 + j * 4] = 27.0 * phi[0 + j * 4] * phi[1 + j * 4] * phi[2 + j * 4];

                dphidx[3 + j * 4] = 27.0 * (
                    dphidx[0 + j * 4] * phi[1 + j * 4] * phi[2 + j * 4]
                    + phi[0 + j * 4] * dphidx[1 + j * 4] * phi[2 + j * 4]
                    + phi[0 + j * 4] * phi[1 + j * 4] * dphidx[2 + j * 4]);

                dphidy[3 + j * 4] = 27.0 * (
                    dphidy[0 + j * 4] * phi[1 + j * 4] * phi[2 + j * 4]
                    + phi[0 + j * 4] * dphidy[1 + j * 4] * phi[2 + j * 4]
                    + phi[0 + j * 4] * phi[1 + j * 4] * dphidy[2 + j * 4]);
            }

            //
            //  Subtract 1/3 of the cubic bubble function from each of the three linears.
            //
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < 3; i++)
                {
                    phi[i + j * 4] = phi[i + j * 4] - phi[3 + j * 4] / 3.0;
                    dphidx[i + j * 4] = dphidx[i + j * 4] - dphidx[3 + j * 4] / 3.0;
                    dphidy[i + j * 4] = dphidy[i + j * 4] - dphidy[3 + j * 4] / 3.0;
                }
            }

            return;
        }

        public static void basis_mn_t4_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BASIS_MN_T4_TEST verifies BASIS_MN_T4.
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
            //  Parameters:
            //
        {
            int NODE_NUM = 4;

            double[] dphidx = new double[NODE_NUM * NODE_NUM];
            double[] dphidy = new double[NODE_NUM * NODE_NUM];
            int i;
            int j;
            double[] phi = new double[NODE_NUM * NODE_NUM];
            double sum_x;
            double sum_y;
            double[] t =
            {
                2.0, 0.0,
                4.0, 2.0,
                0.0, 4.0,
                2.0, 2.0
            };

            Console.WriteLine("");
            Console.WriteLine("BASIS_MN_T4_TEST:");
            Console.WriteLine("  Verify basis functions for element T4.");
            Console.WriteLine("");
            Console.WriteLine("  Number of nodes = " + NODE_NUM + "");

            Console.WriteLine("");
            Console.WriteLine("  Physical Nodes:");
            Console.WriteLine("");
            for (j = 0; j < NODE_NUM; j++)
            {
                Console.WriteLine("  "
                                  + t[0 + j * 2].ToString().PadLeft(10) + "  "
                                  + t[1 + j * 2].ToString().PadLeft(10) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("  The basis function values at basis nodes");
            Console.WriteLine("  should form the identity matrix.");
            Console.WriteLine("");

            basis_mn_t4(t, NODE_NUM, t, ref phi, ref dphidx, ref dphidy);

            for (j = 0; j < NODE_NUM; j++)
            {
                string cout = "";
                for (i = 0; i < NODE_NUM; i++)
                {
                    cout += "  " + phi[i + j * NODE_NUM].ToString().PadLeft(10);
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("  The X and Y derivatives should sum to 0.");
            Console.WriteLine("");
            Console.WriteLine("  dPhidX sum, dPhidY sum:");
            Console.WriteLine("");

            for (j = 0; j < NODE_NUM; j++)
            {
                sum_x = 0.0;
                sum_y = 0.0;
                for (i = 0; i < NODE_NUM; i++)
                {
                    sum_x = sum_x + dphidx[i + j * NODE_NUM];
                    sum_y = sum_y + dphidy[i + j * NODE_NUM];
                }

                Console.WriteLine("  "
                                  + sum_x.ToString().PadLeft(10) + "  "
                                  + sum_y.ToString().PadLeft(10) + "");
            }
        }

        public static void basis_mn_t6(double[] t, int n, double[] p, ref double[] phi,
                ref double[] dphidx, ref double[] dphidy)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BASIS_MN_T6: all bases at N points for a T6 element.
            //
            //  Discussion:
            //
            //    The routine is given the coordinates of the vertices and midside
            //    nodes of a triangle.  It works directly with these coordinates, and does 
            //    not refer to a reference element.
            //
            //    This routine requires that the midside nodes be "in line"
            //    with the vertices, that is, that the sides of the triangle be
            //    straight.  However, the midside nodes do not actually have to
            //    be halfway along the side of the triangle.  
            //
            //  The physical element T6:
            //
            //    This picture indicates the assumed ordering of the six nodes
            //    of the triangle.
            //
            //    |
            //    |   
            //    |        3
            //    |       . .
            //    |      .   .
            //    Y     6     5
            //    |    .       .
            //    |   .         .
            //    |  1-----4-----2
            //    |
            //    +--------X-------->
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 February 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double T[2*6], the nodal oordinates of the element.
            //    It is common to list these points in counter clockwise order.
            //
            //    Input, int N, the number of evaluation points.
            //
            //    Input, double P[2*N], the coordinates of the point where
            //    the basis functions are to be evaluated.
            //
            //    Output, double PHI[6*N], the value of the basis functions at P.
            //
            //    Output, double DPHIDX[6*N], DPHIDY[6*N], the value of the X 
            //    and Y derivatives of the basis functions at P.
            //
        {
            double gn;
            double gx;
            double hn;
            double hx;
            int j;

            for (j = 0; j < n; j++)
            {
                //
                //  Basis function 1: PHI(X,Y) = G(3,2) * H(6,4) / normalization.
                //
                gx = (p[0 + j * 2] - t[0 + 1 * 2]) * (t[1 + 2 * 2] - t[1 + 1 * 2])
                     - (t[0 + 2 * 2] - t[0 + 1 * 2]) * (p[1 + j * 2] - t[1 + 1 * 2]);

                gn = (t[0 + 0 * 2] - t[0 + 1 * 2]) * (t[1 + 2 * 2] - t[1 + 1 * 2])
                     - (t[0 + 2 * 2] - t[0 + 1 * 2]) * (t[1 + 0 * 2] - t[1 + 1 * 2]);

                hx = (p[0 + j * 2] - t[0 + 3 * 2]) * (t[1 + 5 * 2] - t[1 + 3 * 2])
                     - (t[0 + 5 * 2] - t[0 + 3 * 2]) * (p[1 + j * 2] - t[1 + 3 * 2]);

                hn = (t[0 + 0 * 2] - t[0 + 3 * 2]) * (t[1 + 5 * 2] - t[1 + 3 * 2])
                     - (t[0 + 5 * 2] - t[0 + 3 * 2]) * (t[1 + 0 * 2] - t[1 + 3 * 2]);

                phi[0 + j * 6] = (gx * hx) / (gn * hn);
                dphidx[0 + j * 6] = ((t[1 + 2 * 2] - t[1 + 1 * 2]) * hx
                                     + gx * (t[1 + 5 * 2] - t[1 + 3 * 2])) / (gn * hn);
                dphidy[0 + j * 6] = -((t[0 + 2 * 2] - t[0 + 1 * 2]) * hx
                                      + gx * (t[0 + 5 * 2] - t[0 + 3 * 2])) / (gn * hn);
                //
                //  Basis function 2: PHI(X,Y) = G(3,1) * H(4,5) / normalization.
                //
                gx = (p[0 + j * 2] - t[0 + 0 * 2]) * (t[1 + 2 * 2] - t[1 + 0 * 2])
                     - (t[0 + 2 * 2] - t[0 + 0 * 2]) * (p[1 + j * 2] - t[1 + 0 * 2]);

                gn = (t[0 + 1 * 2] - t[0 + 0 * 2]) * (t[1 + 2 * 2] - t[1 + 0 * 2])
                     - (t[0 + 2 * 2] - t[0 + 0 * 2]) * (t[1 + 1 * 2] - t[1 + 0 * 2]);

                hx = (p[0 + j * 2] - t[0 + 4 * 2]) * (t[1 + 3 * 2] - t[1 + 4 * 2])
                     - (t[0 + 3 * 2] - t[0 + 4 * 2]) * (p[1 + j * 2] - t[1 + 4 * 2]);

                hn = (t[0 + 1 * 2] - t[0 + 4 * 2]) * (t[1 + 3 * 2] - t[1 + 4 * 2])
                     - (t[0 + 3 * 2] - t[0 + 4 * 2]) * (t[1 + 1 * 2] - t[1 + 4 * 2]);

                phi[1 + j * 6] = (gx * hx) / (gn * hn);
                dphidx[1 + j * 6] = ((t[1 + 2 * 2] - t[1 + 0 * 2]) * hx
                                     + gx * (t[1 + 3 * 2] - t[1 + 4 * 2])) / (gn * hn);
                dphidy[1 + j * 6] = -((t[0 + 2 * 2] - t[0 + 0 * 2]) * hx
                                      + gx * (t[0 + 3 * 2] - t[0 + 4 * 2])) / (gn * hn);
                //
                //  Basis function 3: PHI(X,Y) = G(1,2) * H(5,6) / normalization.
                //
                gx = (p[0 + j * 2] - t[0 + 1 * 2]) * (t[1 + 0 * 2] - t[1 + 1 * 2])
                     - (t[0 + 0 * 2] - t[0 + 1 * 2]) * (p[1 + j * 2] - t[1 + 1 * 2]);

                gn = (t[0 + 2 * 2] - t[0 + 1 * 2]) * (t[1 + 0 * 2] - t[1 + 1 * 2])
                     - (t[0 + 0 * 2] - t[0 + 1 * 2]) * (t[1 + 2 * 2] - t[1 + 1 * 2]);

                hx = (p[0 + j * 2] - t[0 + 5 * 2]) * (t[1 + 4 * 2] - t[1 + 5 * 2])
                     - (t[0 + 4 * 2] - t[0 + 5 * 2]) * (p[1 + j * 2] - t[1 + 5 * 2]);

                hn = (t[0 + 2 * 2] - t[0 + 5 * 2]) * (t[1 + 4 * 2] - t[1 + 5 * 2])
                     - (t[0 + 4 * 2] - t[0 + 5 * 2]) * (t[1 + 2 * 2] - t[1 + 5 * 2]);

                phi[2 + j * 6] = (gx * hx) / (gn * hn);
                dphidx[2 + j * 6] = ((t[1 + 0 * 2] - t[1 + 1 * 2]) * hx
                                     + gx * (t[1 + 4 * 2] - t[1 + 5 * 2])) / (gn * hn);
                dphidy[2 + j * 6] = -((t[0 + 0 * 2] - t[0 + 1 * 2]) * hx
                                      + gx * (t[0 + 4 * 2] - t[0 + 5 * 2])) / (gn * hn);
                //
                //  Basis function 4: PHI(X,Y) = G(1,3) * H(2,3) / normalization.
                //
                gx = (p[0 + j * 2] - t[0 + 2 * 2]) * (t[1 + 0 * 2] - t[1 + 2 * 2])
                     - (t[0 + 0 * 2] - t[0 + 2 * 2]) * (p[1 + j * 2] - t[1 + 2 * 2]);

                gn = (t[0 + 3 * 2] - t[0 + 2 * 2]) * (t[1 + 0 * 2] - t[1 + 2 * 2])
                     - (t[0 + 0 * 2] - t[0 + 2 * 2]) * (t[1 + 3 * 2] - t[1 + 2 * 2]);

                hx = (p[0 + j * 2] - t[0 + 2 * 2]) * (t[1 + 1 * 2] - t[1 + 2 * 2])
                     - (t[0 + 1 * 2] - t[0 + 2 * 2]) * (p[1 + j * 2] - t[1 + 2 * 2]);

                hn = (t[0 + 3 * 2] - t[0 + 2 * 2]) * (t[1 + 1 * 2] - t[1 + 2 * 2])
                     - (t[0 + 1 * 2] - t[0 + 2 * 2]) * (t[1 + 3 * 2] - t[1 + 2 * 2]);

                phi[3 + j * 6] = (gx * hx) / (gn * hn);
                dphidx[3 + j * 6] = ((t[1 + 0 * 2] - t[1 + 2 * 2]) * hx
                                     + gx * (t[1 + 1 * 2] - t[1 + 2 * 2])) / (gn * hn);
                dphidy[3 + j * 6] = -((t[0 + 0 * 2] - t[0 + 2 * 2]) * hx
                                      + gx * (t[0 + 1 * 2] - t[0 + 2 * 2])) / (gn * hn);
                //
                //  Basis function 5: PHI(X,Y) = G(2,1) * H(3,1) / normalization.
                //
                gx = (p[0 + j * 2] - t[0 + 0 * 2]) * (t[1 + 1 * 2] - t[1 + 0 * 2])
                     - (t[0 + 1 * 2] - t[0 + 0 * 2]) * (p[1 + j * 2] - t[1 + 0 * 2]);

                gn = (t[0 + 4 * 2] - t[0 + 0 * 2]) * (t[1 + 1 * 2] - t[1 + 0 * 2])
                     - (t[0 + 1 * 2] - t[0 + 0 * 2]) * (t[1 + 4 * 2] - t[1 + 0 * 2]);

                hx = (p[0 + j * 2] - t[0 + 0 * 2]) * (t[1 + 2 * 2] - t[1 + 0 * 2])
                     - (t[0 + 2 * 2] - t[0 + 0 * 2]) * (p[1 + j * 2] - t[1 + 0 * 2]);

                hn = (t[0 + 4 * 2] - t[0 + 0 * 2]) * (t[1 + 2 * 2] - t[1 + 0 * 2])
                     - (t[0 + 2 * 2] - t[0 + 0 * 2]) * (t[1 + 4 * 2] - t[1 + 0 * 2]);

                phi[4 + j * 6] = (gx * hx) / (gn * hn);
                dphidx[4 + j * 6] = ((t[1 + 1 * 2] - t[1 + 0 * 2]) * hx
                                     + gx * (t[1 + 2 * 2] - t[1 + 0 * 2])) / (gn * hn);
                dphidy[4 + j * 6] = -((t[0 + 1 * 2] - t[0 + 0 * 2]) * hx
                                      + gx * (t[0 + 2 * 2] - t[0 + 0 * 2])) / (gn * hn);
                //
                //  Basis function 6: PHI(X,Y) = G(1,2) * H(3,2) / normalization.
                //
                gx = (p[0 + j * 2] - t[0 + 1 * 2]) * (t[1 + 0 * 2] - t[1 + 1 * 2])
                     - (t[0 + 0 * 2] - t[0 + 1 * 2]) * (p[1 + j * 2] - t[1 + 1 * 2]);

                gn = (t[0 + 5 * 2] - t[0 + 1 * 2]) * (t[1 + 0 * 2] - t[1 + 1 * 2])
                     - (t[0 + 0 * 2] - t[0 + 1 * 2]) * (t[1 + 5 * 2] - t[1 + 1 * 2]);

                hx = (p[0 + j * 2] - t[0 + 1 * 2]) * (t[1 + 2 * 2] - t[1 + 1 * 2])
                     - (t[0 + 2 * 2] - t[0 + 1 * 2]) * (p[1 + j * 2] - t[1 + 1 * 2]);

                hn = (t[0 + 5 * 2] - t[0 + 1 * 2]) * (t[1 + 2 * 2] - t[1 + 1 * 2])
                     - (t[0 + 2 * 2] - t[0 + 1 * 2]) * (t[1 + 5 * 2] - t[1 + 1 * 2]);

                phi[5 + j * 6] = (gx * hx) / (gn * hn);
                dphidx[5 + j * 6] = ((t[1 + 0 * 2] - t[1 + 1 * 2]) * hx
                                     + gx * (t[1 + 2 * 2] - t[1 + 1 * 2])) / (gn * hn);
                dphidy[5 + j * 6] = -((t[0 + 0 * 2] - t[0 + 1 * 2]) * hx
                                      + gx * (t[0 + 2 * 2] - t[0 + 1 * 2])) / (gn * hn);
            }

            return;
        }

        public static void basis_mn_t6_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BASIS_MN_T6_TEST verifies BASIS_MN_T6.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 February 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    None
            //
        {
            int NODE_NUM = 6;

            double[] dphidx = new double[NODE_NUM * NODE_NUM];
            double[] dphidy = new double[NODE_NUM * NODE_NUM];
            int i;
            int j;
            double[] phi = new double[NODE_NUM * NODE_NUM];
            double sum_x;
            double sum_y;
            double[] t =
            {
                2.0, 0.0,
                4.0, 3.0,
                0.0, 4.0,
                3.0, 1.5,
                2.0, 3.5,
                1.0, 2.0
            };

            Console.WriteLine("");
            Console.WriteLine("BASIS_MN_T6_TEST:");
            Console.WriteLine("  Verify basis functions for element T6.");
            Console.WriteLine("");
            Console.WriteLine("  Number of nodes = " + NODE_NUM + "");

            Console.WriteLine("");
            Console.WriteLine("  Physical Nodes:");
            Console.WriteLine("");
            for (j = 0; j < NODE_NUM; j++)
            {
                Console.WriteLine("  "
                                  + j.ToString().PadLeft(6) + "  "
                                  + t[0 + j * 2].ToString().PadLeft(7) + "  "
                                  + t[1 + j * 2].ToString().PadLeft(7) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("  The basis function values at basis nodes");
            Console.WriteLine("  should form the identity matrix.");
            Console.WriteLine("");

            basis_mn_t6(t, NODE_NUM, t, ref phi, ref dphidx, ref dphidy);

            for (i = 0; i < NODE_NUM; i++)
            {
                string cout = "";
                for (j = 0; j < NODE_NUM; j++)
                {
                    cout += "  " + phi[i + j * NODE_NUM].ToString().PadLeft(7);
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("  The X and Y derivatives should sum to 0.");
            Console.WriteLine("");
            Console.WriteLine("  dPhidX sum, dPhidY sum:");
            Console.WriteLine("");
            for (j = 0; j < NODE_NUM; j++)
            {
                sum_x = 0.0;
                for (i = 0; i < NODE_NUM; i++)
                {
                    sum_x = sum_x + dphidx[i + j * NODE_NUM];
                }

                sum_y = 0.0;
                for (i = 0; i < NODE_NUM; i++)
                {
                    sum_y = sum_y + dphidy[i + j * NODE_NUM];
                }

                Console.WriteLine("  "
                                  + sum_x.ToString().PadLeft(14) + "  "
                                  + sum_y.ToString().PadLeft(14) + "");
            }

        }

    }
}