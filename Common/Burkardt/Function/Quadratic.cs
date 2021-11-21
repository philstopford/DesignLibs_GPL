using System;

namespace Burkardt.Function;

public static class Quadratic
{
    public static void qbf(double x, double y, int element, int inode, double[] node_xy,
            int[] element_node, int element_num, int nnodes,
            int node_num, ref double b, ref double dbdx, ref double dbdy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    QBF evaluates the quadratic basis functions.
        //
        //  Discussion:
        //
        //    This routine assumes that the "midpoint" nodes are, in fact,
        //    exactly the average of the two extreme nodes.  This is NOT true
        //    for a general quadratic triangular element.
        //
        //    Assuming this property of the midpoint nodes makes it easy to
        //    determine the values of (R,S) in the reference element that
        //    correspond to (X,Y) in the physical element.
        //
        //    Once we know the (R,S) coordinates, it's easy to evaluate the
        //    basis functions and derivatives.
        //
        //  The physical element T6:
        //
        //    In this picture, we don't mean to suggest that the bottom of
        //    the physical triangle is horizontal.  However, we do assume that
        //    each of the sides is a straight line, and that the intermediate
        //    points are exactly halfway on each side.
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
        //  Reference element T6:
        //
        //    In this picture of the reference element, we really do assume
        //    that one side is vertical, one horizontal, of length 1.
        //
        //    |
        //    |
        //    1  3
        //    |  ..
        //    |  . .
        //    S  6  5
        //    |  .   .
        //    |  .    .
        //    0  1--4--2
        //    |
        //    +--0--R--1-------->
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 September 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, Y, the (global) coordinates of the point
        //    at which the basis function is to be evaluated.
        //
        //    Input, int ELEMENT, the index of the element which contains the point.
        //
        //    Input, int INODE, the local index (between 1 and 6) that
        //    specifies which basis function is to be evaluated.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the nodes.
        //
        //    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM];
        //    ELEMENT_NODE(I,J) is the global index of local node I in element J.
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        //
        //    Input, int NNODES, the number of nodes used to form one element.
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Output, double *B, *DBDX, *DBDY, the value of the basis function
        //    and its X and Y derivatives at (X,Y).
        //
    {
        double dbdr;
        double dbds;
        int i;
        double[] xn = new double[6];
        double[] yn=  new double[6];

        for (i = 0; i < 6; i++)
        {
            xn[i] = node_xy[0 + (element_node[i + (element - 1) * nnodes] - 1) * 2];
            yn[i] = node_xy[1 + (element_node[i + (element - 1) * nnodes] - 1) * 2];
        }

        //
        //  Determine the (R,S) coordinates corresponding to (X,Y).
        //
        //  What is happening here is that we are solving the linear system:
        //
        //    ( X2-X1  X3-X1 ) * ( R ) = ( X - X1 )
        //    ( Y2-Y1  Y3-Y1 )   ( S )   ( Y - Y1 )
        //
        //  by computing the inverse of the coefficient matrix and multiplying
        //  it by the right hand side to get R and S.
        //
        //  The values of dRdX, dRdY, dSdX and dSdY are easily from the formulas
        //  for R and S.
        //
        double det = (xn[1] - xn[0]) * (yn[2] - yn[0])
                     - (xn[2] - xn[0]) * (yn[1] - yn[0]);

        double r = ((yn[2] - yn[0]) * (x - xn[0])
                    + (xn[0] - xn[2]) * (y - yn[0])) / det;

        double drdx = (yn[2] - yn[0]) / det;
        double drdy = (xn[0] - xn[2]) / det;

        double s = ((yn[0] - yn[1]) * (x - xn[0])
                    + (xn[1] - xn[0]) * (y - yn[0])) / det;

        double dsdx = (yn[0] - yn[1]) / det;
        double dsdy = (xn[1] - xn[0]) / det;
        switch (inode)
        {
            //
            //  The basis functions can now be evaluated in terms of the
            //  reference coordinates R and S.  It's also easy to determine
            //  the values of the derivatives with respect to R and S.
            //
            case 1:
                b = 2.0E+00 * (1.0E+00 - r - s) * (0.5E+00 - r - s);
                dbdr = -3.0E+00 + 4.0E+00 * r + 4.0E+00 * s;
                dbds = -3.0E+00 + 4.0E+00 * r + 4.0E+00 * s;
                break;
            case 2:
                b = 2.0E+00 * r * (r - 0.5E+00);
                dbdr = -1.0E+00 + 4.0E+00 * r;
                dbds = 0.0E+00;
                break;
            case 3:
                b = 2.0E+00 * s * (s - 0.5E+00);
                dbdr = 0.0E+00;
                dbds = -1.0E+00 + 4.0E+00 * s;
                break;
            case 4:
                b = 4.0E+00 * r * (1.0E+00 - r - s);
                dbdr = 4.0E+00 - 8.0E+00 * r - 4.0E+00 * s;
                dbds = -4.0E+00 * r;
                break;
            case 5:
                b = 4.0E+00 * r * s;
                dbdr = 4.0E+00 * s;
                dbds = 4.0E+00 * r;
                break;
            case 6:
                b = 4.0E+00 * s * (1.0E+00 - r - s);
                dbdr = -4.0E+00 * s;
                dbds = 4.0E+00 - 4.0E+00 * r - 8.0E+00 * s;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("QBF - Fatal error!");
                Console.WriteLine("  Request for local basis function INODE = " + inode + "");
                return;
        }

        //
        //  We need to convert the derivative information from (R(X,Y),S(X,Y))
        //  to (X,Y) using the chain rule.
        //
        dbdx = dbdr * drdx + dbds * dsdx;
        dbdy = dbdr * drdy + dbds * dsdy;
    }
}