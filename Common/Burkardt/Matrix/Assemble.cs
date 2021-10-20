using System;
using Burkardt.Function;

namespace Burkardt.MatrixNS
{
    public static partial class Matrix
    {
        public static void assemble(int node_num, double[] node_xy, int nnodes,
                int element_num, int[] element_node, int nq,
                double[] wq, double[] xq, double[] yq, double[] element_area, int[] indx,
                int ib, int nunk, ref double[] a, ref double[] f)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ASSEMBLE assembles the matrix and right-hand side using piecewise quadratics.
            //
            //  Discussion:
            //
            //    The matrix is known to be banded.  A special matrix storage format
            //    is used to reduce the space required.  Details of this format are
            //    discussed in the routine DGB_FA.
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
            //    Output, double A[(3*IB+1)*NUNK], the NUNK by NUNK coefficient matrix,
            //    stored in a compressed format.
            //
            //    Output, double F[NUNK], the right hand side.
            //
            //    Input, int IB, the half-bandwidth of the matrix.
            //
            //    Input, double NODE_XY[2*NODE_NUM], the X and Y coordinates of nodes.
            //
            //    Input, double XQ[NQ*ELEMENT_NUM], YQ[NQ*ELEMENT_NUM], the X and Y
            //    coordinates of the quadrature points in each element.
            //
            //    Input, double WQ[NQ], quadrature weights.
            //
            //    Input, double ELEMENT_AREA[ELEMENT_NUM], the area of each element.
            //
            //    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM]; ELEMENT_NODE(I,J) is the global
            //    index of local node I in element J.
            //
            //    Input, int INDX[NODE_NUM], gives the index of the unknown quantity
            //    associated with the given node.
            //
            //    Input, int NNODES, the number of nodes used to form one element.
            //
            //    Input, int NUNK, the number of unknowns.
            //
            //    Input, int NQ, the number of quadrature points used in assembly.
            //
            //    Input, int ELEMENT_NUM, the number of elements.
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //  Local parameters:
            //
            //    Local, double BB, BX, BY, the value of some basis function
            //    and its first derivatives at a quadrature point.
            //
            //    Local, double BJ, DBJDX, DBJDY, the value of another basis
            //    function and its first derivatives at a quadrature point.
            //
            //    Local, int NODE_NUM, the number of global nodes.
            //
        {
            double aij;
            int basis;
            double bi = 0;
            double bj = 0;
            double dbidx = 0;
            double dbidy = 0;
            double dbjdx = 0;
            double dbjdy = 0;
            int element;
            int i;
            int ip;
            int ipp;
            int j;
            int quad;
            int test;
            double w;
            double x;
            double y;
            //
            //  Initialize the arrays to zero.
            //
            for (i = 1; i <= nunk; i++)
            {
                f[i - 1] = 0.0E+00;
            }

            for (j = 1; j <= nunk; j++)
            {
                for (i = 1; i <= 3 * ib + 1; i++)
                {
                    a[i - 1 + (j - 1) * (3 * ib + 1)] = 0.0E+00;
                }
            }

            //
            //  The actual values of A and F are determined by summing up
            //  contributions from all the elements.
            //
            for (element = 1; element <= element_num; element++)
            {
                for (quad = 1; quad <= nq; quad++)
                {
                    x = xq[quad - 1 + (element - 1) * nq];
                    y = yq[quad - 1 + (element - 1) * nq];
                    w = element_area[element - 1] * wq[quad - 1];

                    for (test = 1; test <= nnodes; test++)
                    {
                        ip = element_node[test - 1 + (element - 1) * nnodes];
                        i = indx[ip - 1];

                        Quadratic.qbf(x, y, element, test, node_xy, element_node,
                            element_num, nnodes, node_num, ref bi, ref dbidx, ref dbidy);

                        f[i - 1] = f[i - 1] + w * rhs(x, y) * bi;
                        //
                        //  We are about to compute a contribution associated with the
                        //  I-th test function and the J-th basis function, and add this
                        //  to the entry A(I,J).
                        //
                        //  Because of the compressed storage of the matrix, the element
                        //  will actually be stored in A(I-J+2*IB+1,J).
                        //
                        //  Two extra complications: we are storing the array as a vector,
                        //  and C uses 0-based indices rather than 1-based indices.
                        //
                        //  Therefore, we ACTUALLY store the entry in A[I-J+2*IB+1-1 + (J-1) * (3*IB+1)];
                        //
                        for (basis = 1; basis <= nnodes; basis++)
                        {
                            ipp = element_node[basis - 1 + (element - 1) * nnodes];
                            j = indx[ipp - 1];

                            Quadratic.qbf(x, y, element, basis, node_xy, element_node,
                                element_num, nnodes, node_num, ref bj, ref dbjdx, ref dbjdy);

                            aij = dbidx * dbjdx + dbidy * dbjdy;

                            a[i - j + 2 * ib + (j - 1) * (3 * ib + 1)] =
                                a[i - j + 2 * ib + (j - 1) * (3 * ib + 1)] + w * aij;
                        }
                    }
                }
            }

        }
        
        public static double rhs ( double x, double y )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    RHS gives the right-hand side of the differential equation.
            //
            //  Discussion:
            //
            //    The function specified here depends on the problem being
            //    solved.  This is one of the routines that a user will
            //    normally want to change.
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
            //    Input, double X, Y, the coordinates of a point
            //    in the region, at which the right hand side of the
            //    differential equation is to be evaluated.
            //
            //    Output, double RHS, the value of the right
            //    hand side of the differential equation at (X,Y).
            //
        {

            double value;

            value = 2.0E+00 * Math.PI * Math.PI * Math.Sin ( Math.PI * x ) * Math.Sin ( Math.PI * y );

            return value;
        }
    }
}