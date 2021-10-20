using System;

namespace Burkardt.FEM
{
    public static class Boundary
    {
        public static void boundary ( int nx, int ny, int node_num, double[] node_xy, int[] indx,
        int ib, int nunk, ref double[] a, ref double[] f, Func<double, double, ExactResult> exact )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BOUNDARY modifies the linear system for boundary conditions.
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
        //    Input, int NX, NY, controls the number of elements along the
        //    X and Y directions.  The number of elements will be
        //    2 * ( NX - 1 ) * ( NY - 1 ).
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the coordinates of nodes.
        //
        //    Input, int INDX[NODE_NUM], gives the index of the unknown quantity
        //    associated with the given node.
        //
        //    Input, int IB, the half-bandwidth of the matrix.
        //
        //    Input, int NUNK, the number of unknowns.
        //
        //    Input/output, double A[(3*IB+1)*NUNK], the NUNK by NUNK
        //    coefficient matrix, stored in a compressed format.
        //    On output, A has been adjusted for boundary conditions.
        //
        //    Input/output, double F[NUNK], the right hand side.
        //    On output, F has been adjusted for boundary conditions.
        //
        {
            int col;
            double dudx;
            double dudy;
            int i;
            int j;
            int jhi;
            int jlo;
            int node;
            int row;
            double u;
            double x;
            double y;
            //
            //  Consider each node.
            //
            node = 0;

            for (row = 1; row <= 2 * ny - 1; row++)
            {
                for (col = 1; col <= 2 * nx - 1; col++)
                {
                    node = node + 1;

                    if (row == 1 ||
                        row == 2 * ny - 1 ||
                        col == 1 ||
                        col == 2 * nx - 1)
                    {
                        i = indx[node - 1];
                        x = node_xy[0 + (node - 1) * 2];
                        y = node_xy[1 + (node - 1) * 2];
                        ExactResult res = exact(x, y);
                        u = res.u;
                        dudx = res.dudx;
                        dudy = res.dudy;

                        jlo = Math.Max(i - ib, 1);
                        jhi = Math.Min(i + ib, nunk);

                        for (j = jlo; j <= jhi; j++)
                        {
                            a[i - j + 2 * ib + (j - 1) * (3 * ib + 1)] = 0.0;
                        }

                        a[i - i + 2 * ib + (i - 1) * (3 * ib + 1)] = 1.0;

                        f[i - 1] = u;
                    }
                }
            }
        }

        public static void indx_set(int nx, int ny, int node_num, ref int[] indx, ref int nunk)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    INDX_SET assigns a boundary value index or unknown value index at each node.
            //
            //  Discussion:
            //
            //    Every finite element node will is assigned an index which
            //    indicates the finite element basis function and its coefficient
            //    which are associated with that node.
            //
            //  Example:
            //
            //    On a simple 5 by 5 grid, where the nodes are numbered starting
            //    at the lower left, and increasing in X first, we would have the
            //    following values of INDX:
            //
            //       21  22  23  24  25
            //       16  17  18  19  20
            //       11  12  13  14  15
            //        6   7   8   9  10
            //        1   2   3   4   5
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
            //    Input, int NX, NY, the number of elements in the X and Y directions.
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Output, int INDX[NODE_NUM], the index of the unknown in the finite
            //    element linear system.
            //
            //    Output, int *NUNK, the number of unknowns.
            //
        {
            int i;
            int in_;
            int j;

            nunk = 0;
            in_ = 0;

            for (j = 1; j <= 2 * ny - 1; j++)
            {
                for (i = 1; i <= 2 * nx - 1; i++)
                {
                    in_ = in_ + 1;
                    nunk = nunk + 1;
                    indx[in_ - 1] = nunk;
                }
            }
        }
    }
}