namespace Burkardt.FEM;

public class Derivative
{
    public static void derivative_average_t3(int node_num, double[] node_xy, int element_num,
            int[] element_node, double[] c, double[] dcdx, double[] dcdy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DERIVATIVE_AVERAGE_T3 averages derivatives at the nodes of a T3 mesh.
        //
        //  Discussion:
        //
        //    This routine can be used to compute an averaged nodal value of any
        //    quantity associated with the finite element function.  At a node
        //    that is shared by several elements, the fundamental function
        //    U will be continuous, but its spatial derivatives, for instance,
        //    will generally be discontinuous.  This routine computes the
        //    value of the spatial derivatives in each element, and averages
        //    them, to make a reasonable assignment of a nodal value.
        //
        //    Note that the ELEMENT_NODE array is assumed to be 1-based, rather
        //    than 0-based.  Thus, entries from this array must be decreased by
        //    1 before being used!
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 June 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        //
        //    Input, int ELEMENT_NODE[3*ELEMENT_NUM], the element->node data.
        //
        //    Input, double C[NODE_NUM], the finite element coefficient vector.
        //
        //    Output, double DCDX[NODE_NUM], DCDY[NODE_NUM], the averaged
        //    values of dCdX and dCdY at the nodes.
        //
    {
        int OFFSET = 1;

        int dim;
        double[] dphidx = new double[3 * 3];
        double[] dphidy = new double[3 * 3];
        int element;
        int j;
        int node;
        int[] node_count = new int[node_num];
        int node_global1;
        int node_global2;
        int node_local1;
        int node_local2;
        double[] phi = new double[3 * 3];
        double[] t = new double[2 * 3];

        for (node = 0; node < node_num; node++)
        {
            node_count[node] = 0;
            dcdx[node] = 0.0;
            dcdy[node] = 0.0;
        }

        //
        //  Consider every element.
        //
        for (element = 0; element < element_num; element++)
        {
            //
            //  Get the coordinates of the nodes of the element.
            //
            for (j = 0; j < 3; j++)
            {
                for (dim = 0; dim < 2; dim++)
                {
                    t[dim + 2 * j] = node_xy[dim + (element_node[j + element * 3] - OFFSET)];
                }
            }

            //
            //  Evaluate the X and Y derivatives of the 3 basis functions at the
            //  3 nodes.
            //
            Basis_mn.basis_mn_t3(t, 3, t, ref phi, ref dphidx, ref dphidy);
            //
            //  Evaluate dCdX and dCdY at each node in the element, and add
            //  them to the running totals.
            //
            for (node_local1 = 0; node_local1 < 3; node_local1++)
            {
                node_global1 = element_node[node_local1 + element * 3] - OFFSET;

                for (node_local2 = 0; node_local2 < 3; node_local2++)
                {
                    node_global2 = element_node[node_local2 + element * 3] - OFFSET;

                    dcdx[node_global1] += c[node_global2] * dphidx[node_local2 + node_local1 * 3];

                    dcdy[node_global1] += c[node_global2] * dphidy[node_local2 + node_local1 * 3];
                }

                node_count[node_global1] += 1;
            }
        }

        //
        //  Average the running totals.
        //
        for (node = 0; node < node_num; node++)
        {
            dcdx[node] /= node_count[node];
            dcdy[node] /= node_count[node];
        }
    }
}