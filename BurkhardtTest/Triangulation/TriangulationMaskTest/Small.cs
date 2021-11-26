﻿namespace TriangulationMaskTest;

public static class Small
{
    public static bool triangle_mask(int dim_num, int triangle_order, int[] nodes,
            double[] coord)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_MASK is a user routine which masks triangles.
        //
        //  Discussion:
        //
        //    The region to be considered is the [0,4]x[0,4] square.
        //
        //    We want to remove the lower left triangular corner,
        //    and part of the upper right triangular corner.
        //
        //    The following diagram of the 25 nodes indicates by "O" the
        //    nodes that should end up being deleted, although the deletion
        //    is actually done by triangles.
        //
        //    Before masking:
        //
        //      X - X - X - X - X
        //      | \ | \ | \ | \ |
        //      X - X - X - X - X
        //      | \ | \ | \ | \ |
        //      X - X - X - X - X
        //      | \ | \ | \ | \ |
        //      X - X - X - X - X
        //      | \ | \ | \ | \ |
        //      X - X - X - X - X
        //
        //    After masking:
        //
        //      X - X   O   O   O
        //      : \ : \.          
        //      X - X - X   O   O
        //      : \ : \ : \.      
        //      X - X - X - X - X
        //        \ : \ : \ : \ :
        //      O   X - X - X - X
        //            \ : \ : \ :
        //      O   O   X - X - X
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int TRIANGLE_ORDER, the number of nodes in the triangle.
        //
        //    Input, int NODES[TRIANGLE_ORDER], the indices of the nodes.
        //
        //    Input, double COORD[DIM_NUM*TRIANGLE_ORDER], the coordinates
        //    of the nodes.
        //
        //    Output, bool TRIANGLE_MASK, is TRUE if the triangle should be discarded,
        //    and FALSE if the triangle should be retained.
        //
    {
        double[] centroid = new double[2];
        int dim;
        bool mask;
        //
        //  Compute the centroid.
        //
        for (dim = 0; dim < 2; dim++)
        {
            centroid[dim] = 0.0;
            int order;
            for (order = 0; order < triangle_order; order++)
            {
                centroid[dim] += coord[dim + order * dim_num];
            }

            centroid[dim] /= triangle_order;
        }

        switch (centroid[0] + centroid[1])
        {
            //
            //  Remove the lower left corner
            //
            case < 2.0:
            //
            //  Remove the upper right section.
            //
            case > 5.0 when 2.0 < centroid[1]:
                mask = true;
                break;
            //
            default:
                mask = false;
                break;
        }

        return mask;
    }
}