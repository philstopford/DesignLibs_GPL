using System;

namespace Burkardt.FEM
{
    public static class Bandwidth
    {
        public static void bandwidth_mesh(int element_order, int element_num, int[] element_node,
                ref int ml, ref int mu, ref int m)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BANDWIDTH_MESH determines the bandwidth of the coefficient matrix.
            //
            //  Discussion:
            //
            //    The quantity computed here is the "geometric" bandwidth determined
            //    by the finite element mesh alone.
            //
            //    If a single finite element variable is associated with each node
            //    of the mesh, and if the nodes and variables are numbered in the
            //    same way, then the geometric bandwidth is the same as the bandwidth
            //    of a typical finite element matrix.
            //
            //    The bandwidth M is defined in terms of the lower and upper bandwidths:
            //
            //      M = ML + 1 + MU
            //
            //    where 
            //
            //      ML = maximum distance from any diagonal entry to a nonzero
            //      entry in the same row, but earlier column,
            //
            //      MU = maximum distance from any diagonal entry to a nonzero
            //      entry in the same row, but later column.
            //
            //    Because the finite element node adjacency relationship is symmetric,
            //    we are guaranteed that ML = MU.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 January 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int ELEMENT_ORDER, the order of the elements.
            //
            //    Input, int ELEMENT_NUM, the number of elements.
            //
            //    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM];
            //    ELEMENT_NODE(I,J) is the global index of local node I in element J.
            //
            //    Output, int *ML, *MU, the lower and upper bandwidths of the matrix.
            //
            //    Output, int *M, the bandwidth of the matrix.
            //
        {
            int element;
            int global_i;
            int global_j;
            int local_i;
            int local_j;

            ml = 0;
            mu = 0;

            for (element = 0; element < element_num; element++)
            {
                for (local_i = 0; local_i < element_order; local_i++)
                {
                    global_i = element_node[local_i + element * element_order];

                    for (local_j = 0; local_j < element_order; local_j++)
                    {
                        global_j = element_node[local_j + element * element_order];

                        mu = Math.Max(mu, global_j - global_i);
                        ml = Math.Max(ml, global_i - global_j);
                    }
                }
            }

            m = ml + 1 + mu;
        }

        public static void bandwidth_var(int element_order, int element_num, int[] element_node,
        int node_num, int[] var_node, int var_num, int[] var, ref int ml, ref int mu,
        ref int m )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BANDWIDTH_VAR determines the bandwidth for finite element variables.
        //
        //  Discussion:
        //
        //    We assume that, attached to each node in the finite element mesh
        //    there are a (possibly zero) number of finite element variables.
        //    We wish to determine the bandwidth necessary to store the stiffness
        //    matrix associated with these variables.
        //
        //    An entry K(I,J) of the stiffness matrix must be zero unless the
        //    variables I and J correspond to nodes N(I) and N(J) which are
        //    common to some element.
        //
        //    In order to determine the bandwidth of the stiffness matrix, we
        //    essentially seek a nonzero entry K(I,J) for which abs ( I - J )
        //    is maximized.
        //
        //    The bandwidth M is defined in terms of the lower and upper bandwidths:
        //
        //      M = ML + 1 + MU
        //
        //    where
        //
        //      ML = maximum distance from any diagonal entry to a nonzero
        //      entry in the same row, but earlier column,
        //
        //      MU = maximum distance from any diagonal entry to a nonzero
        //      entry in the same row, but later column.
        //
        //    We assume the finite element variable adjacency relationship is 
        //    symmetric, so we are guaranteed that ML = MU.
        //
        //    Note that the user is free to number the variables in any way
        //    whatsoever, and to associate variables to nodes in any way,
        //    so that some nodes have no variables, some have one, and some
        //    have several.  
        //
        //    The storage of the indices of the variables is fairly simple.
        //    In VAR, simply list all the variables associated with node 1, 
        //    then all those associated with node 2, and so on.  Then set up
        //    the pointer array VAR_NODE so that we can jump to the section of
        //    VAR where the list begins for any particular node.
        //
        //    The routine does not check that each variable is only associated
        //    with a single node.  This would normally be the case in a finite
        //    element setting.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ELEMENT_ORDER, the order of the elements.
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        //
        //    Input,  ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM];
        //    ELEMENT_NODE(I,J) is the global index of local node I in element J.
        //
        //    Output, int *ML, *MU, the lower and upper bandwidths of the matrix.
        //
        //    Output, int *M, the bandwidth of the matrix.
        //
        {
            int element;
            int node_global_i;
            int node_global_j;
            int node_local_i;
            int node_local_j;
            int var_global_i;
            int var_global_j;
            int var_local_i;
            int var_local_j;

            ml = 0;
            mu = 0;

            for (element = 0; element < element_num; element++)
            {
                for (node_local_i = 0; node_local_i < element_order; node_local_i++)
                {
                    node_global_i = element_node[node_local_i + element * element_order];

                    for (var_local_i = var_node[node_global_i - 1];
                        var_local_i <= var_node[node_global_i] - 1;
                        var_local_i++)
                    {
                        var_global_i = var[var_local_i - 1];

                        for (node_local_j = 0; node_local_j < element_order; node_local_j++)
                        {
                            node_global_j = element_node[node_local_j + element * element_order];

                            for (var_local_j = var_node[node_global_j - 1];
                                var_local_j <= var_node[node_global_j] - 1;
                                var_local_j++)
                            {
                                var_global_j = var[var_local_j - 1];

                                mu = Math.Max(mu, var_global_j - var_global_i);
                                ml = Math.Max(ml, var_global_i - var_global_j);
                            }
                        }
                    }
                }
            }

            m = ml + 1 + mu;
        }
    }
}