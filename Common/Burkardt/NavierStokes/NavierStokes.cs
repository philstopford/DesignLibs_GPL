using Burkardt.Types;

namespace Burkardt.NavierStokesNS
{
    public static class NavierStokes
    {
        public static int ns_t6_var_count(int element_num, int[] element_node, int node_num,
                ref int[] var_node)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NS_T6_VAR_COUNT counts the Navier Stokes variables on a T6 grid.
            //
            //  Discussion:
            //
            //    We are given a mesh of T6 elements, and asked to count, in advance,
            //    the number of Navier-Stokes variables associated with the grid.
            //    In particular, every node has two velocity variables associated with
            //    it, but only a node that is a vertex of the element will also have
            //    an associated pressure variable.
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
            //    Input, int ELEMENT_NUM, the number of elements.
            //
            //    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM]; 
            //    ELEMENT_NODE(I,J) is the global index of local node I in element J.
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Output, int VAR_NODE[NODE_NUM+1], used to find the variables 
            //    associated with a given node, which are in VAR in locations 
            //    VAR_NODE(NODE) to VAR_NODE(NODE+1)-1.  Note that the last entry of
            //    this array points to the location just after the last location in VAR.
            //
            //    Output, int NS_T6_VAR_COUNT, the number of variables.
            //
        {
            int count;
            int element;
            int element_order = 6;
            int node;
            int num;
            int order;
            int var_num;
            //
            //  Our job is easy once we determine which nodes are vertices.
            //  So to begin with, let VAR_NODE count the number of variables
            //  associated with each node.
            //
            for (node = 0; node < node_num; node++)
            {
                var_node[node] = 2;
            }

            for (element = 0; element < element_num; element++)
            {
                for (order = 0; order < 3; order++)
                {
                    node = element_node[order + element * element_order];
                    var_node[node - 1] = 3;
                }
            }

            //
            //  Count them.
            //
            var_num = 0;
            for (node = 0; node < node_num; node++)
            {
                var_num = var_num + var_node[node];
            }

            //
            //  Make pointers.
            //
            count = 1;

            for (node = 0; node < node_num; node++)
            {
                num = var_node[node];
                var_node[node] = count;
                count = count + num;
            }

            var_node[node_num] = count;

            return var_num;
        }

        public static int[] ns_t6_var_set(int element_num, int[] element_node, int node_num,
                int[] var_node, int var_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NS_T6_VAR_SET sets the Navier Stokes variables on a T6 grid.
            //
            //  Discussion:
            //
            //    We are given a mesh of T6 elements, and asked to create the natural
            //    list of indices for Navier-Stokes variables associated with each node.
            //    In particular, every node has two velocity variables associated with
            //    it, but only a node that is a vertex of the element will also have
            //    an associated pressure variable.
            //
            //    The hard work has been done for us alread, because the variables
            //    have been counted, and the pointers to the occurrence of the
            //    first variable associated with each node have been created.
            //
            //    The indexing of the nodes can be arbitrary, although a bad
            //    indexing will result in a miserably large bandwidth (if band
            //    storage is being tried for the stiffness matrix).  Here, we
            //    simply try to natural ordering, that is, the variables are
            //    numbered in order of the node with which they are associated.
            //
            //    For the Navier Stokes problem on a T6 grid, we take it as
            //    understood that each node has either 2 or 3 variables associated
            //    with it, that the first two are always the horizontal and
            //    then vertical velocity coefficients, and that the third, if
            //    present, is a pressure coefficient.
            //
            //    In other settings, it might be necessary not merely to assign
            //    the variables an index, but also to identify them as to type.
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
            //    Input, int ELEMENT_NUM, the number of elements.
            //
            //    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM]; 
            //    ELEMENT_NODE(I,J) is the global index of local node I in element J.
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Input, int VAR_NODE[NODE_NUM+1], used to find the variables 
            //    associated with a given node, which are in VAR in locations 
            //    VAR_NODE(NODE) to VAR_NODE(NODE+1)-1.  Note that the last entry of
            //    this array points to the location just after the last location in VAR.
            //
            //    Input, int VAR_NUM, the number of variables.
            //
            //    Output, int NS_T6_VAR_SET[VAR_NUM], the indexes of the variables, which
            //    are simply 1, 2, 3, ..., VAR_NUM.
            //
        {
            int i;
            int[] var;

            var = new int[var_num];

            for (i = 0; i < var_num; i++)
            {
                var[i] = i + 1;
            }

            return var;
        }

        public static int ns_adj_col_set(int node_num, int triangle_num, int variable_num,
            int[] triangle_node, int[] triangle_neighbor, int[] node_u_variable,
        int[] node_v_variable, int[] node_p_variable, ref int[] adj_col )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NS_ADJ_COL_SET sets the COL array in a Navier Stokes triangulation.
        //
        //  Discussion:
        //
        //    This routine also counts the the value and returns the value of
        //    ADJ_NUM, the number of Navier-Stokes variable adjacencies, which
        //    should be identical to the value that would have been computed
        //    by calling NS_ADJ_COUNT.
        //
        //    This routine is called to set up the ADJ_COL array, which indicates
        //    the number of entries needed to store each column in the sparse
        //    compressed column storage to be used for the adjancency matrix.
        //
        //    The triangulation is assumed to involve 6-node triangles.
        //
        //    Variables for the horizontal and vertical velocities are associated
        //    with every node.  Variables for the pressure are associated only with
        //    the vertex nodes.
        //
        //    We are interested in determining the number of nonzero entries in the
        //    stiffness matrix of the Stokes equations, or the jacobian matrix of
        //    the Navier Stokes equations.  To this end, we will say, somewhat
        //    too broadly, that two variables are "adjacent" if their associated
        //    nodes both occur in some common element.  This adjacency of variables
        //    I and J is taken to be equivalent to the possible nonzeroness of
        //    matrix entries A(I,J) and A(J,I).
        //
        //    A sparse compressed column format is used to store the counts for
        //    the nonzeroes.  In other words, while the value ADJ_NUM reports the
        //    number of adjacencies, the vector ADJ_COL is sufficient to allow us
        //    to properly set up a sparse compressed matrix for the actual storage
        //    of the sparse matrix, if we desire to proceed.
        //
        //  Local Node Numbering:
        //
        //       3
        //    s  |.
        //    i  | .
        //    d  |  .
        //    e  6   5  side 2
        //       |    .
        //    3  |     .
        //       |      .
        //       1---4---2
        //
        //         side 1
        //
        //  Variable Diagram:
        //
        //      UVP
        //       |.
        //       | .
        //       |  .
        //      UV   UV
        //       |    .
        //       |     .
        //       |      .
        //      UVP--UV--UVP
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int TRIANGLE_NUM, the number of triangles.
        //
        //    Input, int VARIABLE_NUM, the number of variables.
        //
        //    Input, int TRIANGLE_NODE[6*TRIANGLE_NUM], lists the nodes that
        //    make up each triangle.  The first three nodes are the vertices,
        //    in counterclockwise order.  The fourth value is the midside
        //    node between nodes 1 and 2; the fifth and sixth values are
        //    the other midside nodes in the logical order.
        //
        //    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], for each side of
        //    a triangle, lists the neighboring triangle, or -1 if there is
        //    no neighbor.
        //
        //    Input, int NODE_U_VARIABLE[NODE_NUM], NODE_V_VARIABLE[NODE_NUM],
        //    NODE_P_VARIABLE[NODE_NUM], the index of the horizontal velocity,
        //    vertical velocity and pressure variables associated with a node,
        //    or -1 if no such variable is associated with the node.
        //
        //    Output, int ADJ_COL[VARIABLE_NUM+1].  Information about variable J
        //    is stored in entries ADJ_COL[J] through ADJ_COL[J+1]-1 of ADJ.
        //
        //    Output, int NS_ADJ_COL_SET, the number of Navier Stokes variable
        //    adjacencies.
        //
        {
            int adj_num;
            int n1;
            int n2;
            int n3;
            int n4;
            int n5;
            int n6;
            int node;
            int p1;
            int p2;
            int p3;
            int triangle;
            int triangle_order = 6;
            int triangle2;
            int u1;
            int u2;
            int u3;
            int u4;
            int u5;
            int u6;
            int v1;
            int v2;
            int v3;
            int v4;
            int v5;
            int v6;
            int variable;

            adj_num = 0;
            //
            //  Set every variable to be adjacent to itself.
            //
            for (variable = 0; variable < variable_num; variable++)
            {
                adj_col[variable] = 1;
            }

            //
            //  Set every variable to be adjacent to the other variables associated with
            //  that node.
            //
            //  U <=> V
            //  U <=> P (if there is a P variable)
            //  V <=> P (if there is a P variable)
            //
            for (node = 0; node < node_num; node++)
            {
                u1 = node_u_variable[node] - 1;
                v1 = node_v_variable[node] - 1;
                p1 = node_p_variable[node] - 1;

                adj_col[u1] = adj_col[u1] + 1;
                adj_col[v1] = adj_col[v1] + 1;

                if (0 <= p1)
                {
                    adj_col[u1] = adj_col[u1] + 1;
                    adj_col[v1] = adj_col[v1] + 1;
                    adj_col[p1] = adj_col[p1] + 2;
                }
            }

            //
            //  Examine each triangle.
            //
            for (triangle = 0; triangle < triangle_num; triangle++)
            {
                n1 = triangle_node[0 + triangle * triangle_order] - 1;
                n2 = triangle_node[1 + triangle * triangle_order] - 1;
                n3 = triangle_node[2 + triangle * triangle_order] - 1;
                n4 = triangle_node[3 + triangle * triangle_order] - 1;
                n5 = triangle_node[4 + triangle * triangle_order] - 1;
                n6 = triangle_node[5 + triangle * triangle_order] - 1;

                u1 = node_u_variable[n1] - 1;
                v1 = node_v_variable[n1] - 1;
                p1 = node_p_variable[n1] - 1;

                u2 = node_u_variable[n2] - 1;
                v2 = node_v_variable[n2] - 1;
                p2 = node_p_variable[n2] - 1;

                u3 = node_u_variable[n3] - 1;
                v3 = node_v_variable[n3] - 1;
                p3 = node_p_variable[n3] - 1;

                u4 = node_u_variable[n4] - 1;
                v4 = node_v_variable[n4] - 1;

                u5 = node_u_variable[n5] - 1;
                v5 = node_v_variable[n5] - 1;

                u6 = node_u_variable[n6] - 1;
                v6 = node_v_variable[n6] - 1;
                //
                //  For sure, we add the new adjacencies:
                //
                //    U5 V5 <=> U1 V1 P1
                //    U6 V6 <=> U2 V2 P2
                //    U4 V4 <=> U3 V3 P3
                //    U5 V5 <=> U4 V4
                //    U6 V6 <=> U4 V4
                //    U6 V6 <=> U5 V5
                //
                adj_col[u1] = adj_col[u1] + 2;
                adj_col[v1] = adj_col[v1] + 2;
                adj_col[p1] = adj_col[p1] + 2;

                adj_col[u2] = adj_col[u2] + 2;
                adj_col[v2] = adj_col[v2] + 2;
                adj_col[p2] = adj_col[p2] + 2;

                adj_col[u3] = adj_col[u3] + 2;
                adj_col[v3] = adj_col[v3] + 2;
                adj_col[p3] = adj_col[p3] + 2;

                adj_col[u4] = adj_col[u4] + 7;
                adj_col[v4] = adj_col[v4] + 7;

                adj_col[u5] = adj_col[u5] + 7;
                adj_col[v5] = adj_col[v5] + 7;

                adj_col[u6] = adj_col[u6] + 7;
                adj_col[v6] = adj_col[v6] + 7;
                //
                //  Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
                //  that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
                //  or if this triangle is the first of the pair in which the edge
                //  occurs (TRIANGLE < TRIANGLE2).
                //
                //  Maybe add
                //
                //    U1 V1 P1 <=> U2 V2 P2
                //    U1 V1 P1 <=> U4 V4
                //    U2 V2 P2 <=> U4 V4
                //
                triangle2 = triangle_neighbor[0 + triangle * 3];

                if (triangle2 < 0 || triangle < triangle2)
                {
                    adj_col[u1] = adj_col[u1] + 5;
                    adj_col[v1] = adj_col[v1] + 5;
                    adj_col[p1] = adj_col[p1] + 5;

                    adj_col[u2] = adj_col[u2] + 5;
                    adj_col[v2] = adj_col[v2] + 5;
                    adj_col[p2] = adj_col[p2] + 5;

                    adj_col[u4] = adj_col[u4] + 6;
                    adj_col[v4] = adj_col[v4] + 6;
                }

                //
                //  Maybe add
                //
                //    U2 V2 P2 <=> U3 V3 P3
                //    U2 V2 P2 <=> U5 V5
                //    U3 V3 P3 <=> U5 V5
                //
                triangle2 = triangle_neighbor[1 + triangle * 3];

                if (triangle2 < 0 || triangle < triangle2)
                {
                    adj_col[u2] = adj_col[u2] + 5;
                    adj_col[v2] = adj_col[v2] + 5;
                    adj_col[p2] = adj_col[p2] + 5;

                    adj_col[u3] = adj_col[u3] + 5;
                    adj_col[v3] = adj_col[v3] + 5;
                    adj_col[p3] = adj_col[p3] + 5;

                    adj_col[u5] = adj_col[u5] + 6;
                    adj_col[v5] = adj_col[v5] + 6;
                }

                //
                //  Maybe add
                //
                //    U1 V1 P1 <=> U3 V3 P3
                //    U1 V1 P1 <=> U6 V6
                //    U3 V3 P3 <=> U6 V6
                //
                triangle2 = triangle_neighbor[2 + triangle * 3];

                if (triangle2 < 0 || triangle < triangle2)
                {
                    adj_col[u1] = adj_col[u1] + 5;
                    adj_col[v1] = adj_col[v1] + 5;
                    adj_col[p1] = adj_col[p1] + 5;

                    adj_col[u3] = adj_col[u3] + 5;
                    adj_col[v3] = adj_col[v3] + 5;
                    adj_col[p3] = adj_col[p3] + 5;

                    adj_col[u6] = adj_col[u6] + 6;
                    adj_col[v6] = adj_col[v6] + 6;
                }

            }

            //
            //  We used ADJ_COL to count the number of entries in each column.
            //  Convert it to pointers into the ADJ array.
            //
            for (variable = variable_num; 0 < variable; variable--)
            {
                adj_col[variable] = adj_col[variable - 1];
            }

            adj_col[0] = 1;
            for (variable = 1; variable <= variable_num; variable++)
            {
                adj_col[variable] = adj_col[variable - 1] + adj_col[variable];
            }

            adj_num = adj_col[variable_num] - 1;

            return adj_num;
        }

        public static int ns_adj_count(int node_num, int triangle_num, int variable_num,
            int[] triangle_node, int[] triangle_neighbor, int[] node_u_variable,
        int[] node_v_variable, int[] node_p_variable )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NS_ADJ_COUNT counts adjacencies in a Navier Stokes triangulation.
        //
        //  Discussion:
        //
        //    This routine is called to count the adjacencies, so that the
        //    appropriate amount of memory can be set aside for storage when
        //    the adjacency structure is created.
        //
        //    The value of ADJ_NUM computed and returned by this routine should
        //    be identical to the value computed by NS_ADJ_COL_SET.
        //
        //    The triangulation is assumed to involve 6-node triangles.
        //
        //    Variables for the horizontal and vertical velocities are associated
        //    with every node.  Variables for the pressure are associated only with
        //    the vertex nodes.
        //
        //    We are interested in determining the number of nonzero entries in the
        //    stiffness matrix of the Stokes equations, or the jacobian matrix of
        //    the Navier Stokes equations.  To this end, we will say, somewhat
        //    too broadly, that two variables are "adjacent" if their associated
        //    nodes both occur in some common element.  This adjacency of variables
        //    I and J is taken to be equivalent to the possible nonzeroness of
        //    matrix entries A(I,J) and A(J,I).
        //
        //    A sparse compressed column format is used to store the counts for
        //    the nonzeroes.  In other words, while the value ADJ_NUM reports the
        //    number of adjacencies, the vector ADJ_COL is sufficient to allow us
        //    to properly set up a sparse compressed matrix for the actual storage
        //    of the sparse matrix, if we desire to proceed.
        //
        //  Local Node Numbering:
        //
        //       3
        //    s  |.
        //    i  | .
        //    d  |  .
        //    e  6   5  side 2
        //       |    .
        //    3  |     .
        //       |      .
        //       1---4---2
        //
        //         side 1
        //
        //  Variable Diagram:
        //
        //      UVP
        //       |.
        //       | .
        //       |  .
        //      UV   UV
        //       |    .
        //       |     .
        //       |      .
        //      UVP--UV--UVP
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int TRIANGLE_NUM, the number of triangles.
        //
        //    Input, int VARIABLE_NUM, the number of variables.
        //
        //    Input, int TRIANGLE_NODE[6*TRIANGLE_NUM], lists the nodes that
        //    make up each triangle.  The first three nodes are the vertices,
        //    in counterclockwise order.  The fourth value is the midside
        //    node between nodes 1 and 2; the fifth and sixth values are
        //    the other midside nodes in the logical order.
        //
        //    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], for each side of
        //    a triangle, lists the neighboring triangle, or -1 if there is
        //    no neighbor.
        //
        //    Input, int NODE_U_VARIABLE[NODE_NUM], NODE_V_VARIABLE[NODE_NUM],
        //    NODE_P_VARIABLE[NODE_NUM], the index of the horizontal velocity,
        //    vertical velocity and pressure variables associated with a node,
        //    or -1 if no such variable is associated with the node.
        //
        //    Output, int NS_ADJ_COUNT, the value of ADJ_NUM, the number of
        //    Navier Stokes variable adjacencies.
        //
        {
            int adj_num;
            int node;
            int p1;
            int triangle;
            int triangle2;

            adj_num = 0;
            //
            //  Set every variable to be adjacent to itself.
            //
            adj_num = variable_num;
            //
            //  Set every variable to be adjacent to the other variables associated with
            //  that node.
            //
            //  U <=> V
            //  U <=> P (if there is a P variable)
            //  V <=> P (if there is a P variable)
            //
            for (node = 0; node < node_num; node++)
            {
                adj_num = adj_num + 2;

                p1 = node_p_variable[node] - 1;

                if (0 <= p1)
                {
                    adj_num = adj_num + 4;
                }
            }

            //
            //  Examine each triangle.
            //
            for (triangle = 0; triangle < triangle_num; triangle++)
            {

                //
                //  For sure, we add the new adjacencies:
                //
                //    U5 V5 <=> U1 V1 P1
                //    U6 V6 <=> U2 V2 P2
                //    U4 V4 <=> U3 V3 P3
                //    U5 V5 <=> U4 V4
                //    U6 V6 <=> U4 V4
                //    U6 V6 <=> U5 V5
                //
                adj_num = adj_num + 60;
                //
                //  Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
                //  that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
                //  or if this triangle is the first of the pair in which the edge
                //  occurs (TRIANGLE < TRIANGLE2).
                //
                //  Maybe add
                //
                //    U1 V1 P1 <=> U2 V2 P2
                //    U1 V1 P1 <=> U4 V4
                //    U2 V2 P2 <=> U4 V4
                //
                triangle2 = triangle_neighbor[0 + triangle * 3];

                if (triangle2 < 0 || triangle < triangle2)
                {
                    adj_num = adj_num + 42;
                }

                //
                //  Maybe add
                //
                //    U2 V2 P2 <=> U3 V3 P3
                //    U2 V2 P2 <=> U5 V5
                //    U3 V3 P3 <=> U5 V5
                //
                triangle2 = triangle_neighbor[1 + triangle * 3];

                if (triangle2 < 0 || triangle < triangle2)
                {
                    adj_num = adj_num + 42;
                }

                //
                //  Maybe add
                //
                //    U1 V1 P1 <=> U3 V3 P3
                //    U1 V1 P1 <=> U6 V6
                //    U3 V3 P3 <=> U6 V6
                //
                triangle2 = triangle_neighbor[2 + triangle * 3];

                if (triangle2 < 0 || triangle < triangle2)
                {
                    adj_num = adj_num + 42;
                }

            }

            return adj_num;
        }

        public static void ns_adj_insert(int v1, int v2, int variable_num, int adj_num,
            ref int[] adj_col_free, ref int[] adj_row )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NS_ADJ_INSERT inserts an adjacency into a compressed column adjacency matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int V1, V2, the indices of two items which are adjacent.
        //
        //    Input, int VARIABLE_NUM, the number of items.
        //
        //    Input, int ADJ_NUM, the number of entries available in ADJ_ROW.
        //
        //    Input/output, int ADJ_COL_FREE[VARIABLE_NUM], contains the next free
        //    location in which an entry for a given column can be stored.  On output,
        //    two pointers have been updated.
        //
        //    Input/output, int ADJ_ROW[ADJ_NUM], the row indices of the Navier Stokes
        //    variable adjacency matrix.  On output, two new entries have been added.
        //
        {
            adj_row[adj_col_free[v1 - 1] - 1] = v2;
            adj_col_free[v1 - 1] = adj_col_free[v1 - 1] + 1;

            if (v1 == v2)
            {
                return;
            }

            adj_row[adj_col_free[v2 - 1] - 1] = v1;
            adj_col_free[v2 - 1] = adj_col_free[v2 - 1] + 1;

            return;
        }

        public static void ns_adj_row_set(int node_num, int triangle_num, int variable_num,
            int[] triangle_node, int[] triangle_neighbor, int[] node_u_variable,
        int[] node_v_variable, int[] node_p_variable, int adj_num, int[] adj_col,
        int[] adj_row )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NS_ADJ_ROW_SET sets the Navier Stokes sparse compressed column row indices.
        //
        //  Discussion:
        //
        //    After NS_ADJ_COUNT has been called to count ADJ_NUM, the number of
        //    variable adjacencies and to set up ADJ_COL, the compressed column pointer,
        //    this routine can be called to assign values to ADJ_ROW, the row
        //    indices of the sparse compressed column adjacency matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int TRIANGLE_NUM, the number of triangles.
        //
        //    Input, int VARIABLE_NUM, the number of variables.
        //
        //    Input, int TRIANGLE_NODE[6*TRIANGLE_NUM], lists the nodes that
        //    make up each triangle.  The first three nodes are the vertices,
        //    in counterclockwise order.  The fourth value is the midside
        //    node between nodes 1 and 2; the fifth and sixth values are
        //    the other midside nodes in the logical order.
        //
        //    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], for each side of
        //    a triangle, lists the neighboring triangle, or -1 if there is
        //    no neighbor.
        //
        //    Input, int NODE_U_VARIABLE[NODE_NUM], NODE_V_VARIABLE[NODE_NUM],
        //    NODE_P_VARIABLE[NODE_NUM], the index of the horizontal velocity,
        //    vertical velocity and pressure variables associated with a node,
        //    or -1 if no such variable is associated with the node.
        //
        //    Input, int ADJ_NUM, the number of Navier Stokes variable adjacencies.
        //
        //    Input, int ADJ_COL[VARIABLE_NUM+1].  Information about variable J
        //    is stored in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
        //
        //    Output, int ADJ_ROW[ADJ_NUM], the row indices of the Navier Stokes
        //    variable adjacency matrix.
        //
        //  Local Parameters:
        //
        //    Local, int ADJ_COL_FREE[VARIABLE_NUM], for each column,
        //    the location in ADJ_ROW which can store the next row index.
        //
        {
            int[] adj_col_free;
            int k1;
            int k2;
            int n1;
            int n2;
            int n3;
            int n4;
            int n5;
            int n6;
            int node;
            int p1;
            int p2;
            int p3;
            int triangle;
            int triangle2;
            int u1;
            int u2;
            int u3;
            int u4;
            int u5;
            int u6;
            int v;
            int v1;
            int v2;
            int v3;
            int v4;
            int v5;
            int v6;

            for (v = 0; v < adj_num; v++)
            {
                adj_row[v] = -1;
            }

            adj_col_free = new int[variable_num];

            for (v = 0; v < variable_num; v++)
            {
                adj_col_free[v] = adj_col[v];
            }

            //
            //  Set every variable to be adjacent to itself.
            //  Here, we have to be careful to start at index 1.
            //
            for (v = 1; v <= variable_num; v++)
            {
                ns_adj_insert(v, v, variable_num, adj_num, ref adj_col_free, ref adj_row);
            }

            //
            //  Set every variable to be adjacent to the other variables associated with
            //  that node.
            //
            //  U <=> V
            //  U <=> P (if there is a P variable)
            //  V <=> P (if there is a P variable)
            //
            for (node = 0; node < node_num; node++)
            {
                u1 = node_u_variable[node];
                v1 = node_v_variable[node];
                p1 = node_p_variable[node];

                ns_adj_insert(u1, v1, variable_num, adj_num, ref adj_col_free, ref adj_row);

                if (0 < p1)
                {
                    ns_adj_insert(u1, p1, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(v1, p1, variable_num, adj_num, ref adj_col_free, ref adj_row);
                }
            }

            //
            //  Examine each triangle.
            //
            for (triangle = 0; triangle < triangle_num; triangle++)
            {
                n1 = triangle_node[0 + triangle * 6];
                n2 = triangle_node[1 + triangle * 6];
                n3 = triangle_node[2 + triangle * 6];
                n4 = triangle_node[3 + triangle * 6];
                n5 = triangle_node[4 + triangle * 6];
                n6 = triangle_node[5 + triangle * 6];

                u1 = node_u_variable[n1 - 1];
                v1 = node_v_variable[n1 - 1];
                p1 = node_p_variable[n1 - 1];

                u2 = node_u_variable[n2 - 1];
                v2 = node_v_variable[n2 - 1];
                p2 = node_p_variable[n2 - 1];

                u3 = node_u_variable[n3 - 1];
                v3 = node_v_variable[n3 - 1];
                p3 = node_p_variable[n3 - 1];

                u4 = node_u_variable[n4 - 1];
                v4 = node_v_variable[n4 - 1];

                u5 = node_u_variable[n5 - 1];
                v5 = node_v_variable[n5 - 1];

                u6 = node_u_variable[n6 - 1];
                v6 = node_v_variable[n6 - 1];
                //
                //  For sure, we add the new adjacencies:
                //
                //    U5 V5 <=> U1 V1 P1
                //    U6 V6 <=> U2 V2 P2
                //    U4 V4 <=> U3 V3 P3
                //    U5 V5 <=> U4 V4
                //    U6 V6 <=> U4 V4
                //    U6 V6 <=> U5 V5
                //
                ns_adj_insert(u5, u1, variable_num, adj_num, ref adj_col_free, ref adj_row);
                ns_adj_insert(u5, v1, variable_num, adj_num, ref adj_col_free, ref adj_row);
                ns_adj_insert(u5, p1, variable_num, adj_num, ref adj_col_free, ref adj_row);
                ns_adj_insert(v5, u1, variable_num, adj_num, ref adj_col_free, ref adj_row);
                ns_adj_insert(v5, v1, variable_num, adj_num, ref adj_col_free, ref adj_row);
                ns_adj_insert(v5, p1, variable_num, adj_num, ref adj_col_free, ref adj_row);

                ns_adj_insert(u6, u2, variable_num, adj_num, ref adj_col_free, ref adj_row);
                ns_adj_insert(u6, v2, variable_num, adj_num, ref adj_col_free, ref adj_row);
                ns_adj_insert(u6, p2, variable_num, adj_num, ref adj_col_free, ref adj_row);
                ns_adj_insert(v6, u2, variable_num, adj_num, ref adj_col_free, ref adj_row);
                ns_adj_insert(v6, v2, variable_num, adj_num, ref adj_col_free, ref adj_row);
                ns_adj_insert(v6, p2, variable_num, adj_num, ref adj_col_free, ref adj_row);

                ns_adj_insert(u4, u3, variable_num, adj_num, ref adj_col_free, ref adj_row);
                ns_adj_insert(u4, v3, variable_num, adj_num, ref adj_col_free, ref adj_row);
                ns_adj_insert(u4, p3, variable_num, adj_num, ref adj_col_free, ref adj_row);
                ns_adj_insert(v4, u3, variable_num, adj_num, ref adj_col_free, ref adj_row);
                ns_adj_insert(v4, v3, variable_num, adj_num, ref adj_col_free, ref adj_row);
                ns_adj_insert(v4, p3, variable_num, adj_num, ref adj_col_free, ref adj_row);

                ns_adj_insert(u5, u4, variable_num, adj_num, ref adj_col_free, ref adj_row);
                ns_adj_insert(u5, v4, variable_num, adj_num, ref adj_col_free, ref adj_row);
                ns_adj_insert(v5, u4, variable_num, adj_num, ref adj_col_free, ref adj_row);
                ns_adj_insert(v5, v4, variable_num, adj_num, ref adj_col_free, ref adj_row);

                ns_adj_insert(u6, u4, variable_num, adj_num, ref adj_col_free, ref adj_row);
                ns_adj_insert(u6, v4, variable_num, adj_num, ref adj_col_free, ref adj_row);
                ns_adj_insert(v6, u4, variable_num, adj_num, ref adj_col_free, ref adj_row);
                ns_adj_insert(v6, v4, variable_num, adj_num, ref adj_col_free, ref adj_row);

                ns_adj_insert(u6, u5, variable_num, adj_num, ref adj_col_free, ref adj_row);
                ns_adj_insert(u6, v5, variable_num, adj_num, ref adj_col_free, ref adj_row);
                ns_adj_insert(v6, u5, variable_num, adj_num, ref adj_col_free, ref adj_row);
                ns_adj_insert(v6, v5, variable_num, adj_num, ref adj_col_free, ref adj_row);
                //
                //  Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
                //  that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
                //  or if this triangle is the first of the pair in which the edge
                //  occurs (TRIANGLE < TRIANGLE2).
                //
                //  Maybe add
                //
                //    U1 V1 P1 <=> U2 V2 P2
                //    U1 V1 P1 <=> U4 V4
                //    U2 V2 P2 <=> U4 V4
                //
                triangle2 = triangle_neighbor[0 + triangle * 3];

                if (triangle2 < 0 || triangle < triangle2)
                {
                    ns_adj_insert(u1, u2, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(u1, v2, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(u1, p2, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(v1, u2, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(v1, v2, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(v1, p2, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(p1, u2, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(p1, v2, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(p1, p2, variable_num, adj_num, ref adj_col_free, ref adj_row);

                    ns_adj_insert(u1, u4, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(u1, v4, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(v1, u4, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(v1, v4, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(p1, u4, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(p1, v4, variable_num, adj_num, ref adj_col_free, ref adj_row);

                    ns_adj_insert(u2, u4, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(u2, v4, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(v2, u4, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(v2, v4, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(p2, u4, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(p2, v4, variable_num, adj_num, ref adj_col_free, ref adj_row);
                }

                //
                //  Maybe add
                //
                //    U2 V2 P2 <=> U3 V3 P3
                //    U2 V2 P2 <=> U5 V5
                //    U3 V3 P3 <=> U5 V5
                //
                triangle2 = triangle_neighbor[1 + triangle * 3];

                if (triangle2 < 0 || triangle < triangle2)
                {
                    ns_adj_insert(u2, u3, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(u2, v3, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(u2, p3, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(v2, u3, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(v2, v3, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(v2, p3, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(p2, u3, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(p2, v3, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(p2, p3, variable_num, adj_num, ref adj_col_free, ref adj_row);

                    ns_adj_insert(u2, u5, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(u2, v5, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(v2, u5, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(v2, v5, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(p2, u5, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(p2, v5, variable_num, adj_num, ref adj_col_free, ref adj_row);

                    ns_adj_insert(u3, u5, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(u3, v5, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(v3, u5, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(v3, v5, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(p3, u5, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(p3, v5, variable_num, adj_num, ref adj_col_free, ref adj_row);
                }

                //
                //  Maybe add
                //
                //    U1 V1 P1 <=> U3 V3 P3
                //    U1 V1 P1 <=> U6 V6
                //    U3 V3 P3 <=> U6 V6
                //
                triangle2 = triangle_neighbor[2 + triangle * 3];

                if (triangle2 < 0 || triangle < triangle2)
                {
                    ns_adj_insert(u1, u3, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(u1, v3, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(u1, p3, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(v1, u3, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(v1, v3, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(v1, p3, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(p1, u3, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(p1, v3, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(p1, p3, variable_num, adj_num, ref adj_col_free, ref adj_row);

                    ns_adj_insert(u1, u6, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(u1, v6, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(v1, u6, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(v1, v6, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(p1, u6, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(p1, v6, variable_num, adj_num, ref adj_col_free, ref adj_row);

                    ns_adj_insert(u3, u6, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(u3, v6, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(v3, u6, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(v3, v6, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(p3, u6, variable_num, adj_num, ref adj_col_free, ref adj_row);
                    ns_adj_insert(p3, v6, variable_num, adj_num, ref adj_col_free, ref adj_row);
                }
            }

            //
            //  Ascending sort the entries for each variable.
            //
            for (v = 0; v < variable_num; v++)
            {
                k1 = adj_col[v];
                k2 = adj_col[v + 1] - 1;

                typeMethods.i4vec_sort_heap_a(k2 + 1 - k1, ref adj_row, aIndex: + k1 - 1);
            }

            return;
        }
    }
}