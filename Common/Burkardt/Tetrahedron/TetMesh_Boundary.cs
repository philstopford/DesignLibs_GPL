using Burkardt.Types;

namespace Burkardt.TetrahedronNS;

public static class TetMesh_Boundary
{
    public static void tet_mesh_boundary_count(int element_order, int element_num,
            int[] element_node, int node_num, ref int boundary_node_num,
            ref int boundary_element_num, ref int[] boundary_node_mask )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TET_MESH_BOUNDARY_COUNT counts boundary faces and nodes in a tet mesh.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 December 2010
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
        //    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the
        //    nodes that make up each element.
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Output, int *BOUNDARY_NODE_NUM, the number of boundary nodes.
        //
        //    Output, int *BOUNDARY_ELEMENT_NUM, the number of boundary elements.
        //
        //    Output, int BOUNDARY_NODE_MASK[NODE_NUM], is 0 for
        //    interior nodes, 1 for boundary nodes.
        //
    {
        int a = 0;
        int b = 0;
        int c= 0;
        int element;
        int f;
        int face;
        int[] faces;
        int i;
        int j;
        int jj;
        int k;
        int l;

        faces = new int[5 * 4 * element_num];
        //
        //  Step 1.
        //  From the list of nodes forming tetrahedron T, of the form: 
        //
        //    (I,J,K,L)
        //
        //  or
        //
        //    (I,J,K,L,I+J,I+K,I+L,J+K,J+L,K+L),
        //    (1,2,3,4, 5,  6,  7,  8,  9, 10 ),
        //
        //  construct the four face relations:
        //
        //    F = 1: (J,K,L,F,T)
        //    F = 2: (I,K,L,F,T)
        //    F = 3: (I,J,L,F,T)
        //    F = 4: (I,J,K,F,T)
        //
        //  If T is actually order 10, we can retrieve the indices of the midside
        //  nodes from the values of F and T.  In that case, the 4 faces are:
        //
        //    F = 1: 2, 3, 4, 8, 10, 9
        //    F = 2: 1, 3, 4, 6, 10, 7
        //    F = 3: 1, 2, 4, 5,  9, 7
        //    F = 4: 1, 2, 3, 5,  8, 6
        //
        //  In order to make matching easier, we reorder each triple of nodes
        //  into ascending order.
        //
        for (element = 0; element < element_num; element++)
        {
            i = element_node[0 + element * element_order];
            j = element_node[1 + element * element_order];
            k = element_node[2 + element * element_order];
            l = element_node[3 + element * element_order];

            typeMethods.i4i4i4_sort_a(j, k, l, ref a, ref b, ref c);

            jj = 4 * element;
            faces[0 + jj * 5] = a;
            faces[1 + jj * 5] = b;
            faces[2 + jj * 5] = c;
            faces[3 + jj * 5] = 0;
            faces[4 + jj * 5] = element;

            typeMethods.i4i4i4_sort_a(i, k, l, ref a, ref b, ref c);

            jj = 4 * element + 1;
            faces[0 + jj * 5] = a;
            faces[1 + jj * 5] = b;
            faces[2 + jj * 5] = c;
            faces[3 + jj * 5] = 1;
            faces[4 + jj * 5] = element;

            typeMethods.i4i4i4_sort_a(i, j, l, ref a, ref b, ref c);

            jj = 4 * element + 2;
            faces[0 + jj * 5] = a;
            faces[1 + jj * 5] = b;
            faces[2 + jj * 5] = c;
            faces[3 + jj * 5] = 2;
            faces[4 + jj * 5] = element;

            typeMethods.i4i4i4_sort_a(i, j, k, ref a, ref b, ref c);

            jj = 4 * element + 3;
            faces[0 + jj * 5] = a;
            faces[1 + jj * 5] = b;
            faces[2 + jj * 5] = c;
            faces[3 + jj * 5] = 3;
            faces[4 + jj * 5] = element;
        }

        //
        //  Step 2. Perform an ascending dictionary sort on the neighbor relations.
        //  We only intend to sort on rows 1:3; the routine we call here
        //  sorts on rows 1 through 5 but that won't hurt us.
        //
        //  What we need is to find cases where two tetrahedrons share a face.
        //  By sorting the columns of the FACES array, we will put shared faces
        //  next to each other.
        //
        typeMethods.i4col_sort_a(5, 4 * element_num, ref faces);
        //
        //  Step 3. Neighboring faces show up as consecutive columns with
        //  identical first three entries.  Count columns which don't have
        //  a following column that matches the first three entries.
        //
        boundary_element_num = 0;
        face = 0;

        for (i = 0; i < node_num; i++)
        {
            boundary_node_mask[i] = 0;
        }

        while (face < 4 * element_num)
        {
            if (face < 4 * element_num - 1)
            {
                if (faces[0 + face * 5] == faces[0 + (face + 1) * 5] &&
                    faces[1 + face * 5] == faces[1 + (face + 1) * 5] &&
                    faces[2 + face * 5] == faces[2 + (face + 1) * 5])
                {
                    face += 2;
                    continue;
                }
            }

            boundary_element_num += 1;
            //
            //  The vertices of the triangle are boundary nodes.
            //
            boundary_node_mask[faces[0 + face * 5]] = 1;
            boundary_node_mask[faces[1 + face * 5]] = 1;
            boundary_node_mask[faces[2 + face * 5]] = 1;
            switch (element_order)
            {
                //
                //  For quadratic tetrahedrons, we need to add three more side nodes.
                //  We retrieve the face index by F = FACES(4,*).
                //  We can determine the local indices from the value of F.
                //  We can determine the global indices from ELEMENT_NODE.
                //
                case 10:
                    f = faces[3 + face * 5];
                    element = faces[4 + face * 5];

                    switch (f)
                    {
                        case 0:
                            boundary_node_mask[element_node[7 + element * element_order]] = 1;
                            boundary_node_mask[element_node[9 + element * element_order]] = 1;
                            boundary_node_mask[element_node[8 + element * element_order]] = 1;
                            break;
                        case 1:
                            boundary_node_mask[element_node[5 + element * element_order]] = 1;
                            boundary_node_mask[element_node[9 + element * element_order]] = 1;
                            boundary_node_mask[element_node[6 + element * element_order]] = 1;
                            break;
                        case 2:
                            boundary_node_mask[element_node[4 + element * element_order]] = 1;
                            boundary_node_mask[element_node[8 + element * element_order]] = 1;
                            boundary_node_mask[element_node[6 + element * element_order]] = 1;
                            break;
                        case 3:
                            boundary_node_mask[element_node[4 + element * element_order]] = 1;
                            boundary_node_mask[element_node[7 + element * element_order]] = 1;
                            boundary_node_mask[element_node[5 + element * element_order]] = 1;
                            break;
                    }

                    break;
            }

            face += 1;
        }

        boundary_node_num = typeMethods.i4vec_sum(node_num, boundary_node_mask);
    }

    public static int[] tet_mesh_boundary_set(int element_order, int element_num,
            int[] element_node, int boundary_element_order, int boundary_element_num )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TET_MESH_BOUNDARY_SET sets the boundary faces in a tet mesh.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 December 2010
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
        //    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the
        //    nodes that make up each element.
        //
        //    Input, int BOUNDARY_ELEMENT_ORDER, the order of the boundary faces.
        //
        //    Input, int BOUNDARY_ELEMENT_NUM, the number of boundary faces.
        //
        //    Output, int 
        //    TET_MESH_BOUNDARY_SET[BOUNDARY_ELEMENT_ORDER*BOUNDARY_ELEMENT_NUM],
        //    the nodes that make up each boundary face.
        //
    {
        int a = 0;
        int b = 0;
        int boundary_element;
        int[] boundary_element_node;
        int c = 0;
        int element;
        int f;
        int face;
        int[] faces;
        int i;
        int j;
        int jj;
        int k;
        int l;

        faces = new int[5 * 4 * element_num];
        //
        //  Step 1.
        //  From the list of nodes forming tetrahedron T, of the form: 
        //
        //    (I,J,K,L)
        //
        //  or
        //
        //    (I,J,K,L,I+J,I+K,I+L,J+K,J+L,K+L),
        //    (1,2,3,4, 5,  6,  7,  8,  9, 10 ),
        //
        //  construct the four face relations:
        //
        //    F = 1: (J,K,L,F,T)
        //    F = 2: (I,K,L,F,T)
        //    F = 3: (I,J,L,F,T)
        //    F = 4: (I,J,K,F,T)
        //
        //  If T is actually order 10, we can retrieve the indices of the midside
        //  nodes from the values of F and T.  In that case, the 4 faces are:
        //
        //    F = 1: 2, 3, 4, 8, 10, 9
        //    F = 2: 1, 3, 4, 6, 10, 7
        //    F = 3: 1, 2, 4, 5,  9, 7
        //    F = 4: 1, 2, 3, 5,  8, 6
        //
        //  In order to make matching easier, we reorder each triple of nodes
        //  into ascending order.
        //
        for (element = 0; element < element_num; element++)
        {
            i = element_node[0 + element * element_order];
            j = element_node[1 + element * element_order];
            k = element_node[2 + element * element_order];
            l = element_node[3 + element * element_order];

            typeMethods.i4i4i4_sort_a(j, k, l, ref a, ref b, ref c);

            jj = 4 * element;
            faces[0 + jj * 5] = a;
            faces[1 + jj * 5] = b;
            faces[2 + jj * 5] = c;
            faces[3 + jj * 5] = 0;
            faces[4 + jj * 5] = element;

            typeMethods.i4i4i4_sort_a(i, k, l, ref a, ref b, ref c);

            jj = 4 * element + 1;
            faces[0 + jj * 5] = a;
            faces[1 + jj * 5] = b;
            faces[2 + jj * 5] = c;
            faces[3 + jj * 5] = 1;
            faces[4 + jj * 5] = element;

            typeMethods.i4i4i4_sort_a(i, j, l, ref a, ref b, ref c);

            jj = 4 * element + 2;
            faces[0 + jj * 5] = a;
            faces[1 + jj * 5] = b;
            faces[2 + jj * 5] = c;
            faces[3 + jj * 5] = 2;
            faces[4 + jj * 5] = element;

            typeMethods.i4i4i4_sort_a(i, j, k, ref a, ref b, ref c);

            jj = 4 * element + 3;
            faces[0 + jj * 5] = a;
            faces[1 + jj * 5] = b;
            faces[2 + jj * 5] = c;
            faces[3 + jj * 5] = 3;
            faces[4 + jj * 5] = element;
        }

        //
        //  Step 2. Perform an ascending dictionary sort on the neighbor relations.
        //  We only intend to sort on rows 1:3; the routine we call here
        //  sorts on rows 1 through 5 but that won't hurt us.
        //
        //  What we need is to find cases where two tetrahedrons share a face.
        //  By sorting the columns of the FACES array, we will put shared faces
        //  next to each other.
        //
        typeMethods.i4col_sort_a(5, 4 * element_num, ref faces);
        //
        //  Step 3. Neighboring faces show up as consecutive columns with
        //  identical first three entries.  Count columns which don't have
        //  a following column that matches the first three entries.
        //
        boundary_element = 0;
        face = 0;

        boundary_element_node = new int[boundary_element_order * boundary_element_num];

        while (face < 4 * element_num)
        {
            if (face < 4 * element_num - 1)
            {
                if (faces[0 + face * 5] == faces[0 + (face + 1) * 5] &&
                    faces[1 + face * 5] == faces[1 + (face + 1) * 5] &&
                    faces[2 + face * 5] == faces[2 + (face + 1) * 5])
                {
                    face += 2;
                    continue;
                }
            }

            f = faces[3 + face * 5];
            element = faces[4 + face * 5];

            switch (f)
            {
                case 0:
                    boundary_element_node[0 + boundary_element * boundary_element_order] =
                        element_node[1 + element * element_order];
                    boundary_element_node[1 + boundary_element * boundary_element_order] =
                        element_node[2 + element * element_order];
                    boundary_element_node[2 + boundary_element * boundary_element_order] =
                        element_node[3 + element * element_order];
                    break;
                case 1:
                    boundary_element_node[0 + boundary_element * boundary_element_order] =
                        element_node[0 + element * element_order];
                    boundary_element_node[1 + boundary_element * boundary_element_order] =
                        element_node[2 + element * element_order];
                    boundary_element_node[2 + boundary_element * boundary_element_order] =
                        element_node[3 + element * element_order];
                    break;
                case 2:
                    boundary_element_node[0 + boundary_element * boundary_element_order] =
                        element_node[0 + element * element_order];
                    boundary_element_node[1 + boundary_element * boundary_element_order] =
                        element_node[1 + element * element_order];
                    boundary_element_node[2 + boundary_element * boundary_element_order] =
                        element_node[3 + element * element_order];
                    break;
                case 3:
                    boundary_element_node[0 + boundary_element * boundary_element_order] =
                        element_node[0 + element * element_order];
                    boundary_element_node[1 + boundary_element * boundary_element_order] =
                        element_node[1 + element * element_order];
                    boundary_element_node[2 + boundary_element * boundary_element_order] =
                        element_node[2 + element * element_order];
                    break;
            }

            switch (element_order)
            {
                //
                //  For quadratic tetrahedrons, we need to add three more side nodes.
                //
                case 10:
                    switch (f)
                    {
                        case 0:
                            boundary_element_node[3 + boundary_element * boundary_element_order] =
                                element_node[7 + element * element_order];
                            boundary_element_node[4 + boundary_element * boundary_element_order] =
                                element_node[9 + element * element_order];
                            boundary_element_node[5 + boundary_element * boundary_element_order] =
                                element_node[8 + element * element_order];
                            break;
                        case 1:
                            boundary_element_node[3 + boundary_element * boundary_element_order] =
                                element_node[5 + element * element_order];
                            boundary_element_node[4 + boundary_element * boundary_element_order] =
                                element_node[9 + element * element_order];
                            boundary_element_node[5 + boundary_element * boundary_element_order] =
                                element_node[6 + element * element_order];
                            break;
                        case 2:
                            boundary_element_node[3 + boundary_element * boundary_element_order] =
                                element_node[4 + element * element_order];
                            boundary_element_node[4 + boundary_element * boundary_element_order] =
                                element_node[8 + element * element_order];
                            boundary_element_node[5 + boundary_element * boundary_element_order] =
                                element_node[6 + element * element_order];
                            break;
                        case 3:
                            boundary_element_node[3 + boundary_element * boundary_element_order] =
                                element_node[4 + element * element_order];
                            boundary_element_node[4 + boundary_element * boundary_element_order] =
                                element_node[7 + element * element_order];
                            boundary_element_node[5 + boundary_element * boundary_element_order] =
                                element_node[5 + element * element_order];
                            break;
                    }

                    break;
            }

            boundary_element += 1;
            face += 1;
        }

        return boundary_element_node;
    }
}