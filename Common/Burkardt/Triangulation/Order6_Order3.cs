namespace Burkardt.TriangulationNS;

public static partial class Conversion
{
    public static int[] triangulation_order6_to_order3(int triangle_num1, int[] triangle_node1)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_ORDER6_TO_ORDER3 linearizes a quadratic triangulation.
        //
        //  Discussion:
        //
        //    A quadratic triangulation is assumed to consist of 6-node triangles,
        //    as in the following:
        //
        //    11-12-13-14-15
        //     |\    |\    |
        //     | \   | \   |
        //     6  7  8  9 10
        //     |   \ |   \ |
        //     |    \|    \|
        //     1--2--3--4--5
        //
        //   This routine rearranges information so as to define the 3-node
        //   triangulation:
        //
        //    11-12-13-14-15
        //     |\ |\ |\ |\ |
        //     | \| \| \| \|
        //     6--7--8--9-10
        //     |\ |\ |\ |\ |
        //     | \| \| \| \|
        //     1--2--3--4--5
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int TRIANGLE_NUM1, the number of triangles in the quadratic
        //    triangulation.
        //
        //    Input, int TRIANGLE_NODE1[6*TRIANGLE_NUM1], the quadratic triangulation.
        //
        //    Output, int TRIANGULATION_ORDER6_TO_ORDER3[3*TRIANGLE_NUM2], the linear
        //    triangulation.  Here, TRIANGLE_NUM2 = 4 * TRIANGLE_NUM1.
        //
    {
        int n1;
        int n2;
        int n3;
        int n4;
        int n5;
        int n6;
        int triangle_num2;
        int tri1;
        int tri2;
        int[] triangle_node2;

        triangle_num2 = 4 * triangle_num1;
        triangle_node2 = new int[3 * triangle_num2];

        tri2 = 0;

        for (tri1 = 0; tri1 < triangle_num1; tri1++)
        {
            n1 = triangle_node1[0 + tri1 * 6];
            n2 = triangle_node1[1 + tri1 * 6];
            n3 = triangle_node1[2 + tri1 * 6];
            n4 = triangle_node1[3 + tri1 * 6];
            n5 = triangle_node1[4 + tri1 * 6];
            n6 = triangle_node1[5 + tri1 * 6];

            triangle_node2[0 + tri2 * 3] = n1;
            triangle_node2[1 + tri2 * 3] = n4;
            triangle_node2[2 + tri2 * 3] = n6;
            tri2 += 1;

            triangle_node2[0 + tri2 * 3] = n2;
            triangle_node2[1 + tri2 * 3] = n5;
            triangle_node2[2 + tri2 * 3] = n4;
            tri2 += 1;

            triangle_node2[0 + tri2 * 3] = n3;
            triangle_node2[1 + tri2 * 3] = n6;
            triangle_node2[2 + tri2 * 3] = n5;
            tri2 += 1;

            triangle_node2[0 + tri2 * 3] = n4;
            triangle_node2[1 + tri2 * 3] = n5;
            triangle_node2[2 + tri2 * 3] = n6;
            tri2 += 1;
        }

        return triangle_node2;
    }
}