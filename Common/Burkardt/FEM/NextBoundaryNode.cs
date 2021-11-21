namespace Burkardt.FEM;

public static class NextBoundaryNode
{
    public static int next_boundary_node(int node, string code)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    NEXT_BOUNDARY_NODE returns the next boundary node in any element.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE, the index of the current node.  An input
        //    value of 0 (or any "unusual" value") indicates that the
        //    first edge node is desired.
        //
        //    Input, string CODE, identifies the element desired.
        //    Legal values include "Q4", "Q8", "Q9", "Q12", "Q16", "QL", 
        //    "T3", "T4", "T6" and "T10".
        //
        //    Output, int NEXT_BOUNDARY_NODE, the index of the next edge node.
        //
    {
        int value = code switch
        {
            "Q4" => next_boundary_node_q4(node),
            "Q8" => next_boundary_node_q8(node),
            "Q9" => next_boundary_node_q9(node),
            "Q12" => next_boundary_node_q12(node),
            "Q16" => next_boundary_node_q16(node),
            "QL" => next_boundary_node_ql(node),
            "T3" => next_boundary_node_t3(node),
            "T4" => next_boundary_node_t4(node),
            "T6" => next_boundary_node_t6(node),
            "T10" => next_boundary_node_t10(node),
            _ => -1
        };

        return value;
    }

    public static int next_boundary_node_q4(int node)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    NEXT_BOUNDARY_NODE_Q4 returns the next boundary node in a Q4 element.
        //
        //  Element Q4:
        //
        //    |
        //    1  4-----3
        //    |  |     |
        //    |  |     |
        //    S  |     |
        //    |  |     |
        //    |  |     |
        //    0  1-----2
        //    |
        //    +--0--R--1-->
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE, the index of the current node.  An input
        //    value of 0 (or any "unusual" value") indicates that the
        //    first edge node is desired.
        //
        //    Output, int NEXT_BOUNDARY_NODE_Q4, the index of the next edge node.
        //
    {
        int value = node switch
        {
            1 => 2,
            2 => 3,
            3 => 4,
            4 => 1,
            _ => -1
        };

        return value;
    }

    public static int next_boundary_node_q8(int node)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    NEXT_BOUNDARY_NODE_Q8 returns the next boundary node in a Q8 element.
        //
        //  Element Q8:
        //
        //    |
        //    1  4--7--3
        //    |  |     |
        //    |  |     |
        //    S  8     6
        //    |  |     |
        //    |  |     |
        //    0  1--5--2
        //    |
        //    +--0--R--1-->
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE, the index of the current node.  An input
        //    value of 0 (or any "unusual" value") indicates that the
        //    first edge node is desired.
        //
        //    Output, int NEXT_BOUNDARY_NODE_Q8, the index of the next edge node.
        //
    {
        int value = node switch
        {
            1 => 5,
            5 => 2,
            2 => 6,
            6 => 3,
            3 => 7,
            7 => 4,
            4 => 8,
            8 => 1,
            _ => -1
        };

        return value;
    }

    public static int next_boundary_node_q9(int node)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    NEXT_BOUNDARY_NODE_Q9 returns the next boundary node in a Q9 element.
        //
        //  Element Q9:
        //
        //    |
        //    1  4--7--3
        //    |  |     |
        //    |  |     |
        //    S  8  9  6
        //    |  |     |
        //    |  |     |
        //    0  1--5--2
        //    |
        //    +--0--R--1-->
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE, the index of the current node.  An input
        //    value of 0 (or any "unusual" value") indicates that the
        //    first edge node is desired.
        //
        //    Output, int NEXT_BOUNDARY_NODE_Q9, the index of the next edge node.
        //
    {
        int value = node switch
        {
            1 => 5,
            5 => 2,
            2 => 6,
            6 => 3,
            3 => 7,
            7 => 4,
            4 => 8,
            _ => 1
        };

        return value;
    }

    public static int next_boundary_node_q12(int node)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    NEXT_BOUNDARY_NODE_Q12 returns the next boundary node in a Q12 element.
        //
        //  Element Q12:
        //
        //    |
        //    1  9-10-11-12
        //    |  |        |
        //    |  7        8
        //    S  |        |
        //    |  5        6
        //    |  |        |
        //    0  1--2--3--4
        //    |
        //    +--0---R---1-->
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE, the index of the current node.  An input
        //    value of 0 (or any "unusual" value") indicates that the
        //    first edge node is desired.
        //
        //    Output, int NEXT_BOUNDARY_NODE_Q12, the index of the next edge node.
        //
    {
        int value = node switch
        {
            1 => 2,
            2 => 3,
            3 => 4,
            4 => 6,
            6 => 8,
            8 => 12,
            12 => 11,
            11 => 10,
            10 => 9,
            9 => 7,
            7 => 5,
            5 => 1,
            _ => -1
        };

        return value;
    }

    public static int next_boundary_node_q16(int node)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    NEXT_BOUNDARY_NODE_Q16 returns the next boundary node in a Q16 element.
        //
        //  Element Q16:
        //
        //    |
        //    1 13--14--15--16
        //    |  |   :   :   |
        //    |  |   :   :   |
        //    |  9..10..11..12
        //    S  |   :   :   |
        //    |  |   :   :   |
        //    |  5...6...7...8
        //    |  |   :   :   |
        //    |  |   :   :   |  
        //    0  1---2---3---4
        //    |
        //    +--0-----R-----1-->
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE, the index of the current node.  An input
        //    value of 0 (or any "unusual" value") indicates that the
        //    first edge node is desired.
        //
        //    Output, int NEXT_BOUNDARY_NODE_Q16, the index of the next edge node.
        //
    {
        int value = node switch
        {
            1 => 2,
            2 => 3,
            3 => 4,
            4 => 8,
            8 => 12,
            12 => 16,
            16 => 15,
            15 => 14,
            14 => 13,
            13 => 9,
            9 => 5,
            5 => 1,
            _ => -1
        };

        return value;
    }

    public static int next_boundary_node_ql(int node)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    NEXT_BOUNDARY_NODE_QL returns the next boundary node in a QL element.
        //
        //  Element QL:
        //
        //    |
        //    1  4---5---6
        //    |  |       |
        //    |  |       |
        //    S  |       |
        //    |  |       |
        //    |  |       |
        //    0  1---2---3
        //    |
        //    +--0---R---1-->
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE, the index of the current node.  An input
        //    value of 0 (or any "unusual" value") indicates that the
        //    first edge node is desired.
        //
        //    Output, int NEXT_BOUNDARY_NODE_QL, the index of the next edge node.
        //
    {
        int value = node switch
        {
            1 => 2,
            2 => 3,
            3 => 6,
            6 => 5,
            5 => 4,
            4 => 1,
            _ => -1
        };

        return value;
    }

    public static int next_boundary_node_t3(int node)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    NEXT_BOUNDARY_NODE_T3 returns the next boundary node in a T3 element.
        //
        //  Element T3:
        //
        //    |
        //    1  3
        //    |  ..
        //    |  . .
        //    S  .  .
        //    |  .   .
        //    |  .    .
        //    0  1-----2
        //    |
        //    +--0--R--1-->
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE, the index of the current node.  An input
        //    value of 0 (or any "unusual" value") indicates that the
        //    first edge node is desired.
        //
        //    Output, int NEXT_BOUNDARY_NODE_T3, the index of the next edge node.
        //
    {
        int value = node switch
        {
            1 => 2,
            2 => 3,
            3 => 1,
            _ => -1
        };

        return value;
    }

    public static int next_boundary_node_t4(int node)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    NEXT_BOUNDARY_NODE_T4 returns the next boundary node in a T4 element.
        //
        //  Element T4:
        //
        //    |
        //    1  3
        //    |  ..
        //    |  . .
        //    S  .  .
        //    |  . 4 .
        //    |  .    .
        //    0  1-----2
        //    |
        //    +--0--R--1-->
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 March 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE, the index of the current node.  An input
        //    value of 0 (or any "unusual" value") indicates that the
        //    first edge node is desired.
        //
        //    Output, int NEXT_BOUNDARY_NODE_T4, the index of the next edge node.
        //
    {
        int value = node switch
        {
            1 => 2,
            2 => 3,
            3 => 1,
            _ => -1
        };

        return value;
    }

    public static int next_boundary_node_t6(int node)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    NEXT_BOUNDARY_NODE_T6 returns the next boundary node in a T6 element.
        //
        //  Element T6:
        //
        //    |
        //    1  3
        //    |  ..
        //    |  . .
        //    S  6  5
        //    |  .   .
        //    |  .    .
        //    0  1--4--2
        //    |
        //    +--0--R--1-->
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE, the index of the current node.  An input
        //    value of 0 (or any "unusual" value") indicates that the
        //    first edge node is desired.
        //
        //    Output, int NEXT_BOUNDARY_NODE_T6, the index of the next edge node.
        //
    {
        int value = node switch
        {
            1 => 4,
            4 => 2,
            2 => 5,
            5 => 3,
            3 => 6,
            6 => 1,
            _ => -1
        };

        return value;
    }

    public static int next_boundary_node_t10(int node)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    NEXT_BOUNDARY_NODE_T10 returns the next boundary node in a T10 element.
        //
        //  Element T10:
        //
        //    |
        //    1  10
        //    |  ..
        //    |  . .
        //    |  8  9
        //    |  .   .
        //    S  .    .
        //    |  5  6  7
        //    |  .      .
        //    |  .       .
        //    0  1--2--3--4
        //    |
        //    +--0----R---1-->
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE, the index of the current node.  An input
        //    value of 0 (or any "unusual" value") indicates that the
        //    first edge node is desired.
        //
        //    Output, int NEXT_BOUNDARY_NODE_T10, the index of the next edge node.
        //
    {
        int value = node switch
        {
            1 => 2,
            2 => 3,
            3 => 4,
            4 => 7,
            7 => 9,
            9 => 10,
            10 => 8,
            8 => 5,
            _ => 1
        };

        return value;
    }

}