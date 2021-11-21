using System;

namespace Burkardt.FEM;

public static class NodeReference
{

    public static void node_reference(string code, ref double[] r, ref double[] s, ref double area )

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    NODE_REFERENCE returns the basis nodes for any available element.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 April 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string CODE, identifies the element desired.
        //    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
        //    'T3', 'T4', 'T6' and 'T10'.
        //
        //    Output, double R[N], S[N], the coordinates of the basis nodes.
        //
        //    Output, ref double AREA, the area of the element.
        //
    {
        switch (code)
        {
            case "Q4":
                node_reference_q4(ref r, ref s, ref area);
                break;
            case "Q8":
                node_reference_q8(ref r, ref s, ref area);
                break;
            case "Q9":
                node_reference_q9(ref r, ref s, ref area);
                break;
            case "Q12":
                node_reference_q12(ref r, ref s, ref area);
                break;
            case "Q16":
                node_reference_q16(ref r, ref s, ref area);
                break;
            case "QL":
                node_reference_ql(ref r, ref s, ref area);
                break;
            case "T3":
                node_reference_t3(ref r, ref s, ref area);
                break;
            case "T4":
                node_reference_t4(ref r, ref s, ref area);
                break;
            case "T6":
                node_reference_t6(ref r, ref s, ref area);
                break;
            case "T10":
                node_reference_t10(ref r, ref s, ref area);
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("NODE_REFERENCE - Fatal error!");
                Console.WriteLine("  Illegal value of CODE = " + code + "");
                break;
        }
    }

    public static void node_reference_q4(ref double[] r, ref double[] s, ref double area )

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    NODE_REFERENCE_Q4 returns the basis nodes for a 4 node quadrilateral.
        //
        //  Element Q4:
        //
        //    |
        //    1  4-------3
        //    |  |       |
        //    |  |       |
        //    S  |       |
        //    |  |       |
        //    |  |       |
        //    0  1-------2
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
        //    Output, double R[4], S[4], the coordinates of the basis nodes.
        //
        //    Output, ref double AREA, the area of the element.
        //
    {
        r[0] = 0.0;
        r[1] = 1.0;
        r[2] = 1.0;
        r[3] = 0.0;

        s[0] = 0.0;
        s[1] = 0.0;
        s[2] = 1.0;
        s[3] = 1.0;

        area = 1.0;
    }

    public static void node_reference_q8(ref double[] r, ref double[] s, ref double area )

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    NODE_REFERENCE_Q8 returns the basis nodes for an 8 node quadrilateral.
        //
        //  Comment:
        //
        //    This element is known as the quadratic "serendipity" element.
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
        //    Output, double R[8], S[8], the coordinates of the basis nodes.
        //
        //    Output, ref double AREA, the area of the element.
        //
    {
        r[0] = 0.0;
        r[1] = 1.0;
        r[2] = 1.0;
        r[3] = 0.0;
        r[4] = 0.5;
        r[5] = 1.0;
        r[6] = 0.5;
        r[7] = 0.0;

        s[0] = 0.0;
        s[1] = 0.0;
        s[2] = 1.0;
        s[3] = 1.0;
        s[4] = 0.0;
        s[5] = 0.5;
        s[6] = 1.0;
        s[7] = 0.5;

        area = 1.0;
    }

    public static void node_reference_q9(ref double[] r, ref double[] s, ref double area )

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    NODE_REFERENCE_Q9 returns the basis nodes for a 9 node quadrilateral.
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
        //    01 April 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R[9], S[9], the coordinates of the basis nodes.
        //
        //    Output, ref double AREA, the area of the element.
        //
    {
        r[0] = 0.0;
        r[1] = 1.0;
        r[2] = 1.0;
        r[3] = 0.0;
        r[4] = 0.5;
        r[5] = 1.0;
        r[6] = 0.5;
        r[7] = 0.0;
        r[8] = 0.5;

        s[0] = 0.0;
        s[1] = 0.0;
        s[2] = 1.0;
        s[3] = 1.0;
        s[4] = 0.0;
        s[5] = 0.5;
        s[6] = 1.0;
        s[7] = 0.5;
        s[8] = 0.5;

        area = 1.0;
    }

    public static void node_reference_q12(ref double[] r, ref double[] s, ref double area )

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    NODE_REFERENCE_Q12 returns the basis nodes for a 12 node quadrilateral.
        //
        //  Discussion:
        //
        //    This element is known as the cubic "serendipity" element.
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
        //    Output, double R[12], S[12], the coordinates of the basis nodes.
        //
        //    Output, ref double AREA, the area of the element.
        //
    {
        double a = 0.0;
        double b = 1.0 / 3.0;
        double c = 2.0 / 3.0;
        double d = 1.0;

        r[0] = a;
        r[1] = b;
        r[2] = c;
        r[3] = d;
        r[4] = a;
        r[5] = d;
        r[6] = a;
        r[7] = d;
        r[8] = a;
        r[9] = b;
        r[10] = c;
        r[11] = d;

        s[0] = a;
        s[1] = a;
        s[2] = a;
        s[3] = a;
        s[4] = b;
        s[5] = b;
        s[6] = c;
        s[7] = c;
        s[8] = d;
        s[9] = d;
        s[10] = d;
        s[11] = d;

        area = 1.0;
    }

    public static void node_reference_q16(ref double[] r, ref double[] s, ref double area )

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    NODE_REFERENCE_Q16 returns the basis nodes for a 16 node quadrilateral.
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
        //    Output, double R[16], S[16], the coordinates of the basis nodes.
        //
        //    Output, ref double AREA, the area of the element.
        //
    {
        int i;

        int k = 0;
        for (i = 0; i <= 3; i++)
        {
            int j;
            for (j = 0; j <= 3; j++)
            {
                r[k] = j / 3.0;
                s[k] = i / 3.0;
                k += 1;
            }
        }

        area = 1.0;
    }

    public static void node_reference_ql(ref double[] r, ref double[] s, ref double area )

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    NODE_REFERENCE_QL returns the basis nodes for a quadratic/linear.
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
        //    Output, double R[6], S[6], the coordinates of the basis nodes.
        //
        //    Output, ref double AREA, the area of the element.
        //
    {
        r[0] = 0.0;
        r[1] = 0.5;
        r[2] = 1.0;
        r[3] = 0.0;
        r[4] = 0.5;
        r[5] = 1.0;

        s[0] = 0.0;
        s[1] = 0.0;
        s[2] = 0.0;
        s[3] = 1.0;
        s[4] = 1.0;
        s[5] = 1.0;

        area = 1.0;
    }

    public static void node_reference_t3(ref double[] r, ref double[] s, ref double area )

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    NODE_REFERENCE_T3 returns the basis nodes for the 3 node triangle.
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
        //    Output, double R[3], S[3], the coordinates of the basis nodes.
        //
        //    Output, ref double AREA, the area of the element.
        //
    {
        r[0] = 0.0;
        r[1] = 1.0;
        r[2] = 0.0;

        s[0] = 0.0;
        s[1] = 0.0;
        s[2] = 1.0;

        area = 0.5;
    }

    public static void node_reference_t4(ref double[] r, ref double[] s, ref double area )

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    NODE_REFERENCE_T4 returns the basis nodes for the 4 node triangle.
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
        //    Output, double R[4], S[4], the coordinates of the basis nodes.
        //
        //    Output, ref double AREA, the area of the element.
        //
    {
        r[0] = 0.0;
        r[1] = 1.0;
        r[2] = 0.0;
        r[3] = 1.0 / 3.0;

        s[0] = 0.0;
        s[1] = 0.0;
        s[2] = 1.0;
        s[3] = 1.0 / 3.0;

        area = 0.5;
    }

    public static void node_reference_t6(ref double[] r, ref double[] s, ref double area )

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    NODE_REFERENCE_T6 returns the basis nodes for a 6 node triangle.
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
        //    Output, double R[6], S[6], the coordinates of the basis nodes.
        //
        //    Output, ref double AREA, the area of the element.
        //
    {
        r[0] = 0.0;
        r[1] = 1.0;
        r[2] = 0.0;
        r[3] = 0.5;
        r[4] = 0.5;
        r[5] = 0.0;

        s[0] = 0.0;
        s[1] = 0.0;
        s[2] = 1.0;
        s[3] = 0.0;
        s[4] = 0.5;
        s[5] = 0.5;

        area = 0.5;
    }

    public static void node_reference_t10(ref double[] r, ref double[] s, ref double area )

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    NODE_REFERENCE_T10 returns the basis nodes for a 10 node triangle.
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
        //    Output, double R[10], S[10], the coordinates of the basis nodes.
        //
        //    Output, ref double AREA, the area of the element.
        //
    {
        r[0] = 0.0;
        s[0] = 0.0;

        r[1] = 1.0 / 3.0;
        s[1] = 0.0;

        r[2] = 2.0 / 3.0;
        s[2] = 0.0;

        r[3] = 1.0;
        s[3] = 0.0;

        r[4] = 0.0;
        s[4] = 1.0 / 3.0;

        r[5] = 1.0 / 3.0;
        s[5] = 1.0 / 3.0;

        r[6] = 2.0 / 3.0;
        s[6] = 1.0 / 3.0;

        r[7] = 0.0;
        s[7] = 2.0 / 3.0;

        r[8] = 1.0 / 3.0;
        s[8] = 2.0 / 3.0;

        r[9] = 0.0;
        s[9] = 1.0;

        area = 0.5;
    }

}