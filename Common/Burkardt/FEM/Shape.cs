using System;
using System.Globalization;
using Burkardt.OrderNS;
using Burkardt.Types;

namespace Burkardt.FEM;

public class Shape
{
    public static void shape(string code, double r, double s, ref double[] t,
            ref double[] dtdr, ref double[] dtds)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    SHAPE evaluates shape functions for any available element.
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
        //    Input, string CODE, identifies the element.
        //    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
        //    'T3', 'T4', 'T6' and 'T10'.
        //
        //    Input, double R, S, the reference coordinates of a point.
        //
        //    Output, double T[N], the basis functions at the point.
        //
        //    Output, double DTDR[N], the R basis derivatives at the point.
        //
        //    Output, double DTDS[N], the S basis derivatives at the point.
        //
    {
        if (typeMethods.s_eqi(code, "Q4"))
        {
            shape_q4(r, s, t, ref dtdr, ref dtds);
        }
        else if (typeMethods.s_eqi(code, "Q8"))
        {
            shape_q8(r, s, t, ref dtdr, ref dtds);
        }
        else if (typeMethods.s_eqi(code, "Q9"))
        {
            shape_q9(r, s, t, ref dtdr, ref dtds);
        }
        else if (typeMethods.s_eqi(code, "Q12"))
        {
            shape_q12(r, s, t, ref dtdr, ref dtds);
        }
        else if (typeMethods.s_eqi(code, "Q16"))
        {
            shape_q16(r, s, t, ref dtdr, ref dtds);
        }
        else if (typeMethods.s_eqi(code, "QL"))
        {
            shape_ql(r, s, t, ref dtdr, ref dtds);
        }
        else if (typeMethods.s_eqi(code, "T3"))
        {
            shape_t3(r, s, t, ref dtdr, ref dtds);
        }
        else if (typeMethods.s_eqi(code, "T4"))
        {
            shape_t4(r, s, t, ref dtdr, ref dtds);
        }
        else if (typeMethods.s_eqi(code, "T6"))
        {
            shape_t6(r, s, t, ref dtdr, ref dtds);
        }
        else if (typeMethods.s_eqi(code, "T10"))
        {
            shape_t10(r, s, t, ref dtdr, ref dtds);
        }
        else
        {
            Console.WriteLine("");
            Console.WriteLine("SHAPE - Fatal error!");
            Console.WriteLine("  Unrecognized code = " + code + "");
        }
    }

    public static void shape_q4(double r, double s, double[] t, ref double[] dtdr,
            ref double[] dtds)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    SHAPE_Q4 evaluates shape functions for a 4 node quadrilateral.
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
        //    Input, double R, S, the reference coordinates of a point.
        //
        //    Output, double T[4], the basis functions at the point.
        //
        //    Output, double DTDR[4], the R basis derivatives at the point.
        //
        //    Output, double DTDS[4], the S basis derivatives at the point.
        //
    {
        t[0] = (1.0 - r) * (1.0 - s);
        t[1] = r * (1.0 - s);
        t[2] = r * s;
        t[3] = (1.0 - r) * s;

        dtdr[0] = -1.0 + s;
        dtdr[1] = 1.0 - s;
        dtdr[2] = s;
        dtdr[3] = -s;

        dtds[0] = -1.0 + r;
        dtds[1] = -r;
        dtds[2] = r;
        dtds[3] = 1.0 - r;
    }

    public static void shape_q8(double r, double s, double[] t, ref double[] dtdr,
            ref double[] dtds)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    SHAPE_Q8 evaluates shape functions for an 8 node quadrilateral.
        //
        //  Comment:
        //
        //    This element is known as the "serendipity" element.
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
        //    Input, double R, S, the reference coordinates of a point.
        //
        //    Output, double T[8], the basis functions at the point.
        //
        //    Output, double DTDR[8], the R basis derivatives at the point.
        //
        //    Output, double DTDS[8], the S basis derivatives at the point.
        //
    {
        t[0] = (r - 1.0) * (s - 1.0)
                         * (1.0 - 2.0 * r - 2.0 * s);
        t[1] = r * (s - 1.0)
                 * (1.0 - 2.0 * r + 2.0 * s);
        t[2] = r * s
                 * (2.0 * r + 2.0 * s - 3.0);
        t[3] = (r - 1.0) * s
                         * (2.0 * r - 2.0 * s + 1.0);
        t[4] = 4.0 * r * (r - 1.0) * (s - 1.0);
        t[5] = -4.0 * r * s * (s - 1.0);
        t[6] = -4.0 * r * (r - 1.0) * s;
        t[7] = 4.0 * (r - 1.0) * s * (s - 1.0);

        dtdr[0] = (s - 1.0) * (-4.0 * r - 2.0 * s + 3.0);
        dtdr[1] = (s - 1.0) * (-4.0 * r + 2.0 * s + 1.0);
        dtdr[2] = s * (4.0 * r + 2.0 * s - 3.0);
        dtdr[3] = s * (4.0 * r - 2.0 * s - 1.0);
        dtdr[4] = 4.0 * (2.0 * r - 1.0) * (s - 1.0);
        dtdr[5] = -4.0 * s * (s - 1.0);
        dtdr[6] = -4.0 * (2.0 * r - 1.0) * s;
        dtdr[7] = 4.0 * s * (s - 1.0);

        dtds[0] = (r - 1.0) * (-4.0 * s - 2.0 * r + 3.0);
        dtds[1] = r * (4.0 * s - 2.0 * r - 1.0);
        dtds[2] = r * (4.0 * s + 2.0 * r - 3.0);
        dtds[3] = (r - 1.0) * (-4.0 * s + 2.0 * r + 1.0);
        dtds[4] = 4.0 * r * (r - 1.0);
        dtds[5] = -4.0 * r * (2.0 * s - 1.0);
        dtds[6] = -4.0 * r * (r - 1.0);
        dtds[7] = 4.0 * (r - 1.0) * (2.0 * s - 1.0);
    }

    public static void shape_q9(double r, double s, double[] t, ref double[] dtdr,
            ref double[] dtds)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    SHAPE_Q9 evaluates shape functions for a 9 node quadrilateral.
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
        //    Input, double R, S, the reference coordinates of a point.
        //
        //    Output, double T[9], the basis functions at the point.
        //
        //    Output, double DTDR[9], the R basis derivatives at the point.
        //
        //    Output, double DTDS[9], the S basis derivatives at the point.
        //
    {
        t[0] = 4.0 * (r - 1.0) * (r - 0.5) * (s - 1.0)
               * (s - 0.5);
        t[1] = 4.0 * r * (r - 0.5) * (s - 1.0) * (s - 0.5);
        t[2] = 4.0 * r * (r - 0.5) * s * (s - 0.5);
        t[3] = 4.0 * (r - 1.0) * (r - 0.5) * s * (s - 0.5);
        t[4] = -8.0 * r * (r - 1.0) * (s - 1.0) * (s - 0.5);
        t[5] = -8.0 * r * (r - 0.5) * s * (s - 1.0);
        t[6] = -8.0 * r * (r - 1.0) * s * (s - 0.5);
        t[7] = -8.0 * (r - 1.0) * (r - 0.5) * s * (s - 1.0);
        t[8] = 16.0 * r * (r - 1.0) * s * (s - 1.0);

        dtdr[0] = 4.0 * (2.0 * r - 1.5) * (s - 1.0)
                  * (s - 0.5);
        dtdr[1] = 4.0 * (2.0 * r - 0.5) * (s - 1.0)
                  * (s - 0.5);
        dtdr[2] = 4.0 * (2.0 * r - 0.5) * s * (s - 0.5);
        dtdr[3] = 4.0 * (2.0 * r - 1.5) * s * (s - 0.5);

        dtdr[4] = -8.0 * (2.0 * r - 1.0) * (s - 1.0)
                  * (s - 0.5);
        dtdr[5] = -8.0 * (2.0 * r - 0.5) * s * (s - 1.0);
        dtdr[6] = -8.0 * (2.0 * r - 1.0) * s * (s - 0.5);
        dtdr[7] = -8.0 * (2.0 * r - 1.5) * s * (s - 1.0);
        dtdr[8] = 16.0 * (2.0 * r - 1.0) * s * (s - 1.0);

        dtds[0] = 4.0 * (r - 1.0) * (r - 0.5)
                  * (2.0 * s - 1.5);
        dtds[1] = 4.0 * r * (r - 0.5) * (2.0 * s - 1.5);
        dtds[2] = 4.0 * r * (r - 0.5) * (2.0 * s - 0.5);
        dtds[3] = 4.0 * (r - 1.0) * (r - 0.5)
                  * (2.0 * s - 0.5);
        dtds[4] = -8.0 * r * (r - 1.0) * (2.0 * s - 1.5);
        dtds[5] = -8.0 * r * (r - 0.5) * (2.0 * s - 1.0);
        dtds[6] = -8.0 * r * (r - 1.0) * (2.0 * s - 0.5);
        dtds[7] = -8.0 * (r - 1.0) * (r - 0.5)
                  * (2.0 * s - 1.0);
        dtds[8] = 16.0 * r * (r - 1.0) * (2.0 * s - 1.0);
    }

    public static void shape_q12(double r, double s, double[] t, ref double[] dtdr,
            ref double[] dtds)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    SHAPE_Q12 evaluates shape functions for a 12 node quadrilateral.
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
        //    Input, double R, S, the reference coordinates of a point.
        //
        //    Output, double T[12], the basis functions at the point.
        //
        //    Output, double DTDR[12], the R basis derivatives at the point.
        //
        //    Output, double DTDS[12], the S basis derivatives at the point.
        //
    {
        double a;
        double b;
        double c;
        double corner;
        double d;
        double dcdr;
        double dcds;

        a = 0.0;
        b = 1.0 / 3.0;
        c = 2.0 / 3.0;
        d = 1.0;

        corner = 9.0 * ((2.0 * r - 1.0) * (2.0 * r - 1.0)
                        + (2.0 * s - 1.0) * (2.0 * s - 1.0))
                 - 10.0;

        t[0] = 0.125 * (r - d) * (s - d) * corner;
        t[1] = -13.5 * (r - a) * (r - c) * (r - d) * (s - d);
        t[2] = 13.5 * (r - a) * (r - b) * (r - d) * (s - d);
        t[3] = -0.125 * (r - a) * (s - d) * corner;
        t[4] = -13.5 * (r - d) * (s - a) * (s - c) * (s - d);
        t[5] = 13.5 * (r - a) * (s - a) * (s - c) * (s - d);
        t[6] = 13.5 * (r - d) * (s - a) * (s - b) * (s - d);
        t[7] = -13.5 * (r - a) * (s - a) * (s - b) * (s - d);
        t[8] = -0.125 * (r - d) * (s - a) * corner;
        t[9] = 13.5 * (r - a) * (r - c) * (r - d) * (s - a);
        t[10] = -13.5 * (r - a) * (r - b) * (r - d) * (s - a);
        t[11] = 0.125 * (r - a) * (s - a) * corner;

        dcdr = 36.0 * (2.0 * r - 1.0);

        dtdr[0] = 0.125 * (s - d) * ((r - d) * dcdr + corner);
        dtdr[1] = -13.5 * (s - d) * (3.0 * r * r
            - 2.0 * (a + c + d) * r + a * c + c * d + d * a);
        dtdr[2] = 13.5 * (s - d) * (3.0 * r * r
            - 2.0 * (a + b + d) * r + a * b + b * d + d * a);
        dtdr[3] = -0.125 * (s - d) * ((r - a) * dcdr + corner);
        dtdr[4] = -13.5 * (s - a) * (s - c) * (s - d);
        dtdr[5] = 13.5 * (s - a) * (s - c) * (s - d);
        dtdr[6] = 13.5 * (s - a) * (s - b) * (s - d);
        dtdr[7] = -13.5 * (s - a) * (s - b) * (s - d);
        dtdr[8] = -0.125 * (s - a) * ((r - d) * dcdr + corner);
        dtdr[9] = 13.5 * (s - a) * (3.0 * r * r
            - 2.0 * (a + c + d) * r + a * c + c * d + d * a);
        dtdr[10] = -13.5 * (s - a) * (3.0 * r * r
            - 2.0 * (a + b + d) * r + a * b + b * d + d * a);
        dtdr[11] = 0.125 * (s - a) * ((r - a) * dcdr + corner);

        dcds = 36.0 * (2.0 * s - 1.0);

        dtds[0] = 0.125 * (r - d) * (corner + (s - d) * dcds);
        dtds[1] = -13.5 * (r - a) * (r - c) * (r - d);
        dtds[2] = 13.5 * (r - a) * (r - b) * (r - d);
        dtds[3] = -0.125 * (r - a) * (corner + (s - d) * dcds);
        dtds[4] = -13.5 * (r - d) * (3.0 * s * s
            - 2.0 * (a + c + d) * s + a * c + c * d + d * a);
        dtds[5] = 13.5 * (r - a) * (3.0 * s * s
            - 2.0 * (a + c + d) * s + a * c + c * d + d * a);
        dtds[6] = 13.5 * (r - d) * (3.0 * s * s
            - 2.0 * (a + b + d) * s + a * b + b * d + d * a);
        dtds[7] = -13.5 * (r - a) * (3.0 * s * s
            - 2.0 * (a + b + d) * s + a * b + b * d + d * a);
        dtds[8] = -0.125 * (r - d) * (corner + (s - a) * dcds);
        dtds[9] = 13.5 * (r - a) * (r - c) * (r - d);
        dtds[10] = -13.5 * (r - a) * (r - b) * (r - d);
        dtds[11] = 0.125 * (r - a) * (corner + (s - a) * dcds);
    }

    public static void shape_q16(double r, double s, double[] t, ref double[] dtdr,
            ref double[] dtds)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    SHAPE_Q16 evaluates shape functions for a 16 node quadrilateral.
        //
        //  Diagram:
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
        //    Input, double R, S, the reference coordinates of a point.
        //
        //    Output, double T[16], the basis functions at the point.
        //
        //    Output, double DTDR[16], the R basis derivatives at the point.
        //
        //    Output, double DTDS[16], the S basis derivatives at the point.
        //
    {
        double dabc;
        double dabd;
        double dacd;
        double dbcd;
        double ra;
        double rb;
        double rc;
        double rd;
        double sa;
        double sb;
        double sc;
        double sd;

        ra = r - 0.0;
        rb = r - 1.0 / 3.0;
        rc = r - 2.0 / 3.0;
        rd = r - 1.0;

        sa = s - 0.0;
        sb = s - 1.0 / 3.0;
        sc = s - 2.0 / 3.0;
        sd = s - 1.0;

        t[0] = 81.0 / 4.0 * rb * rc * rd * sb * sc * sd;
        t[1] = -(243.0 / 4.0) * ra * rc * rd * sb * sc * sd;
        t[2] = 243.0 / 4.0 * ra * rb * rd * sb * sc * sd;
        t[3] = -(81.0 / 4.0) * ra * rb * rc * sb * sc * sd;

        t[4] = -(243.0 / 4.0) * rb * rc * rd * sa * sc * sd;
        t[5] = 729.0 / 4.0 * ra * rc * rd * sa * sc * sd;
        t[6] = -(729.0 / 4.0) * ra * rb * rd * sa * sc * sd;
        t[7] = 243.0 / 4.0 * ra * rb * rc * sa * sc * sd;

        t[8] = 243.0 / 4.0 * rb * rc * rd * sa * sb * sd;
        t[9] = -(729.0 / 4.0) * ra * rc * rd * sa * sb * sd;
        t[10] = 729.0 / 4.0 * ra * rb * rd * sa * sb * sd;
        t[11] = -(243.0 / 4.0) * ra * rb * rc * sa * sb * sd;

        t[12] = -(81.0 / 4.0) * rb * rc * rd * sa * sb * sc;
        t[13] = 243.0 / 4.0 * ra * rc * rd * sa * sb * sc;
        t[14] = -(243.0 / 4.0) * ra * rb * rd * sa * sb * sc;
        t[15] = 81.0 / 4.0 * ra * rb * rc * sa * sb * sc;

        dbcd = 3.0 * r * r - 4.0 * r + 11.0 / 9.0;
        dacd = 3.0 * r * r - 10.0 * r / 3.0 + 2.0 / 3.0;
        dabd = 3.0 * r * r - 8.0 * r / 3.0 + 1.0 / 3.0;
        dabc = 3.0 * r * r - 2.0 * r + 2.0 / 9.0;

        dtdr[0] = 81.0 / 4.0 * dbcd * sb * sc * sd;
        dtdr[1] = -(243.0 / 4.0) * dacd * sb * sc * sd;
        dtdr[2] = 243.0 / 4.0 * dabd * sb * sc * sd;
        dtdr[3] = -(81.0 / 4.0) * dabc * sb * sc * sd;
        dtdr[4] = -(243.0 / 4.0) * dbcd * sa * sc * sd;
        dtdr[5] = 729.0 / 4.0 * dacd * sa * sc * sd;
        dtdr[6] = -(729.0 / 4.0) * dabd * sa * sc * sd;
        dtdr[7] = 243.0 / 4.0 * dabc * sa * sc * sd;
        dtdr[8] = 243.0 / 4.0 * dbcd * sa * sb * sd;
        dtdr[9] = -(729.0 / 4.0) * dacd * sa * sb * sd;
        dtdr[10] = 729.0 / 4.0 * dabd * sa * sb * sd;
        dtdr[11] = -(243.0 / 4.0) * dabc * sa * sb * sd;
        dtdr[12] = -(81.0 / 4.0) * dbcd * sa * sb * sc;
        dtdr[13] = 243.0 / 4.0 * dacd * sa * sb * sc;
        dtdr[14] = -(243.0 / 4.0) * dabd * sa * sb * sc;
        dtdr[15] = 81.0 / 4.0 * dabc * sa * sb * sc;

        dbcd = 3.0 * s * s - 4.0 * s + 11.0 / 9.0;
        dacd = 3.0 * s * s - 10.0 * s / 3.0 + 2.0 / 3.0;
        dabd = 3.0 * s * s - 8.0 * s / 3.0 + 1.0 / 3.0;
        dabc = 3.0 * s * s - 2.0 * s + 2.0 / 9.0;

        dtds[0] = 81.0 / 4.0 * rb * rc * rd * dbcd;
        dtds[1] = -(243.0 / 4.0) * ra * rc * rd * dbcd;
        dtds[2] = 243.0 / 4.0 * ra * rb * rd * dbcd;
        dtds[3] = -(81.0 / 4.0) * ra * rb * rc * dbcd;
        dtds[4] = -(243.0 / 4.0) * rb * rc * rd * dacd;
        dtds[5] = 729.0 / 4.0 * ra * rc * rd * dacd;
        dtds[6] = -(729.0 / 4.0) * ra * rb * rd * dacd;
        dtds[7] = 243.0 / 4.0 * ra * rb * rc * dacd;
        dtds[8] = 243.0 / 4.0 * rb * rc * rd * dabd;
        dtds[9] = -(729.0 / 4.0) * ra * rc * rd * dabd;
        dtds[10] = 729.0 / 4.0 * ra * rb * rd * dabd;
        dtds[11] = -(243.0 / 4.0) * ra * rb * rc * dabd;
        dtds[12] = -(81.0 / 4.0) * rb * rc * rd * dabc;
        dtds[13] = 243.0 / 4.0 * ra * rc * rd * dabc;
        dtds[14] = -(243.0 / 4.0) * ra * rb * rd * dabc;
        dtds[15] = 81.0 / 4.0 * ra * rb * rc * dabc;
    }

    public static void shape_ql(double r, double s, double[] t, ref double[] dtdr,
            ref double[] dtds)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    SHAPE_QL evaluates shape functions for a 6 node quadratic/linear.
        //
        //  Diagram:
        //
        //    |
        //    1  4--5--6
        //    |  |     |
        //    |  |     |
        //    S  |     |
        //    |  |     |
        //    |  |     |
        //    0  1--2--3
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
        //    Input, double R, S, the reference coordinates of a point.
        //
        //    Output, double T[6], the basis functions at the point.
        //
        //    Output, double DTDR[6], the R basis derivatives at the point.
        //
        //    Output, double DTDS[6], the S basis derivatives at the point.
        //
    {
        t[0] = -2.0 * (r - 0.5) * (r - 1.0) * (s - 1.0);
        t[1] = 4.0 * r * (r - 1.0) * (s - 1.0);
        t[2] = -2.0 * r * (r - 0.5) * (s - 1.0);
        t[3] = 2.0 * (r - 0.5) * (r - 1.0) * s;
        t[4] = -4.0 * r * (r - 1.0) * s;
        t[5] = 2.0 * r * (r - 0.5) * s;

        dtdr[0] = 2.0 * (-2.0 * r + 1.5) * (s - 1.0);
        dtdr[1] = 4.0 * (2.0 * r - 1.0) * (s - 1.0);
        dtdr[2] = 2.0 * (-2.0 * r + 0.5) * (s - 1.0);
        dtdr[3] = 2.0 * (2.0 * r - 1.5) * s;
        dtdr[4] = 4.0 * (-2.0 * r + 1.0) * s;
        dtdr[5] = 2.0 * (2.0 * r - 0.5) * s;

        dtds[0] = -2.0 * (r - 0.5) * (r - 1.0);
        dtds[1] = 4.0 * r * (r - 1.0);
        dtds[2] = -2.0 * r * (r - 0.5);
        dtds[3] = 2.0 * (r - 0.5) * (r - 1.0);
        dtds[4] = -4.0 * r * (r - 1.0);
        dtds[5] = 2.0 * r * (r - 0.5);
    }

    public static void shape_t3(double r, double s, double[] t, ref double[] dtdr,
            ref double[] dtds)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    SHAPE_T3 evaluates shape functions for a 3 node triangle.
        //
        //  Element T3:
        //
        //    |
        //    1  3
        //    |  |.
        //    |  | .
        //    S  |  .
        //    |  |   .
        //    |  |    .
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
        //    Input, double R, S, the reference coordinates of a point.
        //
        //    Output, double T[3], the basis functions at the point.
        //
        //    Output, double DTDR[3], the R basis derivatives at the point.
        //
        //    Output, double DTDS[3], the S basis derivatives at the point.
        //
    {
        t[0] = 1.0 - r - s;
        t[1] = r;
        t[2] = s;

        dtdr[0] = -1.0;
        dtdr[1] = 1.0;
        dtdr[2] = 0.0;

        dtds[0] = -1.0;
        dtds[1] = 0.0;
        dtds[2] = 1.0;
    }

    public static void shape_t4(double r, double s, double[] t, ref double[] dtdr,
            ref double[] dtds)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    SHAPE_T4 evaluates shape functions for a T4 triangle.
        //
        //  Element T4:
        //
        //    |
        //    1  3
        //    |  |.
        //    |  | .
        //    S  |  .
        //    |  | 4 .
        //    |  |    .
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
        //    Input, double R, S, the reference coordinates of a point.
        //
        //    Output, double T[4], the basis functions at the point.
        //
        //    Output, double DTDR[3], the R basis derivatives at the point.
        //
        //    Output, double DTDS[3], the S basis derivatives at the point.
        //
    {
        t[0] = (1.0 - 9.0 * r * s) * (1.0 - r - s);
        t[1] = r * (1.0 - 9.0 * (1.0 - r - s) * s);
        t[2] = s * (1.0 - 9.0 * (1.0 - r - s) * r);
        t[3] = 27.0 * (1.0 - r - s) * r * s;

        dtdr[0] = -1.0 + 9.0 * (-s + 2.0 * r * s + s * s);
        dtdr[1] = 1.0 + 9.0 * (-s + 2.0 * r * s + s * s);
        dtdr[2] = 9.0 * (-s + 2.0 * r * s + s * s);
        dtdr[3] = -27.0 * (-s + 2.0 * r * s + s * s);

        dtds[0] = -1.0 + 9.0 * (-r + r * r + 2.0 * r * s);
        dtds[1] = 9.0 * (-r + r * r + 2.0 * r * s);
        dtds[2] = 1.0 + 9.0 * (-r + r * r + 2.0 * r * s);
        dtds[3] = -27.0 * (-r + r * r + 2.0 * r * s);
    }

    public static void shape_t6(double r, double s, double[] t, ref double[] dtdr,
            ref double[] dtds)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    SHAPE_T6 evaluates shape functions for a 6 node triangle.
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
        //    Input, double R, S, the reference coordinates of a point.
        //
        //    Output, double T[6], the basis functions at the point.
        //
        //    Output, double DTDR[6], the R basis derivatives at the point.
        //
        //    Output, double DTDS[6], the S basis derivatives at the point.
        //
    {
        t[0] = 2.0 * (1.0 - r - s) * (0.5 - r - s);
        t[1] = 2.0 * r * (r - 0.5);
        t[2] = 2.0 * s * (s - 0.5);
        t[3] = 4.0 * r * (1.0 - r - s);
        t[4] = 4.0 * r * s;
        t[5] = 4.0 * s * (1.0 - r - s);

        dtdr[0] = -3.0 + 4.0 * r + 4.0 * s;
        dtdr[1] = -1.0 + 4.0 * r;
        dtdr[2] = 0.0;
        dtdr[3] = 4.0 - 8.0 * r - 4.0 * s;
        dtdr[4] = 4.0 * s;
        dtdr[5] = -4.0 * s;

        dtds[0] = -3.0 + 4.0 * r + 4.0 * s;
        dtds[1] = 0.0;
        dtds[2] = -1.0 + 4.0 * s;
        dtds[3] = -4.0 * r;
        dtds[4] = 4.0 * r;
        dtds[5] = 4.0 - 4.0 * r - 8.0 * s;
    }

    public static void shape_t10(double r, double s, double[] t, ref double[] dtdr,
            ref double[] dtds)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    SHAPE_T10 evaluates shape functions for a 10 node triangle.
        //
        //  Diagram:
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
        //    Input, double R, S, the reference coordinates of a point.
        //
        //    Output, double T[10], the basis functions at the point.
        //
        //    Output, double DTDR[10], the R basis derivatives at the point.
        //
        //    Output, double DTDS[10], the S basis derivatives at the point.
        //
    {
        double a;
        double b;
        double c;

        a = 1.0 / 3.0;
        b = 2.0 / 3.0;
        c = 1.0;

        t[0] = 4.5 * (a - r - s) * (b - r - s) * (c - r - s);
        t[1] = 13.5 * r * (b - r - s) * (c - r - s);
        t[2] = -13.5 * r * (a - r) * (c - r - s);
        t[3] = 4.5 * r * (a - r) * (b - r);
        t[4] = 13.5 * s * (b - r - s) * (c - r - s);
        t[5] = 27.0 * r * s * (c - r - s);
        t[6] = 13.5 * r * s * (r - a);
        t[7] = 13.5 * s * (s - a) * (c - r - s);
        t[8] = 13.5 * r * s * (s - a);
        t[9] = 4.5 * s * (a - s) * (b - s);

        dtdr[0] = 4.5 * ((a - s) * (2.0 * r - c - b + 2.0 * s)
                         - (s - b) * (s - c) - 2.0 * (2.0 * s - b - c) * r
                         - 3.0 * r * r);
        dtdr[1] = 13.5 * (
            (s - b) * (s - c) + 2.0 * (2.0 * s - b - c) * r
                              + 3.0 * r * r);
        dtdr[2] = -13.5 * (a * (c - s) + 2.0 * (s - a - c) * r
                                       + 3.0 * r * r);
        dtdr[3] = 4.5 * (a * b - 2.0 * (a + b) * r + 3.0 * r * r);
        dtdr[4] = 13.5 * s * (2.0 * s - b - c + 2.0 * r);
        dtdr[5] = 27.0 * s * (c - s - 2.0 * r);
        dtdr[6] = 13.5 * s * (2.0 * r - a);
        dtdr[7] = -13.5 * s * (s - a);
        dtdr[8] = 13.5 * s * (s - a);
        dtdr[9] = 0.0;

        dtds[0] = 4.5 * ((a - r) * (2.0 * s - c - b + 2.0 * r)
                         - (r - b) * (r - c) - 2.0 * (2.0 * r - b - c) * s
                         - 3.0 * s * s);
        dtds[1] = 13.5 * r * (2.0 * s + 2.0 * r - b - c);
        dtds[2] = 13.5 * r * (a - r);
        dtds[3] = 0.0;
        dtds[4] = 13.5 * ((r - b) * (r - c) +
                          2.0 * (2.0 * r - b - c) * s + 3.0 * s * s);
        dtds[5] = 27.0 * r * (c - r - 2.0 * s);
        dtds[6] = 13.5 * r * (r - a);
        dtds[7] = -13.5 * (a * (c - r) + 2.0 * (r - c - a) * s
                                       + 3.0 * s * s);
        dtds[8] = 13.5 * r * (2.0 * s - a);
        dtds[9] = 4.5 * (a * b - 2.0 * (a + b) * s + 3.0 * s * s);
    }

    public static void shape_test(string code)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    SHAPE_TEST verifies the shape function values at the basis nodes.
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
        //    Input, string CODE, identifies the element to be used.
        //    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
        //    'T3', 'T6' and 'T10'.
        //
    {
        double area = 0;
        double[] dtdr;
        double[] dtds;
        int element_order;
        int i;
        int j;
        double[] r;
        double rsum;
        double[] s;
        double ssum;
        double[] t;

        Console.WriteLine("");
        Console.WriteLine("  SHAPE_TEST: Verify shape functions of type " + code + "");

        element_order = Order.order_code(code);

        dtdr = new double[element_order];
        dtds = new double[element_order];
        r = new double[element_order];
        s = new double[element_order];
        t = new double[element_order];

        NodeReference.node_reference(code, ref r, ref s, ref area);

        Console.WriteLine("");
        Console.WriteLine("  Element order = " + element_order + "");
        Console.WriteLine("  Basis function values at basis nodes");
        Console.WriteLine("  should form the identity matrix.");
        Console.WriteLine("");

        for (i = 0; i < element_order; i++)
        {
            shape(code, r[i], s[i], ref t, ref dtdr, ref dtds);
            string cout = "";
            for (j = 0; j < element_order; j++)
            {
                cout += "  " + t[j].ToString(CultureInfo.InvariantCulture).PadLeft(7);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  The R and S derivatives should sum to 0.");
        Console.WriteLine("");
        Console.WriteLine("  dTdR sum, dTdS sum:");
        Console.WriteLine("");
        for (i = 0; i < element_order; i++)
        {
            shape(code, r[i], s[i], ref t, ref dtdr, ref dtds);
            rsum = 0.0;
            for (j = 0; j < element_order; j++)
            {
                rsum += dtdr[j];
            }

            ssum = 0.0;
            for (j = 0; j < element_order; j++)
            {
                ssum += dtds[j];
            }

            Console.WriteLine("  " + rsum.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + ssum.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

}