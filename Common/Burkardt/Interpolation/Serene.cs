﻿namespace Burkardt.Interpolation
{
    public class Serene
    {
        public static double serene(string type, double ve, double vn, double vne, double vnw,
                double vs, double vse, double vsw, double vw)

            //****************************************************************************80
            //
            //  Purpose: 
            //
            //    SERENE interpolates data using a Q8 element.
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
            //    Input, string TYPE, tells SERENE the geometry of the
            //    finite element that surrounds the point of interest.  The options
            //    are displayed in the following table, which suggests the meaning
            //    of each option by its position:
            //
            //        |   |
            //     NW * N * NE
            //        |   |
            //     -*-*-*-*-*-
            //        |   |
            //      W * C * E
            //        |   |
            //     -*-*-*-*-*-
            //        |   |
            //     SW * S * SE
            //        |   |
            //
            //    Input, double VE, VN, VNE, VNW, VS, VSE, VSW, VW,
            //    are the values of the function at the nodes to the east,
            //    north, northeast, northwest, south, southeast, southwest and
            //    west of the point of interest.  If the finite element is of
            //    type 'C', then all 8 values are needed.  However, if the
            //    finite element is of type 'SE', for instance, then only three
            //    values are needed, namely VE, VN, and VNW, since these are
            //    the only node positions defined in such a finite element.
            //
            //    Output, double SERENE, the interpolated value of the
            //    function at the point of interest.
            //
        {
            double eta;
            double pe;
            double pn;
            double pne;
            double pnw;
            double ps;
            double pse;
            double psw;
            double pw;
            double vterp = 0;
            double xsi;
            //
            //  To make this routine more general, simply pass in the values of XSI
            //  and ETA at which the interpolated value is desired.
            //
            //  By setting XSI = ETA = 0, we are asking for the interpolated value
            //  at the center of the finite element.
            //
            xsi = 0.0;
            eta = 0.0;
            //
            //  8 node center
            //
            //  Polynomial space is spanned by:
            //
            //         1
            //       x    y
            //    x^2  xy  y^2
            //      x^2y xy^2
            //
            //
            //    ^   1    4--7--3
            //    |        !     !
            //    E        !     !
            //    T   0    8  X  6
            //    A        !     !
            //    |        !     !
            //    V  -1    1--5--2
            //
            //            -1  0  1
            //
            //           <---XSI--->
            //
            if (type == "C")
            {
                psw = -0.25 * (1.0 - xsi) * (1.0 - eta)
                      * (1.0 + xsi + eta);
                pse = -0.25 * (1.0 + xsi) * (1.0 - eta)
                      * (1.0 - xsi + eta);
                pne = -0.25 * (1.0 + xsi) * (1.0 + eta)
                      * (1.0 - xsi - eta);
                pnw = -0.25 * (1.0 - xsi) * (1.0 + eta)
                      * (1.0 + xsi - eta);
                ps = 0.50 * (1.0 - xsi) * (1.0 + xsi)
                     * (1.0 - eta);
                pe = 0.50 * (1.0 + xsi) * (1.0 + eta)
                     * (1.0 - eta);
                pn = 0.50 * (1.0 - xsi) * (1.0 + xsi)
                     * (1.0 + eta);
                pw = 0.50 * (1.0 - xsi) * (1.0 + eta)
                     * (1.0 - eta);

                vterp = vsw * psw + vse * pse + vne * pne + vnw * pnw
                        + vs * ps + ve * pe + vn * pn + vw * pw;
            }
            //
            //  5 node side
            //
            //    ^   1
            //    |
            //    E
            //    T   0    8  X  6
            //    A        !     !
            //    |        !     !
            //    V  -1    1--5--2
            //
            //            -1  0  1
            //
            //           <---XSI--->
            //
            else if (type == "N")
            {
                psw = 0.5 * (xsi - 1.0) * (1.0 + xsi + eta);
                pse = -0.5 * (xsi + 1.0) * (1.0 - xsi + eta);
                ps = -(xsi + 1.0) * (xsi - 1.0);
                pe = 0.5 * (xsi + 1.0) * (eta + 1.0);
                pw = -0.5 * (xsi - 1.0) * (eta + 1.0);

                vterp = vsw * psw + vse * pse + vs * ps + ve * pe + vw * pw;
            }
            //
            //    ^   1    4--7
            //    |        !
            //    E        !
            //    T   0    8  X
            //    A        !
            //    |        !
            //    V  -1    1--5
            //
            //            -1  0  1
            //
            //           <---XSI--->
            //
            else if (type == "E")
            {
                pse = 0.5 * (eta - 1.0) * (1.0 + xsi + eta);
                pne = -0.5 * (eta + 1.0) * (1.0 + xsi - eta);
                ps = -0.5 * (xsi + 1.0) * (eta - 1.0);
                pn = 0.5 * (xsi + 1.0) * (eta + 1.0);
                pw = -(eta + 1.0) * (eta - 1.0);

                vterp = vse * pse + vne * pne + vs * ps + vn * pn + vw * pw;
            }
            //
            //  5 node side
            //
            //    ^   1       7--3
            //    |              !
            //    E              !
            //    T   0       X  6
            //    A              !
            //    |              !
            //    V  -1       5--2
            //
            //            -1  0  1
            //
            //           <---XSI--->
            //
            else if (type == "W")
            {
                pse = 0.5 * (eta - 1.0) * (1.0 - xsi + eta);
                pne = -0.5 * (eta + 1.0) * (1.0 - xsi - eta);
                ps = 0.5 * (xsi - 1.0) * (eta - 1.0);
                pe = -(eta - 1.0) * (eta + 1.0);
                pn = -0.5 * (xsi - 1.0) * (eta + 1.0);

                vterp = vse * pse + vne * pne + vs * ps + ve * pe + vn * pn;
            }
            //
            //  5 node side
            //
            //    ^   1    4--7--3
            //    |        !     !
            //    E        !     !
            //    T   0    8  X  6
            //    A
            //    |
            //    V  -1
            //
            //            -1  0  1
            //
            //           <---XSI--->
            //
            else if (type == "S")
            {
                pne = -0.5 * (xsi + 1.0) * (1.0 - xsi - eta);
                pnw = 0.5 * (xsi - 1.0) * (1.0 + xsi - eta);
                pe = -0.5 * (eta - 1.0) * (xsi + 1.0);
                pn = -(xsi + 1.0) * (xsi - 1.0);
                pw = 0.5 * (eta - 1.0) * (xsi - 1.0);

                vterp = vne * pne + vnw * pnw + ve * pe + vn * pn + vw * pw;
            }
            //
            //  3 node corner
            //
            //  Polynomial space is spanned by:
            //
            //         1
            //       x    y
            //
            //
            //    ^   1
            //    |
            //    E
            //    T   0    8  X
            //    A        !
            //    |        !
            //    V  -1    1--5
            //
            //            -1  0  1
            //
            //           <---XSI--->
            //
            else if (type == "NE")
            {
                psw = -1.0 - xsi - eta;
                ps = 1.0 + xsi;
                pw = 1.0 + eta;

                vterp = vsw * psw + vs * ps + vw * pw;
            }
            //
            //  3 node corner
            //
            //  Polynomial space is spanned by:
            //
            //         1
            //       x    y
            //
            //    ^   1
            //    |
            //    E
            //    T   0       X  6
            //    A              !
            //    |              !
            //    V  -1       5--2
            //
            //            -1  0  1
            //
            //           <---XSI--->
            //
            else if (type == "NW")
            {
                pse = 1.0 + xsi - eta;
                ps = 1.0 - xsi;
                pe = 1.0 + eta;

                vterp = vse * pse + vs * ps + ve * pe;
            }
            //
            //  3 node corner
            //
            //  Polynomial space is spanned by:
            //         1
            //       x    y
            //
            //
            //    ^   1    4--7
            //    |        !
            //    E        !
            //    T   0    8  X
            //    A
            //    |
            //    V  -1
            //
            //            -1  0  1
            //
            //           <---XSI--->
            //
            else if (type == "SE")
            {
                pnw = -1.0 - xsi + eta;
                pn = 1.0 + xsi;
                pw = 1.0 - eta;

                vterp = vnw * pnw + vn * pn + vw * pw;
            }
            //
            //  3 node corner
            //
            //  Polynomial space is spanned by:
            //
            //         1
            //       x    y
            //
            //    ^   1       7--3
            //    |              !
            //    E              !
            //    T   0       X  6
            //    A
            //    |
            //    V  -1
            //
            //            -1  0  1
            //
            //           <---XSI--->
            //
            else if (type == "SW")
            {
                pne = -1.0 + xsi + eta;
                pe = 1.0 - eta;
                pn = 1.0 - xsi;

                vterp = vne * pne + ve * pe + vn * pn;
            }

            return vterp;
        }

    }
}