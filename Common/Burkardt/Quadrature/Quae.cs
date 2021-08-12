using System;
using Burkardt.Types;

namespace Burkardt.Quadrature
{
    public static class Quae
    {
        public static void quaecopy(int nf, double[] xs, double[] ys, double[] ws, ref double[] z,
                ref double[] w)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    QUAECOPY copies a quadrature rule into user arrays Z and W.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU GPL license.
            //
            //  Modified:
            //
            //    30 June 2014
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
            //    This C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Hong Xiao, Zydrunas Gimbutas,
            //    A numerical algorithm for the construction of efficient quadrature
            //    rules in two and higher dimensions,
            //    Computers and Mathematics with Applications,
            //    Volume 59, 2010, pages 663-676.
            //
            //  Parameters:
            //
            //    Input, int NF, the number of values to copy.
            //
            //    Input, double XS[NF], YS[NF], the point coordinates to copy.
            //
            //    Input, double WS[NF], the weights to copy.
            //
            //    Output, double Z[2*NF], the copied point coordinates.
            //
            //    Output, double W[NF], the copied weights.
            //
        {
            int j;

            for (j = 0; j < nf; j++)
            {
                z[0 + j * 2] = xs[j];
                z[1 + j * 2] = ys[j];
            }

            typeMethods.r8vec_copy(nf, ws, ref w);
        }

        public static void quaecopy2(double[] xs, double[] ys, double[] ws, ref double[] xnew,
                ref double[] ynew, ref double[] w, int kk)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    QUAECOPY2 copies a quadrature rule into a user arrays X, Y, and W.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU GPL license.
            //
            //  Modified:
            //
            //    30 June 2014
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
            //    This C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Hong Xiao, Zydrunas Gimbutas,
            //    A numerical algorithm for the construction of efficient quadrature
            //    rules in two and higher dimensions,
            //    Computers and Mathematics with Applications,
            //    Volume 59, 2010, pages 663-676.
            //
            //  Parameters:
            //
            //    Input, double XS[KK], YS[KK], the point coordinates to copy.
            //
            //    Input, double WS[KK], the weights to copy.
            //
            //    Output, double XNEW[KK], YNEW[KK], the copied point coordinates.
            //
            //    Output, double W[KK], the copied weights.
            //
            //    Input, int KK, the number of values to copy.
            //
        {
            typeMethods.r8vec_copy(kk, xs, ref xnew);
            typeMethods.r8vec_copy(kk, ys, ref ynew);
            typeMethods.r8vec_copy(kk, ws, ref w);
        }

        public static int quaeinside(int iitype, double xsout, double ysout)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    QUAEINSIDE checks whether a point is inside a triangle.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU GPL license.
            //
            //  Modified:
            //
            //    30 June 2014
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
            //    This C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Hong Xiao, Zydrunas Gimbutas,
            //    A numerical algorithm for the construction of efficient quadrature
            //    rules in two and higher dimensions,
            //    Computers and Mathematics with Applications,
            //    Volume 59, 2010, pages 663-676.
            //
            //  Parameters:
            //
            //    Input, int IITYPE, indicates the check to perform.
            //    * 0 check whether the point is inside the whole triangle
            //    * 1 check whether the point is inside the bottom 1/3 of the triangle
            //    * 2 check whether the point is inside the lower-left 1/6 of the triangle
            //
            //    Input, double XSOUT, YSOUT, the coordinates of the point.
            //
            //    Output, int QUAEINSIDE.
            //    * 1, the point is inside.
            //    * 2, the point is outside.
            //
        {
            int nbool = 0;
            double s;

            s = Math.Sqrt(3.0);
            //
            //  The 1/6 of triangle.
            //
            if (iitype == 2)
            {
                nbool = 1;
                if (xsout < -1.0 || 0.0 < xsout)
                {
                    nbool = 0;
                }

                //
                //  hx
                //
                if (ysout < -1.0 / s - 1.0E-30 || xsout / s < ysout)
                {
                    nbool = 0;
                }
            }
            //
            //  The 1/3 of triangle.
            //
            else if (iitype == 1)
            {
                nbool = 0;

                if (xsout <= 0.0 && -1.0 <= xsout && -1.0 / s <= ysout && ysout <= xsout / s)
                {
                    nbool = 1;
                }

                if (iitype == 1)
                {
                    if (0.0 <= xsout && xsout <= 1.0 &&
                        -1.0 / s <= ysout && ysout <= -xsout / s)
                    {
                        nbool = 1;
                    }
                }
            }
            //
            //  The entire triangle.
            //
            else if (iitype == 0)
            {
                nbool = 1;
                if (ysout < -1.0 / s)
                {
                    nbool = 0;
                }

                if (s * xsout + 2.0 / s < ysout)
                {
                    nbool = 0;
                }

                if (-s * xsout + 2.0 / s < ysout)
                {
                    nbool = 0;
                }
            }

            return nbool;
        }

        public static void quaenodes(int nptsout, double[] xsout, double[] ysout, double[] wsout,
                ref int nptsoutout, ref double[] xs2, ref double[] ys2, ref double[] ws2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    QUAENODES expands nodes to the reference triangle.
            //
            //  Discussion:
            //
            //    This routine expands nodes to the reference triangle
            //    assuming that the points are already in the lower-left 1/6 of the
            //    triangle.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU GPL license.
            //
            //  Modified:
            //
            //    30 June 2014
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
            //    This C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Hong Xiao, Zydrunas Gimbutas,
            //    A numerical algorithm for the construction of efficient quadrature
            //    rules in two and higher dimensions,
            //    Computers and Mathematics with Applications,
            //    Volume 59, 2010, pages 663-676.
            //
            //  Parameters:
            //
            //    Input, int NPTSOUT, the number of points expanded.
            //
            //    Input, double XSOUT[NPTSOUT], YSOUT[NPTSOUT], WSOUT[NPTSOUT],
            //    the coordinates and weights of the points to be expanded.
            //
            //    Output, int NPTSOUTPUT, the number of points in the expanded set.
            //
            //    Output, double XS2[NPTSOUTOUT], YS2[NPTSOUTOUT], WS2[NPTSOUTOUT],
            //    the coordinates and weights of the expanded set of points.
            //
        {
            const double eps = 1.0E-12;
            int i;
            int ntot;
            double w0 = 0;
            double x0 = 0;
            double x1 = 0;
            double x2 = 0;
            double y0 = 0;
            double y1 = 0;
            double y2 = 0;

            ntot = 0;

            for (i = 0; i < nptsout; i++)
            {
                if (Math.Pow(xsout[i], 2) + Math.Pow(ysout[i], 2) < eps)
                {
                    xs2[ntot] = xsout[i];
                    ys2[ntot] = ysout[i];
                    ws2[ntot] = wsout[i];
                    ntot = ntot + 1;
                }
                else if (Math.Pow(xsout[i], 2) < eps ||
                         Math.Abs(ysout[i] - xsout[i] / Math.Sqrt(3.0)) < Math.Sqrt(eps))
                {
                    x0 = xsout[i];
                    y0 = ysout[i];
                    w0 = wsout[i] / 3.0;

                    xs2[ntot] = x0;
                    ys2[ntot] = y0;
                    ws2[ntot] = w0;
                    ntot = ntot + 1;

                    quaerotate(x0, y0, ref x1, ref y1);
                    xs2[ntot] = x1;
                    ys2[ntot] = y1;
                    ws2[ntot] = w0;
                    ntot = ntot + 1;

                    quaerotate(x1, y1, ref x2, ref y2);
                    xs2[ntot] = x2;
                    ys2[ntot] = y2;
                    ws2[ntot] = w0;
                    ntot = ntot + 1;
                }
                else
                {
                    x0 = xsout[i];
                    y0 = ysout[i];
                    w0 = wsout[i] / 6.0;

                    xs2[ntot] = x0;
                    ys2[ntot] = y0;
                    ws2[ntot] = w0;
                    ntot = ntot + 1;

                    quaerotate(x0, y0, ref x1, ref y1);
                    xs2[ntot] = x1;
                    ys2[ntot] = y1;
                    ws2[ntot] = w0;
                    ntot = ntot + 1;

                    quaerotate(x1, y1, ref x2, ref y2);
                    xs2[ntot] = x2;
                    ys2[ntot] = y2;
                    ws2[ntot] = w0;
                    ntot = ntot + 1;

                    xs2[ntot] = -x0;
                    ys2[ntot] = y0;
                    ws2[ntot] = w0;
                    ntot = ntot + 1;

                    quaerotate(-x0, y0, ref x1, ref y1);
                    xs2[ntot] = x1;
                    ys2[ntot] = y1;
                    ws2[ntot] = w0;
                    ntot = ntot + 1;

                    quaerotate(x1, y1, ref x2, ref y2);
                    xs2[ntot] = x2;
                    ys2[ntot] = y2;
                    ws2[ntot] = w0;
                    ntot = ntot + 1;
                }
            }

        }

        public static void quaenodes2(int nptsout, double[] xsout, double[] ysout, double[] wsout,
                ref int nptsoutout, ref double[] xs2, ref double[] ys2, ref double[] ws2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    QUAENODES2 expands nodes from 1/6 to 1/3 of the triangle.
            //
            //  Discussion:
            //
            //    This routine only expands to 1/3 of the triangle, assuming the points are
            //    already in the lower-left 1/6 of the triangle.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU GPL license.
            //
            //  Modified:
            //
            //    30 June 2014
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
            //    This C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Hong Xiao, Zydrunas Gimbutas,
            //    A numerical algorithm for the construction of efficient quadrature
            //    rules in two and higher dimensions,
            //    Computers and Mathematics with Applications,
            //    Volume 59, 2010, pages 663-676.
            //
            //  Parameters:
            //
            //    Input, int NPTSOUT, the number of points expanded.
            //
            //    Input, double XSOUT[NPTSOUT], YSOUT[NPTSOUT], WSOUT[NPTSOUT],
            //    the coordinates and weights of the points to be expanded.
            //
            //    Output, int NPTSOUTPUT, the number of points in the expanded set.
            //
            //    Output, double XS2[NPTSOUTOUT], YS2[NPTSOUTOUT], WS2[NPTSOUTOUT],
            //    the coordinates and weights of the expanded set of points.
            //
        {
            const double eps = 1.0E-12;
            int i;
            int ntot;

            ntot = 0;

            for (i = 0; i < nptsout; i++)
            {
                if (Math.Pow(xsout[i], 2) + Math.Pow(ysout[i], 2) < eps)
                {
                    xs2[ntot] = xsout[i];
                    ys2[ntot] = ysout[i];
                    ws2[ntot] = wsout[i];
                    ntot = ntot + 1;
                }
                else if (Math.Pow(xsout[i], 2) < eps)
                {
                    xs2[ntot] = xsout[i];
                    ys2[ntot] = ysout[i];
                    ws2[ntot] = wsout[i];
                    ntot = ntot + 1;
                }
                else
                {
                    xs2[ntot] = -xsout[i];
                    ys2[ntot] = ysout[i];
                    ws2[ntot] = wsout[i] / 2.0;
                    ntot = ntot + 1;

                    xs2[ntot] = xsout[i];
                    ys2[ntot] = ysout[i];
                    ws2[ntot] = wsout[i] / 2.0;
                    ntot = ntot + 1;
                }
            }

        }

        public static void quaequad(int itype, int mmax, ref double[] zs, ref double[] whts, int numnodes)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    QUAEQUAD returns a symmetric quadrature formula for a reference triangle.
            //
            //  Discussion:
            //
            //    This routine constructs (or rather, retrieves)
            //    D_3 symmetric quadrature formulae for smooth functions
            //    on the triangle with vertices
            //      (-1,-1/sqrt(3)), (1,-1/sqrt(3)), (0,2/sqrt(3)).
            //
            //    All quadratures are obtained to the extended precision.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU GPL license.
            //
            //  Modified:
            //
            //    30 June 2014
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
            //    This C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Hong Xiao, Zydrunas Gimbutas,
            //    A numerical algorithm for the construction of efficient quadrature
            //    rules in two and higher dimensions,
            //    Computers and Mathematics with Applications,
            //    Volume 59, 2010, pages 663-676.
            //
            //  Parameters:
            //
            //    Input, int ITYPE, the configuration of the quadrature nodes.
            //    * 0, all nodes on the reference triangle
            //    * 1, quadrature nodes on the bottom 1/3 of the reference triangle
            //    * 2, quadrature nodes on lower-left 1/6 of the reference triangle
            //
            //    Input, int MMAX, the degree of the quadrature (the maximum degree of
            //    the polynomials of two variables that are integrated
            //    exactly. 1 <= MMAX <= 50.
            //
            //    Output, double ZS[2*NUMNODES], the nodes.
            //
            //    Output, double WHTS[NUMNODES], the weights.
            //
            //    Input, int NUMNODES, the number of nodes in the quadrature rule.
            //
        {
            int nc;
            int nf;
            double[] wc;
            double[] wf;
            double[] xc;
            double[] xf;
            double[] yc;
            double[] yf;

            if (mmax < 1 || 50 < mmax)
            {
                Console.WriteLine("");
                Console.WriteLine("QUAEQUAD - Fatal error!");
                Console.WriteLine("  1 <= MMAX <= 50 is required.");
                return;
            }

            if (itype < 0 || 2 < itype)
            {
                Console.WriteLine("");
                Console.WriteLine("QUAEQUAD - Fatal error!");
                Console.WriteLine("  0 <= ITYPE <= 2 is required.");
                return;
            }

            //
            //  Get the size of the compressed rule.
            //
            nc = QuadratureRule.rule_compressed_size(mmax);
            //
            //  Retrieve the compressed rule.
            //
            xc = new double[nc];
            yc = new double[nc];
            wc = new double[nc];

            quaequad0(mmax, nc, ref xc, ref yc, ref wc);
            //
            //  Expand the nodes to the entire triangle.
            //
            nf = QuadratureRule.rule_full_size(mmax);
            xf = new double[nf];
            yf = new double[nf];
            wf = new double[nf];

            if (itype == 0)
            {
                quaenodes(nc, xc, yc, wc, ref nf, ref xf, ref yf, ref wf);
            }
            //
            //  Expand the nodes to the lower 1/3 of the triangle.
            //
            else if (itype == 1)
            {
                quaenodes2(nc, xc, yc, wc, ref nf, ref xf, ref yf, ref wf);
            }
            //
            //  Simply copy the nodes; they are already in the lower-left
            //  1/6 of the triangle.
            //
            else if (itype == 2)
            {
                typeMethods.r8vec_copy(nc, xc, ref xf);
                typeMethods.r8vec_copy(nc, yc, ref yf);
                typeMethods.r8vec_copy(nc, wc, ref wf);
            }

            quaecopy(nf, xf, yf, wf, ref zs, ref whts);
        }

        public static void quaequad0(int mmax, int kk, ref double[] xnew, ref double[] ynew, ref double[] w)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    QUAEQUAD0 returns the requested quadrature rule.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU GPL license.
            //
            //  Modified:
            //
            //    30 June 2014
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
            //    This C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Hong Xiao, Zydrunas Gimbutas,
            //    A numerical algorithm for the construction of efficient quadrature
            //    rules in two and higher dimensions,
            //    Computers and Mathematics with Applications,
            //    Volume 59, 2010, pages 663-676.
            //
            //  Parameters:
            //
            //    Input, int MMAX, the degree of the quadrature (the 
            //    maximum degree of the polynomials of two variables that are integrated
            //    exactly.  1 <= MMAX <= 50.
            //
            //    Input, int KK, the size of the compressed rule.
            //
            //    Output, double XNEW[KK], YNEW[KK], the coordinates of the nodes.
            //
            //    Output, double W[KK], the weights.
            //
        {
            int i;
            int iitype;
            int nbool0;
            int nbool1;
            int nbool2;
            double x0 = 0;
            double x1 = 0;
            double y0 = 0;
            double y1 = 0;
            //
            //  Copy the arrays defining the compressed rule.
            //
            if (mmax == 1)
            {
                QuadratureRule.rule01(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 2)
            {
                QuadratureRule.rule02(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 3)
            {
                QuadratureRule.rule03(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 4)
            {
                QuadratureRule.rule04(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 5)
            {
                QuadratureRule.rule05(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 6)
            {
                QuadratureRule.rule06(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 7)
            {
                QuadratureRule.rule07(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 8)
            {
                QuadratureRule.rule08(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 9)
            {
                QuadratureRule.rule09(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 10)
            {
                QuadratureRule.rule10(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 11)
            {
                QuadratureRule.rule11(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 12)
            {
                QuadratureRule.rule12(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 13)
            {
                QuadratureRule.rule13(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 14)
            {
                QuadratureRule.rule14(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 15)
            {
                QuadratureRule.rule15(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 16)
            {
                QuadratureRule.rule16(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 17)
            {
                QuadratureRule.rule17(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 18)
            {
                QuadratureRule.rule18(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 19)
            {
                QuadratureRule.rule19(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 20)
            {
                QuadratureRule.rule20(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 21)
            {
                QuadratureRule.rule21(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 22)
            {
                QuadratureRule.rule22(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 23)
            {
                QuadratureRule.rule23(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 24)
            {
                QuadratureRule.rule24(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 25)
            {
                QuadratureRule.rule25(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 26)
            {
                QuadratureRule.rule26(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 27)
            {
                QuadratureRule.rule27(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 28)
            {
                QuadratureRule.rule28(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 29)
            {
                QuadratureRule.rule29(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 30)
            {
                QuadratureRule.rule30(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 31)
            {
                QuadratureRule.rule31(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 32)
            {
                QuadratureRule.rule32(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 33)
            {
                QuadratureRule.rule33(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 34)
            {
                QuadratureRule.rule34(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 35)
            {
                QuadratureRule.rule35(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 36)
            {
                QuadratureRule.rule36(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 37)
            {
                QuadratureRule.rule37(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 38)
            {
                QuadratureRule.rule38(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 39)
            {
                QuadratureRule.rule39(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 40)
            {
                QuadratureRule.rule40(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 41)
            {
                QuadratureRule.rule41(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 42)
            {
                QuadratureRule.rule42(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 43)
            {
                QuadratureRule.rule43(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 44)
            {
                QuadratureRule.rule44(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 45)
            {
                QuadratureRule.rule45(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 46)
            {
                QuadratureRule.rule46(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 47)
            {
                QuadratureRule.rule47(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 48)
            {
                QuadratureRule.rule48(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 49)
            {
                QuadratureRule.rule49(ref xnew, ref ynew, ref w);
            }
            else if (mmax == 50)
            {
                QuadratureRule.rule50(ref xnew, ref ynew, ref w);
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("QUAEQUAD0 - Fatal error!");
                Console.WriteLine("  Illegal input value of MMAX.");
                Console.WriteLine("  1 <= MMAX <= 50 required.");
                return;
            }

            for (i = 0; i < kk; i++)
            {
                //
                //  The lower-left 1/6
                //
                iitype = 2;
                nbool2 = quaeinside(iitype, xnew[i], ynew[i]);
                //
                //  The lower 1/3
                //
                iitype = 1;
                nbool1 = quaeinside(iitype, xnew[i], ynew[i]);
                //
                //  The whole triangle
                //
                iitype = 0;
                nbool0 = quaeinside(iitype, xnew[i], ynew[i]);

                if (nbool2 == 1)
                {
                }
                else if (nbool1 == 1)
                {
                    xnew[i] = -xnew[i];
                    ynew[i] = ynew[i];
                }
                else if (nbool0 == 1)
                {
                    x0 = xnew[i];
                    y0 = ynew[i];
                    quaerotate(x0, y0, ref x1, ref y1);
                    xnew[i] = x1;
                    ynew[i] = y1;
                }
                else
                {
                    Console.WriteLine("");
                    Console.WriteLine("QUAEQUAD0 - Fatal error!");
                    Console.WriteLine("  Point does not lie inside triangle.");
                    return;
                }
            }
        }

        public static void quaerotate(double xin, double yin, ref double xout, ref double yout)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    QUAEROTATE applies a rotation.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU GPL license.
            //
            //  Modified:
            //
            //    30 June 2014
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
            //    This C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Hong Xiao, Zydrunas Gimbutas,
            //    A numerical algorithm for the construction of efficient quadrature
            //    rules in two and higher dimensions,
            //    Computers and Mathematics with Applications,
            //    Volume 59, 2010, pages 663-676.
            //
            //  Parameters:
            //
            //    Input, double XIN, YIN, the coordinates of the point.
            //
            //    Output, double &XOUT, &YOUT, the coordinates of the point
            //    after rotation.
            //
        {
            double a11;
            double a12;
            double a21;
            double a22;
            const double r8_pi = 3.141592653589793;
            double theta;
            //
            //  Initialize the matrix of rotation.
            //
            theta = 2.0 * r8_pi / 3.0;
            a11 = Math.Cos(theta);
            a22 = Math.Cos(theta);
            a12 = -Math.Sin(theta);
            a21 = -a12;
            //
            //  Apply the rotation matrix to the input vector.
            //
            xout = a11 * xin + a12 * yin;
            yout = a21 * xin + a22 * yin;

        }
    }
}