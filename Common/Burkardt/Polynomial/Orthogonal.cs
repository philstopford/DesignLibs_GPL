using System;

namespace Burkardt.PolynomialNS;

public static class Orthogonal
{
    public static double[] ortho2eva(int mmax, double[] z)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ORTHO2EVA evaluates orthogonal polynomials on the reference triangle.
        //
        //  Discussion:
        //
        //    This routine evaluates at the user-supplied point Z
        //    a collection of polynomials (of X,Y) orthogonal on the
        //    reference triangle.
        //
        //    The reference triangle has the vertices
        //      (0,2/Math.Sqrt(3)), (-1,-1/Math.Sqrt(3)), (1,-1/Math.Sqrt(3))
        //
        //    The polynomials evaluated by this routine are all
        //    orthogonal polynomials up to order mmax, arranged in the
        //    increasing order.
        //
        //    This routine is based on the Koornwinder's representation
        //    of the orthogonal polynomials on the right triangle
        //      (-1,-1), (-1,1), (1,-1)
        //    given by
        //      K_mn(x,y) = P_m ((2*x+1+y)/(1-y)) ((1-y)/2)^m P_n^{2m+1,0} (y)
        //    where P_m are the Legendre polynomials or order m
        //    and P_n^{2m+1,0} are the Jacobi polynomials of order n with
        //    the parameters alpha=2*m+1 and beta=0.
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
        //    C++ version by John Burkardt.
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
        //    Input, int MMAX, the maximum order to which the polynomials
        //    are to be evaluated;
        //
        //    Input, double Z[2], the location where the polynomials are
        //    to be evaluated; normally, expected to be inside (including boundary)
        //    the reference triangle.
        //
        //    Output, double POLS[(mmax+1)*(mmax+2)/2], the orthogonal
        //    polynomials evaluated at the point Z.
        //
    {
        double[] pols;

        switch (mmax)
        {
            case 0:
                pols = new double[1];

                pols[0] = 1.0 / Math.Sqrt(3.0) * Math.Sqrt(Math.Sqrt(3.0));
                break;
            case 1:
                pols = new double[3];

                pols[0] = 1.0 / Math.Sqrt(3.0) * Math.Sqrt(Math.Sqrt(3.0));
                pols[1] = z[0] * Math.Sqrt(2.0) * Math.Sqrt(Math.Sqrt(3.0));
                pols[2] = z[1] * Math.Sqrt(2.0) * Math.Sqrt(Math.Sqrt(3.0));
                break;
            default:
                pols = ortho2eva0(mmax, z);
                break;
        }

        return pols;
    }

    public static double[] ortho2eva0(int mmax, double[] z)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ORTHO2EVA0 evaluates the orthonormal polynomials on the triangle.
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
        //    C++ version by John Burkardt.
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
        //    Input, int MMAX, the maximum order to which the polynomials are
        //    to be evaluated.
        //
        //    Input, double Z[2], the coordinates of the evaluation point.
        //
        //    Output, double ORTHO2EVA0[(mmax+1)*(mmax+2)/2], the orthogonal
        //    polynomials evaluated at the point Z.
        //
    {
        int m;

        int npols = (mmax + 1) * (mmax + 2) / 2;
        double[] pols = new double[npols];

        const double zero = 0.0;
        double sqrt3 = Math.Sqrt(3.0);
        const double r11 = -1.0 / 3.0;
        double r12 = -1.0 / sqrt3;
        const double r21 = -1.0 / 3.0;
        double r22 = 2.0 / sqrt3;

        double a = z[0];
        double b = z[1];
        //
        //  Map the reference triangle to the right
        //  triangle with the vertices (-1,-1), (1,-1), (-1,1)
        //
        double x = r11 + r12 * b + a;
        double y = r21 + r22 * b;
        //
        //  Evaluate the Koornwinder's polynomials via the three term recursion.
        //
        double par1 = (2.0 * x + 1.0 + y) / 2.0;
        double par2 = (1.0 - y) / 2.0;

        double[] f1 = LegendreScaled.klegeypols(par1, par2, mmax);

        double[] f2 = new double[(mmax + 1) * (mmax + 1)];

        for (m = 0; m <= mmax; m++)
        {
            par1 = 2 * m + 1;
            Jacobi.kjacopols(y, par1, zero, mmax - m, ref f2, polsIndex: + m * (mmax + 1));
        }

        int kk = 0;
        for (m = 0; m <= mmax; m++)
        {
            int n;
            for (n = 0; n <= m; n++)
            {
                //
                //  Evaluate the polynomial (m-n, n)
                //
                pols[kk] = f1[m - n] * f2[n + (m - n) * (mmax + 1)];
                //
                //  Normalize.
                //
                double scale = Math.Sqrt
                ((1 + (m - n) + n) * (1 + (m - n) + (m - n)) / sqrt3
                );

                pols[kk] *= scale;
                kk += 1;
            }
        }

        return pols;
    }
}