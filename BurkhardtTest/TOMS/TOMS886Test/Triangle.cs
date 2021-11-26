using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Function;
using Burkardt.Interpolation;
using Burkardt.Types;

namespace TOMS886Test;

public static class TriangleTest
{
    public static void triangle()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TRIANGLE.
        //
        //  Discussion:
        //
        //    This driver computes the interpolation of the Franke function
        //    on the triangle T(U,V,W) with vertices U=(U1,U2)=(0,0), 
        //    V=(V1,V2)=(1,0) and W=(W1,W2)=(0,1) (unit triangle) 
        //    at the first family of Padua points. 
        //
        //    The degree of interpolation is DEG = 60 and the number of target 
        //    points is NTG = NTG1 ** 2 - NTG1 + 1, NTG1 = 100.
        //
        //    The maps from the reference square [-1,1]^2 to the triangle
        //    are SIGMA1 and SIGMA2 with inverses ISIGM1 and ISIGM2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //  
        //  Modified:
        //
        //    16 February 2014
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
        //    Marco Vianello.
        //    This C version by John Burkardt.
        //
        //  Reference:
        //
        //    Marco Caliari, Stefano de Marchi, Marco Vianello,
        //    Algorithm 886:
        //    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
        //    ACM Transactions on Mathematical Software,
        //    Volume 35, Number 3, October 2008, Article 21, 11 pages.
        //
        //  Parameters:
        //
        //    Local, int DEGMAX, the maximum degree of interpolation.
        //
        //    Local, int NPDMAX, the maximum number of Padua points
        //    = (DEGMAX + 1) * (DEGMAX + 2) / 2.
        //
        //    Local, int NTG1MX, the maximum value of the parameter determining 
        //    the number of target points.
        //
        //    Local, int NTGMAX, the maximum number of target points,
        //    dependent on NTG1MX.
        //
        //    Local, int DEG, the degree of interpolation.
        //
        //    Local, int NTG1, the parameter determining the number 
        //   of target points.
        //
        //    Local, int NPD, the number of Padua points = (DEG + 1) * (DEG + 2) / 2.
        //
        //    Local, int NTG, the number of target points, dependent on NTG1.
        //
        //    Local, double PD1(NPDMAX), the first coordinates of 
        //    the Padua points.
        //
        //    Local, double PD2(NPDMAX), the second coordinates of the 
        //    Padua points.
        //
        //    Local, double WPD(NPDMAX), the weights.
        //
        //    Local, double FPD(NPDMAX), the function at the Padua points.
        //
        //    Workspace, double RAUX1(DEGMAX+1)*(DEGMAX+2)).
        //
        //    Workspace, double RAUX2(DEGMAX+1)*(DEGMAX+2)).
        //
        //    Local, double C0(0:DEGMAX+1,0:DEGMAX+1), the coefficient matrix.
        //
        //    Local, double TG1(NTGMAX), the first coordinates of the 
        //    target points.
        //
        //    Local, double TG2(NTGMAX), the second coordinates of the 
        //    target points.
        //
        //    Local, double INTFTG(NTGMAX), the values of the 
        //    interpolated function.
        //
        //    Local, double MAXERR, the maximum norm of the error at target 
        //    points.
        //
        //    Local, double ESTERR, the estimated error.
        //
    {
        const int DEGMAX = 60;
        const int NTG1MX = 100;
        const int NPDMAX = (DEGMAX + 1) * (DEGMAX + 2) / 2;
        const int NTGMAX = NTG1MX * NTG1MX - NTG1MX + 1;

        double[] c0 = new double[(DEGMAX + 2) * (DEGMAX + 2)];
        double esterr = 0;
        double[] fpd = new double[NPDMAX];
        int i;
        double[] intftg = new double[NTGMAX];
        double maxdev;
        int npd = 0;
        int ntg = 0;
        List<string> output = new();
        double[] pd1 = new double[NPDMAX];
        double[] pd2 = new double[NPDMAX];
        double[] raux1 = new double[(DEGMAX + 1) * (DEGMAX + 2)];
        double[] raux2 = new double[(DEGMAX + 1) * (DEGMAX + 2)];
        double[] tg1 = new double[NTGMAX];
        double[] tg2 = new double[NTGMAX];
        double[] wpd = new double[NPDMAX];
        double x;
        double y;

        const double u1 = 0.0;
        const double u2 = 0.0;
        const double v1 = 1.0;
        const double v2 = 0.0;
        const double w1 = 0.0;
        const double w2 = 1.0;
        const int deg = 60;
        const int ntg1 = 100;

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE:");
        Console.WriteLine("  Interpolation of the Franke function");
        Console.WriteLine("  on the unit triangle T((0,0),(1,0),(0,1))");
        Console.WriteLine("  at degree = " + deg + "");

        //
        //  Build the first family of Padua points in the square [-1,1]^2
        // 
        Padua.pdpts(deg, ref pd1, ref pd2, ref wpd, ref npd);
        // 
        //  Compute the Franke function at Padua points mapped to T(U,V,W).
        //
        for (i = 0; i < npd; i++)
        {
            x = sigma1(pd1[i], pd2[i], u1, u2, v1, v2, w1, w2);
            y = sigma2(pd1[i], pd2[i], u1, u2, v1, v2, w1, w2);
            fpd[i] = Franke.franke(x, y);
        }

        //
        //  Write X, Y, F(X,Y) to a file.
        //
        string filename = "triangle_fpd.txt";
        for (i = 0; i < npd; i++)
        {
            x = sigma1(pd1[i], pd2[i], u1, u2, v1, v2, w1, w2);
            y = sigma2(pd1[i], pd2[i], u1, u2, v1, v2, w1, w2);
            output.Add(x
                       + "  " + y
                       + "  " + fpd[i] + "");
        }

        File.WriteAllLines(filename, output);
        Console.WriteLine("");
        Console.WriteLine("  Wrote F(x,y) at Padua points in '" + filename + "'");
        //    
        //  Compute the matrix C0 of the coefficients in the bivariate
        //  orthonormal Chebyshev basis
        //
        Padua.padua2(deg, DEGMAX, npd, wpd, fpd, raux1, raux2, ref c0, ref esterr);
        //    
        //  Build the set of target points on T(U,V,W)
        //
        target(u1, u2, v1, v2, w1, w2, ntg1, NTGMAX, ref tg1, ref tg2, ref ntg);
        //    
        //  Evaluate the interpolant at the target points.
        //
        for (i = 0; i < ntg; i++)
        {
            x = isigm1(tg1[i], tg2[i], u1, u2, v1, v2, w1, w2);
            y = isigm2(tg1[i], tg2[i], u1, u2, v1, v2, w1, w2);
            intftg[i] = Padua.pd2val(deg, DEGMAX, c0, x, y);
        }

        //
        //  Write the function value at target points to a file.
        //
        filename = "triangle_ftg.txt";
        for (i = 0; i < ntg; i++)
        {
            output.Add(tg1[i]
                       + "  " + tg2[i]
                       + "  " + Franke.franke(tg1[i], tg2[i]) + "");
        }

        File.WriteAllLines(filename, output);
        Console.WriteLine("  Wrote F(x,y) at target points in '" + filename + "'");
        //
        //  Write the interpolated function value at target points to a file.
        //
        filename = "triangle_itg.txt";
        for (i = 0; i < ntg; i++)
        {
            output.Add(tg1[i]
                       + "  " + tg2[i]
                       + "  " + intftg[i] + "");
        }

        File.WriteAllLines(filename, output);
        Console.WriteLine("  Wrote I(F)(x,y) at target points in '" + filename + "'");
        //
        //  Compute the error relative to the max deviation from the mean.
        //
        double maxerr = 0.0;
        double mean = 0.0;
        double fmax = -typeMethods.r8_huge();
        double fmin = +typeMethods.r8_huge();

        for (i = 0; i < ntg; i++)
        {
            double fxy = Franke.franke(tg1[i], tg2[i]);
            double ixy = intftg[i];
            maxerr = Math.Max(maxerr, Math.Abs(fxy - ixy));
            mean += fxy;
            fmax = Math.Max(fmax, fxy);
            fmin = Math.Min(fmin, fxy);
        }

        if (Math.Abs(fmax - fmin) <= double.Epsilon)
        {
            maxdev = 1.0;
        }
        else
        {
            mean /= ntg;
            maxdev = Math.Max(fmax - mean, mean - fmin);
        }

        //
        //  Print error ratios.
        //
        Console.WriteLine("");
        Console.WriteLine("  Estimated error:  " + esterr / maxdev + "");
        Console.WriteLine("  Actual error:     " + maxerr / maxdev + "");
        Console.WriteLine("  Expected error:   " + 0.1226E-09 + "");
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("TRIANGLE:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static double sigma1(double t1, double t2, double u1, double u2, double v1,
            double v2, double w1, double w2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SIGMA1 maps first coordinate from square to triangle.
        //
        //  Discussion:
        //
        //    This function returns the first component of the map
        //    from the square [-1,1]^2 to the triangle T(U,V,W). 
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //  
        //  Modified:
        //
        //    16 February 2014
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
        //    Marco Vianello.
        //    This C version by John Burkardt.
        //
        //  Reference:
        //
        //    Marco Caliari, Stefano de Marchi, Marco Vianello,
        //    Algorithm 886:
        //    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
        //    ACM Transactions on Mathematical Software,
        //    Volume 35, Number 3, October 2008, Article 21, 11 pages.
        //
        //  Parameters:
        //
        //    Input, double T1, T2, the coordinates of a point in the square.
        //
        //    Input, double U1, U2, V1, V2, W1, W2, the coordinates of the 
        //    vertices of the triangle.
        //
        //    Output, double SIGMA1, the X coordinate of the corresponding
        //    point in the triangle.
        //
    {
        double value = 0;

        value = (v1 - u1) * (1.0 + t1)
                          * (1.0 - t2) / 4.0
                + (w1 - u1) * (1.0 + t2) / 2.0 + u1;

        return value;
    }

    private static double isigm1(double sigma1, double sigma2, double u1, double u2, double v1,
            double v2, double w1, double w2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ISIGM1 maps first coordinate from triangle to the square.
        //
        //  Discussion:
        //
        //    This functions returns the first component of the map
        //    from the the triangle T(U,V,W) to the square [-1,1]^2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //  
        //  Modified:
        //
        //    16 February 2014
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
        //    Marco Vianello.
        //    This C version by John Burkardt.
        //
        //  Reference:
        //
        //    Marco Caliari, Stefano de Marchi, Marco Vianello,
        //    Algorithm 886:
        //    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
        //    ACM Transactions on Mathematical Software,
        //    Volume 35, Number 3, October 2008, Article 21, 11 pages.
        //
        //  Parameters:
        //
        //    Input, double SIGMA1, SIGMA2, the coordinates of a point 
        //    in the triangle.
        //
        //    Input, double U1, U2, V1, V2, W1, W2, the coordinates of the 
        //    vertices of the triangle.
        //
        //    Output, double ISIGM1, the X coordinate of the corresponding
        //    point in the square.
        //
    {
        double rho1 = (sigma1 * (w2 - u2) - sigma2 * (w1 - u1)
                          + (w1 - u1) * u2 - (w2 - u2) * u1) /
                      ((v1 - u1) * (w2 - u2) - (v2 - u2) * (w1 - u1));

        double rho2 = (sigma1 * (v2 - u2) - sigma2 * (v1 - u1)
                          + (v1 - u1) * u2 - (v2 - u2) * u1) /
                      ((w1 - u1) * (v2 - u2) - (w2 - u2) * (v1 - u1));

        double value = rho2 switch
        {
            1.0 => 0.0,
            _ => 2.0 * rho1 / (1.0 - rho2) - 1.0
        };

        return value;
    }

    private static double sigma2(double t1, double t2, double u1, double u2, double v1,
            double v2, double w1, double w2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SIGMA2 maps the second coordinate from square to triangle.
        //
        //  Discussion:
        //
        //    This functions returns the second component of the map
        //    from the square [-1,1]^2 to the triangle T(U,V,W).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //  
        //  Modified:
        //
        //    16 February 2014
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
        //    Marco Vianello.
        //    This C version by John Burkardt.
        //
        //  Reference:
        //
        //    Marco Caliari, Stefano de Marchi, Marco Vianello,
        //    Algorithm 886:
        //    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
        //    ACM Transactions on Mathematical Software,
        //    Volume 35, Number 3, October 2008, Article 21, 11 pages.
        //
        //  Parameters:
        //
        //    Input, double T1, T2, the coordinates of a point in the square.
        //
        //    Input, double U1, U2, V1, V2, W1, W2, the coordinates of the 
        //    vertices of the triangle.
        //
        //    Output, double SIGMA2, the Y coordinate of the corresponding
        //    point in the triangle.
        //
    {
        double value = 0;

        value = (v2 - u2) * (1.0 + t1)
                          * (1.0 - t2) / 4.0 + (w2 - u2)
            * (1.0 + t2) / 2.0 + u2;

        return value;
    }

    private static double isigm2(double sigma1, double sigma2, double u1, double u2, double v1,
            double v2, double w1, double w2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ISIGM2 maps second coordinate from triangle to the square.
        //
        //  Discussion:
        //
        //    This functions returns the second component of the map
        //    from the the triangle T(U,V,W) to the square [-1,1]^2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //  
        //  Modified:
        //
        //    16 February 2014
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
        //    Marco Vianello.
        //    This C version by John Burkardt.
        //
        //  Reference:
        //
        //    Marco Caliari, Stefano de Marchi, Marco Vianello,
        //    Algorithm 886:
        //    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
        //    ACM Transactions on Mathematical Software,
        //    Volume 35, Number 3, October 2008, Article 21, 11 pages.
        //
        //  Parameters:
        //
        //    Input, double SIGMA1, SIGMA2, the coordinates of a point 
        //    in the ellipse.
        //
        //    Input, double U1, U2, V1, V2, W1, W2, the coordinates of the 
        //    vertices of the triangle.
        //
        //    Output, double ISIGM2, the Y coordinate of the corresponding
        //    point in the triangle.
        //
    {
        double rho2 = (sigma1 * (v2 - u2) -
                          sigma2 * (v1 - u1) + (v1 - u1) * u2 - (v2 - u2) * u1) /
                      ((w1 - u1) * (v2 - u2) - (w2 - u2) * (v1 - u1));

        double value = rho2 switch
        {
            1.0 => 1.0,
            _ => 2.0 * rho2 - 1.0
        };

        return value;
    }

    private static void target(double u1, double u2, double v1, double v2, double w1, double w2,
            int ntg1, int ntgmax, ref double[] tg1, ref double[] tg2, ref int ntg)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TARGET returns the target points on the triangle.
        //
        //  Discussion:
        //
        //    Target points on the triangle T(U,V,W).
        //    The number of target points is NTG = NTG1^2 - NGT1 + 1.
        //
        //  Licensing:
        //
        //    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
        //    Marco Vianello.
        //    This C version by John Burkardt.
        //  
        //  Modified:
        //
        //    16 February 2014
        //
        //  Author:
        //
        //    Marco Caliari, Stefano De Marchi, Marco Vianello
        //
        //  Reference:
        //
        //    Marco Caliari, Stefano de Marchi, Marco Vianello,
        //    Algorithm 886:
        //    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
        //    ACM Transactions on Mathematical Software,
        //    Volume 35, Number 3, October 2008, Article 21, 11 pages.
        //
        //  Parameters:
        //
        //    Input, double U1, U2, V1, V2, W1, W2, the coordinates of the 
        //    vertices of the triangle.
        //
        //    Input, int NTG1, a parameter determining the number 
        //    of target points
        //
        //    Input, int NTGMAX, the maximum number of target points.
        //
        //    Output, double TG1[NTG], TG2[NTG], the X and Y coordinates
        //    of the target points.
        //
        //    Output, int &NTG, the number of target points computed.
        //
    {
        int i;
        int j;

        switch (ntg1)
        {
            case < 2:
                Console.WriteLine("");
                Console.WriteLine("TARGET - Fatal error!");
                Console.WriteLine("  NTG1 < 2");
                Console.WriteLine("  NTG1 = " + ntg1 + "");
                return;
        }

        if (ntgmax < ntg1 * ntg1 - ntg1 + 1)
        {
            Console.WriteLine("");
            Console.WriteLine("TARGET - Fatal error!");
            Console.WriteLine("  NTGMAX < NTG1 * NTG1 - NTG1 + 1.");
            Console.WriteLine("  NTG1 = " + ntg1 + "");
            Console.WriteLine("  NTGMAX = " + ntgmax + "");
            return;
        }

        ntg = 0;

        for (i = 1; i <= ntg1 - 1; i++)
        {
            for (j = 1; j <= ntg1; j++)
            {
                tg1[ntg] = (v1 - u1) * (i - 1) / (ntg1 - 1)
                           + (w1 - u1) * ((j - 1) / (double) (ntg1 - 1))
                                       * (1.0 - (i - 1) / (double) (ntg1 - 1)) + u1;

                tg2[ntg] = (v2 - u2) * (i - 1) / (ntg1 - 1)
                           + (w2 - u2) * ((j - 1) / (double) (ntg1 - 1))
                                       * (1.0 - (i - 1) / (double) (ntg1 - 1)) + u2;

                ntg += 1;
            }
        }

        i = ntg1;
        j = 1;

        tg1[ntg] = (v1 - u1) * (i - 1) / (ntg1 - 1)
                   + (w1 - u1) * ((j - 1) / (double) (ntg1 - 1))
                               * (1.0 - (i - 1) / (double) (ntg1 - 1)) + u1;

        tg2[ntg] = (v2 - u2) * (i - 1) / (ntg1 - 1)
                   + (w2 - u2) * ((j - 1) / (double) (ntg1 - 1))
                               * (1.0 - (i - 1) / (double) (ntg1 - 1)) + u2;

        ntg += 1;

    }
}