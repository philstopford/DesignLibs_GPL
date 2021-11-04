using System;
using System.Collections.Generic;
using System.Drawing;
using System.IO;
using Burkardt;
using Burkardt.Function;
using Burkardt.Interpolation;
using Burkardt.Types;

namespace TOMS886Test
{
    class RectangleTest
    {
        public static void rectangle()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for RECTANGLE.
            //
            //  Discussion:
            //
            //    This driver computes the interpolation of the Franke function
            //    on the rectangle R(A,B) = [A1,B1] x [A2,B2] with A=(A1,A2)=(0,0) 
            //    and B=(B1,B2)=(1,1) (unit square) at the FAMILY = 1 of Padua points. 
            //
            //    The degree of interpolation is DEG = 60 and the number of target 
            //    points is NTG = NTG1^2, NTG1 = 100. 
            //
            //    The maps from the reference square [-1,1]^2 to the rectangle
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
            //    C++ version by John Burkardt.
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
            //    Local, int DEG, the degree of interpolation,
            //
            //    Local, int NTG1, the parameter determining the number 
            //    of target points.
            //
            //    Local, int FAMILY, specifies the desired family of Padua points.
            //
            //    Local, int NPD, the number of Padua points = (DEG + 1) * (DEG + 2) / 2.
            //
            //    Local, int NTG, the number of target points, dependent on NTG1.
            //
            //    Local, double PD1[NPDMAX], the first coordinates of 
            //    the Padua points.
            //
            //    Local, double PD2[NPDMAX], the second coordinates of the 
            //    Padua points.
            //
            //    Local, double WPD[NPDMAX], the weights.
            //
            //    Local, double FPD[NPDMAX], the function at the Padua points.
            //
            //    Workspace, double RAUX1[(DEGMAX+1)*(DEGMAX+2)].
            //
            //    Workspace, double RAUX2[(DEGMAX+1)*(DEGMAX+2)].
            //
            //    Local, double C0[(0:DEGMAX+1)*(0:DEGMAX+1)], the coefficient matrix.
            //
            //    Local, double TG1[NTGMAX], the first coordinates of the 
            //    target points.
            //
            //    Local, double TG2[NTGMAX], the second coordinates of the 
            //    target points.
            //
            //    Local, double INTFTG[NTGMAX], the values of the 
            //    interpolated function.
            //
            //    Local, double MAXERR, the maximum norm of the error at target 
            //    points.
            //
            //    Local, double ESTERR, the estimated error.
            //
        {
            int DEGMAX = 60;
            int NTG1MX = 100;
            int NPDMAX = ((DEGMAX + 1) * (DEGMAX + 2) / 2);
            int NTGMAX = (NTG1MX * NTG1MX);

            double a1;
            double a2;
            double b1;
            double b2;
            double[] c0 = new double[(DEGMAX + 2) * (DEGMAX + 2)];
            int deg;
            int degmax = DEGMAX;
            double esterr = 0;
            int family;
            string filename;
            double fmax;
            double fmin;
            double[] fpd = new double[NPDMAX];
            double fxy;
            int i;
            double[] intftg = new double[NTGMAX];
            double ixy;
            double maxdev;
            double maxerr;
            double mean;
            int npd = 0;
            int ntg = 0;
            int ntg1;
            int ntgmax = NTGMAX;
            List<string> output = new List<string>();
            double[] pd1 = new double[NPDMAX];
            double[] pd2 = new double[NPDMAX];
            double[] raux1 = new double[(DEGMAX + 1) * (DEGMAX + 2)];
            double[] raux2 = new double[(DEGMAX + 1) * (DEGMAX + 2)];
            double[] tg1 = new double[NTGMAX];
            double[] tg2 = new double[NTGMAX];
            double[] wpd = new double[NPDMAX];
            double x;
            double y;

            a1 = 0.0;
            a2 = 0.0;
            b1 = 1.0;
            b2 = 1.0;
            family = 1;
            deg = 60;
            ntg1 = 100;

            Console.WriteLine("");
            Console.WriteLine("RECTANGLE:");
            Console.WriteLine("  Interpolation of the Franke function");
            Console.WriteLine("  on the unit square [0,1] x [0,1]");
            Console.WriteLine("  of degree = " + deg + "");

            if (degmax < deg)
            {
                Console.WriteLine("");
                Console.WriteLine("RECTANGLE - Fatal error!");
                Console.WriteLine("  DEGMAX < DEG.");
                Console.WriteLine("  DEG =    " + deg + "");
                Console.WriteLine("  DEGMAX = " + degmax + "");
                return;
            }

            //
            //  Build the first family of Padua points in the square [-1,1]^2
            //
            Padua.pdpts(deg, ref pd1, ref pd2, ref wpd, ref npd);
            //    
            //  Compute the Franke function at Padua points mapped to the region.
            //
            for (i = 0; i < npd; i++)
            {
                x = sigma1(pd1[i], pd2[i], a1, a2, b1, b2, family, deg);
                y = sigma2(pd1[i], pd2[i], a1, a2, b1, b2, family, deg);
                fpd[i] = Franke.franke(x, y);
            }

            //
            //  Write X, Y, F(X,Y) to a file.
            //
            filename = "rectangle_fpd.txt";
            output.Clear();
            for (i = 0; i < npd; i++)
            {
                x = sigma1(pd1[i], pd2[i], a1, a2, b1, b2, family, deg);
                y = sigma2(pd1[i], pd2[i], a1, a2, b1, b2, family, deg);
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
            Padua.padua2(deg, degmax, npd, wpd, fpd, raux1, raux2, ref c0, ref esterr);
            //    
            //  Evaluate the target points in the region.
            //
            target(a1, b1, a2, b2, ntg1, ntgmax, ref tg1, ref tg2, ref ntg);
            //    
            //  Evaluate the interpolant at the target points.
            //
            for (i = 0; i < ntg; i++)
            {
                x = isigm1(tg1[i], tg2[i], a1, a2, b1, b2, family, deg);
                y = isigm2(tg1[i], tg2[i], a1, a2, b1, b2, family, deg);
                intftg[i] = Padua.pd2val(deg, degmax, c0, x, y);
            }

            //
            //  Write the function value at target points to a file.
            //
            filename = "rectangle_ftg.txt";

            output.Clear();
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
            filename = "ellipse_itg.txt";

            output.Clear();
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
            maxerr = 0.0;
            mean = 0.0;
            fmax = -typeMethods.r8_huge();
            fmin = +typeMethods.r8_huge();

            for (i = 0; i < ntg; i++)
            {
                fxy = Franke.franke(tg1[i], tg2[i]);
                ixy = intftg[i];
                maxerr = Math.Max(maxerr, Math.Abs(fxy - ixy));
                mean = mean + fxy;
                fmax = Math.Max(fmax, fxy);
                fmin = Math.Min(fmin, fxy);
            }

            if (fmax == fmin)
            {
                maxdev = 1.0;
            }
            else
            {
                mean = mean / (double) (ntg);
                maxdev = Math.Max(fmax - mean, mean - fmin);
            }

            //
            //  Print error ratios.
            //
            Console.WriteLine("");
            Console.WriteLine("  Estimated error:  " + esterr / maxdev + "");
            Console.WriteLine("  Actual error:     " + maxerr / maxdev + "");
            Console.WriteLine("  Expected error:   " + 0.2468E-10 + "");

            Console.WriteLine("");
            Console.WriteLine("RECTANGLE:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static double sigma1(double t1, double t2, double a1, double a2, double b1,
                double b2, int family, int deg)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIGMA1 maps first coordinate from square to the rectangle.
            //
            //  Discussion:
            //
            //    This function returns the first component of the map 
            //    from the square [-1,1]^2 to the rectangle [A1,B1] x [A2,B2]. 
            //    FAMILY and DEG select the rotation in order to get 
            //    the corresponding FAMILY of Padua points at degree DEG.
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
            //    C++ version by John Burkardt.
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
            //    Input, double A1, B1, A2, B2, the coordinates of the extreme
            //    corners of the rectangle.
            //
            //    Input, int FAMILY, DEG, select the family of Padua points at 
            //    degree DEG.
            //
            //    Output, double SIGMA1, the X coordinate of the corresponding
            //    point in the rectangle.
            //
        {
            double pi = 3.1415926535897931;
            double theta;
            double value;

            theta = (double) (2 * (deg % 2) - 1)
                * (double) (family - 1) * pi / 2.0;
            value = t1 * Math.Cos(theta) - t2 * Math.Sin(theta);
            value = ((b1 - a1) * value + (b1 + a1)) / 2.0;

            return value;
        }

        static double isigm1(double sigma1, double sigma2, double a1, double a2, double b1,
                double b2, int family, int deg)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ISIGM1 maps first coordinate from the rectangle to the square.
            //
            //  Discussion:
            //
            //    This function returns the first component of the map 
            //    from the rectangle [A1,B1] x [A2,B2] to the square [-1,1]^2. 
            //    FAMILY and DEG select the rotation in order to get 
            //    the corresponding FAMILY of Padua points at degree DEG.
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
            //    C++ version by John Burkardt.
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
            //    in the rectangle.
            //
            //    Input, double A1, B1, A2, B2, the coordinates of the extreme
            //    corners of the rectangle.
            //
            //    Input, int FAMILY, DEG, select the family of Padua points at 
            //    degree DEG.
            //
            //    Output, double ISIGM1, the X coordinate of the corresponding
            //    point in the square.
            //
        {
            double isigm2;
            const double pi = 3.1415926535897931;
            double theta;
            double value;

            theta = (double) (2 * (deg % 2) - 1)
                * (double) (family - 1) * pi / 2.0;
            value = (2.0 * sigma1 - (b1 + a1)) / (b1 - a1);
            isigm2 = (2.0 * sigma2 - (b2 + a2)) / (b2 - a2);
            value = value * Math.Cos(theta) + isigm2 * Math.Sin(theta);

            return value;
        }

        static double sigma2(double t1, double t2, double a1, double a2, double b1,
                double b2, int family, int deg)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIGMA2 maps second coordinate from square to the rectangle.
            //
            //  Discussion:
            //
            //    This function returns the second component of the map 
            //    from the square [-1,1]^2 to the rectangle [A1,B1] x [A2,B2]. 
            //    FAMILY and DEG select the rotation in order to get 
            //    the corresponding FAMILY of Padua points at degree DEG.
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
            //    C++ version by John Burkardt.
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
            //    Input, double A1, B1, A2, B2, the coordinates of the extreme
            //    corners of the rectangle.
            //
            //    Input, int FAMILY, DEG, select the family of Padua points at 
            //    degree DEG.
            //
            //    Output, double SIGMA2, the Y coordinate of the corresponding
            //    point in the rectangle.
            //
        {
            const double pi = 3.1415926535897931;
            double theta;
            double value;

            theta = (double) (2 * (deg % 2) - 1)
                * (double) (family - 1) * pi / 2.0;
            value = t1 * Math.Sin(theta) + t2 * Math.Cos(theta);
            value = ((b2 - a2) * value + (b2 + a2)) / 2.0;

            return value;
        }

        static double isigm2(double sigma1, double sigma2, double a1, double a2, double b1,
                double b2, int family, int deg)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ISIGM2 maps the second coordinate from the rectangle to the square.
            //
            //  Discussion:
            //
            //    This function returns the second component of the map 
            //    from the rectangle [A1,B1] x [A2,B2] to the square [-1,1]^2. 
            //
            //    FAMILY and DEG select the rotation in order to get 
            //    the corresponding FAMILY of Padua points at degree DEG.
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
            //    C++ version by John Burkardt.
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
            //    Input, double A1, B1, A2, B2, the coordinates of the extreme
            //    corners of the rectangle.
            //
            //    Input, int FAMILY, DEG, select the family of Padua points at 
            //    degree DEG.
            //
            //    Output, double ISIGM2, the Y coordinate of the corresponding
            //    point in the rectangle.
            //
        {
            double isigm1;
            const double pi = 3.1415926535897931;
            double theta;
            double value;

            theta = (double) (2 * (deg % 2) - 1)
                * (double) (family - 1) * pi / 2.0;
            isigm1 = (2.0 * sigma1 - (b1 + a1)) / (b1 - a1);
            value = (2.0 * sigma2 - (b2 + a2)) / (b2 - a2);
            value = -isigm1 * Math.Sin(theta) + value * Math.Cos(theta);

            return value;
        }

        static void target(double a1, double b1, double a2, double b2, int ntg1, int ntgmax,
                ref double[] tg1, ref double[] tg2, ref int ntg)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TARGET returns the target points on the rectangle.
            //
            //  Discussion:
            //
            //    Target points (uniform grid) on the rectangle [A1,B1] x [A2,B2].
            //    The number of target points is NTG = NTG1^2.
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
            //    C++ version by John Burkardt.
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
            //    Input, double A1, B1, A2, B2, the coordinates of the extreme
            //    corners of the rectangle.
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

            if (ntg1 < 2)
            {
                Console.WriteLine("");
                Console.WriteLine("TARGET - Fatal error!");
                Console.WriteLine("  NTG1 < 2.");
                Console.WriteLine("  NTG1 = " + ntg1 + "");
                return;
            }

            if (ntgmax < ntg1 * ntg1)
            {
                Console.WriteLine("");
                Console.WriteLine("TARGET - Fatal error!");
                Console.WriteLine("  NTGMAX < NTG1 * NTG1.");
                Console.WriteLine("  NTG1 = " + ntg1 + "");
                Console.WriteLine("  NTGMAX = " + ntgmax + "");
                return;
            }

            ntg = 0;

            for (i = 1; i <= ntg1; i++)
            {
                for (j = 1; j <= ntg1; j++)
                {
                    tg1[ntg] = a1 + (double) (j - 1) * (b1 - a1)
                        / (double) (ntg1 - 1);
                    tg2[ntg] = a2 + (double) (i - 1) * (b2 - a2)
                        / (double) (ntg1 - 1);
                    ntg = ntg + 1;
                }
            }

        }
    }
}