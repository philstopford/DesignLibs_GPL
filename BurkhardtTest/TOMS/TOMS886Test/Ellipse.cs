using System;
using System.Collections.Generic;
using System.IO;
using Burkardt;
using Burkardt.Interpolation;
using Burkardt.Types;

namespace TOMS886Test
{
    class EllipseTest
    {
        public static void ellipse()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for ELLIPSE.
            //
            //  Discussion:
            //
            //    This driver computes the interpolation of the Franke function
            //    on the ellipse E((C1,C2),ALPHA,BETA) = E((0.5,0.5),0.5,0.5)  
            //    at the first family of Padua points. 
            //
            //    The ellipse has the equation:
            //
            //      ( ( X - C1 ) / ALPHA )^2 + ( ( Y - C2 ) / BETA )^2 = 1
            //
            //    The degree of interpolation DEG = 60 and the number of target 
            //    points is NTG = NTG1 ^ 2 - 2 * NTG1 + 2, NTG1 = 100.  
            //
            //    The maps from the reference square [-1,1]^2 to the current domain 
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
            //    This C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Marco Caliari, Stefano de Marchi, Marco Vianello,
            //    Algorithm 886:
            //    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
            //    ACM Transactions on Mathematical Software,
            //    Volume 35, Number 3, October 2008, Article 21, 11 pages.
            //
            //  Local Parameters:
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
            //    Local, int NTG1, the parameter determining the number of target points.
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
            int NTGMAX = (NTG1MX * NTG1MX - 2 * NTG1MX + 2);

            double alpha;
            double beta;
            double[] c0 = new double[(DEGMAX + 2) * (DEGMAX + 2)];
            double c1;
            double c2;
            int deg;
            int degmax = DEGMAX;
            double esterr = 0;
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

            alpha = 0.5;
            beta = 0.5;
            c1 = 0.5;
            c2 = 0.5;
            deg = 60;
            ntg1 = 100;

            Console.WriteLine("");
            Console.WriteLine("ELLIPSE:");
            Console.WriteLine("  Interpolation of the Franke function");
            Console.WriteLine("  on the disk with center = (0.5,0.5) and radius = 0.5");
            Console.WriteLine("  of degree = " + deg + "");

            if (degmax < deg)
            {
                Console.WriteLine("");
                Console.WriteLine("ELLIPSE - Fatal error!");
                Console.WriteLine("  DEGMAX < DEG.");
                Console.WriteLine("  DEG =    " + deg + "");
                Console.WriteLine("  DEGMAX = " + degmax + "");
                return;
            }

            //   
            //  Build the first family of Padua points in the square [-1,1]^2.
            //
            Padua.pdpts(deg, ref pd1, ref pd2, ref wpd, ref npd);
            //    
            //  Compute the Franke function at Padua points mapped to the region.
            //
            for (i = 0; i < npd; i++)
            {
                x = sigma1(pd1[i], pd2[i], c1, c2, alpha, beta);
                y = sigma2(pd1[i], pd2[i], c1, c2, alpha, beta);
                fpd[i] = Franke.franke(x, y);
            }

            //
            //  Write X, Y, F(X,Y) to a file.
            //
            filename = "ellipse_fpd.txt";
            for (i = 0; i < npd; i++)
            {
                x = sigma1(pd1[i], pd2[i], c1, c2, alpha, beta);
                y = sigma2(pd1[i], pd2[i], c1, c2, alpha, beta);
                output.Add(x
                           + "  " + y
                           + "  " + fpd[i] + "");
            }

            File.WriteAllLines(filename, output);
            Console.WriteLine("");
            Console.WriteLine("  Wrote F(x,y) at Padua points in '" + filename + "'");
            //
            //  Compute the matrix C0 of the coefficients in the bivariate
            //  orthonormal Chebyshev basis.
            //
            Padua.padua2(deg, degmax, npd, wpd, fpd, raux1, raux2, ref c0, ref esterr);
            //    
            //  Evaluate the target points in the region.
            //
            target(c1, c2, alpha, beta, ntg1, ntgmax, ref tg1, ref tg2, ref ntg);
            //
            //  Evaluate the interpolant at the target points.
            //
            for (i = 0; i < ntg; i++)
            {
                x = isigm1(tg1[i], tg2[i], c1, c2, alpha, beta);
                y = isigm2(tg1[i], tg2[i], c1, c2, alpha, beta);
                intftg[i] = Padua.pd2val(deg, degmax, c0, x, y);
            }

            //
            //  Write the function value at target points to a file.
            //
            output.Clear();
            filename = "ellipse_ftg.txt";
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
            output.Clear();
            filename = "ellipse_itg.txt";
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
            Console.WriteLine("  Expected error:   " + 0.1769E-09 + "");
            //
            //  Terminate.
            //
            Console.WriteLine("");
            Console.WriteLine("ELLIPSE:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static double sigma1(double t1, double t2, double c1, double c2, double alpha,
                double beta)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIGMA1 maps first coordinate from square to ellipse.
            //
            //  Discussion:
            //
            //    This function returns the first component of the map 
            //    from the square [-1,1]^2 to the ellipse E((C1,C2),ALPHA,BETA).
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
            //    This C++ version by John Burkardt.
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
            //    Input, double C1, C2, ALPHA, BETA, the center and scale
            //    parameters of the ellipse.
            //
            //    Output, double SIGMA1, the X coordinate of the corresponding
            //    point in the ellipse.
            //
        {
            double value;

            value = c1 - alpha * t2 * Math.Sin(phi(t1));

            return value;
        }

        static double isigm1(double sigma1, double sigma2, double c1, double c2,
                double alpha, double beta)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ISIGM1 maps the first coordinate from the ellipse to the square.
            //
            //  Discussion:
            //
            //    This function returns the first component of the map 
            //    from the ellipse E((C1,C2),ALPHA,BETA) to the square [-1,1]^2. 
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //  
            //  Modified:
            //
            //    09 February 2014
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
            //    Marco Vianello.
            //    This C++ version by John Burkardt.
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
            //    Input, double C1, C2, ALPHA, BETA, the center and scale
            //    parameters of the ellipse.
            //
            //    Output, double ISIGM1, the X coordinate of the corresponding
            //    point in the square.
            //
        {
            double value;

            if (sigma2 == c2)
            {
                value = 1.0;
            }
            else
            {
                value = iphi(Math.Atan(beta * (c1 - sigma1) /
                                       (alpha * (sigma2 - c2))));
            }

            return value;
        }

        static double sigma2(double t1, double t2, double c1, double c2, double alpha,
                double beta)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIGMA2 maps the second coordinate from square to ellipse.
            //
            //  Discussion:
            //
            //    This function returns the second component of the map 
            //    from the square [-1,1]^2 to the ellipse E((C1,C2),ALPHA,BETA).
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
            //    This C++ version by John Burkardt.
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
            //    Input, double C1, C2, ALPHA, BETA, the center and scale
            //    parameters of the ellipse.
            //
            //    Output, double SIGMA2, the Y coordinate of the corresponding
            //    point in the ellipse.
            //
        {
            double value;

            value = c2 + beta * t2 * Math.Cos(phi(t1));

            return value;
        }

        static double isigm2(double sigma1, double sigma2, double c1, double c2,
                double alpha, double beta)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ISIGM2 maps second coordinate from ellipse to the square.
            //
            //  Discussion:
            //
            //    This function returns the second component of the map 
            //    from the ellipse E((C1,C2),ALPHA,BETA) to the square [-1,1]^2. 
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
            //    This C++ version by John Burkardt.
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
            //    Input, double C1, C2, ALPHA, BETA, the center and scale
            //    parameters of the ellipse.
            //
            //    Output, double ISIGM2, the Y coordinate of the corresponding
            //    point in the square.
            //
        {
            double value;

            if (sigma2 == c2)
            {
                value = (c1 - sigma1) / alpha;
            }
            else
            {
                value = Math.Sqrt(beta * beta * Math.Pow(c1 - sigma1, 2) +
                                  alpha * alpha * Math.Pow(c2 - sigma2, 2))
                    / (alpha * beta) * typeMethods.r8_sign(sigma2 - c2);
            }

            return value;
        }

        static double phi(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PHI maps from [-1,+1] to [-pi/2,+pi/2].
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
            //    This C++ version by John Burkardt.
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
            //    Input, double X, a point in [-1,+1];
            //
            //    Output, double PHI, a corresponding point in [-pi/2,+pi/2].
            //
        {
            const double pi = 3.1415926535897931;
            double value;

            value = pi * x / 2.0;

            return value;
        }

        static double iphi(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    IPHI maps from [-pi/2,+pi/2] to [-1,+1].
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
            //    This C++ version by John Burkardt.
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
            //    Input, double X, a point in [-pi/2,+pi/2].
            //
            //    Output, double IPHI, a corresponding point in [-1,+1].
            //
        {
            const double pi = 3.1415926535897931;
            double value;

            value = 2.0 * x / pi;

            return value;
        }

        static void target(double c1, double c2, double alpha, double beta, int ntg1,
                int ntgmax, ref double[] tg1, ref double[] tg2, ref int ntg)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TARGET returns the target points on the ellipse.
            //
            //  Discussion:
            //
            //    Target points on the ellipse E((C1,C2),ALPHA,BETA).
            //    The number of target points is NTG = NTG1^2 - 2 * NTG1 + 2.
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
            //    This C++ version by John Burkardt.
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
            //    Input, double C1, C2, ALPHA, BETA, the center and scale
            //    parameters of the ellipse.
            //
            //    Input, int NTG1, a parameter determining the number 
            //    of target points.  2 <= NTG1.
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
            double t;

            if (ntg1 < 2)
            {
                Console.WriteLine("");
                Console.WriteLine("TARGET - Fatal error!");
                Console.WriteLine("  NTG1 < 2");
                Console.WriteLine("  NTG1 = " + ntg1 + "");
                return;
            }

            if (ntgmax < ntg1 * ntg1 - 2 * ntg1 + 2)
            {
                Console.WriteLine("");
                Console.WriteLine("TARGET - Fatal error!");
                Console.WriteLine("  NTGMAX < NTG1 * NTG1 - 2 * NTG1 + 2.");
                Console.WriteLine("  NTG1 = " + ntg1 + "");
                Console.WriteLine("  NTGMAX = " + ntgmax + "");
                return;
            }

            i = 1;
            j = 1;
            ntg = 0;

            tg1[ntg] = alpha * (-1.0 + (double) (i - 1) * 2.0
                / (double) (ntg1 - 1)) + c1;

            t = -1.0 + (double) (i - 1) * 2.0 / (double) (ntg1 - 1);

            tg2[ntg] = beta * (-1.0 + (double) (j - 1) * 2.0
                / (double) (ntg1 - 1)) * Math.Sqrt(1.0 - t * t) + c2;

            ntg = ntg + 1;

            for (i = 2; i <= ntg1 - 1; i++)
            {
                for (j = 1; j <= ntg1; j++)
                {
                    tg1[ntg] = alpha * (-1.0 + (double) (i - 1) * 2.0
                        / (double) (ntg1 - 1)) + c1;

                    t = -1.0 + (double) (i - 1) * 2.0 / (double) (ntg1 - 1);

                    tg2[ntg] = beta * (-1.0 + (double) (j - 1) * 2.0
                        / (double) (ntg1 - 1)) * Math.Sqrt(1.0 - t * t) + c2;

                    ntg = ntg + 1;
                }
            }

            i = ntg1;
            j = 1;

            tg1[ntg] = alpha * (-1.0 + (double) (i - 1) * 2.0
                / (double) (ntg1 - 1)) + c1;

            t = -1.0 + (double) (i - 1) * 2.0 / (double) (ntg1 - 1);

            tg2[ntg] = beta * (-1.0 + (double) (j - 1) * 2.0
                / (double) (ntg1 - 1)) * Math.Sqrt(1.0 - t * t) + c2;

            ntg = ntg + 1;

        }
    }
}