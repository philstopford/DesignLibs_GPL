using System;
using Burkardt.BLAS;

namespace Burkardt.Interpolation;

public static class Padua
{
    public static void padua2(int deg, int degmax, int npd, double[] wpd, double[] fpd,
            double[] raux1, double[] raux2, ref double[] c0, ref double esterr)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PADUA2 computes the Padua interpolation coefficient matrix.
        //
        //  Discussion:
        //
        //    This function computes the coefficient matrix C0, in the 
        //    orthonormal Chebyshev basis T_j(x)T_{k-j}(y), 0 <= j <= k <= DEG, 
        //    T_0(x)=1, T_j(x) = sqrt(2) * cos(j * acos(x)), of the 
        //    interpolation polynomial of degree DEG of the function values FPD 
        //    at the set of NPD Padua points (PD1,PD2) in the square [-1,1]^2. 
        //
        //    The interpolant may be evaluated at an arbitrary point by the 
        //    function PD2VAL. PD1, PD2 and WPD are the Padua points and weights 
        //    computed by PDPTS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //  
        //  Modified:
        //
        //    15 February 2014
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
        //    Input, int DEG, the degree of approximation.
        //
        //    Input, int DEGMAX, the maximum degree allowed.
        //
        //    Input, int NPD, the number of Padua points.
        //
        //    Input, double WPD[NPD], the weights.
        //
        //    Input, double FPD[NPD], the value at the Padua points
        //    of the function to be interpolated.
        //
        //    Workspace, double RAUX1[(DEGMAX+1)*(DEG+2)].
        //
        //    Workspace, double RAUX2[(DEGMAX+1)*(DEG+2)].
        //
        //    Output, double C0[(DEGMAX+2)*(DEG+1)], the coefficient matrix.
        //
        //    Output, double &ESTERR, the estimated error.
        //
    {
        double angle;
        int i;
        int j;
        int k;
            
        double pt;
        //
        //  Build the matrix P_2 and store it in RAUX2.
        //
        for (i = 0; i <= deg + 1; i++)
        {
            angle = i * Math.PI / (deg + 1);
            pt = -Math.Cos(angle);
            PolynomialNS.Chebyshev.cheb(deg, pt, ref raux2, +i * (degmax + 1));
        }

        //
        //  Build the matrix G(f) and store it in C0.
        //
        for (j = 0; j <= deg + 1; j++)
        {
            for (i = 0; i <= degmax + 1; i++)
            {
                c0[i + j * (degmax + 2)] = 0.0;
            }
        }

        k = 0;
        for (j = 0; j <= deg + 1; j++)
        {
            for (i = 0; i <= deg; i++)
            {
                switch ((i + j) % 2)
                {
                    case 0:
                        c0[i + j * (degmax + 2)] = fpd[k] * wpd[k];
                        k += 1;
                        break;
                    default:
                        c0[i + j * (degmax + 2)] = 0.0;
                        break;
                }
            }
        }

        //
        //  Compute the matrix-matrix product G(f)*P_2' and store it in RAUX1.
        //
        BLAS3D.dgemm('n', 't', deg + 1, deg + 1, deg + 2, 1.0,
            c0, degmax + 2, raux2, degmax + 1, 0.0, ref raux1, degmax + 1);
        //
        //  Build the matrix P_1 and store it in RAUX2.
        //
        for (i = 0; i <= deg; i++)
        {
            angle = i * Math.PI / deg;
            pt = -Math.Cos(angle);
            PolynomialNS.Chebyshev.cheb(deg, pt, ref raux2, +i * (degmax + 1));
        }

        //
        //  Compute the matrix-matrix product C(f) = P_1 * ( G(f) * P_2' ) 
        //  and store it in C0.
        //
        BLAS3D.dgemm('n', 'n', deg + 1, deg + 1, deg + 1, 1.0,
            raux2, degmax + 1, raux1, degmax + 1, 0.0, ref c0, degmax + 2);

        c0[deg + 0 * (degmax + 2)] /= 2.0;
        //
        //  Estimate the error.
        //
        esterr = 0.0;
        for (j = 0; j <= 2; j++)
        {
            for (i = 0; i <= deg - j; i++)
            {
                esterr += Math.Abs(c0[i + (deg - i - j) * (degmax + 2)]);
            }
        }

        esterr = 2.0 * esterr;
    }

    public static double pd2val(int deg, int degmax, double[] c0, double tg1, double tg2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PD2VAL evaluates the Padua2 interpolant.
        //
        //  Discussion:
        //
        //    This function returns the value of the interpolant at (TG1,TG2).
        //    C0 is the matrix of the coefficients computed by PADUA2.
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
        //    Input, int DEG, the degree of approximation.
        //
        //    Input, int DEGMAX, the maximum degree allowed.         
        //
        //    Input, double C0[(0:DEGMAX+1)*(0:DEG)], the coefficient matrix.
        //
        //    Input, double TG1, TG2, the first and second coordinates of
        //    the target point.
        //
        //    Output, double PD2VAL, the value of the interpolant at
        //    the target point.
        //
    {
        int i;
        int j;
        double t;
        double[] ttg1;
        double[] ttg2;
        double value = 0;
        //
        //  Compute the normalized Chebyshev polynomials at the target point.
        //
        ttg1 = new double[deg + 1];
        PolynomialNS.Chebyshev.cheb(deg, tg1, ref ttg1);

        ttg2 = new double[deg + 1];
        PolynomialNS.Chebyshev.cheb(deg, tg2, ref ttg2);
        //
        //  Evaluate the interpolant
        //
        value = 0.0;
        for (i = deg; 0 <= i; i--)
        {
            t = 0.0;
            for (j = 0; j <= i; j++)
            {
                t += ttg1[j] * c0[j + (deg - i) * (degmax + 2)];
            }

            value += ttg2[deg - i] * t;
        }

        return value;
    }

    public static void pdpts(int deg, ref double[] pd1, ref double[] pd2, ref double[] wpd, ref int npd)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PDPTS returns the points and weights for Padua interpolation.
        //
        //  Discussion:
        //
        //    This subroutine computes the first family of Padua points and 
        //    weights corresponding to degree DEG.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //  
        //  Modified:
        //
        //    14 February 2014
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
        //    Input, int DEG, the degree of approximation.
        //
        //    Output, double PD1[NPD], PD2[NPD], the first and second
        //    coordinates of the Padua points
        //
        //    Output, double WPD[NPD], the weights.
        //
        //    Output, int &NPD, the number of Padua points.
        //    NPD = ( DEG + 1 ) * ( DEG + 2 ) / 2.
        //
    {
        int itemp0;
        int j;
        int k;
            
        double rtemp0;
        switch (deg)
        {
            //
            //  Compute the Padua points of the first family at degree DEG.
            //
            case 0:
                pd1[0] = -1.0;
                pd2[0] = -1.0;
                wpd[0] = 2.0;
                npd = 1;
                return;
        }

        npd = 0;
        itemp0 = deg * (deg + 1);
        rtemp0 = Math.PI / itemp0;

        for (j = 0; j <= deg + 1; j++)
        {
            for (k = j % 2; k <= deg; k += 2)
            {
                pd1[npd] = -Math.Cos((deg + 1) * k * rtemp0);
                pd2[npd] = -Math.Cos(deg * j * rtemp0);
                wpd[npd] = 2.0 / itemp0;

                if (k == 0 || k == deg)
                {
                    wpd[npd] /= 2.0;
                }

                if (j == 0 || j == deg + 1)
                {
                    wpd[npd] /= 2.0;
                }

                npd += 1;
            }
        }

    }
}