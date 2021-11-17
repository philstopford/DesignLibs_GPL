using System;
using Burkardt.PDFLib;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Disk;

public static class MonteCarlo
{
    public static double disk_area(double[] center, double r )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DISK_AREA returns the area of a general disk in 2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 July 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double CENTER[2], coordinates of the center.
        //    This information is not actually needed to compute the area.
        //
        //    Input, double R, the radius of the disk.
        //
        //    Output, double DISK01_AREA, the area of the unit disk.
        //
    {
        double area;
            

        area = Math.PI * r * r;

        return area;
    }

    public static double[] disk_sample(double[] center, double r, int n, ref typeMethods.r8vecNormalData data, ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DISK_SAMPLE uniformly samples the general disk in 2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 July 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double CENTER[2], coordinates of the center.
        //    This information is not actually needed to compute the area.
        //
        //    Input, double R, the radius of the disk.
        //
        //    Input, int N, the number of points.
        //
        //    Input/output, int &SEED, a seed for the random 
        //    number generator.
        //
        //    Output, double X[2*N], the points.
        //
    {
        int i;
        int j;
        double norm;
        double r2;
        double[] v;
        double[] x;

        x = new double[2 * n];

        for (j = 0; j < n; j++)
        {
            v = typeMethods.r8vec_normal_01_new(2, ref data, ref seed);
            //
            //  Compute the length of the vector.
            //
            norm = Math.Sqrt(Math.Pow(v[0], 2) + Math.Pow(v[1], 2));
            //
            //  Normalize the vector.
            //
            for (i = 0; i < 2; i++)
            {
                v[i] /= norm;
            }

            //
            //  Now compute a value to map the point ON the circle INTO the circle.
            //
            r2 = UniformRNG.r8_uniform_01(ref seed);

            for (i = 0; i < 2; i++)
            {
                x[i + j * 2] = center[i] + r * Math.Sqrt(r2) * v[i];
            }
        }

        return x;
    }

    public static double disk01_monomial_integral(int[] e )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DISK01_MONOMIAL_INTEGRAL returns monomial integrals in the unit disk in 2D.
        //
        //  Discussion:
        //
        //    The integration region is 
        //
        //      X^2 + Y^2 <= 1.
        //
        //    The monomial is F(X,Y) = X^E(1) * Y^E(2).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Philip Davis, Philip Rabinowitz,
        //    Methods of Numerical Integration,
        //    Second Edition,
        //    Academic Press, 1984, page 263.
        //
        //  Parameters:
        //
        //    Input, int E[2], the exponents of X and Y in the 
        //    monomial.  Each exponent must be nonnegative.
        //
        //    Output, double DISK01_MONOMIAL_INTEGRAL, the integral.
        //
    {
        double arg;
        int i;
        double integral;
        const double r = 1.0;
        double s;

        if (e[0] < 0 || e[1] < 0)
        {
            Console.WriteLine("");
            Console.WriteLine("DISK01_MONOMIAL_INTEGRAL - Fatal error!");
            Console.WriteLine("  All exponents must be nonnegative.");
            Console.WriteLine("  E[0] = " + e[0] + "");
            Console.WriteLine("  E[1] = " + e[1] + "");
            return 1;
        }

        if (e[0] % 2 == 1 || e[1] % 2 == 1)
        {
            integral = 0.0;
        }
        else
        {
            integral = 2.0;

            for (i = 0; i < 2; i++)
            {
                arg = 0.5 * (e[i] + 1);
                integral *= typeMethods.r8_gamma(arg);
            }

            arg = 0.5 * (e[0] + e[1] + 2);
            integral /= typeMethods.r8_gamma(arg);
        }

        //
        //  Adjust the surface integral to get the volume integral.
        //
        s = e[0] + e[1] + 2;
        integral = integral * Math.Pow(r, s) / s;

        return integral;
    }

}