﻿using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Ball;

public static class MonteCarlo
{
    public static double ball01_monomial_integral(int[] e)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BALL01_MONOMIAL_INTEGRAL returns monomial integrals in the unit ball.
        //
        //  Discussion:
        //
        //    The integration region is 
        //
        //      X^2 + Y^2 + Z^2 <= 1.
        //
        //    The monomial is F(X,Y,Z) = X^E(1) * Y^E(2) * Z^E(3).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 January 2014
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
        //    Input, int E[3], the exponents of X, Y and Z in the 
        //    monomial.  Each exponent must be nonnegative.
        //
        //    Output, double BALL01_MONOMIAL_INTEGRAL, the integral.
        //
    {
        double integral;
        const double r = 1.0;

        if (e[0] < 0 || e[1] < 0 || e[2] < 0)
        {
            Console.WriteLine("");
            Console.WriteLine("BALL01_MONOMIAL_INTEGRAL - Fatal error!");
            Console.WriteLine("  All exponents must be nonnegative.");
            Console.WriteLine("  E[0] = " + e[0] + "");
            Console.WriteLine("  E[1] = " + e[1] + "");
            Console.WriteLine("  E[2] = " + e[2] + "");
            return 1;
        }

        switch (e[0])
        {
            case 0 when e[1] == 0 && e[2] == 0:
                integral = 2.0 * Math.Sqrt(Math.PI * Math.PI * Math.PI) / typeMethods.r8_gamma(1.5);
                break;
            default:
            {
                if (e[0] % 2 == 1 ||
                    e[1] % 2 == 1 ||
                    e[2] % 2 == 1)
                {
                    integral = 0.0;
                }
                else
                {
                    integral = 2.0;

                    int i;
                    for (i = 0; i < 3; i++)
                    {
                        integral *= typeMethods.r8_gamma(0.5 * (e[i] + 1));
                    }

                    integral /= typeMethods.r8_gamma(0.5 * (e[0] + e[1] + e[2] + 3));
                }

                break;
            }
        }

        /*
        The surface integral is now adjusted to give the volume integral.
        */
        double s = e[0] + e[1] + e[2] + 3;

        integral = integral * Math.Pow(r, s) / s;

        return integral;
    }

    public static double[] ball01_sample(int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BALL01_SAMPLE uniformly samples points from the unit ball.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Russell Cheng,
        //    Random Variate Generation,
        //    in Handbook of Simulation,
        //    edited by Jerry Banks,
        //    Wiley, 1998, pages 168.
        //
        //    Reuven Rubinstein,
        //    Monte Carlo Optimization, Simulation, and Sensitivity 
        //    of Queueing Networks,
        //    Krieger, 1992,
        //    ISBN: 0894647644,
        //    LC: QA298.R79.
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input/output, int &SEED, a seed for the random 
        //    number generator.
        //
        //    Output, double X[3*N], the points.
        //
    {
        const double exponent = 1.0 / 3.0;
        int j;
        typeMethods.r8vecNormalData data = new();

        double[] x = typeMethods.r8mat_normal_01_new(3, n, ref data, ref seed);

        for (j = 0; j < n; j++)
        {
            //
            //  Compute the length of the vector.
            //
            double norm = Math.Sqrt(Math.Pow(x[0 + j * 3], 2)
                                    + Math.Pow(x[1 + j * 3], 2)
                                    + Math.Pow(x[2 + j * 3], 2));
            //
            //  Normalize the vector.
            //
            int i;
            for (i = 0; i < 3; i++)
            {
                x[i + j * 3] /= norm;
            }

            //
            //  Now compute a value to map the point ON the sphere INTO the sphere.
            //
            double r = UniformRNG.r8_uniform_01(ref seed);

            for (i = 0; i < 3; i++)
            {
                x[i + j * 3] = Math.Pow(r, exponent) * x[i + j * 3];
            }
        }

        return x;
    }

    public static double ball01_volume()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BALL01_VOLUME returns the volume of the unit ball.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double BALL01_VOLUME, the volume of the unit ball.
        //
    {
        const double r = 1.0;

        double volume = 4.0 * Math.PI * Math.Pow(r, 3) / 3.0;

        return volume;
    }

}