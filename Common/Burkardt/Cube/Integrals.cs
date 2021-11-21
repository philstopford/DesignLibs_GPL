using System;
using System.Globalization;
using Burkardt.Uniform;

namespace Burkardt.Cube;

public static class Integrals
{
    public static double cube_monomial(double[] a, double[] b, int[] expon)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CUBE_MONOMIAL integrates a monomial over a cube in 3D.
        //
        //  Discussion:
        //
        //    This routine integrates a monomial of the form
        //
        //      product ( 1 <= dim <= 3 ) x(dim)^expon(dim)
        //
        //    The combination 0^0 should be treated as 1.
        //
        //    The integration region is:
        //      A(1) <= X <= B(1)
        //      A(2) <= Y <= B(2)
        //      A(3) <= Z <= B(3)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A[3], B[3], the lower and upper limits.
        //
        //    Input, int EXPON[3], the exponents.
        //
        //    Output, double CUBE_MONOMIAL, the integral of the monomial.
        //
    {
        int i;
        double value = 0;

        for (i = 0; i < 3; i++)
        {
            switch (expon[i])
            {
                case -1:
                    Console.WriteLine("");
                    Console.WriteLine("CUBE_MONOMIAL - Fatal error!");
                    Console.WriteLine("  Exponent of -1 encountered.");
                    return 1;
            }
        }

        value = 1.0;

        for (i = 0; i < 3; i++)
        {
            value = value
                    * (Math.Pow(b[i], expon[i] + 1) - Math.Pow(a[i], expon[i] + 1))
                    / (expon[i] + 1);
        }

        return value;
    }

    public static void cube_monomial_test(int degree_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CUBE_MONOMIAL_TEST tests CUBE_MONOMIAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DEGREE_MAX, the maximum total degree of the
        //    monomials to check.
        //
    {
        double[] a =
            {
                -1.0, -1.0, -1.0
            }
            ;
        int alpha;
        double[] b =
            {
                +1.0, +1.0, +1.0
            }
            ;
        int[] expon = new int[3];

        Console.WriteLine("");
        Console.WriteLine("CUBE_MONOMIAL_TEST");
        Console.WriteLine("  For the unit hexahedron,");
        Console.WriteLine("  CUBE_MONOMIAL returns the exact value of the");
        Console.WriteLine("  integral of X^ALPHA Y^BETA Z^GAMMA");
        Console.WriteLine("");
        Console.WriteLine("  Volume = " + cube_volume(a, b) + "");
        Console.WriteLine("");
        Console.WriteLine("     ALPHA      BETA     GAMMA      INTEGRAL");
        Console.WriteLine("");

        for (alpha = 0; alpha <= degree_max; alpha++)
        {
            expon[0] = alpha;
            int beta;
            for (beta = 0; beta <= degree_max - alpha; beta++)
            {
                expon[1] = beta;
                int gamma;
                for (gamma = 0; gamma <= degree_max - alpha - beta; gamma++)
                {
                    expon[2] = gamma;

                    double value = cube_monomial(a, b, expon);

                    Console.WriteLine("  " + expon[0].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                           + "  " + expon[1].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                           + "  " + expon[2].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                           + "  " + value.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                }
            }
        }
    }

    public static double cube_volume(double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CUBE_VOLUME: volume of a cube in 3D.
        //
        //  Discussion:
        //
        //    The integration region is:
        //      A(1) <= X <= B(1)
        //      A(2) <= Y <= B(2)
        //      A(3) <= Z <= B(3)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A[3], B[3], the lower and upper limits.
        //
        //    Output, double CUBE_VOLUME, the volume.
        //
    {
        int i;
        double value = 0;

        value = 1.0;

        for (i = 0; i < 3; i++)
        {
            value *= (b[i] - a[i]);
        }

        return value;
    }

    public static double cube01_monomial_integral(int[] e)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CUBE01_MONOMIAL_INTEGRAL: monomial integrals on the unit cube in 3D.
        //
        //  Discussion:
        //
        //    The integration region is 
        //
        //      0 <= X <= 1,
        //      0 <= Y <= 1,
        //      0 <= Z <= 1.
        //
        //    The monomial is F(X,Y,Z) = X^E(1) * Y^E(2) * Z^E(3).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 January 2014
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
        //    Input, int E[3], the exponents.  
        //    Each exponent must be nonnegative.
        //
        //    Output, double CUBE01_MONOMIAL_INTEGRAL, the integral.
        //
    {
        int i;
        const int m = 3;

        for (i = 0; i < m; i++)
        {
            switch (e[i])
            {
                case < 0:
                    Console.WriteLine("");
                    Console.WriteLine("CUBE01_MONOMIAL_INTEGRAL - Fatal error!");
                    Console.WriteLine("  All exponents must be nonnegative.");
                    Console.WriteLine("  E[" + i + "] = " + e[i] + "");
                    return 1;
            }
        }

        double integral = 1.0;
        for (i = 0; i < m; i++)
        {
            integral /= e[i] + 1;
        }

        return integral;
    }

    public static double[] cube01_sample(int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CUBE01_SAMPLE samples the unit cube in 3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 January 2014
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
        const int m = 3;

        double[] x = UniformRNG.r8mat_uniform_01_new(m, n, ref seed);

        return x;
    }

    public static double cube01_volume()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CUBE01_VOLUME: volume of the unit cube in 3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double CUBE01_VOLUME, the volume.
        //
    {
        const double volume = 1.0;

        return volume;
    }

}