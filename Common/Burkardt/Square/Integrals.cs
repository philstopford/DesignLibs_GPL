using System;
using Burkardt.Uniform;

namespace Burkardt.Square;

public static class Integrals
{
    public static double square01_area ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SQUARE01_AREA: area of the unit square in 2D.
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
        //    Output, double SQUARE01_AREA, the area
        //
    {
        double area;

        area = 1.0;

        return area;
    }
        
    public static double square01_monomial_integral(int[] e)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SQUARE01_MONOMIAL_INTEGRAL: monomial integrals on the unit square in 2D.
        //
        //  Discussion:
        //
        //    The integration region is 
        //
        //      0 <= X <= 1,
        //      0 <= Y <= 1.
        //
        //    The monomial is F(X,Y) = X^E(1) * Y^E(2).
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
        //    Input, int E[2], the exponents.  
        //    Each exponent must be nonnegative.
        //
        //    Output, double SQUARE01_MONOMIAL_INTEGRAL, the integral.
        //
    {
        int i;
        double integral;
        const int m = 2;

        for (i = 0; i < m; i++)
        {
            switch (e[i])
            {
                case < 0:
                    Console.WriteLine("");
                    Console.WriteLine("SQUARE01_MONOMIAL_INTEGRAL - Fatal error!");
                    Console.WriteLine("  All exponents must be nonnegative.");
                    Console.WriteLine("  E[" + i + "] = " + e[i] + "");
                    return 1;
            }
        }

        integral = 1.0;
        for (i = 0; i < m; i++)
        {
            integral /= e[i] + 1;
        }

        return integral;
    }

    public static double[] square01_sample(int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SQUARE01_SAMPLE samples the unit square in 2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 February 2018
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
        //    Output, double SQUARE01_SAMPLE[2*N], the points.
        //
    {
        double a;
        double b;
        const int m = 2;
        double[] x;

        a = 0.0;
        b = 1.0;
        x = UniformRNG.r8mat_uniform_ab_new(m, n, a, b, ref seed);

        return x;
    }

    public static double squaresym_area()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SQUARE01_AREA: area of the symmetric unit square in 2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 February 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double SQUARESYM_AREA, the area
        //
    {
        double area;

        area = 4.0;

        return area;
    }

    public static double squaresym_monomial_integral(int[] e)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SQUARESYM_MONOMIAL_INTEGRAL: monomial integrals on the symmetric unit square in 2D.
        //
        //  Discussion:
        //
        //    The integration region is 
        //
        //      -1 <= X <= 1,
        //      -1 <= Y <= 1.
        //
        //    The monomial is F(X,Y) = X^E(1) * Y^E(2).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 February 2018
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
        //    Input, int E[2], the exponents.  
        //    Each exponent must be nonnegative.
        //
        //    Output, double SQUARESYM_MONOMIAL_INTEGRAL, the integral.
        //
    {
        int i;
        double integral;
        const int m = 2;

        for (i = 0; i < m; i++)
        {
            switch (e[i])
            {
                case < 0:
                    Console.WriteLine("");
                    Console.WriteLine("SQUARESYM_MONOMIAL_INTEGRAL - Fatal error!");
                    Console.WriteLine("  All exponents must be nonnegative.");
                    Console.WriteLine("  E[" + i + "] = " + e[i] + "");
                    return 1;
            }
        }

        if (e[0] % 2 == 1 || e[1] % 2 == 1)
        {
            integral = 0.0;
        }
        else
        {
            integral = 4.0
                       / (e[0] + 1)
                       / (e[1] + 1);
        }

        return integral;
    }

    public static double[] squaresym_sample(int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SQUARESYM_SAMPLE samples the symmetric unit square in 2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 February 2018
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
        //    Output, double SQUARESYM_SAMPLE[2*N], the points.
        //
    {
        double a;
        double b;
        const int m = 2;
        double[] x;

        a = -1.0;
        b = 1.0;
        x = UniformRNG.r8mat_uniform_ab_new(m, n, a, b, ref seed);

        return x;
    }
}