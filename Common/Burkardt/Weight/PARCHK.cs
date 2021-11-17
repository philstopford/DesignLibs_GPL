using System;

namespace Burkardt.Weight;

public static class PARCHK
{
    public static void parchk(int kind, int m, double alpha, double beta)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARCHK checks parameters ALPHA and BETA for classical weight functions. 
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 January 2010
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Sylvan Elhay, Jaroslav Kautsky,
        //    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
        //    Interpolatory Quadrature,
        //    ACM Transactions on Mathematical Software,
        //    Volume 13, Number 4, December 1987, pages 399-415.
        //
        //  Parameters:
        //
        //    Input, int KIND, the rule.
        //    1, Legendre,             (a,b)       1.0
        //    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
        //    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
        //    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
        //    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
        //    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
        //    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
        //    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
        //
        //    Input, int M, the order of the highest moment to
        //    be calculated.  This value is only needed when KIND = 8.
        //
        //    Input, double ALPHA, BETA, the parameters, if required
        //    by the value of KIND.
        //
    {
        switch (kind)
        {
            case <= 0:
                Console.WriteLine("");
                Console.WriteLine("PARCHK - Fatal error!");
                Console.WriteLine("  KIND <= 0.");
                return;
            //
            //  Check ALPHA for Gegenbauer, Jacobi, Laguerre, Hermite, Exponential.
            //
            case >= 3 when alpha <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("PARCHK - Fatal error!");
                Console.WriteLine("  3 <= KIND and ALPHA <= -1.");
                return;
            //
            //  Check BETA for Jacobi.
            //
            case 4 when beta <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("PARCHK - Fatal error!");
                Console.WriteLine("  KIND == 4 and BETA <= -1.0.");
                return;
            //
            //  Check ALPHA and BETA for rational.
            //
            case 8:
            {
                double tmp = alpha + beta + m + 1.0;
                if (0.0 <= tmp || tmp <= beta)
                {
                    Console.WriteLine("");
                    Console.WriteLine("PARCHK - Fatal error!");
                    Console.WriteLine("  KIND == 8 but condition on ALPHA and BETA fails.");
                }

                break;
            }
        }
    }
}