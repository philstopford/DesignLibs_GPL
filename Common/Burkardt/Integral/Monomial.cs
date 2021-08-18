using System;
using Burkardt.Types;

namespace Burkardt.IntegralNS
{
    public static class Monomial
    {
        public static double monomial_int01 ( int dim_num, int[] expon )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MONOMIAL_INT01 returns the integral of a monomial over the [0,1] hypercube.
            //
            //  Discussion:
            //
            //    This routine evaluates a monomial of the form
            //
            //      product ( 1 <= dim <= dim_num ) x(dim)^expon(dim)
            //
            //    where the exponents are nonnegative integers.  Note that
            //    if the combination 0^0 is encountered, it should be treated
            //    as 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 November 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, int EXPON[DIM_NUM], the exponents.
            //
            //    Output, double MONOMIAL_INT01, the value of the integral of the
            //    monomial.
            //
        {
            int dim;
            double value;

            value = 1.0;
            for ( dim = 0; dim < dim_num; dim++ )
            {
                value = value / ( double ) ( expon[dim] + 1 );
            }

            return value;
        }
        public static double monomial_integral_generalized_hermite(int expon, double alpha)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MONOMIAL_INTEGRAL_GENERALIZED_HERMITE evaluates a 1D monomial generalized Hermite integral.
            //
            //  Discussion:
            //
            //    The integral being computed is
            //
            //      integral ( -oo < x < +oo ) x^n |x|^alpha exp(-x*x) dx 
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 October 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int EXPON, the exponent of the monomial.
            //    0 <= EXPON.
            //
            //    Input, double ALPHA, the exponent of |x| in the weight.
            //    -1.0 < ALPHA.
            //
            //    Output, double MONOMIAL_INTEGRAL_GENERALIZED_HERMITE, 
            //    the value of the integral.
            //
        {
            double arg;
            double value;

            if ((expon % 2) == 1)
            {
                value = 0.0;
            }
            else
            {
                arg = alpha + (double) (expon);

                if (arg <= -1.0)
                {
                    value = -typeMethods.r8_huge();
                }
                else
                {
                    arg = (arg + 1.0) / 2.0;
                    value = typeMethods.r8_gamma(arg);
                }
            }

            return value;
        }

        public static double monomial_integral_generalized_laguerre(int expon, double alpha)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MONOMIAL_INTEGRAL_GENERALIZED_LAGUERRE evaluates a 1D monomial generalized Laguerre integral.
            //
            //  Discussion:
            //
            //    The integral being computed is
            //
            //      integral ( 0 <= x < +oo ) x^n x^alpha exp(-x) dx 
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    10 October 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int EXPON, the exponent of the monomial.
            //    0 <= EXPON.
            //
            //    Input, double ALPHA, the exponent of x in the weight.
            //    -1.0 < ALPHA.
            //
            //    Output, double MONOMIAL_INTEGRAL_GENERALIZED_LAGUERRE, 
            //    the value of the integral.
            //
        {
            double arg;
            double value;

            arg = alpha + (double) (expon + 1);

            value = typeMethods.r8_gamma(arg);

            return value;
        }

        public static double monomial_integral_hermite(int expon)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MONOMIAL_INTEGRAL_HERMITE evaluates a 1D monomial Hermite integral.
            //
            //  Discussion:
            //
            //    H(n) = Integral ( -oo < x < +oo ) x^n exp(-x*x) dx
            //
            //    H(n) is 0 for n odd.
            //
            //    H(n) = (n-1)!! * sqrt(pi) / 2^(n/2) for n even.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 August 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int EXPON, the exponent.  
            //    0 <= EXPON.
            //
            //    Output, double MONOMIAL_INTEGRAL_HERMITE, 
            //    the value of the integral.
            //
        {
            double pi = 3.141592653589793;
            double value;

            if (expon < 0)
            {
                value = -typeMethods.r8_huge();
            }
            else if ((expon % 2) == 1)
            {
                value = 0.0;
            }
            else
            {
                value = typeMethods.r8_factorial2(expon - 1) * Math.Sqrt(pi) / Math.Pow(2.0, expon / 2);
            }

            return value;
        }

        public static double monomial_integral_jacobi(int expon, double alpha, double beta)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MONOMIAL_INTEGRAL_JACOBI evaluates the integral of a monomial with Jacobi weight.
            //
            //  Discussion:
            //
            //    VALUE = Integral ( -1 <= X <= +1 ) x^EXPON (1-x)^ALPHA (1+x)^BETA dx
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 October 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int EXPON, the exponent.
            //
            //    Input, double ALPHA, the exponent of (1-X) in the weight factor.
            //
            //    Input, double BETA, the exponent of (1+X) in the weight factor.
            //
            //    Output, double MONOMIAL_INTEGRAL_JACOBI, the value of the integral.
            //
        {
            double arg1;
            double arg2;
            double arg3;
            double arg4;
            double c;
            double s;
            double value;
            double value1;
            double value2;

            c = (double) (expon);

            if ((expon % 2) == 0)
            {
                s = +1.0;
            }
            else
            {
                s = -1.0;
            }

            arg1 = -alpha;
            arg2 = 1.0 + c;
            arg3 = 2.0 + beta + c;
            arg4 = -1.0;

            value1 = typeMethods.r8_hyper_2f1(arg1, arg2, arg3, arg4);

            arg1 = -beta;
            arg2 = 1.0 + c;
            arg3 = 2.0 + alpha + c;
            arg4 = -1.0;

            value2 = typeMethods.r8_hyper_2f1(arg1, arg2, arg3, arg4);

            value = typeMethods.r8_gamma(1.0 + c) * (s * typeMethods.r8_gamma(1.0 + beta) * value1
                                                     / typeMethods.r8_gamma(2.0 + beta + c)
                                                     + typeMethods.r8_gamma(1.0 + alpha) * value2
                                                     / typeMethods.r8_gamma(2.0 + alpha + c));

            return value;
        }

        public static double monomial_integral_laguerre(int expon)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MONOMIAL_INTEGRAL_LAGUERRE evaluates a 1D monomial Laguerre integral.
            //
            //  Discussion:
            //
            //    The integral being computed is
            //
            //      integral ( 0 <= x < +oo ) x^n * exp ( -x ) dx
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 August 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int EXPON, the exponent.
            //    0 <= EXPON.
            //
            //    Output, double MONOMIAL_INTEGRAL_LAGUERRE, 
            //    the value of the integral.
            //
        {
            double value;

            value = typeMethods.r8_factorial(expon);

            return value;
        }

        public static double monomial_integral_legendre(int expon)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MONOMIAL_INTEGRAL_LEGENDRE evaluates a 1D monomial Legendre integral.
            //
            //  Discussion:
            //
            //    The integral being computed is
            //
            //      integral ( -1 <= x < +1 ) x^n dx 
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 August 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int EXPON, the exponent.
            //    0 <= EXPON.
            //
            //    Output, double MONOMIAL_INTEGRAL_LEGENDRE, 
            //    the value of the integral.
            //
        {
            double value = 0;

            if ((expon % 2) == 1)
            {
                value = 0.0;
            }
            else if ((expon % 2) == 0)
            {
                value = 2.0 / (double) (expon + 1);
            }

            return value;
        }

        public static double monomial_integral_mixed(int dim_num, int[] rule, double[] alpha,
                double[] beta, int[] expon)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MONOMIAL_INTEGRAL_MIXED evaluates a multi-D monomial mixed integral.
            //
            //  Discussion:
            //
            //    This routine evaluates a monomial of the form
            //
            //      product ( 1 <= dim <= dim_num ) x(dim)^expon(dim)
            //
            //    where the exponents are nonnegative integers.  Note that
            //    if the combination 0^0 is encountered, it should be treated
            //    as 1.
            //
            //    The integration is carried out in a region that is a direct product
            //    of 1D factors that may be of Legendre, Laguerre or Hermite type,
            //    and the integration includes the weight functions associated with
            //    the 1D factors.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 August 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, int RULE[DIM_NUM], the component rules.
            //    1, Gauss-Legendre rule on [-1,+1];
            //    2, Gauss-Jacobi rule on [-1,+1];
            //    3, Gauss-Laguerre rule on [0,+oo);
            //    4, Generalized Gauss-Laguerre rule on [0,+oo);
            //    5, Gauss-Hermite rule on (-oo,+oo);
            //    6, Generalized Gauss-Hermite rule on (-oo,+oo).
            //
            //    Input, double ALPHA[DIM_NUM], BETA[DIM_NUM], parameters that
            //    may be needed for Jacobi, Generalized-Laguerre, or Generalized Hermite rules.
            //
            //    Input, int EXPON[DIM_NUM], the exponents.
            //
            //    Output, double MONOMIAL_INTEGRAL_MIXED, 
            //    the value of the integral.
            //
        {
            int dim;
            double value;

            value = 1.0;

            for (dim = 0; dim < dim_num; dim++)
            {
                if (rule[dim] == 1)
                {
                    value = value * monomial_integral_legendre(expon[dim]);
                }
                else if (rule[dim] == 2)
                {
                    value = value
                            * monomial_integral_jacobi(expon[dim], alpha[dim], beta[dim]);
                }
                else if (rule[dim] == 3)
                {
                    value = value * monomial_integral_laguerre(expon[dim]);
                }
                else if (rule[dim] == 4)
                {
                    value = value
                            * monomial_integral_generalized_laguerre(expon[dim], alpha[dim]);
                }
                else if (rule[dim] == 5)
                {
                    value = value * monomial_integral_hermite(expon[dim]);
                }
                else if (rule[dim] == 6)
                {
                    value = value
                            * monomial_integral_generalized_hermite(expon[dim], alpha[dim]);
                }
            }

            return value;
        }
    }
}