using Burkardt.Types;

namespace Burkardt.IntegralNS;

public static partial class Integral
{
    public static double gen_laguerre_integral ( int expon, double alpha )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEN_LAGUERRE_INTEGRAL evaluates a monomial generalized Laguerre integral.
        //
        //  Discussion:
        //
        //    L(n,alpha) = Integral ( 0 <= x < +oo ) x^n * x^alpha exp(-x) dx
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 February 2008
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
        //    Input, double ALPHA, the exponent of X in the weight function.
        //    -1.0 < ALPHA.
        //
        //    Output, double GEN_LAGUERRE_INTEGRAL, the value of the integral.
        //
    {
        double arg;
        double value = 0;

        arg = alpha + (expon + 1.0);
        value = typeMethods.r8_gamma ( arg );

        return value;
    }
    public static double laguerre_integral ( int p )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGUERRE_INTEGRAL evaluates a monomial Laguerre integral.
        //
        //  Discussion:
        //
        //    The integral being computed is
        //
        //      integral ( 0 <= x < +oo ) x^p * exp ( -x ) dx
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int P, the exponent.
        //    0 <= P.
        //
        //    Output, double LAGUERRE_INTEGRAL, the value of the integral.
        //
    {
        double s;

        s = typeMethods.r8_factorial ( p );

        return s;
    }
        
    public static double laguerre_integral_nd ( int dim_num, int[] expon )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGUERRE_INTEGRAL_ND evaluates a multidimensional Laguerre polynomial integral.
        //
        //  Discussion:
        //
        //    L(1,n) = Integral ( 0 <= x < Infinity ) x^n exp ( -x ) dx
        //           = n!
        //
        //    L(d,n) = Integral ( 0 <= x(i) < Infinity ) 
        //             x1^n1 * x2^n2...*xd^nd * exp(-x1-x2...-xd ) dx
        //           = Product ( n(i)! ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 October 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        // Parameters:
        //
        //    Input, int DIM_NUM, the dimension of the integral.
        //
        //    Input, int EXPON[DIM_NUM], the order of the integral.  
        //    0 <= EXPON(1:DIM_NUM).
        //
        //    Output, double LAGUERRE_INTEGRAL, the value of the integral.
        //
    { 
        int dim;
        double value = 0;

        value = 1.0;
        for ( dim = 0; dim < dim_num; dim++ )
        {
            value *= typeMethods.r8_factorial ( expon[dim] );
        }

        return value;
    }
}