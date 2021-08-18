using System;
using Burkardt.Types;

namespace Burkardt.IntegralNS
{
    public static partial class Integral
    {
        public static double hermite_integral_nd ( int dim_num, int[] expon )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_INTEGRAL_ND evaluates a multidimensional Hermite polynomial integral.
            //
            //  Discussion:
            //
            //    H(d,n) = Integral ( -Infinity < x < Infinity ) 
            //      x1^n1 * x2^n2...*xd^nd * exp(-x1^2-x2^2...-xd^2 ) dx
            //
            //    H(d,n) is 0 if any n(i) odd.
            //
            //    H(d,n) = product ( 1 <= i <= d ) 
            //      ( (n(i)-1)!! * sqrt(pi) / 2^(n(i)/2) for all n(i) even.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    03 October 2007
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
            //    Output, double HERMITE_INTEGRAL, the value of the integral.
            //
        { 
            int dim;
            double pi = 3.141592653589793;
            double value;

            for ( dim = 0; dim < dim_num; dim++ )
            {
                if ( expon[dim] < 0 )
                {
                    value = - typeMethods.r8_huge ( );
                    return value;
                }
            }

            for ( dim = 0; dim < dim_num; dim++ )
            {
                if ( ( expon[dim] % 2 ) == 1 )
                {
                    value = 0.0;
                    return value;
                }
            }

            value = 1.0;
            for ( dim = 0; dim < dim_num; dim++ )
            {
                value = value * typeMethods.r8_factorial2 ( expon[dim] - 1 ) * Math.Sqrt ( pi ) 
                        / ( double ) Math.Pow ( 2, expon[dim] / 2 );
            }

            return value;
        }
        
        public static double hermite_integral ( int p )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_INTEGRAL evaluates a monomial Hermite integral.
            //
            //  Discussion:
            //
            //    Integral ( -oo < x < oo ) x^p exp(-x^2) dx
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
            //    Input, int P, the exponent of the monomial.  
            //    0 <= P.
            //
            //    Output, double HERMITE_INTEGRAL, the value of the integral.
            //
        {
            
            double value;

            if ( ( p % 2 ) == 0 )
            {
                value = typeMethods.r8_factorial2 ( p - 1 ) * Math.Sqrt ( Math.PI ) / Math.Pow ( 2.0, p / 2 );
            }
            else
            {
                value = 0.0;
            }
            return value;
        }
        
        public static double h_integral ( int n )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    H_INTEGRAL evaluates the integral of H(i,x).
            //
            //  Discussion:
            //
            //    H(i,x) is the physicist's Hermite polynomial of degree I.
            //
            //    The integral computed is:
            //
            //      integral ( -oo < x < +oo ) H(i,x) exp(-x^2) dx
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 March 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the integral.  
            //    0 <= N.
            //
            //    Output, double H_INTEGRAL, the value of the integral.
            //
        {
            const double r8_pi = 3.141592653589793;
            double value;

            if ( ( n % 2 ) == 1 )
            {
                value = 0.0;
            }
            else
            {
                value = typeMethods.r8_factorial2 ( n - 1 ) * Math.Sqrt ( r8_pi ) / Math.Pow ( 2.0, n / 2 );
            }

            return value;
        }
        
        public static double he_double_product_integral ( int i, int j )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HE_DOUBLE_PRODUCT_INTEGRAL: integral of He(i,x)*He(j,x)*e^(-x^2/2).
            //
            //  Discussion:
            //
            //    He(i,x) represents the probabilist's Hermite polynomial.
            //
            //    VALUE = integral ( -oo < x < +oo ) He(i,x)*He(j,x) exp(-x^2/2) dx
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 March 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Dongbin Xiu,
            //    Numerical Methods for Stochastic Computations: A Spectral Method Approach,
            //    Princeton, 2010,
            //    ISBN13: 978-0-691-14212-8,
            //    LC: QA274.23.X58.
            //
            //  Parameters:
            //
            //    Input, int I, J, the polynomial indices.
            //
            //    Output, double HE_DOUBLE_PRODUCT_INTEGRAL, the value of the integral.
            //
        {
            double value;

            if ( i == j )
            {
                value = typeMethods.r8_factorial ( i );
            }
            else
            {
                value = 0.0;
            }
            return value;
        }

        public static double he_integral ( int n )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HE_INTEGRAL evaluates the integral of He(i,x).
            //
            //  Discussion:
            //
            //    He(i,x) represents the probabilist's Hermite polynomial.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 March 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the integral.  
            //    0 <= N.
            //
            //    Output, double HE_INTEGRAL, the value of the integral.
            //
        {
            const double r8_pi = 3.141592653589793;
            double value;

            if ( ( n % 2 ) == 1 )
            {
                value = 0.0;
            }
            else
            {
                value = typeMethods.r8_factorial2 ( n - 1 ) * Math.Sqrt ( 2.0 * r8_pi );
            }

            return value;
        }
        
        public static double he_triple_product_integral ( int i, int j, int k )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HE_TRIPLE_PRODUCT_INTEGRAL: integral of He(i,x)*He(j,x)*He(k,x)*e^(-x^2/2).
            //
            //  Discussion:
            //
            //    He(i,x) represents the probabilist's Hermite polynomial.
            //
            //    VALUE = integral ( -oo < x < +oo ) He(i,x)*He(j,x)*He(k,x) exp(-x^2/2) dx
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 March 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Dongbin Xiu,
            //    Numerical Methods for Stochastic Computations: A Spectral Method Approach,
            //    Princeton, 2010,
            //    ISBN13: 978-0-691-14212-8,
            //    LC: QA274.23.X58.
            //
            //  Parameters:
            //
            //    Input, int I, J, K, the polynomial indices.
            //
            //    Output, double HE_TRIPLE_PRODUCT_INTEGRAL, the value of the integral.
            //
        {
            int s;
            double value;

            s = ( i + j + k ) / 2;

            if ( s < i || s < j || s < k )
            {
                value = 0.0;
            }
            else if ( ( ( i + j + k ) % 2 ) != 0 )
            {
                value = 0.0;
            }
            else
            {
                value = typeMethods.r8_factorial ( i ) / typeMethods.r8_factorial ( s - i ) 
                    * typeMethods.r8_factorial ( j ) / typeMethods.r8_factorial ( s - j ) 
                    * typeMethods.r8_factorial ( k ) / typeMethods.r8_factorial ( s - k );
            }

            return value;
        }
    }
}