﻿using System;
using Burkardt.Types;

namespace Burkardt.IntegralNS;

public static partial class Integral
{
    public static double gen_hermite_integral ( int expon, double alpha )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEN_HERMITE_INTEGRAL evaluates a monomial generalized Hermite integral.
        //
        //  Discussion:
        //
        //    H(n,alpha) = Integral ( -oo < x < +oo ) x^n |x|^alpha exp(-x^2) dx
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 February 2008
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
        //    Input, double ALPHA, the exponent of |X| in the weight function.
        //    -1.0 < ALPHA.
        //
        //    Output, double GEN_HERMITE_INTEGRAL, the value of the integral.
        //
    {
        double value;

        switch (expon % 2)
        {
            case 1:
                value = 0.0;
                break;
            default:
            {
                double a = alpha + expon;
                switch (a)
                {
                    case <= - 1.0:
                        value = - typeMethods.r8_huge ( );
                        break;
                    default:
                        double arg = ( a + 1.0 ) / 2.0;
                        value = typeMethods.r8_gamma ( arg );
                        break;
                }

                break;
            }
        }
        return value;
    }
    public static double fn_integral ( int d )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FN_INTEGRAL is the integral of the Hermite test function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 December 2012
        //
        //  Author:
        //
        //    John Burkardt.
        //
        //  Parameters:
        //
        //    Input, int D, the spatial dimension.
        //
        //    Output, double FN_INTEGRAL, the integral value.
        //
    {
        const int exponent = 6;

        double value = typeMethods.i4_factorial2 ( exponent - 1 );

        return value;
    }

    public static double[] fn_value ( int d, int n, double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FN_VALUE is a Hermite test function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 May 2012
        //
        //  Author:
        //
        //    John Burkardt.
        //
        //  Parameters:
        //
        //    Input, int D, the spatial dimension.
        //
        //    Input, int N, the number of points.
        //
        //    Input, double X[D*N], the points.
        //
        //    Output, double FN_VALUE[N], the function values.
        //
    {
        const int exponent = 6;
        int i;

        double[] fx = new double[n];

        for ( i = 0; i < n; i++ )
        {
            fx[i] = Math.Pow ( x[0+i*d], exponent );
        }

        return fx;
    }
        
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
            
        double value = 0;

        for ( dim = 0; dim < dim_num; dim++ )
        {
            switch (expon[dim])
            {
                case < 0:
                    value = - typeMethods.r8_huge ( );
                    return value;
            }
        }

        for ( dim = 0; dim < dim_num; dim++ )
        {
            switch (expon[dim] % 2)
            {
                case 1:
                    value = 0.0;
                    return value;
            }
        }

        value = 1.0;
        for ( dim = 0; dim < dim_num; dim++ )
        {
            value = value * typeMethods.r8_factorial2 ( expon[dim] - 1 ) * Math.Sqrt ( Math.PI ) 
                    / Math.Pow ( 2, (double)expon[dim] / 2 );
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
        double value = (p % 2) switch
        {
            0 => typeMethods.r8_factorial2(p - 1) * Math.Sqrt(Math.PI) / Math.Pow(2.0, (double)p / 2),
            _ => 0.0
        };
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
        double value = (n % 2) switch
        {
            1 => 0.0,
            _ => typeMethods.r8_factorial2(n - 1) * Math.Sqrt(Math.PI) / Math.Pow(2.0, (double)n / 2)
        };

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
        double value = i == j ? typeMethods.r8_factorial ( i ) : 0.0;
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
        double value = (n % 2) switch
        {
            1 => 0.0,
            _ => typeMethods.r8_factorial2(n - 1) * Math.Sqrt(2.0 * Math.PI)
        };

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
        double value = 0;

        int s = ( i + j + k ) / 2;

        if ( s < i || s < j || s < k )
        {
            value = 0.0;
        }
        else if ( ( i + j + k ) % 2 != 0 )
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