using Burkardt.Types;

namespace Burkardt.IntegralNS
{
    public static partial class Integral
    {
        public static double jacobi_integral ( int expon, double alpha, double beta )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    JACOBI_INTEGRAL evaluates the integral of a monomial with Jacobi weight.
            //
            //  Discussion:
            //
            //    VALUE = Integral ( -1 <= X <= +1 ) x^EXPON (1-x)^ALPHA (1+x)^BETA dx
            //
            //  Modified:
            //
            //    08 September 2007
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
            //    Output, double JACOBI_INTEGRAL, the value of the integral.
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

            c = ( double ) ( expon );

            if ( ( expon % 2 ) == 0 )
            {
                s = +1.0;
            }
            else
            {
                s = -1.0;
            }

            arg1 = - alpha;
            arg2 =   1.0 + c;
            arg3 =   2.0 + beta + c;
            arg4 = - 1.0;

            value1 = typeMethods.r8_hyper_2f1 ( arg1, arg2, arg3, arg4 );

            arg1 = - beta;
            arg2 =   1.0 + c;
            arg3 =   2.0 + alpha + c;
            arg4 = - 1.0;

            value2 = typeMethods.r8_hyper_2f1 ( arg1, arg2, arg3, arg4 );

            value = typeMethods.r8_gamma ( 1.0 + c ) * ( 
                s * typeMethods.r8_gamma ( 1.0 + beta  ) * value1 
                / typeMethods.r8_gamma ( 2.0 + beta  + c ) 
                +     typeMethods.r8_gamma ( 1.0 + alpha ) * value2 
                / typeMethods.r8_gamma ( 2.0 + alpha + c ) );

            return value;
        }
    }
}