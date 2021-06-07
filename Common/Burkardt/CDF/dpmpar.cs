using System;

namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        public static double dpmpar ( int i )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DPMPAR provides machine constants for double precision arithmetic.
            //
            //  Discussion:
            //
            //     DPMPAR PROVIDES THE double PRECISION MACHINE CONSTANTS FOR
            //     THE COMPUTER BEING USED.   It is assumed that tHE
            //     double PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND
            //     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 February 2021
            //
            //  Author:
            //
            //    ALFRED H. MORRIS, JR
            //
            //  Input:
            //
            //    int I, a value of 1, 2 or 3.
            //
            //  Output:
            //
            //    double DPMPAR, the value requested by I.
            //    DPMPAR(1) = B^(1 - M), THE MACHINE PRECISION,
            //    DPMPAR(2) = B^(EMIN - 1), THE SMALLEST MAGNITUDE,
            //    DPMPAR(3) = B^EMAX*(1 - B^(-M)), THE LARGEST MAGNITUDE.
        {
            double b;
            double binv;
            double bm1;
            int emax;
            int emin;
            int ibeta;
            int K1 = 4;
            int K2 = 8;
            int K3 = 9;
            int K4 = 10;
            int m;
            double one;
            double value;
            double w;
            double z;

            if ( i == 1 )
            {
                b = ipmpar ( K1 );
                m = ipmpar ( K2 );
                value = Math.Pow ( b, ( double ) ( 1 - m ) );
            }
            else if ( i == 2 )
            {
                b = ipmpar(K1);
                emin = ipmpar(K3);
                one = 1.0;
                binv = one/b;
                w = Math.Pow(b,(double)(emin+2));
                value = w*binv*binv*binv;
            }
            else 
            {
                ibeta = ipmpar(K1);
                m = ipmpar(K2);
                emax = ipmpar(K4);
                b = ibeta;
                bm1 = ibeta-1;
                one = 1.0;
                z = Math.Pow(b,(double)(m-1));
                w = ((z-one)*b+bm1)/(b*z);
                z = Math.Pow(b,(double)(emax-2));
                value = w*z*b*b;
            }

            return value;
        }
    }
}
