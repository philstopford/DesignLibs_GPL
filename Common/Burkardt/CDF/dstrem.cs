using System;

namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static double dstrem(double z)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DSTREM computes the Sterling remainder ln ( Gamma ( Z ) ) - Sterling ( Z ).
        //
        //  Discussion:
        //
        //    This routine returns
        //
        //      ln ( Gamma ( Z ) ) - Sterling ( Z )
        //
        //    where Sterling(Z) is Sterling's approximation to ln ( Gamma ( Z ) ).
        //
        //    Sterling(Z) = ln ( sqrt ( 2 * PI ) ) + ( Z - 0.5 ) * ln ( Z ) - Z
        //
        //    If 6 <= Z, the routine uses 9 terms of a series in Bernoulli numbers,
        //    with values calculated using Maple.
        //
        //    Otherwise, the difference is computed explicitly.
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
        //    Barry Brown, James Lovato, Kathy Russell.
        //
        //  Parameters:
        //
        //    Input, double *Z, the value at which the Sterling
        //    remainder is to be calculated.  Z must be positive.
        //
        //    Output, double DSTREM, the Sterling remainder.
        //
    {
        const double hln2pi = 0.91893853320467274178e0;

        double[] coef =  {
                0.0e0,0.0833333333333333333333333333333e0,
                -0.00277777777777777777777777777778e0,0.000793650793650793650793650793651e0,
                -0.000595238095238095238095238095238e0,
                0.000841750841750841750841750841751e0,-0.00191752691752691752691752691753e0,
                0.00641025641025641025641025641026e0,-0.0295506535947712418300653594771e0,
                0.179644372368830573164938490016e0
            }
            ;
        const int K1 = 10;
        double dstrem;
        switch (z)
        {
            //
            //    For information, here are the next 11 coefficients of the
            //    remainder term in Sterling's formula
            //            -1.39243221690590111642743221691
            //            13.4028640441683919944789510007
            //            -156.848284626002017306365132452
            //            2193.10333333333333333333333333
            //            -36108.7712537249893571732652192
            //            691472.268851313067108395250776
            //            -0.152382215394074161922833649589D8
            //            0.382900751391414141414141414141D9
            //            -0.108822660357843910890151491655D11
            //            0.347320283765002252252252252252D12
            //            -0.123696021422692744542517103493D14
            //
            case <= 0.0e0:
                throw new Exception("Zero or negative argument in DSTREM");
        }

        switch (z > 6.0e0)
        {
            case false:
                goto S10;
        }
        double T2 = 1.0e0 / Math.Pow(z, 2.0);
        dstrem = eval_pol(coef, K1, T2) * z;
        goto S20;
        S10:
        double sterl = hln2pi + (z - 0.5e0) * Math.Log(z) - z;
        dstrem = gamma_log(z) - sterl;
        S20:
        return dstrem;
    }
}