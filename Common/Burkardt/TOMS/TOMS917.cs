using System;
using System.Numerics;
using Burkardt.Types;

namespace Burkardt;

public static class WrightOmega
{
    public static Complex wrightomega(Complex z)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WRIGHTOMEGA is the simple routine for evaluating the Wright Omega function. 
        //
        //  Discussion:
        //
        //    This function is called by:
        //
        //      w = wrightomega ( z )
        //
        //    This function makes a call to the more powerful wrightomega_ext() function.
        //
        //  Modified:
        //
        //    14 May 2016
        //
        //  Author:
        //
        //    Piers Lawrence, Robert Corless, David Jeffrey
        //
        //  Reference:
        //
        //    Piers Lawrence, Robert Corless, David Jeffrey,
        //    Algorithm 917: Complex Double-Precision Evaluation of the Wright Omega 
        //    Function,
        //    ACM Transactions on Mathematical Software,
        //    Volume 38, Number 3, Article 20, April 2012, 17 pages.
        //
        //  Parameters:
        //
        //    Input, Complex Z, the argument.
        //
        //    Output, Complex WRIGHTOMEGA, the value of the Wright Omega
        //    function of Z.
        //
    {
        Complex cond = new();
        Complex e = new();
        Complex r = new();
        Complex w = new();

        wrightomega_ext(z, ref w, ref e, ref r, ref cond);

        return w;
    }

    public static int wrightomega_ext(Complex z, ref Complex w,
            ref Complex e, ref Complex r, ref Complex cond)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WRIGHTOMEGA_EXT computes the Wright Omega function with extra information.
        //
        //  Discussion:
        //
        //    WRIGHTOMEGA_EXT is the extended routine for evaluating the Wright
        //    Omega function with the option of extracting the last update step,
        //    the penultimate residual and the condition number estimate.
        //
        //  Modified:
        //
        //    14 May 2016
        //
        //  Author:
        //
        //    Piers Lawrence, Robert Corless, David Jeffrey
        //
        //  Reference:
        //
        //    Piers Lawrence, Robert Corless, David Jeffrey,
        //    Algorithm 917: Complex Double-Precision Evaluation of the Wright Omega 
        //    Function,
        //    ACM Transactions on Mathematical Software,
        //    Volume 38, Number 3, Article 20, April 2012, 17 pages.
        //
        //  Parameters:
        //
        //    Input, Complex Z, value at which to evaluate Wrightomega().
        //
        //    Output, Complex &W, the value of Wrightomega(z).
        //
        //    Output, Complex &E, the last update step in the iterative scheme.
        //
        //    Output, Complex &R, the penultimate residual,
        //    r_k = z - w_k - log(w_k)
        //
        //    Output, Complex &COND, the condition number estimate. 
        //
        //    Output, int WRIGHTOMEGA_EXT, error flag;
        //    0, successful computation.
        //    nonzero, the computation failed.        
        //
    {
        double near;
        Complex pz;
        double s = 1.0;
        Complex t;
        Complex wp1;
        double x;
        double y;
        double ympi;
        double yppi;
        // 
        //  Extract real and imaginary parts of Z. 
        //
        x = z.Real;
        y = z.Imaginary;
        // 
        //  Compute if we are near the branch cuts.
        //
        ympi = y - Math.PI;
        yppi = y + Math.PI;
        near = 0.01;
        // 
        //  Test for floating point exceptions:
        //

        //
        //  NaN output for NaN input.
        //
        if (x == double.NaN || y == double.NaN)
        {
            w = new Complex(0.0 / 0.0, 0.0 / 0.0);
            e = new Complex(0.0, 0.0);
            r = new Complex(0.0, 0.0);
            return 0;
        }
        //
        //  Signed zeros between branches.
        //

        switch (x)
        {
            case (double.NegativeInfinity or double.PositiveInfinity) and < 0.0 when -Math.PI < y && y <= Math.PI:
            {
                w = Math.Abs(y) switch
                {
                    <= Math.PI / 2.0 => +0.0,
                    _ => -0.0
                };

                w += y switch
                {
                    >= 0.0 => new Complex(0.0, 0.0),
                    _ => new Complex(0.0, -1.0 * 0.0)
                };

                e = new Complex(0.0, 0.0);
                r = new Complex(0.0, 0.0);
                return 0;
            }
        }
        //
        //  Asymptotic for large z.
        //
        if (x == double.NegativeInfinity || x == double.PositiveInfinity || y == double.NegativeInfinity || y == double.PositiveInfinity)
        {
            w = new Complex(x, y);
            e = new Complex(0.0, 0.0);
            r = new Complex(0.0, 0.0);
            return 0;
        }

        switch (x)
        {
            //
            //  Test if exactly on the singular points.
            //
            case -1.0 when Complex.Abs(y) == Math.PI:
                w = new Complex(-1.0, 0.0);
                e = new Complex(0.0, 0.0);
                r = new Complex(0.0, 0.0);
                return 0;
            // 
            //  Choose approximation based on region.
            //
            //
            //  Region 1: upper branch point.
            //  Series about z=-1+Pi*I.
            //
            case > -2.0 and <= 1.0 when 1.0 < y && y < 2.0 * Math.PI:
                pz = Complex.Conjugate(Complex.Sqrt(Complex.Conjugate(2.0 * (z + new Complex(1.0, -Math.PI)))));

                w = -1.0
                    + (new Complex(0.0, 1.0)
                       + (1.0 / 3.0
                          + (-1.0 / 36.0 * new Complex(0.0, 1.0)
                             + (1.0 / 270.0 + 1.0 / 4320.0 * new Complex(0.0, 1.0) * pz)
                             * pz) * pz) * pz) * pz;
                break;
            //
            //  Region 2: lower branch point.
            //  Series about z=-1-Pi*I.
            //
            case > -2.0 and <= 1.0 when -2.0 * Math.PI < y && y < -1.0:
                pz = Complex.Conjugate(Complex.Sqrt(Complex.Conjugate(2.0 * (z + 1.0 + new Complex(0.0, Math.PI)))));

                w = -1.0
                    + (-new Complex(0.0, 1.0) + (1.0 / 3.0
                                                 + (1.0 / 36.0 * new Complex(0.0, 1.0)
                                                    + (1.0 / 270.0 - 1.0 / 4320.0 * new Complex(0.0, 1.0) * pz)
                                                    * pz) * pz) * pz) * pz;
                break;
            //
            //  Region 3: between branch cuts.
            //  Series: About -infinity.
            //
            case <= -2.0 when -Math.PI < y && y <= Math.PI:
                pz = Complex.Exp(z);
                w = (1.0
                     + (-1.0
                        + (3.0 / 2.0
                           + (-8.0 / 3.0
                              + 125.0 / 24.0 * pz) * pz) * pz) * pz) * pz;
                break;
            //
            //  Region 4: Mushroom.
            //  Series about z=1.
            //
            case > -2.0 and <= 1.0 when -1.0 <= y && y <= 1.0:
            case > -2.0 when (x - 1.0) * (x - 1.0) + y * y <= Math.PI * Math.PI:
                pz = z - 1.0;
                w = 1.0 / 2.0 + 1.0 / 2.0 * z
                              + (1.0 / 16.0
                                 + (-1.0 / 192.0
                                    + (-1.0 / 3072.0 + 13.0 / 61440.0 * pz) * pz) * pz) * pz * pz;
                break;
            //
            //  Region 5: Top wing.
            //  Negative log series.
            //
            case <= -1.05 when Math.PI < y && y - Math.PI <= -0.75 * (x + 1.0):
                t = z - new Complex(0.0, Math.PI);
                pz = Complex.Log(-t);
                w = ((1.0 + (-3.0 / 2.0 + 1.0 / 3.0 * pz) * pz) * pz
                     + ((-1.0 + 1.0 / 2.0 * pz) * pz + (pz + (-pz + t) * t) * t) * t)
                    / (t * t * t);
                break;
            //
            //  Region 6: Bottom wing.
            //  Negative log series.
            //
            case <= -1.05 when 0.75 * (x + 1.0) < y + Math.PI && y + Math.PI <= 0.0:
                t = z + new Complex(0.0, Math.PI);
                pz = Complex.Log(-t);
                w = ((1.0 + (-3.0 / 2.0 + 1.0 / 3.0 * pz) * pz) * pz
                     + ((-1.0 + 1.0 / 2.0 * pz) * pz + (pz + (-pz + t) * t) * t) * t)
                    / (t * t * t);
                break;
            //
            default:
                pz = Complex.Log(z);
                w = ((1.0 + (-3.0 / 2.0 + 1.0 / 3.0 * pz) * pz) * pz
                     + ((-1.0 + 1.0 / 2.0 * pz) * pz + (pz + (-pz + z) * z) * z) * z)
                    / (z * z * z);
                break;
        }

        //
        //  Regularize if near branch cuts.
        ///
        if (x <= -1.0 + near && (Complex.Abs(ympi) <= near || Complex.Abs(yppi) <= near))
        {
            s = -1.0;
            if (Complex.Abs(ympi) <= near)
            {
                ympi = ympi switch
                {
                    <= 0.0 =>
                        // fesetround ( FE_DOWNWARD );
                        y - Math.PI,
                    //
                    //  Recompute ympi with directed rounding.
                    //
                    // fesetround ( FE_UPWARD );
                    _ => y - Math.PI
                };

                z = new Complex(x, ympi);
                // 
                //  Return rounding to default.
                //
                // fesetround ( FE_TONEAREST );
            }
            else
            {
                yppi = yppi switch
                {
                    <= 0.0 =>
                        // fesetround ( FE_DOWNWARD );
                        y + Math.PI,
                    //
                    //  Recompute yppi with directed rounding.
                    //
                    // fesetround ( FE_UPWARD );
                    _ => y + Math.PI
                };

                z = new Complex(x, yppi);
                // 
                //  Return rounding to default.
                //
                // fesetround ( FE_TONEAREST );
            }
        }

        //
        //  Iteration one.
        //
        w = s * w;
        r = z - s * w - Complex.Log(w);
        wp1 = s * w + 1.0;
        e = r / wp1 * (2.0 * wp1 * (wp1 + 2.0 / 3.0 * r) - r)
            / (2.0 * wp1 * (wp1 + 2.0 / 3.0 * r) - 2.0 * r);
        w *= (1.0 + e);
        //
        //  Iteration two.
        //

        if (Complex.Abs((2.0 * w * w - 8.0 * w - 1.0) * Complex.Pow(Complex.Abs(r), 4.0))
            >= typeMethods.r8_epsilon() * 72.0 * Math.Pow(Complex.Abs(wp1), 6.0))
        {
            r = z - s * w - Complex.Log(w);
            wp1 = s * w + 1.0;
            e = r / wp1 * (2.0 * wp1 * (wp1 + 2.0 / 3.0 * r) - r)
                / (2.0 * wp1 * (wp1 + 2.0 / 3.0 * r) - 2.0 * r);
            w *= (1.0 + e);
        }

        //
        //  Undo regularization.
        //
        w = s * w;
        //
        //  Provide condition number estimate.
        //
        cond = z / (1.0 + w);

        return 0;
    }
}