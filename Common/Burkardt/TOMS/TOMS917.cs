using System;
using System.Numerics;

namespace Burkardt
{
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
            Complex cond = new Complex();
            Complex e = new Complex();
            Complex r = new Complex();
            Complex w = new Complex();

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
            double pi = Math.PI;
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
            x = (z.Real);
            y = (z.Imaginary);
            // 
            //  Compute if we are near the branch cuts.
            //
            ympi = y - pi;
            yppi = y + pi;
            near = 0.01;
            // 
            //  Test for floating point exceptions:
            //

            //
            //  NaN output for NaN input.
            //
            if ((x == Double.NaN) || (y == Double.NaN))
            {
                w = new Complex((0.0 / 0.0), (0.0 / 0.0));
                e = new Complex(0.0, 0.0);
                r = new Complex(0.0, 0.0);
                return 0;
            }
            //
            //  Signed zeros between branches.
            //
            else if (((x == Double.NegativeInfinity) || (x == Double.PositiveInfinity)) && (x < 0.0) && (-pi < y) &&
                     (y <= pi))
            {
                if (Math.Abs(y) <= pi / 2.0)
                {
                    w = +0.0;
                }
                else
                {
                    w = -0.0;
                }

                if (0.0 <= y)
                {
                    w = w + new Complex(0.0, 0.0);
                }
                else
                {
                    w = w + new Complex(0.0, -1.0 * 0.0);
                }

                e = new Complex(0.0, 0.0);
                r = new Complex(0.0, 0.0);
                return 0;
            }
            //
            //  Asymptotic for large z.
            //
            else if (((x == Double.NegativeInfinity) || (x == Double.PositiveInfinity)) ||
                     (((y == Double.NegativeInfinity) || (y == Double.PositiveInfinity))))
            {
                w = new Complex(x, y);
                e = new Complex(0.0, 0.0);
                r = new Complex(0.0, 0.0);
                return 0;
            }

            //
            //  Test if exactly on the singular points.
            //
            if ((x == -1.0) && (Complex.Abs(y) == pi))
            {
                w = new Complex(-1.0, 0.0);
                e = new Complex(0.0, 0.0);
                r = new Complex(0.0, 0.0);
                return 0;
            }
            // 
            //  Choose approximation based on region.
            //

            //
            //  Region 1: upper branch point.
            //  Series about z=-1+Pi*I.
            //
            if ((-2.0 < x && x <= 1.0 && 1.0 < y && y < 2.0 * pi))
            {
                pz = Complex.Conjugate(Complex.Sqrt(Complex.Conjugate(2.0 * (z + new Complex(1.0, -pi)))));

                w = -1.0
                    + (new Complex(0.0, 1.0)
                       + (1.0 / 3.0
                          + (-1.0 / 36.0 * new Complex(0.0, 1.0)
                             + (1.0 / 270.0 + 1.0 / 4320.0 * new Complex(0.0, 1.0) * pz)
                             * pz) * pz) * pz) * pz;
            }
            //
            //  Region 2: lower branch point.
            //  Series about z=-1-Pi*I.
            //
            else if ((-2.0 < x && x <= 1.0 && -2.0 * pi < y && y < -1.0))
            {
                pz = Complex.Conjugate(Complex.Sqrt(Complex.Conjugate(2.0 * (z + 1.0 + new Complex(0.0, pi)))));

                w = -1.0
                    + (-new Complex(0.0, 1.0) + (1.0 / 3.0
                                                 + (1.0 / 36.0 * new Complex(0.0, 1.0)
                                                    + (1.0 / 270.0 - 1.0 / 4320.0 * new Complex(0.0, 1.0) * pz)
                                                    * pz) * pz) * pz) * pz;
            }
            //
            //  Region 3: between branch cuts.
            //  Series: About -infinity.
            //
            else if (x <= -2.0 && -pi < y && y <= pi)
            {
                pz = Complex.Exp(z);
                w = (1.0
                     + (-1.0
                        + (3.0 / 2.0
                           + (-8.0 / 3.0
                              + 125.0 / 24.0 * pz) * pz) * pz) * pz) * pz;
            }
            //
            //  Region 4: Mushroom.
            //  Series about z=1.
            //
            else if (((-2.0 < x) && (x <= 1.0) && (-1.0 <= y) && (y <= 1.0))
                     || ((-2.0 < x) && (x - 1.0) * (x - 1.0) + y * y <= pi * pi))
            {
                pz = z - 1.0;
                w = 1.0 / 2.0 + 1.0 / 2.0 * z
                              + (1.0 / 16.0
                                 + (-1.0 / 192.0
                                    + (-1.0 / 3072.0 + 13.0 / 61440.0 * pz) * pz) * pz) * pz * pz;
            }
            //
            //  Region 5: Top wing.
            //  Negative log series.
            //
            else if (x <= -1.05 && pi < y && y - pi <= -0.75 * (x + 1.0))
            {
                t = z - new Complex(0.0, pi);
                pz = Complex.Log(-t);
                w = ((1.0 + (-3.0 / 2.0 + 1.0 / 3.0 * pz) * pz) * pz
                     + ((-1.0 + 1.0 / 2.0 * pz) * pz + (pz + (-pz + t) * t) * t) * t)
                    / (t * t * t);
            }
            //
            //  Region 6: Bottom wing.
            //  Negative log series.
            //
            else if (x <= -1.05 && 0.75 * (x + 1.0) < y + pi && y + pi <= 0.0)
            {
                t = z + new Complex(0.0, pi);
                pz = Complex.Log(-t);
                w = ((1.0 + (-3.0 / 2.0 + 1.0 / 3.0 * pz) * pz) * pz
                     + ((-1.0 + 1.0 / 2.0 * pz) * pz + (pz + (-pz + t) * t) * t) * t)
                    / (t * t * t);
            }
            //
            //  Region 7: Everywhere else.
            //  Series solution about infinity.
            //
            else
            {
                pz = Complex.Log(z);
                w = ((1.0 + (-3.0 / 2.0 + 1.0 / 3.0 * pz) * pz) * pz
                     + ((-1.0 + 1.0 / 2.0 * pz) * pz + (pz + (-pz + z) * z) * z) * z)
                    / (z * z * z);
            }

            //
            //  Regularize if near branch cuts.
            ///
            if (x <= -1.0 + near && (Complex.Abs(ympi) <= near || Complex.Abs(yppi) <= near))
            {
                s = -1.0;
                if (Complex.Abs(ympi) <= near)
                {
                    //
                    //  Recompute ympi with directed rounding.
                    //
                    // fesetround ( FE_UPWARD );
                    ympi = y - pi;

                    if (ympi <= 0.0)
                    {
                        // fesetround ( FE_DOWNWARD );
                        ympi = y - pi;
                    }

                    z = new Complex(x, ympi);
                    // 
                    //  Return rounding to default.
                    //
                    // fesetround ( FE_TONEAREST );
                }
                else
                {
                    //
                    //  Recompute yppi with directed rounding.
                    //
                    // fesetround ( FE_UPWARD );
                    yppi = y + pi;

                    if (yppi <= 0.0)
                    {
                        // fesetround ( FE_DOWNWARD );
                        yppi = y + pi;
                    }

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
            w = w * (1.0 + e);
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
                w = w * (1.0 + e);
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
}