using System;
using Burkardt.Types;

namespace Burkardt.Praxis
{
    public static class MINNY
    {
        public static void minny(int n, int jsearch, int nits, ref double d2, ref double x1, ref double f1,
                bool fk, Func<double[], int, double> f, double[] x, double t, double h,
                double[] v, double[] q0, double[] q1, ref int nl, ref int nf, double dmin,
                double ldt, ref double fx, ref double qa, ref double qb, ref double qc, ref double qd0,
                ref double qd1)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MINNY minimizes a scalar function of N variables along a line.
            //
            //  Discussion:
            //
            //    MINNY minimizes F along the line from X in the direction V(*,JSEARCH) 
            //    or else using a quadratic search in the plane defined by Q0, Q1 and X.
            //
            //    If FK = true, then F1 is FLIN(X1).  Otherwise X1 and F1 are ignored
            //    on entry unless final FX is greater than F1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 August 2016
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Richard Brent.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Richard Brent,
            //    Algorithms for Minimization with Derivatives,
            //    Prentice Hall, 1973,
            //    Reprinted by Dover, 2002.
            //
            //  Parameters:
            //
            //    Input, int N, the number of variables.
            //
            //    Input, int JSEARCH, indicates the kind of search.
            //    If J is a legal columnindex, linear search in the direction of V(*,JSEARCH).
            //    Otherwise, the search is parabolic, based on X, Q0 and Q1.
            //
            //    Input, int NITS, the maximum number of times the interval 
            //    may be halved to retry the calculation.
            //
            //    Input/output, ref double D2, is either zero, or an approximation to 
            //    the value of (1/2) times the second derivative of F.
            //
            //    Input/output, ref double X1, on entry, an estimate of the 
            //    distance from X to the minimum along V(*,JSEARCH), or a curve.  
            //    On output, the distance between X and the minimizer that was found.
            //
            //    Input/output, ref double F1, ?
            //
            //    Input, bool FK; if FK is TRUE, then on input F1 contains 
            //    the value FLIN(X1).
            //
            //    Input, double F ( double X[], int N ), is the name of the function to 
            //    be minimized.
            //
            //    Input/output, double X[N], ?
            //
            //    Input, double T, ?
            //
            //    Input, double H, ?
            //
            //    Input, double V[N,N], a matrix whose columns are direction
            //    vectors along which the function may be minimized.
            //
            //    ?, double Q0[N], ?
            //
            //    ?, double Q1[N], ?
            //
            //    Input/output, int &NL, the number of linear searches.
            //
            //    Input/output, int &NF, the number of function evaluations.
            //
            //    Input, double DMIN, an estimate for the smallest eigenvalue.
            //
            //    Input, double LDT, the length of the step.
            //
            //    Input/output, ref double FX, the value of F(X,N).
            //
            //    Input/output, ref double QA, &QB, &QC;
            //
            //    Input/output, ref double QD0, &QD1, ?.
            //
        {
            double d1;
            bool dz;
            double f0;
            double f2;
            double fm;
            int i;
            int k;
            double m2;
            double m4;
            double machep;
            bool ok;
            double s;
            double sf1;
            double small;
            double sx1;
            double t2;
            double temp;
            double x2;
            double xm;

            machep = typeMethods.r8_epsilon();
            small = machep * machep;
            m2 = Math.Sqrt(machep);
            m4 = Math.Sqrt(m2);
            sf1 = f1;
            sx1 = x1;
            k = 0;
            xm = 0.0;
            fm = fx;
            f0 = fx;
            dz = (d2 < machep);
            //
            //  Find the step size.
            //
            s = typeMethods.r8vec_norm(n, x);

            if (dz)
            {
                temp = dmin;
            }
            else
            {
                temp = d2;
            }

            t2 = m4 * Math.Sqrt(Math.Abs(fx) / temp + s * ldt) + m2 * ldt;
            s = m4 * s + t;
            if (dz && s < t2)
            {
                t2 = s;
            }

            t2 = Math.Max(t2, small);
            t2 = Math.Min(t2, 0.01 * h);

            if (fk && f1 <= fm)
            {
                xm = x1;
                fm = f1;
            }

            if ((!fk) || Math.Abs(x1) < t2)
            {
                if (0.0 <= x1)
                {
                    temp = 1.0;
                }
                else
                {
                    temp = -1.0;
                }

                x1 = temp * t2;
                f1 = FLIN.flin(n, jsearch, x1, f, x, ref nf, v, q0, q1, ref qd0, ref qd1, ref qa, ref qb, ref qc);
            }

            if (f1 <= fm)
            {
                xm = x1;
                fm = f1;
            }

            //
            //  Evaluate FLIN at another point and estimate the second derivative.
            //
            for (;;)
            {
                if (dz)
                {
                    if (f1 <= f0)
                    {
                        x2 = 2.0 * x1;
                    }
                    else
                    {
                        x2 = -x1;
                    }

                    f2 = FLIN.flin(n, jsearch, x2, f, x, ref nf, v, q0, q1, ref qd0, ref qd1, ref qa, ref qb, ref qc);

                    if (f2 <= fm)
                    {
                        xm = x2;
                        fm = f2;
                    }

                    d2 = (x2 * (f1 - f0) - x1 * (f2 - f0))
                         / ((x1 * x2) * (x1 - x2));
                }

                //
                //  Estimate the first derivative at 0.
                //
                d1 = (f1 - f0) / x1 - x1 * d2;
                dz = true;
                //
                //  Predict the minimum.
                //
                if (d2 <= small)
                {
                    if (0.0 <= d1)
                    {
                        x2 = -h;
                    }
                    else
                    {
                        x2 = h;
                    }
                }
                else
                {
                    x2 = (-0.5 * d1) / d2;
                }

                if (h < Math.Abs(x2))
                {
                    if (x2 <= 0.0)
                    {
                        x2 = -h;
                    }
                    else
                    {
                        x2 = h;
                    }
                }

                //
                //  Evaluate F at the predicted minimum.
                //
                ok = true;

                for (;;)
                {
                    f2 = FLIN.flin(n, jsearch, x2, f, x, ref nf, v, q0, q1, ref qd0, ref qd1, ref qa, ref qb, ref qc);

                    if (nits <= k || f2 <= f0)
                    {
                        break;
                    }

                    k = k + 1;

                    if (f0 < f1 && 0.0 < x1 * x2)
                    {
                        ok = false;
                        break;
                    }

                    x2 = 0.5 * x2;
                }

                if (ok)
                {
                    break;
                }
            }

            //
            //  Increment the one-dimensional search counter.
            //
            nl = nl + 1;

            if (fm < f2)
            {
                x2 = xm;
            }
            else
            {
                fm = f2;
            }

            //
            //  Get a new estimate of the second derivative.
            //
            if (small < Math.Abs(x2 * (x2 - x1)))
            {
                d2 = (x2 * (f1 - f0) - x1 * (fm - f0))
                     / ((x1 * x2) * (x1 - x2));
            }
            else
            {
                if (0 < k)
                {
                    d2 = 0.0;
                }
            }

            d2 = Math.Max(d2, small);

            x1 = x2;
            fx = fm;

            if (sf1 < fx)
            {
                fx = sf1;
                x1 = sx1;
            }

            //
            //  Update X for linear search.
            //
            if (0 <= jsearch)
            {
                for (i = 0; i < n; i++)
                {
                    x[i] = x[i] + x1 * v[i + jsearch * n];
                }
            }

        }
    }
}