using System;

namespace Burkardt.Praxis;

public static class QUAD
{
    public static void quad(int n, Func < double[], int, double > f, ref double[] x, double t,
            double h, double[] v, ref double[] q0, ref double[] q1, ref int  nl, ref int  nf, double dmin,
            double ldt, ref double fx, ref double qf1, ref double qa, ref double qb, ref double qc,
            ref double qd0, ref double qd1 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    QUAD seeks to minimize the scalar function F along a particular curve.
        //
        //  Discussion:
        //
        //    The minimizer to be sought is required to lie on a curve defined
        //    by Q0, Q1 and X.
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
        //    Input, double F ( double X[], int N ), the name of the function to 
        //    be minimized.
        //
        //    Input/output, double X[N], ?
        //
        //    Input, double T, ?
        //
        //    Input, double H, ?
        //
        //    Input, double V[N,N], the matrix of search directions.
        //
        //    Input/output, double Q0[N], Q1[N], two auxiliary points used to define
        //    a curve through X.
        //
        //    Input/output, ref int  NL, the number of linear searches.
        //
        //    Input/output, ref int  NF, the number of function evaluations.
        //
        //    Input, double DMIN, an estimate for the smallest eigenvalue.
        //
        //    Input, double LDT, the length of the step.
        //
        //    Input/output, ref double FX, the value of F(X,N).
        //
        //    Input/output, ref double QF1, &QA, &QB, &QC, &QD0, &QD1 ?
        //
    {
        bool fk;
        int i;
        int jsearch;
        double l;
        int nits;
        double s;
        double temp;
        double value = 0;

        temp = fx;
        fx = qf1;
        qf1 = temp;

        for (i = 0; i < n; i++)
        {
            temp = x[i];
            x[i] = q1[i];
            q1[i] = temp;
        }

        qd1 = 0.0;
        for (i = 0; i < n; i++)
        {
            qd1 += (x[i] - q1[i]) * (x[i] - q1[i]);
        }

        qd1 = Math.Sqrt(qd1);

        if (qd0 <= 0.0 || qd1 <= 0.0 || nl < 3 * n * n)
        {
            fx = qf1;
            qa = 0.0;
            qb = 0.0;
            qc = 1.0;
            s = 0.0;
        }
        else
        {
            jsearch = -1;
            nits = 2;
            s = 0.0;
            l = qd1;
            value = qf1;
            fk = true;

            MINNY.minny(n, jsearch, nits, ref s, ref l, ref value, fk, f, x, t,
                h, v, q0, q1, ref nl, ref nf, dmin, ldt, ref fx, ref qa, ref qb, ref qc, ref qd0, ref qd1);

            qa = l * (l - qd1) / (qd0 + qd1) / qd0;
            qb = -(l + qd0) * (l - qd1) / qd1 / qd0;
            qc = (l + qd0) * l / qd1 / (qd0 + qd1);
        }

        qd0 = qd1;

        for (i = 0; i < n; i++)
        {
            s = q0[i];
            q0[i] = x[i];
            x[i] = qa * s + qb * x[i] + qc * q1[i];
        }

    }
}