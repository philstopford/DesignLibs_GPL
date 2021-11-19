using System;

namespace Burkardt.Praxis;

public static class FLIN
{
    public static double flin(int n, int jsearch, double l, Func<double[], int, double> f,
            double[] x, ref int nf, double[] v, double[] q0, double[] q1, ref double qd0,
            ref double qd1, ref double qa, ref double qb, ref double qc)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FLIN is the function of one variable to be minimized by MINNY.
        //
        //  Discussion:
        //
        //    F(X) is a scalar function of a vector argument X.
        //
        //    A minimizer of F(X) is sought along a line or parabola.
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
        //    If JSEARCH is a legal column index, linear search along V(*,JSEARCH).
        //    If JSEARCH is -1, then the search is parabolic, based on X, Q0 and Q1.
        //
        //    Input, double L, is the parameter determining the particular
        //    point at which F is to be evaluated.  
        //    For a linear search, L is the step size.
        //    For a quadratic search, L is a parameter which specifies
        //    a point in the plane of X, Q0 and Q1.
        //
        //    Input, double F ( double X[], int N ), the function to be minimized.
        //
        //    Input, double X[N], the base point of the search.
        //
        //    Input/output, int &NF, the function evaluation counter.
        //
        //    Input, double V[N,N], a matrix whose columns constitute 
        //    search directions.
        //
        //    Input, double Q0[N], Q1[N], two auxiliary points used to
        //    determine the plane when a quadratic search is performed.
        //
        //    Input, double &QD0, &QD1, values needed to compute the 
        //    coefficients QA, QB, QC.
        //
        //    Output, double &QA, &QB, &QC, coefficients used to combine
        //    Q0, X, and A1 if a quadratic search is used.
        //
        //    Output, double FLIN, the value of the function at the 
        //    minimizing point.
        //
    {
        int i;

        double[] t = new double[n];
        switch (jsearch)
        {
            //
            //  The search is linear.
            //
            case >= 0:
            {
                for (i = 0; i < n; i++)
                {
                    t[i] = x[i] + l * v[i + jsearch * n];
                }

                break;
            }
            //
            default:
            {
                qa = l * (l - qd1) / (qd0 + qd1) / qd0;
                qb = -(l + qd0) * (l - qd1) / qd1 / qd0;
                qc = (l + qd0) * l / qd1 / (qd0 + qd1);

                for (i = 0; i < n; i++)
                {
                    t[i] = qa * q0[i] + qb * x[i] + qc * q1[i];
                }

                break;
            }
        }

        //
        //  The function evaluation counter NF is incremented.
        //
        nf += 1;
        //
        //  Evaluate the function.
        //
        double value = f(t, n);

        return value;
    }
}