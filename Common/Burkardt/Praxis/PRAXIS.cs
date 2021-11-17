using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Praxis;

public static class PRAXIS
{
    public static double praxis(double t0, double h0, int n, int prin, ref double[] x,
            Func<double[], int, double> f)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PRAXIS seeks an N-dimensional minimizer X of a scalar function F(X).
        //
        //  Discussion:
        //
        //    PRAXIS returns the minimum of the function F(X,N) of N variables
        //    using the principal axis method.  The gradient of the function is
        //    not required.
        //
        //    The approximating quadratic form is
        //
        //      Q(x") = F(x,n) + (1/2) * (x"-x)" * A * (x"-x)
        //
        //    where X is the best estimate of the minimum and 
        //
        //      A = inverse(V") * D * inverse(V)
        //
        //    V(*,*) is the matrix of search directions; 
        //    D(*) is the array of second differences.  
        //
        //    If F(X) has continuous second derivatives near X0, then A will tend 
        //    to the hessian of F at X0 as X approaches X0.
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
        //    Input, double T0, is a tolerance.  PRAXIS attempts to return 
        //    praxis = f(x) such that if X0 is the true local minimum near X, then
        //    norm ( x - x0 ) < T0 + Math.Sqrt ( EPSILON ( X ) ) * norm ( X ),
        //    where EPSILON ( X ) is the machine precision for X.
        //
        //    Input, double H0, is the maximum step size.  H0 should be 
        //    set to about the maximum distance from the initial guess to the minimum.
        //    If H0 is set too large or too small, the initial rate of
        //    convergence may be slow.
        //
        //    Input, int N, the number of variables.
        //
        //    Input, int PRIN, controls printing intermediate results.
        //    0, nothing is printed.
        //    1, F is printed after every n+1 or n+2 linear minimizations.  
        //       final X is printed, but intermediate X is printed only 
        //       if N is at most 4.
        //    2, the scale factors and the principal values of the approximating 
        //       quadratic form are also printed.
        //    3, X is also printed after every few linear minimizations.
        //    4, the principal vectors of the approximating quadratic form are 
        //       also printed.
        //
        //    Input/output, double X[N], is an array containing on entry a
        //    guess of the point of minimum, on return the estimated point of minimum.
        //
        //    Input, double F ( double X[], int N ), is the name of the function to be
        //    minimized.
        //
        //    Output, double PRAXIS, the function value at the minimizer.
        //
        //  Local parameters:
        //
        //    Local, double DMIN, an estimate for the smallest eigenvalue.
        //
        //    Local, double FX, the value of F(X,N).
        //
        //    Local, bool ILLC, is TRUE if the system is ill-conditioned.
        //
        //    Local, double LDT, the length of the step.
        //
        //    Local, int NF, the number of function evaluations.
        //
        //    Local, int NL, the number of linear searches.
        //
    {
        double[] d;
        double d2;
        double df;
        double dmin;
        double dn;
        double dni;
        double f1;
        bool fk;
        double fx;
        double h;
        int i;
        bool illc;
        int j;
        int jsearch;
        int k;
        int k2;
        int kl;
        int kt;
        int ktm;
        double large;
        double ldfac;
        double lds;
        double ldt;
        double m2;
        double m4;
        double machep;
        int nits;
        int nl;
        int nf;
        double[] q0;
        double[] q1;
        double qa;
        double qb;
        double qc;
        double qd0;
        double qd1;
        double qf1;
        double r;
        double s;
        double scbd;
        int seed;
        double sf;
        double sl;
        double small;
        double t;
        double temp;
        double t2;
        double[] v;
        double value = 0;
        double vlarge;
        double vsmall;
        double[] y;
        double[] z;
        //
        //  Allocation.
        //
        d = new double[n];
        q0 = new double[n];
        q1 = new double[n];
        v = new double[n * n];
        y = new double[n];
        z = new double[n];
        //
        //  Initialization.
        //
        machep = typeMethods.r8_epsilon();
        small = machep * machep;
        vsmall = small * small;
        large = 1.0 / small;
        vlarge = 1.0 / vsmall;
        m2 = Math.Sqrt(machep);
        m4 = Math.Sqrt(m2);
        seed = 123456789;
        //
        //  Heuristic numbers:
        //
        //  If the axes may be badly scaled (which is to be avoided if
        //  possible), then set SCBD = 10.  Otherwise set SCBD = 1.
        //
        //  If the problem is known to be ill-conditioned, initialize ILLC = true.
        //
        //  KTM is the number of iterations without improvement before the
        //  algorithm terminates.  KTM = 4 is very cautious; usually KTM = 1
        //  is satisfactory.
        //
        scbd = 1.0;
        illc = false;
        ktm = 1;

        ldfac = illc switch
        {
            true => 0.1,
            _ => 0.01
        };

        kt = 0;
        nl = 0;
        nf = 1;
        fx = f(x, n);
        qf1 = fx;
        t = small + Math.Abs(t0);
        t2 = t;
        dmin = small;
        h = h0;
        h = Math.Max(h, 100.0 * t);
        ldt = h;
        //
        //  The initial set of search directions V is the identity matrix.
        //
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                v[i + j * n] = 0.0;
            }

            v[j + j * n] = 1.0;
        }

        for (i = 0; i < n; i++)
        {
            d[i] = 0.0;
        }

        qa = 0.0;
        qb = 0.0;
        qc = 0.0;
        qd0 = 0.0;
        qd1 = 0.0;
        typeMethods.r8vec_copy(n, x, ref q0);
        typeMethods.r8vec_copy(n, x, ref q1);

        switch (prin)
        {
            case > 0:
                print2(n, x, prin, fx, nf, nl);
                break;
        }

        //
        //  The main loop starts here.
        //
        for (;;)
        {
            sf = d[0];
            d[0] = 0.0;
            //
            //  Minimize along the first direction V(*,1).
            //
            jsearch = 0;
            nits = 2;
            d2 = d[0];
            s = 0.0;
            value = fx;
            fk = false;

            MINNY.minny(n, jsearch, nits, ref d2, ref s, ref value, fk, f, x, t,
                h, v, q0, q1, ref nl, ref nf, dmin, ldt, ref fx, ref qa, ref qb, ref qc, ref qd0, ref qd1);

            d[0] = d2;

            switch (s)
            {
                case <= 0.0:
                {
                    for (i = 0; i < n; i++)
                    {
                        v[i + 0 * n] = -v[i + 0 * n];
                    }

                    break;
                }
            }

            if (sf <= 0.9 * d[0] || d[0] <= 0.9 * sf)
            {
                for (i = 1; i < n; i++)
                {
                    d[i] = 0.0;
                }
            }

            //
            //  The inner loop starts here.
            //
            for (k = 2; k <= n; k++)
            {
                typeMethods.r8vec_copy(n, x, ref y);

                sf = fx;

                illc = kt switch
                {
                    > 0 => true,
                    _ => illc
                };

                for (;;)
                {
                    kl = k;
                    df = 0.0;
                    switch (illc)
                    {
                        //
                        //  A random step follows, to avoid resolution valleys.
                        //
                        case true:
                        {
                            for (j = 0; j < n; j++)
                            {
                                r = UniformRNG.r8_uniform_01(ref seed);
                                s = (0.1 * ldt + t2 * Math.Pow(10.0, kt)) * (r - 0.5);
                                z[j] = s;
                                for (i = 0; i < n; i++)
                                {
                                    x[i] += s * v[i + j * n];
                                }
                            }

                            fx = f(x, n);
                            nf += 1;
                            break;
                        }
                    }

                    //
                    //  Minimize along the "non-conjugate" directions V(*,K),...,V(*,N).
                    //
                    for (k2 = k; k2 <= n; k2++)
                    {
                        sl = fx;

                        jsearch = k2 - 1;
                        nits = 2;
                        d2 = d[k2 - 1];
                        s = 0.0;
                        value = fx;
                        fk = false;

                        MINNY.minny(n, jsearch, nits, ref d2, ref s, ref value, fk, f, x, t,
                            h, v, q0, q1, ref nl, ref nf, dmin, ldt, ref fx, ref qa, ref qb, ref qc, ref qd0,
                            ref qd1);

                        d[k2 - 1] = d2;

                        s = illc switch
                        {
                            true => d[k2 - 1] * Math.Pow(s + z[k2 - 1], 2),
                            _ => sl - fx
                        };

                        if (df <= s)
                        {
                            df = s;
                            kl = k2;
                        }
                    }

                    //
                    //  If there was not much improvement on the first try, set
                    //  ILLC = true and start the inner loop again.
                    //
                    if (illc)
                    {
                        break;
                    }

                    if (Math.Abs(100.0 * machep * fx) <= df)
                    {
                        break;
                    }

                    illc = true;
                }

                switch (k)
                {
                    case 2 when 1 < prin:
                        typeMethods.r8vec_print(n, d, "  The second difference array:");
                        break;
                }

                //
                //  Minimize along the "conjugate" directions V(*,1),...,V(*,K-1).
                //
                for (k2 = 1; k2 < k; k2++)
                {
                    jsearch = k2 - 1;
                    nits = 2;
                    d2 = d[k2 - 1];
                    s = 0.0;
                    value = fx;
                    fk = false;

                    MINNY.minny(n, jsearch, nits, ref d2, ref s, ref value, fk, f, x, t,
                        h, v, q0, q1, ref nl, ref nf, dmin, ldt, ref fx, ref qa, ref qb, ref qc, ref qd0, ref qd1);

                    d[k2 - 1] = d2;
                }

                f1 = fx;
                fx = sf;

                for (i = 0; i < n; i++)
                {
                    temp = x[i];
                    x[i] = y[i];
                    y[i] = temp - y[i];
                }

                lds = typeMethods.r8vec_norm(n, y);
                //
                //  Discard direction V(*,kl).
                //
                //  If no random step was taken, V(*,KL) is the "non-conjugate"
                //  direction along which the greatest improvement was made.
                //
                if (small < lds)
                {
                    for (j = kl - 1; k <= j; j--)
                    {
                        for (i = 1; i <= n; i++)
                        {
                            v[i - 1 + j * n] = v[i - 1 + (j - 1) * n];
                        }

                        d[j] = d[j - 1];
                    }

                    d[k - 1] = 0.0;

                    for (i = 1; i <= n; i++)
                    {
                        v[i - 1 + (k - 1) * n] = y[i - 1] / lds;
                    }

                    //
                    //  Minimize along the new "conjugate" direction V(*,k), which is
                    //  the normalized vector:  (new x) - (old x).
                    //
                    jsearch = k - 1;
                    nits = 4;
                    d2 = d[k - 1];
                    value = f1;
                    fk = true;

                    MINNY.minny(n, jsearch, nits, ref d2, ref lds, ref value, fk, f, x, t,
                        h, v, q0, q1, ref nl, ref nf, dmin, ldt, ref fx, ref qa, ref qb, ref qc, ref qd0, ref qd1);

                    d[k - 1] = d2;

                    switch (lds)
                    {
                        case <= 0.0:
                        {
                            lds = -lds;
                            for (i = 1; i <= n; i++)
                            {
                                v[i - 1 + (k - 1) * n] = -v[i - 1 + (k - 1) * n];
                            }

                            break;
                        }
                    }
                }

                ldt = ldfac * ldt;
                ldt = Math.Max(ldt, lds);

                switch (prin)
                {
                    case > 0:
                        print2(n, x, prin, fx, nf, nl);
                        break;
                }

                t2 = typeMethods.r8vec_norm(n, x);

                t2 = m2 * t2 + t;
                //
                //  See whether the length of the step taken since starting the
                //  inner loop exceeds half the tolerance.
                //
                if (0.5 * t2 < ldt)
                {
                    kt = -1;
                }

                kt += 1;

                if (ktm < kt)
                {
                    switch (prin)
                    {
                        case > 0:
                            typeMethods.r8vec_print(n, x, "  X:");
                            break;
                    }

                    return fx;
                }
            }

            //
            //  The inner loop ends here.
            //
            //  Try quadratic extrapolation in case we are in a curved valley.
            //
            QUAD.quad(n, f, ref x, t, h, v, ref q0, ref q1, ref nl, ref nf, dmin, ldt, ref fx, ref qf1,
                ref qa, ref qb, ref qc, ref qd0, ref qd1);

            for (j = 0; j < n; j++)
            {
                d[j] = 1.0 / Math.Sqrt(d[j]);
            }

            dn = typeMethods.r8vec_max(n, d);

            switch (prin)
            {
                case > 3:
                    typeMethods.r8mat_print(n, n, v, "  The new direction vectors:");
                    break;
            }

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    v[i + j * n] = d[j] / dn * v[i + j * n];
                }
            }

            switch (scbd)
            {
                //
                //  Scale the axes to try to reduce the condition number.
                //
                case > 1.0:
                {
                    for (i = 0; i < n; i++)
                    {
                        s = 0.0;
                        for (j = 0; j < n; j++)
                        {
                            s += v[i + j * n] * v[i + j * n];
                        }

                        s = Math.Sqrt(s);
                        z[i] = Math.Max(m4, s);
                    }

                    s = typeMethods.r8vec_min(n, z);

                    for (i = 0; i < n; i++)
                    {
                        sl = s / z[i];
                        z[i] = 1.0 / sl;

                        if (scbd < z[i])
                        {
                            sl = 1.0 / scbd;
                            z[i] = scbd;
                        }

                        for (j = 0; j < n; j++)
                        {
                            v[i + j * n] = sl * v[i + j * n];
                        }
                    }

                    break;
                }
            }

            //
            //  Calculate a new set of orthogonal directions before repeating
            //  the main loop.
            //
            //  Transpose V for MINFIT:
            //
            typeMethods.r8mat_transpose_in_place(n, ref v);
            //
            //  MINFIT finds the singular value decomposition of V.
            //
            //  This gives the principal values and principal directions of the
            //  approximating quadratic form without squaring the condition number.
            //
            MINFIT.minfit(n, vsmall, ref v, ref d);
            switch (scbd)
            {
                //
                //  Unscale the axes.
                //
                case > 1.0:
                {
                    for (i = 0; i < n; i++)
                    {
                        for (j = 0; j < n; j++)
                        {
                            v[i + j * n] = z[i] * v[i + j * n];
                        }
                    }

                    for (j = 0; j < n; j++)
                    {
                        s = 0.0;
                        for (i = 0; i < n; i++)
                        {
                            s += v[i + j * n] * v[i + j * n];
                        }

                        s = Math.Sqrt(s);

                        d[j] = s * d[j];
                        for (i = 0; i < n; i++)
                        {
                            v[i + j * n] /= s;
                        }
                    }

                    break;
                }
            }

            for (i = 0; i < n; i++)
            {
                dni = dn * d[i];

                if (large < dni)
                {
                    d[i] = vsmall;
                }
                else if (dni < small)
                {
                    d[i] = vlarge;
                }
                else
                {
                    d[i] = 1.0 / dni / dni;
                }
            }

            //
            //  Sort the eigenvalues and eigenvectors.
            //
            SVSORT.svsort(n, ref d, ref v);
            //
            //  Determine the smallest eigenvalue.
            //
            dmin = Math.Max(d[n - 1], small);
            //
            //  The ratio of the smallest to largest eigenvalue determines whether
            //  the system is ill conditioned.
            //
            if (dmin < m2 * d[0])
            {
                illc = true;
            }
            else
            {
                illc = false;
            }

            switch (prin)
            {
                case > 1:
                {
                    switch (scbd)
                    {
                        case > 1.0:
                            typeMethods.r8vec_print(n, z, "  The scale factors:");
                            break;
                    }

                    typeMethods.r8vec_print(n, d, "  Principal values of the quadratic form:");
                    break;
                }
            }

            switch (prin)
            {
                case > 3:
                    typeMethods.r8mat_print(n, n, v, "  The principal axes:");
                    break;
            }
            //
            //  The main loop ends here.
            //
        }

        switch (prin)
        {
            case > 0:
                typeMethods.r8vec_print(n, x, "  X:");
                break;
        }

        return fx;
    }
        
    public static void print2 ( int n, double[] x, int prin, double fx, int nf, int nl )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PRINT2 prints certain data about the progress of the iteration.
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
        //    Input, double X[N], the current estimate of the minimizer.
        //
        //    Input, int PRIN, the user-specifed print level.
        //    0, nothing is printed.
        //    1, F is printed after every n+1 or n+2 linear minimizations.  
        //       final X is printed, but intermediate X is printed only 
        //       if N is at most 4.
        //    2, the scale factors and the principal values of the approximating 
        //       quadratic form are also printed.
        //    3, X is also printed after every few linear minimizations.
        //    4, the principal vectors of the approximating quadratic form are 
        //       also printed.
        //
        //    Input, double FX, the smallest value of F(X) found so far.
        //
        //    Input, int NF, the number of function evaluations.
        //
        //    Input, int NL, the number of linear searches.
        //
    {
        Console.WriteLine("");
        Console.WriteLine("  Linear searches = " + nl + "");
        Console.WriteLine("  Function evaluations " + nf + "");
        Console.WriteLine("  Function value FX = " + fx + "");

        if ( n <= 4 || 2 < prin )
        {
            typeMethods.r8vec_print ( n, x, "  X:" );
        }
    }
}