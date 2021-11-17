using System;
using Burkardt.Types;

namespace Burkardt.Grid;

public static class Poisson
{
    public static void monogrid_poisson_1d(int n, double a, double b, double ua, double ub,
            Func<double, double> force, Func<double, double> exact, ref int it_num,
            ref double[] u)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MONOGRID_POISSON_1D solves a 1D PDE, using the Gauss-Seidel method.
        //
        //  Discussion:
        //
        //    This routine solves a 1D boundary value problem of the form
        //
        //      - U''(X) = F(X) for A < X < B,
        //
        //    with boundary conditions U(A) = UA, U(B) = UB.
        //
        //    The Gauss-Seidel method is used. 
        //
        //    This routine is provided primarily for comparison with the
        //    multigrid solver.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    William Hager,
        //    Applied Numerical Linear Algebra,
        //    Prentice-Hall, 1988,
        //    ISBN13: 978-0130412942,
        //    LC: QA184.H33.
        //
        //  Parameters:
        //
        //    Input, int N, the number of intervals.
        //
        //    Input, double A, B, the endpoints.
        //
        //    Input, double UA, UB, the boundary values at the endpoints.
        //
        //    Input, double FORCE ( double x ), the name of the function 
        //    which evaluates the right hand side.
        //
        //    Input, double EXACT ( double x ), the name of the function 
        //    which evaluates the exact solution.
        //
        //    Output, int &IT_NUM, the number of iterations.
        //
        //    Output, double U[N+1], the computed solution.
        //
    {
        double d1 = 0;
        double h;
        int i;
        double[] r;
        double tol;
        double[] x;
        //
        //  Initialization.
        //
        tol = 0.0001;
        //
        //  Set the nodes.
        //
        x = typeMethods.r8vec_linspace_new(n + 1, a, b);
        //
        //  Set the right hand side.
        //
        r = new double[n + 1];

        r[0] = ua;
        h = (b - a) / n;
        for (i = 1; i < n; i++)
        {
            r[i] = h * h * force(x[i]);
        }

        r[n] = ub;

        for (i = 0; i <= n; i++)
        {
            u[i] = 0.0;
        }

        it_num = 0;
        //
        //  Gauss-Seidel iteration.
        //
        for (;;)
        {
            it_num += 1;

            GaussSeidel.gauss_seidel(n + 1, r, ref u, ref d1);

            if (d1 <= tol)
            {
                break;
            }
        }
    }

    public static void multigrid_poisson_1d(int n, double a, double b, double ua, double ub,
            Func<double, double> force, Func<double, double> exact, ref int it_num,
            ref double[] u)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MULTIGRID_POISSON_1D solves a 1D PDE using the multigrid method.
        //
        //  Discussion:
        //
        //    This routine solves a 1D boundary value problem of the form
        //
        //      - U''(X) = F(X) for A < X < B,
        //
        //    with boundary conditions U(A) = UA, U(B) = UB.
        //
        //    The multigrid method is used. 
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 July 2014
        //
        //  Author:
        //
        //    Original FORTRAN77 version by William Hager.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    William Hager,
        //    Applied Numerical Linear Algebra,
        //    Prentice-Hall, 1988,
        //    ISBN13: 978-0130412942,
        //    LC: QA184.H33.
        //
        //  Parameters:
        //
        //    Input, int N, the number of intervals.
        //    N must be a power of 2.
        //
        //    Input, double A, B, the endpoints.
        //
        //    Input, double UA, UB, the boundary values at the endpoints.
        //
        //    Input, double FORCE ( double x ), the name of the function 
        //    which evaluates the right hand side.
        //
        //    Input, double EXACT ( double x ), the name of the function 
        //    which evaluates the exact solution.
        //
        //    Output, int &IT_NUM, the number of iterations.
        //
        //    Output, double U[N+1], the computed solution.
        //
    {
        double d0;
        double d1;
        double h;
        int i;
        int it;
        int j;
        int k;
        int l;
        int ll;
        int m;
        int nl;
        double[] r;
        double tol;
        double utol;
        double[] uu;
        double[] x;
        //
        //  Determine if we have enough storage.
        //
        k = (int) Math.Log2(n);

        if (n != Math.Pow(2, k))
        {
            Console.WriteLine("");
            Console.WriteLine("MULTIGRID_POISSON_1D - Fatal error!");
            Console.WriteLine("  N is not a power of 2.");
            return;
        }

        nl = n + n + k - 2;
        //
        //  Initialization.
        //
        it = 4;
        it_num = 0;
        tol = 0.0001;
        utol = 0.7;
        m = n;
        //
        //  Set the nodes.
        //
        x = typeMethods.r8vec_linspace_new(n + 1, a, b);
        //
        //  Set the right hand side.
        //
        r = new double[nl];
        r[0] = ua;
        h = (b - a) / n;
        for (i = 1; i < n; i++)
        {
            r[i] = h * h * force(x[i]);
        }

        r[n] = ub;

        uu = new double[nl];

        for (i = 0; i < nl; i++)
        {
            uu[i] = 0.0;
        }

        //
        //  L points to first entry of solution
        //  LL points to penultimate entry.
        //
        l = 0;
        ll = n - 1;
        //
        //  Gauss-Seidel iteration
        //
        d1 = 0.0;
        j = 0;

        for (;;)
        {
            d0 = d1;
            j += 1;
            GaussSeidel.gauss_seidel(n + 1, r, ref uu, ref d1, rIndex: +l, uIndex: +l);
            it_num += 1;
            //
            //  Do at least 4 iterations at each level.
            //
            if (j < it)
            {
            }
            //
            //  Enough iterations, satisfactory decrease, on finest grid, exit.
            //
            else if (d1 < tol && n == m)
            {
                break;
            }
            //
            //  Enough iterations, satisfactory convergence, go finer.
            //
            else if (d1 < tol)
            {
                Transfer.ctof(n + 1, uu, n + n + 1, ref uu, ucIndex: +l, ufIndex: +(l - 1 - n - n));

                n += n;
                ll = l - 2;
                l = l - 1 - n;
                j = 0;
            }
            //
            //  Enough iterations, slow convergence, 2 < N, go coarser.
            //
            else if (utol * d0 <= d1 && 2 < n)
            {
                Transfer.ftoc(n + 1, uu, r, n / 2 + 1, ref uu, ref r, ufIndex: +l, rfIndex: +l,
                    ucIndex: +(l + n + 1), rcIndex: +(l + n + 1));

                n /= 2;
                l = ll + 2;
                ll = ll + n + 1;
                j = 0;
            }
        }

        for (i = 0; i < n + 1; i++)
        {
            u[i] = uu[i];
        }
    }
}