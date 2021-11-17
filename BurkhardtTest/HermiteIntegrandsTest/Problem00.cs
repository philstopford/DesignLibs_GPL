using System;
using Burkardt;
using Burkardt.PDFLib;
using Burkardt.Probability;
using Burkardt.Quadrature;
using Burkardt.Types;

namespace HermiteIntegrandsTest;

public static class Problem00
{
    public class p00Data
    {
        public Problem06.p06Data p6data = new();
    }
    public static double p00_exact(ref p00Data data, int problem)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P00_EXACT returns the exact integral for any problem.
        //
        //  Discussion:
        //
        //    This routine provides a "generic" interface to the exact integral
        //    routines for the various problems, and allows a problem to be called
        //    by index (PROBLEM) rather than by name.
        //
        //    In most cases, the "exact" value of the integral is not given;
        //    instead a "respectable" approximation is available.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 Julyl 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int PROBLEM, the index of the problem.
        //
        //    Output, double P00_EXACT, the exact value of the integral.
        //
    {
        double exact;

        switch (problem)
        {
            case 1:
                exact = Problem01.p01_exact();
                break;
            case 2:
                exact = Problem02.p02_exact();
                break;
            case 3:
                exact = Problem03.p03_exact();
                break;
            case 4:
                exact = Problem04.p04_exact();
                break;
            case 5:
                exact = Problem05.p05_exact();
                break;
            case 6:
                exact = Problem06.p06_exact(ref data.p6data);
                break;
            case 7:
                exact = Problem07.p07_exact();
                break;
            case 8:
                exact = Problem08.p08_exact();
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("P00_EXACT - Fatal error!");
                Console.WriteLine("  Illegal problem number = " + problem + "");
                return 1;
        }

        return exact;
    }

    public static void p00_fun(ref p00Data data, int problem, int option, int n, double[] x, ref double[] f )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P00_FUN evaluates the integrand for any problem.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int PROBLEM, the index of the problem.
        //
        //    Input, int OPTION:
        //    0, integrand is f(x).
        //    1, integrand is exp(-x*x) * f(x);
        //    2, integrand is exp(-x*x/2) * f(x);
        //
        //    Input, int N, the number of points.
        //
        //    Input, double X[N], the evaluation points.
        //
        //    Output, double F[N], the function values.
        //
    {
        switch (problem)
        {
            case 1:
                Problem01.p01_fun(option, n, x, ref f);
                break;
            case 2:
                Problem02.p02_fun(option, n, x, ref f);
                break;
            case 3:
                Problem03.p03_fun(option, n, x, ref f);
                break;
            case 4:
                Problem04.p04_fun(option, n, x, ref f);
                break;
            case 5:
                Problem05.p05_fun(option, n, x, ref f);
                break;
            case 6:
                Problem06.p06_fun(ref data.p6data, option, n, x, ref f);
                break;
            case 7:
                Problem07.p07_fun(option, n, x, ref f);
                break;
            case 8:
                Problem08.p08_fun(option, n, x, ref f);
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("P00_FUN - Fatal error!");
                Console.WriteLine("  Illegal problem number = " + problem + "");
                break;
        }
    }

    public static double p00_gauss_hermite(ref p00Data data, int problem, int order)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P00_GAUSS_HERMITE applies a Gauss-Hermite quadrature rule.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int PROBLEM, the index of the problem.
        //
        //    Input, int ORDER, the order of the rule to apply.
        //
        //    Output, double P00_GAUSS_HERMITE, the approximate integral.
        //
    {
        double[] f_vec;
        int option;
        double result;
        double[] weight;
        double[] xtab;

        f_vec = new double[order];
        weight = new double[order];
        xtab = new double[order];

        GaussHermite.hermite_compute(order, ref xtab, ref weight);

        option = 1;
        p00_fun(ref data, problem, option, order, xtab, ref f_vec);

        result = typeMethods.r8vec_dot_product(order, weight, f_vec);

        return result;
    }

    public static double p00_monte_carlo(ref p00Data data, int problem, int order)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P00_MONTE_CARLO applies a Monte Carlo procedure to Hermite integrals.
        //
        //  Discussion:
        //
        //    We wish to estimate the integral:
        //
        //      I(f) = integral ( -oo < x < +oo ) f(x) exp ( - x * x ) dx
        //
        //    We do this by a Monte Carlo sampling procedure, in which 
        //    we select N points X(1:N) from a standard normal distribution,
        //    and estimate
        //
        //      Q(f) = sum ( 1 <= I <= N ) f(x(i)) / sqrt ( pi )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int PROBLEM, the index of the problem.
        //
        //    Input, int ORDER, the order of the Gauss-Laguerre rule
        //    to apply.
        //
        //    Output, double P00_MONTE_CARLO, the approximate integral.
        //
    {
        double[] f_vec;
        int option;
        const double r8_pi = 3.141592653589793;
        double result;
        int seed;
        double weight;
        double[] x_vec;

        seed = 123456789;
        typeMethods.r8vecNormalData ndata = new();
        x_vec = typeMethods.r8vec_normal_01_new(order, ref ndata, ref seed);

        option = 2;
        f_vec = new double[order];

        p00_fun(ref data, problem, option, order, x_vec, ref f_vec);

        weight = order / Math.Sqrt(r8_pi) / Math.Sqrt(2.0);

        result = typeMethods.r8vec_sum(order, f_vec) / weight;

        return result;
    }

    public static int p00_problem_num()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P00_PROBLEM_NUM returns the number of test integration problems.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int P00_PROBLEM_NUM, the number of test problems.
        //
    {
        int problem_num;

        problem_num = 8;

        return problem_num;
    }

    public static string p00_title(int problem)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P00_TITLE returns the title for any problem.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 February 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int PROBLEM, the index of the problem.
        //
        //    Output, string P00_TITLE, the title of the problem.
        //
    {
        string title;

        switch (problem)
        {
            case 1:
                title = Problem01.p01_title();
                break;
            case 2:
                title = Problem02.p02_title();
                break;
            case 3:
                title = Problem03.p03_title();
                break;
            case 4:
                title = Problem04.p04_title();
                break;
            case 5:
                title = Problem05.p05_title();
                break;
            case 6:
                title = Problem06.p06_title();
                break;
            case 7:
                title = Problem07.p07_title();
                break;
            case 8:
                title = Problem08.p08_title();
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("P00_TITLE - Fatal error!");
                Console.WriteLine("  Illegal problem number = " + problem + "");
                return "";
        }

        return title;
    }

    public static double p00_turing(ref p00Data data, int problem, double h, double tol, ref int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P00_TURING applies the Turing quadrature rule.
        //
        //  Discussion
        //
        //    We consider the approximation:
        //
        //      Integral ( -oo < x < +oo ) f(x) dx
        //
        //      = h * Sum ( -oo < i < +oo ) f(nh) + error term
        //
        //    Given H and a tolerance TOL, we start summing at I = 0, and
        //    adding one more term in the positive and negative I directions,
        //    until the absolute value of the next two terms being added 
        //    is less than TOL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Alan Turing,
        //    A Method for the Calculation of the Zeta Function,
        //    Proceedings of the London Mathematical Society,
        //    Volume 48, 1943, pages 180-197.
        //
        //  Parameters:
        //
        //    Input, int PROBLEM, the index of the problem.
        //
        //    Input, double H, the spacing to use.
        //
        //    Input, double TOL, the tolerance.  
        //
        //    Output, int N, the number of pairs of steps taken.
        //    The actual number of function evaluations is 2*N+1.
        //
        //    Output, double P00_TURING, the approximate integral.
        //
    {
        double[] f_vec = new double[2];
        int n_too_many = 100000;
        int option;
        int order;
        double result;
        double[] xtab = new double[2];

        option = 0;
        n = 0;

        result = 0.0;
        order = 1;
        xtab[0] = 0.0;
        p00_fun(ref data, problem, option, order, xtab, ref f_vec);
        result += h * f_vec[0];

        for (;;)
        {
            n += 1;

            xtab[0] = n * h;
            xtab[1] = -(double) n * h;

            order = 2;
            p00_fun(ref data, problem, option, order, xtab, ref f_vec);

            result += h * (f_vec[0] + f_vec[1]);
            //
            //  Just do a simple-minded absolute error tolerance check to start with.
            //
            if (Math.Abs(f_vec[0]) + Math.Abs(f_vec[1]) <= tol)
            {
                break;
            }

            //
            //  Just in case things go crazy.
            //
            if (n_too_many <= n)
            {
                Console.WriteLine("");
                Console.WriteLine("P00_TURING - Warning!");
                Console.WriteLine("  Number of steps exceeded N_TOO_MANY = " + n_too_many + "");
                break;
            }

        }

        return result;
    }
}