using System;
using Burkardt.Types;

namespace Burkardt.Laguerre;

public static partial class Integrands
{
    public class p00Data
    {
        public p17Data p17data = new();
        public p18Data p18data = new();
        public p19Data p19data = new();
        public p20Data p20data = new();
    }
    public static double p00_alpha(int problem)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P00_ALPHA returns the value of ALPHA for any problem.
        //
        //  Discussion:
        //
        //    ALPHA is the lower, finite limit of integration in the integral.
        //
        //    The typical or default value is 0.0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int PROBLEM, the index of the problem.
        //
        //    Output, double P00_ALPHA, the value of ALPHA.
        //
    {
        double alpha;

        switch (problem)
        {
            case 1:
                alpha = p01_alpha();
                break;
            case 2:
                alpha = p02_alpha();
                break;
            case 3:
                alpha = p03_alpha();
                break;
            case 4:
                alpha = p04_alpha();
                break;
            case 5:
                alpha = p05_alpha();
                break;
            case 6:
                alpha = p06_alpha();
                break;
            case 7:
                alpha = p07_alpha();
                break;
            case 8:
                alpha = p08_alpha();
                break;
            case 9:
                alpha = p09_alpha();
                break;
            case 10:
                alpha = p10_alpha();
                break;
            case 11:
                alpha = p11_alpha();
                break;
            case 12:
                alpha = p12_alpha();
                break;
            case 13:
                alpha = p13_alpha();
                break;
            case 14:
                alpha = p14_alpha();
                break;
            case 15:
                alpha = p15_alpha();
                break;
            case 16:
                alpha = p16_alpha();
                break;
            case 17:
                alpha = p17_alpha();
                break;
            case 18:
                alpha = p18_alpha();
                break;
            case 19:
                alpha = p19_alpha();
                break;
            case 20:
                alpha = p20_alpha();
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("P00_ALPHA - Fatal error!");
                Console.WriteLine("  Illegal problem number = " + problem + "");
                return 1;
        }

        return alpha;
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
        //    27 December 2011
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
                exact = p01_exact();
                break;
            case 2:
                exact = p02_exact();
                break;
            case 3:
                exact = p03_exact();
                break;
            case 4:
                exact = p04_exact();
                break;
            case 5:
                exact = p05_exact();
                break;
            case 6:
                exact = p06_exact();
                break;
            case 7:
                exact = p07_exact();
                break;
            case 8:
                exact = p08_exact();
                break;
            case 9:
                exact = p09_exact();
                break;
            case 10:
                exact = p10_exact();
                break;
            case 11:
                exact = p11_exact();
                break;
            case 12:
                exact = p12_exact();
                break;
            case 13:
                exact = p13_exact();
                break;
            case 14:
                exact = p14_exact();
                break;
            case 15:
                exact = p15_exact();
                break;
            case 16:
                exact = p16_exact();
                break;
            case 17:
                exact = p17_exact();
                break;
            case 18:
                exact = p18_exact(ref data.p18data);
                break;
            case 19:
                exact = p19_exact();
                break;
            case 20:
                exact = p20_exact(ref data.p20data);
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("P00_EXACT - Fatal error!");
                Console.WriteLine("  Illegal problem number = " + problem + "");
                return 1;
        }

        return exact;
    }

    public static double p00_exp_transform(ref p00Data data, int problem, int order)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P00_EXP_TRANSFORM applies an exponential transform and Gauss-Legendre rule.
        //
        //  Discussion:
        //
        //    To approximate:
        //
        //      Integral ( alpha <= x < +oo ) f(x) dx
        //
        //    Transform:
        //
        //      u = exp ( -x )
        //      du = - exp ( -x ) dx
        //
        //      x = - log ( u )
        //      dx = - du / u
        //
        //      x = alpha    => u = exp ( -alpha )
        //      x = Infinity => u = 0
        //
        //    Transformed integral:
        //
        //      Integral ( 0 <= u <= exp ( -alpha ) ) f ( -log(u) ) du / u
        //
        //    We apply a Gauss-Legendre rule here, but we could easily use any rule
        //    that avoids evaluation at U = 0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Arthur Stroud, Don Secrest,
        //    Gaussian Quadrature Formulas,
        //    Prentice Hall, 1966,
        //    LC: QA299.4G3S7.
        //
        //  Parameters:
        //
        //    Input, int PROBLEM, the index of the problem.
        //
        //    Input, int ORDER, the order of the Gauss-Legendre rule 
        //    to apply.
        //
        //    Output, double P00_EXP_TRANSFORM, the approximate integral.
        //
    {
        int i;

        double[] u = new double[order];
        double[] u_log = new double[order];
        double[] weight = new double[order];

        double alpha = p00_alpha(problem);
        //
        //  Get the abscissas and weights for Gauss-Legendre quadrature.
        //
        Legendre.QuadratureRule.legendre_compute(order, ref u, ref weight);
        //
        //  Modify the weights from [-1,1] to [0,exp(-alpha)].
        //
        for (i = 0; i < order; i++)
        {
            weight[i] = Math.Exp(-alpha) * weight[i] / 2.0;
        }

        //
        //  Linear transform of abscissas from [-1,1] to [0,exp(-alpha)].
        //
        for (i = 0; i < order; i++)
        {
            u[i] = ((1.0 + u[i]) * Math.Exp(-alpha)
                    + (1.0 - u[i]) * 0.0)
                   / 2.0;
        }

        //
        //  Define U_LOG = - log ( U )
        //
        for (i = 0; i < order; i++)
        {
            u_log[i] = -Math.Log(u[i]);
        }

        //
        //  Evaluate F ( -LOG(U) ).
        //
        double[] fu = p00_fun(ref data, problem, order, u_log);
        //
        //  The integrand is F ( -LOG(U) ) / U
        //
        for (i = 0; i < order; i++)
        {
            fu[i] /= u[i];
        }

        //
        //  Sum.
        //
        double result = typeMethods.r8vec_dot_product(order, weight, fu);

        return result;
    }

    public static double[] p00_fun(ref p00Data data, int problem, int n, double[] x)

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
        //    27 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int PROBLEM, the index of the problem.
        //
        //    Input, int N, the number of points.
        //
        //    Input, double X[N], the evaluation points.
        //
        //    Output, double P00_FUN[N], the function values.
        //
    {
        double[] f;

        switch (problem)
        {
            case 1:
                f = p01_fun(n, x);
                break;
            case 2:
                f = p02_fun(n, x);
                break;
            case 3:
                f = p03_fun(n, x);
                break;
            case 4:
                f = p04_fun(n, x);
                break;
            case 5:
                f = p05_fun(n, x);
                break;
            case 6:
                f = p06_fun(n, x);
                break;
            case 7:
                f = p07_fun(n, x);
                break;
            case 8:
                f = p08_fun(n, x);
                break;
            case 9:
                f = p09_fun(n, x);
                break;
            case 10:
                f = p10_fun(n, x);
                break;
            case 11:
                f = p11_fun(n, x);
                break;
            case 12:
                f = p12_fun(n, x);
                break;
            case 13:
                f = p13_fun(n, x);
                break;
            case 14:
                f = p14_fun(n, x);
                break;
            case 15:
                f = p15_fun(n, x);
                break;
            case 16:
                f = p16_fun(n, x);
                break;
            case 17:
                f = p17_fun(ref data.p17data, n, x);
                break;
            case 18:
                f = p18_fun(ref data.p18data, n, x);
                break;
            case 19:
                f = p19_fun(ref data.p19data, n, x);
                break;
            case 20:
                f = p20_fun(ref data.p20data, n, x);
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("P00_FUN - Fatal error!");
                Console.WriteLine("  Illegal problem number = " + problem + "");
                return null;
        }

        return f;
    }

    public static double p00_gauss_laguerre(ref p00Data data, int problem, int order)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P00_GAUSS_LAGUERRE applies a Gauss-Laguerre rule.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 July 2007
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
        //    Output, double P00_GAUSS_LAGUERRE, the approximate integral.
        //
    {
        int i;

        double[] weight = new double[order];
        double[] xtab = new double[order];

        double alpha = p00_alpha(problem);

        double alpha2 = 0.0;
        QuadratureRule.laguerre_compute(order, ref xtab, ref weight, alpha2);

        for (i = 0; i < order; i++)
        {
            xtab[i] += alpha;
        }

        double[] fx = p00_fun(ref data, problem, order, xtab);
        //
        //  The Gauss-Laguerre rule assumes a weight of EXP(-X).
        //
        //  We need to multiply each F(X) by EXP(X) to implicitly 
        //  adjust for this weight.
        //
        for (i = 0; i < order; i++)
        {
            fx[i] *= Math.Exp(xtab[i]);
        }

        double result = Math.Exp(-alpha) * typeMethods.r8vec_dot_product(order, weight, fx);

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
        //    27 December 2011
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
        const int problem_num = 20;

        return problem_num;
    }

    public static double p00_rat_transform(ref p00Data data, int problem, int order)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P00_RAT_TRANSFORM applies a rational transform and Gauss-Legendre rule.
        //
        //  Discussion:
        //
        //    To approximate:
        //
        //      Integral ( alpha <= x < +oo ) f(x) dx
        //
        //    Transform:
        //
        //      u = 1 / ( 1 + x )
        //      du = - dx / ( 1 + x )^2
        //
        //      x = ( 1 - u ) / u
        //      dx = - du / u^2
        //
        //      x = alpha    => u = 1 / ( 1 + alpha )
        //      x = Infinity => u = 0
        //
        //    Transformed integral:
        //
        //      Integral ( 0 < u <= 1 / ( 1 + alpha ) ) f ( ( 1 - u ) / u ) du / u^2
        //
        //    We apply a Gauss-Legendre rule here, but we could easily use any rule
        //    that avoids evaluation at U = 0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Arthur Stroud, Don Secrest,
        //    Gaussian Quadrature Formulas,
        //    Prentice Hall, 1966,
        //    LC: QA299.4G3S7.
        //
        //  Parameters:
        //
        //    Input, int PROBLEM, the index of the problem.
        //
        //    Input, int ORDER, the order of the Gauss-Legendre rule 
        //    to apply.
        //
        //    Output, double P00_RAT_TRANSFORM, the approximate integral.
        //
    {
        int i;

        double[] u = new double[order];
        double[] u_rat = new double[order];
        double[] weight = new double[order];

        double alpha = p00_alpha(problem);
        //
        //  Get the abscissas and weights for Gauss-Legendre quadrature.
        //
        Legendre.QuadratureRule.legendre_compute(order, ref u, ref weight);
        //
        //  Modify the weights from [-1,1] to [0,1/(1+alpha)].
        //
        for (i = 0; i < order; i++)
        {
            weight[i] = weight[i] / 2.0 / (1.0 + alpha);
        }

        //
        //  Linear transform of abscissas from [-1,1] to [0,exp(-alpha)].
        //
        for (i = 0; i < order; i++)
        {
            u[i] = ((1.0 + u[i]) / (1.0 + alpha)
                    + (1.0 - u[i]) * 0.0)
                   / 2.0;
        }

        //
        //  Define U_RAT = ( 1 - U ) / U.
        //
        for (i = 0; i < order; i++)
        {
            u_rat[i] = (1.0 - u[i]) / u[i];
        }

        //
        //  Evaluate F ( ( 1 - U ) / U ).
        //
        double[] fu = p00_fun(ref data, problem, order, u_rat);
        //
        //  The integrand is F ( ( 1 - U ) / U ) / U^2
        //
        for (i = 0; i < order; i++)
        {
            fu[i] = fu[i] / u[i] / u[i];
        }

        //
        //  Sum.
        //
        double result = typeMethods.r8vec_dot_product(order, weight, fu);

        return result;
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
        //    27 December 2011
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
                title = p01_title();
                break;
            case 2:
                title = p02_title();
                break;
            case 3:
                title = p03_title();
                break;
            case 4:
                title = p04_title();
                break;
            case 5:
                title = p05_title();
                break;
            case 6:
                title = p06_title();
                break;
            case 7:
                title = p07_title();
                break;
            case 8:
                title = p08_title();
                break;
            case 9:
                title = p09_title();
                break;
            case 10:
                title = p10_title();
                break;
            case 11:
                title = p11_title();
                break;
            case 12:
                title = p12_title();
                break;
            case 13:
                title = p13_title();
                break;
            case 14:
                title = p14_title();
                break;
            case 15:
                title = p15_title();
                break;
            case 16:
                title = p16_title();
                break;
            case 17:
                title = p17_title();
                break;
            case 18:
                title = p18_title();
                break;
            case 19:
                title = p19_title();
                break;
            case 20:
                title = p20_title();
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("P00_TITLE - Fatal error!");
                Console.WriteLine("  Illegal problem number = " + problem + "");
                return "";
        }

        return title;
    }
}