using System;

namespace InterpTest;

public static class Data_1D
{
    public static double[] p00_f(int prob, int n, double[] x)

//****************************************************************************80
//
//  Purpose:
//
//    P00_F evaluates the function for any problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the number of the desired test problem.
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P00_F[N], the function values.
//
    {
        double[] f;

        switch (prob)
        {
            case 1:
                f = p01_f(n, x);
                break;
            case 2:
                f = p02_f(n, x);
                break;
            case 3:
                f = p03_f(n, x);
                break;
            case 4:
                f = p04_f(n, x);
                break;
            case 5:
                f = p05_f(n, x);
                break;
            case 6:
                f = p06_f(n, x);
                break;
            case 7:
                f = p07_f(n, x);
                break;
            case 8:
                f = p08_f(n, x);
                break;
            default:
                Console.WriteLine();
                Console.WriteLine("P00_F - Fatal error!");
                Console.WriteLine("  Illegal problem number = " + prob );
                return new double[1];
        }

        return f;
    }

    public static int p00_prob_num()
//****************************************************************************80
//
//  Purpose:
//
//    P00_PROB_NUM returns the number of problems.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P00_PROB_NUM, the number of problems.
//
    {
        const int prob_num = 8;

        return prob_num;
    }

    public static string p00_title(int prob)

//****************************************************************************80
//
//  Purpose:
//
//    P00_TITLE returns the title of any problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the number of the desired test problem.
//
//    Output, string P00_TITLE, the title of the problem.
//
    {
        string title;

        switch (prob)
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
            default:
                Console.WriteLine();
                Console.WriteLine("P00_TITLE - Fatal error!");
                Console.WriteLine("  Illegal problem number = " + prob );
                return "";
        }

        return title;
    }

    public static double[] p01_f(int n, double[] x)
//****************************************************************************80
//
//  Purpose:
//
//    P01_F evaluates the function for problem p01.
//
//  Discussion:
//
//    Step function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P01_F[N], the function values.
//
    {
        double[] f = new double[n];

        for (int i = 0; i < n; i++)
        {
            f[i] = x[i] switch
            {
                <= 1.0 / 3.0 => -1.0,
                <= 4.0 / 5.0 => 2.0,
                _ => 1.0
            };
        }

        return f;
    }

    public static string p01_title()
//****************************************************************************80
//
//  Purpose:
//
//    P01_TITLE returns the title of problem p01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P01_TITLE, the title of the problem.
//
    {
        const string title = "f(x) = steps -1/2/1 at [0,1/3], [1/3,4/5], [4/5,1].";

        return title;
    }

    public static double[] p02_f(int n, double[] x)
//****************************************************************************80
//
//  Purpose:
//
//    P02_F evaluates the function for problem p01.
//
//  Discussion:
//
//    Nondifferentiable function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P02_F[N], the function values.
//
    {
        double[] f = new double[n];

        for (int i = 0; i < n; i++)
        {
            f[i] = x[i] switch
            {
                <= 1.0 / 3.0 => 1.0 - 3.0 * x[i],
                _ => 6.0 * x[i] - 2.0
            };
        }

        return f;
    }

    public static string p02_title()
//****************************************************************************80
//
//  Purpose:
//
//    P02_TITLE returns the title of problem p02.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P02_TITLE, the title of the problem.
//
    {
        const string title = "f(x) = (1-3x), x < 1/3 (6x-2) if 1/3 < x";

        return title;
    }

    public static double[] p03_f(int n, double[] x)
//****************************************************************************80
//
//  Purpose:
//
//    P03_F evaluates the function for problem p03.
//
//  Discussion:
//
//    Step function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P03_F[N], the function values.
//
    {
        double[] f = new double[n];

        for (int i = 0; i < n; i++)
        {
            f[i] = x[i] * (10.0 * x[i] - 1.0)
                        * (5.0 * x[i] - 2.0) * (5.0 * x[i] - 2.0)
                        * (4.0 * x[i] - 3.4) * (x[i] - 1.0);
        }

        return f;
    }

    public static string p03_title()
//****************************************************************************80
//
//  Purpose:
//
//    P03_TITLE returns the title of problem p03.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P03_TITLE, the title of the problem.
//
    {
        const string title = "f(x) = x (10*x-1) (5x-2) (5x-2) (4x-3.4) (x-1)";

        return title;
    }

    public static double[] p04_f(int n, double[] x)
//****************************************************************************80
//
//  Purpose:
//
//    P04_F evaluates the function for problem p04.
//
//  Discussion:
//
//    Step function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P04_F[N], the function values.
//
    {
        double[] f = new double[n];

        for (int i = 0; i < n; i++)
        {
            f[i] = Math.Atan(40.0 * x[i] - 15.0);
        }

        return f;
    }

    public static string p04_title()

//****************************************************************************80
//
//  Purpose:
//
//    P04_TITLE returns the title of problem p04.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P04_TITLE, the title of the problem.
//
    {
        const string title = "f(x) = atan ( 40 * x - 15 )";

        return title;
    }

    public static double[] p05_f(int n, double[] x)
//****************************************************************************80
//
//  Purpose:
//
//    P05_F evaluates the function for problem p05.
//
//  Discussion:
//
//    Step function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P05_F[N], the function values.
//
    {
        double[] f = new double[n];

        for (int i = 0; i < n; i++)
        {
            f[i] = Math.Cos(7.0 * x[i])
                   + 5.0 * Math.Cos(11.2 * x[i])
                   - 2.0 * Math.Cos(14.0 * x[i])
                   + 5.0 * Math.Cos(31.5 * x[i])
                   + 7.0 * Math.Cos(63.0 * x[i]);
        }

        return f;
    }

    public static string p05_title()
//****************************************************************************80
//
//  Purpose:
//
//    P05_TITLE returns the title of problem p05.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P05_TITLE, the title of the problem.
//
    {
        const string title = "f(x) = cos(7*x)+5*cos(11.2*x)-2*cos(14*x)+5*cos(31.5*x)+7*cos(63*x).";

        return title;
    }

    public static double[] p06_f(int n, double[] x)
//****************************************************************************80
//
//  Purpose:
//
//    P06_F evaluates the function for problem p06.
//
//  Discussion:
//
//    f(x) = exp ( - (4 * x - 1)^2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Genz,
//    A Package for Testing Multiple Integration Subroutines,
//    in Numerical Integration: Recent Developments, Software
//    and Applications,
//    edited by Patrick Keast and Graeme Fairweather,
//    D Reidel, 1987, pages 337-340,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P06_F[N], the function values.
//
    {
        double[] f = new double[n];

        for (int i = 0; i < n; i++)
        {
            f[i] = Math.Exp(-Math.Pow(4.0 * x[i] - 1.0, 2));
        }

        return f;
    }

    public static string p06_title()
//****************************************************************************80
//
//  Purpose:
//
//    P06_TITLE returns the title of problem p06.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P06_TITLE, the title of the problem.
//
    {
        const string title = "f(x) = exp ( - ( 4*x-1 )^2 )";

        return title;
    }

    public static double[] p07_f(int n, double[] x)
//****************************************************************************80
//
//  Purpose:
//
//    P07_F evaluates the function for problem p07.
//
//  Discussion:
//
//    f(x) = exp ( 4 * x ) if x <= 1/2
//           0                  otherwise
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Genz,
//    A Package for Testing Multiple Integration Subroutines,
//    in Numerical Integration: Recent Developments, Software
//    and Applications,
//    edited by Patrick Keast and Graeme Fairweather,
//    D Reidel, 1987, pages 337-340,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P07_F[N], the function values.
//
    {
        double[] f = new double[n];

        for (int i = 0; i < n; i++)
        {
            f[i] = x[i] switch
            {
                < 0.5 => Math.Exp(4.0 * x[i]),
                _ => 0.0
            };
        }

        return f;
    }

    public static string p07_title()
//****************************************************************************80
//
//  Purpose:
//
//    P07_TITLE returns the title of problem p07.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P07_TITLE, the title of the problem.
//
    {
        const string title = "f(x) = exp ( 2 x ) if x < 0.5, 0 otherwise";

        return title;
    }

    public static double[] p08_f(int n, double[] x)
//****************************************************************************80
//
//  Purpose:
//
//    P08_F evaluates the function for problem p08.
//
//  Discussion:
//
//    This is a famous example, due to Runge.  If equally spaced
//    abscissas are used, the sequence of interpolating polynomials Pn(X)
//    diverges, in the sense that the max norm of the difference
//    between Pn(X) and F(X) becomes arbitrarily large as N increases.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P08_F[N], the function values.
//
    {
        double[] f = new double[n];

        for (int i = 0; i < n; i++)
        {
            f[i] = 1.0 / (Math.Pow(10.0 * (x[i] - 0.5), 2) + 1.0);
        }

        return f;
    }

    public static string p08_title()
//****************************************************************************80
//
//  Purpose:
//
//    P08_TITLE returns the title of problem p08.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P09_TITLE, the title of the problem.
//
    {
        const string title = "f(x) = 1 / ( 1 + ( 10 * (x-1/2) )^2 )";

        return title;
    }
}