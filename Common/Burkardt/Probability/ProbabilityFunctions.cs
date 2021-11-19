using System;
using Burkardt.Types;

namespace Burkardt.Probability;

public static class ProbabilityFunctions
{
    public static double[] p00_a(int prob, int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P00_A returns the matrix A for any least squares problem.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int PROB, the problem index.
        //
        //    Input, int M, the number of equations.
        //
        //    Input, int N, the number of variables.
        //
        //    Output, double P00_A[M*N], the matrix.
        //
    {
        double[] a;

        switch (prob)
        {
            case 1:
                a = p01_a(m, n);
                break;
            case 2:
                a = p02_a(m, n);
                break;
            case 3:
                a = p03_a(m, n);
                break;
            case 4:
                a = p04_a(m, n);
                break;
            case 5:
                a = p05_a(m, n);
                break;
            case 6:
                a = p06_a(m, n);
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("P00_A - Fatal error!");
                Console.WriteLine("  Illegal value of PROB = " + prob + "");
                return null;
        }

        return a;
    }

    public static double[] p00_b(int prob, int m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P00_B returns the right hand side B for any least squares problem.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int PROB, the problem index.
        //
        //    Input, int M, the number of equations.
        //
        //    Output, double P00_B[M], the right hand side.
        //
    {
        double[] b;

        switch (prob)
        {
            case 1:
                b = p01_b(m);
                break;
            case 2:
                b = p02_b(m);
                break;
            case 3:
                b = p03_b(m);
                break;
            case 4:
                b = p04_b(m);
                break;
            case 5:
                b = p05_b(m);
                break;
            case 6:
                b = p06_b(m);
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("P00_B - Fatal error!");
                Console.WriteLine("  Illegal value of PROB = " + prob + "");
                return null;
        }

        return b;
    }

    public static int p00_m(int prob)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P00_M returns the number of equations M for any least squares problem.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int PROB, the problem index.
        //
        //    Output, int P00_M, the number of equations.
        //
    {
        int m;

        switch (prob)
        {
            case 1:
                m = p01_m();
                break;
            case 2:
                m = p02_m();
                break;
            case 3:
                m = p03_m();
                break;
            case 4:
                m = p04_m();
                break;
            case 5:
                m = p05_m();
                break;
            case 6:
                m = p06_m();
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("P00_M - Fatal error!");
                Console.WriteLine("  Illegal value of PROB = " + prob + "");
                return 0;
        }

        return m;
    }

    public static int p00_n(int prob)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P00_N returns the number of variables N for any least squares problem.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int PROB, the problem index.
        //
        //    Output, int P00_N, the number of variables.
        //
    {
        int n;

        switch (prob)
        {
            case 1:
                n = p01_n();
                break;
            case 2:
                n = p02_n();
                break;
            case 3:
                n = p03_n();
                break;
            case 4:
                n = p04_n();
                break;
            case 5:
                n = p05_n();
                break;
            case 6:
                n = p06_n();
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("P00_N - Fatal error!");
                Console.WriteLine("  Illegal value of PROB = " + prob + "");
                return 0;
        }

        return n;
    }

    public static int p00_prob_num()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P00_PROB_NUM returns the number of least squares problems.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double P00_PROB_NUM, the number of problems.
        //
    {
        const int prob_num = 6;

        return prob_num;
    }

    public static double[] p00_x(int prob, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P00_X returns the least squares solution X for any least squares problem.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int PROB, the problem index.
        //
        //    Input, int N, the number of variables.
        //
        //    Output, double P00_X[N], the least squares solution.
        //
    {
        double[] x;

        switch (prob)
        {
            case 1:
                x = p01_x(n);
                break;
            case 2:
                x = p02_x(n);
                break;
            case 3:
                x = p03_x(n);
                break;
            case 4:
                x = p04_x(n);
                break;
            case 5:
                x = p05_x(n);
                break;
            case 6:
                x = p06_x(n);
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("P00_X - Fatal error!");
                Console.WriteLine("  Illegal value of PROB = " + prob + "");
                return null;
        }

        return x;
    }

    public static double[] p01_a(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P01_A returns the matrix A for problem 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of equations.
        //
        //    Input, int N, the number of variables.
        //
        //    Output, double P01_A[M*N], the matrix.
        //
    {
        double[] a = new double[m * n];

        for (int i = 0; i < m; i++)
        {
            a[i + 0 * m] = 1.0;
            for (int j = 1; j < n; j++)
            {
                a[i + j * m] = a[i + (j - 1) * m] * (i + 1);
            }
        }

        return a;
    }

    public static double[] p01_b(int m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P01_B returns the right hand side B for problem 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of equations.
        //
        //    Output, double P01_B[M], the right hand side.
        //
    {
        double[] b_save =  {
                1.0, 2.3, 4.6, 3.1, 1.2
            }
            ;

        double[] b = typeMethods.r8vec_copy_new(m, b_save);

        return b;
    }

    public static int p01_m()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P01_M returns the number of equations M for problem 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int P01_M, the number of equations.
        //
    {
        const int m = 5;

        return m;
    }

    public static int p01_n()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P01_N returns the number of variables N for problem 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int P01_N, the number of variables.
        //
    {
        const int n = 3;

        return n;
    }

    public static double[] p01_x(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P01_X returns the least squares solution X for problem 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of variables.
        //
        //    Output, double P01_X[N], the least squares solution.
        //
    {
        double[] x_save =  {
                -3.0200000, 4.4914286, -0.72857143
            }
            ;

        double[] x = typeMethods.r8vec_copy_new(n, x_save);

        return x;
    }

    public static double[] p02_a(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P02_A returns the matrix A for problem 2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Cleve Moler,
        //    Numerical Computing with MATLAB,
        //    SIAM, 2004,
        //    ISBN13: 978-0-898716-60-3,
        //    LC: QA297.M625,
        //    ebook: http://www.mathworks.com/moler/chapters.html
        //
        //  Parameters:
        //
        //    Input, int M, the number of equations.
        //
        //    Input, int N, the number of variables.
        //
        //    Output, double P02_A[M*N], the matrix.
        //
    {
        double[] a = new double[m * n];

        for (int i = 0; i < m; i++)
        {
            a[i + (n - 1) * m] = 1.0;
            for (int j = n - 2; 0 <= j; j--)
            {
                a[i + j * m] = a[i + (j + 1) * m] * i / 5.0;
            }
        }

        return a;
    }

    public static double[] p02_b(int m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P02_B returns the right hand side B for problem 2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of equations.
        //
        //    Output, double P02_B[M], the right hand side.
        //
    {
        double[] b_save =  {
                150.697, 179.323, 203.212, 226.505, 249.633, 281.422
            }
            ;

        double[] b = typeMethods.r8vec_copy_new(m, b_save);

        return b;
    }

    public static int p02_m()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P02_M returns the number of equations M for problem 2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int P02_M, the number of equations.
        //
    {
        const int m = 6;

        return m;
    }

    public static int p02_n()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P02_N returns the number of variables N for problem 2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int P02_N, the number of variables.
        //
    {
        const int n = 3;

        return n;
    }

    public static double[] p02_x(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P02_X returns the least squares solution X for problem 2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of variables.
        //
        //    Output, double P02_X[N], the least squares solution.
        //
    {
        double[] x_save =  {
                5.7013, 121.1341, 152.4745
            }
            ;

        double[] x = typeMethods.r8vec_copy_new(n, x_save);

        return x;
    }

    public static double[] p03_a(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P03_A returns the matrix A for problem 3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Cleve Moler,
        //    Numerical Computing with MATLAB,
        //    SIAM, 2004,
        //    ISBN13: 978-0-898716-60-3,
        //    LC: QA297.M625,
        //    ebook: http://www.mathworks.com/moler/chapters.html
        //
        //  Parameters:
        //
        //    Input, int M, the number of equations.
        //
        //    Input, int N, the number of variables.
        //
        //    Output, double P03_A[M*N], the matrix.
        //
    {
        double[] a_save =  {
                1.0, 4.0, 7.0, 10.0, 13.0,
                2.0, 5.0, 8.0, 11.0, 14.0,
                3.0, 6.0, 9.0, 12.0, 15.0
            }
            ;

        double[] a = typeMethods.r8mat_copy_new(m, n, a_save);

        return a;
    }

    public static double[] p03_b(int m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P03_B returns the right hand side B for problem 3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of equations.
        //
        //    Output, double P03_B[M], the right hand side.
        //
    {
        double[] b_save =  {
                16.0, 17.0, 18.0, 19.0, 20.0
            }
            ;

        double[] b = typeMethods.r8vec_copy_new(m, b_save);

        return b;
    }

    public static int p03_m()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P03_M returns the number of equations M for problem 3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int P03_M, the number of equations.
        //
    {
        const int m = 5;

        return m;
    }

    public static int p03_n()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P03_N returns the number of variables N for problem 3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int P03_N, the number of variables.
        //
    {
        const int n = 3;

        return n;
    }

    public static double[] p03_x(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P03_X returns the least squares solution X for problem 3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of variables.
        //
        //    Output, double P03_X[N], the least squares solution.
        //
    {
        double[] x_save =  {
                -7.5555556, 0.1111111, 7.7777778
            }
            ;

        double[] x = typeMethods.r8vec_copy_new(n, x_save);

        return x;
    }

    public static double[] p04_a(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P04_A returns the matrix A for problem 4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of equations.
        //
        //    Input, int N, the number of variables.
        //
        //    Output, double P04_A[M*N], the matrix.
        //
    {
        double[] a = new double[m * n];

        for (int j = 0; j < n; j++)
        {
            for (int i = 0; i < m; i++)
            {
                a[i + j * m] = Math.Pow(j + 1, i);
            }
        }

        return a;
    }

    public static double[] p04_b(int m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P04_B returns the right hand side B for problem 4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of equations.
        //
        //    Output, double P04_B[M], the right hand side.
        //
    {
        double[] b_save =  {
                15.0, 55.0, 225.0
            }
            ;

        double[] b = typeMethods.r8vec_copy_new(m, b_save);

        return b;
    }

    public static int p04_m()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P04_M returns the number of equations M for problem 4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int P04_M, the number of equations.
        //
    {
        int m = 3;

        return m;
    }

    public static int p04_n()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P04_N returns the number of variables N for problem 4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int P04_N, the number of variables.
        //
    {
        const int n = 5;

        return n;
    }

    public static double[] p04_x(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P04_X returns the least squares solution X for problem 4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of variables.
        //
        //    Output, double P04_X[N], the least squares solution.
        //
    {
        double[] x_save =  {
                1.0, 2.0, 3.0, 4.0, 5.0
            }
            ;

        double[] x = typeMethods.r8vec_copy_new(n, x_save);

        return x;
    }

    public static double[] p05_a(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P05_A returns the matrix A for problem 5.
        //
        //  Discussion:
        //
        //    A is the Hilbert matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of equations.
        //
        //    Input, int N, the number of variables.
        //
        //    Output, double P05_A[M*N], the matrix.
        //
    {
        double[] a = new double[m * n];

        for (int j = 0; j < n; j++)
        {
            for (int i = 0; i < m; i++)
            {
                a[i + j * m] = 1.0 / (i + j + 1);
            }
        }

        return a;
    }

    public static double[] p05_b(int m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P05_B returns the right hand side B for problem 5.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of equations.
        //
        //    Output, double P05_B[M], the right hand side.
        //
    {
        double[] b = new double[m];

        b[0] = 1.0;
        for (int i = 1; i < m; i++)
        {
            b[i] = 0.0;
        }

        return b;
    }

    public static int p05_m()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P05_M returns the number of equations M for problem 5.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int P05_M, the number of equations.
        //
    {
        const int m = 10;

        return m;
    }

    public static int p05_n()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P05_N returns the number of variables N for problem 5.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int P05_N, the number of variables.
        //
    {
        const int n = 10;

        return n;
    }

    public static double[] p05_x(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P05_X returns the least squares solution X for problem 5.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of variables.
        //
        //    Output, double P05_X[N], the least squares solution.
        //
    {
        double[] x = new double[n];

        for (int i = 0; i < n; i++)
        {
            x[i] = typeMethods.r8_mop(i + 2) * (i + 1)
                                             * typeMethods.r8_choose(n + i, n - 1) * typeMethods.r8_choose(n, n - i - 1);
        }

        return x;
    }

    public static double[] p06_a(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P06_A returns the matrix A for problem 6.
        //
        //  Discussion:
        //
        //    A is a symmetric, orthogonal matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of equations.
        //
        //    Input, int N, the number of variables.
        //
        //    Output, double P06_A[M*N], the matrix.
        //
    {
        double[] a = new double[m * n];

        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                double angle = (i + 1) * (j + 1) * Math.PI / (n + 1);
                a[i + j * m] = Math.Sin(angle) * Math.Sqrt(2.0 / (n + 1));
            }
        }

        return a;
    }

    public static double[] p06_b(int m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P06_B returns the right hand side B for problem 6.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of equations.
        //
        //    Output, double P06_B[M], the right hand side.
        //
    {
        double[] b = new double[m];

        b[0] = 1.0;
        for (int i = 1; i < m; i++)
        {
            b[i] = 0.0;
        }

        return b;
    }

    public static int p06_m()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P06_M returns the number of equations M for problem 6.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int P06_M, the number of equations.
        //
    {
        const int m = 10;

        return m;
    }

    public static int p06_n()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P06_N returns the number of variables N for problem 6.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int P06_N, the number of variables.
        //
    {
        const int n = 10;

        return n;
    }

    public static double[] p06_x(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P06_X returns the least squares solution X for problem 6.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of variables.
        //
        //    Output, double P06_X[N], the least squares solution.
        //
    {
        double[] x = new double[n];

        for (int i = 0; i < n; i++)
        {
            double angle = (i + 1) * Math.PI / (n + 1);
            x[i] = Math.Sin(angle) * Math.Sqrt(2.0 / (n + 1));
        }

        return x;
    }
}