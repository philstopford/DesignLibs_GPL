using System;
using Burkardt.FullertonFnLib;
using Burkardt.Types;
using Burkardt.Uniform;

namespace InterpTest
{
    public static class Data_nD
    {
        public static double csevl(double x, double[] a, int n )
//****************************************************************************80
//
//  Purpose:
//
//    CSEVL evaluates a Chebyshev series.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Roger Broucke,
//    Algorithm 446:
//    Ten Subroutines for the Manipulation of Chebyshev Series,
//    Communications of the ACM,
//    Volume 16, Number 4, April 1973, pages 254-256.
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Input, double CS[N], the Chebyshev coefficients.
//
//    Input, int N, the number of Chebyshev coefficients.
//
//    Output, double CSEVL, the Chebyshev series evaluated at X.
//
        {
            double b2 = 0;

            if (n < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("CSEVL - Fatal error!");
                Console.WriteLine("  Number of terms <= 0.");
                return (1);
            }

            if (1000 < n)
            {
                Console.WriteLine("");
                Console.WriteLine("CSEVL - Fatal error!");
                Console.WriteLine("  Number of terms greater than 1000.");
                return (1);
            }

            if (x < -1.1 || 1.1 < x)
            {
                Console.WriteLine("");
                Console.WriteLine("CSEVL - Fatal error!");
                Console.WriteLine("  X outside (-1,+1).");
                return (1);
            }

            double twox = 2.0 * x;
            double b1 = 0.0;
            double b0 = 0.0;

            for (int i = n - 1; 0 <= i; i--)
            {
                b2 = b1;
                b1 = b0;
                b0 = twox * b1 - b2 + a[i];
            }

            double value = 0.5 * (b0 - b2);

            return value;
        }
        
        public static int inits(double[] dos, int nos, double eta )
//****************************************************************************80
//
//  Purpose:
//
//    INITS initializes a Chebyshev series.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Roger Broucke,
//    Algorithm 446:
//    Ten Subroutines for the Manipulation of Chebyshev Series,
//    Communications of the ACM,
//    Volume 16, Number 4, April 1973, pages 254-256.
//
//  Parameters:
//
//    Input, double DOS[NOS], the Chebyshev coefficients.
//
//    Input, int NOS, the number of coefficients.
//
//    Input, double ETA, the desired accuracy.
//
//    Output, int INITS, the number of terms of the series needed
//    to ensure the requested accuracy.
//
        {
            int i;
            int value;

            if (nos < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("INITS - Fatal error!");
                Console.WriteLine("  Number of coefficients < 1.");
                return (1);
            }

            double err = 0.0;

            for (i = nos - 1; 0 <= i; i--)
            {
                err = err + Math.Abs(dos[i]);
                if (eta < err)
                {
                    value = i + 1;
                    return value;
                }
            }

            value = i;
            Console.WriteLine("");
            Console.WriteLine("INITS - Warning!");
            Console.WriteLine("  ETA may be too small.");

            return value;
        }

        public static double[] p00_c(int prob, int m, ref int seed)
//****************************************************************************80
//
//  Purpose:
//
//    P00_CW computes a random C parameter vector for any problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem number.
//
//    Input, int M, the spatial dimension.
//
//    Input/output, int &SEED, a seed for the random 
//    number generator.
//
//    Output, double P00_C[M], the parameter vector.
//
        {
            double[] b =  {
                1.5, 0.0, 1.85, 7.03, 20.4, 4.3
            }
            ;

            b[2 - 1] = (double) (m);

            double[] c = UniformRNG.r8vec_uniform_01_new(m, ref seed);
            double c_sum = typeMethods.r8vec_sum(m, c);

            for (int i = 0; i < m; i++)
            {
                c[i] = b[prob - 1] * c[i] / c_sum;
            }

            return c;
        }

        public static double[] p00_d(int prob, int m, int id, double[] c, double[] w, int n,
        double[] x )
//****************************************************************************80
//
//  Purpose:
//
//    P00_D returns a derivative component of any function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the index of the function.
//
//    Input, int M, the spatial dimension.
//
//    Input, int ID, the spatial coordinate to differentiate.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[M*N], the evalution points.
//
//    Output, double P00_D[N], the ID-th derivative component.
//
        {
            double[] d;

            if (id < 0 || m < id)
            {
                Console.WriteLine("");
                Console.WriteLine("P00_D - Fatal error!");
                Console.WriteLine("  Illegal spatial coordinate ID = " + id + "");
                return new double[1];
            }

            if (prob == 1)
            {
                d = p01_d(m, id, c, w, n, x);
            }
            else if (prob == 2)
            {
                d = p02_d(m, id, c, w, n, x);
            }
            else if (prob == 3)
            {
                d = p03_d(m, id, c, w, n, x);
            }
            else if (prob == 4)
            {
                d = p04_d(m, id, c, w, n, x);
            }
            else if (prob == 5)
            {
                d = p05_d(m, id, c, w, n, x);
            }
            else if (prob == 6)
            {
                d = p06_d(m, id, c, w, n, x);
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("P00_D - Fatal error!");
                Console.WriteLine("  Illegal function index PROB = " + prob + "");
                return new double[1];
            }

            return d;
        }

        public static double[] p00_f(int prob, int m, double[] c, double[] w, int n, double[] x )
//****************************************************************************80
//
//  Purpose:
//
//    P00_F returns the value of any function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the index of the function.
//
//    Input, int M, the spatial dimension.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[M*N], the evalution points.
//
//    Output, double P00_F[N], the function values.
//
        {
            double[] f;

            if (prob == 1)
            {
                f = p01_f(m, c, w, n, x);
            }
            else if (prob == 2)
            {
                f = p02_f(m, c, w, n, x);
            }
            else if (prob == 3)
            {
                f = p03_f(m, c, w, n, x);
            }
            else if (prob == 4)
            {
                f = p04_f(m, c, w, n, x);
            }
            else if (prob == 5)
            {
                f = p05_f(m, c, w, n, x);
            }
            else if (prob == 6)
            {
                f = p06_f(m, c, w, n, x);
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("P00_F - Fatal error!");
                Console.WriteLine("  Illegal function index PROB = " + prob + "");
                return new double[1];
            }

            return f;
        }

        public static int p00_prob_num()
//****************************************************************************80
//
//  Purpose:
//
//    P00_PROB_NUM returns the number of test functions available.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//   Output, int P00_PROB_NUM, the number of test functions.
//
        {
            int prob_num = 6;

            return prob_num;
        }

        public static double p00_q(ref typeMethods.r8ErrorData data, ref typeMethods.r8ErrorcData cdata, int prob, int m, double[] c, double[] w )
//****************************************************************************80
//
//  Purpose:
//
//    P00_Q returns the integral of any function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the index of the function.
//
//    Input, int M, the spatial dimension.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Output, double P00_Q, the integral.
//
        {
            double q;

            if (prob == 1)
            {
                q = p01_q(m, c, w);
            }
            else if (prob == 2)
            {
                q = p02_q(m, c, w);
            }
            else if (prob == 3)
            {
                q = p03_q(m, c, w);
            }
            else if (prob == 4)
            {
                q = p04_q(ref data, ref cdata, m, c, w);
            }
            else if (prob == 5)
            {
                q = p05_q(m, c, w);
            }
            else if (prob == 6)
            {
                q = p06_q(m, c, w);
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("P00_Q - Fatal error!");
                Console.WriteLine("  Illegal function index PROB = " + prob + "");
                return (1);
            }

            return q;
        }

        public static string p00_title(int prob)
//****************************************************************************80
//
//  Purpose:
//
//    P00_TITLE returns the title for any function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the index of the function.
//
//    Output, string P00_TITLE, the function title.
//
        {
            string title;

            if (prob == 1)
            {
                title = p01_title();
            }
            else if (prob == 2)
            {
                title = p02_title();
            }
            else if (prob == 3)
            {
                title = p03_title();
            }
            else if (prob == 4)
            {
                title = p04_title();
            }
            else if (prob == 5)
            {
                title = p05_title();
            }
            else if (prob == 6)
            {
                title = p06_title();
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("P00_TITLE - Fatal error!");
                Console.WriteLine("  Illegal function index PROB = " + prob + "");
                return "";
            }

            return title;
        }

        public static double[] p00_w(int prob, int m, ref int seed)
//****************************************************************************80
//
//  Purpose:
//
//    P00_W computes a random W parameter vector for any problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem number.
//
//    Input, int M, the spatial dimension.
//
//    Input/output, int &SEED, a seed for the random 
//    number generator.
//
//    Output, double P00_W[M], the parameter vector.
//
        {
            double[] w = UniformRNG.r8vec_uniform_01_new(m, ref seed);

            return w;
        }

        public static double[] p01_d(int m, int id, double[] c, double[] w, int n, double[] x )
//****************************************************************************80
//
//  Purpose:
//
//    P01_D evaluates any derivative component for problem p01.
//
//  Discussion:
//
//    f(x) = cos ( 2 * pi * w(1) + sum ( c(1:m) * x(1:m) ) )
//
//    Default values are:
//
//    c(1:m) = 1/m
//    w(1) = 0.3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
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
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, int ID, the spatial coordinate to differentiate.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Input, int N, the number of points.
//
//    Input, double X[M*N], the evaluation points.
//
//    Output, double P01_D[N], the ID-th derivative component.
//
        {
            double[] d = new double[n];

            for (int j = 0; j < n; j++)
            {
                d[j] = 2.0 * Math.PI * w[0];
            }

            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    d[j] = d[j] + c[i] * x[i + j * m];
                }
            }

            for (int j = 0; j < n; j++)
            {
                d[j] = -c[id] * Math.Sin(d[j]);
            }

            return d;
        }

        public static double[] p01_f(int m, double[] c, double[] w, int n, double[] x )
//****************************************************************************80
//
//  Purpose:
//
//    P01_F evaluates the function for problem p01.
//
//  Discussion:
//
//    f(x) = cos ( 2 * pi * w(1) + sum ( c(1:m) * x(1:m) ) )
//
//    Default values are:
//
//    c(1:m) = 1/m
//    w(1) = 0.3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
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
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Input, int N, the number of points.
//
//    Input, double X[M*N], the evaluation points.
//
//    Output, double P01_F[N], the function values.
//
        {
            double[] f = new double[n];

            for (int j = 0; j < n; j++)
            {
                f[j] = 2.0 * Math.PI * w[0];
            }

            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    f[j] = f[j] + c[i] * x[i + j * m];
                }
            }

            for (int j = 0; j < n; j++)
            {
                f[j] = Math.Cos(f[j]);
            }

            return f;
        }

        public static double p01_q(int m, double[] c, double[] w )
//****************************************************************************80
//
//  Purpose:
//
//    P01_Q evaluates the integral for problem p01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2012
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
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Output, double P01_Q, the integral.
//
        {
            double c_sum = typeMethods.r8vec_sum(m, c);

            double c_prod = 1.0;
            for (int i = 0; i < m; i++)
            {
                c_prod = c_prod * Math.Sin(0.5 * c[i]) / c[i];
            }

            double q = Math.Pow(2.0, m) * Math.Cos(2.0 * Math.PI * w[0] + 0.5 * c_sum) * c_prod;

            return q;
        }

        public static string p01_title()
//****************************************************************************80
//
//  Purpose:
//
//    P01_TITLE returns the name of problem p01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
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
            string title = "Oscillatory";

            return title;
        }

        public static double[] p02_d(int m, int id, double[] c, double[] w, int n, double[] x )

//****************************************************************************80
//
//  Purpose:
//
//    P02_D evaluates an derivative component for problem p02.
//
//  Discussion:
//
//    f(x) = 1 / product ( c(1:m)^(-2) + ( x(1:m) - w(1:m) )^2 )
//
//    Default values are:
//
//    c(1:m) = 1
//    w(1:m) = 0.5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
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
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, int ID, the spatial coordinate to differentiate.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Input, int N, the number of points.
//
//    Input, double X[M*N], the evaluation points.
//
//    Output, double P02_D[N], the ID-th derivative component.
//
        {
            double[] d = new double[n];

            for (int j = 0; j < n; j++)
            {
                d[j] = 1.0;
            }

            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    d[j] = d[j] * (Math.Pow(c[i], -2) + Math.Pow(x[i + j * m] - w[i], 2));
                }
            }

            for (int j = 0; j < n; j++)
            {
                d[j] = -2.0 / d[j] * (x[id + j * m] - w[id]) /
                       (Math.Pow(c[id], -2) + Math.Pow(x[id + j * m] - w[id], 2));
            }

            return d;
        }

        public static double[] p02_f(int m, double[] c, double[] w, int n, double[] x )
//****************************************************************************80
//
//  Purpose:
//
//    P02_F evaluates the function for problem p02.
//
//  Discussion:
//
//    f(x) = 1 / product ( c(1:m)^(-2) + ( x(1:m) - w(1:m) )^2 )
//
//    Default values are:
//
//    c(1:m) = 1
//    w(1:m) = 0.5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
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
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Input, int N, the number of points.
//
//    Input, double X[M*N], the evaluation points.
//
//    Output, double P02_F[N], the function values.
//
        {
            double[] f = new double[n];

            for (int j = 0; j < n; j++)
            {
                f[j] = 1.0;
            }

            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    f[j] = f[j] * (Math.Pow(c[i], -2) + Math.Pow(x[i + j * m] - w[i], 2));
                }
            }

            for (int j = 0; j < n; j++)
            {
                f[j] = 1.0 / f[j];
            }

            return f;
        }

        public static double p02_q(int m, double[] c, double[] w )
//****************************************************************************80
//
//  Purpose:
//
//    P02_Q evaluates the integral for problem p02.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2012
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
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Output, double P02_Q, the integral.
//
        {
            double q = 1.0;

            for (int i = 0; i < m; i++)
            {
                q = q *
                    (Math.Atan((1.0 - w[i]) * c[i])
                     + Math.Atan(w[i] * c[i])
                    ) * c[i];
            }

            return q;
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
//    27 August 2012
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
            string title = "Product Peak";

            return title;
        }

        public static double[] p03_d(int m, int id, double[] c, double[] w, int n, double[] x )
//****************************************************************************80
//
//  Purpose:
//
//    P03_D evaluates any derivative component for problem p03.
//
//  Discussion:
//
//    f(x) = 1 / ( 1 + sum ( c(1:m) * x(1:m) ) ) ^ ( m + 1 )
//
//    Default values are:
//
//    c(1:m) = 1/m
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
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
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, int ID, the spatial coordinate to differentiate.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Input, int N, the number of points.
//
//    Input, double X[M*N], the evaluation points.
//
//    Output, double P03_D[N], the ID-th derivative component.
//
        {
            double[] d = new double[n];

            for (int j = 0; j < n; j++)
            {
                d[j] = 1.0;
            }

            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    d[j] = d[j] + c[i] * x[i + j * m];
                }
            }

            for (int j = 0; j < n; j++)
            {
                d[j] = -c[id] * (double) (m + 1) / Math.Pow(d[j], m + 2);
            }

            return d;
        }

        public static double[] p03_f(int m, double[] c, double[] w, int n, double[] x )
//****************************************************************************80
//
//  Purpose:
//
//    P03_F evaluates the function for problem p03.
//
//  Discussion:
//
//    f(x) = 1 / ( 1 + sum ( c(1:m) * x(1:m) ) ) ^ ( m + 1 )
//
//    Default values are:
//
//    c(1:m) = 1/m
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
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
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Input, int N, the number of points.
//
//    Input, double X[M*N], the evaluation points.
//
//    Output, double P03_F[N], the function values.
//
        {
            double[] f = new double[n];

            for (int j = 0; j < n; j++)
            {
                f[j] = 1.0;
            }

            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    f[j] = f[j] + c[i] * x[i + j * m];
                }
            }

            for (int j = 0; j < n; j++)
            {
                f[j] = 1.0 / Math.Pow(f[j], m + 1);
            }

            return f;
        }

        public static double p03_q(int m, double[] c, double[] w )
//****************************************************************************80
//
//  Purpose:
//
//    P03_Q evaluates the integral for problem p03.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2012
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
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Output, double P03_Q, the integral.
//
        {
            //
//  Here, we need to generate all possible DIM_NUM tuples with
//  values of 0 or 1.
//
            int[] ivec = new int[m];

            double q = 0.0;
            int rank = 0;

            while (true)
            {
                tuple_next(0, 1, m, ref rank, ivec);

                if (rank == 0)
                {
                    break;
                }

                int s = typeMethods.i4vec_sum(m, ivec);

                q = q + typeMethods.r8_mop(s) / (1.0 + typeMethods.r8vec_i4vec_dot_product(m, c, ivec));
            }

            q = q / (typeMethods.r8_factorial(m) * typeMethods.r8vec_product(m, c));

            return q;
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
//    27 August 2012
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
            string title = "Corner Peak";

            return title;
        }

        public static double[] p04_d(int m, int id, double[] c, double[] w, int n, double[] x )
//****************************************************************************80
//
//  Purpose:
//
//    P04_D evaluates any derivative component for problem p04.
//
//  Discussion:
//
//    f(x) = exp ( - sum ( c(1:m)^2 * ( x(1:m) - w(1:m) )^2 )
//
//    Default values are:
//
//    c(1:m) = 1 / m
//    w(1:m) = 0.5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
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
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, int ID, the spatial coordinate to differentiate.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Input, int N, the number of points.
//
//    Input, double X[M*N], the evaluation points.
//
//    Output, double P04_D[N], the ID-th derivative component.
//
        {
            double[] d = new double[n];

            for (int j = 0; j < n; j++)
            {
                d[j] = 0.0;
            }

            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    d[j] = d[j] + Math.Pow(c[i] * (x[i + j * m] - w[i]), 2);
                }
            }

            for (int j = 0; j < n; j++)
            {
                d[j] = Math.Exp(-d[j]) * Math.Pow(c[id], 2) * (-2.0) * (x[id + j * m] - w[id]);
            }

            return d;
        }

        public static double[] p04_f(int m, double[] c, double[] w, int n, double[] x )
//****************************************************************************80
//
//  Purpose:
//
//    P04_F evaluates the function for problem p04.
//
//  Discussion:
//
//    f(x) = exp ( - sum ( c(1:m)^2 * ( x(1:m) - w(1:m) )^2 )
//
//    Default values are:
//
//    c(1:m) = 1 / m
//    w(1:m) = 0.5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
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
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Input, int N, the number of points.
//
//    Input, double X[M*N], the evaluation points.
//
//    Output, double P04_F[N], the function values.
//
        {
            double[] f = new double[n];

            for (int j = 0; j < n; j++)
            {
                f[j] = 0.0;
            }

            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    f[j] = f[j] + Math.Pow(c[i] * (x[i + j * m] - w[i]), 2);
                }
            }

            for (int j = 0; j < n; j++)
            {
                f[j] = Math.Exp(-f[j]);
            }

            return f;
        }

        public static double p04_q(ref typeMethods.r8ErrorData data, ref typeMethods.r8ErrorcData cdata, int m, double[] c, double[] w )
//****************************************************************************80
//
//  Purpose:
//
//    P04_Q evaluates the integral for problem p04.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2012
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
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Output, double P04_Q, the integral.
//
        {
            double q = 1.0;

            for (int i = 0; i < m; i++)
            {
                q = q * Math.Sqrt(Math.PI)
                      * (typeMethods.r8_error(ref data, ref cdata, c[i] * (1.0 - w[i]))
                         + typeMethods.r8_error( ref data, ref cdata,c[i] * w[i]))
                    / (2.0 * c[i]);
            }

            return q;
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
//    27 August 2012
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
            string title = "Gaussian";

            return title;
        }

        public static double[] p05_d(int m, int id, double[] c, double[] w, int n, double[] x )
//****************************************************************************80
//
//  Purpose:
//
//    P05_D evaluates any derivative component for problem p05.
//
//  Discussion:
//
//    f(x) = exp ( - sum ( c(1:m) * abs ( x(1:m) - w(1:m) ) ) )
//
//    Default values are:
//
//    c(1:m) = 2.0
//    w(1:m) = 0.5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
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
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, int ID, the spatial coordinate to differentiate.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Input, int N, the number of points.
//
//    Input, double X[M*N], the evaluation points.
//
//    Output, double P05_D[N], the ID-th derivative component.
//
        {
            double[] d = new double[n];

            for (int j = 0; j < n; j++)
            {
                d[j] = 0.0;
            }

            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    d[j] = d[j] + c[i] * Math.Abs(x[i + j * m] - w[i]);
                }
            }

            for (int j = 0; j < n; j++)
            {
                d[j] = Math.Exp(-d[j]);
            }

            for (int j = 0; j < n; j++)
            {
                if (x[id + j * m] - w[id] <= 0.0)
                {
                    d[j] = d[j] * c[id];
                }
                else
                {
                    d[j] = -d[j] * c[id];
                }
            }

            return d;
        }

        public static double[] p05_f(int m, double[] c, double[] w, int n, double[] x )
//****************************************************************************80
//
//  Purpose:
//
//    P05_F evaluates the function for problem p05.
//
//  Discussion:
//
//    f(x) = exp ( - sum ( c(1:m) * abs ( x(1:m) - w(1:m) ) ) )
//
//    Default values are:
//
//    c(1:m) = 2.0
//    w(1:m) = 0.5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
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
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Input, int N, the number of points.
//
//    Input, double X[M*N], the evaluation points.
//
//    Output, double P05_F[N], the function values.
//
        {
            double[] f = new double[n];

            for (int j = 0; j < n; j++)
            {
                f[j] = 0.0;
            }

            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    f[j] = f[j] + c[i] * Math.Abs(x[i + j * m] - w[i]);
                }
            }

            for (int j = 0; j < n; j++)
            {
                f[j] = Math.Exp(-f[j]);
            }

            return f;
        }

        public static double p05_q(int m, double[] c, double[] w )
//****************************************************************************80
//
//  Purpose:
//
//    P05_Q evaluates the integral for problem p05.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2012
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
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Output, double P05_Q, the integral.
//
        {
            double q = 1.0;

            for (int i = 0; i < m; i++)
            {
//
//  W < 0 < 1
//
//  | X - W | = X - W from 0 to 1.
//
                if (w[i] < 0.0)
                {
                    q = q *
                        (Math.Exp(-c[i] * (-w[i]))
                         - Math.Exp(-c[i] * (1.0 - w[i]))) / c[i];
                }
//
//  0 < W < 1
//
//  | X - W | = W - X from 0 to Z, 
//            = X - W from      Z to 1.
//
                else if (w[i] < 1.0)
                {
                    q = q * (2.0
                             - Math.Exp(-c[i] * (w[i]))
                             - Math.Exp(-c[i] * (1.0 - w[i]))) / c[i];
                }
//
//  0 < 1 < W
//
//  | X - W | = W - X from 0 to 1.
//
                else
                {
                    q = q *
                        (Math.Exp(-c[i] * (w[i] - 1.0))
                         - Math.Exp(-c[i] * (w[i]))) / c[i];

                }
            }

            return q;
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
//    27 August 2012
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
            string title = "Continuous";

            return title;
        }

        public static double[] p06_d(int m, int id, double[] c, double[] w, int n, double[] x )
//****************************************************************************80
//
//  Purpose:
//
//    P06_D evaluates any derivative component for problem p06.
//
//  Discussion:
//
//    f(x) = exp ( c(1:m) * x(1:m) ) if x(1) <= w(1) and x(2) <= w(2).
//           0                          otherwise
//
//    Default values are:
//
//    c(1:m) = 0.5^(1/m)
//    w(1:2) = 0.5^(1/m)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
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
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, int ID, the spatial coordinate to differentiate.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Input, int N, the number of points.
//
//    Input, double X[M*N], the evaluation points.
//
//    Output, double P06_D[N], the ID-th derivative component.
//
        {
            double[] d = new double[n];

            if (m == 1)
            {
                for (int j = 0; j < n; j++)
                {
                    d[j] = c[0] * Math.Exp(c[0] * x[0 + j * m]);
                }

                for (int j = 0; j < n; j++)
                {
                    if (w[0] < x[0 + j * m])
                    {
                        d[j] = 0.0;
                    }
                }
            }
            else
            {
                for (int j = 0; j < n; j++)
                {
                    d[j] = 0.0;
                }

                for (int i = 0; i < m; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        d[j] = d[j] + c[i] * x[i + j * m];
                    }
                }

                for (int j = 0; j < n; j++)
                {
                    d[j] = c[id] * Math.Exp(d[j]);
                }

                for (int j = 0; j < n; j++)
                {
                    if (w[0] < x[0 + j * m] || w[1] < x[1 + j * m])
                    {
                        d[j] = 0.0;
                    }
                }
            }

            return d;
        }

        public static double[] p06_f(int m, double[] c, double[] w, int n, double[] x )
//****************************************************************************80
//
//  Purpose:
//
//    P06_F evaluates the function for problem p06.
//
//  Discussion:
//
//    f(x) = exp ( c(1:m) * x(1:m) ) if x(1) <= w(1) and x(2) <= w(2).
//           0                          otherwise
//
//    Default values are:
//
//    c(1:m) = 0.5^(1/m)
//    w(1:2) = 0.5^(1/m)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
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
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Input, int N, the number of points.
//
//    Input, double X[M*N], the evaluation points.
//
//    Output, double P06_F[N], the function values.
//
        {
            double[] f = new double[n];

            if (m == 1)
            {
                for (int j = 0; j < n; j++)
                {
                    f[j] = Math.Exp(c[0] * x[0 + j * m]);
                }

                for (int j = 0; j < n; j++)
                {
                    if (w[0] < x[0 + j * m])
                    {
                        f[j] = 0.0;
                    }
                }
            }
            else
            {
                for (int j = 0; j < n; j++)
                {
                    f[j] = 0.0;
                }

                for (int i = 0; i < m; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        f[j] = f[j] + c[i] * x[i + j * m];
                    }
                }

                for (int j = 0; j < n; j++)
                {
                    f[j] = Math.Exp(f[j]);
                }

                for (int j = 0; j < n; j++)
                {
                    if (w[0] < x[0 + j * m] || w[1] < x[1 + j * m])
                    {
                        f[j] = 0.0;
                    }
                }
            }

            return f;
        }

        public static double p06_q(int m, double[] c, double[] w )
//****************************************************************************80
//
//  Purpose:
//
//    P06_Q evaluates the integral for problem p06.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2012
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
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Output, double P06_Q, the integral.
//
        {
            //
//  To simplify the calculation, force W(3:M) to be at least 1.0.
//
            for (int i = 2; i < m; i++)
            {
                w[i] = 1.0;
            }

            double q = 1.0;

            for (int i = 0; i < m; i++)
            {
                if (w[i] <= 0.0)
                {
                    q = q * 0.0;
                }
                else if (w[i] <= 1.0)
                {
                    if (c[i] == 0.0)
                    {
                        q = q * w[i];
                    }
                    else
                    {
                        q = q * (Math.Exp(c[i] * w[i]) - 1.0) / c[i];
                    }
                }
                else
                {
                    if (c[i] != 0.0)
                    {
                        q = q * (Math.Exp(c[i] * w[i]) - 1.0) / c[i];
                    }
                }
            }

            return q;
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
//    27 August 2012
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
            string title = "Discontinuous";

            return title;
        }
        
        public static void tuple_next(int m1, int m2, int n, ref int rank, int[] x )

//****************************************************************************80
//
//  Purpose:
//
//    TUPLE_NEXT computes the next element of a tuple space.
//
//  Discussion:
//
//    The elements are N vectors.  Each entry is constrained to lie
//    between M1 and M2.  The elements are produced one at a time.
//    The first element is
//      (M1,M1,...,M1),
//    the second element is
//      (M1,M1,...,M1+1),
//    and the last element is
//      (M2,M2,...,M2)
//    Intermediate elements are produced in lexicographic order.
//
//  Example:
//
//    N = 2, M1 = 1, M2 = 3
//
//    INPUT        OUTPUT
//    -------      -------
//    Rank  X      Rank   X
//    ----  ---    -----  ---
//    0     * *    1      1 1
//    1     1 1    2      1 2
//    2     1 2    3      1 3
//    3     1 3    4      2 1
//    4     2 1    5      2 2
//    5     2 2    6      2 3
//    6     2 3    7      3 1
//    7     3 1    8      3 2
//    8     3 2    9      3 3
//    9     3 3    0      0 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M1, M2, the minimum and maximum entries.
//
//    Input, int N, the number of components.
//
//    Input/output, int &RANK, counts the elements.
//    On first call, set RANK to 0.  Thereafter, the output value of RANK
//    will indicate the order of the element returned.  When there are no
//    more elements, RANK will be returned as 0.
//
//    Input/output, int X[N], on input the previous tuple.
//    On output, the next tuple.
//
        {
            int i;
            int j;

            if (m2 < m1)
            {
                rank = 0;
                return;
            }

            if (rank <= 0)
            {
                for (i = 0; i < n; i++)
                {
                    x[i] = m1;
                }

                rank = 1;
            }
            else
            {
                rank = rank + 1;
                i = n - 1;

                for (;;)
                {

                    if (x[i] < m2)
                    {
                        x[i] = x[i] + 1;
                        break;
                    }

                    x[i] = m1;

                    if (i == 0)
                    {
                        rank = 0;
                        for (j = 0; j < n; j++)
                        {
                            x[j] = m1;
                        }

                        break;
                    }

                    i = i - 1;
                }
            }
        }
    }
}