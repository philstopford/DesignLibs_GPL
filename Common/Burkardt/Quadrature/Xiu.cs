using System;
using Burkardt.IntegralNS;
using Burkardt.Types;

namespace Burkardt.Quadrature;

public static class Xiu
{
    public static void cn_leg_02_xiu(int n, int o, ref double[] x, ref double[] w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CN_LEG_02_XIU implements the Xiu precision 2 rule for region CN_LEG.
        //
        //  Discussion:
        //
        //    The rule has order 
        //
        //      O = N + 1.
        //
        //    The rule has precision P = 2.
        //
        //    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
        //
        //      w(x) = 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 March 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dongbin Xiu,
        //    Numerical integration formulas of degree two,
        //    Applied Numerical Mathematics,
        //    Volume 58, 2008, pages 1515-1520.
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, int O, the order.
        //
        //    Output, double X[N*O], the abscissas.
        //
        //    Output, double W[O], the weights.
        //
    {
        double arg;
        double c1;
        double delta0;
        int expon;
        double gamma0;
        int i;
        int j;
            
        int r;
        double volume;
        double volume_1d;

        for (j = 0; j < o; j++)
        {
            i = 0;
            for (r = 1; r <= n / 2; r++)
            {
                arg = 2 * r * j * Math.PI / (n + 1);

                x[i + j * n] = Math.Sqrt(2.0) * Math.Cos(arg);
                i += 1;
                x[i + j * n] = Math.Sqrt(2.0) * Math.Sin(arg);
                i += 1;
            }

            if (i < n)
            {
                x[i + j * n] = typeMethods.r8_mop(j);
                i += 1;
            }
        }

        gamma0 = 1.0;
        delta0 = 0.0;
        c1 = 1.0 / 3.0;

        for (j = 0; j < o; j++)
        {
            for (i = 0; i < n; i++)
            {
                x[i + j * n] = (Math.Sqrt(gamma0 * c1) * x[i + j * n] - delta0) / gamma0;
            }
        }

        expon = 0;
        volume_1d = C1.c1_leg_monomial_integral(expon);
        volume = Math.Pow(volume_1d, n);

        for (j = 0; j < o; j++)
        {
            w[j] = volume / o;
        }
    }

    public static int cn_leg_02_xiu_size(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CN_LEG_02_XIU_SIZE sizes the Xiu rule for region CN_LEG.
        //
        //  Discussion:
        //
        //    The rule has order 
        //
        //      O = N + 1.
        //
        //    The rule has precision P = 2.
        //
        //    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
        //
        //      w(x) = 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 February 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dongbin Xiu,
        //    Numerical integration formulas of degree two,
        //    Applied Numerical Mathematics,
        //    Volume 58, 2008, pages 1515-1520.
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Output, int CN_LEG_02_XIU_SIZE, the order.
        //
    {
        int o;

        o = n + 1;

        return o;
    }



    public static void cn_leg_03_xiu(int n, int o, ref double[] x, ref double[] w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CN_LEG_03_XIU implements the Xiu precision 3 rule for region CN_LEG.
        //
        //  Discussion:
        //
        //    The rule has order 
        //
        //      O = 2 * N.
        //
        //    The rule has precision P = 3.
        //
        //    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
        //
        //      w(x) = 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 February 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dongbin Xiu,
        //    Numerical integration formulas of degree two,
        //    Applied Numerical Mathematics,
        //    Volume 58, 2008, pages 1515-1520.
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, int O, the order.
        //
        //    Output, double X[N*O], the abscissas.
        //
        //    Output, double W[O], the weights.
        //
    {
        double arg;
        int expon;
        int i;
        int j;
            
        int r;
        double volume;

        expon = 0;
        volume = C1.c1_leg_monomial_integral(expon);
        volume = Math.Pow(volume, n);

        for (j = 0; j < o; j++)
        {
            i = 0;
            for (r = 1; r <= n / 2; r++)
            {
                arg = (2 * r - 1) * (j + 1) * Math.PI / n;

                x[i + j * n] = Math.Sqrt(2.0) * Math.Cos(arg) / Math.Sqrt(3.0);
                i += 1;
                x[i + j * n] = Math.Sqrt(2.0) * Math.Sin(arg) / Math.Sqrt(3.0);
                i += 1;
            }

            if (i < n)
            {
                x[i + j * n] = n switch
                {
                    1 => typeMethods.r8_mop(j + 1) / Math.Sqrt(3.0),
                    _ => Math.Sqrt(2.0) * typeMethods.r8_mop(j + 1) / Math.Sqrt(3.0)
                };

                i += 1;
            }
        }

        for (j = 0; j < o; j++)
        {
            w[j] = volume / o;
        }
    }

    public static int cn_leg_03_xiu_size(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CN_LEG_03_XIU_SIZE sizes the Xiu precision 3 rule for region CN_LEG.
        //
        //  Discussion:
        //
        //    The rule has order 
        //
        //      O = 2 * N.
        //
        //    The rule has precision P = 3.
        //
        //    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
        //
        //      w(x) = 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 February 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dongbin Xiu,
        //    Numerical integration formulas of degree two,
        //    Applied Numerical Mathematics,
        //    Volume 58, 2008, pages 1515-1520.
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Output, int CN_LEG_03_XIU_SIZE, the order.
        //
    {
        int o;

        o = 2 * n;

        return o;
    }



    public static void cn_geg_02_xiu(int n, double alpha, int o, ref double[] x, ref double[] w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CN_GEG_02_XIU implements the Xiu precision 2 rule for region CN_GEG.
        //
        //  Discussion:
        //
        //    The rule has order 
        //
        //      O = N + 1.
        //
        //    The rule has precision P = 2.
        //
        //    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
        //
        //      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
        //
        //    with -1.0 < alpha.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 March 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dongbin Xiu,
        //    Numerical integration formulas of degree two,
        //    Applied Numerical Mathematics,
        //    Volume 58, 2008, pages 1515-1520.
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, double ALPHA, the parameter.
        //    -1.0 < ALPHA.
        //
        //    Input, int O, the order.
        //
        //    Output, double X[N*O], the abscissas.
        //
        //    Output, double W[O], the weights.
        //
    {
        double arg;
        double c1;
        double delta0;
        int expon;
        double gamma0;
        int i;
        int j;
            
        int r;
        double volume;
        double volume_1d;

        switch (alpha)
        {
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("CN_GEG_02_XIU - Fatal error!");
                Console.WriteLine("  ALPHA <= -1.0");
                return;
        }

        for (j = 0; j < o; j++)
        {
            i = 0;
            for (r = 1; r <= n / 2; r++)
            {
                arg = 2 * r * j * Math.PI / (n + 1);

                x[i + j * n] = Math.Sqrt(2.0) * Math.Cos(arg);
                i += 1;
                x[i + j * n] = Math.Sqrt(2.0) * Math.Sin(arg);
                i += 1;
            }

            if (i < n)
            {
                x[i + j * n] = typeMethods.r8_mop(j);
                i += 1;
            }
        }

        gamma0 = 1.0;
        delta0 = 0.0;
        c1 = 1.0 / (2.0 * alpha + 3.0);

        for (j = 0; j < o; j++)
        {
            for (i = 0; i < n; i++)
            {
                x[i + j * n] = (Math.Sqrt(gamma0 * c1) * x[i + j * n] - delta0) / gamma0;
            }
        }

        expon = 0;
        volume_1d = C1.c1_geg_monomial_integral(alpha, expon);
        volume = Math.Pow(volume_1d, n);

        for (j = 0; j < o; j++)
        {
            w[j] = volume / o;
        }
    }

    public static int cn_geg_02_xiu_size(int n, double alpha)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CN_GEG_02_XIU_SIZE sizes the Xiu rule for region CN_GEG.
        //
        //  Discussion:
        //
        //    The rule has order 
        //
        //      O = N + 1.
        //
        //    The rule has precision P = 2.
        //
        //    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
        //
        //      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
        //
        //    with -1.0 < alpha.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dongbin Xiu,
        //    Numerical integration formulas of degree two,
        //    Applied Numerical Mathematics,
        //    Volume 58, 2008, pages 1515-1520.
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, double ALPHA, the parameter.
        //    -1.0 < ALPHA.
        //
        //    Output, int CN_GEG_02_XIU_SIZE, the order.
        //
    {
        int o;

        switch (alpha)
        {
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("CN_GEG_02_XIU_SIZE - Fatal error!");
                Console.WriteLine("  ALPHA <= -1.0");
                return 1;
            default:
                o = n + 1;

                return o;
        }
    }

    public static void cn_geg_03_xiu(int n, double alpha, int o, ref double[] x, ref double[] w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CN_GEG_03_XIU implements the Xiu precision 3 rule for region CN_GEG.
        //
        //  Discussion:
        //
        //    The rule has order 
        //
        //      O = 2 * N.
        //
        //    The rule has precision P = 3.
        //
        //    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
        //
        //      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
        //
        //    with -1.0 < alpha.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dongbin Xiu,
        //    Numerical integration formulas of degree two,
        //    Applied Numerical Mathematics,
        //    Volume 58, 2008, pages 1515-1520.
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, double ALPHA, the parameter.
        //    -1.0 < ALPHA.
        //
        //    Input, int O, the order.
        //
        //    Output, double X[N*O], the abscissas.
        //
        //    Output, double W[O], the weights.
        //
    {
        double arg;
        int expon;
        int i;
        int j;
            
        int r;
        double volume;

        switch (alpha)
        {
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("CN_GEG_03_XIU - Fatal error!");
                Console.WriteLine("  ALPHA <= -1.0");
                return;
        }

        expon = 0;
        volume = C1.c1_geg_monomial_integral(alpha, expon);
        volume = Math.Pow(volume, n);

        for (j = 0; j < o; j++)
        {
            i = 0;
            for (r = 1; r <= n / 2; r++)
            {
                arg = (2 * r - 1) * j * Math.PI / n;

                x[i + j * n] = Math.Sqrt(2.0) * Math.Cos(arg) / Math.Sqrt(2.0 * alpha + 3.0);
                i += 1;
                x[i + j * n] = Math.Sqrt(2.0) * Math.Sin(arg) / Math.Sqrt(2.0 * alpha + 3.0);
                i += 1;
            }

            if (i < n)
            {
                x[i + j * n] = Math.Sqrt(2.0) * typeMethods.r8_mop(j) / Math.Sqrt(2.0 * alpha + 3.0);
                switch (n)
                {
                    case 1:
                        x[i + j * n] /= Math.Sqrt(2.0);
                        break;
                }

                i += 1;
            }
        }

        for (j = 0; j < o; j++)
        {
            w[j] = volume / o;
        }
    }

    public static int cn_geg_03_xiu_size(int n, double alpha)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CN_GEG_03_XIU_SIZE sizes the Xiu precision 3 rule for region CN_GEG.
        //
        //  Discussion:
        //
        //    The rule has order 
        //
        //      O = 2 * N.
        //
        //    The rule has precision P = 3.
        //
        //    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
        //
        //      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
        //
        //    with -1.0 < alpha.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dongbin Xiu,
        //    Numerical integration formulas of degree two,
        //    Applied Numerical Mathematics,
        //    Volume 58, 2008, pages 1515-1520.
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, double ALPHA, the parameter.
        //    -1.0 < ALPHA.
        //
        //    Output, int CN_GEG_03_XIU_SIZE, the order.
        //
    {
        int o;

        switch (alpha)
        {
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("CN_GEG_03_XIU_SIZE - Fatal error!");
                Console.WriteLine("  ALPHA <= -1.0");
                return 1;
            default:
                o = 2 * n;

                return o;
        }
    }
        
    public static void en_her_02_xiu ( int n, int o, ref double[] x, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EN_HER_02_XIU implements the Xiu precision 2 rule for region EN_HER.
        //
        //  Discussion:
        //
        //    The rule has order 
        //
        //      O = N + 1.
        //
        //    The rule has precision P = 2.
        //
        //    EN_HER is the entire N-dimensional space with weight function
        //
        //      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 March 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dongbin Xiu,
        //    Numerical integration formulas of degree two,
        //    Applied Numerical Mathematics,
        //    Volume 58, 2008, pages 1515-1520.
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, int O, the order.
        //
        //    Output, double X[N*O], the abscissas.
        //
        //    Output, double W[O], the weights.
        //
    {
        double arg;
        double c1;
        double delta0;
        double gamma0;
        int i;
        int j;
            
        int r;
        double volume;
        double volume_1d;

        for (j = 0; j < o; j++)
        {
            i = 0;
            for (r = 1; r <= n / 2; r++)
            {
                arg = 2 * r * j * Math.PI / (n + 1);

                x[i + j * n] = Math.Sqrt(2.0) * Math.Cos(arg);
                i += 1;
                x[i + j * n] = Math.Sqrt(2.0) * Math.Sin(arg);
                i += 1;
            }

            if (i < n)
            {
                x[i + j * n] = typeMethods.r8_mop(j);
                i += 1;
            }
        }

        gamma0 = 2.0;
        delta0 = 0.0;
        c1 = 1.0;

        for (j = 0; j < o; j++)
        {
            for (i = 0; i < n; i++)
            {
                x[i + j * n] = (Math.Sqrt(gamma0 * c1) * x[i + j * n] - delta0) / gamma0;
            }
        }

        volume_1d = Math.Sqrt(Math.PI);
        volume = Math.Pow(volume_1d, n);

        for (j = 0; j < o; j++)
        {
            w[j] = volume / o;
        }
    }

    public static int en_her_02_xiu_size(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EN_HER_02_XIU_SIZE sizes the Xiu precision 2 rule for region EN_HER.
        //
        //  Discussion:
        //
        //    The rule has order O = N + 1;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dongbin Xiu,
        //    Numerical integration formulas of degree two,
        //    Applied Numerical Mathematics,
        //    Volume 58, 2008, pages 1515-1520.
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Output, int EN_HER_01_1_SIZE, the order.
        //
    {
        int o;

        o = n + 1;

        return o;
    }

    public static void en_her_03_xiu ( int n, int o, ref double[] x, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EN_HER_03_XIU implements the Xiu precision 3 rule for region EN_HER.
        //
        //  Discussion:
        //
        //    The rule has order 
        //
        //      O = 2 * N.
        //
        //    The rule has precision P = 3.
        //
        //    EN_HER is the entire N-dimensional space with weight function
        //
        //      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dongbin Xiu,
        //    Numerical integration formulas of degree two,
        //    Applied Numerical Mathematics,
        //    Volume 58, 2008, pages 1515-1520.
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, int O, the order.
        //
        //    Output, double X[N*O], the abscissas.
        //
        //    Output, double W[O], the weights.
        //
    {
        double arg;
        int i;
        int j;
            
        int r;
        double volume;
        double volume_1d;

        for (j = 0; j < o; j++)
        {
            i = 0;
            for (r = 1; r <= n / 2; r++)
            {
                arg = (2 * r - 1) * (j + 1) * Math.PI / n;
                x[i + j * n] = Math.Cos(arg);
                i += 1;
                x[i + j * n] = Math.Sin(arg);
                i += 1;
            }

            if (i < n)
            {
                x[i + j * n] = typeMethods.r8_mop(j + 1);
                switch (n)
                {
                    case 1:
                        x[i + j * n] /= Math.Sqrt(2.0);
                        break;
                }

                i += 1;
            }
        }

        volume_1d = Math.Sqrt(Math.PI);
        volume = Math.Pow(volume_1d, n);

        for (j = 0; j < o; j++)
        {
            w[j] = volume / o;
        }
    }

    public static int en_her_03_xiu_size(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EN_HER_03_XIU_SIZE sizes the Xiu precision 3 rule for region EN_HER.
        //
        //  Discussion:
        //
        //    The rule has order 
        //
        //      O = 2 * N.
        //
        //    The rule has precision P = 3.
        //
        //    EN_HER is the entire N-dimensional space with weight function
        //
        //      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dongbin Xiu,
        //    Numerical integration formulas of degree two,
        //    Applied Numerical Mathematics,
        //    Volume 58, 2008, pages 1515-1520.
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Output, int EN_HER_XIU_SIZE, the order.
        //
    {
        int o;

        o = 2 * n;

        return o;
    }


    public static void en_r2_02_xiu(int n, int o, ref double[] x, ref double[] w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EN_R2_02_XIU implements the Xiu rule for region EN_R2.
        //
        //  Discussion:
        //
        //    The rule has order 
        //
        //      O = N + 1.
        //
        //    The rule has precision P = 2.
        //
        //    EN_R2 is the entire N-dimensional space with weight function
        //
        //      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 March 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dongbin Xiu,
        //    Numerical integration formulas of degree two,
        //    Applied Numerical Mathematics,
        //    Volume 58, 2008, pages 1515-1520.
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, int O, the order.
        //
        //    Output, double X[N*O], the abscissas.
        //
        //    Output, double W[O], the weights.
        //
    {
        double arg;
        double c1;
        double delta0;
        double gamma0;
        int i;
        int j;
            
        int r;
        double volume;
        double volume_1d;

        for (j = 0; j < o; j++)
        {
            i = 0;
            for (r = 1; r <= n / 2; r++)
            {
                arg = 2 * r * j * Math.PI / (n + 1);

                x[i + j * n] = Math.Sqrt(2.0) * Math.Cos(arg);
                i += 1;
                x[i + j * n] = Math.Sqrt(2.0) * Math.Sin(arg);
                i += 1;
            }

            if (i < n)
            {
                x[i + j * n] = typeMethods.r8_mop(j);
                i += 1;
            }
        }

        gamma0 = 2.0;
        delta0 = 0.0;
        c1 = 1.0;

        for (j = 0; j < o; j++)
        {
            for (i = 0; i < n; i++)
            {
                x[i + j * n] = (Math.Sqrt(gamma0 * c1) * x[i + j * n] - delta0) / gamma0;
            }
        }

        volume_1d = Math.Sqrt(Math.PI);
        volume = Math.Pow(volume_1d, n);

        for (j = 0; j < o; j++)
        {
            w[j] = volume / o;
        }
    }

    public static int en_r2_02_xiu_size(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EN_R2_02_XIU_SIZE sizes the Xiu for region EN_R2.
        //
        //  Discussion:
        //
        //    The rule has order O = N + 1;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dongbin Xiu,
        //    Numerical integration formulas of degree two,
        //    Applied Numerical Mathematics,
        //    Volume 58, 2008, pages 1515-1520.
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Output, int EN_R2_01_1_SIZE, the order.
        //
    {
        int o;

        o = n + 1;

        return o;
    }

    public static void en_r2_03_xiu(int n, int o, ref double[] x, ref double[] w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EN_R2_03_XIU implements the Xiu precision 3 rule for region EN_R2.
        //
        //  Discussion:
        //
        //    The rule has order 
        //
        //      O = 2 * N.
        //
        //    The rule has precision P = 3.
        //
        //    EN_R2 is the entire N-dimensional space with weight function
        //
        //      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dongbin Xiu,
        //    Numerical integration formulas of degree two,
        //    Applied Numerical Mathematics,
        //    Volume 58, 2008, pages 1515-1520.
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, int O, the order.
        //
        //    Output, double X[N*O], the abscissas.
        //
        //    Output, double W[O], the weights.
        //
    {
        double arg;
        int i;
        int j;
            
        int r;
        double volume;

        volume = Math.Sqrt(Math.Pow(Math.PI, n));

        for (j = 0; j < o; j++)
        {
            i = 0;
            for (r = 1; r <= n / 2; r++)
            {
                arg = (2 * r - 1) * j * Math.PI / n;
                x[i + j * n] = Math.Cos(arg);
                i += 1;
                x[i + j * n] = Math.Sin(arg);
                i += 1;
            }

            if (i < n)
            {
                x[i + j * n] = typeMethods.r8_mop(j);
                switch (n)
                {
                    case 1:
                        x[i + j * n] /= Math.Sqrt(2.0);
                        break;
                }

                i += 1;
            }
        }

        for (j = 0; j < o; j++)
        {
            w[j] = volume / o;
        }
    }

    public static int en_r2_03_xiu_size(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EN_R2_03_XIU_SIZE sizes the Xiu precision 3 rule for region EN_R2.
        //
        //  Discussion:
        //
        //    The rule has order 
        //
        //      O = 2 * N.
        //
        //    The rule has precision P = 3.
        //
        //    EN_R2 is the entire N-dimensional space with weight function
        //
        //      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dongbin Xiu,
        //    Numerical integration formulas of degree two,
        //    Applied Numerical Mathematics,
        //    Volume 58, 2008, pages 1515-1520.
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Output, int EN_R2_XIU_SIZE, the order.
        //
    {
        int o;

        o = 2 * n;

        return o;
    }
}