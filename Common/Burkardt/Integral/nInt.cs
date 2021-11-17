using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.IntegralNS;

public static class nInt
{
    public static double box_nd(Func<int, double[], double> func, int dim_num,
            int order, double[] xtab, double[] weight, ref int eval_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BOX_ND estimates a multidimensional integral using a product rule.
        //
        //  Discussion:
        //
        //    The routine creates a DIM_NUM-dimensional product rule from a 1D rule
        //    supplied by the user.  The routine is fairly inflexible.  If
        //    you supply a rule for integration from -1 to 1, then your product
        //    box must be a product of DIM_NUM copies of the interval [-1,1].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Philip Davis, Philip Rabinowitz,
        //    Methods of Numerical Integration,
        //    Second Edition,
        //    Dover, 2007,
        //    ISBN: 0486453391,
        //    LC: QA299.3.D28.
        //
        //  Parameters:
        //
        //    Input, double FUNC ( int dim_num, double x[] ), evaluates
        //    the function to be integrated.
        //      
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int ORDER, the number of points used in the 1D rule.
        //
        //    Input, double XTAB[ORDER], the abscissas of the 1D rule.
        //
        //    Input, double WEIGHT[ORDER], the weights of the 1D rule.
        //
        //    Output, int *EVAL_NUM, the number of function evaluations.
        //
        //    Output, double BOX_ND, the approximate value of the integral.
        //
    {
        int dim;
        int[] indx;
        int k;
        double result;
        double w;
        double[] x;

        eval_num = 0;

        switch (dim_num)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("BOX_ND - Fatal error!");
                Console.WriteLine("  DIM_NUM < 1.");
                Console.WriteLine("  DIM_NUM = " + dim_num + "");
                return 1;
        }

        switch (order)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("BOX_ND - Fatal error!");
                Console.WriteLine("  ORDER < 1.");
                Console.WriteLine("  ORDER = " + order + "");
                return 1;
        }

        k = 0;
        result = 0.0;

        indx = new int[dim_num];
        x = new double[dim_num];

        for (;;)
        {
            BTuple.tuple_next(1, order, dim_num, ref k, ref indx);

            if (k == 0)
            {
                break;
            }

            w = 1.0;
            for (dim = 0; dim < dim_num; dim++)
            {
                w *= weight[indx[dim] - 1];
            }

            for (dim = 0; dim < dim_num; dim++)
            {
                x[dim] = xtab[indx[dim] - 1];
            }

            result += w * func(dim_num, x);
            eval_num += 1;
        }

        return result;
    }

    public static double monte_carlo_nd(Func<int, double[], double> func, int dim_num,
            double[] a, double[] b, int eval_num, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MONTE_CARLO_ND estimates a multidimensional integral using Monte Carlo.
        //
        //  Discussion:
        //
        //    Unlike the other routines, this routine requires the user to specify
        //    the number of function evaluations as an INPUT quantity.  
        //    No attempt at error estimation is made.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Philip Davis, Philip Rabinowitz,
        //    Methods of Numerical Integration,
        //    Second Edition,
        //    Dover, 2007,
        //    ISBN: 0486453391,
        //    LC: QA299.3.D28.
        //
        //  Parameters:
        //
        //    Input, double FUNC ( int dim_num, double x[] ), evaluates
        //    the function to be integrated.
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, double A[DIM_NUM], B[DIM_NUM], the integration limits.
        //
        //    Input, int EVAL_NUM, the number of function evaluations.
        //
        //    Input/output, int *SEED, a seed for the random number generator.
        //
        //    Output, double MONTE_CARLO_ND, the approximate value of the integral.
        //
    {
        int dim;
        int i;
        double result;
        double volume;
        double[] x;

        result = 0.0;

        for (i = 0; i < eval_num; i++)
        {
            x = UniformRNG.r8vec_uniform_01_new(dim_num, ref seed);

            result += func(dim_num, x);
        }

        volume = 1.0;
        for (dim = 0; dim < dim_num; dim++)
        {
            volume *= (b[dim] - a[dim]);
        }

        result = result * volume / eval_num;

        return result;
    }

    public static double p5_nd(Func<int, double[], double> func, int dim_num,
            double[] a, double[] b, ref int eval_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P5_ND estimates a multidimensional integral with a formula of exactness 5.
        //
        //  Discussion:
        //
        //    The routine uses a method which is exact for polynomials of total 
        //    degree 5 or less.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 February 2007
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
        //    C++ version by John Burkardt
        //
        //  Reference:
        //
        //    Philip Davis, Philip Rabinowitz,
        //    Methods of Numerical Integration,
        //    Second Edition,
        //    Dover, 2007,
        //    ISBN: 0486453391,
        //    LC: QA299.3.D28.
        //
        //  Parameters:
        //
        //    Input, double FUNC ( int dim_num, double x[] ), evaluates
        //    the function to be integrated.
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, double A[DIM_NUM], B[DIM_NUM], the integration limits.
        //
        //    Output, int *EVAL_NUM, the number of function evaluations.
        //
        //    Output, double P5_ND, the approximate value of the integral.
        //
    {
        double a0;
        double a1;
        double a2;
        double a3;
        double a4;
        double a5;
        int dim;
        double en;
        int i;
        int j;
        double result;
        double sum1;
        double sum2;
        double sum3;
        double volume;
        double[] work;

        eval_num = 0;

        switch (dim_num)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("P5_ND - Fatal error!");
                Console.WriteLine("  DIM_NUM < 1, DIM_NUM = " + dim_num + "");
                return 1;
        }

        a2 = 25.0 / 324.0;
        a3 = Math.Sqrt(0.6);
        en = dim_num;
        a0 = (25.0 * en * en - 115.0 * en + 162.0) / 162.0;
        a1 = (70.0 - 25.0 * en) / 162.0;

        volume = 1.0;
        for (dim = 0; dim < dim_num; dim++)
        {
            volume *= (b[dim] - a[dim]);
        }

        work = new double[dim_num];
        for (dim = 0; dim < dim_num; dim++)
        {
            work[dim] = 0.5 * (a[dim] + b[dim]);
        }

        result = 0.0;
        switch (volume)
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("P5_ND - Warning!");
                Console.WriteLine("  Volume = 0, integral = 0.");
                return result;
        }

        sum1 = a0 * func(dim_num, work);
        eval_num += 1;

        sum2 = 0.0;
        sum3 = 0.0;

        for (i = 0; i < dim_num; i++)
        {
            work[i] = 0.5 * (a[i] + b[i] + a3 * (b[i] - a[i]));
            sum2 += func(dim_num, work);
            eval_num += 1;

            work[i] = 0.5 * (a[i] + b[i] - a3 * (b[i] - a[i]));
            sum2 += func(dim_num, work);
            eval_num += 1;

            work[i] = 0.5 * (a[i] + b[i]);
        }

        switch (dim_num)
        {
            case > 1:
            {
                a4 = a3;

                for (;;)
                {
                    for (i = 0; i < dim_num - 1; i++)
                    {
                        work[i] = 0.5 * (a[i] + b[i] + a4 * (b[i] - a[i]));
                        a5 = a3;

                        for (;;)
                        {
                            for (j = i + 1; j < dim_num; j++)
                            {
                                work[j] = 0.5 * (a[j] + b[j] + a5 * (b[j] - a[j]));
                                sum3 += func(dim_num, work);
                                eval_num += 1;
                                work[j] = 0.5 * (a[j] + b[j]);
                            }

                            a5 = -a5;

                            if (0.0 <= a5)
                            {
                                break;
                            }
                        }

                        work[i] = 0.5 * (a[i] + b[i]);
                    }

                    a4 = -a4;

                    if (0.0 <= a4)
                    {
                        break;
                    }
                }

                break;
            }
        }

        result = volume * (sum1 + a1 * sum2 + a2 * sum3);

        return result;
    }

    public static double romberg_nd(Func<int, double[], double> func, double[] a,
            double[] b, int dim_num, int[] sub_num, int it_max, double tol, ref int ind,
            ref int eval_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ROMBERG_ND estimates a multidimensional integral using Romberg integration.
        //
        //  Discussion:
        //
        //    The routine uses a Romberg method based on the midpoint rule.
        //
        //    In the reference, this routine is called "NDIMRI".
        //
        //    Thanks to Barak Bringoltz for pointing out problems in a previous
        //    FORTRAN90 implementation of this routine.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 February 2007
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
        //    C++ version by John Burkardt
        //
        //  Reference:
        //
        //    Philip Davis, Philip Rabinowitz,
        //    Methods of Numerical Integration,
        //    Second Edition,
        //    Dover, 2007,
        //    ISBN: 0486453391,
        //    LC: QA299.3.D28.
        //
        //  Parameters:
        //
        //    Input, double FUNC ( int dim_num, double x[] ), evaluates
        //    the function to be integrated.
        //
        //    Input, double A[DIM_NUM], B[DIM_NUM], the integration limits.
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int SUB_NUM[DIM_NUM], the number of subintervals into
        //    which the I-th integration interval (A(I), B(I)) is
        //    initially subdivided.  SUB_NUM(I) must be greater than 0.
        //
        //    Input, int IT_MAX, the maximum number of iterations to
        //    be performed.  The number of function evaluations on
        //    iteration J is at least J**DIM_NUM, which grows very rapidly.
        //    IT_MAX should be small!
        //
        //    Input, double TOL, an error tolerance for the approximation
        //    of the integral.
        //
        //    Output, int *IND, error return flag.
        //    IND = -1 if the error tolerance could not be achieved.
        //    IND = 1 if the error tolerance was achieved.
        //
        //    Output, int *EVAL_NUM, the number of function evaluations.
        //
        //    Output, double ROMBERG_ND, the approximate value of the integral.
        //
        //  Local Parameters:
        //
        //    Local, int IWORK[DIM_NUM], a pointer used to generate all the
        //    points X in the product region.
        //
        //    Local, int IWORK2[IT_MAX], a counter of the number of points
        //    used at each step of the Romberg iteration.
        //
        //    Local, int SUB_NUM2[DIM_NUM], the number of subintervals used
        //    in each direction, a refinement of the user's input SUB_NUM.
        //
        //    Local, double TABLE[IT_MAX], the difference table.
        //
        //    Local, double X[DIM_NUM], an evaluation point.
        //
    {
        int dim;
        double factor;
        int i;
        int it;
        int[] iwork;
        int[] iwork2;
        int kdim;
        int ll;
        double result;
        double result_old = 0;
        double rnderr;
        int[] sub_num2;
        double sum1;
        double weight;
        double[] table;
        double[] x;

        eval_num = 0;

        switch (dim_num)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("ROMBERG_ND - Fatal error!");
                Console.WriteLine("  DIM_NUM is less than 1.  DIM_NUM = " + dim_num + "");
                return 1;
        }

        switch (it_max)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("ROMBERG_ND - Fatal error!");
                Console.WriteLine("  IT_MAX is less than 1.  IT_MAX = " + it_max + "");
                return 1;
        }

        for (i = 0; i < dim_num; i++)
        {
            switch (sub_num[i])
            {
                case <= 0:
                    Console.WriteLine("");
                    Console.WriteLine("ROMBERG_ND - Fatal error!");
                    Console.WriteLine("  SUB_NUM(I) is less than 1.");
                    Console.WriteLine("  for I = " + i + "");
                    Console.WriteLine("  SUB_NUM(I) = " + sub_num[i] + "");
                    return 1;
            }
        }

        iwork = new int[dim_num];
        iwork2 = new int[it_max];
        sub_num2 = new int[dim_num];
        table = new double[it_max];
        x = new double[dim_num];

        ind = 0;
        rnderr = typeMethods.r8_epsilon();
        iwork2[0] = 1;

        for (dim = 0; dim < dim_num; dim++)
        {
            sub_num2[dim] = sub_num[dim];
        }

        iwork2[1] = it_max switch
        {
            > 1 => 2,
            _ => iwork2[1]
        };

        it = 1;

        for (;;)
        {
            sum1 = 0.0;

            weight = 1.0;
            for (dim = 0; dim < dim_num; dim++)
            {
                weight = weight * (b[dim] - a[dim]) / sub_num2[dim];
            }

            //
            //  Generate every point X in the product region, and evaluate F(X).
            //
            for (dim = 0; dim < dim_num; dim++)
            {
                iwork[dim] = 1;
            }

            for (;;)
            {
                for (dim = 0; dim < dim_num; dim++)
                {
                    x[dim] =
                        ((2 * sub_num2[dim] - 2 * iwork[dim] + 1) * a[dim]
                         + (+2 * iwork[dim] - 1) * b[dim])
                        / (2 * sub_num2[dim]);
                }

                sum1 += func(dim_num, x);
                eval_num += 1;

                kdim = dim_num;

                while (0 < kdim)
                {
                    if (iwork[kdim - 1] < sub_num2[kdim - 1])
                    {
                        iwork[kdim - 1] += 1;
                        break;
                    }

                    iwork[kdim - 1] = 1;
                    kdim -= 1;
                }

                if (kdim == 0)
                {
                    break;
                }
            }

            //
            //  Done with summing.
            //
            table[it - 1] = weight * sum1;

            if (it <= 1)
            {
                result = table[0];
                result_old = result;

                if (it_max <= it)
                {
                    ind = 1;
                    break;
                }

                it += 1;
                for (dim = 0; dim < dim_num; dim++)
                {
                    sub_num2[dim] = iwork2[it - 1] * sub_num2[dim];
                }

                continue;
            }

            //
            //  Compute the difference table for Richardson extrapolation.
            // 
            for (ll = 2; ll <= it; ll++)
            {
                i = it + 1 - ll;
                factor = iwork2[i - 1] * iwork2[i - 1]
                         / (double) (iwork2[it - 1] * iwork2[it - 1] - iwork2[i - 1] * iwork2[i - 1]);
                table[i] += (table[i] - table[i - 1]) * factor;
            }

            result = table[0];
            //
            //  Terminate successfully if the estimated error is acceptable.
            //
            if (Math.Abs(result - result_old) <= Math.Abs(result * (tol + rnderr)))
            {
                ind = 1;
                break;
            }

            //
            //  Terminate unsuccessfully if the iteration limit has been reached.
            //
            if (it_max <= it)
            {
                ind = -1;
                break;
            }

            //
            //  Prepare for another step.
            //
            result_old = result;

            it += 1;

            iwork2[it - 1] = (int) (1.5 * iwork2[it - 2]);

            for (dim = 0; dim < dim_num; dim++)
            {
                sub_num2[dim] = (int) (1.5 * sub_num2[dim]);
            }
        }

        return result;
    }

    public static void sample_nd(Func<int, double[], double> func, int k1, int k2,
            int dim_num, double[] est1, double[] err1, double[] dev1, double[] est2,
            double[] err2, double[] dev2, ref int eval_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SAMPLE_ND estimates a multidimensional integral using sampling.
        //
        //  Discussion:
        //
        //    This routine computes two sequences of integral estimates, EST1 
        //    and EST2, for indices K going from K1 to K2.  These estimates are 
        //    produced by the generation of 'random' abscissas in the region.  
        //    The process can become very expensive if high accuracy is needed.
        //
        //    The total number of function evaluations is
        //    4*(K1^DIM_NUM+(K1+1)^DIM_NUM+...+(K2-1)^DIM_NUM+K2^DIM_NUM), and K2
        //    should be chosen so as to make this quantity reasonable.
        //    In most situations, EST2(K) are much better estimates than EST1(K).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 March 2007
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
        //    C++ version by John Burkardt
        //
        //  Reference:
        //
        //    Philip Davis, Philip Rabinowitz,
        //    Methods of Numerical Integration,
        //    Second Edition,
        //    Dover, 2007,
        //    ISBN: 0486453391,
        //    LC: QA299.3.D28.
        //
        //  Parameters:
        //
        //    Input, double FUNC ( int dim_num, double x[] ), evaluates
        //    the function to be integrated.
        //
        //    Input, int K1, the beginning index for the iteration.
        //    1 <= K1 <= K2.
        //
        //    Input, int K2, the final index for the iteration.  K1 <= K2.
        //    Increasing K2 increases the accuracy of the calculation,
        //    but vastly increases the work and running time of the code.
        //
        //    Input, int DIM_NUM, the spatial dimension.  1 <= DIM_NUM <= 10.
        //
        //    Output, double EST1[K2].  Entries K1 through K2 contain
        //    successively better estimates of the integral.
        //
        //    Output, double ERR1[K2].  Entries K1 through K2 contain
        //    the corresponding estimates of the integration errors.
        //
        //    Output, double DEV1[K2].  Entries K1 through K2 contain
        //    estimates of the reliability of the the integration.
        //    If consecutive values DEV1(K) and DEV1(K+1) do not differ
        //    by more than 10 percent, then ERR1(K) can be taken as
        //    a reliable upper bound on the difference between EST1(K)
        //    and the true value of the integral.
        //
        //    Output, double EST2[K2].  Entries K2 through K2 contain
        //    successively better estimates of the integral.
        //
        //    Output, double ERR2[K2].  Entries K2 through K2 contain
        //    the corresponding estimates of the integration errors.
        //
        //    Output, double DEV2[K2].  Entries K2 through K2 contain
        //    estimates of the reliability of the the integration.
        //    If consecutive values DEV2(K) and DEV2(K+2) do not differ
        //    by more than 10 percent, then ERR2(K) can be taken as
        //    a reliable upper bound on the difference between EST2(K)
        //    and the true value of the integral.
        //
        //    Output, int *EVAL_NUM, the number of function evaluations.
        //
    {
        int DIM_MAX = 10;

        double ak;
        double ak1;
        double akn;
        double[] al =
        {
            0.4142135623730950,
            0.7320508075688773,
            0.2360679774997897,
            0.6457513110645906,
            0.3166247903553998,
            0.6055512754639893,
            0.1231056256176605,
            0.3589989435406736,
            0.7958315233127195,
            0.3851648071345040
        };
        double b;
        double[] be;
        double bk;
        double d1;
        double d2;
        double[] dex;
        int dim;
        double g;
        double[] ga;
        int i;
        int j;
        int k;
        int key;
        bool more;
        double[] p1;
        double[] p2;
        double[] p3;
        double[] p4;
        double s1;
        double s2;
        double t;
        double y1;
        double y2;
        double y3;
        double y4;

        eval_num = 0;
        switch (dim_num)
        {
            //
            //  Check input.
            //
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("SAMPLE_ND - Fatal error!");
                Console.WriteLine("  DIM_NUM must be at least 1,");
                Console.WriteLine("  but DIM_NUM = " + dim_num + "");
                return;
        }

        if (DIM_MAX < dim_num)
        {
            Console.WriteLine("");
            Console.WriteLine("SAMPLE_ND - Fatal error!");
            Console.WriteLine("  DIM_NUM must be no more than DIM_MAX = " + DIM_MAX + "");
            Console.WriteLine("  but DIM_NUM = " + dim_num + "");
            return;
        }

        switch (k1)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("SAMPLE_ND - Fatal error!");
                Console.WriteLine("  K1 must be at least 1, but K1 = " + k1 + "");
                return;
        }

        if (k2 < k1)
        {
            Console.WriteLine("");
            Console.WriteLine("SAMPLE_ND - Fatal error!");
            Console.WriteLine("  K1 may not be greater than K2, but ");
            Console.WriteLine("  K1 = " + k1 + "");
            Console.WriteLine("  K2 = " + k2 + "");
            return;
        }

        be = new double[dim_num];
        dex = new double[dim_num];
        ga = new double[dim_num];
        p1 = new double[dim_num];
        p2 = new double[dim_num];
        p3 = new double[dim_num];
        p4 = new double[dim_num];

        for (dim = 0; dim < dim_num; dim++)
        {
            be[dim] = al[dim];
        }

        for (dim = 0; dim < dim_num; dim++)
        {
            ga[dim] = al[dim];
        }

        for (dim = 0; dim < dim_num; dim++)
        {
            dex[dim] = 0.0;
        }

        for (k = k1; k <= k2; k++)
        {
            ak = k;
            key = 0;
            ak1 = ak - 1.1;
            s1 = 0.0;
            d1 = 0.0;
            s2 = 0.0;
            d2 = 0.0;
            akn = Math.Pow(ak, dim_num);
            t = Math.Sqrt(Math.Pow(ak, dim_num)) * ak;
            bk = 1.0 / ak;

            for (;;)
            {
                key += 1;

                if (key != 1)
                {
                    key -= 1;
                    more = false;
                    for (j = 0; j < dim_num; j++)
                    {
                        if (dex[j] <= ak1)
                        {
                            dex[j] += 1.0;
                            more = true;
                            break;
                        }

                        dex[j] = 0.0;
                    }

                    if (!more)
                    {
                        break;
                    }
                }

                for (i = 0; i < dim_num; i++)
                {
                    b = be[i] + al[i];
                    switch (b)
                    {
                        case > 1.0:
                            b -= 1.0;
                            break;
                    }

                    g = ga[i] + b;
                    switch (g)
                    {
                        case > 1.0:
                            g -= 1.0;
                            break;
                    }

                    be[i] = b + al[i];
                    switch (be[i])
                    {
                        case > 1.0:
                            be[i] -= 1.0;
                            break;
                    }

                    ga[i] = be[i] + g;
                    switch (ga[i])
                    {
                        case > 1.0:
                            ga[i] -= 1.0;
                            break;
                    }

                    p1[i] = (dex[i] + g) * bk;
                    p2[i] = (dex[i] + 1.0 - g) * bk;
                    p3[i] = (dex[i] + ga[i]) * bk;
                    p4[i] = (dex[i] + 1.0 - ga[i]) * bk;
                }

                y1 = func(dim_num, p1);
                eval_num += 1;
                //
                //  There may be an error in the next two lines,
                //  but oddly enough, that is how the original reads
                //
                y3 = func(dim_num, p2);
                eval_num += 1;
                y2 = func(dim_num, p3);
                eval_num += 1;
                y4 = func(dim_num, p4);
                eval_num += 1;

                s1 = s1 + y1 + y2;
                d1 += (y1 - y2) * (y1 - y2);
                s2 = s2 + y3 + y4;
                d2 += (y1 + y3 - y2 - y4) * (y1 + y3 - y2 - y4);
            }

            est1[k - 1] = 0.5 * s1 / akn;
            err1[k - 1] = 1.5 * Math.Sqrt(d1) / akn;
            dev1[k - 1] = err1[k - 1] * t;
            est2[k - 1] = 0.25 * (s1 + s2) / akn;
            err2[k - 1] = 0.75 * Math.Sqrt(d2) / akn;
            dev2[k - 1] = err2[k - 1] * t * ak;
        }
    }

    public static double sum2_nd(Func<int, double[], double> func, double[] xtab,
            double[] weight, int[] order, int dim_num, ref int eval_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SUM2_ND estimates a multidimensional integral using a product rule.
        //
        //  Discussion:
        //
        //    The routine uses a product rule supplied by the user.
        //
        //    The region may be a product of any combination of finite,
        //    semi-infinite, or infinite intervals.
        //
        //    For each factor in the region, it is assumed that an integration
        //    rule is given, and hence, the region is defined implicitly by
        //    the integration rule chosen.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 February 2007
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
        //    C++ version by John Burkardt
        //
        //  Reference:
        //
        //    Philip Davis, Philip Rabinowitz,
        //    Methods of Numerical Integration,
        //    Second Edition,
        //    Dover, 2007,
        //    ISBN: 0486453391,
        //    LC: QA299.3.D28.
        //
        //  Parameters:
        //
        //    Input, double FUNC ( int dim_num, double x[] ), evaluates
        //    the function to be integrated.
        //
        //    Input, double XTAB[DIM_NUM*ORDER_MAX].  XTAB(I,J) is the 
        //    I-th abscissa of the J-th rule.
        //
        //    Input, double WEIGHT[DIM_NUM*ORDER_MAX].  WEIGHT(I,J) is the 
        //    I-th weight for the J-th rule.
        //
        //    Input, int ORDER[DIM_NUM].  ORDER(I) is the number of
        //    abscissas to be used in the J-th rule.  ORDER(I) must be
        //    greater than 0 and less than or equal to ORDER_MAX.
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Output, int EVAL_NUM, the number of function evaluations.
        //
        //    Output, double SUM2_ND, the approximate value of the integral.
        //
    {
        int dim;
        int i;
        int[] iwork;
        int k;
        int m1;
        double result;
        double w1;
        double[] work;
        //
        //  Default values.
        //
        result = 0.0;
        eval_num = 0;

        switch (dim_num)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("SUM2_ND - Fatal error!");
                Console.WriteLine("  DIM_NUM < 1");
                Console.WriteLine("  DIM_NUM = " + dim_num + "");
                return 1;
        }

        for (i = 0; i < dim_num; i++)
        {
            switch (order[i])
            {
                case < 1:
                    Console.WriteLine("");
                    Console.WriteLine("SUM2_ND - Fatal error!");
                    Console.WriteLine("  ORDER(I) < 1.");
                    Console.WriteLine("  For I = " + i + "");
                    Console.WriteLine("  ORDER(I) = " + order[i] + "");
                    return 1;
            }
        }

        iwork = new int[dim_num];
        work = new double[dim_num];

        for (dim = 0; dim < dim_num; dim++)
        {
            iwork[dim] = 1;
        }

        for (;;)
        {
            k = 1;

            w1 = 1.0;
            for (i = 0; i < dim_num; i++)
            {
                m1 = iwork[i];
                work[i] = xtab[i + (m1 - 1) * dim_num];
                w1 *= weight[i + (m1 - 1) * dim_num];
            }

            result += w1 * func(dim_num, work);
            eval_num += 1;

            while (iwork[k - 1] == order[k - 1])
            {
                iwork[k - 1] = 1;
                k += 1;

                if (dim_num < k)
                {
                    return result;
                }
            }

            iwork[k - 1] += 1;
        }

        return result;
    }
}