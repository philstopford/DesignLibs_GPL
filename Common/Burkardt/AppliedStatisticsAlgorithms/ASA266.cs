using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.AppliedStatistics;

public static partial class Algorithms
{

    public static double alogam(double x, ref int ifault)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ALOGAM computes the logarithm of the Gamma function.
        //
        //  Modified:
        //
        //    22 January 2008
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Malcolm Pike, David Hill.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Malcolm Pike, David Hill,
        //    Algorithm 291:
        //    Logarithm of Gamma Function,
        //    Communications of the ACM,
        //    Volume 9, Number 9, September 1966, page 684.
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the Gamma function.
        //    X should be greater than 0.
        //
        //    Output, int *IFAULT, error flag.
        //    0, no error.
        //    1, X <= 0.
        //
        //    Output, double ALOGAM, the logarithm of the Gamma
        //    function of X.
        //
    {
        double f;
        double value;
        double z;

        switch (x)
        {
            case <= 0.0:
                ifault = 1;
                value = 0.0;
                return value;
        }

        ifault = 0;
        double y = x;

        switch (x)
        {
            case < 7.0:
            {
                f = 1.0;
                z = y;

                while (z < 7.0)
                {
                    f *= z;
                    z += 1.0;
                }

                y = z;
                f = -Math.Log(f);
                break;
            }
            default:
                f = 0.0;
                break;
        }

        z = 1.0 / y / y;

        value = f + (y - 0.5) * Math.Log(y) - y
                + 0.918938533204673 +
                (((
                      -0.000595238095238 * z
                      + 0.000793650793651) * z
                  - 0.002777777777778) * z
                 + 0.083333333333333) / y;

        return value;
    }

    public static double digamma(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIGAMMA calculates DIGAMMA ( X ) = d ( LOG ( GAMMA ( X ) ) ) / dX
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 June 2013
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Jose Bernardo.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Jose Bernardo,
        //    Algorithm AS 103:
        //    Psi ( Digamma ) Function,
        //    Applied Statistics,
        //    Volume 25, Number 3, 1976, pages 315-317.
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the digamma function.
        //    0 < X.
        //
        //    Output, double DIGAMMA, the value of the digamma function at X.
        //
    {
        double euler_mascheroni = 0.57721566490153286060;
        switch (x)
        {
            //
            //  Check the input.
            //
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("DIGAMMA - Fatal error!");
                Console.WriteLine("  X <= 0.0.");
                return 1;
        }

        //
        //  Initialize.
        //
        double x2 = x;
        double value = 0.0;
        switch (x2)
        {
            //
            //  Use approximation for small argument.
            //
            case <= 0.00001:
                value = -euler_mascheroni - 1.0 / x2;
                return value;
        }

        //
        //  Reduce to DIGAMMA(X + N).
        //
        while (x2 < 8.5)
        {
            value -= 1.0 / x2;
            x2 += 1.0;
        }

        //
        //  Use Stirling's (actually de Moivre's) expansion.
        //
        double r = 1.0 / x2;
        value = value + Math.Log(x2) - 0.5 * r;
        r *= r;
        value -= r * (1.0 / 12.0
                      - r * (1.0 / 120.0
                             - r * 1.0 / 252.0));

        return value;
    }

    public static void dirichlet_check(int n, double[] a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIRICHLET_CHECK checks the parameters of the Dirichlet PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of components.
        //
        //    Input, double A[N], the probabilities for each component.
        //    Each A(I) should be nonnegative, and at least one should be positive.
        //
    {
        int i;

        int positive = 0;

        for (i = 0; i < n; i++)
        {
            switch (a[i])
            {
                case < 0.0:
                    Console.WriteLine("");
                    Console.WriteLine("DIRICHLET_CHECK - Fatal error!");
                    Console.WriteLine("  A(I) < 0.");
                    Console.WriteLine("  For I = " + i + "");
                    Console.WriteLine("  A[I] = " + a[i] + "");
                    break;
                case > 0.0:
                    positive = 1;
                    break;
            }
        }

        switch (positive)
        {
            case 0:
                Console.WriteLine("");
                Console.WriteLine("DIRICHLET_CHECK - Fatal error!");
                Console.WriteLine("  All parameters are zero!");
                break;
        }
    }

    public static void dirichlet_estimate(int k, int n, double[] x, int ix, int init,
            ref double[] alpha, ref double rlogl, ref double[] v, ref double[] g, ref int niter,
            ref double s, ref double eps, ref int ifault)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIRICHLET_ESTIMATE estimates the parameters of a Dirichlet distribution.
        //
        //  Discussion:
        //
        //    This routine requires several auxilliary routines:
        //
        //      ALOGAM (CACM algorithm 291 or AS 245),
        //      DIGAMA (AS 103),
        //      GAMMAD (AS 239),
        //      PPCHI2 (AS 91),
        //      TRIGAM (AS 121).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2013
        //
        //  Author:
        //
        //    Original FORTRAN77 version by A Naryanan.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    A. Naryanan,
        //    Algorithm AS 266:
        //    Maximum Likelihood Estimation of the Parameters of the
        //    Dirichlet Distribution,
        //    Applied Statistics,
        //    Volume 40, Number 2, 1991, pages 365-374.
        //
        //  Parameters:
        //
        //    Input, int K, the number of parameters.
        //    2 <= K.
        //
        //    Input, int N, the number of observations.
        //    K < N.
        //
        //    Input, double X[IX*K], contains the N by K array of samples
        //    from the distribution.  X(I,J) is the J-th component of
        //    the I-th sample.
        //
        //    Input, int IX, the leading dimension of the array X.
        //    N <= IX.
        //
        //    Input, int INIT, specifies how the parameter estimates
        //    are to be initialized:
        //    1, use the method of moments;
        //    2, initialize each ALPHA to the minimum of X;
        //    otherwise, the input values of ALPHA already contain estimates.
        //
        //    Input/output, double ALPHA[K].
        //    On input, if INIT is not 1 or 2, then ALPHA must contain
        //    initial estimates for the parameters.
        //    On output, with IFAULT = 0, ALPHA contains the computed
        //    estimates for the parameters.
        //
        //    Output, double &RLOGL, the value of the log-likelihood function
        //    at the solution point.
        //
        //    Output, double V[K*K]; the covariance between ALPHA(I) and ALPHA(J).
        //
        //    Output, double G[K], contains an estimate of the derivative of
        //    the log-likelihood with respect to each component of ALPHA.
        //
        //    Output, int &NITER, contains the number of Newton-Raphson
        //    iterations performed.
        //
        //    Output, double &S, the value of the chi-squared statistic.
        //
        //    Output, double &EPS, contains the probability that the 
        //    chi-squared statistic is less than S.
        //
        //    Output, int &IFAULT, error indicator.
        //    0, no error, the results were computed successfully;
        //    1, K < 2;
        //    2, N <= K;
        //    3, IX < N;
        //    4, if X(I,J) <= 0 for any I or J, or if
        //       ABS ( Sum ( 1 <= J <= K ) X(I,J) - 1 ) >= GAMMA = 0.001;
        //    5, if IFAULT is returned nonzero from the chi-square
        //       routine PPCHI2;
        //    6, if ALPHA(J) <= 0 for any J during any step of the iteration;
        //    7, if MAXIT iterations were carried out but convergence
        //       was not achieved.
        //
    {
        const double alpha_min = 0.00001;
        const double gamma = 0.0001;
        int ifault2 = 0;
        const int it_max = 100;
        double sum1;
        double sum2;

        ifault = 0;
        switch (k)
        {
            //
            //  Check the input arguments.
            //
            case < 2:
                ifault = 1;
                Console.WriteLine("");
                Console.WriteLine("DIRICHLET_ESTIMATE - Fatal error!");
                Console.WriteLine("  K < 2.");
                break;
        }

        if (n <= k)
        {
            ifault = 2;
            Console.WriteLine("");
            Console.WriteLine("DIRICHLET_ESTIMATE - Fatal error!");
            Console.WriteLine("  N <= K.");
        }

        if (ix < n)
        {
            ifault = 3;
            Console.WriteLine("");
            Console.WriteLine("DIRICHLET_ESTIMATE - Fatal error!");
            Console.WriteLine("  IX < N.");
        }

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < k; j++)
            {
                switch (x[i + j * ix])
                {
                    case <= 0.0:
                        niter = i;
                        ifault = 4;
                        Console.WriteLine("");
                        Console.WriteLine("DIRICHLET_ESTIMATE - Fatal error!");
                        Console.WriteLine("  X(I,J) <= 0.");
                        break;
                }
            }

            sum2 = 0.0;
            for (int j = 0; j < k; j++)
            {
                sum2 += x[i + j * ix];
            }

            if (gamma <= Math.Abs(sum2 - 1.0))
            {
                ifault = 4;
                niter = i;
                Console.WriteLine("");
                Console.WriteLine("DIRICHLET_ESTIMATE - Fatal error!");
                Console.WriteLine("  Row I does not sum to 1.");
            }
        }

        ifault = 0;

        double an = n;
        double rk = k;
        niter = 0;
        switch (init)
        {
            //
            //  Calculate initial estimates using the method of moments.
            //
            case 1:
            {
                sum2 = 0.0;
                for (int j = 0; j < k - 1; j++)
                {
                    sum1 = 0.0;
                    for (int i = 0; i < n; i++)
                    {
                        sum1 += x[i + j * ix];
                    }

                    alpha[j] = sum1 / an;
                    sum2 += alpha[j];
                }

                alpha[k - 1] = 1.0 - sum2;

                double x12 = 0.0;
                for (int i = 0; i < n; i++)
                {
                    x12 += Math.Pow(x[i + 0 * ix], 2);
                }

                x12 /= an;
                double varp1 = x12 - Math.Pow(alpha[0], 2);

                double x11 = (alpha[0] - x12) / varp1;
                for (int j = 0; j < k; j++)
                {
                    alpha[j] = x11 * alpha[j];
                }

                break;
            }
            //
            //  Calculate initial estimates using Ronning's suggestion.
            //
            case 2:
            {
                double x_min = x[0 + 0 * ix];
                for (int j = 0; j < k; j++)
                {
                    for (int i = 0; i < n; i++)
                    {
                        x_min = Math.Min(x_min, x[i + j * ix]);
                    }
                }

                x_min = Math.Max(x_min, alpha_min);

                for (int j = 0; j < k; j++)
                {
                    alpha[j] = x_min;
                }

                break;
            }
        }

        //
        //  Check whether any ALPHA's are negative or zero.
        //
        for (int j = 0; j < k; j++)
        {
            switch (alpha[j])
            {
                case <= 0.0:
                    ifault = 6;
                    Console.WriteLine("");
                    Console.WriteLine("DIRICHLET_ESTIMATE - Fatal error!");
                    Console.WriteLine("  For J = " + j + "");
                    Console.WriteLine("  ALPHA(J) = " + alpha[j] + "");
                    Console.WriteLine("  but ALPHA(J) must be positive.");
                    break;
            }
        }

        //
        //  Calculate N * log ( G(J) ) for J = 1,2,...,K and store in WORK array.
        //
        double[] work = new double[k];

        for (int j = 0; j < k; j++)
        {
            work[j] = 0.0;
            for (int i = 0; i < n; i++)
            {
                work[j] += Math.Log(x[i + j * ix]);
            }
        }

        //
        //  Call Algorithm AS 91 to compute CHI2, the chi-squared value.
        //
        double gg = Helpers.LogGamma(rk / 2.0);

        double chi2 = ppchi2(gamma, rk, gg, ref ifault2);

        if (ifault2 != 0)
        {
            ifault = 5;
            Console.WriteLine("");
            Console.WriteLine("DIRICHLET_ESTIMATE - Fatal error!");
            Console.WriteLine("  PPCHI2 returns error code.");
        }

        //
        //  Carry out the Newton iteration.
        //
        double[] work2 = new double[k];

        for (int it_num = 1; it_num <= it_max; it_num++)
        {
            sum2 = typeMethods.r8vec_sum(k, alpha);

            sum1 = 0.0;
            for (int j = 0; j < k; j++)
            {
                work2[j] = trigamma(alpha[j], ref ifault2);
                sum1 += 1.0 / work2[j];
            }

            double beta = trigamma(sum2, ref ifault2);
            beta = an * beta / (1.0 - beta * sum1);

            double temp = digamma(sum2);

            for (int j = 0; j < k; j++)
            {
                g[j] = an * (temp - digamma(alpha[j])) + work[j];
            }

            //
            //  Calculate the lower triangle of the Variance-Covariance matrix V.
            //
            sum2 = beta / an / an;
            for (int i = 0; i < k; i++)
            {
                for (int j = 0; j < k; j++)
                {
                    v[i + j * k] = sum2 / (work2[i] * work2[j]);
                    if (i == j)
                    {
                        v[i + j * k] += 1.0 / (an * work2[j]);
                    }
                }
            }

            //
            //  Postmultiply the Variance-Covariance matrix V by G and store in WORK2.
            //
            work2 = typeMethods.r8mat_mv_new(k, k, v, g);
            //
            //  Update the ALPHA's.
            //
            niter = it_num;

            for (int j = 0; j < k; j++)
            {
                alpha[j] += work2[j];
                alpha[j] = Math.Max(alpha[j], alpha_min);
            }

            for (int j = 0; j < k; j++)
            {
                switch (alpha[j])
                {
                    case <= 0.0:
                        ifault = 6;
                        Console.WriteLine("");
                        Console.WriteLine("DIRICHLET_ESTIMATE - Fatal error!");
                        Console.WriteLine("  Newton iteration " + it_num + "");
                        Console.WriteLine("  Computed ALPHA[J] <= 0.");
                        Console.WriteLine("  J = " + j + "");
                        Console.WriteLine("  ALPHA[J] = " + alpha[j] + "");
                        break;
                }
            }

            //
            //  Test for convergence.
            //
            s = typeMethods.r8vec_dot_product(k, g, work2);

            if (s < chi2)
            {
                eps = gammad(s / 2.0, rk / 2.0, ref ifault2);

                sum2 = typeMethods.r8vec_sum(k, alpha);

                rlogl = 0.0;
                for (int j = 0; j < k; j++)
                {
                    rlogl = rlogl + (alpha[j] - 1.0) * work[j] - an *
                        Helpers.LogGamma(alpha[j]);
                }

                rlogl += an * Helpers.LogGamma(sum2);

                return;
            }
        }

        ifault = 7;

        Console.WriteLine("");
        Console.WriteLine("DIRICHLET_ESTIMATE - Fatal error!");
        Console.WriteLine("  No convergence.");
    }

    public static double[] dirichlet_mean(int n, double[] a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIRICHLET_MEAN returns the means of the Dirichlet PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of components.
        //
        //    Input, double A[N], the probabilities for each component.
        //    Each A[I] should be nonnegative, and at least one should be positive.
        //
        //    Output, double DIRICHLET_MEAN[N], the means of the PDF.
        //
    {
        dirichlet_check(n, a);

        double[] mean = new double[n];

        double a_sum = typeMethods.r8vec_sum(n, a);

        for (int i = 0; i < n; i++)
        {
            mean[i] = a[i] / a_sum;
        }

        return mean;
    }

    public static double[] dirichlet_mix_mean(int comp_max, int comp_num, int elem_num,
            double[] a, double[] comp_weight)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIRICHLET_MIX_MEAN returns the means of a Dirichlet mixture PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int COMP_MAX, the leading dimension of A, which
        //    must be at least COMP_NUM.
        //
        //    Input, int COMP_NUM, the number of components in the
        //    Dirichlet mixture density, that is, the number of distinct Dirichlet PDF's
        //    that are mixed together.
        //
        //    Input, int ELEM_NUM, the number of elements of an
        //    observation.
        //
        //    Input, double A[COMP_MAX*ELEM_NUM], the probabilities for
        //    element ELEM_NUM in component COMP_NUM.
        //    Each A(I,J) should be greater than or equal to 0.0.
        //
        //    Input, double COMP_WEIGHT[COMP_NUM], the mixture weights of
        //    the densities.
        //    These do not need to be normalized.  The weight of a given component is
        //    the relative probability that that component will be used to generate
        //    the sample.
        //
        //    Output, double DIRICHLET_MIX_MEAN[ELEM_NUM], the means for each element.
        //
    {
        //
        //  Check.
        //
        for (int comp_i = 0; comp_i < comp_num; comp_i++)
        {
            for (int elem_i = 0; elem_i < elem_num; elem_i++)
            {
                switch (a[comp_i + elem_i * comp_max])
                {
                    case < 0.0:
                        Console.WriteLine("");
                        Console.WriteLine("DIRICHLET_MIX_MEAN - Fatal error!");
                        Console.WriteLine("  A(COMP,ELEM) < 0.");
                        Console.WriteLine("  COMP = " + comp_i + "");
                        Console.WriteLine("  ELEM = " + elem_i + "");
                        Console.WriteLine("  A[COMP,ELEM] = " + a[comp_i + elem_i * comp_max] + "");
                        break;
                }
            }
        }

        double comp_weight_sum = typeMethods.r8vec_sum(comp_num, comp_weight);

        double[] a_sum = new double[comp_num];

        for (int comp_i = 0; comp_i < comp_num; comp_i++)
        {
            a_sum[comp_i] = 0.0;
            for (int elem_i = 0; elem_i < elem_num; elem_i++)
            {
                a_sum[comp_i] += a[comp_i + elem_i * comp_max];
            }
        }

        double[] mean = new double[elem_num];

        for (int elem_i = 0; elem_i < elem_num; elem_i++)
        {
            mean[elem_i] = 0.0;
            for (int comp_i = 0; comp_i < comp_num; comp_i++)
            {
                mean[elem_i] += comp_weight[comp_i] * a[comp_i + elem_i * comp_max] / a_sum[comp_i];
            }

            mean[elem_i] /= comp_weight_sum;
        }

        return mean;
    }

    public static double[] dirichlet_mix_sample(int comp_max, int comp_num, int elem_num,
            double[] a, double[] comp_weight, ref typeMethods.r8NormalData data, ref int seed, ref int comp)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIRICHLET_MIX_SAMPLE samples a Dirichlet mixture PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int COMP_MAX, the leading dimension of A, which 
        //    must be at least COMP_NUM.
        //
        //    Input, int COMP_NUM, the number of components in the
        //    Dirichlet mixture density, that is, the number of distinct Dirichlet PDF's
        //    that are mixed together.
        //
        //    Input, int ELEM_NUM, the number of elements of an
        //    observation.
        //
        //    Input, double A[COMP_MAX*ELEM_NUM], the probabilities for
        //    element ELEM_NUM in component COMP_NUM.
        //    Each A(I,J) should be greater than or equal to 0.0.
        //
        //    Input, double COMP_WEIGHT[COMP_NUM], the mixture weights of
        //    the densities.
        //    These do not need to be normalized.  The weight of a given component is
        //    the relative probability that that component will be used to generate
        //    the sample.
        //
        //    Input/output, int &SEED, a seed for the random 
        //    number generator.
        //
        //    Output, int &COMP, the index of the component of the
        //    Dirichlet mixture that was chosen to generate the sample.
        //
        //    Output, double DIRICHLET_MIX_SAMPLE[ELEM_NUM], a sample of the PDF.
        //
    {
        //
        //  Check.
        //
        for (int comp_i = 0; comp_i < comp_num; comp_i++)
        {
            for (int elem_i = 0; elem_i < elem_num; elem_i++)
            {
                switch (a[comp_i + elem_i * comp_max])
                {
                    case < 0.0:
                        Console.WriteLine("");
                        Console.WriteLine("DIRICHLET_MIX_SAMPLE - Fatal error!");
                        Console.WriteLine("  A(COMP,ELEM) < 0.'");
                        Console.WriteLine("  COMP = " + comp_i + "");
                        Console.WriteLine("  ELEM = " + elem_i + "");
                        Console.WriteLine("  A(COMP,ELEM) = " + a[comp_i + elem_i * comp_max] + "");
                        break;
                }
            }
        }

        //
        //  Choose a particular density MIX.
        //
        double comp_weight_sum = typeMethods.r8vec_sum(comp_num, comp_weight);

        double r = UniformRNG.r8_uniform_ab(0.0, comp_weight_sum, ref seed);

        comp = 0;
        double sum2 = 0.0;

        while (comp < comp_num)
        {
            sum2 += comp_weight[comp];

            if (r <= sum2)
            {
                break;
            }

            comp += 1;
        }

        //
        //  Sample density COMP.
        //
        double[] x = new double[elem_num];

        for (int elem_i = 0; elem_i < elem_num; elem_i++)
        {
            x[elem_i] = gamma_sample(a[comp + elem_i * comp_max], 1.0, ref data, ref seed);
        }

        //
        //  Normalize the result.
        //
        double x_sum = typeMethods.r8vec_sum(elem_num, x);

        for (int elem_i = 0; elem_i < elem_num; elem_i++)
        {
            x[elem_i] /= x_sum;
        }

        return x;
    }

    public static double[] dirichlet_sample(int n, double[] a, ref typeMethods.r8NormalData data, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIRICHLET_SAMPLE samples the Dirichlet PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Jerry Banks, editor,
        //    Handbook of Simulation,
        //    Engineering and Management Press Books, 1998, page 169.
        //
        //  Parameters:
        //
        //    Input, int N, the number of components.
        //
        //    Input, double A[N], the probabilities for each component.
        //    Each A(I) should be nonnegative, and at least one should be
        //    positive.
        //
        //    Input/output, int &SEED, a seed for the random 
        //    number generator.
        //
        //    Output, double DIRICHLET_SAMPLE[N], a sample of the PDF.  The entries 
        //    of X should sum to 1.
        //
    {
        dirichlet_check(n, a);

        double[] x = new double[n];

        for (int i = 0; i < n; i++)
        {
            x[i] = gamma_sample(a[i], 1.0, ref data, ref seed);
        }

        //
        //  Normalize the result.
        //
        double x_sum = typeMethods.r8vec_sum(n, x);

        for (int i = 0; i < n; i++)
        {
            x[i] /= x_sum;
        }

        return x;
    }

    public static double[] dirichlet_variance(int n, double[] a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIRICHLET_VARIANCE returns the variances of the Dirichlet PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of components.
        //
        //    Input, double A[N], the probabilities for each component.
        //    Each A[I] should be nonnegative, and at least one should be positive.
        //
        //    Output, double DIRICHLET_VARIANCE[N], the variances of the PDF.
        //
    {
        dirichlet_check(n, a);

        double a_sum = typeMethods.r8vec_sum(n, a);

        double[] variance = new double[n];

        for (int i = 0; i < n; i++)
        {
            variance[i] = a[i] * (a_sum - a[i]) / (a_sum * a_sum * (a_sum + 1.0));
        }

        return variance;
    }

    public static double exponential_01_sample(ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXPONENTIAL_01_SAMPLE samples the Exponential PDF with parameters 0, 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, int &SEED, a seed for the random 
        //    number generator.
        //
        //    Output, double EXPONENTIAL_01_SAMPLE, a sample of the PDF.
        //
    {
        double cdf = UniformRNG.r8_uniform_01(ref seed);

        double a = 0.0;
        double b = 1.0;
        double x = exponential_cdf_inv(cdf, a, b);

        return x;
    }

    public static double exponential_cdf_inv(double cdf, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXPONENTIAL_CDF_INV inverts the Exponential CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double CDF, the value of the CDF.
        //    0.0 <= CDF <= 1.0.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Output, double EXPONENTIAL_CDF_INV, the corresponding argument.
        //
    {
        switch (b)
        {
            //
            //  Check.
            //
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("EXPONENTIAL_CDF_INV - Fatal error!");
                Console.WriteLine("  B <= 0.0");
                break;
        }

        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine("");
                Console.WriteLine("EXPONENTIAL_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                break;
        }

        double x = a - b * Math.Log(1.0 - cdf);

        return x;
    }

    public static double gamma_sample(double a, double b, ref typeMethods.r8NormalData data, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GAMMA_SAMPLE samples the Gamma PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2013
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Joachim Ahrens, Ulrich Dieter.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Joachim Ahrens, Ulrich Dieter,
        //    Computer Methods for Sampling from Gamma, Beta, Poisson and
        //    Binomial Distributions,
        //    Computing,
        //    Volume 12, Number 3, September 1974, pages 223-246.
        //
        //    Joachim Ahrens, Ulrich Dieter,
        //    Generating Gamma Variates by a Modified Rejection Technique,
        //    Communications of the ACM,
        //    Volume 25, Number 1, January 1982, pages 47-54.
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < A,
        //    0.0 < B.
        //
        //    Input/output, int &SEED, a seed for the random 
        //    number generator.
        //
        //    Output, double GAMMA_SAMPLE, a sample of the PDF.
        //
    {
        const double a1 = 0.3333333;
        const double a2 = -0.2500030;
        const double a3 = 0.2000062;
        const double a4 = -0.1662921;
        const double a5 = 0.1423657;
        const double a6 = -0.1367177;
        const double a7 = 0.1233795;
        const double e1 = 1.0;
        const double e2 = 0.4999897;
        const double e3 = 0.1668290;
        const double e4 = 0.0407753;
        const double e5 = 0.0102930;
        const double q1 = 0.04166669;
        const double q2 = 0.02083148;
        const double q3 = 0.00801191;
        const double q4 = 0.00144121;
        const double q5 = -0.00007388;
        const double q6 = 0.00024511;
        const double q7 = 0.00024240;

        double e;
        double x;
        switch (a)
        {
            //
            //  Allow A = 0.
            //
            case 0.0:
                x = 0.0;
                return x;
            //
            //  A < 1.
            //
            case < 1.0:
            {
                for (;;)
                {
                    double p = UniformRNG.r8_uniform_01(ref seed);
                    p = (1.0 + 0.3678794 * a) * p;

                    e = exponential_01_sample(ref seed);

                    switch (p)
                    {
                        case >= 1.0:
                        {
                            x = -Math.Log((1.0 + 0.3678794 * a - p) / a);

                            if ((1.0 - a) * Math.Log(x) <= e)
                            {
                                x /= b;
                                return x;
                            }

                            break;
                        }
                        default:
                        {
                            x = Math.Exp(Math.Log(p) / a);

                            if (x <= e)
                            {
                                x /= b;
                                return x;
                            }

                            break;
                        }
                    }
                }
            }
            //
            default:
            {
                double s2 = a - 0.5;
                double s = Math.Sqrt(a - 0.5);
                double d = Math.Sqrt(32.0) - 12.0 * Math.Sqrt(a - 0.5);

                double t = typeMethods.r8_normal_01(ref data, ref seed);
                x = Math.Pow(Math.Sqrt(a - 0.5) + 0.5 * t, 2);

                switch (t)
                {
                    case >= 0.0:
                        x /= b;
                        return x;
                }

                double u = UniformRNG.r8_uniform_01(ref seed);

                if (d * u <= t * t * t)
                {
                    x /= b;
                    return x;
                }

                double r = 1.0 / a;
                double q0 = ((((((
                                     q7 * r
                                     + q6) * r
                                 + q5) * r
                                + q4) * r
                               + q3) * r
                              + q2) * r
                             + q1) * r;

                double bcoef;
                double c;
                double si;
                switch (a)
                {
                    case <= 3.686:
                        bcoef = 0.463 + s - 0.178 * s2;
                        si = 1.235;
                        c = 0.195 / s - 0.079 + 0.016 * s;
                        break;
                    case <= 13.022:
                        bcoef = 1.654 + 0.0076 * s2;
                        si = 1.68 / s + 0.275;
                        c = 0.062 / s + 0.024;
                        break;
                    default:
                        bcoef = 1.77;
                        si = 0.75;
                        c = 0.1515 / s;
                        break;
                }

                double q;
                double v;
                switch (Math.Sqrt(a - 0.5) + 0.5 * t)
                {
                    case > 0.0:
                    {
                        v = 0.5 * t / s;

                        q = Math.Abs(v) switch
                        {
                            > 0.25 => q0 - s * t + 0.25 * t * t + 2.0 * s2 * Math.Log(1.0 + v),
                            _ => q0 + 0.5 * t * t *
                                ((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v + a2) * v + a1) * v
                        };

                        if (Math.Log(1.0 - u) <= q)
                        {
                            x /= b;
                            return x;
                        }

                        break;
                    }
                }

                for (;;)
                {
                    e = exponential_01_sample(ref seed);

                    u = UniformRNG.r8_uniform_01(ref seed);
                    t = u switch
                    {
                        < 0.5 => bcoef - Math.Abs(si * e),
                        _ => bcoef + Math.Abs(si * e)
                    };

                    switch (t)
                    {
                        case >= -0.7187449:
                        {
                            v = 0.5 * t / s;

                            q = Math.Abs(v) switch
                            {
                                > 0.25 => q0 - s * t + 0.25 * t * t + 2.0 * s2 * Math.Log(1.0 + v),
                                _ => q0 + 0.5 * t * t *
                                    ((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v + a2) * v + a1) * v
                            };

                            switch (q)
                            {
                                case > 0.0:
                                {
                                    double w = q switch
                                    {
                                        > 0.5 => Math.Exp(q) - 1.0,
                                        _ => ((((e5 * q + e4) * q + e3) * q + e2) * q + e1) * q
                                    };

                                    if (c * Math.Abs(u) <= w * Math.Exp(e - 0.5 * t * t))
                                    {
                                        x = Math.Pow(s + 0.5 * t, 2) / b;
                                        return x;
                                    }

                                    break;
                                }
                            }

                            break;
                        }
                    }
                }
            }
        }
    }

    public static double ppnd16(double p, ref int ifault)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PPND16 inverts the standard normal CDF.
        //
        //  Discussion:
        //
        //    The result is accurate to about 1 part in 10**16.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 March 2010
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Michael Wichura.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Michael Wichura,
        //    The Percentage Points of the Normal Distribution,
        //    Algorithm AS 241,
        //    Applied Statistics,
        //    Volume 37, Number 3, pages 477-484, 1988.
        //
        //  Parameters:
        //
        //    Input, double P, the value of the cumulative probability 
        //    densitity function.  0 < P < 1.  If P is outside this range, an "infinite"
        //    value is returned.
        //
        //    Output, integer *IFAULT, error flag.
        //    0, no error.
        //    1, P <= 0 or P >= 1.  PPND is returned as 0.
        //
        //    Output, double PPND16, the normal deviate value 
        //    with the property that the probability of a standard normal deviate being 
        //    less than or equal to this value is P.
        //
    {
        double[] a =
        {
            3.3871328727963666080, 1.3314166789178437745e+2,
            1.9715909503065514427e+3, 1.3731693765509461125e+4,
            4.5921953931549871457e+4, 6.7265770927008700853e+4,
            3.3430575583588128105e+4, 2.5090809287301226727e+3
        };
        double[] b =
        {
            1.0, 4.2313330701600911252e+1,
            6.8718700749205790830e+2, 5.3941960214247511077e+3,
            2.1213794301586595867e+4, 3.9307895800092710610e+4,
            2.8729085735721942674e+4, 5.2264952788528545610e+3
        };
        double[] c =
        {
            1.42343711074968357734, 4.63033784615654529590,
            5.76949722146069140550, 3.64784832476320460504,
            1.27045825245236838258, 2.41780725177450611770e-1,
            2.27238449892691845833e-2, 7.74545014278341407640e-4
        };
        const double const1 = 0.180625;
        const double const2 = 1.6;
        double[] d =
        {
            1.0, 2.05319162663775882187,
            1.67638483018380384940, 6.89767334985100004550e-1,
            1.48103976427480074590e-1, 1.51986665636164571966e-2,
            5.47593808499534494600e-4, 1.05075007164441684324e-9
        };
        double[] e =
        {
            6.65790464350110377720, 5.46378491116411436990,
            1.78482653991729133580, 2.96560571828504891230e-1,
            2.65321895265761230930e-2, 1.24266094738807843860e-3,
            2.71155556874348757815e-5, 2.01033439929228813265e-7
        };
        double[] f =
        {
            1.0, 5.99832206555887937690e-1,
            1.36929880922735805310e-1, 1.48753612908506148525e-2,
            7.86869131145613259100e-4, 1.84631831751005468180e-5,
            1.42151175831644588870e-7, 2.04426310338993978564e-15
        };
        double r;
        const double split1 = 0.425;
        const double split2 = 5.0;
        double value = 0;

        ifault = 0;

        switch (p)
        {
            case <= 0.0:
                ifault = 1;
                value = -typeMethods.r8_huge();
                return value;
            case >= 1.0:
                ifault = 1;
                value = typeMethods.r8_huge();
                return value;
        }

        double q = p - 0.5;

        if (Math.Abs(q) <= split1)
        {
            r = const1 - q * q;
            value = q * typeMethods.r8poly_value(8, a, r) / typeMethods.r8poly_value(8, b, r);
        }
        else
        {
            r = q switch
            {
                < 0.0 => p,
                _ => 1.0 - p
            };

            value = r switch
            {
                <= 0.0 => -1.0,
                _ => value
            };

            r = Math.Sqrt(-Math.Log(r));

            if (r <= split2)
            {
                r -= const2;
                value = typeMethods.r8poly_value(8, c, r) / typeMethods.r8poly_value(8, d, r);
            }
            else
            {
                r -= split2;
                value = typeMethods.r8poly_value(8, e, r) / typeMethods.r8poly_value(8, f, r);
            }

            value = q switch
            {
                < 0.0 => -value,
                _ => value
            };
        }

        return value;
    }
}