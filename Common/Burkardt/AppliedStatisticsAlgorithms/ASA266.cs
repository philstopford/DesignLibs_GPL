﻿using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.AppliedStatistics
{
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
            double y;
            double z;

            if (x <= 0.0)
            {
                ifault = 1;
                value = 0.0;
                return value;
            }

            ifault = 0;
            y = x;

            if (x < 7.0)
            {
                f = 1.0;
                z = y;

                while (z < 7.0)
                {
                    f = f * z;
                    z = z + 1.0;
                }

                y = z;
                f = -Math.Log(f);
            }
            else
            {
                f = 0.0;
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

        static double digamma(double x)
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
            //
            //  Check the input.
            //
            if (x <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("DIGAMMA - Fatal error!");
                Console.WriteLine("  X <= 0.0.");
                return (1);
            }

            //
            //  Initialize.
            //
            double x2 = x;
            double value = 0.0;
            //
            //  Use approximation for small argument.
            //
            if (x2 <= 0.00001)
            {
                value = -euler_mascheroni - 1.0 / x2;
                return value;
            }

            //
            //  Reduce to DIGAMMA(X + N).
            //
            while (x2 < 8.5)
            {
                value = value - 1.0 / x2;
                x2 = x2 + 1.0;
            }

            //
            //  Use Stirling's (actually de Moivre's) expansion.
            //
            double r = 1.0 / x2;
            value = value + Math.Log(x2) - 0.5 * r;
            r = r * r;
            value = value
                    - r * (1.0 / 12.0
                           - r * (1.0 / 120.0
                                  - r * 1.0 / 252.0));

            return value;
        }

        static void dirichlet_check(int n, double[] a)
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
                if (a[i] < 0.0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("DIRICHLET_CHECK - Fatal error!");
                    Console.WriteLine("  A(I) < 0.");
                    Console.WriteLine("  For I = " + i + "");
                    Console.WriteLine("  A[I] = " + a[i] + "");
                }
                else if (0.0 < a[i])
                {
                    positive = 1;
                }
            }

            if (positive == 0)
            {
                Console.WriteLine("");
                Console.WriteLine("DIRICHLET_CHECK - Fatal error!");
                Console.WriteLine("  All parameters are zero!");
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
            double alpha_min = 0.00001;
            double gamma = 0.0001;
            double gg;
            int ifault2 = 0;
            int it_max = 100;
            double sum1;
            double sum2;

            ifault = 0;
            //
            //  Check the input arguments.
            //
            if (k < 2)
            {
                ifault = 1;
                Console.WriteLine("");
                Console.WriteLine("DIRICHLET_ESTIMATE - Fatal error!");
                Console.WriteLine("  K < 2.");
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
                    if (x[i + j * ix] <= 0.0)
                    {
                        niter = i;
                        ifault = 4;
                        Console.WriteLine("");
                        Console.WriteLine("DIRICHLET_ESTIMATE - Fatal error!");
                        Console.WriteLine("  X(I,J) <= 0.");
                    }
                }

                sum2 = 0.0;
                for (int j = 0; j < k; j++)
                {
                    sum2 = sum2 + x[i + j * ix];
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

            double an = (double) (n);
            double rk = (double) (k);
            niter = 0;
            //
            //  Calculate initial estimates using the method of moments.
            //
            if (init == 1)
            {
                sum2 = 0.0;
                for (int j = 0; j < k - 1; j++)
                {
                    sum1 = 0.0;
                    for (int i = 0; i < n; i++)
                    {
                        sum1 = sum1 + x[i + j * ix];
                    }

                    alpha[j] = sum1 / an;
                    sum2 = sum2 + alpha[j];
                }

                alpha[k - 1] = 1.0 - sum2;

                double x12 = 0.0;
                for (int i = 0; i < n; i++)
                {
                    x12 = x12 + Math.Pow(x[i + 0 * ix], 2);
                }

                x12 = x12 / an;
                double varp1 = x12 - Math.Pow(alpha[0], 2);

                double x11 = (alpha[0] - x12) / varp1;
                for (int j = 0; j < k; j++)
                {
                    alpha[j] = x11 * alpha[j];
                }
            }
            //
            //  Calculate initial estimates using Ronning's suggestion.
            //
            else if (init == 2)
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
            }

            //
            //  Check whether any ALPHA's are negative or zero.
            //
            for (int j = 0; j < k; j++)
            {
                if (alpha[j] <= 0.0)
                {
                    ifault = 6;
                    Console.WriteLine("");
                    Console.WriteLine("DIRICHLET_ESTIMATE - Fatal error!");
                    Console.WriteLine("  For J = " + j + "");
                    Console.WriteLine("  ALPHA(J) = " + alpha[j] + "");
                    Console.WriteLine("  but ALPHA(J) must be positive.");
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
                    work[j] = work[j] + Math.Log(x[i + j * ix]);
                }
            }

            //
            //  Call Algorithm AS 91 to compute CHI2, the chi-squared value.
            //
            gg = Helpers.LogGamma(rk / 2.0);

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
                    sum1 = sum1 + 1.0 / work2[j];
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
                            v[i + j * k] = v[i + j * k] + 1.0 / (an * work2[j]);
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
                    alpha[j] = alpha[j] + work2[j];
                    alpha[j] = Math.Max(alpha[j], alpha_min);
                }

                for (int j = 0; j < k; j++)
                {
                    if (alpha[j] <= 0.0)
                    {
                        ifault = 6;
                        Console.WriteLine("");
                        Console.WriteLine("DIRICHLET_ESTIMATE - Fatal error!");
                        Console.WriteLine("  Newton iteration " + it_num + "");
                        Console.WriteLine("  Computed ALPHA[J] <= 0.");
                        Console.WriteLine("  J = " + j + "");
                        Console.WriteLine("  ALPHA[J] = " + alpha[j] + "");
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

                    rlogl = rlogl + an * Helpers.LogGamma(sum2);

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
            double[] a_sum;
            double comp_weight_sum;
            double[] mean;
            //
            //  Check.
            //
            for (int comp_i = 0; comp_i < comp_num; comp_i++)
            {
                for (int elem_i = 0; elem_i < elem_num; elem_i++)
                {
                    if (a[comp_i + elem_i * comp_max] < 0.0)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("DIRICHLET_MIX_MEAN - Fatal error!");
                        Console.WriteLine("  A(COMP,ELEM) < 0.");
                        Console.WriteLine("  COMP = " + comp_i + "");
                        Console.WriteLine("  ELEM = " + elem_i + "");
                        Console.WriteLine("  A[COMP,ELEM] = " + a[comp_i + elem_i * comp_max] + "");
                    }
                }
            }

            comp_weight_sum = typeMethods.r8vec_sum(comp_num, comp_weight);

            a_sum = new double[comp_num];

            for (int comp_i = 0; comp_i < comp_num; comp_i++)
            {
                a_sum[comp_i] = 0.0;
                for (int elem_i = 0; elem_i < elem_num; elem_i++)
                {
                    a_sum[comp_i] = a_sum[comp_i] + a[comp_i + elem_i * comp_max];
                }
            }

            mean = new double[elem_num];

            for (int elem_i = 0; elem_i < elem_num; elem_i++)
            {
                mean[elem_i] = 0.0;
                for (int comp_i = 0; comp_i < comp_num; comp_i++)
                {
                    mean[elem_i] = mean[elem_i]
                                   + comp_weight[comp_i] * a[comp_i + elem_i * comp_max] / a_sum[comp_i];
                }

                mean[elem_i] = mean[elem_i] / comp_weight_sum;
            }

            return mean;
        }

        public static double[] dirichlet_mix_sample(int comp_max, int comp_num, int elem_num,
            double[] a, double[] comp_weight, ref int seed, ref int comp)
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
            double comp_weight_sum;
            double r;
            double sum2;
            double[] x;
            double x_sum;
            //
            //  Check.
            //
            for (int comp_i = 0; comp_i < comp_num; comp_i++)
            {
                for (int elem_i = 0; elem_i < elem_num; elem_i++)
                {
                    if (a[comp_i + elem_i * comp_max] < 0.0)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("DIRICHLET_MIX_SAMPLE - Fatal error!");
                        Console.WriteLine("  A(COMP,ELEM) < 0.'");
                        Console.WriteLine("  COMP = " + comp_i + "");
                        Console.WriteLine("  ELEM = " + elem_i + "");
                        Console.WriteLine("  A(COMP,ELEM) = " + a[comp_i + elem_i * comp_max] + "");
                    }
                }
            }

            //
            //  Choose a particular density MIX.
            //
            comp_weight_sum = typeMethods.r8vec_sum(comp_num, comp_weight);

            r = UniformRNG.r8_uniform_ab(0.0, comp_weight_sum, ref seed);

            comp = 0;
            sum2 = 0.0;

            while (comp < comp_num)
            {
                sum2 = sum2 + comp_weight[comp];

                if (r <= sum2)
                {
                    break;
                }

                comp = comp + 1;
            }

            //
            //  Sample density COMP.
            //
            x = new double[elem_num];

            for (int elem_i = 0; elem_i < elem_num; elem_i++)
            {
                x[elem_i] = gamma_sample(a[comp + elem_i * comp_max], 1.0, ref seed);
            }

            //
            //  Normalize the result.
            //
            x_sum = typeMethods.r8vec_sum(elem_num, x);

            for (int elem_i = 0; elem_i < elem_num; elem_i++)
            {
                x[elem_i] = x[elem_i] / x_sum;
            }

            return x;
        }

        public static double[] dirichlet_sample(int n, double[] a, ref int seed)
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
                x[i] = gamma_sample(a[i], 1.0, ref seed);
            }

            //
            //  Normalize the result.
            //
            double x_sum = typeMethods.r8vec_sum(n, x);

            for (int i = 0; i < n; i++)
            {
                x[i] = x[i] / x_sum;
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

        static double exponential_01_sample(ref int seed)
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

        static double exponential_cdf_inv(double cdf, double a, double b)
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
            //
            //  Check.
            //
            if (b <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("EXPONENTIAL_CDF_INV - Fatal error!");
                Console.WriteLine("  B <= 0.0");
            }

            if (cdf < 0.0 || 1.0 < cdf)
            {
                Console.WriteLine("");
                Console.WriteLine("EXPONENTIAL_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
            }

            double x = a - b * Math.Log(1.0 - cdf);

            return x;
        }

        public static double gamma_sample(double a, double b, ref int seed)
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
            //
            //  Allow A = 0.
            //
            if (a == 0.0)
            {
                x = 0.0;
                return x;
            }

            //
            //  A < 1.
            //
            if (a < 1.0)
            {
                for (;;)
                {
                    double p = UniformRNG.r8_uniform_01(ref seed);
                    p = (1.0 + 0.3678794 * a) * p;

                    e = exponential_01_sample(ref seed);

                    if (1.0 <= p)
                    {
                        x = -Math.Log((1.0 + 0.3678794 * a - p) / a);

                        if ((1.0 - a) * Math.Log(x) <= e)
                        {
                            x = x / b;
                            return x;
                        }
                    }
                    else
                    {
                        x = Math.Exp(Math.Log(p) / a);

                        if (x <= e)
                        {
                            x = x / b;
                            return x;
                        }
                    }
                }
            }
            //
            //  1 <= A.
            //
            else
            {
                double s2 = a - 0.5;
                double s = Math.Sqrt(a - 0.5);
                double d = Math.Sqrt(32.0) - 12.0 * Math.Sqrt(a - 0.5);

                double t = r8_normal_01(ref seed);
                x = Math.Pow(Math.Sqrt(a - 0.5) + 0.5 * t, 2);

                if (0.0 <= t)
                {
                    x = x / b;
                    return x;
                }

                double u = UniformRNG.r8_uniform_01(ref seed);

                if (d * u <= t * t * t)
                {
                    x = x / b;
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
                if (a <= 3.686)
                {
                    bcoef = 0.463 + s - 0.178 * s2;
                    si = 1.235;
                    c = 0.195 / s - 0.079 + 0.016 * s;
                }
                else if (a <= 13.022)
                {
                    bcoef = 1.654 + 0.0076 * s2;
                    si = 1.68 / s + 0.275;
                    c = 0.062 / s + 0.024;
                }
                else
                {
                    bcoef = 1.77;
                    si = 0.75;
                    c = 0.1515 / s;
                }

                double q;
                double v;
                if (0.0 < Math.Sqrt(a - 0.5) + 0.5 * t)
                {
                    v = 0.5 * t / s;

                    if (0.25 < Math.Abs(v))
                    {
                        q = q0 - s * t + 0.25 * t * t + 2.0 * s2 * Math.Log(1.0 + v);
                    }
                    else
                    {
                        q = q0 + 0.5 * t * t * ((((((
                                                        a7 * v
                                                        + a6) * v
                                                    + a5) * v
                                                   + a4) * v
                                                  + a3) * v
                                                 + a2) * v
                                                + a1) * v;
                    }

                    if (Math.Log(1.0 - u) <= q)
                    {
                        x = x / b;
                        return x;
                    }
                }

                for (;;)
                {
                    e = exponential_01_sample(ref seed);

                    u = UniformRNG.r8_uniform_01(ref seed);
                    if (u < 0.5)
                    {
                        t = bcoef - Math.Abs(si * e);
                    }
                    else
                    {
                        t = bcoef + Math.Abs(si * e);
                    }

                    if (-0.7187449 <= t)
                    {
                        v = 0.5 * t / s;

                        if (0.25 < Math.Abs(v))
                        {
                            q = q0 - s * t + 0.25 * t * t + 2.0 * s2 * Math.Log(1.0 + v);
                        }
                        else
                        {
                            q = q0 + 0.5 * t * t * ((((((
                                                            a7 * v
                                                            + a6) * v
                                                        + a5) * v
                                                       + a4) * v
                                                      + a3) * v
                                                     + a2) * v
                                                    + a1) * v;
                        }

                        if (0.0 < q)
                        {
                            double w;
                            if (0.5 < q)
                            {
                                w = Math.Exp(q) - 1.0;
                            }
                            else
                            {
                                w = ((((
                                           e5 * q
                                           + e4) * q
                                       + e3) * q
                                      + e2) * q
                                     + e1) * q;
                            }

                            if (c * Math.Abs(u) <= w * Math.Exp(e - 0.5 * t * t))
                            {
                                x = Math.Pow(s + 0.5 * t, 2) / b;
                                return x;
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
            double const1 = 0.180625;
            double const2 = 1.6;
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
            double split1 = 0.425;
            double split2 = 5.0;
            double value;

            ifault = 0;

            if (p <= 0.0)
            {
                ifault = 1;
                value = -typeMethods.r8_huge();
                return value;
            }

            if (1.0 <= p)
            {
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
                if (q < 0.0)
                {
                    r = p;
                }
                else
                {
                    r = 1.0 - p;
                }

                if (r <= 0.0)
                {
                    value = -1.0;
                }

                r = Math.Sqrt(-Math.Log(r));

                if (r <= split2)
                {
                    r = r - const2;
                    value = typeMethods.r8poly_value(8, c, r) / typeMethods.r8poly_value(8, d, r);
                }
                else
                {
                    r = r - split2;
                    value = typeMethods.r8poly_value(8, e, r) / typeMethods.r8poly_value(8, f, r);
                }

                if (q < 0.0)
                {
                    value = -value;
                }

            }

            return value;
        }
        
        public static double r8_gamma_log(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_GAMMA_LOG evaluates the logarithm of the gamma function.
        //
        //  Discussion:
        //
        //    This routine calculates the LOG(GAMMA) function for a positive real
        //    argument X.  Computation is based on an algorithm outlined in
        //    references 1 and 2.  The program uses rational functions that
        //    theoretically approximate LOG(GAMMA) to at least 18 significant
        //    decimal digits.  The approximation for X > 12 is from reference
        //    3, while approximations for X < 12.0 are similar to those in
        //    reference 1, but are unpublished.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 April 2013
        //
        //  Author:
        //
        //    Original FORTRAN77 version by William Cody, Laura Stoltz.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    William Cody, Kenneth Hillstrom,
        //    Chebyshev Approximations for the Natural Logarithm of the
        //    Gamma Function,
        //    Mathematics of Computation,
        //    Volume 21, Number 98, April 1967, pages 198-203.
        //
        //    Kenneth Hillstrom,
        //    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
        //    May 1969.
        //
        //    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
        //    Charles Mesztenyi, John Rice, Henry Thatcher,
        //    Christoph Witzgall,
        //    Computer Approximations,
        //    Wiley, 1968,
        //    LC: QA297.C64.
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the function.
        //
        //    Output, double R8_GAMMA_LOG, the value of the function.
        //
        {
            double[] c =
            {
                -1.910444077728E-03,
                8.4171387781295E-04,
                -5.952379913043012E-04,
                7.93650793500350248E-04,
                -2.777777777777681622553E-03,
                8.333333333333333331554247E-02,
                5.7083835261E-03
            };
            double d1 = -5.772156649015328605195174E-01;
            double d2 = 4.227843350984671393993777E-01;
            double d4 = 1.791759469228055000094023;
            double frtbig = 2.25E+76;
            int i;
            double[] p1 =
            {
                4.945235359296727046734888,
                2.018112620856775083915565E+02,
                2.290838373831346393026739E+03,
                1.131967205903380828685045E+04,
                2.855724635671635335736389E+04,
                3.848496228443793359990269E+04,
                2.637748787624195437963534E+04,
                7.225813979700288197698961E+03
            };
            double[] p2 =
            {
                4.974607845568932035012064,
                5.424138599891070494101986E+02,
                1.550693864978364947665077E+04,
                1.847932904445632425417223E+05,
                1.088204769468828767498470E+06,
                3.338152967987029735917223E+06,
                5.106661678927352456275255E+06,
                3.074109054850539556250927E+06
            };
            double[] p4 =
            {
                1.474502166059939948905062E+04,
                2.426813369486704502836312E+06,
                1.214755574045093227939592E+08,
                2.663432449630976949898078E+09,
                2.940378956634553899906876E+10,
                1.702665737765398868392998E+11,
                4.926125793377430887588120E+11,
                5.606251856223951465078242E+11
            };
            double[] q1 =
            {
                6.748212550303777196073036E+01,
                1.113332393857199323513008E+03,
                7.738757056935398733233834E+03,
                2.763987074403340708898585E+04,
                5.499310206226157329794414E+04,
                6.161122180066002127833352E+04,
                3.635127591501940507276287E+04,
                8.785536302431013170870835E+03
            };
            double[] q2 =
            {
                1.830328399370592604055942E+02,
                7.765049321445005871323047E+03,
                1.331903827966074194402448E+05,
                1.136705821321969608938755E+06,
                5.267964117437946917577538E+06,
                1.346701454311101692290052E+07,
                1.782736530353274213975932E+07,
                9.533095591844353613395747E+06
            };
            double[] q4 =
            {
                2.690530175870899333379843E+03,
                6.393885654300092398984238E+05,
                4.135599930241388052042842E+07,
                1.120872109616147941376570E+09,
                1.488613728678813811542398E+10,
                1.016803586272438228077304E+11,
                3.417476345507377132798597E+11,
                4.463158187419713286462081E+11
            };
            double res;
            double sqrtpi = 0.9189385332046727417803297;
            double xbig = 2.55E+305;
            double xinf = 1.79E+308;

            double y = x;

            if (0.0 < y && y <= xbig)
            {
                if (y <= double.Epsilon)
                {
                    res = -Math.Log(y);
                }
                //
                //  EPS < X <= 1.5.
                //
                else
                {
                    double corr;
                    double xden;
                    double xm2;
                    double xnum;
                    if (y <= 1.5)
                    {
                        double xm1;
                        if (y < 0.6796875)
                        {
                            corr = -Math.Log(y);
                            xm1 = y;
                        }
                        else
                        {
                            corr = 0.0;
                            xm1 = (y - 0.5) - 0.5;
                        }

                        if (y <= 0.5 || 0.6796875 <= y)
                        {
                            xden = 1.0;
                            xnum = 0.0;
                            for (i = 0; i < 8; i++)
                            {
                                xnum = xnum * xm1 + p1[i];
                                xden = xden * xm1 + q1[i];
                            }

                            res = corr + (xm1 * (d1 + xm1 * (xnum / xden)));
                        }
                        else
                        {
                            xm2 = (y - 0.5) - 0.5;
                            xden = 1.0;
                            xnum = 0.0;
                            for (i = 0; i < 8; i++)
                            {
                                xnum = xnum * xm2 + p2[i];
                                xden = xden * xm2 + q2[i];
                            }

                            res = corr + xm2 * (d2 + xm2 * (xnum / xden));
                        }
                    }
                    //
                    //  1.5 < X <= 4.0.
                    //
                    else if (y <= 4.0)
                    {
                        xm2 = y - 2.0;
                        xden = 1.0;
                        xnum = 0.0;
                        for (i = 0; i < 8; i++)
                        {
                            xnum = xnum * xm2 + p2[i];
                            xden = xden * xm2 + q2[i];
                        }

                        res = xm2 * (d2 + xm2 * (xnum / xden));
                    }
                    //
                    //  4.0 < X <= 12.0.
                    //
                    else if (y <= 12.0)
                    {
                        double xm4 = y - 4.0;
                        xden = -1.0;
                        xnum = 0.0;
                        for (i = 0; i < 8; i++)
                        {
                            xnum = xnum * xm4 + p4[i];
                            xden = xden * xm4 + q4[i];
                        }

                        res = d4 + xm4 * (xnum / xden);
                    }
                    //
                    //  Evaluate for 12 <= argument.
                    //
                    else
                    {
                        res = 0.0;

                        if (y <= frtbig)
                        {
                            res = c[6];
                            double ysq = y * y;
                            for (i = 0; i < 6; i++)
                            {
                                res = res / ysq + c[i];
                            }
                        }

                        res = res / y;
                        corr = Math.Log(y);
                        res = res + sqrtpi - 0.5 * corr;
                        res = res + y * (corr - 1.0);
                    }
                }
            }
            //
            //  Return for bad arguments.
            //
            else
            {
                res = xinf;
            }

            //
            //  Final adjustments and return.
            //
            return res;
        }


        static double r8_normal_01(ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_NORMAL_01 samples the standard normal probability distribution.
        //
        //  Discussion:
        //
        //    The standard normal probability distribution function (PDF) has
        //    mean 0 and standard deviation 1.
        //
        //    The Box-Muller method is used, which is efficient, but
        //    generates two values at a time.
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
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double R8_NORMAL_01, a normally distributed random value.
        //
        {
            double pi = 3.141592653589793;
            int seed2 = 0;
            int used = 0;
            double x;
            double y = 0.0;
            //
            //  If we've used an even number of values so far, generate two more, 
            //  return one, and save one.
            //
            if ((used % 2) == 0)
            {
                double r1 = UniformRNG.r8_uniform_01(ref seed);
                seed2 = seed;
                double r2 = UniformRNG.r8_uniform_01(ref seed2);

                x = Math.Sqrt(-2.0 * Math.Log(r1)) * Math.Cos(2.0 * pi * r2);
                y = Math.Sqrt(-2.0 * Math.Log(r1)) * Math.Sin(2.0 * pi * r2);
            }
            else
            {
                seed = seed2;
                x = y;
            }

            used = used + 1;

            return x;
        }

        public static double r8_psi(double xx)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_PSI evaluates the psi function.
        //
        //  Discussion:
        //
        //    This routine evaluates the logarithmic derivative of the
        //    Gamma function,
        //
        //      PSI(X) = d/dX ( GAMMA(X) ) / GAMMA(X)
        //             = d/dX LN ( GAMMA(X) )
        //
        //    for real X, where either
        //
        //      - XMAX1 < X < - XMIN, and X is not a negative integer,
        //
        //    or
        //
        //      XMIN < X.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 June 2013
        //
        //  Author:
        //
        //    Original FORTRAN77 version by William Cody, Anthony Strecok, Henry Thacher.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    William Cody, Anthony Strecok, Henry Thacher,
        //    Chebyshev Approximations for the Psi Function,
        //    Mathematics of Computation,
        //    Volume 27, Number 121, January 1973, pages 123-127.
        //
        //  Parameters:
        //
        //    Input, double XX, the argument of the psi function.
        //
        //    Output, double R8_PSI, the value of the psi function.  
        //    PSI is assigned the value 0 when the psi function is undefined.
        //
        {
            double den;
            double fourth = 0.25;
            double[] p1 =
            {
                4.5104681245762934160E-03,
                5.4932855833000385356E+00,
                3.7646693175929276856E+02,
                7.9525490849151998065E+03,
                7.1451595818951933210E+04,
                3.0655976301987365674E+05,
                6.3606997788964458797E+05,
                5.8041312783537569993E+05,
                1.6585695029761022321E+05
            };
            double[] p2 =
            {
                -2.7103228277757834192E+00,
                -1.5166271776896121383E+01,
                -1.9784554148719218667E+01,
                -8.8100958828312219821E+00,
                -1.4479614616899842986E+00,
                -7.3689600332394549911E-02,
                -6.5135387732718171306E-21
            };
            double piov4 = 0.78539816339744830962;
            double[] q1 =
            {
                9.6141654774222358525E+01,
                2.6287715790581193330E+03,
                2.9862497022250277920E+04,
                1.6206566091533671639E+05,
                4.3487880712768329037E+05,
                5.4256384537269993733E+05,
                2.4242185002017985252E+05,
                6.4155223783576225996E-08
            };
            double[] q2 =
            {
                4.4992760373789365846E+01,
                2.0240955312679931159E+02,
                2.4736979003315290057E+02,
                1.0742543875702278326E+02,
                1.7463965060678569906E+01,
                8.8427520398873480342E-01
            };
            double sgn;
            double three = 3.0;
            double upper;
            double value;
            double w;
            double x;
            double x01 = 187.0E+00;
            double x01d = 128.0E+00;
            double x02 = 6.9464496836234126266E-04;
            double xlarge = 2.04E+15;
            double xmax1 = 3.60E+16;
            double xmin1 = 5.89E-39;
            double xsmall = 2.05E-09;
            double z;

            x = xx;
            w = Math.Abs(x);
            double aug = 0.0;
            //
            //  Check for valid arguments, then branch to appropriate algorithm.
            //
            if (xmax1 <= -x || w < xmin1)
            {
                if (0.0 < x)
                {
                    value = -typeMethods.r8_huge();
                }
                else
                {
                    value = typeMethods.r8_huge();
                }

                return value;
            }

            //
            //  X < 0.5, use reflection formula: psi(1-x) = psi(x) + pi * cot(pi*x)
            //  Use 1/X for PI*COTAN(PI*X)  when  XMIN1 < |X| <= XSMALL.
            //
            if (x < 0.5)
            {
                if (w <= xsmall)
                {
                    aug = -1.0 / x;
                }
                //
                //  Argument reduction for cotangent.
                //
                else
                {
                    if (x < 0.0)
                    {
                        sgn = piov4;
                    }
                    else
                    {
                        sgn = -piov4;
                    }

                    w = w - (double) ((int) (w));
                    int nq = (int) (w * 4.0);
                    w = 4.0 * (w - (double) (nq) * fourth);
                    //
                    //  W is now related to the fractional part of 4.0 * X.
                    //  Adjust argument to correspond to values in the first
                    //  quadrant and determine the sign.
                    //
                    int n = nq / 2;

                    if (n + n != nq)
                    {
                        w = 1.0 - w;
                    }

                    z = piov4 * w;

                    if ((n % 2) != 0)
                    {
                        sgn = -sgn;
                    }

                    //
                    //  Determine the final value for  -pi * cotan(pi*x).
                    //
                    n = (nq + 1) / 2;
                    if ((n % 2) == 0)
                    {
                        //
                        //  Check for singularity.
                        //
                        if (z == 0.0)
                        {
                            if (0.0 < x)
                            {
                                value = -typeMethods.r8_huge();
                            }
                            else
                            {
                                value = typeMethods.r8_huge();
                            }

                            return value;
                        }

                        aug = sgn * (4.0 / Math.Tan(z));
                    }
                    else
                    {
                        aug = sgn * (4.0 * Math.Tan(z));
                    }
                }

                x = 1.0 - x;
            }

            //
            //  0.5 <= X <= 3.0.
            //
            if (x <= three)
            {
                den = x;
                upper = p1[0] * x;
                for (int i = 0; i < 7; i++)
                {
                    den = (den + q1[i]) * x;
                    upper = (upper + p1[i + 1]) * x;
                }

                den = (upper + p1[8]) / (den + q1[7]);
                x = (x - x01 / x01d) - x02;
                value = den * x + aug;
                return value;
            }

            //
            //  3.0 < X.
            //
            if (x < xlarge)
            {
                w = 1.0 / (x * x);
                den = w;
                upper = p2[0] * w;
                for (int i = 0; i < 5; i++)
                {
                    den = (den + q2[i]) * w;
                    upper = (upper + p2[i + 1]) * w;
                }

                aug = (upper + p2[6]) / (den + q2[5]) - 0.5 / x + aug;
            }

            value = aug + Math.Log(x);

            return value;
        }
        
   }
}