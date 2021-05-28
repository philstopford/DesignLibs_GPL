using System;
using Burkardt.Types;

namespace Burkardt.Probability
{
    public static class Dirichlet
    {
        public static bool dirichlet_check(int n, double[] a)
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
        //    30 October 2004
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
        //    Each A[I] should be positive.
        //
        //    Output, bool DIRICHLET_CHECK, is true if the parameters are legal.
        //
        {
            int i;

            bool positive = false;

            for (i = 0; i < n; i++)
            {
                if (a[i] <= 0.0)
                {
                    Console.WriteLine(" ");
                    Console.WriteLine("DIRICHLET_CHECK - Warning!");
                    Console.WriteLine("  A[" + i + "] <= 0.");
                    Console.WriteLine("  A[" + i + "] = " + a[i] + ".");
                    return false;
                }
                else if (0.0 < a[i])
                {
                    positive = true;
                }
            }

            if (!positive)
            {
                Console.WriteLine(" ");
                Console.WriteLine("DIRICHLET_CHECK - Warning!");
                Console.WriteLine("  All parameters are zero!");
                return false;
            }

            return true;
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
        //    14 October 2004
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
        //    Each A[I] should be positive.
        //
        //    Output, double DIRICHLET_MEAN[N], the means of the PDF.
        //
        {
            double[] mean = new double[n];

            for (int i = 0; i < n; i++)
            {
                mean[i] = a[i];
            }

            typeMethods.r8vec_unit_sum(n, ref mean);

            return mean;
        }

        public static double[] dirichlet_moment2(int n, double[] a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIRICHLET_MOMENT2 returns the second moments of the Dirichlet PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 October 2004
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
        //    Each A(I) should be positive.
        //
        //    Output, double DIRICHLET_MOMENT[N*N], the second moments of the PDF.
        //
        {
            double[] m2 = new double[n * n];

            double a_sum = typeMethods.r8vec_sum(n, a);

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (i == j)
                    {
                        m2[i + j * n] = a[i] * (a[i] + 1.0) / (a_sum * (a_sum + 1.0));
                    }
                    else
                    {
                        m2[i + j * n] = a[i] * a[j] / (a_sum * (a_sum + 1.0));
                    }
                }
            }

            return m2;
        }

        public static double dirichlet_pdf(double[] x, int n, double[] a )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIRICHLET_PDF evaluates the Dirichlet PDF.
        //
        //  Discussion:
        //
        //    PDF(N,A;X) = Product ( 1 <= I <= N ) X(I)^( A(I) - 1 )
        //      * Gamma ( A_SUM ) / A_PROD
        //
        //    where
        //
        //      0 < A(I) for all I;
        //      0 <= X(I) for all I;
        //      Sum ( 1 <= I <= N ) X(I) = 1;
        //      A_SUM = Sum ( 1 <= I <= N ) A(I).
        //      A_PROD = Product ( 1 <= I <= N ) Gamma ( A(I) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X(N), the argument of the PDF.  Each X(I) should
        //    be greater than 0.0, and the X(I)'s must add up to 1.0.
        //
        //    Input, int N, the number of components.
        //
        //    Input, double A(N), the probabilities for each component.
        //    Each A(I) should be positive.
        //
        //    Output, double PDF, the value of the PDF.
        //
        {
            double a_prod;
            double a_sum;
            int i;
            double pdf;
            double tol = 0.0001;
            double x_sum;

            for (i = 0; i < n; i++)
            {
                if (x[i] <= 0.0)
                {
                    Console.WriteLine(" ");
                    Console.WriteLine("DIRICHLET_PDF - Fatal error!");
                    Console.WriteLine("  X(I) <= 0.");
                    return (1);
                }
            }

            x_sum = typeMethods.r8vec_sum(n, x);

            if (tol < Math.Abs(x_sum - 1.0))
            {
                Console.WriteLine(" ");
                Console.WriteLine("DIRICHLET_PDF - Fatal error!");
                Console.WriteLine("  SUM X(I) =/= 1.");
                return (1);
            }

            a_sum = typeMethods.r8vec_sum(n, a);

            a_prod = 1.0;
            for (i = 0; i < n; i++)
            {
                a_prod = a_prod * Helpers.Gamma(a[i]);
            }

            pdf = Helpers.Gamma(a_sum) / a_prod;

            for (i = 0; i < n; i++)
            {
                pdf = pdf * Math.Pow(x[i], a[i] - 1.0);
            }

            return pdf;
        }

        public static double[] dirichlet_sample(int n, double[] a, ref int seed )
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
        //    18 October 2004
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
        //    Input, double A(N), the probabilities for each component.
        //    Each A(I) should be positive.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double DIRICHLET_SAMPLE[N], a sample of the PDF.  The entries
        //    of X should sum to 1.
        //
        {
            double[] x = new double[n];

            double a2 = 0.0;
            double b2 = 1.0;

            for (int i = 0; i < n; i++)
            {
                double c2 = a[i];
                x[i] = gamma_sample(a2, b2, c2, seed);
            }

            //
            //  Rescale the vector to have unit sum.
            //
            typeMethods.r8vec_unit_sum(n, ref x);

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
        //    18 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of components.
        //
        //    Input, double A(N), the probabilities for each component.
        //    Each A(I) should be positive.
        //
        //    Output, double DIRICHLET_VARIANCE(N), the variances of the PDF.
        //
        {
            double[] variance = new double[n];

            double a_sum = typeMethods.r8vec_sum(n, a);

            for (int i = 0; i < n; i++)
            {
                variance[i] = a[i] * (a_sum - a[i]) / (a_sum * a_sum * (a_sum + 1.0));
            }

            return variance;
        }

        public static bool dirichlet_mix_check(int comp_num, int elem_num, double[] a,
        double[] comp_weight )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIRICHLET_MIX_CHECK checks the parameters of a Dirichlet mixture PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int COMP_NUM, the number of components in the Dirichlet
        //    mixture density, that is, the number of distinct Dirichlet PDF's
        //    that are mixed together.
        //
        //    Input, int ELEM_NUM, the number of elements of an observation.
        //
        //    Input, double A[ELEM_NUM*COMP_NUM], the probabilities
        //    for element ELEM_NUM in component COMP_NUM.
        //    Each A[I,J] should be positive.
        //
        //    Input, double COMP_WEIGHT[COMP_NUM], the mixture weights of the densities.
        //    These do not need to be normalized.  The weight of a given component is
        //    the relative probability that that component will be used to generate
        //    the sample.
        //
        //    Output, bool DIRICHLET_MIX_CHECK, is true if the parameters are legal.
        //
        {
            for (int comp_i = 0; comp_i < comp_num; comp_i++)
            {
                for (int elem_i = 0; elem_i < elem_num; elem_i++)
                {
                    if (a[elem_i + comp_i * elem_num] <= 0.0)
                    {
                        Console.WriteLine(" ");
                        Console.WriteLine("DIRICHLET_MIX_CHECK - Warning!");
                        Console.WriteLine("  A(ELEM,COMP) <= 0.");
                        Console.WriteLine("  COMP = " + comp_i + "");
                        Console.WriteLine("  ELEM = " + elem_i + "");
                        Console.WriteLine("  A(COMP,ELEM) = " + a[elem_i + comp_i * elem_num] + "");
                        return false;
                    }
                }
            }

            bool positive = false;

            for (int comp_i = 0; comp_i < comp_num; comp_i++)
            {
                if (comp_weight[comp_i] < 0.0)
                {
                    Console.WriteLine(" ");
                    Console.WriteLine("DIRICHLET_MIX_CHECK - Warning!");
                    Console.WriteLine("  COMP_WEIGHT(COMP) < 0.");
                    Console.WriteLine("  COMP = " + comp_i + "");
                    Console.WriteLine("  COMP_WEIGHT(COMP) = " + comp_weight[comp_i] + "");
                    return false;
                }
                else if (0.0 < comp_weight[comp_i])
                {
                    positive = true;
                }
            }

            if (!positive)
            {
                Console.WriteLine(" ");
                Console.WriteLine("DIRICHLET_MIX_CHECK - Warning!");
                Console.WriteLine("  All component weights are zero.");
                return false;
            }

            return true;
        }

        public static double[] dirichlet_mix_mean(int comp_num, int elem_num, double[] a,
        double[] comp_weight )
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
        //    22 November 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int COMP_NUM, the number of components in the Dirichlet
        //    mixture density, that is, the number of distinct Dirichlet PDF's
        //    that are mixed together.
        //
        //    Input, int ELEM_NUM, the number of elements of an observation.
        //
        //    Input, double A[ELEM_NUM*COMP_NUM], the probabilities for
        //    element ELEM_NUM in component COMP_NUM.
        //    Each A[I,J] should be positive.
        //
        //    Input, double COMP_WEIGHT[COMP_NUM], the mixture weights of the densities.
        //    These do not need to be normalized.  The weight of a given component is
        //    the relative probability that that component will be used to generate
        //    the sample.
        //
        //    Output, double DIRICHLET_MIX_MEAN[ELEM_NUM], the means for each element.
        //
        {
            double[] a_vec = new double[elem_num];
            double[] mean = new double[elem_num];

            double comp_weight_sum = 0.0;
            for (int comp_i = 0; comp_i < comp_num; comp_i++)
            {
                comp_weight_sum = comp_weight_sum + comp_weight[comp_i];
            }

            for (int elem_i = 0; elem_i < elem_num; elem_i++)
            {
                mean[elem_i] = 0.0;
            }

            for (int comp_i = 0; comp_i < comp_num; comp_i++)
            {
                for (int elem_i = 0; elem_i < elem_num; elem_i++)
                {
                    a_vec[elem_i] = a[elem_i + comp_i * elem_num];
                }

                double[] comp_mean = dirichlet_mean(elem_num, a_vec);
                for (int elem_i = 0; elem_i < elem_num; elem_i++)
                {
                    mean[elem_i] = mean[elem_i] + comp_weight[comp_i] * comp_mean[elem_i];
                }
            }

            for (int elem_i = 0; elem_i < elem_num; elem_i++)
            {
                mean[elem_i] = mean[elem_i] / comp_weight_sum;
            }
            return mean;
        }

        public static double dirichlet_mix_pdf(double[] x, int comp_num, int elem_num, double[] a,
        double[] comp_weight )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIRICHLET_MIX_PDF evaluates a Dirichlet mixture PDF.
        //
        //  Discussion:
        //
        //    The PDF is a weighted sum of Dirichlet PDF's.  Each PDF is a
        //    "component", with an associated weight.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X[ELEM_NUM], the argument of the PDF.
        //
        //    Input, int COMP_NUM, the number of components in the Dirichlet
        //    mixture density, that is, the number of distinct Dirichlet PDF's
        //    that are mixed together.
        //
        //    Input, int ELEM_NUM, the number of elements of an observation.
        //
        //    Input, double A[ELEM_NUM*COMP_NUM], the probabilities for
        //    element ELEM_NUM in component COMP_NUM.
        //    Each A[I,J] should be positive.
        //
        //    Input, double COMP_WEIGHT[COMP_NUM], the mixture weights of the densities.
        //    These do not need to be normalized.  The weight of a given component is
        //    the relative probability that that component will be used to generate
        //    the sample.
        //
        //    Output, double DIRICHLET_MIX_PDF, the value of the PDF.
        //
        {
            double[] a_vec = new double[elem_num];

            double comp_weight_sum = 0.0;
            for (int j = 0; j < comp_num; j++)
            {
                comp_weight_sum = comp_weight_sum + comp_weight[j];
            }

            double pdf = 0.0;
            for (int j = 0; j < comp_num; j++)
            {
                for (int i = 0; i < elem_num; i++)
                {
                    a_vec[i] = a[i + j * elem_num];
                }

                double comp_pdf = dirichlet_pdf(x, elem_num, a_vec);

                pdf = pdf + comp_weight[j] * comp_pdf / comp_weight_sum;
            }
            return pdf;
        }

        public static double[] dirichlet_mix_sample(int comp_num, int elem_num, double[] a,
        double[] comp_weight, ref int seed, ref int comp )
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
        //    30 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int COMP_NUM, the number of components in the Dirichlet
        //    mixture density, that is, the number of distinct Dirichlet PDF's
        //    that are mixed together.
        //
        //    Input, int ELEM_NUM, the number of elements of an observation.
        //
        //    Input, double A[ELEM_NUM*COMP_NUM], the probabilities for
        //    element ELEM_NUM in component COMP_NUM.
        //    Each A[I,J] should be positive.
        //
        //    Input, double COMP_WEIGHT[COMP_NUM], the mixture weights of the densities.
        //    These do not need to be normalized.  The weight of a given component is
        //    the relative probability that that component will be used to generate
        //    the sample.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, int *COMP, the index of the component of the Dirichlet
        //    mixture that was chosen to generate the sample.
        //
        //    Output, double DIRICHLET_MIX_SAMPLE[ELEM_NUM], a sample of the PDF.
        //
        {
            double[] a_vec = new double[elem_num];
            //
            //  Choose a particular density component COMP.
            //
            comp = Discrete.discrete_sample(comp_num, comp_weight, ref seed);
            //
            //  Sample the density number COMP.
            //
            for (int elem_i = 0; elem_i < elem_num; elem_i++)
            {
                a_vec[elem_i] = a[elem_i + (comp - 1) * elem_num];
            }

            double[] x = dirichlet_sample(elem_num, a_vec, ref seed);
            
            return x;
        }

        public static double dirichlet_multinomial_pdf(int[] x, int a, int b, double[] c )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIRICHLET_MULTINOMIAL_PDF evaluates a Dirichlet Multinomial PDF.
        //
        //  Formula:
        //
        //    PDF(A,B,C;X) = Comb(A,B,X) * ( Gamma(C_Sum) / Gamma(C_Sum+A) )
        //      Product ( 1 <= I <= B ) Gamma(C(I)+X(I)) / Gamma(C(I))
        //
        //    where:
        //
        //      Comb(A,B,X) is the multinomial coefficient C( A; X(1), X(2), ..., X(B) ),
        //      C_Sum = Sum ( 1 <= I <= B ) C(I)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Kenneth Lange,
        //    Mathematical and Statistical Methods for Genetic Analysis,
        //    Springer, 1997, page 45.
        //
        //  Parameters:
        //
        //    Input, int X[B]; X[I] counts the number of occurrences of
        //    outcome I, out of the total of A trials.
        //
        //    Input, int A, the total number of trials.
        //
        //    Input, int B, the number of different possible outcomes on
        //    one trial.
        //
        //    Input, double C[B]; C[I] is the Dirichlet parameter associated
        //    with outcome I.
        //
        //    Output, double DIRICHLET_MULTINOMIAL_PDF, the value of the Dirichlet
        //     multinomial PDF.
        //
        {
            double c_sum;
            int i;
            double pdf;
            double pdf_log;

            c_sum = typeMethods.r8vec_sum(b, c);

            pdf_log = -Helpers.LogGamma(c_sum + (double) (a)) + Helpers.LogGamma(c_sum)
                                                    + Helpers.LogGamma((double) (a + 1));

            for (i = 0; i < b; i++)
            {
                pdf_log = pdf_log + Helpers.LogGamma(c[i] + (double) (x[i]))
                          - Helpers.LogGamma(c[i]) - Helpers.LogGamma((double) (x[i] + 1));
            }

            pdf = Math.Exp(pdf_log);

            return pdf;
        }
    }
}