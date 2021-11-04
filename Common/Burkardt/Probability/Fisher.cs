using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class Fisher
    {
        public static double fisher_pdf(double[] x, double kappa, double[] mu )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FISHER_PDF evaluates the Fisher PDF.
        //
        //  Discussion:
        //
        //    The formula for the PDF is:
        //
        //      PDF(KAPPA,MU;X) = C(KAPPA) * exp ( KAPPA * MU' * X )
        //
        //    where:
        //
        //      0 <= KAPPA is the concentration parameter,
        //      MU is a point on the unit sphere, the mean direction,
        //      X is any point on the unit sphere,
        //      and C(KAPPA) is a normalization factor:
        //
        //      C(KAPPA) = sqrt ( KAPPA ) / ( ( 2 * Math.PI )^(3/2) * I(0.5,KAPPA) )
        //
        //    where
        //
        //      I(nu,X) is the Bessel function of order NU and argument X.
        //
        //    For a fixed value of MU, the value of KAPPA determines the
        //    tendency of sample points to tend to be near MU.  In particular,
        //    KAPPA = 0 corresponds to a uniform distribution of points on the
        //    sphere, but as KAPPA increases, the sample points will tend to
        //    cluster more closely to MU.
        //
        //    The Fisher distribution for points on the unit sphere is
        //    analogous to the normal distribution of points on a line,
        //    and, more precisely, to the von Mises distribution on a circle.
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
        //    John Burkardt
        //
        //  Reference:
        //
        //    Kanti Mardia, Peter Jupp,
        //    Directional Statistics,
        //    Wiley, 2000,
        //    LC: QA276.M335
        //
        //  Parameters:
        //
        //    Input, double X[3], the argument of the PDF.
        //    X should have unit Euclidean norm, but this routine will
        //    automatically work with a normalized version of X.
        //
        //    Input, double KAPPA, the concentration parameter.
        //    0 <= KAPPA is required.
        //
        //    Input, double MU[3], the mean direction.
        //    MU should have unit Euclidean norm, but this routine will
        //    automatically work with a normalized version of MU.
        //
        //    Output, double FISHER_PDF, the value of the PDF.
        //
        {
            int NB = 1;

            double[] b = new double[NB];
            double pdf;
            

            if (kappa < 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("FISHER_PDF - Fatal error!");
                Console.WriteLine("  KAPPA must be nonnegative.");
                Console.WriteLine("  Input KAPPA = " + kappa + "");
                return (1);
            }

            if (kappa == 0.0)
            {
                pdf = 1.0 / (4.0 * Math.PI);
                return pdf;
            }

            double alpha = 0.5;
            int ize = 1;

            Ribesl.ribesl(kappa, alpha, NB, ize, ref b);

            double cf = Math.Sqrt(kappa) / (Math.Sqrt(Math.Pow(2.0 * Math.PI, 3)) * b[0]);

            double mu_norm = typeMethods.r8vec_length(3, mu);

            if (mu_norm == 0.0)
            {
                pdf = cf;
                return pdf;
            }

            double x_norm = typeMethods.r8vec_length(3, x);

            if (x_norm == 0.0)
            {
                pdf = cf;
                return pdf;
            }

            double arg = kappa * typeMethods.r8vec_dot(3, x, mu) / (x_norm * mu_norm);

            pdf = cf * Math.Exp(arg);

            return pdf;
        }

        public static double[] fisher_sample(double kappa, double[] mu, int n, ref int seed )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FISHER_SAMPLE samples the Fisher distribution.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 March 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Nicholas Fisher, Toby Lewis, Brian Embleton,
        //    Statistical Analysis of Spherical Data,
        //    Cambridge, 2003,
        //    ISBN13: 978-0521456999,
        //    LC: QA276.F489.
        //
        //  Parameters:
        //
        //    Input, double KAPPA, the concentration parameter.
        //
        //    Input, double MU[3], the mean direction.
        //    MU should have unit Euclidean norm, but this routine will
        //    automatically work with a normalized version of MU.
        //
        //    Input, int N, the number of samples to choose.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double FISHER_SAMPLE[3*N], a sample of the Fisher distribution.
        //
        //  Local Parameters:
        //
        //    Local, double ALPHA, BETA, the colatitude (theta) and
        //    longitude (phi) of the mean direction.
        //
        {
            double[] a = new double[3*3];
            double lambda;
            double mu_norm;
            double[] phi;
            
            double[] rst = new double[3];
            double[] theta;
            double[] xyz;

            mu_norm = typeMethods.r8vec_length(3, mu);

            if (mu_norm == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("FISHER_SAMPLE - Fatal error!");
                Console.WriteLine("  MU = 0.");
                return new double[0];
            }

            double alpha = -Math.Acos(mu[2] / mu_norm);
            double beta = Math.Atan2(mu[1], mu[0]);

            lambda = Math.Exp(-2.0 * kappa);

            theta = UniformRNG.r8vec_uniform_01_new(n, ref seed);

            for (int k = 0; k < n; k++)
            {
                if (kappa == 0.0)
                {
                    theta[k] = 2.0 * Math.Asin(Math.Sqrt(1.0 - theta[k]));
                }
                else
                {
                    theta[k] = 2.0 * Math.Asin(Math.Sqrt(
                        -Math.Log(theta[k] * (1.0 - lambda) + lambda)
                        / (2.0 * kappa)));
                }
            }

            phi = UniformRNG.r8vec_uniform_01_new(n, ref seed);

            for (int k = 0; k < n; k++)
            {
                phi[k] = 2.0 * Math.PI * phi[k];
            }

            //
            //  Compute the rotation matrix.
            //
            a[0 + 0 * 3] = Math.Cos(alpha) * Math.Cos(beta);
            a[1 + 0 * 3] = -Math.Sin(beta);
            a[2 + 0 * 3] = Math.Sin(alpha) * Math.Cos(beta);

            a[0 + 1 * 3] = Math.Cos(alpha) * Math.Sin(beta);
            a[1 + 1 * 3] = +Math.Cos(beta);
            a[2 + 1 * 3] = Math.Sin(alpha) * Math.Sin(beta);

            a[0 + 2 * 3] = -Math.Sin(alpha);
            a[1 + 2 * 3] = 0.0;
            a[2 + 2 * 3] = Math.Cos(alpha);
            //
            //  Compute the unrotated points.
            //
            xyz = new double[3 * n];

            for (int k = 0; k < n; k++)
            {
                rst[0] = Math.Sin(theta[k]) * Math.Cos(phi[k]);
                rst[1] = Math.Sin(theta[k]) * Math.Sin(phi[k]);
                rst[2] = Math.Cos(theta[k]);
                //
                //  Rotate the points.
                //
                for (int i = 0; i < 3; i++)
                {
                    xyz[i + k * 3] = 0.0;
                    for (int j = 0; j < 3; j++)
                    {
                        xyz[i + k * 3] = xyz[i + k * 3] + a[i + j * 3] * rst[j];
                    }
                }
            }

            return xyz;
        }
    }
}