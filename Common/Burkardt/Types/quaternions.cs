using System;
using Burkardt.Uniform;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double[] q8_conjugate(double[] q1)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    Q8_CONJUGATE conjugates a quaternion.
            //
            //  Discussion:
            //
            //    A quaternion is a quadruplet (A,B,C,D) of real numbers, which
            //    may be written as
            //
            //      Q = A + Bi + Cj + Dk.
            //
            //    The conjugate of Q is
            //
            //      conj ( Q ) = A - Bi - Cj - Dk.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 May 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double Q1[4], the quaternion to be conjugated.
            //
            //    Output, double Q8_CONJUGATE[4], the conjugated quaternion.
            //
        {
            double[] q2;

            q2 = new double[4];

            q2[0] = q1[0];
            q2[1] = -q1[1];
            q2[2] = -q1[2];
            q2[3] = -q1[3];

            return q2;
        }

        public static double[] q8_exponentiate(double[] q1)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    Q8_EXPONENTIATE exponentiates a quaternion.
            //
            //  Discussion:
            //
            //    A quaternion is a quadruplet (A,B,C,D) of real numbers, which
            //    may be written as
            //
            //      Q = A + Bi + Cj + Dk.
            //    
            //    The exponential of Q can be set by
            //      V = Math.Sqrt ( B^2 + C^2 + D^2 )
            //      e^Q = e^A * ( Math.Cos ( ||V|| ) + V/||V|| Math.Sin ||V|| )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 August 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double Q1[4], the quaternions to exponentiate.
            //
            //    Output, double Q8_EXPONENTIATE[4], the exponential of the quaternion.
            //
        {
            int i;
            double[] q2;
            double v_norm;

            v_norm = Math.Sqrt(q1[1] * q1[1] + q1[2] * q1[2] + q1[3] * q1[3]);

            q2 = new double[4];
            q2[0] = Math.Cos(v_norm);
            if (v_norm != 0.0)
            {
                for (i = 1; i < 4; i++)
                {
                    q2[i] = Math.Sin(v_norm) * q1[i] / v_norm;
                }
            }
            else
            {
                for (i = 1; i < 4; i++)
                {
                    q2[i] = 0.0;
                }
            }

            for (i = 0; i < 4; i++)
            {
                q2[i] = Math.Exp(q1[0]) * q2[i];
            }

            return q2;
        }

        public static double[] q8_inverse(double[] q)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    Q8_INVERSE inverts a quaternion.
            //
            //  Discussion:
            //
            //    A quaternion is a quadruplet (A,B,C,D) of real numbers, which
            //    may be written as
            //
            //      Q = A + Bi + Cj + Dk.
            //
            //    The inverse of Q is
            //
            //      inverse ( Q ) = conjugate ( Q ) / ( norm ( Q ) )^2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 May 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double Q[4], the quaternion to be inverted.
            //
            //    Output, double QUAT_INVERSE[4], the inverse of the input quaternion.
            //
        {
            double norm;
            double[] q2;

            q2 = new double[4];

            norm = q[0] * q[0]
                   + q[1] * q[1]
                   + q[2] * q[2]
                   + q[3] * q[3];

            q2[0] = q[0] / norm;
            q2[1] = -q[1] / norm;
            q2[2] = -q[2] / norm;
            q2[3] = -q[3] / norm;

            return q2;
        }

        public static double[] q8_multiply(double[] q1, double[] q2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    Q8_MULTIPLY multiplies two quaternions.
            //
            //  Discussion:
            //
            //    A quaternion is a quadruplet (A,B,C,D) of real numbers, which
            //    may be written as
            //
            //      Q = A + Bi + Cj + Dk.
            //
            //    To multiply two quaternions, use the relationships:
            //
            //      i * j = -j * i = k
            //      j * k = -k * j = i
            //      k * i = -i * k = j
            //      i * i =  j * j = k * k = -1
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 May 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double Q1[4], Q2[4], the two quaternions to be multiplied.
            //
            //    Output, double Q8_MULTIPLY[4], the product of the two quaternions.
            //
        {
            double[] q3;

            q3 = new double[4];

            q3[0] = q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3];
            q3[1] = q1[0] * q2[1] + q1[1] * q2[0] + q1[2] * q2[3] - q1[3] * q2[2];
            q3[2] = q1[0] * q2[2] - q1[1] * q2[3] + q1[2] * q2[0] + q1[3] * q2[1];
            q3[3] = q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1] + q1[3] * q2[0];

            return q3;
        }

        public static double[] q8_multiply2(double[] q1, double[] q2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    Q8_MULTIPLY2 multiplies two quaternions using a matrix format.
            //
            //  Discussion:
            //
            //    A quaternion is a quadruplet (A,B,C,D) of real numbers, which
            //    may be written as
            //
            //      Q = A + Bi + Cj + Dk.
            //
            //    To multiply two quaternions, use the relationships:
            //
            //      i * j = -j * i = k
            //      j * k = -k * j = i
            //      k * i = -i * k = j
            //      i * i =  j * j = k * k = -1
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 August 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double Q1[4], Q2[4], the quaternions to be multiplied.
            //
            //    Output, double Q8_MULTIPLY2[4], the product of the two quaternions.
            //
        {
            double[] q3;
            //
            //  The matrix entries are listed by column, not row.
            //
            double[] qm =
            {
                q1[0], q1[1], q1[2], q1[3],
                -q1[1], +q1[0], +q1[3], -q1[2],
                -q1[2], -q1[3], +q1[0], +q1[1],
                -q1[3], +q1[2], -q1[1], +q1[0]
            };

            q3 = r8mat_mv_new(4, 4, qm, q2);

            return q3;
        }

        public static double q8_norm(double[] q)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    Q8_NORM computes the norm of a quaternion.
            //
            //  Discussion:
            //
            //    A quaternion is a quadruplet (A,B,C,D) of real numbers, which
            //    may be written as
            //
            //      Q = A + Bi + Cj + Dk.
            //
            //    The norm of Q is
            //
            //      norm(Q) = Math.Sqrt ( A * A + B * B + C * C + D * D ).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 May 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double Q[4], the quaternion.
            //
            //    Output, double Q8_NORM, the norm of the quaternion.
            //
        {
            double norm;

            norm = r8vec_norm(4, q);

            return norm;
        }

        public static double[] q8_normal_01(ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    Q8_NORMAL_01 returns a normally distributed quaternion.
            //
            //  Discussion:
            //
            //    The normal distribution with mean 0 and variance 1 is used.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 August 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input/output, int &SEED, a seed for the random number 
            //    generator.
            //
            //    Output, double[] Q8_NORMAL_01[4], the sampled quaternion.
            //
        {
            double[] q;
            double[] r;

            r = UniformRNG.r8vec_uniform_01_new(4, ref seed);

            q = new double[4];

            q[0] = Math.Sqrt(-2.0 * Math.Log(r[0])) * Math.Cos(2.0 * Math.PI * r[1]);
            q[1] = Math.Sqrt(-2.0 * Math.Log(r[0])) * Math.Sin(2.0 * Math.PI * r[1]);
            q[2] = Math.Sqrt(-2.0 * Math.Log(r[2])) * Math.Cos(2.0 * Math.PI * r[3]);
            q[3] = Math.Sqrt(-2.0 * Math.Log(r[2])) * Math.Sin(2.0 * Math.PI * r[3]);

            return q;
        }

        public static void q8_transpose_print(double[] q, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    Q8_TRANSPOSE_PRINT prints a Q8 "transposed".
            //
            //  Discussion:
            //
            //    A Q8 is a quaternion using R8 arithmetic.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 August 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double Q[4], the quaternion to be printed.
            //
            //    Input, string TITLE, a title.
            //
        {
            Console.WriteLine(title
                              + "  " + q[0]
                              + "  " + q[1]
                              + "  " + q[2]
                              + "  " + q[3] + "");

            return;
        }
    }
}