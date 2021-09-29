using System;
using System.Linq;
using Burkardt.BLAS;
using Burkardt.Uniform;

namespace BLAS1DTest
{
    class Program
    {
        static void Main(string[] args)
//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BLAS1_D_TEST.
//
//  Discussion:
//
//    BLAS1_D_TEST tests the BLAS1_D library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2006
//
//  Author:
//
//    John Burkardt
//
        {
            Console.WriteLine("");
            Console.WriteLine("BLAS1_D_TEST:");

            Console.WriteLine("  Test the BLAS1_D library.");

            dasum_test();
            daxpy_test();
            dcopy_test();
            ddot_test();
            dnrm2_test();
            drot_test();
            drotg_test();
            dscal_test();
            dswap_test();
            idamax_test();

            Console.WriteLine("");
            Console.WriteLine("BLAS1_D_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void dasum_test()

//****************************************************************************80
//
//  Purpose:
//
//    DASUM_TEST tests DASUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2006
//
//  Author:
//
//    John Burkardt
//
        {
            int LDA = 6;
            int MA = 5;
            int NA = 4;
            int NX = 10;

            double[] a = new double[LDA * NA];
            double[] x = new double[NX];

            for (int i = 0; i < NX; i++)
            {
                x[i] = Math.Pow(-1.0, i + 1) * (double) (2 * (i + 1));
            }

            Console.WriteLine("");
            Console.WriteLine("DASUM_TEST");
            Console.WriteLine("  DASUM adds the absolute values of elements of a vector.");
            Console.WriteLine("");
            Console.WriteLine("  X = ");
            Console.WriteLine("");
            for (int i = 0; i < NX; i++)
            {
                Console.WriteLine("  "
                     + (i + 1).ToString().PadLeft(6) + "  "
                     + x[i].ToString().PadLeft(14) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("  DASUM ( NX,   X, 1 ) =    " + BLAS1D.dasum(NX, x, 1,0) + "");
            Console.WriteLine("  DASUM ( NX/2, X, 2 ) =    " + BLAS1D.dasum(NX / 2, x, 2,0) + "");
            Console.WriteLine("  DASUM ( 2,    X, NX/2 ) = " + BLAS1D.dasum(2, x, NX / 2,0) + "");

            for (int i = 0; i < MA; i++)
            {
                for (int j = 0; j < NA; j++)
                {
                    a[i + j * LDA] = Math.Pow(-1.0, i + 1 + j + 1)
                                     * (double) (10 * (i + 1) + j + 1);
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  Demonstrate with a matrix A:");
            Console.WriteLine("");
            for (int i = 0; i < MA; i++)
            {
                string cout = "";
                for (int j = 0; j < NA; j++)
                {
                    cout += "  " + a[i + j * LDA].ToString().PadLeft(14);
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("  DASUM(MA,A(1,2),1) =   " + BLAS1D.dasum(MA, a, 1, ( 0 + 1 * LDA)) + "");
            Console.WriteLine("  DASUM(NA,A(2,1),LDA) = " + BLAS1D.dasum(NA, a, LDA, (1 + 0 * LDA)) + "");

        }

        static void daxpy_test()

//****************************************************************************80
//
//  Purpose:
//
//    DAXPY_TEST tests DAXPY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2006
//
//  Author:
//
//    John Burkardt
//
        {
            int N = 6;

            double da;
            double[] x = new double[N];
            double[] y = new double[N];

            for (int i = 0; i < N; i++)
            {
                x[i] = (double) (i + 1);
            }

            for (int i = 0; i < N; i++)
            {
                y[i] = (double) (100 * (i + 1));
            }

            Console.WriteLine("");
            Console.WriteLine("DAXPY_TEST");
            Console.WriteLine("  DAXPY adds a multiple of vector X to vector Y.");
            Console.WriteLine("");
            Console.WriteLine("  X =");
            Console.WriteLine("");
            for (int i = 0; i < N; i++)
            {
                Console.WriteLine("  "
                     + (i + 1).ToString().PadLeft(6) + "  "
                     + x[i].ToString().PadLeft(14) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("  Y =");
            Console.WriteLine("");
            for (int i = 0; i < N; i++)
            {
                Console.WriteLine("  "
                     + (i + 1).ToString().PadLeft(6) + "  "
                     + y[i].ToString().PadLeft(14) + "");
            }

            da = 1.0;
            BLAS1D.daxpy(N, da, x,1, ref y, 1);
            Console.WriteLine("");
            Console.WriteLine("  DAXPY ( N, " + da + ", X, 1, Y, 1 )");
            Console.WriteLine("");
            Console.WriteLine("  Y =");
            Console.WriteLine("");
            for (int i = 0; i < N; i++)
            {
                Console.WriteLine("  "
                     + (i + 1).ToString().PadLeft(6) + "  "
                     + y[i].ToString().PadLeft(14) + "");
            }

            for (int i = 0; i < N; i++)
            {
                y[i] = (double) (100 * (i + 1));
            }

            da = -2.0;
            BLAS1D.daxpy(N, da, x, 1, ref y, 1);
            Console.WriteLine("");
            Console.WriteLine("  DAXPY ( N, " + da + ", X, 1, Y, 1 )");
            Console.WriteLine("");
            Console.WriteLine("  Y =");
            Console.WriteLine("");
            for (int i = 0; i < N; i++)
            {
                Console.WriteLine("  "
                     + (i + 1).ToString().PadLeft(6) + "  "
                     + y[i].ToString().PadLeft(14) + "");
            }

            for (int i = 0; i < N; i++)
            {
                y[i] = (double) (100 * (i + 1));
            }

            da = 3.0;
            BLAS1D.daxpy(3, da, x, 2, ref y, 1);
            Console.WriteLine("");
            Console.WriteLine("  DAXPY ( 3, " + da + ", X, 2, Y, 1 )");
            Console.WriteLine("");
            Console.WriteLine("  Y =");
            Console.WriteLine("");
            for (int i = 0; i < N; i++)
            {
                Console.WriteLine("  "
                     + (i + 1).ToString().PadLeft(6) + "  "
                     + y[i].ToString().PadLeft(14) + "");
            }

            for (int i = 0; i < N; i++)
            {
                y[i] = (double) (100 * (i + 1));
            }

            da = -4.0;
            BLAS1D.daxpy(3, da, x, 1, ref y, 2);
            Console.WriteLine("");
            Console.WriteLine("  DAXPY ( 3, " + da + ", X, 1, Y, 2 )");
            Console.WriteLine("");
            Console.WriteLine("  Y =");
            Console.WriteLine("");
            for (int i = 0; i < N; i++)
            {
                Console.WriteLine("  "
                     + (i + 1).ToString().PadLeft(6) + "  "
                     + y[i].ToString().PadLeft(14) + "");
            }
        }

        static void dcopy_test()

//****************************************************************************80
//
//  Purpose:
//
//    DCOPY_TEST demonstrates DCOPY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2006
//
//  Author:
//
//    John Burkardt
//
        {
            double[] a = new double[5 * 5];
            double[] x = new double[10];
            double[] y = new double[10];

            Console.WriteLine("");
            Console.WriteLine("DCOPY_TEST");
            Console.WriteLine("  DCOPY copies one vector into another.");
            Console.WriteLine("");

            for (int i = 0; i < 10; i++)
            {
                x[i] = (double) (i + 1);
            }

            for (int i = 0; i < 10; i++)
            {
                y[i] = (double) (10 * (i + 1));
            }

            for (int i = 0; i < 5; i++)
            {
                for (int j = 0; j < 5; j++)
                {
                    a[i + j * 5] = (double) (10 * (i + 1) + j + 1);
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  X =");
            Console.WriteLine("");
            for (int i = 0; i < 10; i++)
            {
                Console.WriteLine("  "
                     + (i + 1).ToString().PadLeft(6) + "  "
                     + x[i].ToString().PadLeft(14) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("  Y =");
            Console.WriteLine("");
            for (int i = 0; i < 10; i++)
            {
                Console.WriteLine("  "
                     + (i + 1).ToString().PadLeft(6) + "  "
                     + y[i].ToString().PadLeft(14) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("  A =");
            Console.WriteLine("");
            for (int i = 0; i < 5; i++)
            {
                string cout = "";
                for (int j = 0; j < 5; j++)
                {
                    cout += "  " + a[i + j * 5].ToString().PadLeft(14);
                }

                Console.WriteLine(cout);
            }

            BLAS1D.dcopy(5, x, 1, ref y, 1);
            Console.WriteLine("");
            Console.WriteLine("  DCOPY ( 5, X, 1, Y, 1 )");
            Console.WriteLine("");
            for (int i = 0; i < 10; i++)
            {
                Console.WriteLine("  "
                     + (i + 1).ToString().PadLeft(6) + "  "
                     + y[i].ToString().PadLeft(14) + "");
            }

            for (int i = 0; i < 10; i++)
            {
                y[i] = (double) (10 * (i + 1));
            }

            BLAS1D.dcopy(3, x, 2, ref y, 3);
            Console.WriteLine("");
            Console.WriteLine("  DCOPY ( 3, X, 2, Y, 3 )");
            Console.WriteLine("");
            for (int i = 0; i < 10; i++)
            {
                Console.WriteLine("  "
                     + (i + 1).ToString().PadLeft(6) + "  "
                     + y[i].ToString().PadLeft(14) + "");
            }

            BLAS1D.dcopy(5, x, 1, ref a, 1);
            Console.WriteLine("");
            Console.WriteLine("  DCOPY ( 5, X, 1, A, 1 )");
            Console.WriteLine("");
            Console.WriteLine("  A =");
            Console.WriteLine("");
            for (int i = 0; i < 5; i++)
            {
                string cout = "";
                for (int j = 0; j < 5; j++)
                {
                    cout += "  " + a[i + j * 5].ToString().PadLeft(14);
                }

                Console.WriteLine(cout);
            }

            for (int i = 0; i < 5; i++)
            {
                for (int j = 0; j < 5; j++)
                {
                    a[i + j * 5] = (double) (10 * (i + 1) + j + 1);
                }
            }

            BLAS1D.dcopy(5, x, 2, ref a, 5);
            Console.WriteLine("");
            Console.WriteLine("  DCOPY ( 5, X, 2, A, 5 )");
            Console.WriteLine("");
            Console.WriteLine("  A =");
            Console.WriteLine("");
            for (int i = 0; i < 5; i++)
            {
                string cout = "";
                for (int j = 0; j < 5; j++)
                {
                    cout += "  " + a[i + j * 5].ToString().PadLeft(14);
                }

                Console.WriteLine(cout);
            }
        }

        static void ddot_test()

//****************************************************************************80
//
//  Purpose:
//
//    DDOT_TEST demonstrates DDOT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2006
//
//  Author:
//
//    John Burkardt
//
        {
            int N = 5;
            int LDA = 10;
            int LDB = 7;
            int LDC = 6;

            double[] a = new double[LDA * LDA];
            double[] b = new double[LDB * LDB];
            double[] c = new double[LDC * LDC];
            double sum1;
            double[] x = new double[N];
            double[] y = new double[N];

            Console.WriteLine("");
            Console.WriteLine("DDOT_TEST");
            Console.WriteLine("  DDOT computes the dot product of vectors.");
            Console.WriteLine("");

            for (int i = 0; i < N; i++)
            {
                x[i] = (double) (i + 1);
            }

            for (int i = 0; i < N; i++)
            {
                y[i] = -(double) (i + 1);
            }

            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    a[i + j * LDA] = (double) (i + 1 + j + 1);
                }
            }

            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    b[i + j * LDB] = (double) ((i + 1) - (j + 1));
                }
            }

            sum1 = BLAS1D.ddot(N, x, 1, y, 1 );

            Console.WriteLine("");
            Console.WriteLine("  Dot product of X and Y is " + sum1 + "");
//
//  To multiply a ROW of a matrix A times a vector X, we need to
//  specify the increment between successive entries of the row of A:
//
            sum1 = BLAS1D.ddot(N, a, LDA, x,  1, (1 + 0 * LDA), 0);

            Console.WriteLine("");
            Console.WriteLine("  Product of row 2 of A and X is " + sum1 + "");
//
//  Product of a column of A and a vector is simpler:
//
            sum1 = BLAS1D.ddot(N, a, 1, x, 1, ( + 0 + 1 * LDA), 0);

            Console.WriteLine("");
            Console.WriteLine("  Product of column 2 of A and X is " + sum1 + "");
//
//  Here's how matrix multiplication, c = a*b, could be done
//  with DDOT:
//
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    c[i + j * LDC] = BLAS1D.ddot(N, a, LDA, b, 1, (+ i), ( + 0 + j * LDB));
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  Matrix product computed with DDOT:");
            Console.WriteLine("");
            for (int i = 0; i < N; i++)
            {
                string cout = "";
                for (int j = 0; j < N; j++)
                {
                    cout += "  " + c[i + j * LDC].ToString().PadLeft(14);
                }

                Console.WriteLine(cout);
            }
        }


        static void dnrm2_test()

//****************************************************************************80
//
//  Purpose:
//
//    DNRM2_TEST demonstrates DNRM2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2006
//
//  Author:
//
//    John Burkardt
//
        {
            int N = 5;
            int LDA = 10;
//
//  These parameters illustrate the fact that matrices are typically
//  dimensioned with more space than the user requires.
//
            double[] a = new double[LDA * LDA];
            double[] x = new double[N];

            Console.WriteLine("");
            Console.WriteLine("DNRM2_TEST");
            Console.WriteLine("  DNRM2 computes the Euclidean norm of a vector.");
            Console.WriteLine("");
//
//  Compute the euclidean norm of a vector:
//
            for (int i = 0; i < N; i++)
            {
                x[i] = (double) (i + 1);
            }

            Console.WriteLine("");
            Console.WriteLine("  X =");
            Console.WriteLine("");
            for (int i = 0; i < N; i++)
            {
                Console.WriteLine("  "
                     + (i + 1).ToString().PadLeft(6) + "  "
                     + x[i].ToString().PadLeft(14) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("  The 2-norm of X is " + BLAS1D.dnrm2(N, x, 1, 0) + "");
//
//  Compute the euclidean norm of a row or column of a matrix:
//
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    a[i + j * LDA] = (double) (i + 1 + j + 1);
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The 2-norm of row 2 of A is "
                 + BLAS1D.dnrm2(N, a, LDA,( + 1 + 0 * LDA)) + "");

            Console.WriteLine("");
            Console.WriteLine("  The 2-norm of column 2 of A is "
                 + BLAS1D.dnrm2(N, a, 1, ( + 0 + 1 * LDA)) + "");

        }

        static void drot_test()

//****************************************************************************80
//
//  Purpose:
//
//    DROT_TEST tests DROT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2006
//
//  Author:
//
//    John Burkardt
//
        {
            int N = 6;

            double c;
            double s;
            double[] x = new double[N];
            double[] y = new double[N];

            for (int i = 0; i < N; i++)
            {
                x[i] = (double) (i + 1);
            }

            for (int i = 0; i < N; i++)
            {
                y[i] = (double) ((i + 1) * (i + 1) - 12);
            }

            Console.WriteLine("");
            Console.WriteLine("DROT_TEST");
            Console.WriteLine("  DROT carries out a Givens rotation.");
            Console.WriteLine("");
            Console.WriteLine("  X and Y");
            Console.WriteLine("");
            for (int i = 0; i < N; i++)
            {
                Console.WriteLine("  "
                     + (i + 1).ToString().PadLeft(6) + "  "
                     + x[i].ToString().PadLeft(14) + "  "
                     + y[i].ToString().PadLeft(14) + "");
            }

            c = 0.5;
            s = Math.Sqrt(1.0 - c * c);
            BLAS1D.drot(N, ref x,  1, ref y, 1, c, s);
            Console.WriteLine("");
            Console.WriteLine("  DROT ( N, X, 1, Y, 1, " + c + "," + s + " )");
            Console.WriteLine("");
            for (int i = 0; i < N; i++)
            {
                Console.WriteLine("  "
                     + (i + 1).ToString().PadLeft(6) + "  "
                     + x[i].ToString().PadLeft(14) + "  "
                     + y[i].ToString().PadLeft(14) + "");
            }

            for (int i = 0; i < N; i++)
            {
                x[i] = (double) (i + 1);
            }

            for (int i = 0; i < N; i++)
            {
                y[i] = (double) ((i + 1) * (i + 1) - 12);
            }

            c = x[0] / Math.Sqrt(x[0] * x[0] + y[0] * y[0]);
            s = y[0] / Math.Sqrt(x[0] * x[0] + y[0] * y[0]);
            BLAS1D.drot(N, ref x, 1, ref y, 1, c, s);
            Console.WriteLine("");
            Console.WriteLine("  DROT ( N, X, 1, Y, 1, " + c + "," + s + " )");
            Console.WriteLine("");
            for (int i = 0; i < N; i++)
            {
                Console.WriteLine("  "
                     + (i + 1).ToString().PadLeft(6) + "  "
                     + x[i].ToString().PadLeft(14) + "  "
                     + y[i].ToString().PadLeft(14) + "");
            }
        }

        static void drotg_test()

//****************************************************************************80
//
//  Purpose:
//
//    DROTG_TEST tests DROTG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2006
//
//  Author:
//
//    John Burkardt
//
        {
            double c = 0;
            double s = 0;
            int test_num = 5;

            Console.WriteLine("");
            Console.WriteLine("DROTG_TEST");
            Console.WriteLine("  DROTG generates a real Givens rotation");
            Console.WriteLine("    (  C  S ) * ( A ) = ( R )");
            Console.WriteLine("    ( -S  C )   ( B )   ( 0 )");
            Console.WriteLine("");

            int seed = 123456789;

            for (int test = 1; test <= test_num; test++)
            {
                double a = UniformRNG.r8_uniform_01(ref seed);
                double b = UniformRNG.r8_uniform_01(ref seed);

                double sa = a;
                double sb = b;

                BLAS1D.drotg(ref sa, ref sb, ref c, ref s);

                double r = sa;
                double z = sb;

                Console.WriteLine("");
                Console.WriteLine("  A =  " + a + "  B =  " + b + "");
                Console.WriteLine("  C =  " + c + "  S =  " + s + "");
                Console.WriteLine("  R =  " + r + "  Z =  " + z + "");
                Console.WriteLine("   C*A+S*B = " + (c * a + s * b) + "");
                Console.WriteLine("  -S*A+C*B = " + (-s * a + c * b) + "");
            }
        }

        static void dscal_test()

//****************************************************************************80
//
//  Purpose:
//
//    DSCAL_TEST tests DSCAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2006
//
//  Author:
//
//    John Burkardt
//
        {
            int N = 6;

            double[] x = new double[N];

            for (int i = 0; i < N; i++)
            {
                x[i] = (double) (i + 1);
            }

            Console.WriteLine("");
            Console.WriteLine("DSCAL_TEST");
            Console.WriteLine("  DSCAL multiplies a vector by a scalar.");
            Console.WriteLine("");
            Console.WriteLine("  X =");
            Console.WriteLine("");
            for (int i = 0; i < N; i++)
            {
                Console.WriteLine("  "
                     + (i + 1).ToString().PadLeft(6) + "  "
                     + x[i].ToString().PadLeft(14) + "");
            }

            double da = 5.0;
            BLAS1D.dscal(N, da, ref x,1, 0);
            Console.WriteLine("");
            Console.WriteLine("  DSCAL ( N, " + da + ", X, 1 )");
            Console.WriteLine("");
            for (int i = 0; i < N; i++)
            {
                Console.WriteLine("  "
                     + (i + 1).ToString().PadLeft(6) + "  "
                     + x[i].ToString().PadLeft(14) + "");
            }

            for (int i = 0; i < N; i++)
            {
                x[i] = (double) (i + 1);
            }

            da = -2.0;
            BLAS1D.dscal(3, da, ref x, 2, 0);
            Console.WriteLine("");
            Console.WriteLine("  DSCAL ( 3, " + da + ", X, 2 )");
            Console.WriteLine("");
            for (int i = 0; i < N; i++)
            {
                Console.WriteLine("  "
                     + (i + 1).ToString().PadLeft(6) + "  "
                     + x[i].ToString().PadLeft(14) + "");
            }
        }

        static void dswap_test()

//****************************************************************************80
//
//  Purpose:
//
//    DSWAP_TEST tests DSWAP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2006
//
//  Author:
//
//    John Burkardt
//
        {
            int N = 6;

            double[] x = new double[N];
            double[] y = new double[N];

            for (int i = 0; i < N; i++)
            {
                x[i] = (double) (i + 1);
            }

            for (int i = 0; i < N; i++)
            {
                y[i] = (double) (100 * (i + 1));
            }

            Console.WriteLine("");
            Console.WriteLine("DSWAP_TEST");
            Console.WriteLine("  DSWAP swaps two vectors.");
            Console.WriteLine("");
            Console.WriteLine("  X and Y:");
            Console.WriteLine("");
            for (int i = 0; i < N; i++)
            {
                Console.WriteLine("  "
                     + (i + 1).ToString().PadLeft(6) + "  "
                     + x[i].ToString().PadLeft(14) + "  "
                     + y[i].ToString().PadLeft(14) + "");
            }

            BLAS1D.dswap(N, ref x, 1, ref y, 1);
            Console.WriteLine("");
            Console.WriteLine("  DSWAP ( N, X, 1, Y, 1 )");
            Console.WriteLine("");
            Console.WriteLine("  X and Y:");
            Console.WriteLine("");
            for (int i = 0; i < N; i++)
            {
                Console.WriteLine("  "
                     + (i + 1).ToString().PadLeft(6) + "  "
                     + x[i].ToString().PadLeft(14) + "  "
                     + y[i].ToString().PadLeft(14) + "");
            }

            for (int i = 0; i < N; i++)
            {
                x[i] = (double) (i + 1);
            }

            for (int i = 0; i < N; i++)
            {
                y[i] = (double) (100 * (i + 1));
            }

            BLAS1D.dswap(3, ref x, 2, ref y, 1);
            Console.WriteLine("");
            Console.WriteLine("  DSWAP ( 3, X, 2, Y, 1 )");

            Console.WriteLine("");
            Console.WriteLine("  X and Y:");
            Console.WriteLine("");
            for (int i = 0; i < N; i++)
            {
                Console.WriteLine("  "
                     + (i + 1).ToString().PadLeft(6) + "  "
                     + x[i].ToString().PadLeft(14) + "  "
                     + y[i].ToString().PadLeft(14) + "");
            }
        }

        static void idamax_test()

//****************************************************************************80
//
//  Purpose:
//
//    IDAMAX_TEST demonstrates IDAMAX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2006
//
//  Author:
//
//    John Burkardt
//
        {
            int N = 11;

            double[] x = new double[N];

            Console.WriteLine("");
            Console.WriteLine("IDAMAX_TEST");
            Console.WriteLine("  IDAMAX returns the index of maximum magnitude;");

            for (int i = 1; i <= N; i++)
            {
                x[i - 1] = (double) ((7 * i) % 11) - (double) (N / 2);
            }

            Console.WriteLine("");
            Console.WriteLine("  The vector X:");
            Console.WriteLine("");

            for (int i = 1; i <= N; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(6)
                    + "  " + x[i - 1].ToString().PadLeft(8) + "");
            }

            int incx = 1;

            int i1 = BLAS1D.idamax(N, x, incx);

            Console.WriteLine("");
            Console.WriteLine("  The index of maximum magnitude = " + i1 + "");

        }
    }
}