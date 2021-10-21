using System;

namespace Burkardt.FourierTransform
{
    public static class Serial
    {
        public static void ccopy(int n, double[] x, ref double[] y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CCOPY copies a complex vector.
            //
            //  Discussion:
            //
            //    The "complex" vector A[N] is actually stored as a double vector B[2*N].
            //
            //    The "complex" vector entry A[I] is stored as:
            //
            //      B[I*2+0], the real part,
            //      B[I*2+1], the imaginary part.
            //
            //  Modified:
            //
            //    23 March 2009
            //
            //  Author:
            //
            //    Original C version by Wesley Petersen.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Wesley Petersen, Peter Arbenz, 
            //    Introduction to Parallel Computing - A practical guide with examples in C,
            //    Oxford University Press,
            //    ISBN: 0-19-851576-6,
            //    LC: QA76.58.P47.
            //
            //  Parameters:
            //
            //    Input, int N, the length of the "complex" array.
            //
            //    Input, double X[2*N], the array to be copied.
            //
            //    Output, double Y[2*N], a copy of X.
            //
        {
            int i;

            for (i = 0; i < n; i++)
            {
                y[i * 2 + 0] = x[i * 2 + 0];
                y[i * 2 + 1] = x[i * 2 + 1];
            }
        }

        public static void cfft2(int n, ref double[] x, ref double[] y, double[] w, double sgn)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CFFT2 performs a complex Fast Fourier Transform.
            //
            //  Modified:
            //
            //    23 March 2009
            //
            //  Author:
            //
            //    Original C version by Wesley Petersen.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Wesley Petersen, Peter Arbenz, 
            //    Introduction to Parallel Computing - A practical guide with examples in C,
            //    Oxford University Press,
            //    ISBN: 0-19-851576-6,
            //    LC: QA76.58.P47.
            //
            //  Parameters:
            //
            //    Input, int N, the size of the array to be transformed.
            //
            //    Input/output, double X[2*N], the data to be transformed.  
            //    On output, the contents of X have been overwritten by work information.
            //
            //    Output, double Y[2*N], the forward or backward FFT of X.
            //
            //    Input, double W[N], a table of sines and cosines.
            //
            //    Input, double SGN, is +1 for a "forward" FFT and -1 for a "backward" FFT.
            //
        {
            int j;
            int m;
            int mj;
            bool tgle = false;

            m = (int) (Math.Log((double) n) / Math.Log(1.99));
            mj = 1;
            //
            //  Toggling switch for work array.
            //
            tgle = true;
            step(n, mj, x, 0 * 2 + 0, x, (n / 2) * 2 + 0, ref y, 0 * 2 + 0, ref y, mj * 2 + 0, w, 0, sgn);

            if (n == 2)
            {
                return;
            }

            for (j = 0; j < m - 2; j++)
            {
                mj = mj * 2;
                if (tgle)
                {
                    step(n, mj, y, 0 * 2 + 0, y, (n / 2) * 2 + 0, ref x, 0 * 2 + 0, ref x, mj * 2 + 0, w, 0, sgn);
                    tgle = false;
                }
                else
                {
                    step(n, mj, x, 0 * 2 + 0, x, (n / 2) * 2 + 0, ref y, 0 * 2 + 0, ref y, mj * 2 + 0, w, 0, sgn);
                    tgle = true;
                }
            }

            //
            //  Last pass thru data: move y to x if needed 
            //
            if (tgle)
            {
                ccopy(n, y, ref x);
            }

            mj = n / 2;
            step(n, mj, x, 0 * 2 + 0, x, (n / 2) * 2 + 0, ref y, 0 * 2 + 0, ref y, mj * 2 + 0, w, 0, sgn);

        }

        public static void cffti(int n, ref double[] w)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CFFTI sets up sine and cosine tables needed for the FFT calculation.
            //
            //  Modified:
            //
            //    23 March 2009
            //
            //  Author:
            //
            //    Original C version by Wesley Petersen.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Wesley Petersen, Peter Arbenz, 
            //    Introduction to Parallel Computing - A practical guide with examples in C,
            //    Oxford University Press,
            //    ISBN: 0-19-851576-6,
            //    LC: QA76.58.P47.
            //
            //  Parameters:
            //
            //    Input, int N, the size of the array to be transformed.
            //
            //    Output, double W[N], a table of sines and cosines.
            //
        {
            double arg;
            double aw;
            int i;
            int n2;
            const double pi = 3.141592653589793;

            n2 = n / 2;
            aw = 2.0 * pi / ((double) n);

            for (i = 0; i < n2; i++)
            {
                arg = aw * ((double) i);
                w[i * 2 + 0] = Math.Cos(arg);
                w[i * 2 + 1] = Math.Sin(arg);
            }

            return;
        }

        public static double ggl(ref double seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GGL generates uniformly distributed pseudorandom numbers. 
            //
            //  Modified:
            //
            //    23 March 2009
            //
            //  Author:
            //
            //    Original C version by Wesley Petersen, M Troyer, I Vattulainen.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Wesley Petersen, Peter Arbenz, 
            //    Introduction to Parallel Computing - A practical guide with examples in C,
            //    Oxford University Press,
            //    ISBN: 0-19-851576-6,
            //    LC: QA76.58.P47.
            //
            //  Parameters:
            //
            //    Input/output, double *SEED, used as a seed for the sequence.
            //
            //    Output, double GGL, the next pseudorandom value.
            //
        {
            double d2 = 0.2147483647e10;
            double t;
            double value;

            t = seed;
            t = (16807.0 * t) % d2;
            seed = t;
            value = (t - 1.0) / (d2 - 1.0);

            return value;
        }

        public static void step(int n, int mj, double[] a, int aIndex, double[] b, int bIndex, ref double[] c,
                int cIndex,
                ref double[] d, int dIndex, double[] w, int wIndex, double sgn)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    STEP carries out one step of the workspace version of CFFT2.
            //
            //  Modified:
            //
            //    23 March 2009
            //
            //  Author:
            //
            //    Original C version by Wesley Petersen.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Wesley Petersen, Peter Arbenz, 
            //    Introduction to Parallel Computing - A practical guide with examples in C,
            //    Oxford University Press,
            //    ISBN: 0-19-851576-6,
            //    LC: QA76.58.P47.
            //
            //  Parameters:
            //
        {
            double ambr;
            double ambu;
            int j;
            int ja;
            int jb;
            int jc;
            int jd;
            int jw;
            int k;
            int lj;
            int mj2;
            double[] wjw = new double[2];

            mj2 = 2 * mj;
            lj = n / mj2;

            for (j = 0; j < lj; j++)
            {
                jw = j * mj;
                ja = jw;
                jb = ja;
                jc = j * mj2;
                jd = jc;

                wjw[0] = w[(wIndex + (jw * 2 + 0)) % w.Length];
                wjw[1] = w[(wIndex + (jw * 2 + 1)) % w.Length];

                if (sgn < 0.0)
                {
                    wjw[1] = -wjw[1];
                }

                for (k = 0; k < mj; k++)
                {
                    c[(cIndex + ((jc + k) * 2 + 0)) % c.Length] = a[(aIndex + ((ja + k) * 2 + 0)) % a.Length] +
                                                                  b[(bIndex + ((jb + k) * 2 + 0)) % b.Length];
                    c[(cIndex + ((jc + k) * 2 + 1)) % c.Length] = a[(aIndex + ((ja + k) * 2 + 1)) % a.Length] +
                                                                  b[(bIndex + ((jb + k) * 2 + 1)) % b.Length];

                    ambr = a[(aIndex + ((ja + k) * 2 + 0)) % a.Length] - b[(bIndex + ((jb + k) * 2 + 0)) % b.Length];
                    ambu = a[(aIndex + ((ja + k) * 2 + 1)) % a.Length] - b[(bIndex + ((jb + k) * 2 + 1)) % b.Length];

                    d[(dIndex + ((jd + k) * 2 + 0)) % d.Length] = wjw[0] * ambr - wjw[1] * ambu;
                    d[(dIndex + ((jd + k) * 2 + 1)) % d.Length] = wjw[1] * ambr + wjw[0] * ambu;
                }
            }
        }
    }
}