using System;
using Burkardt.FourierTransform;

namespace FFTSerialTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for FFT_SERIAL.
            //
            //  Discussion:
            //
            //    The complex data in an N vector is stored as pairs of values in a
            //    real vector of length 2*N.
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
        {
            TimeSpan ctime;
            DateTime ctime1;
            DateTime ctime2;
            double error;
            bool first;
            double flops;
            double fnm1;
            int i;
            int icase;
            int it;
            int ln2;
            double mflops;
            int n;
            int nits = 10000;
            double seed;
            double sgn;
            double[] w;
            double[] x;
            double[] y;
            double[] z;
            double z0;
            double z1;

            Console.WriteLine("");
            Console.WriteLine("FFT_SERIAL");
            
            Console.WriteLine("");
            Console.WriteLine("  Demonstrate an implementation of the Fast Fourier Transform");
            Console.WriteLine("  of a complex data vector.");
            //
            //  Prepare for tests.
            //
            Console.WriteLine("");
            Console.WriteLine("  Accuracy check:");
            Console.WriteLine("");
            Console.WriteLine("    FFT ( FFT ( X(1:N) ) ) == N * X(1:N)");
            Console.WriteLine("");
            Console.WriteLine("             N      NITS    Error         Time          Time/Call     MFLOPS");
            Console.WriteLine("");

            seed = 331.0;
            string cout = "";
            n = 1;
            //
            //  LN2 is the log base 2 of N.  Each increase of LN2 doubles N.
            //
            for (ln2 = 1; ln2 <= 20; ln2++)
            {
                n = 2 * n;
                //
                //  Allocate storage for the complex arrays W, X, Y, Z.  
                //
                //  We handle the complex arithmetic,
                //  and store a complex number as a pair of doubles, a complex vector as a doubly
                //  dimensioned array whose second dimension is 2. 
                //
                w = new double[n];
                x = new double[2 * n];
                y = new double[2 * n];
                z = new double[2 * n];

                first = true;

                for (icase = 0; icase < 2; icase++)
                {

                    if (first)
                    {
                        for (i = 0; i < 2 * n; i = i + 2)
                        {
                            z0 = Serial.ggl(ref seed);
                            z1 = Serial.ggl(ref seed);
                            x[i] = z0;
                            z[i] = z0;
                            x[i + 1] = z1;
                            z[i + 1] = z1;
                        }
                    }
                    else
                    {
                        for (i = 0; i < 2 * n; i = i + 2)
                        {
                            z0 = 0.0;
                            z1 = 0.0;
                            x[i] = z0;
                            z[i] = z0;
                            x[i + 1] = z1;
                            z[i + 1] = z1;
                        }
                    }

                    //
                    //  Initialize the sine and cosine tables.
                    //
                    Serial.cffti(n, ref w);
                    //
                    //  Transform forward, back 
                    //
                    if (first)
                    {
                        sgn = +1.0;
                        Serial.cfft2(n, ref x, ref y, w, sgn);
                        sgn = -1.0;
                        Serial.cfft2(n, ref y, ref x, w, sgn);
                        // 
                        //  Results should be same as initial multiplied by N.
                        //
                        fnm1 = 1.0 / (double) n;
                        error = 0.0;
                        for (i = 0; i < 2 * n; i = i + 2)
                        {
                            error = error
                                    + Math.Pow(z[i] - fnm1 * x[i], 2)
                                    + Math.Pow(z[i + 1] - fnm1 * x[i + 1], 2);
                        }

                        error = Math.Sqrt(fnm1 * error);
                        cout += "  " + n.ToString().PadLeft(12)
                                               + "  " + nits.ToString().PadLeft(8)
                                               + "  " + error.ToString().PadLeft(12);
                        first = false;
                    }
                    else
                    {
                        ctime1 = DateTime.Now;
                        for (it = 0; it < nits; it++)
                        {
                            sgn = +1.0;
                            Serial.cfft2(n, ref x, ref y, w, sgn);
                            sgn = -1.0;
                            Serial.cfft2(n, ref y, ref x, w, sgn);
                        }

                        ctime2 = DateTime.Now;
                        ctime = ctime2 - ctime1;

                        flops = 2.0 * (double) nits * (5.0 * (double) n * (double) ln2);

                        mflops = flops / 1.0E+06 / ctime.TotalSeconds;

                        Console.WriteLine(cout + "  " + ctime.ToString().PadLeft(12)
                                               + "  " + (ctime / (double) (2 * nits)).ToString().PadLeft(12)
                                               + "  " + mflops.ToString().PadLeft(12) + "");
                        cout = "";
                    }
                }

                if ((ln2 % 4) == 0)
                {
                    nits = nits / 10;
                }

                if (nits < 1)
                {
                    nits = 1;
                }
            }

            Console.WriteLine("");
            Console.WriteLine("FFT_SERIAL:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }
    }
}