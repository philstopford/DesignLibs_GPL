using System;
using Burkardt;
using Burkardt.IO;

namespace MandelbrotPMMASCIITest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for MANDELBROT.
        //
        //  Discussion:
        //
        //    MANDELBROT computes an image of the Mandelbrot set.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Local Parameters:
        //
        //    Local, integer COUNT_MAX, the maximum number of iterations taken
        //    for a particular pixel.
        //
    {
        int[] b;
        int c;
        int c_max;
        int[] count;
        int count_max = 400;
        string filename = "mandelbrot.ppm";
        int[] g;
        int i;
        int j;
        int k;
        int n = 501;
        int[] r;
        double x;
        double x_max = 1.25;
        double x_min = -2.25;
        double x1;
        double x2;
        double y;
        double y_max = 1.75;
        double y_min = -1.75;
        double y1;
        double y2;

        Console.WriteLine("");
        Console.WriteLine("MANDELBROT");
        Console.WriteLine("");
        Console.WriteLine("  Create an ASCII PPM image of the Mandelbrot set.");
        Console.WriteLine("");
        Console.WriteLine("  For each point C = X + i*Y");
        Console.WriteLine("  with X range [" + x_min + "," + x_max + "]");
        Console.WriteLine("  and  Y range [" + y_min + "," + y_max + "]");
        Console.WriteLine("  carry out " + count_max + " iterations of the map");
        Console.WriteLine("  Z(n+1) = Z(n)^2 + C.");
        Console.WriteLine("  If the iterates stay bounded (norm less than 2)");
        Console.WriteLine("  then C is taken to be a member of the set.");
        Console.WriteLine("");
        Console.WriteLine("  An ASCII PPM image of the set is created using");
        Console.WriteLine("    N = " + n + " pixels in the X direction and");
        Console.WriteLine("    N = " + n + " pixels in the Y direction.");
        //
        //  Carry out the iteration for each pixel, determining COUNT.
        //
        count = new int[n * n];

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                x = (j * x_max
                     + (n - j - 1) * x_min)
                    / (n - 1);

                y = (i * y_max
                     + (n - i - 1) * y_min)
                    / (n - 1);

                count[i + j * n] = 0;

                x1 = x;
                y1 = y;

                for (k = 1; k <= count_max; k++)
                {
                    x2 = x1 * x1 - y1 * y1 + x;
                    y2 = 2 * x1 * y1 + y;

                    if (x2 < -2.0 || 2.0 < x2 || y2 < -2.0 || 2.0 < y2)
                    {
                        count[i + j * n] = k;
                        break;
                    }

                    x1 = x2;
                    y1 = y2;
                }
            }
        }

        //
        //  Determine the coloring of each pixel.
        //
        c_max = 0;
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                if (c_max < count[i + j * n])
                {
                    c_max = count[i + j * n];
                }
            }
        }

        //
        //  Set the image data.
        //
        r = new int[n * n];
        g = new int[n * n];
        b = new int[n * n];
        int index = 0;

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                index = i + j * n;
                switch (count[i + j * n] % 2)
                {
                    case 1:
                        r[index] = 255;
                        g[index] = 255;
                        b[index] = 255;
                        break;
                    default:
                    {
                        c = (int) (255.0 * Math.Sqrt(Math.Sqrt(Math.Sqrt(
                            count[i + j * n] / (double) c_max))));
                        int v = 3 * c / 5;
                        r[index] = v;
                        g[index] = v;
                        b[index] = c;
                        break;
                    }
                }
            }
        }

        //
        //  Write an image file.
        //
        PPM_ASCII.ppma_write(filename, n, n, r, g, b);

        Console.WriteLine("");
        Console.WriteLine("  ASCII PPM image data stored in \"" + filename + "\".");

        Console.WriteLine("");
        Console.WriteLine("MANDELBROT");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}