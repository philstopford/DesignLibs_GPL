using System;
using Burkardt.Uniform;

namespace Burkardt.Probability;

public static class Buffon
{
    public static double buffon_box_pdf(double a, double b, double l)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BUFFON_BOX_PDF evaluates the Buffon Box PDF.
        //
        //  Discussion:
        //
        //    In the Buffon-Laplace needle experiment, we suppose that the plane has
        //    been tiled into a grid of rectangles of width A and height B, and that a
        //    needle of length L is dropped "at random" onto this grid.
        //
        //    We may assume that one end, the "eye" of the needle falls at the point
        //    (X1,Y1), taken uniformly at random in the cell [0,A]x[0,B].
        //
        //    ANGLE, the angle that the needle makes is taken to be uniformly random.
        //    The point of the needle, (X2,Y2), therefore lies at
        //
        //      (X2,Y2) = ( X1+L*cos(ANGLE), Y1+L*sin(ANGLE) )
        //
        //    The needle will have crossed at least one grid line if any of the
        //    following are true:
        //
        //      X2 <= 0, A <= X2, Y2 <= 0, B <= Y2.
        //
        //    If L is larger than sqrt ( A*A + B*B ), then the needle will
        //    cross every time, and the computation is uninteresting.  However, if
        //    L is smaller than this limit, then the probability of a crossing on
        //    a single trial is
        //
        //      P(L,A,B) = ( 2 * L * ( A + B ) - L * L ) / ( PI * A * B )
        //
        //    and therefore, a record of the number of hits for a given number of
        //    trials can be used as a very roundabout way of estimating PI.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Sudarshan Raghunathan,
        //    Making a Supercomputer Do What You Want: High Level Tools for
        //    Parallel Programming,
        //    Computing in Science and Engineering,
        //    Volume 8, Number 5, September/October 2006, pages 70-80.
        //
        //  Parameters:
        //
        //    Input, double A, B, the horizontal and vertical dimensions
        //    of each cell of the grid.  0 <= A, 0 <= B.
        //
        //    Input, double L, the length of the needle.
        //    0 <= L <= min ( A, B ).
        //
        //    Output, double BUFFON_BOX_PDF, the PDF.
        //
    {
        double pdf;

        switch (a)
        {
            case < 0.0:
                Console.WriteLine("");
                Console.WriteLine("BUFFON_BOX_PDF - Fatal error!");
                Console.WriteLine("  Input A < 0.");
                return 1;
            case 0.0:
                pdf = 1.0;
                return pdf;
        }

        switch (b)
        {
            case < 0.0:
                Console.WriteLine("");
                Console.WriteLine("BUFFON_BOX_PDF - Fatal error!");
                Console.WriteLine("  Input B < 0.");
                return 1;
            case 0.0:
                pdf = 1.0;
                return pdf;
        }

        switch (l)
        {
            case < 0.0:
                Console.WriteLine("");
                Console.WriteLine("BUFFON_BOX_PDF - Fatal error!");
                Console.WriteLine("  Input L < 0.");
                return 1;
            case 0.0:
                pdf = 0.0;
                return pdf;
            default:
            {
                if (Math.Min(a, b) < l)
                {
                    Console.WriteLine("");
                    Console.WriteLine("BUFFON_BOX_PDF - Fatal error!");
                    Console.WriteLine("  min ( A, B ) < L.");
                    return 1;
                }

                break;
            }
        }

        pdf = (2.0 * l * (a + b) - l * l) / (Math.PI * a * b);

        return pdf;
    }

    public static int buffon_box_sample(double a, double b, double l, int trial_num)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BUFFON_BOX_SAMPLE samples the Buffon Box distribution.
        //
        //  Discussion:
        //
        //    In the Buffon-Laplace needle experiment, we suppose that the plane has
        //    been tiled into a grid of rectangles of width A and height B, and that a
        //    needle of length L is dropped "at random" onto this grid.
        //
        //    We may assume that one end, the "eye" of the needle falls at the point
        //    (X1,Y1), taken uniformly at random in the cell [0,A]x[0,B].
        //
        //    ANGLE, the angle that the needle makes is taken to be uniformly random.
        //    The point of the needle, (X2,Y2), therefore lies at
        //
        //      (X2,Y2) = ( X1+L*cos(ANGLE), Y1+L*sin(ANGLE) )
        //
        //    The needle will have crossed at least one grid line if any of the
        //    following are true:
        //
        //      X2 <= 0, A <= X2, Y2 <= 0, B <= Y2.
        //
        //    This routine simulates the tossing of the needle, and returns the number
        //    of times that the needle crossed at least one grid line.
        //
        //    If L is larger than sqrt ( A*A + B*B ), then the needle will
        //    cross every time, and the computation is uninteresting.  However, if
        //    L is smaller than this limit, then the probability of a crossing on
        //    a single trial is
        //
        //      P(L,A,B) = ( 2 * L * ( A + B ) - L * L ) / ( PI * A * B )
        //
        //    and therefore, a record of the number of hits for a given number of
        //    trials can be used as a very roundabout way of estimating PI.
        //    (Particularly roundabout, since we actually will use a good value of
        //    PI in order to pick the random angles.)
        //
        //    Since this routine invokes the C++ random number generator,
        //    the user should initialize the random number generator, particularly
        //    if it is desired to control whether the sequence is to be varied
        //    or repeated.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Sudarshan Raghunathan,
        //    Making a Supercomputer Do What You Want: High Level Tools for
        //    Parallel Programming,
        //    Computing in Science and Engineering,
        //    Volume 8, Number 5, September/October 2006, pages 70-80.
        //
        //  Parameters:
        //
        //    Input, double A, B, the horizontal and vertical dimensions
        //    of each cell of the grid.  0 <= A, 0 <= B.
        //
        //    Input, double L, the length of the needle.
        //    0 <= L <= min ( A, B ).
        //
        //    Input, int TRIAL_NUM, the number of times the needle is
        //    to be dropped onto the grid.
        //
        //    Output, int BUFFON_BOX_SAMPLE, the number of times the needle
        //    crossed at least one line of the grid of cells.
        //
    {

        int hits = 0;

        int seed = 12345678;

        for (int trial = 1; trial <= trial_num; trial++)
        {
            //
            //  Randomly choose the location of the eye of the needle in [0,0]x[A,B],
            //  and the angle the needle makes.
            //
            double x1 = a * UniformRNG.r8_uniform_01(ref seed);
            double y1 = b * UniformRNG.r8_uniform_01(ref seed);
            double angle = 2.0 * Math.PI * UniformRNG.r8_uniform_01(ref seed);
            //
            //  Compute the location of the point of the needle.
            //
            double x2 = x1 + l * Math.Cos(angle);
            double y2 = y1 + l * Math.Sin(angle);
            //
            //  Count the end locations that lie outside the cell.
            //
            if (x2 <= 0.0 || a <= x2 || y2 <= 0.0 || b <= y2)
            {
                hits += 1;
            }
        }

        return hits;
    }

    public static double buffon_pdf(double a, double l)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BUFFON_PDF evaluates the Buffon PDF.
        //
        //  Discussion:
        //
        //    In the Buffon needle experiment, we suppose that the plane has been
        //    ruled by vertical lines with a spacing of A units, and that a
        //    needle of length L is dropped "at random" onto this grid.
        //
        //    Because of the various symmetries, we may assume that this eye of
        //    this needle lands in the first infinite strip, and we may further
        //    assume that its Y coordinate is 0.  Thus, we have
        //    the eye as (X1,Y1) with 0 <= X1 <= A and Y1 = 0.
        //
        //    ANGLE, the angle that the needle makes is taken to be uniformly random.
        //    The point of the needle, (X2,Y2), therefore lies at
        //
        //      (X2,Y2) = ( X1+L*cos(ANGLE), Y1+L*sin(ANGLE) )
        //
        //    The needle will have crossed at least one grid line if any of the
        //    following are true:
        //
        //      X2 <= 0, A <= X2.
        //
        //    The probability of a crossing on a single trial is
        //
        //      P(A,L) = ( 2 * L ) / ( PI * A )
        //
        //    and therefore, a record of the number of hits for a given number of
        //    trials can be used as a very roundabout way of estimating PI.
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
        //  Parameters:
        //
        //    Input, double A, the horizontal spacing between the
        //    vertical grid lines.  0 <= A.
        //
        //    Input, double L, the length of the needle.
        //
        //    Output, double BUFFON_PDF, the Buffon PDF.
        //
    {
        double pdf;

        switch (a)
        {
            case < 0.0:
                Console.WriteLine("");
                Console.WriteLine("BUFFON_PDF - Fatal error!");
                Console.WriteLine("  Input A < 0.");
                return 1.0;
            case 0.0:
                pdf = 1.0;
                return pdf;
        }

        switch (l)
        {
            case < 0.0:
                Console.WriteLine("");
                Console.WriteLine("BUFFON_PDF - Fatal error!");
                Console.WriteLine("  Input L < 0.");
                return 1.0;
            case 0.0:
                pdf = 0.0;
                return pdf;
            default:
                pdf = 2.0 * l / (Math.PI * a);

                return pdf;
        }
    }

    public static int buffon_sample(double a, double l, int trial_num)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BUFFON_SAMPLE samples the Buffon PDF.
        //
        //  Discussion:
        //
        //    In the Buffon needle experiment, we suppose that the plane has been
        //    ruled by vertical lines with a spacing of A units, and that a
        //    needle of length L is dropped "at random" onto this grid.
        //
        //    Because of the various symmetries, we may assume that this eye of
        //    this needle lands in the first infinite strip, and we may further
        //    assume that its Y coordinate is 0.  Thus, we have
        //    the eye as (X1,Y1) with 0 <= X1 <= A and Y1 = 0.
        //
        //    ANGLE, the angle that the needle makes is taken to be uniformly random.
        //    The point of the needle, (X2,Y2), therefore lies at
        //
        //      (X2,Y2) = ( X1+L*cos(ANGLE), Y1+L*sin(ANGLE) )
        //
        //    The needle will have crossed at least one grid line if any of the
        //    following are true:
        //
        //      X2 <= 0, A <= X2.
        //
        //    The probability of a crossing on a single trial is
        //
        //      P(A,L) = ( 2 * L ) / ( PI * A )
        //
        //    and therefore, a record of the number of hits for a given number of
        //    trials can be used as a very roundabout way of estimating PI.
        //
        //    Since this routine invokes the C++ random number generator,
        //    the user should initialize the random number generator, particularly
        //    if it is desired to control whether the sequence is to be varied
        //    or repeated.
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
        //  Parameters:
        //
        //    Input, double A, the horizontal spacing between the
        //    vertical grid lines.  0 <= A.
        //
        //    Input, double L, the length of the needle.
        //
        //    Input, int TRIAL_NUM, the number of times the needle is
        //    to be dropped onto the grid.
        //
        //    Output, int BUFFON_SAMPLE, the number of times the needle
        //    crossed at least one line of the grid of cells.
        //
    {

        int hits = 0;

        int seed = 12345678;
            
        for (int trial = 1; trial <= trial_num; trial++)
        {
            //
            //  Randomly choose the location (X1,Y1) of the eye of the needle
            //  in [0,0]x[A,0], and the angle the needle makes.
            //
            double x1 = a * UniformRNG.r8_uniform_01(ref seed); // (double) rand() / (double) RAND_MAX;
            double angle = 2.0 * Math.PI * UniformRNG.r8_uniform_01(ref seed); // (double) rand() / (double) RAND_MAX;
            //
            //  Compute the location of the point of the needle.
            //  We only need to know the value of X2, not Y2!
            //
            double x2 = x1 + l * Math.Cos(angle);
            //
            //  Count the end locations that lie outside the cell.
            //
            if (x2 <= 0.0 || a <= x2)
            {
                hits += 1;
            }
        }

        return hits;
    }
}