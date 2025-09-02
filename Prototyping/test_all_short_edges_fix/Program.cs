using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

using PathD = System.Collections.Generic.List<PointD>;
using PathsD = System.Collections.Generic.List<System.Collections.Generic.List<PointD>>;

/// <summary>
/// Test and fix for the all short edges scenario where the algorithm fails.
/// This demonstrates the issue and implements the solution to connect midpoints.
/// </summary>
public struct PointD
{
    public double x, y;
    
    public PointD(double x_, double y_) { x = x_; y = y_; }
    
    public static PointD operator +(PointD a, PointD b) => new PointD(a.x + b.x, a.y + b.y);
    public static PointD operator -(PointD a, PointD b) => new PointD(a.x - b.x, a.y - b.y);
    public static PointD operator *(PointD a, double s) => new PointD(a.x * s, a.y * s);
    
    public double Length() => Math.Sqrt(x * x + y * y);
    
    public PointD Normalized()
    {
        double len = Length();
        return len > 0 ? new PointD(x / len, y / len) : new PointD(0, 0);
    }
}

enum CornerType { Convex = 0, Concave = 1, ShortEdge = 2 }

class Program
{
    static void Main()
    {
        Console.WriteLine("Testing all short edges scenario...");
        
        // Create a square where all edges will be considered short
        List<PointD> square = new List<PointD>
        {
            new PointD(0, 0),
            new PointD(10, 0),     // Short horizontal edge (10 units)
            new PointD(10, 10),    // Short vertical edge (10 units)  
            new PointD(0, 10),     // Short horizontal edge (10 units)
            new PointD(0, 0)       // Short vertical edge (10 units) - close the path
        };

        double short_edge_threshold = 15.0; // All edges (10 units) are below this threshold
        
        Console.WriteLine("Original square vertices:");
        for (int i = 0; i < square.Count; i++)
        {
            Console.WriteLine($"  {i}: ({square[i].x}, {square[i].y})");
        }

        // Classify corners
        int[] corner_types = CategorizeCorners(square, short_edge_threshold);
        
        Console.WriteLine("\nCorner classifications:");
        for (int i = 0; i < corner_types.Length; i++)
        {
            string type = corner_types[i] == (int)CornerType.ShortEdge ? "ShortEdge" : 
                         corner_types[i] == (int)CornerType.Convex ? "Convex" : "Concave";
            Console.WriteLine($"  Corner {i}: {type}");
        }

        // Create dummy processed corners for testing
        PathsD processedCorners = new PathsD();
        for (int i = 0; i < square.Count - 1; i++)
        {
            PathD corner = new PathD();
            corner.Add(square[i]);
            corner.Add(square[i + 1]);
            processedCorners.Add(corner);
        }

        Console.WriteLine("\nTesting original AssembleWithDiagonals function...");
        try
        {
            // This should cause an infinite loop in the original implementation
            PathD result = AssembleWithDiagonals_Original(processedCorners, corner_types);
            Console.WriteLine("Original function completed (unexpected!)");
        }
        catch (Exception ex)
        {
            Console.WriteLine($"Original function failed as expected: {ex.Message}");
        }

        Console.WriteLine("\nTesting fixed AssembleWithDiagonals function...");
        PathD fixedResult = AssembleWithDiagonals_Fixed(processedCorners, corner_types, square);
        
        Console.WriteLine("Fixed function result (diamond from square midpoints):");
        for (int i = 0; i < fixedResult.Count; i++)
        {
            Console.WriteLine($"  {i}: ({fixedResult[i].x}, {fixedResult[i].y})");
        }

        // Expected diamond points: (5,0), (10,5), (5,10), (0,5)
        var expectedDiamond = new List<PointD>
        {
            new PointD(5, 0),   // Midpoint of bottom edge
            new PointD(10, 5),  // Midpoint of right edge  
            new PointD(5, 10),  // Midpoint of top edge
            new PointD(0, 5)    // Midpoint of left edge
        };

        Console.WriteLine("\nExpected diamond vertices:");
        for (int i = 0; i < expectedDiamond.Count; i++)
        {
            Console.WriteLine($"  {i}: ({expectedDiamond[i].x}, {expectedDiamond[i].y})");
        }

        // Create SVG visualization
        string svg = CreateSvg(square, fixedResult, expectedDiamond);
        File.WriteAllText("all_short_edges_test.svg", svg, Encoding.UTF8);
        Console.WriteLine("\nSVG visualization saved to all_short_edges_test.svg");
    }

    /// <summary>
    /// Original implementation that will fail when all edges are short
    /// </summary>
    static PathD AssembleWithDiagonals_Original(PathsD processedCorners, int[] corner_types)
    {
        int n = processedCorners.Count;
        PathD outPts = new PathD();

        // Mark short corners
        bool[] isShort = new bool[n];
        for (int k = 0; k < n; k++)
            isShort[k] = (k < corner_types.Length && corner_types[k] == (int)CornerType.ShortEdge);

        // Check if all corners are short - this is the problematic case
        bool allShort = isShort.All(x => x);
        if (allShort)
        {
            throw new InvalidOperationException("All edges are short - original algorithm will loop infinitely!");
        }

        int i = 0;
        while (i < n)
        {
            if (!isShort[i])
            {
                // append corner polyline, avoiding duplicating a point
                PathD poly = processedCorners[i];
                if (poly.Count > 0)
                {
                    if (outPts.Count == 0)
                        outPts.AddRange(poly);
                    else
                    {
                        PointD last = outPts[^1];
                        PointD first = poly[0];
                        if (Math.Abs(last.x - first.x) < 1e-9 && Math.Abs(last.y - first.y) < 1e-9)
                            outPts.AddRange(poly.Skip(1));
                        else
                            outPts.AddRange(poly);
                    }
                }
                i++;
                continue;
            }

            // Start of short run
            int runStart = i;
            while (i < n && isShort[i]) i++;
            int runEnd = i - 1;

            // Find previous non-short index (wrap) - INFINITE LOOP when all are short!
            int prevIdx = (runStart - 1 + n) % n;
            while (isShort[prevIdx]) prevIdx = (prevIdx - 1 + n) % n;

            // Find next non-short index (wrap) - INFINITE LOOP when all are short!
            int nextIdx = (runEnd + 1) % n;
            while (isShort[nextIdx]) nextIdx = (nextIdx + 1) % n;

            PointD startPt = processedCorners[prevIdx].Last();
            PointD endPt = processedCorners[nextIdx].First();

            // Add diagonal
            if (outPts.Count == 0)
            {
                outPts.Add(startPt);
            }
            else
            {
                PointD last = outPts[^1];
                if (!(Math.Abs(last.x - startPt.x) < 1e-9 && Math.Abs(last.y - startPt.y) < 1e-9))
                    outPts.Add(startPt);
            }

            outPts.Add(endPt);
        }

        return outPts;
    }

    /// <summary>
    /// Fixed implementation that handles the all-short-edges case
    /// </summary>
    static PathD AssembleWithDiagonals_Fixed(PathsD processedCorners, int[] corner_types, List<PointD> originalPath)
    {
        int n = processedCorners.Count;
        PathD outPts = new PathD();

        // Mark short corners
        bool[] isShort = new bool[n];
        for (int k = 0; k < n; k++)
            isShort[k] = (k < corner_types.Length && corner_types[k] == (int)CornerType.ShortEdge);

        // Special case: all edges are short - connect midpoints
        bool allShort = isShort.All(x => x);
        if (allShort)
        {
            Console.WriteLine("All edges are short - connecting midpoints...");
            return ConnectMidpoints(originalPath);
        }

        // Original logic for mixed edge cases
        int i = 0;
        while (i < n)
        {
            if (!isShort[i])
            {
                // append corner polyline, avoiding duplicating a point
                PathD poly = processedCorners[i];
                if (poly.Count > 0)
                {
                    if (outPts.Count == 0)
                        outPts.AddRange(poly);
                    else
                    {
                        PointD last = outPts[^1];
                        PointD first = poly[0];
                        if (Math.Abs(last.x - first.x) < 1e-9 && Math.Abs(last.y - first.y) < 1e-9)
                            outPts.AddRange(poly.Skip(1));
                        else
                            outPts.AddRange(poly);
                    }
                }
                i++;
                continue;
            }

            // Start of short run
            int runStart = i;
            while (i < n && isShort[i]) i++;
            int runEnd = i - 1;

            // Find previous non-short index (wrap)
            int prevIdx = (runStart - 1 + n) % n;
            while (isShort[prevIdx]) prevIdx = (prevIdx - 1 + n) % n;

            // Find next non-short index (wrap)
            int nextIdx = (runEnd + 1) % n;
            while (isShort[nextIdx]) nextIdx = (nextIdx + 1) % n;

            PointD startPt = processedCorners[prevIdx].Last();
            PointD endPt = processedCorners[nextIdx].First();

            // Add diagonal
            if (outPts.Count == 0)
            {
                outPts.Add(startPt);
            }
            else
            {
                PointD last = outPts[^1];
                if (!(Math.Abs(last.x - startPt.x) < 1e-9 && Math.Abs(last.y - startPt.y) < 1e-9))
                    outPts.Add(startPt);
            }

            outPts.Add(endPt);
        }

        return outPts;
    }

    /// <summary>
    /// Connect the midpoints of edges to form the output polygon.
    /// For a square, this creates a diamond.
    /// </summary>
    static PathD ConnectMidpoints(List<PointD> originalPath)
    {
        PathD midpoints = new PathD();
        
        // Calculate midpoint of each edge (excluding the last point if it's a duplicate of the first)
        int edgeCount = originalPath.Count - 1;
        if (originalPath.Count > 1 && 
            Math.Abs(originalPath[0].x - originalPath[^1].x) < 1e-9 && 
            Math.Abs(originalPath[0].y - originalPath[^1].y) < 1e-9)
        {
            // Path is closed, don't count the duplicate last point
            edgeCount = originalPath.Count - 1;
        }
        else
        {
            edgeCount = originalPath.Count - 1;
        }

        for (int i = 0; i < edgeCount; i++)
        {
            PointD p1 = originalPath[i];
            PointD p2 = originalPath[i + 1];
            PointD midpoint = new PointD((p1.x + p2.x) / 2.0, (p1.y + p2.y) / 2.0);
            midpoints.Add(midpoint);
        }

        return midpoints;
    }

    /// <summary>
    /// Classify corners based on edge lengths
    /// </summary>
    static int[] CategorizeCorners(List<PointD> path_, double short_edge_length)
    {
        int[] status = new int[path_.Count - 1]; // Don't include duplicate last point

        List<PointD> path = new List<PointD>(path_);
        bool trimmed_path = false;
        if (path.Count > 1 && path[0].x == path[^1].x && path[0].y == path[^1].y)
        {
            trimmed_path = true;
            path.RemoveAt(path.Count - 1);
        }

        double area2 = 0;
        for (int i = 0; i < path.Count; i++)
        {
            PointD p1 = path[i];
            PointD p2 = path[(i + 1) % path.Count];
            area2 += p1.x * p2.y - p2.x * p1.y;
        }

        bool isCCW = area2 > 0;

        for (int i = 0; i < path.Count; i++)
        {
            PointD prev = path[(i - 1 + path.Count) % path.Count];
            PointD curr = path[i];
            PointD next = path[(i + 1) % path.Count];

            double vx1 = curr.x - prev.x;
            double vy1 = curr.y - prev.y;
            double vx2 = next.x - curr.x;
            double vy2 = next.y - curr.y;

            double len1 = Math.Sqrt(vx1 * vx1 + vy1 * vy1);
            double len2 = Math.Sqrt(vx2 * vx2 + vy2 * vy2);

            if (len1 <= short_edge_length && len2 <= short_edge_length)
            {
                status[i] = (int)CornerType.ShortEdge;
                continue;
            }

            double crossZ = vx1 * vy2 - vy1 * vx2;
            bool isVertexConvex = isCCW ? (crossZ > 0) : (crossZ < 0);
            status[i] = (int)(isVertexConvex ? CornerType.Convex : CornerType.Concave);
        }

        return status;
    }

    /// <summary>
    /// Create SVG visualization showing original square, result, and expected diamond
    /// </summary>
    static string CreateSvg(List<PointD> original, PathD result, List<PointD> expected)
    {
        var allPoints = original.Concat(result).Concat(expected).ToList();
        double minX = allPoints.Min(p => p.x) - 5;
        double maxX = allPoints.Max(p => p.x) + 5;
        double minY = allPoints.Min(p => p.y) - 5;
        double maxY = allPoints.Max(p => p.y) + 5;

        StringBuilder sb = new StringBuilder();
        sb.AppendLine($"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"400\" height=\"400\" viewBox=\"{minX} {-maxY} {maxX - minX} {maxY - minY}\">");
        
        // Grid
        sb.AppendLine("  <defs>");
        sb.AppendLine("    <pattern id=\"grid\" width=\"5\" height=\"5\" patternUnits=\"userSpaceOnUse\">");
        sb.AppendLine("      <path d=\"M 5 0 L 0 0 0 5\" fill=\"none\" stroke=\"#e0e0e0\" stroke-width=\"0.5\"/>");
        sb.AppendLine("    </pattern>");
        sb.AppendLine("  </defs>");
        sb.AppendLine($"  <rect width=\"100%\" height=\"100%\" fill=\"url(#grid)\" />");

        // Original square (blue)
        sb.Append("  <polyline fill=\"none\" stroke=\"blue\" stroke-width=\"2\" points=\"");
        for (int i = 0; i < original.Count; i++)
        {
            sb.Append($"{original[i].x},{-original[i].y} ");
        }
        sb.AppendLine("\"/>");

        // Result (red)
        sb.Append("  <polyline fill=\"none\" stroke=\"red\" stroke-width=\"3\" points=\"");
        foreach (var p in result)
        {
            sb.Append($"{p.x},{-p.y} ");
        }
        // Close the path
        if (result.Count > 0)
        {
            sb.Append($"{result[0].x},{-result[0].y} ");
        }
        sb.AppendLine("\"/>");

        // Expected diamond (green, dashed)
        sb.Append("  <polyline fill=\"none\" stroke=\"green\" stroke-width=\"2\" stroke-dasharray=\"3,3\" points=\"");
        foreach (var p in expected)
        {
            sb.Append($"{p.x},{-p.y} ");
        }
        // Close the path
        if (expected.Count > 0)
        {
            sb.Append($"{expected[0].x},{-expected[0].y} ");
        }
        sb.AppendLine("\"/>");

        // Legend
        sb.AppendLine("  <text x=\"2\" y=\"-2\" font-size=\"12\" fill=\"blue\">Original Square</text>");
        sb.AppendLine("  <text x=\"2\" y=\"15\" font-size=\"12\" fill=\"red\">Result (Fixed)</text>");
        sb.AppendLine("  <text x=\"2\" y=\"32\" font-size=\"12\" fill=\"green\">Expected Diamond</text>");

        sb.AppendLine("</svg>");
        return sb.ToString();
    }
}