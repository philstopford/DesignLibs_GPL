using System;
using System.Collections.Generic;
using System.Linq;
using Clipper2Lib;

namespace geoWrangler;

public static partial class GeoWrangler
{
    // A vertex and its closest point on an opposing segment include the shortest
    // straight connection between two polygon boundaries. Obstructed candidates
    // are rejected. Kruskal minimizes total length on this candidate graph;
    // finite-width topology conflicts are rejected with alternate attachments
    // tried in length order. This does not optimize arbitrary Steiner junctions.
    private sealed record KeyholeLink(int A, int B, PathD Edge, double Length);

    private static PathsD pPlanKeyholes(PathsD input, double sizing)
    {
        ClipperD clipper = new(Constants.roundingDecimalPrecision);
        clipper.AddSubject(input);
        PolyTreeD tree = new();
        clipper.Execute(ClipType.Union, FillRule.NonZero, tree);
        PathsD result = [];
        void Visit(PolyPathD parent)
        {
            for (int i = 0; i < parent.Count; i++)
            {
                PolyPathD child = parent[i];
                if (!child.IsHole)
                {
                    PathsD component = [pStripCollinear(new PathD(child.Polygon))];
                    for (int h = 0; h < child.Count; h++)
                        component.Add(pStripCollinear(new PathD(child[h].Polygon)));
                    result.AddRange(pPlanKeyholeComponent(component, sizing));
                }
                Visit(child);
            }
        }
        Visit(tree);
        return pClose(pClockwiseAndReorderXY(result));
    }

    private static PathsD pPlanKeyholeComponent(PathsD boundaries, double sizing)
    {
        if (boundaries.Count == 1) return boundaries;
        // Canonicalize both traversal and node order for deterministic tie breaks.
        PathD outer = boundaries[0];
        var holes = boundaries.Skip(1).Select(p => pReorderXY(p)).OrderBy(p => p[0].x)
            .ThenBy(p => p[0].y).ThenBy(p => Math.Abs(Clipper.Area(p))).ToList();
        boundaries = [pReorderXY(outer)];
        boundaries.AddRange(holes);
        PathsD material = Clipper.Union(boundaries, new PathsD(), FillRule.NonZero, Constants.roundingDecimalPrecision);
        List<KeyholeLink> links = [];
        for (int a = 0; a < boundaries.Count; a++)
        for (int b = a + 1; b < boundaries.Count; b++)
        {
            List<KeyholeLink> candidates = [];
            void Project(PathD vertices, PathD segments)
            {
                foreach (PointD point in vertices)
                for (int j = 0; j < segments.Count; j++)
                {
                    PointD start = segments[j], end = segments[(j + 1) % segments.Count];
                    double dx = end.x - start.x, dy = end.y - start.y;
                    double square = dx * dx + dy * dy;
                    if (square <= 1E-16) continue;
                    double t = Math.Clamp(((point.x - start.x) * dx + (point.y - start.y) * dy) / square, 0, 1);
                    PointD target = new(start.x + t * dx, start.y + t * dy);
                    double length = distanceBetweenPoints(point, target);
                    if (length <= Constants.tolerance) continue;
                    PathD edge = point.x < target.x || (point.x == target.x && point.y <= target.y)
                        ? [point, target] : [target, point];
                    candidates.Add(new(a, b, edge, length));
                }
            }
            Project(boundaries[a], boundaries[b]);
            Project(boundaries[b], boundaries[a]);
            // Retain alternate attachments for finite-width conflicts. Test
            // them lazily in Kruskal order, skipping already connected nodes.
            links.AddRange(candidates.DistinctBy(c =>
                (c.Edge[0].x, c.Edge[0].y, c.Edge[1].x, c.Edge[1].y)));
        }

        int[] parents = Enumerable.Range(0, boundaries.Count).ToArray();
        int Root(int n) { while (parents[n] != n) n = parents[n]; return n; }
        PathsD corridors = [];
        int accepted = 0;
        foreach (KeyholeLink link in links.OrderBy(c => c.Length).ThenBy(c => c.Edge[0].x)
                     .ThenBy(c => c.Edge[0].y).ThenBy(c => c.Edge[1].x).ThenBy(c => c.Edge[1].y))
        {
            int a = Root(link.A), b = Root(link.B);
            if (a == b) continue;
            ClipperD lineClipper = new(Constants.roundingDecimalPrecision);
            lineClipper.AddOpenSubject(link.Edge);
            lineClipper.AddClip(material);
            PathsD inside = [], closed = [];
            lineClipper.Execute(ClipType.Intersection, FillRule.NonZero, closed, inside);
            double insideLength = inside.Sum(p => Enumerable.Range(1, p.Count - 1)
                .Sum(i => distanceBetweenPoints(p[i - 1], p[i])));
            if (Math.Abs(insideLength - link.Length) > 0.0005) continue;
            PathsD proposed = new(corridors);
            proposed.AddRange(pInflateEdge(link.Edge, sizing));
            // Union the corridors, rather than XOR: overlapping caps must not
            // leave a material island at a junction.
            proposed = Clipper.Union(proposed, new PathsD(), FillRule.NonZero, Constants.roundingDecimalPrecision);
            PathsD cut = Clipper.Difference(material, proposed, FillRule.NonZero, Constants.roundingDecimalPrecision);
            if (cut.Count(p => Clipper.IsPositive(p)) != 1 ||
                cut.Count(p => !Clipper.IsPositive(p)) != boundaries.Count - accepted - 2) continue;
            parents[a] = b;
            accepted++;
            corridors = proposed;
            if (accepted == boundaries.Count - 1) return cut;
        }
        // Never silently drop an unresolved physical hole or a disconnected
        // material island. The caller needs an explicit failure for this case.
        throw new InvalidOperationException("Could not construct a topology-preserving keyhole tree at the requested corridor width.");
    }
}
