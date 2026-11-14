// Port of the spatial indexes in PolyClipping/Clipper2 (Boost Software License 1.0).
#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;

namespace Clipper2Lib;

/// <summary>A snapshot of polygon edges, partitioned by Y for repeated point queries.</summary>
public sealed class StripeIndex
{
  private (Point64 a, Point64 b)[] _edges = [];
  private List<int>[] _stripes = [];
  private long _yMin, _height = 1;
  public bool Empty => _edges.Length == 0;

  private int YToStripe(long y) => (int)Math.Clamp((y - _yMin) / _height, 0, _stripes.Length - 1);

  public void BuildFromPath(Path64 path)
  {
    var edges = new List<(Point64, Point64)>();
    if (path.Count >= 3)
      for (int i = 0, j = path.Count - 1; i < path.Count; j = i++)
        if (path[j] != path[i]) edges.Add((path[j], path[i]));
    Build(edges);
  }

  internal void BuildFromRing(OutPt op)
  {
    var edges = new List<(Point64, Point64)>();
    if (op.next != op && op.prev != op.next)
    {
      OutPt cur = op;
      do
      {
        if (cur.pt != cur.next!.pt) edges.Add((cur.pt, cur.next!.pt));
        cur = cur.next!;
      } while (cur != op);
    }
    Build(edges);
  }

  private void Build(List<(Point64 a, Point64 b)> edges)
  {
    _edges = [];
    _stripes = [];
    if (edges.Count == 0) return;
    _yMin = edges.Min(e => Math.Min(e.a.Y, e.b.Y));
    long yMax = edges.Max(e => Math.Max(e.a.Y, e.b.Y));
    // Zero-height rings use the brute-force result, which is always outside.
    if (yMax == _yMin) return;
    int count = Math.Clamp((int)Math.Sqrt(edges.Count), 1, 65536);
    _height = (yMax - _yMin + count - 1) / count;
    _stripes = new List<int>[count];
    long refs = 0;
    foreach (var e in edges)
      refs += YToStripe(Math.Max(e.a.Y, e.b.Y)) - YToStripe(Math.Min(e.a.Y, e.b.Y)) + 1;
    if (refs > (long)edges.Count * 8) { _stripes = []; return; }
    for (int s = 0; s < count; s++) _stripes[s] = [];
    _edges = edges.ToArray();
    for (int i = 0; i < edges.Count; i++)
    {
      var e = edges[i];
      int max = YToStripe(Math.Max(e.a.Y, e.b.Y));
      for (int s = YToStripe(Math.Min(e.a.Y, e.b.Y)); s <= max; s++) _stripes[s].Add(i);
    }
  }

  public PointInPolygonResult PointInPolygon(Point64 pt)
  {
    if (Empty) return PointInPolygonResult.IsOutside;
    bool inside = false;
    foreach (int i in _stripes[YToStripe(pt.Y)])
    {
      var (a, b) = _edges[i];
      int cross = InternalClipper.CrossProductSign(a, b, pt);
      if (cross == 0 && pt.X >= Math.Min(a.X, b.X) && pt.X <= Math.Max(a.X, b.X) &&
        pt.Y >= Math.Min(a.Y, b.Y) && pt.Y <= Math.Max(a.Y, b.Y)) return PointInPolygonResult.IsOn;
      if ((a.Y > pt.Y) != (b.Y > pt.Y) && (b.Y > pt.Y ? cross > 0 : cross < 0)) inside = !inside;
    }
    return inside ? PointInPolygonResult.IsInside : PointInPolygonResult.IsOutside;
  }
}

// Packed R-tree used by offset orientation and containment queries.
internal sealed class BoundsTree
{
  private sealed class Node(Rect64 bounds, int first, int count)
  {
    internal Rect64 Bounds = bounds;
    internal int First = first, Count = count, ActiveCount;
  }
  private sealed class Level
  {
    internal int[] Items = [], Parents = [];
    internal List<Node> Nodes = [];
  }
  private readonly Rect64[] _boxes;
  private readonly bool[] _active;
  private readonly List<Level> _levels = [];
  internal BoundsTree(Rect64[] boxes)
  {
    _boxes = boxes;
    _active = new bool[boxes.Length];
    Rect64[] rects = boxes;
    while (rects.Length != 0)
    {
      var level = new Level { Items = TileOrder(rects), Parents = new int[rects.Length] };
      for (int first = 0; first < rects.Length; first += 16)
      {
        int count = Math.Min(16, rects.Length - first);
        var r = new Rect64(long.MaxValue, long.MaxValue, long.MinValue, long.MinValue);
        for (int i = first; i < first + count; i++)
        {
          Rect64 b = rects[level.Items[i]];
          r.left = Math.Min(r.left, b.left); r.top = Math.Min(r.top, b.top);
          r.right = Math.Max(r.right, b.right); r.bottom = Math.Max(r.bottom, b.bottom);
          level.Parents[level.Items[i]] = level.Nodes.Count;
        }
        level.Nodes.Add(new Node(r, first, count));
      }
      _levels.Add(level);
      if (level.Nodes.Count == 1) break;
      rects = level.Nodes.Select(n => n.Bounds).ToArray();
    }
  }
  private static int[] TileOrder(Rect64[] rects)
  {
    int[] order = Enumerable.Range(0, rects.Length).ToArray();
    Array.Sort(order, (a, b) => (rects[a].left / 2.0 + rects[a].right / 2.0)
      .CompareTo(rects[b].left / 2.0 + rects[b].right / 2.0));
    int tiles = (rects.Length + 15) / 16;
    int slices = (int)Math.Ceiling(Math.Sqrt(tiles));
    int len = ((tiles + slices - 1) / slices) * 16;
    for (int first = 0; first < order.Length; first += len)
      Array.Sort(order, first, Math.Min(len, order.Length - first), Comparer<int>.Create((a, b) =>
        (rects[a].top / 2.0 + rects[a].bottom / 2.0).CompareTo(rects[b].top / 2.0 + rects[b].bottom / 2.0)));
    return order;
  }
  internal void SetActive(int box, bool active)
  {
    if (_active[box] == active) return;
    _active[box] = active;
    int item = box;
    foreach (Level level in _levels)
    {
      item = level.Parents[item];
      level.Nodes[item].ActiveCount += active ? 1 : -1;
    }
  }
  internal void Query(Func<Rect64, bool> meets, Action<int> visit, bool activeOnly = false)
  {
    if (_levels.Count == 0) return;
    var stack = new Stack<(int level, int node)>();
    stack.Push((_levels.Count - 1, 0));
    while (stack.Count != 0)
    {
      var (idx, n) = stack.Pop();
      Level level = _levels[idx];
      Node node = level.Nodes[n];
      if ((activeOnly && node.ActiveCount == 0) || !meets(node.Bounds)) continue;
      for (int i = node.First; i < node.First + node.Count; i++)
      {
        int item = level.Items[i];
        if (idx > 0) stack.Push((idx - 1, item));
        else if ((!activeOnly || _active[item]) && meets(_boxes[item])) visit(item);
      }
    }
  }
}

internal sealed class PolygonIndex
{
  private readonly Path64 _path;
  private readonly BoundsTree _edges;
  internal PolygonIndex(Path64 path)
  {
    _path = path;
    var bounds = new Rect64[path.Count];
    for (int i = 0; i < path.Count; i++)
    {
      Point64 a = path[i], b = path[(i + 1) % path.Count];
      bounds[i] = new Rect64(Math.Min(a.X, b.X), Math.Min(a.Y, b.Y), Math.Max(a.X, b.X), Math.Max(a.Y, b.Y));
    }
    _edges = new BoundsTree(bounds);
  }
  internal PointInPolygonResult Locate(Point64 pt)
  {
    bool on = false, inside = false;
    _edges.Query(r => !on && r.left <= pt.X && r.right >= pt.X && r.top <= pt.Y && r.bottom >= pt.Y,
      i => on = InternalClipper.CrossProductSign(_path[i], _path[(i + 1) % _path.Count], pt) == 0);
    if (on) return PointInPolygonResult.IsOn;
    _edges.Query(r => r.left <= pt.X && r.top <= pt.Y && r.bottom > pt.Y, i =>
    {
      Point64 a = _path[i], b = _path[(i + 1) % _path.Count];
      int d = InternalClipper.CrossProductSign(a, b, pt);
      if ((a.Y > pt.Y) != (b.Y > pt.Y) && (a.Y < b.Y ? d < 0 : d > 0)) inside = !inside;
    });
    return inside ? PointInPolygonResult.IsInside : PointInPolygonResult.IsOutside;
  }
}
