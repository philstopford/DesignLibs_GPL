using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using Clipper2Lib;

namespace geoCoreLib;

public class Repetition
{
    public enum RepetitionType
    {
        None,
        Rectangular,
        Regular,
        ExplicitX,
        ExplicitY,
        Explicit
    }

    public RepetitionType type;

    public int columns; // along x or v1
    public int rows; // along y or v2

    public Point64 spacing; // rectangular spacing
    public Point64 rowVector; // regular axis 1
    public Point64 colVector; // regular axis 2

    public Path64 offsets; // explicit
    public List<double> coords; // explicitx and explicity

    public Repetition()
    {
        reset();
    }

    public void reset()
    {
        type = RepetitionType.None;
        columns = 0;
        rows = 0;
        spacing = new(0,0);
        rowVector = new(0,0);
        colVector = new(0,0);
        offsets = new();
        coords = new();
    }
    
    public Repetition(Repetition source)
    {
        copy_from(source);
    }

    public void copy_from(Repetition source)
    {
        type = source.type;
        columns = source.columns;
        rows = source.rows;
        spacing = new(source.spacing);
        rowVector = new(source.rowVector);
        colVector = new(source.colVector);
        offsets = new(source.offsets);
        coords = new(source.coords);
    }

    public void clear()
    {
        if (type == RepetitionType.Explicit)
        {
            offsets.Clear();
        }

        if ((type == RepetitionType.ExplicitX) || (type == RepetitionType.ExplicitY))
        {
            coords.Clear();
        }
    }


    // Not sure whether these methods are needed, but added for now ahead of investigations.
    public int get_count()
    {
        switch (type)
        {
            case RepetitionType.Rectangular:
            case RepetitionType.Regular:
                return columns * rows;
            case RepetitionType.Explicit:
                return offsets.Count + 1; // Assume (0, 0) is not included.
            case RepetitionType.ExplicitX:
            case RepetitionType.ExplicitY:
                return coords.Count + 1; // Assume 0 is not included.
            case RepetitionType.None:
                return 0;
        }

        return 0;
    }

    public Path64 get_offsets()
    {
        Path64 result = new();
        int count = get_count();
        switch (type)
        {
            case RepetitionType.Rectangular:
                for (int i = 0; i < columns; i++)
                {
                    Int64 cx = i * spacing.X;
                    for (int j = 0; j < rows; j++)
                    {
                        result.Append(new(cx, j * spacing.Y));
                    }
                }

                break;
            case RepetitionType.Regular:
                for (int i = 0; i < columns; i++)
                {
                    Point64 vi = new(rowVector.X * i, rowVector.Y * i);
                    for (int j = 0; j < rows; j++)
                    {
                        result.Append(new(vi.X + j * colVector.X, vi.Y + j * colVector.Y));
                    }
                }

                break;
            case RepetitionType.ExplicitX:
                result.Append(new(0, 0));
                for (int j = 1; j < count; j++)
                {
                    result.Append(new(coords[j], 0));
                }

                break;
            case RepetitionType.ExplicitY:
                result.Append(new(0, 0));
                for (int j = 1; j < count; j++)
                {
                    result.Append(new(0, coords[j]));
                }

                break;
            case RepetitionType.Explicit:
                result.Append(new(0, 0));
                result.AddRange(offsets);
                break;
            case RepetitionType.None:
                break;
        }

        return result;
    }

    public Path64 get_extrema()
    {
        Path64 result = new();

        switch (type)
        {
            case RepetitionType.Rectangular:
                if (columns == 0 || rows == 0)
                {
                    return result;
                }

                if (columns == 1)
                {
                    if (rows == 1)
                    {
                        result.Append(new (0, 0));
                    }
                    else
                    {
                        result.Append(new (0, 0));
                        result.Append(new (0, (rows - 1) * spacing.Y));
                    }
                }
                else
                {
                    if (rows == 1)
                    {
                        result.Append(new (0, 0));
                        result.Append(new ((columns - 1) * spacing.X, 0));
                    }
                    else
                    {
                        result.Append(new (0, 0));
                        result.Append(new (0, (rows - 1) * spacing.Y));
                        result.Append(new ((columns - 1) * spacing.X, 0));
                        result.Append(new ((columns - 1) * spacing.X, (rows - 1) * spacing.Y));
                    }
                }

                break;
            case RepetitionType.Regular:
                if (columns == 0 || rows == 0)
                {
                    return result;
                }

                if (columns == 1)
                {
                    if (rows == 1)
                    {
                        result.Append(new (0, 0));
                    }
                    else
                    {
                        result.Append(new (0, 0));
                        result.Append(new((rows - 1) * colVector.X, (rows - 1) * colVector.Y));
                    }
                }
                else
                {
                    if (rows == 1)
                    {
                        result.Append(new (0, 0));
                        result.Append(new((columns - 1) * rowVector.X, (columns - 1) * rowVector.Y));
                    }
                    else
                    {
                        PointD vi = new((columns - 1) * rowVector.X, (columns - 1) * rowVector.Y);
                        PointD vj = new((rows - 1) * colVector.X, (rows - 1) * colVector.Y);
                        result.Append(new (0, 0));
                        result.Append(new (vi));
                        result.Append(new (vj));
                        result.Append(new (vi.x + vj.x, vi.y + vj.y));
                    }
                }

                break;
            case RepetitionType.ExplicitX:
                if (coords.Count == 0)
                {
                    return result;
                }

                double xmin = coords.Min();
                double xmax = coords.Max();
                result.Append(new(xmin, 0));
                if (xmin != xmax)
                {
                    result.Append(new(xmax, 0));
                }

                break;
            case RepetitionType.ExplicitY:
                if (coords.Count == 0)
                {
                    return result;
                }

                double ymin = coords.Min();
                double ymax = coords.Max();
                result.Append(new(0, ymin));
                if (ymin != ymax)
                {
                    result.Append(new(0, ymax));
                }

                break;
            case RepetitionType.Explicit:
                if (coords.Count == 0)
                {
                    return result;
                }

                Point64 vxmin = new(offsets[0]);
                Point64 vxmax = new(offsets[0]);
                Point64 vymin = new(offsets[0]);
                Point64 vymax = new(offsets[0]);
                foreach (Point64 v in offsets)
                {
                    if (v.X < vxmin.X)
                    {
                        vxmin = new(v);
                    }
                    else if (v.X > vxmax.X)
                    {
                        vxmax = new(v);
                    }

                    if (v.Y < vymin.Y)
                    {
                        vymin = new(v);
                    }
                    else if (v.Y > vymax.Y)
                    {
                        vymax = new(v);
                    }
                }

                result.Append(vxmin);
                result.Append(vxmax);
                result.Append(vymin);
                result.Append(vymax);
                break;
            case RepetitionType.None:
                return result;
        }

        return result;
    }

    public List<GCPolygon> transform(List<GCPolygon> geo, double magnification, bool x_reflection, double rotation)
    {
        bool doRotationInFunction = true;
        if (type == RepetitionType.None)
        {
            return geo.ToList();
        }
        
        switch (type)
        {
            case RepetitionType.Rectangular:
                if (magnification != 1)
                {
                    spacing = new (spacing.X * magnification, spacing.Y * magnification);
                }
                if (x_reflection || (rotation != 0 && doRotationInFunction))
                {
                    Point64 v = new(spacing);
                    if (x_reflection) v.Y = -v.Y;
                    double ca = Math.Cos(rotation);
                    double sa = Math.Sin(rotation);
                    type = RepetitionType.Regular;
                    rowVector.X = (Int64)(v.Y * ca);
                    rowVector.Y = (Int64)(v.Y * sa);
                    colVector.X = (Int64)(-v.Y * sa);
                    colVector.Y = (Int64)(v.Y * ca);
                }
                break;
            case RepetitionType.Regular:
                /*
                if (magnification != 1)
                {
                    rowVector = new (rowVector.X * magnification, rowVector.Y * magnification);
                    rowVector = new (colVector.X * magnification, colVector.Y * magnification);
                }

                if (x_reflection)
                {
                    rowVector.Y = -rowVector.Y;
                    colVector.Y = -colVector.Y;
                }

                if (rotation != 0  && doRotationInFunction)
                {
                    Complex z2x = Complex.Cos(rotation);
                    Complex z2y = Complex.Sin(rotation);

                    // Based on
                    / *
                        Vec2 r = {cos(rotation), sin(rotation)};
                        v1 = cplx_mul(v1, r);
                        v2 = cplx_mul(v2, r);

                        inline Vec2 cplx_mul(const Vec2& z1, const Vec2& z2) {
                           return Vec2{z1.re * z2.re - z1.im * z2.im, z1.re * z2.im + z1.im * z2.re};
                        }
                     * /
                    rowVector = new (rowVector.X * z2x.Real - rowVector.X * z2x.Imaginary, rowVector.Y * z2y.Real + rowVector.Y * z2y.Imaginary);
                    colVector = new (colVector.X * z2x.Real - colVector.X * z2x.Imaginary, colVector.Y * z2y.Real + colVector.Y * z2y.Imaginary);

                    //colVector = GeoWrangler.Rotate(new Point64(0.0, 0.0), colVector, rotation+180);
                    //rowVector = GeoWrangler.Rotate(new Point64(0.0, 0.0), rowVector, rotation+180);
                }
                */
                break;
            case RepetitionType.ExplicitX:
                if (rotation != 0  && doRotationInFunction)
                {
                    double ca = magnification * Math.Cos(rotation);
                    double sa = magnification * Math.Sin(rotation);
                    Path64 temp = new ();
                    foreach (double c in coords)
                    {
                        temp.Add(new(c * ca, c * sa));
                    }
                    type = RepetitionType.Explicit;
                    offsets = new(temp);
                }
                else if (magnification != 1)
                {
                    for (int i = 0; i < coords.Count; i++)
                    {
                        coords[i] *= magnification;
                    }
                }
                break;
            case RepetitionType.ExplicitY:
                if (rotation != 0 && doRotationInFunction)
                {
                    double ca = magnification * Math.Cos(rotation);
                    double sa = -magnification * Math.Sin(rotation);
                    if (x_reflection)
                    {
                        ca = -ca;
                        sa = -sa;
                    }

                    Path64 temp = new ();
                    foreach (double c in coords)
                    {
                        temp.Add(new(c * sa, c * ca));
                    }
                    type = RepetitionType.Explicit;
                    type = RepetitionType.Explicit;
                    offsets = new(temp);
                }
                else if (x_reflection || magnification != 1)
                {
                    if (x_reflection)
                    {
                        magnification = -magnification;
                    }
                    for (int i = 0; i < coords.Count; i++)
                    {
                        coords[i] *= magnification;
                    }
                }
                break;
            case RepetitionType.Explicit:
            {
                if (rotation != 0 && doRotationInFunction)
                {
                    Complex z2x = Complex.Cos(rotation);
                    Complex z2y = Complex.Sin(rotation);

                    // Based on 
                    /*
                        Vec2 r = {cos(rotation), sin(rotation)};
                        v1 = cplx_mul(v1, r);
                        v2 = cplx_mul(v2, r);

                        inline Vec2 cplx_mul(const Vec2& z1, const Vec2& z2) {
                           return Vec2{z1.re * z2.re - z1.im * z2.im, z1.re * z2.im + z1.im * z2.re};
                        }
                     */
                    // v1 = new PointD(v1.x * z2x.Real - v1.x * z2x.Imaginary, v1.y * z2y.Imaginary + v1.y * z2y.Imaginary);
                    // v2 = new PointD(v2.x * z2x.Real - v2.x * z2x.Imaginary, v2.y * z2y.Imaginary + v2.y * z2y.Imaginary);

                    if (x_reflection)
                    {
                        z2y = Complex.Conjugate(z2y);
                        for (int i = 0; i < offsets.Count; i++)
                        {
                            offsets[i] = new (offsets[i].X * z2x.Real - offsets[i].X * z2x.Imaginary, offsets[i].Y * z2y.Imaginary + offsets[i].Y * z2y.Imaginary);
                        }
                    }
                    else
                    {
                        for (int i = 0; i < offsets.Count; i++)
                        {
                            offsets[i] = new (offsets[i].X * z2x.Real - offsets[i].X * z2x.Imaginary, offsets[i].Y * z2y.Imaginary + offsets[i].Y * z2y.Imaginary);
                        }
                    }
                }
                else if (x_reflection && magnification != 1)
                {
                    for (int i = 0; i < offsets.Count; i++)
                    {
                        offsets[i] = new(offsets[i].X * magnification, -offsets[i].Y * magnification);
                    }
                }
                else if (x_reflection)
                {
                    for (int i = 0; i < offsets.Count; i++)
                    {
                        offsets[i] = new (offsets[i].X, -offsets[i].Y);
                    }
                }
                else if (magnification != 1)
                {
                    for (int i = 0; i < offsets.Count; i++)
                    {
                        offsets[i] = new (offsets[i].X * magnification, offsets[i].Y * magnification);
                    }
                }
            }
                break;
            default:
                return geo.ToList();
        }
        
        // Scale and rotate.
        List<GCPolygon> scaled_and_rotated = new();
        foreach (GCPolygon poly in geo)
        {
            GCPolygon temp = new(poly);
            temp.rotate(rotation, new(0,0));
            temp.scale(magnification);
            scaled_and_rotated.Add(temp);
        }

        List<GCPolygon> ret = new();

        if (offsets.Count > 0)
        {
            foreach (Point64 offset in offsets)
            {
                foreach (GCPolygon poly in scaled_and_rotated)
                {
                    GCPolygon temp = new(poly);
                    temp.move(offset);
                    ret.Add(temp);
                }
            }
        }

        if (coords.Count > 0)
        {
            for (int idx = 0; idx < coords.Count; idx += 2)
            {
                foreach (GCPolygon poly in scaled_and_rotated)
                {
                    GCPolygon temp = new(poly);
                    temp.move(new Point64(coords[idx], coords[idx+1]));
                    ret.Add(temp);
                }
            }
        }

        for (int y = 0; y < rows; y++)
        {
            for (int x = 0; x < columns; x++)
            {
                foreach (GCPolygon poly in scaled_and_rotated)
                {
                    GCPolygon temp = new(poly);
                    temp.move(new (x * (rowVector.X + colVector.X), y * (rowVector.Y + colVector.Y)));
                    ret.Add(temp);
                }
            }
        }
        
        if (x_reflection)
        {
            foreach (GCPolygon poly in ret)
            {
                Path64 flipped = new();
                foreach (Point64 pt in poly.pointarray)
                {
                    flipped.Add(new (pt.X, -pt.Y));
                }
                poly.pointarray.Clear();
                poly.pointarray.AddRange(flipped);
            }
        }
        
        /*
#if !GCSINGLETHREADED
        Parallel.For(0, ret.Count, (poly, loopstate) =>
#else
            for (int poly = 0; poly < ret.Count; poly++)
#endif
            {
                ret[poly].rotate(trans.angle, new(0,0));
                ret[poly].move(point);
            }
#if !GCSINGLETHREADED
        );
#endif
        */

        return ret;
    }
}