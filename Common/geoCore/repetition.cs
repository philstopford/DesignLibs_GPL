using System;
using System.Collections.Generic;
using System.Linq;
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
        spacing = new Point64(0,0);
        rowVector = new Point64(0,0);
        colVector = new Point64(0,0);
        offsets = [];
        coords = [];
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
        spacing = new Point64(source.spacing);
        rowVector = new Point64(source.rowVector);
        colVector = new Point64(source.colVector);
        offsets = new Path64(source.offsets);
        coords = [..source.coords];
    }

    public void clear()
    {
        switch (type)
        {
            case RepetitionType.Explicit:
                offsets.Clear();
                break;
            case RepetitionType.ExplicitX or RepetitionType.ExplicitY:
                coords.Clear();
                break;
        }
    }


    // Not sure whether these methods are needed, but added for now ahead of investigations.
    public int get_count()
    {
        return type switch
        {
            RepetitionType.Rectangular or RepetitionType.Regular => columns * rows,
            RepetitionType.Explicit => offsets.Count // Original is not included.
            ,
            RepetitionType.ExplicitX or RepetitionType.ExplicitY => coords.Count // Original is not included.
            ,
            RepetitionType.None => 0,
            _ => 0
        };
    }

    public Path64 get_offsets()
    {
        Path64 result = [];
        int count = get_count();
        switch (type)
        {
            case RepetitionType.Rectangular:
                for (int i = 0; i < rows; i++)
                {
                    long cy = i * spacing.Y;
                    for (int j = 0; j < columns; j++)
                    {
                        if ((i == 0) && (j == 0))
                        {
                            continue;
                        }
                        result.Add(new Point64(j * spacing.X, cy));
                    }
                }

                break;
            case RepetitionType.Regular:
                for (int i = 0; i < rows; i++)
                {
                    Point64 vi = new(rowVector.X * i, rowVector.Y * i);
                    for (int j = 1; j < columns; j++)
                    {
                        result.Add(new Point64(vi.X + j * colVector.X, vi.Y + j * colVector.Y));
                    }
                }

                break;
            case RepetitionType.ExplicitX:
                // result.Add(new(0, 0));
                for (int j = 1; j < count; j++)
                {
                    result.Add(new Point64(coords[j], 0));
                }

                break;
            case RepetitionType.ExplicitY:
                // result.Add(new(0, 0));
                for (int j = 1; j < count; j++)
                {
                    result.Add(new Point64(0, coords[j]));
                }

                break;
            case RepetitionType.Explicit:
                result.Add(new Point64(offsets[0]));
                // The explicit offset is made up of delta values, each one a delta from the previous value....
                for (int i = 1; i < offsets.Count; i++)
                {
                    result.Add(new Point64(offsets[i].X + result[^1].X, offsets[i].Y + result[^1].Y));
                }
                break;
            case RepetitionType.None:
                break;
        }

        return result;
    }

    public Path64 get_extrema()
    {
        Path64 result = [];

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
                        result.Add(new Point64(0, 0));
                    }
                    else
                    {
                        result.Add(new Point64(0, 0));
                        result.Add(new Point64(0, (rows - 1) * spacing.Y));
                    }
                }
                else
                {
                    if (rows == 1)
                    {
                        result.Add(new Point64(0, 0));
                        result.Add(new Point64((columns - 1) * spacing.X, 0));
                    }
                    else
                    {
                        result.Add(new Point64(0, 0));
                        result.Add(new Point64(0, (rows - 1) * spacing.Y));
                        result.Add(new Point64((columns - 1) * spacing.X, 0));
                        result.Add(new Point64((columns - 1) * spacing.X, (rows - 1) * spacing.Y));
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
                        result.Add(new Point64(0, 0));
                    }
                    else
                    {
                        result.Add(new Point64(0, 0));
                        result.Add(new Point64((rows - 1) * colVector.X, (rows - 1) * colVector.Y));
                    }
                }
                else
                {
                    if (rows == 1)
                    {
                        result.Add(new Point64(0, 0));
                        result.Add(new Point64((columns - 1) * rowVector.X, (columns - 1) * rowVector.Y));
                    }
                    else
                    {
                        PointD vi = new((columns - 1) * rowVector.X, (columns - 1) * rowVector.Y);
                        PointD vj = new((rows - 1) * colVector.X, (rows - 1) * colVector.Y);
                        result.Add(new Point64(0, 0));
                        result.Add(new Point64(vi));
                        result.Add(new Point64(vj));
                        result.Add(new Point64(vi.x + vj.x, vi.y + vj.y));
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
                result.Add(new Point64(xmin, 0));
                if (Math.Abs(xmin - xmax) > GeoCore.tolerance)
                {
                    result.Add(new Point64(xmax, 0));
                }

                break;
            case RepetitionType.ExplicitY:
                if (coords.Count == 0)
                {
                    return result;
                }

                double ymin = coords.Min();
                double ymax = coords.Max();
                result.Add(new Point64(0, ymin));
                if (Math.Abs(ymin - ymax) > GeoCore.tolerance)
                {
                    result.Add(new Point64(0, ymax));
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
                        vxmin = new Point64(v);
                    }
                    else if (v.X > vxmax.X)
                    {
                        vxmax = new Point64(v);
                    }

                    if (v.Y < vymin.Y)
                    {
                        vymin = new Point64(v);
                    }
                    else if (v.Y > vymax.Y)
                    {
                        vymax = new Point64(v);
                    }
                }

                result.Add(vxmin);
                result.Add(vxmax);
                result.Add(vymin);
                result.Add(vymax);
                break;
            case RepetitionType.None:
                return result;
        }

        return result;
    }

    public List<GCPolygon> transform(List<GCPolygon> geo, double magnification, bool x_reflection, double rotation)
    {
        const bool doRotationInFunction = true;
        bool bypass_regular_array = false; // used for specific cases.
        if (type == RepetitionType.None)
        {
            return geo.ToList();
        }
        
        switch (type)
        {
            case RepetitionType.Rectangular:
                if (Math.Abs(magnification - 1) > GeoCore.tolerance)
                {
                    spacing = new Point64(spacing.X * magnification, spacing.Y * magnification);
                }
                if (x_reflection || (rotation != 0 && doRotationInFunction))
                {
                    Point64 v = new(spacing);
                    if (x_reflection)
                    {
                        v.Y = -v.Y;
                    }
                    double ca = Math.Cos(rotation);
                    double sa = Math.Sin(rotation);
                    type = RepetitionType.Regular;
                    rowVector.X = (long)(v.Y * ca);
                    rowVector.Y = (long)(v.Y * sa);
                    colVector.X = (long)(-v.Y * sa);
                    colVector.Y = (long)(v.Y * ca);
                }
                else
                {
                    rowVector = new Point64(spacing.X, 0);
                    colVector = new Point64(0, spacing.Y);
                }
                break;
            case RepetitionType.Regular:
                if (columns == 1)
                {
                    if (offsets.Count == 0)
                    {
                        offsets.Add(new Point64(0, 0));
                        for (int r = 1; r < rows; r++)
                        {
                            offsets.Add(rowVector);
                        }
                    }

                    bypass_regular_array = true;
                }
                else
                {
                    if (rows == 1)
                    {
                        if (offsets.Count == 0)
                        {
                            offsets.Add(new Point64(0, 0));
                            for (int c = 1; c < columns; c++)
                            {
                                offsets.Add(colVector);
                            }
                        }

                        bypass_regular_array = true;
                    }
                }
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
                    Path64 temp = [];
                    temp.AddRange(coords.Select(c => new Point64(c * ca, c * sa)));
                    type = RepetitionType.Explicit;
                    offsets = new Path64(temp);
                }
                else if (Math.Abs(magnification - 1) > GeoCore.tolerance)
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

                    Path64 temp = [];
                    temp.AddRange(coords.Select(c => new Point64(c * sa, c * ca)));
                    type = RepetitionType.Explicit;
                    offsets = new Path64(temp);
                }
                else if (x_reflection || Math.Abs(magnification - 1) > GeoCore.tolerance)
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
                bypass_regular_array = true;
                // Process offsets
                Point64 _colVector = new(colVector.X, Math.Ceiling((double)colVector.Y / columns));
                Point64 _rowVector = new(Math.Ceiling((double)rowVector.X / rows), rowVector.Y);
                Path64 temp = [];
                for (int c = 1; c < columns; c++)
                {
                    temp.Add(new Point64(_colVector.X * c, _colVector.Y));
                }
                for (int r = 0; r < rows; r++)
                {
                    for (int t = 0; t < temp.Count; t++)
                    {
                        offsets.Add(new Point64(temp[t].X + (_rowVector.X * r), temp[t].Y + (_rowVector.Y * r)));
                    }
                }
                offsets.Add(new Point64(_rowVector));
                /*
                if (rotation != 0 && doRotationInFunction)
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
                */
            }
                break;
            default:
                return geo.ToList();
        }
        
        // Scale and rotate.
        List<GCPolygon> scaled_and_rotated = [];
        foreach (var temp in geo.Select(poly => new GCPolygon(poly)))
        {
            if (x_reflection)
            {
                for (int pt = 0; pt < temp.pointarray.Count; pt++)
                {
                    temp.pointarray[pt] = new Point64(temp.pointarray[pt].X, -temp.pointarray[pt].Y);
                }
            }
            temp.rotate(rotation, new Point64(0,0));
            temp.scale(magnification);
            scaled_and_rotated.Add(temp);
        }

        List<GCPolygon> ret = [];

        if (offsets.Count > 0)
        {
            foreach (Point64 offset in offsets)
            {
                foreach (var temp in scaled_and_rotated.Select(poly => new GCPolygon(poly)))
                {
                    temp.move(offset);
                    ret.Add(temp);
                }
            }
        }

        if (coords.Count > 0)
        {
            for (int idx = 0; idx < coords.Count; idx += 2)
            {
                foreach (var temp in scaled_and_rotated.Select(poly => new GCPolygon(poly)))
                {
                    temp.move(new Point64(coords[idx], coords[idx+1]));
                    ret.Add(temp);
                }
            }
        }

        if (!bypass_regular_array)
        {
            for (int y = 0; y < rows; y++)
            {
                for (int x = 0; x < columns; x++)
                {
                    foreach (var temp in scaled_and_rotated.Select(poly => new GCPolygon(poly)))
                    {
                        temp.move(new Point64(x * (rowVector.X + colVector.X), y * (rowVector.Y + colVector.Y)));
                        ret.Add(temp);
                    }
                }
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