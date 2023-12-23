using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Runtime.InteropServices.JavaScript;
using Clipper2Lib;
using geoWrangler;

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

    public PointD spacing; // rectangular spacing
    public PointD rowVector; // regular axis 1
    public PointD colVector; // regular axis 2

    private List<PointD> offsets; // explicit
    public List<double> coords; // explicitx and explicity

    public Repetition()
    {
        type = RepetitionType.None;
        columns = 0;
        rows = 0;
        spacing = new();
        rowVector = new();
        colVector = new();
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

    // Not sure whether these methods are needed, but added for now ahead of investigations.
    public PathD get_offsets()
    {
        PathD result = new();
        int count = get_count();
        switch (type)
        {
            case RepetitionType.Rectangular:
                for (int i = 0; i < columns; i++)
                {
                    double cx = i * spacing.x;
                    for (int j = 0; j < rows; j++)
                    {
                        result.Append(new(cx, j * spacing.y));
                    }
                }

                break;
            case RepetitionType.Regular:
                for (int i = 0; i < columns; i++)
                {
                    PointD vi = new(rowVector.x * i, rowVector.y * i);
                    for (int j = 0; j < rows; j++)
                    {
                        result.Append(new(vi.x + j * colVector.x, vi.y + j * colVector.y));
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

    public PathD get_extrema()
    {
        PathD result = new();

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
                        result.Append(new PointD(0, 0));
                    }
                    else
                    {
                        result.Append(new PointD(0, 0));
                        result.Append(new PointD(0, (rows - 1) * spacing.y));
                    }
                }
                else
                {
                    if (rows == 1)
                    {
                        result.Append(new PointD(0, 0));
                        result.Append(new PointD((columns - 1) * spacing.x, 0));
                    }
                    else
                    {
                        result.Append(new PointD(0, 0));
                        result.Append(new PointD(0, (rows - 1) * spacing.y));
                        result.Append(new PointD((columns - 1) * spacing.x, 0));
                        result.Append(new PointD((columns - 1) * spacing.x, (rows - 1) * spacing.y));
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
                        result.Append(new PointD(0, 0));
                    }
                    else
                    {
                        result.Append(new PointD(0, 0));
                        result.Append(new((rows - 1) * colVector.x, (rows - 1) * colVector.y));
                    }
                }
                else
                {
                    if (rows == 1)
                    {
                        result.Append(new PointD(0, 0));
                        result.Append(new((columns - 1) * rowVector.x, (columns - 1) * rowVector.y));
                    }
                    else
                    {
                        PointD vi = new((columns - 1) * rowVector.x, (columns - 1) * rowVector.y);
                        PointD vj = new((rows - 1) * colVector.x, (rows - 1) * colVector.y);
                        result.Append(new PointD(0, 0));
                        result.Append(new PointD(vi));
                        result.Append(new PointD(vj));
                        result.Append(new PointD(vi.x + vj.x, vi.y + vj.y));
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

                PointD vxmin = new(offsets[0]);
                PointD vxmax = new(offsets[0]);
                PointD vymin = new(offsets[0]);
                PointD vymax = new(offsets[0]);
                foreach (PointD v in offsets)
                {
                    if (v.x < vxmin.x)
                    {
                        vxmin = new(v);
                    }
                    else if (v.x > vxmax.x)
                    {
                        vxmax = new(v);
                    }

                    if (v.y < vymin.y)
                    {
                        vymin = new(v);
                    }
                    else if (v.y > vymax.y)
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

    public void transform(double magnification, bool x_reflection, double rotation)
    {
        if (type == RepetitionType.None)
        {
            return;
        }
        switch (type)
        {
            case RepetitionType.Rectangular:
                if (magnification != 1)
                {
                    spacing = new (spacing.x * magnification, spacing.y * magnification);
                }
                if (x_reflection || rotation != 0)
                {
                    PointD v = new(spacing);
                    if (x_reflection) v.y = -v.y;
                    double ca = Math.Cos(rotation);
                    double sa = Math.Sin(rotation);
                    type = RepetitionType.Regular;
                    rowVector.x = v.x * ca;
                    rowVector.y = v.x * sa;
                    colVector.x = -v.y * sa;
                    colVector.y = v.y * ca;
                }
                break;
            case RepetitionType.Regular:
                if (magnification != 1)
                {
                    rowVector = new (rowVector.x * magnification, rowVector.y * magnification);
                    rowVector = new (colVector.x * magnification, colVector.y * magnification);
                }

                if (x_reflection)
                {
                    rowVector.y = -rowVector.y;
                    colVector.y = -colVector.y;
                }

                if (rotation != 0)
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
                    rowVector = new PointD(rowVector.x * z2x.Real - rowVector.x * z2x.Imaginary, rowVector.y * z2y.Imaginary + rowVector.y * z2y.Imaginary);
                    colVector = new PointD(colVector.x * z2x.Real - colVector.x * z2x.Imaginary, colVector.y * z2y.Imaginary + colVector.y * z2y.Imaginary);
                }
                break;
            case RepetitionType.ExplicitX:
                if (rotation != 0)
                {
                    double ca = magnification * Math.Cos(rotation);
                    double sa = magnification * Math.Sin(rotation);
                    PathD temp = new ();
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
                if (rotation != 0)
                {
                    double ca = magnification * Math.Cos(rotation);
                    double sa = -magnification * Math.Sin(rotation);
                    if (x_reflection)
                    {
                        ca = -ca;
                        sa = -sa;
                    }

                    PathD temp = new ();
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
                if (rotation != 0)
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
                            offsets[i] = new PointD(offsets[i].x * z2x.Real - offsets[i].x * z2x.Imaginary, offsets[i].y * z2y.Imaginary + offsets[i].y * z2y.Imaginary);
                        }
                    }
                    else
                    {
                        for (int i = 0; i < offsets.Count; i++)
                        {
                            offsets[i] = new PointD(offsets[i].x * z2x.Real - offsets[i].x * z2x.Imaginary, offsets[i].y * z2y.Imaginary + offsets[i].y * z2y.Imaginary);
                        }
                    }
                }
                else if (x_reflection && magnification != 1)
                {
                    for (int i = 0; i < offsets.Count; i++)
                    {
                        offsets[i] = new(offsets[i].x * magnification, -offsets[i].y * magnification);
                    }
                }
                else if (x_reflection)
                {
                    for (int i = 0; i < offsets.Count; i++)
                    {
                        offsets[i] = new (offsets[i].x, -offsets[i].y);
                    }
                }
                else if (magnification != 1)
                {
                    for (int i = 0; i < offsets.Count; i++)
                    {
                        offsets[i] = new (offsets[i].x * magnification, offsets[i].y * magnification);
                    }
                }
            }
                break;
            default:
                return;
        }
    }
}