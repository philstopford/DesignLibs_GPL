using gds;
using geoWrangler;
using oasis;
using System;
using System.Collections.Generic;
using Clipper2Lib;

namespace geoCoreLib;

public class GCPolygon : GCElement
{
    public bool text;
    public string name;
    public Path64 pointarray { get; set; }

    public GCPolygon()
    {
        pGCPolygon();
    }

    private void pGCPolygon()
    {
        pointarray = new();
    }

    public GCPolygon(GCPolygon source)
    {
        pGCPolygon(source);
    }

    private void pGCPolygon(GCPolygon source)
    {
        pointarray = new (source.pointarray);
        layer_nr = source.layer_nr;
        datatype_nr = source.datatype_nr;
    }

    public GCPolygon(Path64 points, int layer, int datatype)
    {
        pGCPolygon(points, layer, datatype);
    }

    private void pGCPolygon(Path64 points, int layer, int datatype)
    {
        pointarray = new(points);
        layer_nr = layer;
        datatype_nr = datatype;
    }

    public override void rotate(double angleDegree, Point64 pos)
    {
        pRotate(angleDegree, pos);
    }

    private void pRotate(double angleDegree, Point64 pos)
    {
        switch (angleDegree)
        {
            case 0:
                return;
            default:
                pointarray = GeoWrangler.Rotate(pos, pointarray, angleDegree);
                break;
        }
    }

    public override void scale(double size)
    {
        pScale(size);
    }

    private void pScale(double size)
    {
        pointarray = GeoWrangler.resize(pointarray, size);
    }

    public override void resize(double factor)
    {
        for (int i = 0; i < pointarray.Count; i++)
        {
            pointarray[i] = new (pointarray[i].X * factor, pointarray[i].Y * factor);
        }
    }

    public override void scale(Point64 origin, double size)
    {
        pScale(origin, size);
    }

    private void pScale(Point64 origin, double size)
    {
        pointarray = GeoWrangler.resize(origin, pointarray, size);
    }

    public void deletePoint(int pos)
    {
        pDeletePoint(pos);
    }

    private void pDeletePoint(int pos)
    {
        Path64 newPointArray = Helper.initedPath64(pointarray.Count - 1);
        for (int i = pos; i < pointarray.Count - 1; i++)
        {
            newPointArray[i] = pointarray[i + 1];
        }
        pointarray = newPointArray;
    }

    public void addPoint(int pos)
    {
        pAddPoint(pos);
    }

    private void pAddPoint(int pos)
    {
        Path64 newPointArray = Helper.initedPath64(pointarray.Count + 1);
        for (int pt = newPointArray.Count - 1; pt > pos; pt--)
        {
            newPointArray[pt] = pointarray[pt - 1];
        }
        pointarray = newPointArray;
    }

    public override void mapSelect(GCStrans m)
    {
        pMapSelect(m);
    }

    private void pMapSelect(GCStrans m)
    {
        switch (select)
        {
            case true:
                m.matrix.TransformPoints(pointarray);
                break;
        }
    }

    public override void map(GCStrans m)
    {
        pMap(m);
    }

    private void pMap(GCStrans m)
    {
        m.matrix.TransformPoints(pointarray);
    }

    public override void minimum(ref Point64 pos)
    {
        pMinimum(ref pos);
    }

    private void pMinimum(ref Point64 pos)
    {
        Point64 t = GeoWrangler.getMinimumPoint(pointarray);

        pos.X = Math.Min(pos.X, t.X);
        pos.Y = Math.Min(pos.Y, t.Y);
    }

    public override void maximum(ref Point64 pos)
    {
        pMaximum(ref pos);
    }

    private void pMaximum(ref Point64 pos)
    {
        Point64 t = GeoWrangler.getMaximumPoint(pointarray);

        pos.X = Math.Max(pos.X, t.X);
        pos.Y = Math.Max(pos.Y, t.Y);
    }

    public override void moveSelect(Point64 pos)
    {
        pMoveSelect(pos);
    }

    private void pMoveSelect(Point64 pos)
    {
        switch (select)
        {
            case true:
            {
                pointarray = GeoWrangler.move(pointarray, pos.X, pos.Y);
                break;
            }
        }
    }

    public override void move(Point64 pos)
    {
        pMove(pos);
    }

    private void pMove(Point64 pos)
    {
        pointarray = GeoWrangler.move(pointarray, pos.X, pos.Y);
    }

    public override List<GCPolygon> convertToPolygons()
    {
        return pConvertToPolygons();
    }

    private List<GCPolygon> pConvertToPolygons()
    {
        List<GCPolygon> ret = new() {new GCPolygon(this)};
        return ret;
    }

    public override void clean()
    {
        pClean();
    }

    private void pClean()
    {
        Point64 p = new(0, 0);
        double a;
        int anz = 0;
        for (a = 0; a is < 350 or > 370 && anz < 3;)
        {
            a = 0;
            for (int i = 0; i < pointarray.Count - 1; i++)
            {
                switch (pointarray.Count)
                {
                    case < 4:
                        return; //no area
                }
                if (pointarray[i] == pointarray[i + 1])
                {
                    deletePoint(i + 1);
                    if (pointarray.Count == i + 1)
                    {
                        pointarray[0] = pointarray[i];
                    }
                    i--;
                }
                if (i < pointarray.Count - 2)
                {
                    if (nearlyParallel(pointarray[i], pointarray[i + 1], pointarray[i + 1], pointarray[i + 2]))
                    {
                        deletePoint(i + 1);
                        i = 0;
                    }
                }
                switch (pointarray.Count)
                {
                    case > 3:
                    {
                        while (pointarray.Count > 3 && parallel(pointarray[0], pointarray[1], pointarray[^1], pointarray[^2]))
                        {
                            deletePoint(pointarray.Count - 1);
                            pointarray[0] = pointarray[^1];
                        }

                        break;
                    }
                }
            }
            // sort points
            for (int i = 0; i < pointarray.Count - 2; i++)
            {
                a += angle(pointarray[i], pointarray[i + 1], pointarray[i + 2]);
            }
            switch (pointarray.Count)
            {
                case > 3:
                    a += angle(pointarray[^2], pointarray[0], pointarray[1]);
                    break;
            }

            switch (a)
            {
                case < -185:
                {
                    for (int i = 1; i < pointarray.Count / 2; i++)
                    {
                        p = pointarray[pointarray.Count - i - 1];
                        pointarray[pointarray.Count - i - 1] = pointarray[i];
                        pointarray[i] = p;
                    }

                    break;
                }
            }

            switch (a)
            {
                // (simple self-intersecting)
                case > 370:
                case < 350 and > -350:
                {
                    for (int i = 0; i < pointarray.Count - 1; i++)
                    {
                        for (int j = i + 2; j < pointarray.Count - 1; j++)
                        {
                            bool b = cutPoint2(pointarray[i], pointarray[i + 1], pointarray[j], pointarray[j + 1], p);
                            switch (b)
                            {
                                case true:
                                {
                                    switch (i == 0 && j == pointarray.Count - 2)
                                    {
                                        case false:
                                        {
                                            addPoint(j + 1);
                                            j++;
                                            pointarray[j] = p;
                                            addPoint(i + 1);
                                            i++;
                                            pointarray[i] = p;
                                            for (int k = i + 1; k <= (i + j) / 2; k++)
                                            {
                                                p = pointarray[j - k + i + 1];
                                                pointarray[j - k + i + 1] = pointarray[k];
                                                pointarray[k] = p;
                                            }
                                            j = pointarray.Count;
                                            i = pointarray.Count;
                                            break;
                                        }
                                    }

                                    break;
                                }
                            }
                        }
                    }

                    break;
                }
            }

            switch (pointarray.Count)
            {
                // remove self-intersection 
                case >= 8:
                {
                    bool ende2 = false;
                    for (int i = 0; i < pointarray.Count - 1; i++)
                    {
                        bool ende1 = false;
                        for (int j = i + 2; j < pointarray.Count - 1; j++)
                        {
                            if (identical(pointarray[i], pointarray[i + 1], pointarray[j + 1], pointarray[j]))
                            {
                                int h1 = 0, h2;
                                bool side = true;
                                bool change;
                                if (pointarray[j + 1] == pointarray[i])
                                {
                                    change = true;
                                    h1 = h1 switch
                                    {
                                        -1 => pointarray.Count - 2,
                                        _ => i - 1
                                    };
                                    h2 = j + 2;
                                    if (h2 == pointarray.Count)
                                    {
                                        h2 = 1;
                                    }

                                    double d1 = angle(pointarray[h1], pointarray[i], pointarray[i + 1]);
                                    double d2 = angle(pointarray[h2], pointarray[j + 1], pointarray[j]);
                                    switch (d1)
                                    {
                                        case < 0:
                                            d1 += 360;
                                            break;
                                    }
                                    switch (d2)
                                    {
                                        case < 0:
                                            d2 += 360;
                                            break;
                                    }
                                    if (d1 < d2)
                                    {
                                        side = false;
                                    }
                                }
                                else if (onLine2(pointarray[i], pointarray[i + 1], pointarray[j + 1]))
                                {
                                    change = true;
                                    h2 = j + 2;
                                    if (h2 == pointarray.Count)
                                    {
                                        h2 = 1;
                                    }
                                    if (distance(pointarray[i], pointarray[i + 1], pointarray[h2]) < 0)
                                    {
                                        side = false;
                                    }
                                }
                                else if (onLine2(pointarray[j], pointarray[j + 1], pointarray[i]))
                                {
                                    change = true;
                                    h1 = h1 switch
                                    {
                                        -1 => pointarray.Count - 2,
                                        _ => i - 1
                                    };
                                    if (distance(pointarray[j], pointarray[j + 1], pointarray[h1]) < 0)
                                    {
                                        side = false;
                                    }
                                }
                                else
                                {
                                    break;
                                }
                                while (pointarray[j] == pointarray[i + 1])
                                {
                                    j--; i++;
                                    switch (j)
                                    {
                                        case -1:
                                            j = pointarray.Count - 2; ende1 = true;
                                            break;
                                    }

                                    if (i != pointarray.Count)
                                    {
                                        continue;
                                    }

                                    i = 1; ende2 = true;
                                }
                                h2 = i + 1;
                                if (h2 == pointarray.Count)
                                {
                                    h2 = 1;
                                }
                                int h3 = h2 + 1;
                                if (h3 == pointarray.Count)
                                {
                                    h3 = 1;
                                }

                                h1 = h1 switch
                                {
                                    -1 => pointarray.Count - 2,
                                    _ => j - 1
                                };
                                int h4 = j + 1;
                                if (h4 == pointarray.Count)
                                {
                                    h4 = 1;
                                }
                                int h5 = h4 + 1;
                                if (h5 == pointarray.Count)
                                {
                                    h5 = 1;
                                }
                                int h6 = i - 1;
                                h6 = h6 switch
                                {
                                    -1 => pointarray.Count - 2,
                                    _ => h6
                                };
                                if (onLine2(pointarray[i], pointarray[h2], pointarray[j]))
                                {
                                    if (distance(pointarray[i], pointarray[h2], pointarray[h1]) < 0 != side)
                                    {
                                        change = false;
                                    }
                                }
                                else if (onLine2(pointarray[j], pointarray[h4], pointarray[h2]))
                                {
                                    if (distance(pointarray[j], pointarray[h4], pointarray[h3]) < 0 != side)
                                    {
                                        change = false;
                                    }
                                }
                                else
                                {
                                    double d1 = angle(pointarray[h6], pointarray[i], pointarray[h2]);
                                    double d2 = angle(pointarray[h5], pointarray[h4], pointarray[j]);
                                    switch (d1)
                                    {
                                        case < 0:
                                            d1 += 360;
                                            break;
                                    }
                                    switch (d2)
                                    {
                                        case < 0:
                                            d2 += 360;
                                            break;
                                    }
                                    if (d1 < d2 != side)
                                    {
                                        change = false;
                                    }
                                }
                                switch (change)
                                {
                                    case true:
                                    {
                                        for (int k = i + 1; k < (i + j + 2) / 2; k++)
                                        {
                                            p = pointarray[j + i - k + 1];
                                            pointarray[j + i - k + 1] = pointarray[k];
                                            pointarray[k] = p;
                                        }

                                        break;
                                    }
                                }
                            }
                            if (ende1 || ende2)
                            {
                                break;
                            }
                        }
                        if (ende2)
                        {
                            break;
                        }
                    }

                    break;
                }
            }
            anz++;
        }
    }

    public override void saveGDS(gdsWriter gw)
    {
        pSaveGDS(gw);
    }

    private void pSaveGDS(gdsWriter gw)
    {
        Point64 pc = new();
        CircleResult res = isCircle();
        if (res.circle)
        {
            GCPath tPath = new(new () { pc }, layer_nr, datatype_nr);
            tPath.setWidth((int)(res.radius * 2));
            tPath.setCap(1); // set cap to 'round' type. Also forced in the saver for a round path type, as a safety.
            tPath.setRound(true); // this is the important thing to set - it forces the correct handling (including cap value)
            tPath.saveGDS(gw);
        }
        else
        {
            // polygon (boundary)
            gw.bw.Write((ushort)4);
            gw.bw.Write(gdsValues.sBOUNDARY);
            //layer
            gw.bw.Write((ushort)6);
            gw.bw.Write(gdsValues.sLAYER);
            gw.bw.Write((short)layer_nr);
            //datatype
            gw.bw.Write((ushort)6);
            gw.bw.Write(gdsValues.sDATATYPE);
            gw.bw.Write((short)datatype_nr);

            Path64 cleaned = GeoWrangler.removeDuplicates(pointarray);
            int i = cleaned.Count;
            i = i switch
            {
                > 8191 => 8191,
                _ => i
            };
            //xy 
            bool closedGeometry = ((cleaned[0].X == cleaned[^1].X) && (cleaned[0].Y == cleaned[^1].Y));
                
            // Add one to the point-list length value (i) to mark the closed shape, if needed.
            int val = (closedGeometry ? i: i + 1) * 2 * 4 + 4;
            gw.bw.Write((ushort)val);
            gw.bw.Write(gdsValues.sXY);

            for (int k = 0; k < i; k++)
            {
                gw.bw.Write((int)cleaned[k].X);
                gw.bw.Write((int)cleaned[k].Y);
            }

            // close the polygon.
            if (!closedGeometry)
            {
                gw.bw.Write((int)cleaned[0].X);
                gw.bw.Write((int)cleaned[0].Y);
            }

            // endel
            gw.bw.Write((ushort)4);
            gw.bw.Write(gdsValues.sENDEL);
        }
    }

    public override void saveOASIS(oasWriter ow)
    {
        pSaveOASIS(ow);
    }

    class CircleResult
    {
        public bool circle { get; set; }
        public double radius { get; set; }
    }
    private CircleResult isCircle()
    {
        CircleResult result = new() {circle = false, radius = 0};
        switch (pointarray.Count)
        {
            case < 10:
                return result;
        }

        PointD mid = GeoWrangler.midPoint(pointarray);

        double min_distance = Math.Abs(GeoWrangler.distanceBetweenPoints(mid, pointarray[0]));
        double max_distance = Math.Abs(GeoWrangler.distanceBetweenPoints(mid, pointarray[0]));

        for (int i = 1; i < pointarray.Count; i++)
        {
            double distance = Math.Abs(GeoWrangler.distanceBetweenPoints(mid, pointarray[i]));
            min_distance = Math.Min(min_distance, distance);
            max_distance = Math.Max(max_distance, distance);
        }

        double min_area = Math.PI * min_distance * min_distance;
        double max_area = Math.PI * max_distance * max_distance;
        double actual_area = Clipper.Area(pointarray);

        if ((actual_area >= min_area) && (actual_area <= max_area))
        {
            result.circle = true;
            // 100 to account for 0.01 nm grid
            result.radius = 100 * Math.Sqrt(actual_area / Math.PI);
        }

        return result;
    }

    private void pSaveOASIS(oasWriter ow)
    {
        byte info_byte;
        switch (GCSetup.oasisSaveCircle)
        {
            case true:
            {
                // pc and r are manipulated in isCircle()
                Point64 pc = new();
                CircleResult res = isCircle();
                if (res.circle)
                {
                    ow.modal.absoluteMode = ow.modal.absoluteMode switch
                    {
                        // save as circle
                        false => true,
                        _ => ow.modal.absoluteMode
                    };
                    info_byte = 0;
                    if (layer_nr != ow.modal.layer)
                    {
                        info_byte += 1;
                    }
                    if (datatype_nr != ow.modal.datatype)
                    {
                        info_byte += 2;
                    }
                    if (pc.X != ow.modal.geometry_x)
                    {
                        info_byte += 16;
                    }
                    if (pc.Y != ow.modal.geometry_y)
                    {
                        info_byte += 8;
                    }
                    switch (Math.Abs(res.radius - ow.modal.circle_radius))
                    {
                        case > double.Epsilon:
                            info_byte += 32;
                            break;
                    }
                    ow.writeUnsignedInteger(oasValues.CIRCLE);
                    ow.writeRaw(info_byte);
                    switch (info_byte & 1)
                    {
                        case > 0:
                            ow.modal.layer = layer_nr;
                            ow.writeUnsignedInteger((uint)ow.modal.layer);
                            break;
                    }
                    switch (info_byte & 2)
                    {
                        case > 0:
                            ow.modal.datatype = datatype_nr;
                            ow.writeUnsignedInteger((uint)datatype_nr);
                            break;
                    }
                    switch (info_byte & 32)
                    {
                        case > 0:
                            ow.modal.circle_radius = res.radius;
                            ow.writeUnsignedInteger((uint)ow.modal.circle_radius);
                            break;
                    }
                    switch (info_byte & 16)
                    {
                        case > 0:
                            ow.modal.geometry_x = (int)pc.X;
                            ow.writeSignedInteger(ow.modal.geometry_x);
                            break;
                    }
                    switch (info_byte & 8)
                    {
                        case > 0:
                            ow.modal.geometry_y = (int)pc.Y;
                            ow.writeSignedInteger(ow.modal.geometry_y);
                            break;
                    }
                    return;

                }

                break;
            }
        }

        switch (GCSetup.oasisSaveCtrapezoid)
        {
            // check if ctrapezoid
            case true when pointarray.Count is 4 or 5:
            {
                int form = 0;
                int off = 0;
                for (int i = 1; i < pointarray.Count; i++)
                {
                    PointD pd = GeoWrangler.distanceBetweenPoints_point(pointarray[i], pointarray[i - 1]);
                    switch (Math.Abs(pd.x - pd.y))
                    {
                        case <= double.Epsilon when pd.x > 0:
                            form = form * 10 + 9;
                            break;
                        case <= double.Epsilon when pd.x < 0:
                            form = form * 10 + 1;
                            break;
                        default:
                        {
                            switch (Math.Abs(pd.x - -pd.y))
                            {
                                case <= double.Epsilon when pd.x < 0:
                                    form = form * 10 + 7;
                                    break;
                                case <= double.Epsilon when pd.x > 0:
                                    form = form * 10 + 3;
                                    break;
                                default:
                                {
                                    form = pd.x switch
                                    {
                                        0 when pd.y > 0 => form * 10 + 8,
                                        0 when pd.y < 0 => form * 10 + 2,
                                        _ => pd.y switch
                                        {
                                            0 when pd.x > 0 => form * 10 + 6,
                                            0 when pd.x < 0 => form * 10 + 4,
                                            _ => -1
                                        }
                                    };

                                    break;
                                }
                            }

                            break;
                        }
                    }
                    if (form == -1)
                    {
                        break;
                    }
                }
                if (form != -1)
                {
                    int w;
                    int h;
                    switch (pointarray.Count)
                    {
                        case 4:
                            switch (form)
                            {
                                case 726: off = 2; form = 672; break;
                                case 267: off = 1; form = 672; break;
                                case 429: off = 2; form = 942; break;
                                case 294: off = 1; form = 942; break;
                                case 816: off = 2; form = 681; break;
                                case 168: off = 1; form = 681; break;
                                case 384: off = 2; form = 438; break;
                                case 843: off = 1; form = 438; break;
                                case 167: off = 1; form = 671; break;
                                case 716: off = 2; form = 671; break;
                                case 394: off = 1; form = 943; break;
                                case 439: off = 2; form = 943; break;
                                case 297: off = 1; form = 972; break;
                                case 729: off = 2; form = 972; break;
                                case 381: off = 1; form = 813; break;
                                case 138: off = 2; form = 813; break;
                            }
                            switch (form)
                            {
                                case 672: ow.writeCtrapezoid(layer_nr, 16, (int)pointarray[0 + off].X, (int)pointarray[0 + off].Y, (int)(pointarray[1 + off].X - pointarray[0 + off].X), 0, datatype_nr); return;
                                case 942: ow.writeCtrapezoid(layer_nr, 17, (int)pointarray[0 + off].X, (int)pointarray[0 + off].Y, (int)(pointarray[1 + off].X - pointarray[0 + off].X), 0, datatype_nr); return;
                                case 681: ow.writeCtrapezoid(layer_nr, 18, (int)pointarray[0 + off].X, (int)pointarray[0 + off].Y, (int)(pointarray[1 + off].X - pointarray[0 + off].X), 0, datatype_nr); return;
                                case 438: w = (int)(pointarray[0 + off].X - pointarray[1 + off].X); ow.writeCtrapezoid(layer_nr, 19, (int)pointarray[0 + off].X - w, (int)pointarray[0 + off].Y - w, w, 0, datatype_nr); return;
                                case 671: h = (int)((pointarray[1 + off].X - pointarray[0 + off].X) / 2); ow.writeCtrapezoid(layer_nr, 20, (int)pointarray[0 + off].X, (int)pointarray[0 + off].Y, 0, h, datatype_nr); return;
                                case 943: h = (int)(pointarray[1 + off].Y - pointarray[0 + off].Y); ow.writeCtrapezoid(layer_nr, 21, (int)pointarray[0 + off].X - h, (int)pointarray[0 + off].Y, 0, h, datatype_nr); return;
                                case 972: w = (int)(pointarray[1 + off].X - pointarray[0 + off].X); ow.writeCtrapezoid(layer_nr, 22, (int)pointarray[0 + off].X, (int)pointarray[0 + off].Y, w, 0, datatype_nr); return;
                                case 813: w = (int)((pointarray[1 + off].Y - pointarray[0 + off].Y) / 2); ow.writeCtrapezoid(layer_nr, 23, (int)pointarray[0 + off].X - w, (int)pointarray[0 + off].Y, w, 0, datatype_nr); return;
                            }

                            break;
                        case 5:
                            switch (form)
                            {
                                //type 24
                                case 2684: off = 1; form = 6842; break;
                                case 4268: off = 2; form = 6842; break;
                                case 8426: off = 3; form = 6842; break;
                                //type 0
                                case 2674: off = 1; form = 6742; break;
                                case 4267: off = 2; form = 6742; break;
                                case 7426: off = 3; form = 6742; break;
                                //type 1
                                case 2694: off = 1; form = 6942; break;
                                case 4269: off = 2; form = 6942; break;
                                case 9426: off = 3; form = 6942; break;
                                //type 2
                                case 1684: off = 1; form = 6841; break;
                                case 4168: off = 2; form = 6841; break;
                                case 8416: off = 3; form = 6841; break;
                                //type 3
                                case 3684: off = 1; form = 6843; break;
                                case 4368: off = 2; form = 6843; break;
                                case 8436: off = 3; form = 6843; break;
                                //type 4
                                case 1674: off = 1; form = 6741; break;
                                case 4167: off = 2; form = 6741; break;
                                case 7416: off = 3; form = 6741; break;
                                //type 5
                                case 3694: off = 1; form = 6943; break;
                                case 4369: off = 2; form = 6943; break;
                                case 9436: off = 3; form = 6943; break;
                                //type 6
                                case 1694: off = 1; form = 6941; break;
                                case 4169: off = 2; form = 6941; break;
                                case 9416: off = 3; form = 6941; break;
                                //type 7
                                case 3674: off = 1; form = 6743; break;
                                case 4367: off = 2; form = 6743; break;
                                case 7436: off = 3; form = 6743; break;
                                //type 8
                                case 2687: off = 1; form = 6872; break;
                                case 7268: off = 2; form = 6872; break;
                                case 8726: off = 3; form = 6872; break;
                                //type 9
                                case 2681: off = 1; form = 6812; break;
                                case 1268: off = 2; form = 6812; break;
                                case 8126: off = 3; form = 6812; break;
                                //type 10
                                case 2984: off = 1; form = 9842; break;
                                case 4298: off = 2; form = 9842; break;
                                case 8429: off = 3; form = 9842; break;
                                //type 11
                                case 3842: off = 1; form = 8423; break;
                                case 2384: off = 2; form = 8423; break;
                                case 4238: off = 3; form = 8423; break;
                                //type 12
                                case 2987: off = 1; form = 9872; break;
                                case 7298: off = 2; form = 9872; break;
                                case 8729: off = 3; form = 9872; break;
                                //type 13
                                case 3812: off = 1; form = 8123; break;
                                case 2381: off = 2; form = 8123; break;
                                case 1238: off = 3; form = 8123; break;
                                //type 14
                                case 2981: off = 1; form = 9812; break;
                                case 1298: off = 2; form = 9812; break;
                                case 8129: off = 3; form = 9812; break;
                                //type 15
                                case 3872: off = 1; form = 8723; break;
                                case 2387: off = 2; form = 8723; break;
                                case 7238: off = 3; form = 8723; break;

                            }
                            switch (form)
                            {
                                case 6842:
                                    w = (int)(pointarray[1 + off].X - pointarray[0 + off].X);
                                    h = off switch
                                    {
                                        3 => (int)(pointarray[1].Y - pointarray[2].Y),
                                        _ => (int)(pointarray[2 + off].Y - pointarray[1 + off].Y)
                                    };
                                    if (h == w)
                                    {
                                        ow.writeCtrapezoid(layer_nr, 25, (int)pointarray[0 + off].X, (int)pointarray[0 + off].Y, w, 0, datatype_nr);
                                    }
                                    else
                                    {
                                        ow.writeCtrapezoid(layer_nr, 24, (int)pointarray[0 + off].X, (int)pointarray[0 + off].Y, w, h, datatype_nr);
                                    }
                                    return;
                                case 6742:
                                    w = (int)(pointarray[1 + off].X - pointarray[0 + off].X);
                                    h = off switch
                                    {
                                        0 => (int)(pointarray[2].Y - pointarray[1].Y),
                                        _ => (int)(pointarray[-1 + off].Y - pointarray[0 + off].Y)
                                    };
                                    ow.writeCtrapezoid(layer_nr, 0, (int)pointarray[0 + off].X, (int)pointarray[0 + off].Y, w, h, datatype_nr);
                                    return;
                                case 6942:
                                    h = off switch
                                    {
                                        0 => (int)(pointarray[2].Y - pointarray[1].Y),
                                        _ => (int)(pointarray[-1 + off].Y - pointarray[0 + off].Y)
                                    };
                                    w = (int)(pointarray[1 + off].X - pointarray[0 + off].X + h);
                                    ow.writeCtrapezoid(layer_nr, 1, (int)pointarray[0 + off].X, (int)pointarray[0 + off].Y, w, h, datatype_nr);
                                    return;
                                case 6841:
                                    w = (int)(pointarray[1 + off].X - pointarray[0 + off].X);
                                    h = off switch
                                    {
                                        0 => (int)(pointarray[2].Y - pointarray[1].Y),
                                        _ => (int)(pointarray[-1 + off].Y - pointarray[0 + off].Y)
                                    };
                                    ow.writeCtrapezoid(layer_nr, 2, (int)pointarray[0 + off].X, (int)pointarray[0 + off].Y, w, h, datatype_nr);
                                    return;
                                case 6843:
                                    h = off switch
                                    {
                                        0 => (int)(pointarray[2].Y - pointarray[1].Y),
                                        _ => (int)(pointarray[-1 + off].Y - pointarray[0 + off].Y)
                                    };
                                    w = (int)(pointarray[1 + off].X - pointarray[0 + off].X + h);
                                    ow.writeCtrapezoid(layer_nr, 3, (int)(pointarray[0 + off].X - h), (int)pointarray[0 + off].Y, w, h, datatype_nr);
                                    return;
                                case 6741:
                                    w = (int)(pointarray[1 + off].X - pointarray[0 + off].X);
                                    h = off switch
                                    {
                                        0 => (int)(pointarray[2].Y - pointarray[1].Y),
                                        _ => (int)(pointarray[-1 + off].Y - pointarray[0 + off].Y)
                                    };
                                    ow.writeCtrapezoid(layer_nr, 4, (int)pointarray[0 + off].X, (int)pointarray[0 + off].Y, w, h, datatype_nr);
                                    return;
                                case 6943:
                                    h = off switch
                                    {
                                        0 => (int)(pointarray[2].Y - pointarray[1].Y),
                                        _ => (int)(pointarray[-1 + off].Y - pointarray[0 + off].Y)
                                    };
                                    w = (int)(pointarray[1 + off].X - pointarray[0 + off].X + h + h);
                                    ow.writeCtrapezoid(layer_nr, 5, (int)(pointarray[0 + off].X - h), (int)pointarray[0 + off].Y, w, h, datatype_nr);
                                    return;
                                case 6941:
                                    h = off switch
                                    {
                                        0 => (int)(pointarray[2].Y - pointarray[1].Y),
                                        _ => (int)(pointarray[-1 + off].Y - pointarray[0 + off].Y)
                                    };
                                    w = (int)(pointarray[1 + off].X - pointarray[0 + off].X + h);
                                    ow.writeCtrapezoid(layer_nr, 6, (int)pointarray[0 + off].X, (int)pointarray[0 + off].Y, w, h, datatype_nr);
                                    return;
                                case 6743:
                                    h = off switch
                                    {
                                        0 => (int)(pointarray[2].Y - pointarray[1].Y),
                                        _ => (int)(pointarray[-1 + off].Y - pointarray[0 + off].Y)
                                    };
                                    w = (int)(pointarray[1 + off].X - pointarray[0 + off].X + h);
                                    ow.writeCtrapezoid(layer_nr, 7, (int)(pointarray[0 + off].X - h), (int)pointarray[0 + off].Y, w, h, datatype_nr);
                                    return;
                                case 6872:
                                    w = (int)(pointarray[1 + off].X - pointarray[0 + off].X);
                                    h = off switch
                                    {
                                        0 => (int)(pointarray[3].Y - pointarray[0].Y),
                                        _ => (int)(pointarray[-1 + off].Y - pointarray[0 + off].Y)
                                    };
                                    ow.writeCtrapezoid(layer_nr, 8, (int)pointarray[0 + off].X, (int)pointarray[0 + off].Y, w, h, datatype_nr);
                                    return;
                                case 6812:
                                    w = (int)(pointarray[1 + off].X - pointarray[0 + off].X);
                                    h = off switch
                                    {
                                        0 => (int)(pointarray[2].Y - pointarray[0].Y),
                                        _ => (int)(pointarray[-1 + off].Y - pointarray[0 + off].Y + w)
                                    };
                                    ow.writeCtrapezoid(layer_nr, 9, (int)pointarray[0 + off].X, (int)pointarray[0 + off].Y, w, h, datatype_nr);
                                    return;
                                case 9842:
                                    w = (int)(pointarray[1 + off].X - pointarray[0 + off].X);
                                    h = off switch
                                    {
                                        0 => (int)(pointarray[3].Y - pointarray[0].Y),
                                        _ => (int)(pointarray[-1 + off].Y - pointarray[0 + off].Y)
                                    };
                                    ow.writeCtrapezoid(layer_nr, 10, (int)pointarray[0 + off].X, (int)pointarray[0 + off].Y, w, h, datatype_nr);
                                    return;
                                case 8423:
                                    h = (int)(pointarray[1 + off].Y - pointarray[0 + off].Y);
                                    w = off switch
                                    {
                                        0 => (int)(pointarray[0].X - pointarray[2].X),
                                        _ => (int)(pointarray[0 + off].X - pointarray[-1 + off].X)
                                    };
                                    ow.writeCtrapezoid(layer_nr, 11, (int)(pointarray[0 + off].X - w), (int)pointarray[0 + off].Y, w, h, datatype_nr);
                                    return;
                                case 9872:
                                    w = (int)(pointarray[1 + off].X - pointarray[0 + off].X);
                                    h = off switch
                                    {
                                        0 => (int)(pointarray[3].Y - pointarray[0].Y),
                                        _ => (int)(pointarray[-1 + off].Y - pointarray[0 + off].Y)
                                    };
                                    ow.writeCtrapezoid(layer_nr, 12, (int)pointarray[0 + off].X, (int)pointarray[0 + off].Y, w, h, datatype_nr);
                                    return;
                                case 8123:
                                    h = (int)(pointarray[1 + off].Y - pointarray[0 + off].Y);
                                    w = off switch
                                    {
                                        0 => (int)(pointarray[0].X - pointarray[2].X),
                                        _ => (int)(pointarray[0 + off].X - pointarray[-1 + off].X)
                                    };
                                    ow.writeCtrapezoid(layer_nr, 13, (int)(pointarray[0 + off].X - w), (int)pointarray[0 + off].Y, w, h, datatype_nr);
                                    return;
                                case 9812:
                                    w = (int)(pointarray[1 + off].X - pointarray[0 + off].X);
                                    h = off switch
                                    {
                                        0 => (int)(pointarray[2].Y - pointarray[0].Y),
                                        _ => (int)(pointarray[-1 + off].Y - pointarray[0 + off].Y + w)
                                    };
                                    ow.writeCtrapezoid(layer_nr, 14, (int)pointarray[0 + off].X, (int)pointarray[0 + off].Y, w, h, datatype_nr);
                                    return;
                                case 8723:
                                    w = off switch
                                    {
                                        0 => (int)(pointarray[0].X - pointarray[2].X),
                                        _ => (int)(pointarray[0 + off].X - pointarray[-1 + off].X)
                                    };
                                    h = (int)(pointarray[1 + off].Y - pointarray[0 + off].Y + w);
                                    ow.writeCtrapezoid(layer_nr, 15, (int)(pointarray[0 + off].X - w), (int)pointarray[0 + off].Y, w, h, datatype_nr);
                                    return;
                            }

                            break;
                    }
                }

                break;
            }
        }

        switch (GCSetup.oasisSaveTrapezoid)
        {
            //trapezoid
            case true when pointarray.Count == 5:
            {
                int form = 0;
                int da, db;
                int minX = (int)pointarray[0].X;
                int maxX = (int)pointarray[0].X;
                int minY = (int)pointarray[0].Y;
                int maxY = (int)pointarray[0].Y;
                for (int i = 1; i < pointarray.Count; i++)
                {
                    if (pointarray[i].X > maxX)
                    {
                        maxX = (int)pointarray[i].X;
                    }
                    if (pointarray[i].Y > maxY)
                    {
                        maxY = (int)pointarray[i].Y;
                    }
                    if (pointarray[i].X < minX)
                    {
                        minX = (int)pointarray[i].X;
                    }
                    if (pointarray[i].Y < minY)
                    {
                        minY = (int)pointarray[i].Y;
                    }
                    PointD pd = GeoWrangler.distanceBetweenPoints_point(pointarray[i], pointarray[i - 1]);
                    switch (pd.x)
                    {
                        case 0:
                            form = form * 10 + 1;
                            break;
                        default:
                        {
                            switch (pd.y)
                            {
                                case 0:
                                    form = form * 10 + 2;
                                    break;
                                default:
                                    form *= 10;
                                    break;
                            }

                            break;
                        }
                    }
                }
                int w = maxX - minX;
                int h = maxY - minY;
                switch (form)
                {
                    case 1010:
                    case 1210:
                    case 1012:
                        if (pointarray[0].X > pointarray[2].X)
                        {
                            da = (int)(pointarray[0].Y - pointarray[3].Y);
                            db = (int)(pointarray[1].Y - pointarray[2].Y);
                        }
                        else
                        {
                            db = (int)(pointarray[0].Y - pointarray[3].Y);
                            da = (int)(pointarray[1].Y - pointarray[2].Y);
                        }
                        ow.writeTrapezoid(layer_nr, 1, minX, minY, w, h, -da, -db, datatype_nr);
                        return;
                    case 101:
                    case 121:
                    case 2101:
                        if (pointarray[1].X > pointarray[3].X)
                        {
                            da = (int)(pointarray[1].Y - pointarray[4].Y);
                            db = (int)(pointarray[2].Y - pointarray[3].Y);
                        }
                        else
                        {
                            db = (int)(pointarray[1].Y - pointarray[4].Y);
                            da = (int)(pointarray[2].Y - pointarray[3].Y);
                        }
                        ow.writeTrapezoid(layer_nr, 1, minX, minY, w, h, -da, -db, datatype_nr);
                        //printf("minX %d, maxX %d, minY %d, maxY %d, w %d, h %d, da %d, db %d\n",minX,maxX,minY,maxY,w,h,da,db);
                        return;
                    case 2020:
                    case 2120:
                    case 2021:
                        if (pointarray[0].Y < pointarray[2].Y)
                        {
                            da = (int)(pointarray[0].X - pointarray[3].X);
                            db = (int)(pointarray[1].X - pointarray[2].X);
                        }
                        else
                        {
                            db = (int)(-pointarray[0].X + pointarray[3].X);
                            da = (int)(-pointarray[1].X + pointarray[2].X);
                        }
                        ow.writeTrapezoid(layer_nr, 0, minX, minY, w, h, -da, -db, datatype_nr);
                        return;
                    case 202:
                    case 1202:
                    case 212:
                        if (pointarray[1].Y < pointarray[3].Y)
                        {
                            da = (int)(pointarray[1].X - pointarray[4].X);
                            db = (int)(pointarray[2].X - pointarray[3].X);
                        }
                        else
                        {
                            db = (int)(-pointarray[1].X + pointarray[4].X);
                            da = (int)(-pointarray[2].X + pointarray[3].X);
                        }
                        ow.writeTrapezoid(layer_nr, 0, minX, minY, w, h, -da, -db, datatype_nr);
                        return;
                }

                break;
            }
        }

        ow.modal.absoluteMode = ow.modal.absoluteMode switch
        {
            // save as polygon
            false => true,
            _ => ow.modal.absoluteMode
        };
        info_byte = 32;  //write point-list;
        if (layer_nr != ow.modal.layer)
        {
            info_byte += 1;
        }
        if (datatype_nr != ow.modal.datatype)
        {
            info_byte += 2;
        }
        if (pointarray[0].X != ow.modal.geometry_x)
        {
            info_byte += 16;
        }
        if (pointarray[0].Y != ow.modal.geometry_y)
        {
            info_byte += 8;
        }
        ow.writeUnsignedInteger(oasValues.POLYGON);
        ow.writeRaw(info_byte);
        switch (info_byte & 1)
        {
            case > 0:
                ow.modal.layer = layer_nr;
                ow.writeUnsignedInteger((uint)ow.modal.layer);
                break;
        }
        switch (info_byte & 2)
        {
            case > 0:
                ow.modal.datatype = datatype_nr;
                ow.writeUnsignedInteger((uint)datatype_nr);
                break;
        }
        ow.writePointArray(pointarray, false);
        switch (info_byte & 16)
        {
            case > 0:
                ow.modal.geometry_x = (int)pointarray[0].X;
                ow.writeSignedInteger(ow.modal.geometry_x);
                break;
        }
        switch (info_byte & 8)
        {
            case > 0:
                ow.modal.geometry_y = (int)pointarray[0].Y;
                ow.writeSignedInteger(ow.modal.geometry_y);
                break;
        }


    }
    public override bool isText()
    {
        return pIsText();
    }

    private bool pIsText()
    {
        return text;
    }

    public override string getName()
    {
        return pGetName();
    }

    private string pGetName()
    {
        return name;
    }

    public override bool isPolygon()
    {
        return true;
    }
}