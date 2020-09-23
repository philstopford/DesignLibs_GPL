using gds;
using geoLib;
using geoWrangler;
using oasis;
using System;
using System.Linq;

namespace geoCoreLib
{
    public class GCPolygon : GCElement
    {
        public GeoLibPoint[] pointarray { get; set; }

        public GCPolygon()
        {
            pGCPolygon();
        }

        void pGCPolygon()
        {
            pointarray = new GeoLibPoint[0];
        }

        public GCPolygon(GCPolygon source)
        {
            pGCPolygon(source);
        }

        void pGCPolygon(GCPolygon source)
        {
            pointarray = source.pointarray.ToArray();
            layer_nr = source.layer_nr;
            datatype_nr = source.datatype_nr;
        }

        public GCPolygon(GeoLibPoint[] points, Int32 layer, Int32 datatype)
        {
            pGCPolygon(points, layer, datatype);
        }

        void pGCPolygon(GeoLibPoint[] points, Int32 layer, Int32 datatype)
        {
            pointarray = points.ToArray();
            layer_nr = layer;
            datatype_nr = datatype;
        }

        public override void rotate(double angleDegree, GeoLibPoint pos)
        {
            pRotate(angleDegree, pos);
        }

        void pRotate(double angleDegree, GeoLibPoint pos)
        {
            if (angleDegree == 0)
            {
                return;
            }

            for (Int32 i = 0; i < pointarray.Length; i++)
            {
                pointarray[i] = pRotate(pointarray[i], pos, angleDegree);
            }
        }

        GeoLibPoint pRotate(GeoLibPoint point, GeoLibPoint pivot, double angleDegree)
        {
            if (angleDegree == 0)
            {
                return point;
            }

            double angle = angleDegree * Math.PI / 180;
            double x = pivot.X + ((point.X - pivot.X) * Math.Cos(angle) -
                (point.Y - pivot.Y) * Math.Sin(angle));
            double y = pivot.Y + ((point.X - pivot.X) * Math.Sin(angle) +
                (point.Y - pivot.Y) * Math.Cos(angle));

            return new GeoLibPoint(x, y);
        }

        public override void scale(double size)
        {
            pScale(size);
        }

        void pScale(double size)
        {
            for (int i = 0; i < pointarray.Length; i++)
            {
                pointarray[i] = new GeoLibPoint(pointarray[i].X * size, pointarray[i].Y * size);
            }
        }

        public override void scale(GeoLibPoint origin, double size)
        {
            pScale(origin, size);
        }

        void pScale(GeoLibPoint origin, double size)
        {
            for (int i = 0; i < pointarray.Length; i++)
            {
                pointarray[i] = new GeoLibPoint(origin.X + ((pointarray[i].X - origin.X) * size), origin.Y + ((pointarray[i].Y - origin.Y) * size));
            }
        }

        public void deletePoint(Int32 pos)
        {
            pDeletePoint(pos);
        }

        void pDeletePoint(Int32 pos)
        {
            GeoLibPoint[] newPointArray = new GeoLibPoint[pointarray.Length - 1];
            for (Int32 i = pos; i < pointarray.Length - 1; i++)
            {
                newPointArray[i] = pointarray[i + 1];
            }
            pointarray = newPointArray;
        }

        public void addPoint(Int32 pos)
        {
            pAddPoint(pos);
        }

        void pAddPoint(Int32 pos)
        {
            GeoLibPoint[] newPointArray = new GeoLibPoint[pointarray.Length + 1];
            for (Int32 pt = newPointArray.Length - 1; pt > pos; pt--)
            {
                newPointArray[pt] = pointarray[pt - 1];
            }
            pointarray = newPointArray;
        }

        public override void mapSelect(GCStrans m)
        {
            pMapSelect(m);
        }

        void pMapSelect(GCStrans m)
        {
            if (select)
            {
                m.matrix.TransformPoints(pointarray);
            }
            /*
			if (select)
			{
				for (i = 0; i < pointarray.Count(); i++)
				{
					Point p = pointarray[i];
					p = m.matrix.map(p);
					pointarray[i] = p;
				}
			}
			else
			{
				for (i = 0; i < pointarray.size(); i++)
				{
					if (p_select.testBit(i))
					{
						p = pointarray[i);
						p = m.matrix.TransformPoints(p);
						pointarray.setPoint(i, p);
					}
				}
			}
			*/
        }

        public override void map(GCStrans m)
        {
            pMap(m);
        }

        void pMap(GCStrans m)
        {
            m.matrix.TransformPoints(pointarray);
            /*
			for (Int32 i = 0; i < pointarray.size(); i++)
			{
				Point p = pointarray[i];
				p = m.matrix.map(p);
				pointarray.setPoint(i, p);
			}
			*/
        }

        public override void minimum(GeoLibPoint pos)
        {
            pMinimum(pos);
        }

        void pMinimum(GeoLibPoint pos)
        {
            for (Int32 i = 0; i < pointarray.Count(); i++)
            {
                GeoLibPoint p = pointarray[i];
                if (p.X < pos.X)
                {
                    pos.X = p.X;
                }
                if (p.Y < pos.Y)
                {
                    pos.Y = p.Y;
                }
            }
        }

        public override void maximum(GeoLibPoint pos)
        {
            pMaximum(pos);
        }

        void pMaximum(GeoLibPoint pos)
        {
            for (Int32 i = 0; i < pointarray.Count(); i++)
            {
                GeoLibPoint p = pointarray[i];
                if (p.X > pos.X)
                {
                    pos.X = p.X;
                }
                if (p.Y > pos.Y)
                {
                    pos.Y = p.Y;
                }
            }
        }

        public override void moveSelect(GeoLibPoint pos)
        {
            pMoveSelect(pos);
        }

        void pMoveSelect(GeoLibPoint pos)
        {
            if (select)
            {
                for (Int32 i = 0; i < pointarray.Count(); i++)
                {
                    pointarray[i].Offset(pos);
                }
            }
        }

        public override void move(GeoLibPoint pos)
        {
            pMove(pos);
        }

        void pMove(GeoLibPoint pos)
        {
            for (Int32 i = 0; i < pointarray.Count(); i++)
            {
                pointarray[i] = new GeoLibPoint(pointarray[i].X + pos.X, pointarray[i].Y + pos.Y);
            }
        }

        public override GCPolygon convertToPolygon()
        {
            return pConvertToPolygon();
        }

        GCPolygon pConvertToPolygon()
        {
            return new GCPolygon(this);
        }

        public override void clean()
        {
            pClean();
        }

        void pClean()
        {
            GeoLibPoint p = new GeoLibPoint(0, 0);
            double a;
            bool b;
            int anz = 0;
            for (a = 0; ((a < 350) || (a > 370)) && (anz < 3);)
            {
                a = 0;
                for (int i = 0; i < pointarray.Length - 1; i++)
                {
                    if (pointarray.Length < 4)
                    {
                        return; //no area
                    }
                    if (pointarray[i] == pointarray[i + 1])
                    {
                        deletePoint(i + 1);
                        if (pointarray.Length == (i + 1))
                        {
                            pointarray[0] = pointarray[i];
                        }
                        i--;
                    }
                    if (i < pointarray.Length - 2)
                    {
                        if (nearlyParallel(pointarray[i], pointarray[i + 1], pointarray[i + 1], pointarray[i + 2]))
                        {
                            deletePoint(i + 1);
                            i = 0;
                        }
                    }
                    if (pointarray.Length > 3)
                    {
                        while ((pointarray.Length > 3) && (parallel(pointarray[0], pointarray[1], pointarray[pointarray.Length - 1], pointarray[pointarray.Length - 2])))
                        {
                            deletePoint(pointarray.Length - 1);
                            pointarray[0] = pointarray[pointarray.Length - 1];
                        }
                    }
                }
                // sort points
                for (Int32 i = 0; i < pointarray.Length - 2; i++)
                {
                    a += angle(pointarray[i], pointarray[i + 1], pointarray[i + 2]);
                }
                if (pointarray.Length > 3)
                {
                    a += angle(pointarray[pointarray.Length - 2], pointarray[0], pointarray[1]);
                }
                // Punkte linksdrehend
                if (a < -185)
                {
                    for (Int32 i = 1; i < (pointarray.Length) / 2; i++)
                    {
                        p = pointarray[pointarray.Length - i - 1];
                        pointarray[pointarray.Length - i - 1] = pointarray[i];
                        pointarray[i] = p;
                    }
                }
                // (simple selfintersecting)
                if ((a > 370) || ((a < 350) && (a > -350)))
                {
                    for (Int32 i = 0; i < pointarray.Length - 1; i++)
                    {
                        for (Int32 j = i + 2; j < pointarray.Length - 1; j++)
                        {
                            b = cutPoint2(pointarray[i], pointarray[i + 1], pointarray[j], pointarray[j + 1], p);
                            if (b)
                            {
                                if (!((i == 0) && (j == pointarray.Length - 2)))
                                {
                                    addPoint(j + 1);
                                    j++;
                                    pointarray[j] = p;
                                    addPoint(i + 1);
                                    i++;
                                    pointarray[i] = p;
                                    for (Int32 k = i + 1; k <= (i + j) / 2; k++)
                                    {
                                        p = pointarray[j - k + i + 1];
                                        pointarray[j - k + i + 1] = pointarray[k];
                                        pointarray[k] = p;
                                    }
                                    j = pointarray.Length;
                                    i = pointarray.Length;
                                }
                            }
                        }
                    }
                }
                // remove selfintersection 
                if (pointarray.Length >= 8)
                {
                    bool ende2 = false;
                    for (int i = 0; i < pointarray.Length - 1; i++)
                    {
                        bool ende1 = false;
                        for (int j = i + 2; j < pointarray.Length - 1; j++)
                        {
                            if (identical(pointarray[i], pointarray[i + 1], pointarray[j + 1], pointarray[j]))
                            {
                                int h1, h2, h3, h4, h5, h6;
                                bool side = true;
                                bool change = false;
                                if (pointarray[j + 1] == pointarray[i])
                                {
                                    change = true;
                                    h1 = i - 1;
                                    if (h1 == (-1))
                                    {
                                        h1 = pointarray.Length - 2;
                                    }
                                    h2 = j + 2;
                                    if (h2 == pointarray.Length)
                                    {
                                        h2 = 1;
                                    }
                                    double d1, d2;
                                    d1 = angle(pointarray[h1], pointarray[i], pointarray[i + 1]);
                                    d2 = angle(pointarray[h2], pointarray[j + 1], pointarray[j]);
                                    if (d1 < 0)
                                    {
                                        d1 += 360;
                                    }
                                    if (d2 < 0)
                                    {
                                        d2 += 360;
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
                                    if (h2 == pointarray.Length)
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
                                    h1 = i - 1;
                                    if (h1 == -1)
                                    {
                                        h1 = pointarray.Length - 2;
                                    }
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
                                    if (j == (-1))
                                    {
                                        j = pointarray.Length - 2; ende1 = true;
                                    }
                                    if (i == pointarray.Length)
                                    {
                                        i = 1; ende2 = true;
                                    }
                                }
                                h2 = i + 1;
                                if (h2 == pointarray.Length)
                                {
                                    h2 = 1;
                                }
                                h3 = h2 + 1;
                                if (h3 == pointarray.Length)
                                {
                                    h3 = 1;
                                }
                                h1 = j - 1;
                                if (h1 == -1)
                                {
                                    h1 = pointarray.Length - 2;
                                }
                                h4 = j + 1;
                                if (h4 == pointarray.Length)
                                {
                                    h4 = 1;
                                }
                                h5 = h4 + 1;
                                if (h5 == pointarray.Length)
                                {
                                    h5 = 1;
                                }
                                h6 = i - 1;
                                if (h6 == -1)
                                {
                                    h6 = pointarray.Length - 2;
                                }
                                if (onLine2(pointarray[i], pointarray[h2], pointarray[j]))
                                {
                                    if ((distance(pointarray[i], pointarray[h2], pointarray[h1]) < 0) != side)
                                    {
                                        change = false;
                                    }
                                }
                                else if (onLine2(pointarray[j], pointarray[h4], pointarray[h2]))
                                {
                                    if ((distance(pointarray[j], pointarray[h4], pointarray[h3]) < 0) != side)
                                    {
                                        change = false;
                                    }
                                }
                                else
                                {
                                    double d1, d2;
                                    d1 = angle(pointarray[h6], pointarray[i], pointarray[h2]);
                                    d2 = angle(pointarray[h5], pointarray[h4], pointarray[j]);
                                    if (d1 < 0)
                                    {
                                        d1 += 360;
                                    }
                                    if (d2 < 0)
                                    {
                                        d2 += 360;
                                    }
                                    if ((d1 < d2) != side)
                                    {
                                        change = false;
                                    }
                                }
                                if (change)
                                {
                                    for (Int32 k = i + 1; k < (i + j + 2) / 2; k++)
                                    {
                                        p = pointarray[j + i - k + 1];
                                        pointarray[j + i - k + 1] = pointarray[k];
                                        pointarray[k] = p;
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
                }
                anz++;
            }
        }

        public override void saveGDS(gdsWriter gw)
        {
            pSaveGDS(gw);
        }

        void pSaveGDS(gdsWriter gw)
        {
            GeoLibPoint pc = new GeoLibPoint();
            int r = 0;
            if (isCircle(ref pc, ref r))
            {
                GCPath tPath = new GCPath(new GeoLibPoint[] { pc }, layer_nr, datatype_nr);
                tPath.setWidth(r * 2);
                tPath.setCap(1); // set cap to 'round' type. Also forced in the saver for a round path type, as a safety.
                tPath.setRound(true); // this is the important thing to set - it forces the correct handling (including cap value)
                tPath.saveGDS(gw);
            }
            else
            {
                // polygon (boundary)
                gw.bw.Write((UInt16)4);
                gw.bw.Write((byte)8);
                gw.bw.Write((byte)0);
                //layer
                gw.bw.Write((UInt16)6);
                gw.bw.Write((byte)0x0D);
                gw.bw.Write((byte)2);
                gw.bw.Write((Int16)layer_nr);
                //datatype
                gw.bw.Write((UInt16)6);
                gw.bw.Write((byte)0x0E);
                gw.bw.Write((byte)2);
                gw.bw.Write((Int16)(datatype_nr));
                int i = pointarray.Length;
                if (i > 8191)
                {
                    i = 8191;
                    // "Polygon with more than 8191 points. Data is lost."
                }
                //xy 
                // Add one to the pointlist length value (i) to mark the closed shape.
                int val = ((i + 1) * 2 * 4) + 4;
                gw.bw.Write((UInt16)val);
                gw.bw.Write((byte)0x10);
                gw.bw.Write((byte)3);

                for (int k = 0; k < i; k++)
                {
                    gw.bw.Write(pointarray[k].X);
                    gw.bw.Write(pointarray[k].Y);
                }

                // close the polygon.
                gw.bw.Write(pointarray[0].X);
                gw.bw.Write(pointarray[0].Y);

                // endel
                gw.bw.Write((UInt16)4);
                gw.bw.Write((byte)0x11);
                gw.bw.Write((byte)0);
            }
        }

        public override void saveOASIS(oasWriter ow)
        {
            pSaveOASIS(ow);
        }

        bool isCircle(ref GeoLibPoint p, ref int radius)
        {
            if (pointarray.Length < 10)
            {
                return false;
            }

            GeoLibPointF mid = GeoWrangler.midPoint(pointarray);

            double min_distance = Math.Abs(GeoWrangler.distanceBetweenPoints(mid, new GeoLibPointF(pointarray[0])));
            double max_distance = Math.Abs(GeoWrangler.distanceBetweenPoints(mid, new GeoLibPointF(pointarray[0])));


            for (int i = 1; i < pointarray.Length; i++)
            {
                double distance = Math.Abs(GeoWrangler.distanceBetweenPoints(mid, new GeoLibPointF(pointarray[i])));
                min_distance = Math.Min(min_distance, distance);
                max_distance = Math.Max(max_distance, distance);
            }

            double delta = max_distance - min_distance;

            // Tolerance value - one DB unit.
            double tol = 1.0f * 1000;

            if (delta < tol)
            {
                return true;
            }

            return false;

            long x, y;
            x = 0; y = 0;
            for (int i = 1; i < pointarray.Length; i++)
            {
                x += pointarray[i].X;
                y += pointarray[i].Y;
            }
            x = round((double)x / (pointarray.Length - 1));
            y = round((double)y / (pointarray.Length - 1));
            GeoLibPoint pc = new GeoLibPoint(x, y);
            double d = distance(pc, pointarray[1]);
            double sum = d;
            double min = d;
            double max = d;
            int offset = (int)(0.02 * distance(pc, pointarray[0])) + 10;
            for (int i = 2; i < pointarray.Length; i++)
            {
                double d1 = distance(pc, pointarray[i]);
                sum += d1;
                if (d1 > max)
                {
                    max = d1;
                    if (max > d + offset)
                    {
                        return false;
                    }
                }
                if (d1 < min)
                {
                    min = d1;
                    if (min < d - offset)
                    {
                        return false;
                    }
                }
            }
            p = pc;
            radius = round(sum / (pointarray.Length - 1));
            offset = (int)(0.001 * radius + 3);
            if (radius < round(max - offset))
            {
                return false;
            }
            if (radius > round(min + offset))
            {
                return false;
            }
            return true;
        }

        void pSaveOASIS(oasWriter ow)
        {
            byte info_byte = 0;
            if (GCSetup.oasisSaveCircle)
            {
                // pc and r are manipulated in isCircle()
                GeoLibPoint pc = new GeoLibPoint();
                int r = 0;
                if (isCircle(ref pc, ref r))
                {
                    // save as circle
                    if (!ow.modal.absoluteMode)
                    {
                        ow.modal.absoluteMode = true;
                    }
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
                    if (r != ow.modal.circle_radius)
                    {
                        info_byte += 32;
                    }
                    ow.writeUnsignedInteger(27);
                    ow.writeRaw(info_byte);
                    if ((info_byte & 1) > 0)
                    {
                        ow.modal.layer = layer_nr;
                        ow.writeUnsignedInteger((uint)ow.modal.layer);
                    }
                    if ((info_byte & 2) > 0)
                    {
                        ow.modal.datatype = datatype_nr;
                        ow.writeUnsignedInteger((uint)datatype_nr);
                    }
                    if ((info_byte & 32) > 0)
                    {
                        ow.modal.circle_radius = r;
                        ow.writeUnsignedInteger((uint)ow.modal.circle_radius);
                    }
                    if ((info_byte & 16) > 0)
                    {
                        ow.modal.geometry_x = pc.X;
                        ow.writeSignedInteger(ow.modal.geometry_x);
                    }
                    if ((info_byte & 8) > 0)
                    {
                        ow.modal.geometry_y = pc.Y;
                        ow.writeSignedInteger(ow.modal.geometry_y);
                    }
                    return;

                }
            }
            // check if ctrapezoid
            if ((GCSetup.oasisSaveCtrapezoid) && ((pointarray.Length == 4) || (pointarray.Length == 5)))
            {
                int form = 0;
                int off = 0;
                int w, h;
                for (int i = 1; i < pointarray.Length; i++)
                {
                    GeoLibPointF pd = GeoWrangler.distanceBetweenPoints_point(pointarray[i], pointarray[i - 1]);
                    if ((pd.X == pd.Y) && (pd.X > 0))
                    {
                        form = form * 10 + 9;
                    }
                    else if ((pd.X == pd.Y) && (pd.X < 0))
                    {
                        form = form * 10 + 1;
                    }
                    else if ((pd.X == -pd.Y) && (pd.X < 0))
                    {
                        form = form * 10 + 7;
                    }
                    else if ((pd.X == -pd.Y) && (pd.X > 0))
                    {
                        form = form * 10 + 3;
                    }
                    else if ((pd.X == 0) && (pd.Y > 0))
                    {
                        form = form * 10 + 8;
                    }
                    else if ((pd.X == 0) && (pd.Y < 0))
                    {
                        form = form * 10 + 2;
                    }
                    else if ((pd.Y == 0) && (pd.X > 0))
                    {
                        form = form * 10 + 6;
                    }
                    else if ((pd.Y == 0) && (pd.X < 0))
                    {
                        form = form * 10 + 4;
                    }
                    else
                    {
                        form = -1;
                    }
                    if (form == -1)
                    {
                        break;
                    }
                }
                if (form != -1)
                {
                    if (pointarray.Length == 4)
                    {
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
                            case 672: ow.writeCtrapezoid(layer_nr, 16, pointarray[0 + off].X, pointarray[0 + off].Y, pointarray[1 + off].X - pointarray[0 + off].X, 0, datatype_nr); return;
                            case 942: ow.writeCtrapezoid(layer_nr, 17, pointarray[0 + off].X, pointarray[0 + off].Y, pointarray[1 + off].X - pointarray[0 + off].X, 0, datatype_nr); return;
                            case 681: ow.writeCtrapezoid(layer_nr, 18, pointarray[0 + off].X, pointarray[0 + off].Y, pointarray[1 + off].X - pointarray[0 + off].X, 0, datatype_nr); return;
                            case 438: w = pointarray[0 + off].X - pointarray[1 + off].X; ow.writeCtrapezoid(layer_nr, 19, pointarray[0 + off].X - w, pointarray[0 + off].Y - w, w, 0, datatype_nr); return;
                            case 671: h = (pointarray[1 + off].X - pointarray[0 + off].X) / 2; ow.writeCtrapezoid(layer_nr, 20, pointarray[0 + off].X, pointarray[0 + off].Y, 0, h, datatype_nr); return;
                            case 943: h = pointarray[1 + off].Y - pointarray[0 + off].Y; ow.writeCtrapezoid(layer_nr, 21, pointarray[0 + off].X - h, pointarray[0 + off].Y, 0, h, datatype_nr); return;
                            case 972: w = pointarray[1 + off].X - pointarray[0 + off].X; ow.writeCtrapezoid(layer_nr, 22, pointarray[0 + off].X, pointarray[0 + off].Y, w, 0, datatype_nr); return;
                            case 813: w = (pointarray[1 + off].Y - pointarray[0 + off].Y) / 2; ow.writeCtrapezoid(layer_nr, 23, pointarray[0 + off].X - w, pointarray[0 + off].Y, w, 0, datatype_nr); return;
                        }
                    }
                    else if (pointarray.Length == 5)
                    {
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
                                w = pointarray[1 + off].X - pointarray[0 + off].X;
                                if (off == 3)
                                {
                                    h = pointarray[1].Y - pointarray[2].Y;
                                }
                                else
                                {
                                    h = pointarray[2 + off].Y - pointarray[1 + off].Y;
                                }
                                if (h == w)
                                {
                                    ow.writeCtrapezoid(layer_nr, 25, pointarray[0 + off].X, pointarray[0 + off].Y, w, 0, datatype_nr);
                                }
                                else
                                {
                                    ow.writeCtrapezoid(layer_nr, 24, pointarray[0 + off].X, pointarray[0 + off].Y, w, h, datatype_nr);
                                }
                                return;
                            case 6742:
                                w = pointarray[1 + off].X - pointarray[0 + off].X;
                                if (off == 0)
                                {
                                    h = pointarray[2].Y - pointarray[1].Y;
                                }
                                else
                                {
                                    h = pointarray[-1 + off].Y - pointarray[0 + off].Y;
                                }
                                ow.writeCtrapezoid(layer_nr, 0, pointarray[0 + off].X, pointarray[0 + off].Y, w, h, datatype_nr);
                                return;
                            case 6942:
                                if (off == 0)
                                {
                                    h = pointarray[2].Y - pointarray[1].Y;
                                }
                                else
                                {
                                    h = pointarray[-1 + off].Y - pointarray[0 + off].Y;
                                }
                                w = pointarray[1 + off].X - pointarray[0 + off].X + h;
                                ow.writeCtrapezoid(layer_nr, 1, pointarray[0 + off].X, pointarray[0 + off].Y, w, h, datatype_nr);
                                return;
                            case 6841:
                                w = pointarray[1 + off].X - pointarray[0 + off].X;
                                if (off == 0)
                                {
                                    h = pointarray[2].Y - pointarray[1].Y;
                                }
                                else
                                {
                                    h = pointarray[-1 + off].Y - pointarray[0 + off].Y;
                                }
                                ow.writeCtrapezoid(layer_nr, 2, pointarray[0 + off].X, pointarray[0 + off].Y, w, h, datatype_nr);
                                return;
                            case 6843:
                                if (off == 0)
                                {
                                    h = pointarray[2].Y - pointarray[1].Y;
                                }
                                else
                                {
                                    h = pointarray[-1 + off].Y - pointarray[0 + off].Y;
                                }
                                w = pointarray[1 + off].X - pointarray[0 + off].X + h;
                                ow.writeCtrapezoid(layer_nr, 3, pointarray[0 + off].X - h, pointarray[0 + off].Y, w, h, datatype_nr);
                                return;
                            case 6741:
                                w = pointarray[1 + off].X - pointarray[0 + off].X;
                                if (off == 0)
                                {
                                    h = pointarray[2].Y - pointarray[1].Y;
                                }
                                else
                                {
                                    h = pointarray[-1 + off].Y - pointarray[0 + off].Y;
                                }
                                ow.writeCtrapezoid(layer_nr, 4, pointarray[0 + off].X, pointarray[0 + off].Y, w, h, datatype_nr);
                                return;
                            case 6943:
                                if (off == 0)
                                {
                                    h = pointarray[2].Y - pointarray[1].Y;
                                }
                                else
                                {
                                    h = pointarray[-1 + off].Y - pointarray[0 + off].Y;
                                }
                                w = pointarray[1 + off].X - pointarray[0 + off].X + h + h;
                                ow.writeCtrapezoid(layer_nr, 5, pointarray[0 + off].X - h, pointarray[0 + off].Y, w, h, datatype_nr);
                                return;
                            case 6941:
                                if (off == 0)
                                {
                                    h = pointarray[2].Y - pointarray[1].Y;
                                }
                                else
                                {
                                    h = pointarray[-1 + off].Y - pointarray[0 + off].Y;
                                }
                                w = pointarray[1 + off].X - pointarray[0 + off].X + h;
                                ow.writeCtrapezoid(layer_nr, 6, pointarray[0 + off].X, pointarray[0 + off].Y, w, h, datatype_nr);
                                return;
                            case 6743:
                                if (off == 0)
                                {
                                    h = pointarray[2].Y - pointarray[1].Y;
                                }
                                else
                                {
                                    h = pointarray[-1 + off].Y - pointarray[0 + off].Y;
                                }
                                w = pointarray[1 + off].X - pointarray[0 + off].X + h;
                                ow.writeCtrapezoid(layer_nr, 7, pointarray[0 + off].X - h, pointarray[0 + off].Y, w, h, datatype_nr);
                                return;
                            case 6872:
                                w = pointarray[1 + off].X - pointarray[0 + off].X;
                                if (off == 0)
                                {
                                    h = pointarray[3].Y - pointarray[0].Y;
                                }
                                else
                                {
                                    h = pointarray[-1 + off].Y - pointarray[0 + off].Y;
                                }
                                ow.writeCtrapezoid(layer_nr, 8, pointarray[0 + off].X, pointarray[0 + off].Y, w, h, datatype_nr);
                                return;
                            case 6812:
                                w = pointarray[1 + off].X - pointarray[0 + off].X;
                                if (off == 0)
                                {
                                    h = pointarray[2].Y - pointarray[0].Y;
                                }
                                else
                                {
                                    h = pointarray[-1 + off].Y - pointarray[0 + off].Y + w;
                                }
                                ow.writeCtrapezoid(layer_nr, 9, pointarray[0 + off].X, pointarray[0 + off].Y, w, h, datatype_nr);
                                return;
                            case 9842:
                                w = pointarray[1 + off].X - pointarray[0 + off].X;
                                if (off == 0)
                                {
                                    h = pointarray[3].Y - pointarray[0].Y;
                                }
                                else
                                {
                                    h = pointarray[-1 + off].Y - pointarray[0 + off].Y;
                                }
                                ow.writeCtrapezoid(layer_nr, 10, pointarray[0 + off].X, pointarray[0 + off].Y, w, h, datatype_nr);
                                return;
                            case 8423:
                                h = pointarray[1 + off].Y - pointarray[0 + off].Y;
                                if (off == 0)
                                {
                                    w = pointarray[0].X - pointarray[2].X;
                                }
                                else
                                {
                                    w = pointarray[0 + off].X - pointarray[-1 + off].X;
                                }
                                ow.writeCtrapezoid(layer_nr, 11, pointarray[0 + off].X - w, pointarray[0 + off].Y, w, h, datatype_nr);
                                return;
                            case 9872:
                                w = pointarray[1 + off].X - pointarray[0 + off].X;
                                if (off == 0)
                                {
                                    h = pointarray[3].Y - pointarray[0].Y;
                                }
                                else
                                {
                                    h = pointarray[-1 + off].Y - pointarray[0 + off].Y;
                                }
                                ow.writeCtrapezoid(layer_nr, 12, pointarray[0 + off].X, pointarray[0 + off].Y, w, h, datatype_nr);
                                return;
                            case 8123:
                                h = pointarray[1 + off].Y - pointarray[0 + off].Y;
                                if (off == 0)
                                {
                                    w = pointarray[0].X - pointarray[2].X;
                                }
                                else
                                {
                                    w = pointarray[0 + off].X - pointarray[-1 + off].X;
                                }
                                ow.writeCtrapezoid(layer_nr, 13, pointarray[0 + off].X - w, pointarray[0 + off].Y, w, h, datatype_nr);
                                return;
                            case 9812:
                                w = pointarray[1 + off].X - pointarray[0 + off].X;
                                if (off == 0)
                                {
                                    h = pointarray[2].Y - pointarray[0].Y;
                                }
                                else
                                {
                                    h = pointarray[-1 + off].Y - pointarray[0 + off].Y + w;
                                }
                                ow.writeCtrapezoid(layer_nr, 14, pointarray[0 + off].X, pointarray[0 + off].Y, w, h, datatype_nr);
                                return;
                            case 8723:
                                if (off == 0)
                                {
                                    w = pointarray[0].X - pointarray[2].X;
                                }
                                else
                                {
                                    w = pointarray[0 + off].X - pointarray[-1 + off].X;
                                }
                                h = pointarray[1 + off].Y - pointarray[0 + off].Y + w;
                                ow.writeCtrapezoid(layer_nr, 15, pointarray[0 + off].X - w, pointarray[0 + off].Y, w, h, datatype_nr);
                                return;
                        }
                    }
                }
            }
            //trapezoid
            if ((GCSetup.oasisSaveTrapezoid) && (pointarray.Length == 5))
            {
                int form = 0;
                int minX, minY, maxX, maxY;
                int w, h, da, db;
                minX = pointarray[0].X;
                maxX = pointarray[0].X;
                minY = pointarray[0].Y;
                maxY = pointarray[0].Y;
                for (int i = 1; i < pointarray.Length; i++)
                {
                    if (pointarray[i].X > maxX)
                    {
                        maxX = pointarray[i].X;
                    }
                    if (pointarray[i].Y > maxY)
                    {
                        maxY = pointarray[i].Y;
                    }
                    if (pointarray[i].X < minX)
                    {
                        minX = pointarray[i].X;
                    }
                    if (pointarray[i].Y < minY)
                    {
                        minY = pointarray[i].Y;
                    }
                    GeoLibPointF pd = GeoWrangler.distanceBetweenPoints_point(pointarray[i], pointarray[i - 1]);
                    if (pd.X == 0)
                    {
                        form = form * 10 + 1;
                    }
                    else if (pd.Y == 0)
                    {
                        form = form * 10 + 2;
                    }
                    else
                    {
                        form = form * 10;
                    }
                }
                w = maxX - minX;
                h = maxY - minY;
                switch (form)
                {
                    case 1010:
                    case 1210:
                    case 1012:
                        if (pointarray[0].X > pointarray[2].X)
                        {
                            da = pointarray[0].Y - pointarray[3].Y;
                            db = pointarray[1].Y - pointarray[2].Y;
                        }
                        else
                        {
                            db = pointarray[0].Y - pointarray[3].Y;
                            da = pointarray[1].Y - pointarray[2].Y;
                        }
                        ow.writeTrapezoid(layer_nr, 1, minX, minY, w, h, -da, -db, datatype_nr);
                        return;
                    case 101:
                    case 121:
                    case 2101:
                        if (pointarray[1].X > pointarray[3].X)
                        {
                            da = pointarray[1].Y - pointarray[4].Y;
                            db = pointarray[2].Y - pointarray[3].Y;
                        }
                        else
                        {
                            db = pointarray[1].Y - pointarray[4].Y;
                            da = pointarray[2].Y - pointarray[3].Y;
                        }
                        ow.writeTrapezoid(layer_nr, 1, minX, minY, w, h, -da, -db, datatype_nr);
                        //printf("minX %d, maxX %d, minY %d, maxY %d, w %d, h %d, da %d, db %d\n",minX,maxX,minY,maxY,w,h,da,db);
                        return;
                    case 2020:
                    case 2120:
                    case 2021:
                        if (pointarray[0].Y < pointarray[2].Y)
                        {
                            da = pointarray[0].X - pointarray[3].X;
                            db = pointarray[1].X - pointarray[2].X;
                        }
                        else
                        {
                            db = -pointarray[0].X + pointarray[3].X;
                            da = -pointarray[1].X + pointarray[2].X;
                        }
                        ow.writeTrapezoid(layer_nr, 0, minX, minY, w, h, -da, -db, datatype_nr);
                        return;
                    case 202:
                    case 1202:
                    case 212:
                        if (pointarray[1].Y < pointarray[3].Y)
                        {
                            da = pointarray[1].X - pointarray[4].X;
                            db = pointarray[2].X - pointarray[3].X;
                        }
                        else
                        {
                            db = -pointarray[1].X + pointarray[4].X;
                            da = -pointarray[2].X + pointarray[3].X;
                        }
                        ow.writeTrapezoid(layer_nr, 0, minX, minY, w, h, -da, -db, datatype_nr);
                        return;
                }
            }
            // save as polygon
            if (!ow.modal.absoluteMode)
            {
                ow.modal.absoluteMode = true;
            }
            info_byte = 32;  //write pointlist;
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
            ow.writeUnsignedInteger(21);
            ow.writeRaw(info_byte);
            if ((info_byte & 1) > 0)
            {
                ow.modal.layer = layer_nr;
                ow.writeUnsignedInteger((uint)ow.modal.layer);
            }
            if ((info_byte & 2) > 0)
            {
                ow.modal.datatype = datatype_nr;
                ow.writeUnsignedInteger((uint)datatype_nr);
            }
            ow.writePointArray(pointarray, false);
            if ((info_byte & 16) > 0)
            {
                ow.modal.geometry_x = pointarray[0].X;
                ow.writeSignedInteger(ow.modal.geometry_x);
            }
            if ((info_byte & 8) > 0)
            {
                ow.modal.geometry_y = pointarray[0].Y;
                ow.writeSignedInteger(ow.modal.geometry_y);
            }


        }
    }
}
