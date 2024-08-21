using gds;
using oasis;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Clipper2Lib;
using geoWrangler;

namespace geoCoreLib;

public class GCCellRefArray : GCElement
{
    private Point64 point;
    public GCCell cell_ref { get; set; }
    public GCStrans trans { get; set; }
    
    public override string getName()
    {
        return cell_ref.getName();
    }
    
    // Used for OASIS
    public Repetition repetition;

    public GCCellRefArray(GCCell c, double angle, double mag, bool mirror_x, Repetition r)
    {
        repetition = new(r);
        cell_ref = c;
        trans = new();
        trans.reset();
        trans.angle = angle;
        trans.mag = mag;
        trans.mirror_x = mirror_x;
    }
    
    public GCCellRefArray(GCCell c, Path64 array, int xCount, int yCount)
    {
        repetition = new Repetition();
        cell_ref = c;
        point = array[0];
        // Shims.
        repetition.columns = xCount;
        repetition.rows = yCount;
        if (xCount < 1)
        {
            repetition.columns = 1;
        }
        if (yCount < 1)
        {
            repetition.rows = 1;
        }
        repetition.rowVector = new ((array[1].X - array[0].X) / repetition.rows, (array[1].Y - array[0].Y) / repetition.rows);
        repetition.colVector = new((array[2].X - array[0].X) / repetition.columns, (array[2].Y - array[0].Y) / repetition.columns);
        repetition.type = Repetition.RepetitionType.Regular;
        if (((repetition.rowVector.X != 0) && (repetition.rowVector.Y != 0)) || ((repetition.colVector.X != 0) && (repetition.colVector.Y != 0)))
        {
            repetition.type = Repetition.RepetitionType.Rectangular;
        }
        // Tag layer and datatype to allow this element to be filtered out from LD and geo lists.
        layer_nr = -1;
        datatype_nr = -1;
        trans = new GCStrans();
        trans.reset();
    }

    public GCCellRefArray(GCCell c, Path64 explicitArray)
    {
        if (explicitArray.Count < 1)
        {
            throw new Exception("No valid array provided");
        }

        repetition = new();
        trans = new();
        trans.reset();

        cell_ref = c;
        point = explicitArray[0];
        
        // Figure out the nature of our repetition in the array.
        Rect64 bounds = Clipper.GetBounds(explicitArray);
        
        bool explicitY = bounds.Width == 0;
        bool explicitX = bounds.Height == 0;
        bool explicitXY = !explicitX && !explicitY;

        if (explicitXY)
        {
            repetition.type = Repetition.RepetitionType.Explicit;
        }
        else
        {
            if (explicitX)
            {
                repetition.type = Repetition.RepetitionType.ExplicitX;
            }

            if (explicitY)
            {
                repetition.type = Repetition.RepetitionType.ExplicitY;
            }
        }

        repetition.offsets = new(explicitArray);
    }

    public GCCellRefArray(GCCell c, Point64 pos1, Point64 pos2, int xCount, int yCount)
    {
        repetition = new();
        cell_ref = c;
        point = pos1;
        repetition.rowVector = new(pos1);
        repetition.colVector = new(pos2);
        repetition.columns = xCount;
        repetition.rows = yCount;
        // Tag layer and datatype to allow this element to be filtered out from LD and geo lists.
        layer_nr = -1;
        datatype_nr = -1;
        trans = new GCStrans();
        trans.reset();
    }

    public GCCellRefArray()
    {
        cell_ref = null;
        // Tag layer and datatype to allow this element to be filtered out from LD and geo lists.
        layer_nr = -1;
        datatype_nr = -1;
        trans = new ();
        trans.reset();
    }

    public override GCCell depend()
    {
        return pDepend();
    }

    private GCCell pDepend()
    {
        return cell_ref;
    }

    public override void scale(double factor)
    {
        trans.mag = factor;
    }

    public override void rotate(double angle)
    {
        trans.angle = angle;
    }

    public override void minimum(ref Point64 p)
    {
        pMinimum(ref p);
    }

    private void pMinimum(ref Point64 p)
    {
        for (int x = 0; x < 2; x++)
        {
            for (int y = 0; y < 2; y++)
            {
                Point64 pos3 = GeoWrangler.move(point, (repetition.rowVector.X + repetition.colVector.X) * x * (repetition.columns - 1), (repetition.rowVector.Y + repetition.colVector.Y) * y * (repetition.rows - 1));
                Point64 pos1 = new(p.X - pos3.X, p.Y - pos3.Y);
                pos1.Y = trans.mirror_x switch
                {
                    true => -pos1.Y,
                    _ => pos1.Y
                };
                Point64 pos2 = pos1;
                cell_ref.maximum(ref pos1);
                cell_ref.minimum(ref pos2);
                switch (trans.mirror_x)
                {
                    case true:
                        pos1.Y = -pos1.Y;
                        pos2.Y = -pos2.Y;
                        break;
                }
                pos1 = GeoWrangler.move(pos1, pos3.X, pos3.Y);
                pos2 = GeoWrangler.move(pos2, pos3.X, pos3.Y);
                p.X = Math.Min(p.X, pos2.X);
                p.X = Math.Min(p.X, pos1.X);
                p.Y = Math.Min(p.Y, pos2.Y);
                p.Y = Math.Min(p.Y, pos1.Y);
            }
        }
    }

    public override void maximum(ref Point64 p)
    {
        pMaximum(ref p);
    }

    private void pMaximum(ref Point64 p)
    {
        for (int x = 0; x < 2; x++)
        {
            for (int y = 0; y < 2; y++)
            {
                Point64 pos3 = GeoWrangler.move(point, (repetition.rowVector.X + repetition.colVector.X) * x * (repetition.columns - 1), (repetition.rowVector.Y + repetition.colVector.Y) * y * (repetition.rows - 1));
                Point64 pos1 = new(p.X - pos3.X, p.Y - pos3.Y);
                pos1.Y = trans.mirror_x switch
                {
                    true => -pos1.Y,
                    _ => pos1.Y
                };
                Point64 pos2 = pos1;
                cell_ref.maximum(ref pos1);
                cell_ref.minimum(ref pos2);
                switch (trans.mirror_x)
                {
                    case true:
                        pos1.Y = -pos1.Y;
                        pos2.Y = -pos2.Y;
                        break;
                }
                pos1 = GeoWrangler.move(pos1, pos3.X, pos3.Y);
                pos2 = GeoWrangler.move(pos2, pos3.X, pos3.Y);
                p.X = Math.Max(p.X, pos1.X);
                p.X = Math.Max(p.X, pos2.X);
                p.Y = Math.Max(p.Y, pos1.Y);
                p.Y = Math.Max(p.Y, pos2.Y);
            }
        }
    }

    public override void move(Point64 p)
    {
        pMove(p);
    }

    private void pMove(Point64 p)
    {
        point = new (point.X + p.X, point.Y + p.Y);
    }

    public override void moveSelect(Point64 p)
    {
        pMoveSelect(p);
    }

    private void pMoveSelect(Point64 p)
    {
        switch (select)
        {
            case true:
                pMove(p);
                break;
        }
    }

    public override void setMirrorx()
    {
        pSetMirrors();
    }

    private void pSetMirrors()
    {
        trans.setMirror_x();
    }

    public override void resize(double factor)
    {
        pResize(factor);
    }

    private void pResize(double factor)
    {
        point.X = (int)(point.X * factor);
        point.Y = (int)(point.Y * factor);
        repetition.rowVector.X = (int)(repetition.rowVector.X * factor);
        repetition.rowVector.Y = (int)(repetition.rowVector.Y * factor);
        repetition.colVector.X = (int)(repetition.colVector.X * factor);
        repetition.colVector.Y = (int)(repetition.colVector.Y * factor);
    }

    public override void setPos(Point64 p)
    {
        pSetPos(p);
    }

    private void pSetPos(Point64 p)
    {
        point = new(p);
    }

    public override double getScale()
    {
        return pGetScale();
    }

    private double pGetScale()
    {
        return trans.mag;
    }

    public override double getAngle()
    {
        return pGetAngle();
    }

    private double pGetAngle()
    {
        return trans.angle;
    }

    public override bool getMirrorX()
    {
        return pGetMirrorX();
    }

    private bool pGetMirrorX()
    {
        return trans.mirror_x;
    }
    public override void setCellRef(GCCell cellRef)
    {
        pSetCellRef(cellRef);
    }

    private void pSetCellRef(GCCell cellRef)
    {
        cell_ref = cellRef;
    }

    public override void setCount(Point64 pt)
    {
        repetition.columns = (int)pt.X;
        repetition.rows = (int)pt.Y;
    }

    public override void setRowPitch(Point64 pt)
    {
        repetition.rowVector = new(pt);
    }
    public override Point64 getRowPitch()
    {
        return new(repetition.rowVector);
    }
    
    public override void setColPitch(Point64 pt)
    {
        repetition.colVector = new(pt);
    }

    public override Point64 getColPitch()
    {
        return new(repetition.colVector);
    }

    public override Point64 getCount()
    {
        return new(repetition.columns, repetition.rows);
    }

    public override void saveGDS(gdsWriter gw)
    {
        pSaveGDS(gw);
    }

    // Cellrefarrays also have to resolve to integer placement.
    // Placement errors will occur if the x, y instance counts do not divide the array X, Y values cleanly.
    private void pSaveGDS(gdsWriter gw)
    {
        //SRef
        gw.bw.Write((ushort)4);
        gw.bw.Write(gdsValues.sAREF);
        gw.writeString(cell_ref.cellName, 0x12);
        int strans_ = 0;
        switch (trans.mirror_x)
        {
            case true:
                strans_ |= 0x8000;
                break;
        }
        //STRANS
        gw.bw.Write((ushort)6);
        gw.bw.Write(gdsValues.sSTRANS);
        gw.bw.Write((short)strans_);
        //mag
        gw.bw.Write((ushort)12);
        gw.bw.Write(gdsValues.sMAG);
        gw.write8ByteReal(trans.mag);
        //angle
        gw.bw.Write((ushort)12);
        gw.bw.Write(gdsValues.sANGLE);
        /*
        if (trans.mirror_x && trans.angle != 0)
        {
            gw.write8ByteReal(360 - trans.angle);
        }
        else
        */
        {
            gw.write8ByteReal(trans.angle);
        }

        //colrow
        gw.bw.Write((ushort)8);
        gw.bw.Write(gdsValues.sCOLROW);
        gw.bw.Write((short)repetition.columns);
        gw.bw.Write((short)repetition.rows);
        //xy
        gw.bw.Write((ushort)(3 * 2 * 4 + 4));
        gw.bw.Write(gdsValues.sXY);
        gw.bw.Write((int)point.X);
        gw.bw.Write((int)point.Y);
        Point64 pos = new(repetition.rowVector.X + point.X, (repetition.rowVector.Y * repetition.rows) + point.Y);
        gw.bw.Write((int)pos.X);
        gw.bw.Write((int)pos.Y);
        pos = new ((repetition.colVector.X * repetition.columns) + point.X, repetition.colVector.Y + point.Y);
        gw.bw.Write((int)pos.X);
        gw.bw.Write((int)pos.Y);
        // endel
        gw.bw.Write((ushort)4);
        gw.bw.Write(gdsValues.sENDEL);
    }

    public override void saveOASIS(oasWriter ow)
    {
        pSaveOASIS(ow);
    }

    // Cellrefarrays also have to resolve to integer placement.
    // Placement errors will occur if the x, y instance counts do not divide the array X, Y values cleanly.
    private void pSaveOASIS(oasWriter ow)
    {
        ow.modal.absoluteMode = ow.modal.absoluteMode switch
        {
            false => true,
            _ => ow.modal.absoluteMode
        };
        byte info_byte = 128;  //explicid cellname,repition;
        if (repetition.columns > 1 || repetition.rows > 1)
        {
            info_byte += 8;
        }

        switch (trans.mirror_x)
        {
            case true:
                info_byte += 1;
                break;
        }

        switch (Math.Abs(trans.mag - 1))
        {
            case > double.Epsilon:
                info_byte += 4;
                break;
        }

        if (trans.angle != 0)
        {
            info_byte += 2;
        }

        if (point.X != ow.modal.placement_x)
        {
            info_byte += 32;
        }

        if (point.Y != ow.modal.placement_y)
        {
            info_byte += 16;
        }

        ow.writeUnsignedInteger(18);
        ow.writeRaw(info_byte);
        ow.writeString(cell_ref.cellName);

        switch (info_byte & 4)
        {
            case > 0:
                ow.writeReal(trans.mag);
                break;
        }

        switch (info_byte & 2)
        {
            case > 0:
                ow.writeReal(trans.angle);
                break;
        }
        switch (info_byte & 32)
        {
            case > 0:
                ow.modal.placement_x = (int)point.X;
                ow.writeSignedInteger(ow.modal.placement_x);
                break;
        }
        switch (info_byte & 16)
        {
            case > 0:
                ow.modal.placement_y = (int)point.Y;
                ow.writeSignedInteger(ow.modal.placement_y);
                break;
        }

        switch (info_byte & 8)
        {
            default:
                switch (repetition.type)
                {
                    case Repetition.RepetitionType.Rectangular:
                        if ((repetition.columns > 1) && (repetition.rows > 1))
                        {
                            ow.modal.x_dimension = repetition.columns;
                            ow.modal.y_dimension = repetition.rows;
                            ow.modal.x_space = (int)(repetition.colVector.X + repetition.rowVector.X);
                            ow.modal.y_space = (int)(repetition.colVector.Y + repetition.rowVector.Y);
                            if (((repetition.colVector.X + repetition.rowVector.X) >= 0) &&
                                ((repetition.colVector.Y + repetition.rowVector.Y) >= 0))
                            {
                                ow.writeUnsignedInteger(1);
                                ow.writeUnsignedInteger((uint)(ow.modal.x_dimension - 2));
                                ow.writeUnsignedInteger((uint)(ow.modal.y_dimension - 2));
                                ow.writeUnsignedInteger((uint)ow.modal.x_space);
                                ow.writeUnsignedInteger((uint)ow.modal.y_space);
                            }
                            else
                            {
                                ow.writeUnsignedInteger(8);
                                ow.writeUnsignedInteger((uint)(ow.modal.x_dimension - 2));
                                ow.writeUnsignedInteger((uint)(ow.modal.y_dimension - 2));
                                ow.writeGDelta(new(ow.modal.x_space, 0));
                                ow.writeGDelta(new(0, ow.modal.y_space));
                            }
                        }
                        else if (repetition.columns > 1)
                        {
                            ow.modal.x_dimension = repetition.columns;
                            ow.modal.x_space = (int)(repetition.colVector.X + repetition.rowVector.X);
                            if ((repetition.colVector.X + repetition.rowVector.X) >= 0)
                            {
                                ow.writeUnsignedInteger(2);
                                ow.writeUnsignedInteger((uint)(ow.modal.x_dimension - 2));
                                ow.writeUnsignedInteger((uint)ow.modal.x_space);
                            }
                            else
                            {
                                ow.writeUnsignedInteger(9);
                                ow.writeUnsignedInteger((uint)(ow.modal.x_dimension - 2));
                                ow.writeGDelta(new(ow.modal.x_space, 0));
                            }
                        }
                        else
                        {
                            ow.modal.y_dimension = repetition.rows;
                            ow.modal.y_space = (int)repetition.rowVector.Y;
                            if ((repetition.colVector.Y + repetition.rowVector.Y) >= 0)
                            {
                                ow.writeUnsignedInteger(3);
                                ow.writeUnsignedInteger((uint)(ow.modal.y_dimension - 2));
                                ow.writeUnsignedInteger((uint)ow.modal.y_space);
                            }
                            else
                            {
                                ow.writeUnsignedInteger(9);
                                ow.writeUnsignedInteger((uint)(ow.modal.y_dimension - 2));
                                ow.writeGDelta(new(0, ow.modal.y_space));
                            }
                        }

                        break;
                    case Repetition.RepetitionType.Regular:
                        if ((repetition.columns > 1) && (repetition.rows > 1))
                        {
                            ow.modal.x_dimension = repetition.columns;
                            ow.modal.y_dimension = repetition.rows;
                            ow.writeUnsignedInteger(8);
                            ow.writeUnsignedInteger((uint)(ow.modal.x_dimension - 2));
                            ow.writeUnsignedInteger((uint)(ow.modal.y_dimension - 2));
                            ow.writeGDelta(new(repetition.colVector));
                            ow.writeGDelta(new(repetition.rowVector));
                        }
                        else if (repetition.columns > 1)
                        {
                            ow.modal.x_dimension = repetition.columns;
                            ow.writeUnsignedInteger(9);
                            ow.writeUnsignedInteger((uint)(ow.modal.x_dimension - 2));
                            ow.writeGDelta(new(repetition.colVector));
                        }
                        else
                        {
                            ow.modal.y_dimension = repetition.rows;
                            ow.writeUnsignedInteger(9);
                            ow.writeUnsignedInteger((uint)(ow.modal.y_dimension - 2));
                            ow.writeGDelta(new(repetition.rowVector));
                        }

                        break;
                    case Repetition.RepetitionType.ExplicitX:
                        if (repetition.coords.Count > 0)
                        {
                            ow.writeUnsignedInteger(4);
                            ow.writeUnsignedInteger((uint)(repetition.coords.Count - 1));
                            int c0_index = 0;
                            int c1_index = 1;
                            ow.writeUnsignedInteger((uint)Math.Round(repetition.coords[c0_index],
                                MidpointRounding.ToEven));
                            for (int i = repetition.coords.Count - 1; i > 0; --i)
                            {
                                c0_index = (c0_index + 1) % repetition.coords.Count;
                                c1_index = (c1_index + 1) % repetition.coords.Count;
                                ow.writeUnsignedInteger(
                                    (uint)Math.Round(repetition.coords[c1_index] - repetition.coords[c0_index]));
                            }
                        }

                        break;
                    case Repetition.RepetitionType.ExplicitY:
                        if (repetition.coords.Count > 0)
                        {
                            ow.writeUnsignedInteger(6);
                            ow.writeUnsignedInteger((uint)(repetition.coords.Count - 1));
                            int c0_index = 0;
                            int c1_index = 1;
                            ow.writeUnsignedInteger((uint)Math.Round(repetition.coords[c0_index],
                                MidpointRounding.ToEven));
                            for (int i = repetition.coords.Count - 1; i > 0; --i)
                            {
                                c0_index = (c0_index + 1) % repetition.coords.Count;
                                c1_index = (c1_index + 1) % repetition.coords.Count;
                                ow.writeUnsignedInteger(
                                    (uint)Math.Round(repetition.coords[c1_index] - repetition.coords[c0_index]));
                            }
                        }

                        break;
                    case Repetition.RepetitionType.Explicit:
                        if (repetition.offsets.Count > 0)
                        {
                            ow.writeUnsignedInteger(10);
                            ow.writeUnsignedInteger((uint)(repetition.offsets.Count - 1));
                            int v0_index = 0;
                            int v1_index = 1;
                            ow.writeGDelta(new(repetition.offsets[v0_index]));
                            for (int i = repetition.offsets.Count - 1; i > 0; --i)
                            {
                                v0_index = (v0_index + 1) % repetition.offsets.Count;
                                v1_index = (v1_index + 1) % repetition.offsets.Count;
                                ow.writeGDelta(new(repetition.offsets[v1_index].X - repetition.offsets[v0_index].X,
                                    repetition.offsets[v1_index].Y - repetition.offsets[v0_index].Y));
                            }
                        }

                        break;
                }

                break;
        }
    }

    public override bool isCellrefArray()
    {
        return pIsCellRefArray();
    }

    private bool pIsCellRefArray()
    {
        return true;
    }

    public override GCCell getCellref() {
        return pGetCellRef();
    }

    private GCCell pGetCellRef()
    {
        return cell_ref;
    }

    public override Point64 getPos()
    {
        return pGetPos();
    }

    private Point64 pGetPos()
    {
        return point;
    }

    public override List<GCPolygon> convertToPolygons()
    {
        return pConvertToPolygons();
    }


    private List<GCPolygon> pConvertToPolygons()
    {
        /*
        // This returns the contents of the reference.
        List<GCPolygon> tmp = new List<GCPolygon>();
        for (int element = 0; element < cell_ref.elementList.Count; element++)
        {
            List<GCPolygon> lp = cell_ref.elementList[element].convertToPolygons().ToList();
            // Above should work with nested arrays due to the recursion of the convertToPolygons override.
            tmp.AddRange(lp);
        }
        */
        
        List<GCPolygon> tmp = cell_ref.convertToPolygons();

        List<GCPolygon> ret = repetition.transform(tmp, trans.mag, trans.mirror_x, trans.angle);

        Parallel.For(0, ret.Count, (p) =>
        {
            ret[p].move(point);
        });

        return ret;
    }

}