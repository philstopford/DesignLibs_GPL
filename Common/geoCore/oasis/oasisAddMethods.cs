using geoCoreLib;
using System;
using Clipper2Lib;
using geoWrangler;

namespace oasis;

internal partial class oasReader
{
    private void registerLayerDatatype()
    {
        if (!layerNames.ContainsKey("L" + modal.layer + "D" + modal.datatype))
        {
            layerNames.Add("L" + modal.layer + "D" + modal.datatype, "L" + modal.layer + "D" + modal.datatype);
        }
    }

    private void addElement(elementType e, Point64 p)
    {
        int x;
        int y;
        switch (e)
        {
            case elementType.boxElement:
                x = modal.geometry_x;
                y = modal.geometry_y;
                modal.geometry_x += (int)p.X;
                modal.geometry_y += (int)p.Y;
                addBox();
                modal.geometry_x = x;
                modal.geometry_y = y;
                break;
            case elementType.polygonElement:
                x = modal.geometry_x;
                y = modal.geometry_y;
                modal.geometry_x += (int)p.X;
                modal.geometry_y += (int)p.Y;
                addPolygon();
                modal.geometry_x = x;
                modal.geometry_y = y;
                break;
            case elementType.pathElement:
                x = modal.geometry_x;
                y = modal.geometry_y;
                modal.geometry_x += (int)p.X;
                modal.geometry_y += (int)p.Y;
                addPath();
                modal.geometry_x = x;
                modal.geometry_y = y;
                break;
            case elementType.cellrefElement:
                x = modal.placement_x;
                y = modal.placement_y;
                modal.placement_x += (int)p.X;
                modal.placement_y += (int)p.Y;
                addCellref();
                modal.placement_x = x;
                modal.placement_y = y;
                break;
            case elementType.textElement:
                x = modal.text_x;
                y = modal.text_y;
                modal.text_x += (int)p.X;
                modal.text_y += (int)p.Y;
                addText();
                modal.text_x = x;
                modal.text_y = y;
                break;
            case elementType.circleElement:
                x = modal.text_x;
                y = modal.text_y;
                modal.text_x += (int)p.X;
                modal.text_y += (int)p.Y;
                addCircle();
                modal.text_x = x;
                modal.text_y = y;
                break;
            case elementType.trapezoidElement:
                x = modal.text_x;
                y = modal.text_y;
                modal.text_x += (int)p.X;
                modal.text_y += (int)p.Y;
                addTrapezoid();
                modal.text_x = x;
                modal.text_y = y;
                break;
            case elementType.ctrapezoidElement:
                x = modal.text_x;
                y = modal.text_y;
                modal.text_x += (int)p.X;
                modal.text_y += (int)p.Y;
                addCtrapezoid();
                modal.text_x = x;
                modal.text_y = y;
                break;
        }
    }

    private void addBox()
    {
        cell_.addBox(modal.geometry_x, modal.geometry_y, modal.geometry_w, modal.geometry_h, modal.layer, modal.datatype);
        registerLayerDatatype();
    }

    private void addCellref()
    {
        cell_.addCellref();
        cell_.elementList[^1].setPos(new (modal.placement_x, modal.placement_y));
        cell_.elementList[^1].setCellRef(drawing_.findCell(modal.placement_cell));
        cell_.elementList[^1].setName(modal.placement_cell);
        cell_.elementList[^1].rotate(modal.angle);
        cell_.elementList[^1].scale(modal.mag);
        if (modal.mirror_x)
        {
            cell_.elementList[^1].setMirrorx();
        }
    }

    private void addCircle()
    {
        cell_.addCircle(modal.layer, modal.datatype, new (modal.geometry_x, modal.geometry_y), modal.circle_radius);
        registerLayerDatatype();
    }

    private void addPath()
    {
        Path64 pa = Helper.initedPath64(modal.polygon_point_list.Count);
        Point64 p = new(modal.geometry_x, modal.geometry_y);
        for (int i = 0; i < modal.polygon_point_list.Count; i++)
        {
            pa[i] = modal.polygon_point_list[i];
            pa[i] = GeoWrangler.move(pa[i], p.X, p.Y);
        }
        cell_.addPath(pa, modal.layer, modal.datatype);
        cell_.elementList[^1].setWidth(2 * modal.geometry_w);
        switch (modal.path_start_extension)
        {
            case 2 when modal.path_end_extension == 2:
                cell_.elementList[^1].setCap(2);
                break;
            case 0 when modal.path_end_extension == 0:
                cell_.elementList[^1].setCap(0);
                break;
            default:
            {
                int start = modal.path_start_extension_value;
                int ende = modal.path_end_extension_value;
                start = modal.path_start_extension switch
                {
                    2 => modal.geometry_w,
                    _ => start
                };
                ende = modal.path_end_extension switch
                {
                    2 => modal.geometry_w,
                    _ => ende
                };
                cell_.elementList[^1].expandCaps(start, ende);
                break;
            }
        }
        registerLayerDatatype();
    }

    private void addPolygon()
    {
        int polySize = modal.polygon_point_list.Count;
        // Let's see if we need to trim a bit based on situation.
        if (polySize >= 5)
        {
            string hash0 = utility.Utils.GetMD5Hash(modal.polygon_point_list[0]);
            string hash1 = utility.Utils.GetMD5Hash(modal.polygon_point_list[^1]);
            string hash2 = utility.Utils.GetMD5Hash(modal.polygon_point_list[^2]);
            string hash3 = utility.Utils.GetMD5Hash(modal.polygon_point_list[^3]);
            string hash4 = utility.Utils.GetMD5Hash(modal.polygon_point_list[^4]);
            if (hash1 == hash0 && hash2 == hash0 && hash3 == hash0 && hash4 == hash0)
            {
                polySize -= 2;
            }
        }

        Path64 pa = Helper.initedPath64(polySize);
        Point64 p = new(modal.geometry_x, modal.geometry_y);
        for (int i = 0; i < polySize; i++)
        {
            pa[i] = new (modal.polygon_point_list[i]);
            pa[i] = GeoWrangler.move(pa[i], p.X, p.Y);
        }
        cell_.addPolygon(pa, modal.layer, modal.datatype);
        registerLayerDatatype();
    }

    private void addText()
    {
        cell_.addText(modal.textlayer, modal.datatype, new (modal.text_x, modal.text_y), modal.text_string);
        cell_.elementList[^1].setWidth(GCSetup.defaultTextWidth);
        cell_.elementList[^1].setPresentation(GCSetup.defaultTextPresentation);
        registerLayerDatatype();
    }

    private void addTrapezoid()
    {
        Path64 pa = Helper.initedPath64(5);

        switch (modal.trapezoid_orientation)
        {
            // (m & 0x80)
            case true:
                //  vertically
                pa[0] = new (modal.geometry_x, modal.geometry_y + Math.Max(modal.trapezoid_delta_a, 0));
                pa[1] = new (modal.geometry_x, modal.geometry_y + modal.geometry_h + Math.Min(modal.trapezoid_delta_b, 0));
                pa[2] = new (modal.geometry_x + modal.geometry_w, modal.geometry_y + modal.geometry_h - Math.Max(modal.trapezoid_delta_b, 0));
                pa[3] = new (modal.geometry_x + modal.geometry_w, modal.geometry_y - Math.Min(modal.trapezoid_delta_a, 0));
                break;
            default:
                //  horizontally
                pa[0] = new (modal.geometry_x + Math.Max(modal.trapezoid_delta_a, 0), modal.geometry_y + modal.geometry_h);
                pa[1] = new (modal.geometry_x + modal.geometry_w + Math.Min(modal.trapezoid_delta_b, 0), modal.geometry_y + modal.geometry_h);
                pa[2] = new (modal.geometry_x + modal.geometry_w - Math.Max(modal.trapezoid_delta_b, 0), modal.geometry_y);
                pa[3] = new (modal.geometry_x - Math.Min(modal.trapezoid_delta_a, 0), modal.geometry_y);
                break;
        }

        pa[4] = new (pa[0]);

        cell_.addPolygon(pa, modal.layer, modal.datatype);
    }

    private void addCtrapezoid()
    {
        Path64 pa = Helper.initedPath64(5);

        int[,] coords = modal.ctrapezoid_type switch
        {
            0 => new[,]
            {
                {0, 0, 0, 0}, // x=0*w+0*h, y=0*w+0*h ...
                {0, 0, 0, 1}, {1, -1, 0, 1}, {1, 0, 0, 0}
            },
            1 => new[,] {{0, 0, 0, 0}, {0, 0, 0, 1}, {1, 0, 0, 1}, {1, -1, 0, 0}},
            2 => new[,] {{0, 0, 0, 0}, {0, 1, 0, 1}, {1, 0, 0, 1}, {1, 0, 0, 0}},
            3 => new[,] {{0, 1, 0, 0}, {0, 0, 0, 1}, {1, 0, 0, 1}, {1, 0, 0, 0}},
            4 => new[,] {{0, 0, 0, 0}, {0, 1, 0, 1}, {1, -1, 0, 1}, {1, 0, 0, 0}},
            5 => new[,] {{0, 1, 0, 0}, {0, 0, 0, 1}, {1, 0, 0, 1}, {1, -1, 0, 0}},
            6 => new[,] {{0, 0, 0, 0}, {0, 1, 0, 1}, {1, 0, 0, 1}, {1, -1, 0, 0}},
            7 => new[,] {{0, 1, 0, 0}, {0, 0, 0, 1}, {1, -1, 0, 1}, {1, 0, 0, 0}},
            8 => new[,] {{0, 0, 0, 0}, {0, 0, 0, 1}, {1, 0, -1, 1}, {1, 0, 0, 0}},
            9 => new[,] {{0, 0, 0, 0}, {0, 0, -1, 1}, {1, 0, 0, 1}, {1, 0, 0, 0}},
            10 => new[,] {{0, 0, 0, 0}, {0, 0, 0, 1}, {1, 0, 0, 1}, {1, 0, 1, 0}},
            11 => new[,] {{0, 0, 1, 0}, {0, 0, 0, 1}, {1, 0, 0, 1}, {1, 0, 0, 0}},
            12 => new[,] {{0, 0, 0, 0}, {0, 0, 0, 1}, {1, 0, -1, 1}, {1, 0, 1, 0}},
            13 => new[,] {{0, 0, 1, 0}, {0, 0, -1, 1}, {1, 0, 0, 1}, {1, 0, 0, 0}},
            14 => new[,] {{0, 0, 0, 0}, {0, 0, -1, 1}, {1, 0, 0, 1}, {1, 0, 1, 0}},
            15 => new[,] {{0, 0, 1, 0}, {0, 0, 0, 1}, {1, 0, -1, 1}, {1, 0, 0, 0}},
            16 => new[,] {{0, 0, 0, 0}, {0, 0, 1, 0}, {1, 0, 0, 0}, {0, 0, 0, 0}},
            17 => new[,] {{0, 0, 0, 0}, {0, 0, 1, 0}, {1, 0, 1, 0}, {0, 0, 0, 0}},
            18 => new[,] {{0, 0, 0, 0}, {1, 0, 1, 0}, {1, 0, 0, 0}, {0, 0, 0, 0}},
            19 => new[,] {{0, 0, 1, 0}, {1, 0, 1, 0}, {1, 0, 0, 0}, {0, 0, 1, 0}},
            20 => new[,] {{0, 0, 0, 0}, {0, 1, 0, 1}, {0, 2, 0, 0}, {0, 0, 0, 0}},
            21 => new[,] {{0, 0, 0, 1}, {0, 2, 0, 1}, {0, 1, 0, 0}, {0, 0, 0, 1}},
            22 => new[,] {{0, 0, 0, 0}, {0, 0, 2, 0}, {1, 0, 1, 0}, {0, 0, 0, 0}},
            23 => new[,] {{1, 0, 0, 0}, {0, 0, 1, 0}, {1, 0, 2, 0}, {1, 0, 0, 0}},
            24 => new[,] {{0, 0, 0, 0}, {0, 0, 0, 1}, {1, 0, 0, 1}, {1, 0, 0, 0}},
            25 => new[,] {{0, 0, 0, 0}, {0, 0, 1, 0}, {1, 0, 1, 0}, {1, 0, 0, 0}},
            _ => new int[4, 4]
        };

        for (int pt = 0; pt < 4; pt++)
        {
            int x = 0;
            if (coords[pt, 0] != 0)
            {
                x += coords[pt, 0] * modal.geometry_w;
            }

            if (coords[pt, 1] != 0)
            {
                x += coords[pt, 1] * modal.geometry_h;
            }

            int y = 0;
            if (coords[pt, 2] != 0)
            {
                y += coords[pt, 2] * modal.geometry_w;
            }

            if (coords[pt, 3] != 0)
            {
                y += coords[pt, 3] * modal.geometry_h;
            }

            pa[pt] = new (modal.geometry_x + x, modal.geometry_y + y);

            if (x > modal.geometry_w)
            {
                modal.geometry_w = x;
            }
            if (y > modal.geometry_h)
            {
                modal.geometry_h = y;
            }
        }

        pa[^1] = new (pa[0]);

        cell_.addPolygon(pa, modal.layer, modal.datatype);
    }

    private void processRepetition(elementType e)
    {
        // The aim here is to ensure we have a consistent representation with the input.
        // If the input is an element, but not a reference, we should just create duplicates of the element.
        Path64 offsets = modal.repetition.get_offsets();
        if (e == elementType.cellrefElement)
        {
            if ((modal.repetition.type == Repetition.RepetitionType.Regular) || (modal.repetition.type == Repetition.RepetitionType.Rectangular))
            {
                cell_.addCellrefArray(drawing_.findCell(modal.placement_cell), new (modal.placement_x, modal.placement_y), modal.angle, modal.mag, modal.mirror_x, modal.repetition);
            }
            else
            {
                // Original is not included in the offsets.
                cell_.addCellref();
                cell_.elementList[^1].setPos(new(modal.placement_x, modal.placement_y));
                cell_.elementList[^1].setCellRef(drawing_.findCell(modal.placement_cell));
                cell_.elementList[^1].setName(modal.placement_cell);
                cell_.elementList[^1].rotate(modal.angle);
                cell_.elementList[^1].scale(modal.mag);
                if (modal.mirror_x)
                {
                    cell_.elementList[^1].setMirrorx();
                }

                foreach (Point64 offset in offsets)
                {
                    cell_.addCellref();
                    cell_.elementList[^1].setPos(new(modal.placement_x + offset.X, modal.placement_y + offset.Y));
                    cell_.elementList[^1].setCellRef(drawing_.findCell(modal.placement_cell));
                    cell_.elementList[^1].setName(modal.placement_cell);
                    cell_.elementList[^1].rotate(modal.angle);
                    cell_.elementList[^1].scale(modal.mag);
                    if (modal.mirror_x)
                    {
                        cell_.elementList[^1].setMirrorx();
                    }
                }
            }
        }
        else
        {
            // To be checked - only add the base version if offset is missing it. Not sure if this is a good idea.
            if ((offsets[0].X != 0) || (offsets[0].Y != 0))
            {
                addElement(e, new(0, 0));
            }
            foreach (Point64 offset in offsets)
            {
                addElement(e, offset);
            }
        }
    }
}