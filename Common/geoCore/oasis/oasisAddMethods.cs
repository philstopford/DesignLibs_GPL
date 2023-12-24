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
        switch (modal.mirror_x)
        {
            case true:
                cell_.elementList[^1].setMirrorx();
                break;
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
        Path64 pa = Helper.initedPath64(modal.polygon_point_list.Count);
        Point64 p = new(modal.geometry_x, modal.geometry_y);
        for (int i = 0; i < modal.polygon_point_list.Count; i++)
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
        Path64 offsets = modal.repetition.get_offsets();
        foreach (Point64 offset in offsets)
        {
            addElement(e, offset);
        }

        /*
        if (e == elementType.cellrefElement)
        {
            switch (modal.repetition)
            {
                // Rectangular for m columns and n rows both >= 1, with spacing in X and Y both >= 0
                case 1:
                    cell_.addCellrefArray(drawing_.findCell(modal.placement_cell),
                        new(modal.placement_x, modal.placement_y),
                        new(modal.x_space + modal.placement_x, modal.y_space + modal.placement_y), modal.x_dimension,
                        modal.y_dimension);
                    cell_.elementList[^1].setName(modal.placement_cell);
                    cell_.elementList[^1].rotate(modal.angle);
                    cell_.elementList[^1].scale(modal.mag);
                    switch (modal.mirror_x)
                    {
                        case true:
                            cell_.elementList[^1].setMirrorx();
                            break;
                    }
                    break;
                // Rectangular m columns for X >= 0
                case 2:
                    cell_.addCellrefArray(drawing_.findCell(modal.placement_cell),
                        new(modal.placement_x, modal.placement_y),
                        new(modal.x_space + modal.placement_x, modal.placement_y), modal.x_dimension, 1);
                    cell_.elementList[^1].setName(modal.placement_cell);
                    cell_.elementList[^1].rotate(modal.angle);
                    cell_.elementList[^1].scale(modal.mag);
                    switch (modal.mirror_x)
                    {
                        case true:
                            cell_.elementList[^1].setMirrorx();
                            break;
                    }
                    break;
                // Rectangular n rows for Y >= 0
                case 3:
                    cell_.addCellrefArray(drawing_.findCell(modal.placement_cell),
                        new(modal.placement_x, modal.placement_y),
                        new(modal.placement_x, modal.y_space + modal.placement_y), 1, modal.y_dimension);
                    cell_.elementList[^1].setName(modal.placement_cell);
                    cell_.elementList[^1].rotate(modal.angle);
                    cell_.elementList[^1].scale(modal.mag);
                    switch (modal.mirror_x)
                    {
                        case true:
                            cell_.elementList[^1].setMirrorx();
                            break;
                    }
                    break;
                // ExplicitX with coords > 0
                case 4:
                    break;
                // ExplicitY with coords > 0
                case 6:
                    break;
                // Regular with m cols and n rows both > 1
                case 8:
                    break;
                default:
                    // We flatten the placements. It might be better to keep the
                    // information via the repetition system.
                    // Need to understand our repetition value to do this properly.
                    // This maps to ExplicitX, ExplicitY, ExplicitXY handling, I think.
                    // 
                    // cell_.addCellrefArray(drawing_.findCell(modal.placement_cell), modal.repArray);
                    // cell_.elementList[^1].setName(modal.placement_cell);
                    // cell_.elementList[^1].rotate(modal.angle);
                    // cell_.elementList[^1].scale(modal.mag);
                    // switch (modal.mirror_x)
                    // {
                    //     case true:
                    //         cell_.elementList[^1].setMirrorx();
                    //         break;
                    // }
                    // 
                    foreach (Point64 t in modal.repArray)
                    {
                        addElement(e, t);
                    }
                    // Throwing to allow this to be investigated....
                    throw new Exception("CellRef with unsupported repetition mode in processRepetition");
                    break;
            }
        }
        else
        {
            switch (modal.repetition)
            {
                case 1:
                    for (int x = 0; x < modal.x_dimension; x++)
                    for (int y = 0; y < modal.y_dimension; y++)
                    {
                        addElement(e, new (x * modal.x_space, y * modal.y_space));
                    }
                    break;
                case 2:
                    for (int x = 0; x < modal.x_dimension; x++)
                    {
                        addElement(e, new (x * modal.x_space, 0));
                    }
                    break;
                case 3:
                    for (int y = 0; y < modal.y_dimension; y++)
                    {
                        addElement(e, new (0, y * modal.y_space));
                    }
                    break;
                default:
                    // We flatten the placements. It might be better to keep the
                    // information via the repetition system.
                    // Need to understand our repetition value to do this properly.
                    // This maps to ExplicitX, ExplicitY, ExplicitXY handling, I think.
                    foreach (Point64 t in modal.repArray)
                    {
                        addElement(e, t);
                    }
                    // Throwing to allow this to be investigated....
                    throw new Exception("Found unsupported repetition mode in processRepetition");
                    break;
            }
        }
        */
    }
}