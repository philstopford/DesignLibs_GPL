using geoCoreLib;
using geoLib;
using System;

namespace oasis
{
    partial class oasReader
    {
        void addElement(elementType e, GeoLibPoint p)
        {
            Int32 x;
            Int32 y;
            switch (e)
            {
                case elementType.boxElement:
                    x = modal.geometry_x;
                    y = modal.geometry_y;
                    modal.geometry_x += p.X;
                    modal.geometry_y += p.Y;
                    addBox();
                    modal.geometry_x = x;
                    modal.geometry_y = y;
                    break;
                case elementType.polygonElement:
                    x = modal.geometry_x;
                    y = modal.geometry_y;
                    modal.geometry_x += p.X;
                    modal.geometry_y += p.Y;
                    addPolygon();
                    modal.geometry_x = x;
                    modal.geometry_y = y;
                    break;
                case elementType.pathElement:
                    x = modal.geometry_x;
                    y = modal.geometry_y;
                    modal.geometry_x += p.X;
                    modal.geometry_y += p.Y;
                    addPath();
                    modal.geometry_x = x;
                    modal.geometry_y = y;
                    break;
                case elementType.cellrefElement:
                    x = modal.placement_x;
                    y = modal.placement_y;
                    modal.placement_x += p.X;
                    modal.placement_y += p.Y;
                    addCellref();
                    modal.placement_x = x;
                    modal.placement_y = y;
                    break;
                case elementType.textElement:
                    x = modal.text_x;
                    y = modal.text_y;
                    modal.text_x += p.X;
                    modal.text_y += p.Y;
                    addText();
                    modal.text_x = x;
                    modal.text_y = y;
                    break;
                case elementType.circleElement:
                    x = modal.text_x;
                    y = modal.text_y;
                    modal.text_x += p.X;
                    modal.text_y += p.Y;
                    addCircle();
                    modal.text_x = x;
                    modal.text_y = y;
                    break;
                case elementType.trapezoidElement:
                    x = modal.text_x;
                    y = modal.text_y;
                    modal.text_x += p.X;
                    modal.text_y += p.Y;
                    addTrapezoid();
                    modal.text_x = x;
                    modal.text_y = y;
                    break;
                case elementType.ctrapezoidElement:
                    x = modal.text_x;
                    y = modal.text_y;
                    modal.text_x += p.X;
                    modal.text_y += p.Y;
                    addCtrapezoid();
                    modal.text_x = x;
                    modal.text_y = y;
                    break;
            }
        }

        void addBox()
        {
            cell_.addBox(modal.geometry_x, modal.geometry_y, modal.geometry_w, modal.geometry_h, modal.layer, modal.datatype);
        }

        void addCellref()
        {
            cell_.addCellref();
            cell_.elementList[^1].setPos(new GeoLibPoint(modal.placement_x, modal.placement_y));
            cell_.elementList[^1].setCellRef(drawing_.findCell(modal.placement_cell));
            cell_.elementList[^1].setName(modal.placement_cell);
            cell_.elementList[^1].rotate(modal.angle);
            cell_.elementList[^1].scale(modal.mag);
            if (modal.mirror_x)
            {
                cell_.elementList[^1].setMirrorx();
            }
        }

        void addCircle()
        {
            cell_.addCircle(modal.layer, modal.datatype, new GeoLibPoint(modal.geometry_x, modal.geometry_y), modal.circle_radius);
        }

        void addPath()
        {
            GeoLibPoint[] pa = new GeoLibPoint[modal.polygon_point_list.Count];
            GeoLibPoint p = new GeoLibPoint(modal.geometry_x, modal.geometry_y);
            for (Int32 i = 0; i < modal.polygon_point_list.Count; i++)
            {
                pa[i] = modal.polygon_point_list[i];
                pa[i].Offset(p.X, p.Y);
            }
            cell_.addPath(pa, modal.layer, modal.datatype);
            cell_.elementList[^1].setWidth(2 * modal.geometry_w);
            if ((modal.path_start_extension == 2) && (modal.path_end_extension == 2))
            {
                cell_.elementList[^1].setCap(2);
            }
            else if ((modal.path_start_extension == 0) && (modal.path_end_extension == 0))
            {
                cell_.elementList[^1].setCap(0);
            }
            else
            {
                Int32 start = modal.path_start_extension_value;
                Int32 ende = modal.path_end_extension_value;
                if (modal.path_start_extension == 2)
                {
                    start = modal.geometry_w;
                }
                if (modal.path_end_extension == 2)
                {
                    ende = modal.geometry_w;
                }
                cell_.elementList[^1].expandCaps(start, ende);
            }
        }

        void addPolygon()
        {
            GeoLibPoint[] pa = new GeoLibPoint[modal.polygon_point_list.Count];
            GeoLibPoint p = new GeoLibPoint(modal.geometry_x, modal.geometry_y);
            for (Int32 i = 0; i < modal.polygon_point_list.Count; i++)
            {
                pa[i] = new GeoLibPoint(modal.polygon_point_list[i]);
                pa[i].Offset(p.X, p.Y);
            }
            cell_.addPolygon(pa, modal.layer, modal.datatype);
        }

        void addText()
        {
            cell_.addText(modal.textlayer, modal.datatype, new GeoLibPoint(modal.text_x, modal.text_y), modal.text_string);
            cell_.elementList[^1].setWidth(GCSetup.defaultTextWidth);
            cell_.elementList[^1].setPresentation(GCSetup.defaultTextPresentation);
        }

        void addTrapezoid()
        {
            GeoLibPoint[] pa = new GeoLibPoint[5];

            if (modal.trapezoid_orientation) // (m & 0x80)
            {
                //  vertically
                pa[0] = new GeoLibPoint(modal.geometry_x, modal.geometry_y + Math.Max(modal.trapezoid_delta_a, 0));
                pa[1] = new GeoLibPoint(modal.geometry_x, modal.geometry_y + modal.geometry_h + Math.Min(modal.trapezoid_delta_b, 0));
                pa[2] = new GeoLibPoint(modal.geometry_x + modal.geometry_w, modal.geometry_y + modal.geometry_h - Math.Max(modal.trapezoid_delta_b, 0));
                pa[3] = new GeoLibPoint(modal.geometry_x + modal.geometry_w, modal.geometry_y - Math.Min(modal.trapezoid_delta_a, 0));
            }
            else
            {
                //  horizontally
                pa[0] = new GeoLibPoint(modal.geometry_x + Math.Max(modal.trapezoid_delta_a, 0), modal.geometry_y + modal.geometry_h);
                pa[1] = new GeoLibPoint(modal.geometry_x + modal.geometry_w + Math.Min(modal.trapezoid_delta_b, 0), modal.geometry_y + modal.geometry_h);
                pa[2] = new GeoLibPoint(modal.geometry_x + modal.geometry_w - Math.Max(modal.trapezoid_delta_b, 0), modal.geometry_y);
                pa[3] = new GeoLibPoint(modal.geometry_x - Math.Min(modal.trapezoid_delta_a, 0), modal.geometry_y);
            }

            pa[4] = new GeoLibPoint(pa[0]);

            cell_.addPolygon(pa, modal.layer, modal.datatype);
        }

        void addCtrapezoid()
        {
            GeoLibPoint[] pa = new GeoLibPoint[5];
            Int32[,] coords = new Int32[4, 4];

            switch (modal.ctrapezoid_type)
            {
                case 0:
                    coords = new [,] {
                                                { 0, 0, 0, 0 },  // x=0*w+0*h, y=0*w+0*h ...
                                                { 0, 0, 0, 1 },
                                                { 1, -1, 0, 1 },
                                                { 1, 0, 0, 0 }
                        };
                    break;

                case 1:
                    coords = new [,] {
                                                { 0, 0, 0, 0 },
                                                { 0, 0, 0, 1 },
                                                { 1, 0, 0, 1 },
                                                { 1, -1, 0, 0 }
                        };
                    break;

                case 2:
                    coords = new [,] {
                                                { 0, 0, 0, 0 },
                                                { 0, 1, 0, 1 },
                                                { 1, 0, 0, 1 },
                                                { 1, 0, 0, 0 }
                        };
                    break;

                case 3:
                    coords = new [,] {
                                                { 0, 1, 0, 0 },
                                                { 0, 0, 0, 1 },
                                                { 1, 0, 0, 1 },
                                                { 1, 0, 0, 0 }
                        };
                    break;

                case 4:
                    coords = new [,] {
                                                { 0, 0, 0, 0 },
                                                { 0, 1, 0, 1 },
                                                { 1, -1, 0, 1 },
                                                { 1, 0, 0, 0 }
                        };
                    break;

                case 5:
                    coords = new [,] {
                                                { 0, 1, 0, 0 },
                                                { 0, 0, 0, 1 },
                                                { 1, 0, 0, 1 },
                                                { 1, -1, 0, 0 }
                        };
                    break;

                case 6:
                    coords = new [,] {
                                                { 0, 0, 0, 0 },
                                                { 0, 1, 0, 1 },
                                                { 1, 0, 0, 1 },
                                                { 1, -1, 0, 0 }
                        };
                    break;

                case 7:
                    coords = new [,] {
                                                { 0, 1, 0, 0 },
                                                { 0, 0, 0, 1 },
                                                { 1, -1, 0, 1 },
                                                { 1, 0, 0, 0 }
                        };
                    break;

                case 8:
                    coords = new [,] {
                                                { 0, 0, 0, 0 },
                                                { 0, 0, 0, 1 },
                                                { 1, 0, -1, 1 },
                                                { 1, 0, 0, 0 }
                        };
                    break;

                case 9:
                    coords = new [,] {
                                                { 0, 0, 0, 0 },
                                                { 0, 0, -1, 1 },
                                                { 1, 0, 0, 1 },
                                                { 1, 0, 0, 0 }
                        };
                    break;

                case 10:
                    coords = new [,] {
                                                { 0, 0, 0, 0 },
                                                { 0, 0, 0, 1 },
                                                { 1, 0, 0, 1 },
                                                { 1, 0, 1, 0 }
                        };
                    break;

                case 11:
                    coords = new [,] {
                                                { 0, 0, 1, 0 },
                                                { 0, 0, 0, 1 },
                                                { 1, 0, 0, 1 },
                                                { 1, 0, 0, 0 }
                        };
                    break;

                case 12:
                    coords = new [,] {
                                                { 0, 0, 0, 0 },
                                                { 0, 0, 0, 1 },
                                                { 1, 0, -1, 1 },
                                                { 1, 0, 1, 0 }
                        };
                    break;

                case 13:
                    coords = new [,] {
                                                { 0, 0, 1, 0 },
                                                { 0, 0, -1, 1 },
                                                { 1, 0, 0, 1 },
                                                { 1, 0, 0, 0 }
                        };
                    break;

                case 14:
                    coords = new [,] {
                                                { 0, 0, 0, 0 },
                                                { 0, 0, -1, 1 },
                                                { 1, 0, 0, 1 },
                                                { 1, 0, 1, 0 }
                        };
                    break;

                case 15:
                    coords = new [,] {
                                                { 0, 0, 1, 0 },
                                                { 0, 0, 0, 1 },
                                                { 1, 0, -1, 1 },
                                                { 1, 0, 0, 0 }
                        };
                    break;

                case 16:
                    coords = new [,] {
                                                { 0, 0, 0, 0 },
                                                { 0, 0, 1, 0 },
                                                { 1, 0, 0, 0 },
                                                { 0, 0, 0, 0 }
                        };
                    break;

                case 17:
                    coords = new [,] {
                                                { 0, 0, 0, 0 },
                                                { 0, 0, 1, 0 },
                                                { 1, 0, 1, 0 },
                                                { 0, 0, 0, 0 }
                        };
                    break;

                case 18:
                    coords = new [,] {
                                                { 0, 0, 0, 0 },
                                                { 1, 0, 1, 0 },
                                                { 1, 0, 0, 0 },
                                                { 0, 0, 0, 0 }
                        };
                    break;

                case 19:
                    coords = new [,] {
                                                { 0, 0, 1, 0 },
                                                { 1, 0, 1, 0 },
                                                { 1, 0, 0, 0 },
                                                { 0, 0, 1, 0 }
                        };
                    break;

                case 20:
                    coords = new [,] {
                                                { 0, 0, 0, 0 },
                                                { 0, 1, 0, 1 },
                                                { 0, 2, 0, 0 },
                                                { 0, 0, 0, 0 }
                        };
                    break;

                case 21:
                    coords = new [,] {
                                                { 0, 0, 0, 1 },
                                                { 0, 2, 0, 1 },
                                                { 0, 1, 0, 0 },
                                                { 0, 0, 0, 1 }
                        };
                    break;

                case 22:
                    coords = new [,] {
                                                { 0, 0, 0, 0 },
                                                { 0, 0, 2, 0 },
                                                { 1, 0, 1, 0 },
                                                { 0, 0, 0, 0 }
                        };
                    break;

                case 23:
                    coords = new [,] {
                                                { 1, 0, 0, 0 },
                                                { 0, 0, 1, 0 },
                                                { 1, 0, 2, 0 },
                                                { 1, 0, 0, 0 }
                        };
                    break;

                case 24:
                    coords = new [,] {
                                                { 0, 0, 0, 0 },
                                                { 0, 0, 0, 1 },
                                                { 1, 0, 0, 1 },
                                                { 1, 0, 0, 0 }
                        };
                    break;

                case 25:
                    coords = new [,] {
                                                { 0, 0, 0, 0 },
                                                { 0, 0, 1, 0 },
                                                { 1, 0, 1, 0 },
                                                { 1, 0, 0, 0 }
                        };
                    break;

            }

            Int32 x, y;

            for (Int32 pt = 0; pt < 4; pt++)
            {
                x = 0;
                if (coords[pt, 0] != 0)
                {
                    x += coords[pt, 0] * modal.geometry_w;
                }

                if (coords[pt, 1] != 0)
                {
                    x += coords[pt, 1] * modal.geometry_h;
                }

                y = 0;
                if (coords[pt, 2] != 0)
                {
                    y += coords[pt, 2] * modal.geometry_w;
                }

                if (coords[pt, 3] != 0)
                {
                    y += coords[pt, 3] * modal.geometry_h;
                }

                pa[pt] = new GeoLibPoint(modal.geometry_x + x, modal.geometry_y + y);

                if (x > modal.geometry_w)
                {
                    modal.geometry_w = x;
                }
                if (y > modal.geometry_h)
                {
                    modal.geometry_h = y;
                }
            }

            pa[^1] = new GeoLibPoint(pa[0]);

            cell_.addPolygon(pa, modal.layer, modal.datatype);
        }

        void processRepetition(elementType e)
        {
            if ((modal.repetition <= 3) && (e == elementType.cellrefElement))
            {
                switch (modal.repetition)
                {
                    case 1:
                        {
                            cell_.addCellrefArray(drawing_.findCell(modal.placement_cell), new GeoLibPoint(modal.placement_x, modal.placement_y),
                                   new GeoLibPoint(modal.x_space + modal.placement_x, modal.y_space + modal.placement_y), modal.x_dimension, modal.y_dimension);
                            cell_.elementList[^1].setName(modal.placement_cell);
                            cell_.elementList[^1].rotate(modal.angle);
                            cell_.elementList[^1].scale(modal.mag);
                            if (modal.mirror_x)
                            {
                                cell_.elementList[^1].setMirrorx();
                            }
                        }
                        break;
                    case 2:
                        {
                            cell_.addCellrefArray(drawing_.findCell(modal.placement_cell), new GeoLibPoint(modal.placement_x, modal.placement_y),
                                   new GeoLibPoint(modal.x_space + modal.placement_x, modal.placement_y), modal.x_dimension, 1);
                            cell_.elementList[^1].setName(modal.placement_cell);
                            cell_.elementList[^1].rotate(modal.angle);
                            cell_.elementList[^1].scale(modal.mag);
                            if (modal.mirror_x)
                            {
                                cell_.elementList[^1].setMirrorx();
                            }
                        }
                        break;
                    case 3:
                        {
                            cell_.addCellrefArray(drawing_.findCell(modal.placement_cell), new GeoLibPoint(modal.placement_x, modal.placement_y),
                                   new GeoLibPoint(modal.placement_x, modal.y_space + modal.placement_y), 1, modal.y_dimension);
                            cell_.elementList[^1].setName(modal.placement_cell);
                            cell_.elementList[^1].rotate(modal.angle);
                            cell_.elementList[^1].scale(modal.mag);
                            if (modal.mirror_x)
                            {
                                cell_.elementList[^1].setMirrorx();
                            }
                        }
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
                                addElement(e, new GeoLibPoint(x * modal.x_space, y * modal.y_space));
                            }
                        break;
                    case 2:
                        for (int x = 0; x < modal.x_dimension; x++)
                        {
                            addElement(e, new GeoLibPoint(x * modal.x_space, 0));
                        }
                        break;
                    case 3:
                        for (int y = 0; y < modal.y_dimension; y++)
                        {
                            addElement(e, new GeoLibPoint(0, y * modal.y_space));
                        }
                        break;
                    default:
                        for (int x = 0; x < modal.repArray.Count; x++)
                        {
                            addElement(e, modal.repArray[x]);
                        }
                        break;
                }
            }
        }
    }
}
