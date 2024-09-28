using Clipper2Lib;

namespace shapeEngine;

public static class shapeOffsets
{
    private static PointD rectangle_offset(ShapeSettings shapeSettings)
    {
        double xOffset = 0;
        double yOffset = 0;
        string posInSubShapeString = ((ShapeSettings.subShapeLocations)shapeSettings.getInt(ShapeSettings.properties_i.posInSubShapeIndex)).ToString();
        double tmp_xOffset = 0;
        double tmp_yOffset = 0;

        if (posInSubShapeString is "TL" or "TR" or "TS" or "RS" or "LS" or "C")
        {
            // Vertical offset needed to put reference corner at world center
            tmp_yOffset = -Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 0));

            // Half the value for a vertical centering requirement
            if (posInSubShapeString is "RS" or "LS" or "C")
            {
                tmp_yOffset = Convert.ToDouble(tmp_yOffset / 2);
            }
        }
        yOffset -= tmp_yOffset;

        if (posInSubShapeString is "TR" or "BR" or "TS" or "RS" or "BS" or "C")
        {
            tmp_xOffset -= Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0));

            // Half the value for horizontal centering conditions
            if (posInSubShapeString is "TS" or "BS" or "C")
            {
                tmp_xOffset = Convert.ToDouble(tmp_xOffset / 2);
            }
        }
        xOffset += tmp_xOffset;

        return new PointD(xOffset, yOffset);
    }

    private static PointD lShape_offset(ShapeSettings shapeSettings)
    {
        double xOffset = 0;
        double yOffset = 0;
        string posInSubShapeString = ((ShapeSettings.subShapeLocations)shapeSettings.getInt(ShapeSettings.properties_i.posInSubShapeIndex)).ToString();
        double tmp_xOffset = 0;
        double tmp_yOffset = 0;

        if (posInSubShapeString is "TL" or "TR" or "TS" or "RS" or "LS" or "C")
        {
            // Vertical offset needed to put reference corner at world center
            if (shapeSettings.getInt(ShapeSettings.properties_i.subShapeRefIndex) == 0)
            {
                tmp_yOffset = -Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 0));
            }
            else
            {
                tmp_yOffset = -Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 1));
            }

            // Half the value for a vertical centering requirement
            if (posInSubShapeString is "RS" or "LS" or "C")
            {
                tmp_yOffset = Convert.ToDouble(tmp_yOffset / 2);
            }
        }
        yOffset -= tmp_yOffset;

        if (shapeSettings.getInt(ShapeSettings.properties_i.subShapeRefIndex) == 1 && posInSubShapeString is "LS" or "BL" or "TL")
        {
            tmp_xOffset -= Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0)); // essentially the same in X as the RS for subshape 1.
        }
        else
        {
            if (posInSubShapeString is "TR" or "BR" or "TS" or "RS" or "BS" or "C")
            {
                if (shapeSettings.getInt(ShapeSettings.properties_i.subShapeRefIndex) == 0)
                {
                    tmp_xOffset -= Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0));
                }
                else
                {
                    tmp_xOffset -= Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0));
                    tmp_xOffset -= Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 1));
                }

                // Half the value for horizontal centering conditions
                if (posInSubShapeString is "TS" or "BS" or "C")
                {
                    if (shapeSettings.getInt(ShapeSettings.properties_i.subShapeRefIndex) == 0)
                    {
                        tmp_xOffset = Convert.ToDouble(tmp_xOffset / 2);
                    }
                    else
                    {
                        tmp_xOffset += Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 1) / 2);
                    }
                }
            }
        }

        xOffset += tmp_xOffset;
        
        return new PointD(xOffset, yOffset);
    }

    private static PointD tShape_offset(ShapeSettings shapeSettings)
    {
        double xOffset = 0;
        double yOffset = 0;
        string posInSubShapeString = ((ShapeSettings.subShapeLocations)shapeSettings.getInt(ShapeSettings.properties_i.posInSubShapeIndex)).ToString();
        double tmp_xOffset = 0;
        double tmp_yOffset = 0;

        if (shapeSettings.getInt(ShapeSettings.properties_i.subShapeRefIndex) == 1 && posInSubShapeString is "BR" or "BL" or "BS")
        {
            tmp_yOffset -= Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.verOffset, 1));
        }
        else
        {
            if (posInSubShapeString is "TL" or "TR" or "TS" or "RS" or "LS" or "C")
            {
                if (shapeSettings.getInt(ShapeSettings.properties_i.subShapeRefIndex) == 0)
                {
                    tmp_yOffset = -Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 0));
                    // Half the value for a vertical centering requirement
                    if (posInSubShapeString is "RS" or "LS" or "C")
                    {
                        tmp_yOffset = Convert.ToDouble(tmp_yOffset / 2);
                    }
                }
                else
                {
                    tmp_yOffset = -Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 1));
                    // Half the value for a vertical centering requirement
                    if (posInSubShapeString is "RS" or "LS" or "C")
                    {
                        tmp_yOffset = Convert.ToDouble(tmp_yOffset / 2);
                    }
                    tmp_yOffset -= Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.verOffset, 1));
                }

            }
        }
        yOffset -= tmp_yOffset;

        if (shapeSettings.getInt(ShapeSettings.properties_i.subShapeRefIndex) == 1 && posInSubShapeString is "LS" or "BL" or "TL")
        {
            tmp_xOffset -= Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0)); // essentially the same in X as the RS for subshape 1.
        }
        else
        {
            if (posInSubShapeString is "TR" or "BR" or "TS" or "RS" or "BS" or "C")
            {
                if (shapeSettings.getInt(ShapeSettings.properties_i.subShapeRefIndex) == 0)
                {
                    tmp_xOffset -= Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0));
                }
                else
                {
                    tmp_xOffset -= Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0));
                    tmp_xOffset -= Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 1));
                }

                // Half the value for horizontal centering conditions
                if (posInSubShapeString is "TS" or "BS" or "C")
                {
                    if (shapeSettings.getInt(ShapeSettings.properties_i.subShapeRefIndex) == 0)
                    {
                        tmp_xOffset = Convert.ToDouble(tmp_xOffset / 2);
                    }
                    else
                    {
                        tmp_xOffset += Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 1) / 2);
                    }
                }
            }
        }

        xOffset += tmp_xOffset;
        return new PointD(xOffset, yOffset);
    }

    private static PointD xShape_offset(ShapeSettings shapeSettings)
    {
        double xOffset = 0;
        double yOffset = 0;
        string posInSubShapeString = ((ShapeSettings.subShapeLocations)shapeSettings.getInt(ShapeSettings.properties_i.posInSubShapeIndex)).ToString();
        double tmp_xOffset = 0;
        double tmp_yOffset = 0;

        if (shapeSettings.getInt(ShapeSettings.properties_i.subShapeRefIndex) == 1 && posInSubShapeString is "BR" or "BL" or "BS")
        {
            tmp_yOffset -= Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.verOffset, 1));
        }
        else
        {
            if (posInSubShapeString is "TL" or "TR" or "TS" or "RS" or "LS" or "C")
            {
                // Vertical offset needed to put reference corner at world center
                if (shapeSettings.getInt(ShapeSettings.properties_i.subShapeRefIndex) == 0)
                {
                    tmp_yOffset = -Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 0));
                    // Half the value for a vertical centering requirement
                    if (posInSubShapeString is "RS" or "LS" or "C")
                    {
                        tmp_yOffset = Convert.ToDouble(tmp_yOffset / 2);
                    }
                }
                else
                {
                    tmp_yOffset = -Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 1));
                    // Half the value for a vertical centering requirement
                    if (posInSubShapeString is "RS" or "LS" or "C")
                    {
                        tmp_yOffset = Convert.ToDouble(tmp_yOffset / 2);
                    }
                    tmp_yOffset -= Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.verOffset, 1));
                }

            }
        }
        yOffset -= tmp_yOffset;

        if (shapeSettings.getInt(ShapeSettings.properties_i.subShapeRefIndex) == 1 && posInSubShapeString is "LS" or "BL" or "TL")
        {
            tmp_xOffset -= Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.horOffset, 1));
        }
        else
        {
            if (posInSubShapeString is "TR" or "BR" or "TS" or "RS" or "BS" or "C")
            {
                if (shapeSettings.getInt(ShapeSettings.properties_i.subShapeRefIndex) == 0)
                {
                    tmp_xOffset -= Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0));
                }
                else
                {
                    tmp_xOffset -= Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 1));
                }

                // Half the value for horizontal centering conditions
                if (posInSubShapeString is "TS" or "BS" or "C")
                {
                    if (shapeSettings.getInt(ShapeSettings.properties_i.subShapeRefIndex) == 0)
                    {
                        tmp_xOffset = Convert.ToDouble(tmp_xOffset / 2);
                    }
                    else
                    {
                        tmp_xOffset += Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 1) / 2);
                    }
                }

                if (shapeSettings.getInt(ShapeSettings.properties_i.subShapeRefIndex) == 1)
                {
                    tmp_xOffset -= Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.horOffset, 1));
                }
            }
        }

        xOffset += tmp_xOffset;
        return new PointD(xOffset, yOffset);
    }

    private static PointD uShape_offset(ShapeSettings shapeSettings)
    {
        double xOffset = 0;
        double yOffset = 0;
        string posInSubShapeString = ((ShapeSettings.subShapeLocations)shapeSettings.getInt(ShapeSettings.properties_i.posInSubShapeIndex)).ToString();
        double tmp_xOffset = 0;
        double tmp_yOffset = 0;

        if (shapeSettings.getInt(ShapeSettings.properties_i.subShapeRefIndex) == 0)
        {
            if (posInSubShapeString is "TL" or "TR" or "TS" or "RS" or "LS" or "C")
            {
                tmp_yOffset = -Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 0));

                // Half the value for a vertical centering requirement
                if (posInSubShapeString is "RS" or "LS" or "C")
                {
                    tmp_yOffset = Convert.ToDouble(tmp_yOffset / 2);
                }
            }
            yOffset -= tmp_yOffset;

            if (posInSubShapeString is "TR" or "BR" or "TS" or "RS" or "BS" or "C")
            {
                tmp_xOffset -= Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0));

                // Half the value for horizontal centering conditions
                if (posInSubShapeString is "TS" or "BS" or "C")
                {
                    tmp_xOffset = Convert.ToDouble(tmp_xOffset / 2);
                }
            }
        }
        else
        {
            // Subshape 2 is always docked against top edge of subshape 1 in U.
            if (posInSubShapeString is "TL" or "TR" or "TS" or "RS" or "LS" or "BL" or "BR" or "BS" or "C")
            {
                tmp_yOffset = -Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 0));

                switch (posInSubShapeString)
                {
                    // Half the value for a vertical centering requirement
                    case "RS" or "LS" or "C":
                        tmp_yOffset += Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 1) / 2);
                        break;
                    // Subtract the value for a subshape 2 bottom edge requirement
                    case "BL" or "BR" or "BS":
                        tmp_yOffset += Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 1));
                        break;
                }
            }
            yOffset -= tmp_yOffset;

            // Subshape 2 is always H-centered in U. Makes it easy.
            tmp_xOffset -= Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0) / 2);

            switch (posInSubShapeString)
            {
                case "TR" or "BR" or "RS":
                    tmp_xOffset -= Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 1) / 2);
                    break;
                case "TL" or "BL" or "LS":
                    tmp_xOffset += Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 1) / 2);
                    break;
            }
        }
        xOffset += tmp_xOffset;
        return new PointD(xOffset, yOffset);
    }

    private static PointD sShape_offset(ShapeSettings shapeSettings)
    {
        double xOffset = 0;
        double yOffset = 0;
        string posInSubShapeString = ((ShapeSettings.subShapeLocations)shapeSettings.getInt(ShapeSettings.properties_i.posInSubShapeIndex)).ToString();
        double tmp_xOffset = 0;
        double tmp_yOffset = 0;

        switch (shapeSettings.getInt(ShapeSettings.properties_i.subShapeRefIndex))
        {
            case 0:
                if (posInSubShapeString is "TL" or "TR" or "TS" or "RS" or "LS" or "C")
                {
                    tmp_yOffset = -Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 0));

                    // Half the value for a vertical centering requirement
                    if (posInSubShapeString is "RS" or "LS" or "C")
                    {
                        tmp_yOffset = Convert.ToDouble(tmp_yOffset / 2);
                    }
                }

                if (posInSubShapeString is "TR" or "BR" or "TS" or "RS" or "BS" or "C")
                {
                    tmp_xOffset -= Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0));

                    // Half the value for horizontal centering conditions
                    if (posInSubShapeString is "TS" or "BS" or "C")
                    {
                        tmp_xOffset = Convert.ToDouble(tmp_xOffset / 2);
                    }
                }
                break;

            case 1:
                // Subshape 2 is always vertically offset relative to bottom edge of subshape 1 in S.
                if (posInSubShapeString is "TL" or "TR" or "TS" or "RS" or "LS" or "BL" or "BR" or "BS" or "C")
                {
                    tmp_yOffset -= Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.verOffset, 1));

                    switch (posInSubShapeString)
                    {
                        // Half the value for a vertical centering requirement
                        case "RS" or "LS" or "C":
                            tmp_yOffset -= Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 1) / 2);
                            break;
                        // Subtract the value for a subshape 2 bottom edge requirement
                        case "TL" or "TR" or "TS":
                            tmp_yOffset -= Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 1));
                            break;
                    }
                }

                // Subshape 2 is always pinned to left edge in S. Makes it easy.

                if (posInSubShapeString is "TR" or "BR" or "RS" or "TS" or "BS" or "C")
                {
                    tmp_xOffset -= Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 1));
                    if (posInSubShapeString is "TS" or "C" or "BS")
                    {
                        tmp_xOffset /= 2;
                    }
                }

                break;

            case 2:
                tmp_yOffset -= Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 0));
                // Subshape 3 is always offset relative to top edge of subshape 1 in S.
                if (posInSubShapeString is "TL" or "TR" or "TS" or "RS" or "LS" or "BL" or "BR" or "BS" or "C")
                {
                    tmp_yOffset += Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.verOffset, 2));

                    switch (posInSubShapeString)
                    {
                        // Half the value for a vertical centering requirement
                        case "RS" or "LS" or "C":
                            tmp_yOffset += Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 2) / 2);
                            break;
                        // Subtract the value for a subshape 2 bottom edge requirement
                        case "BL" or "BR" or "BS":
                            tmp_yOffset += Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 2));
                            break;
                    }
                }

                // Subshape 3 is always pinned to right edge in S. Makes it easy.
                tmp_xOffset -= Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0));

                switch (posInSubShapeString)
                {
                    case "TL" or "BL" or "LS":
                        tmp_xOffset += Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 2));
                        break;
                    case "TS" or "BS" or "C":
                        tmp_xOffset += Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 2) / 2);
                        break;
                }

                break;
        }

        yOffset -= tmp_yOffset;
        xOffset += tmp_xOffset;
        return new PointD(xOffset, yOffset);
    }

    public static PointD doOffsets(int mode, ShapeSettings shapeSettings)
    {
        // Use our shape-specific offset calculation methods :
        PointD offset = new(0, 0);

        switch (mode)
        {
            case 0:
                switch (shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex))
                {
                    case (int)ShapeSettings.typeShapes_mode0.rectangle:
                        offset = rectangle_offset(shapeSettings);
                        break;
                    case (int)ShapeSettings.typeShapes_mode0.L:
                        offset = lShape_offset(shapeSettings);
                        break;
                    case (int)ShapeSettings.typeShapes_mode0.T:
                        offset = tShape_offset(shapeSettings);
                        break;
                    case (int)ShapeSettings.typeShapes_mode0.X:
                        offset = xShape_offset(shapeSettings);
                        break;
                    case (int)ShapeSettings.typeShapes_mode0.U:
                        offset = uShape_offset(shapeSettings);
                        break;
                    case (int)ShapeSettings.typeShapes_mode0.S:
                        offset = sShape_offset(shapeSettings);
                        break;
                    /*
                    case (int)ShapeSettings.typeShapes_mode0.BOOLEAN:
                    case (int)ShapeSettings.typeShapes_mode0.GEOCORE:
                    */
                    default:
                        // customShape_offset(entropyLayerSettings);
                        break;
                }

                break;
            case 1:
                break;
        }

        // Now for global offset.
        offset.x += Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.gHorOffset));
        offset.y -= Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.gVerOffset));

        return offset;
    }

}