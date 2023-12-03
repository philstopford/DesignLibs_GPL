using geoWrangler;

namespace shapeEngine;

[Serializable]
public class ShapeSettings
{
    public static List<string> getAvailableTipsLocations()
    {
        return pGetAvailableTipsLocations();
    }

    private static List<string> pGetAvailableTipsLocations()
    {
        return availableTipsLocations;
    }

    public enum tipLocations { none, L, R, LR, T, B, TB, TL, TR, TLR, BL, BR, BLR, TBL, TBR, all }
    
    public static List<string> getAvailableSubShapePositions()
    {
        return pGetAvailableSubShapePositions();
    }

    private static List<string> pGetAvailableSubShapePositions()
    {
        return availableSubShapePositions;
    }

    public static List<string> getAvailableHorShapePositions()
    {
        return pGetAvailableHorShapePositions();
    }

    private static List<string> pGetAvailableHorShapePositions()
    {
        return availableHorShapePositions;
    }

    public static List<string> getAvailableVerShapePositions()
    {
        return pGetAvailableVerShapePositions();
    }

    private static List<string> pGetAvailableVerShapePositions()
    {
        return availableVerShapePositions;
    }

    public enum subShapeLocations { TL, TR, BL, BR, TS, RS, BS, LS, C }
    
    public static List<string> getPolyFillTypes()
    {
        return pGetPolyFillTypes();
    }

    private static List<string> pGetPolyFillTypes()
    {
        return polyFillTypes;
    }

    private static readonly List<string> availableHorShapePositions = new() { "Left", "Middle", "Right" };
    private static readonly List<string> availableVerShapePositions = new() { "Bottom", "Middle", "Top" };
    public enum subShapeHorLocs { L, M, R }
    public enum subShapeVerLocs { B, M, T }

    
    private static readonly List<string> availableSubShapePositions = new()
    { "Corner: Top Left", "Corner: Top Right", "Corner: Bottom Left", "Corner: Bottom Right",
        "Middle: Top Side", "Middle: Right Side", "Middle: Bottom Side", "Middle: Left Side",
        "Center"};

    private static readonly List<string> availableTipsLocations = new()
    { "None", "Left", "Right", "Left & Right",
        "Top", "Bottom", "Top & Bottom",
        "Top & Left", "Top & Right", "Top & Left & Right",
        "Bottom & Left", "Bottom & Right", "Bottom & Left & Right",
        "Top & Bottom & Left", "Top & Bottom & Right",
        "All"};

    private static readonly List<string> polyFillTypes = new() { "Even/Odd", "Non-zero", "Positive", "Negative" };

    public enum PolyFill { pftEvenOdd, pftNonZero, pftPositive, pftNegative }

    public enum typeShapes_mode0 { none, rectangle, L, T, X, U, S, GEOCORE, BOOLEAN }
    public enum typeShapes_mode1 { none, rectangle, L, T, X, U, S, text, bounding, complex }

    private static readonly int default_subShapeTipLocIndex = 0;
    private static readonly int default_subShape2TipLocIndex = 0;
    private static readonly int default_subShape3TipLocIndex = 0;
    
    private static readonly int default_enabled = 0;

    private static readonly int default_shapeIndex = (int)typeShapes_mode0.none;
    private static readonly int default_subShapeRefIndex = 0;
    private static readonly int default_posInSubShapeIndex = (int)subShapeLocations.BL;

    private static readonly int default_edgeSlide = 0;
    private static readonly int default_proximitySideRays = 2;
    private static readonly int default_proximitySideRaysFallOff = 0;

    private static readonly decimal default_LWR = 0;
    private static readonly int default_LWRNoiseType = (int)NoiseC.noiseIndex.perlin;
    private static readonly decimal default_LWRNoiseFreq = 0.2m;

    private static readonly int default_flipH = 0;
    private static readonly int default_flipV = 0;
    private static readonly int default_alignGeom = 0;

    private int enabled = default_enabled;
    private int shapeIndex = default_shapeIndex;
    private int subShapeTipLocIndex = default_subShapeTipLocIndex;
    private int subShape2TipLocIndex = default_subShape2TipLocIndex;
    private int subShape3TipLocIndex = default_subShape3TipLocIndex;
    private int subShapeRefIndex = default_subShapeRefIndex;
    private int posInSubShapeIndex = default_posInSubShapeIndex;
    private int edgeSlide = default_edgeSlide;
    private int proximitySideRays = default_proximitySideRays;
    private int proxSideRaysFallOff = default_proximitySideRaysFallOff;
    private int LWRNoiseType = default_LWRNoiseType;
    private int LWR2NoiseType = default_LWRNoiseType;
    private int flipH = default_flipH;
    private int flipV = default_flipV;
    private int alignGeomX = default_alignGeom;
    private int alignGeomY = default_alignGeom;

    private static string default_layerName = "";
    
    public enum properties_i
    {
        enabled,
        shapeIndex,
        subShapeTipLocIndex, subShape2TipLocIndex, subShape3TipLocIndex,
        subShapeRefIndex,posInSubShapeIndex,
        edgeSlide,
        proxRays,proxSideRaysFallOff,
        lwrType, lwr2Type,
        flipH, flipV, alignX, alignY,
    }

    public int getInt(properties_i p, int _subShapeRef = -1)
    {
        return pGetInt(p, _subShapeRef);
    }

    int pGetInt(properties_i p, int _subShapeRef)
    {
        int ret = 0;
        switch (p)
        {
            case properties_i.enabled:
                ret = enabled;
                break;
            case properties_i.shapeIndex:
                ret = shapeIndex;
                break;
            case properties_i.subShapeTipLocIndex:
                ret = subShapeTipLocIndex;
                break;
            case properties_i.subShape2TipLocIndex:
                ret = subShape2TipLocIndex;
                break;
            case properties_i.subShape3TipLocIndex:
                ret = subShape3TipLocIndex;
                break;
            case properties_i.subShapeRefIndex:
                ret = subShapeRefIndex;
                break;
            case properties_i.posInSubShapeIndex:
                ret = posInSubShapeIndex;
                break;
            case properties_i.edgeSlide:
                ret = edgeSlide;
                break;
            case properties_i.proxRays:
                ret = proximitySideRays;
                break;
            case properties_i.proxSideRaysFallOff:
                ret = proxSideRaysFallOff;
                break;
            case properties_i.lwrType:
                ret = LWRNoiseType;
                break;
            case properties_i.lwr2Type:
                ret = LWR2NoiseType;
                break;
            case properties_i.flipH:
                ret = flipH;
                break;
            case properties_i.flipV:
                ret = flipV;
                break;
            case properties_i.alignX:
                ret = alignGeomX;
                break;
            case properties_i.alignY:
                ret = alignGeomY;
                break;
        }

        return ret;
    }

    public void setInt(properties_i p, int val, int _subShapeRef = -1)
    {
        pSetInt(p, val, _subShapeRef);
    }

    void pSetInt(properties_i p, int val, int _subShapeRef)
    {
        switch (p)
        {
            case properties_i.enabled:
                enabled = val;
                break;
            case properties_i.shapeIndex:
                shapeIndex = val;
                break;
            case properties_i.subShapeTipLocIndex:
                subShapeTipLocIndex = val;
                break;
            case properties_i.subShape2TipLocIndex:
                subShape2TipLocIndex = val;
                break;
            case properties_i.subShape3TipLocIndex:
                subShape3TipLocIndex = val;
                break;
            case properties_i.subShapeRefIndex:
                subShapeRefIndex = val;
                break;
            case properties_i.posInSubShapeIndex:
                posInSubShapeIndex = val;
                break;
            case properties_i.edgeSlide:
                edgeSlide = val;
                break;
            case properties_i.proxRays:
                proximitySideRays = val;
                break;
            case properties_i.proxSideRaysFallOff:
                proxSideRaysFallOff = val;
                break;
            case properties_i.lwrType:
                LWRNoiseType = val;
                break;
            case properties_i.lwr2Type:
                LWR2NoiseType = val;
                break;
            case properties_i.flipH:
                flipH = val;
                break;
            case properties_i.flipV:
                flipV = val;
                break;
            case properties_i.alignX:
                alignGeomX = val;
                break;
            case properties_i.alignY:
                alignGeomY = val;
                break;
        }
    }

    public void defaultInt(properties_i p)
    {
        pDefaultInt(p);
    }

    private void pDefaultInt(properties_i p)
    {
        switch (p)
        {
            case properties_i.enabled:
                enabled = default_enabled;
                break;
            case properties_i.shapeIndex:
                shapeIndex = default_shapeIndex;
                break;
            case properties_i.subShapeTipLocIndex:
                subShapeTipLocIndex = default_subShapeTipLocIndex;
                break;
            case properties_i.subShape2TipLocIndex:
                subShape2TipLocIndex = default_subShape2TipLocIndex;
                break;
            case properties_i.subShape3TipLocIndex:
                subShape3TipLocIndex = default_subShape3TipLocIndex;
                break;
            case properties_i.subShapeRefIndex:
                subShapeRefIndex = default_subShapeRefIndex;
                break;
            case properties_i.posInSubShapeIndex:
                posInSubShapeIndex = default_posInSubShapeIndex;
                break;
            case properties_i.edgeSlide:
                edgeSlide = default_edgeSlide;
                break;
            case properties_i.proxRays:
                proximitySideRays = default_proximitySideRays;
                break;
            case properties_i.proxSideRaysFallOff:
                proxSideRaysFallOff = default_proximitySideRaysFallOff;
                break;
            case properties_i.lwrType:
                LWRNoiseType = default_LWRNoiseType;
                break;
            case properties_i.lwr2Type:
                LWR2NoiseType = default_LWRNoiseType;
                break;
            case properties_i.flipH:
                flipH = default_flipH;
                break;
            case properties_i.flipV:
                flipV = default_flipV;
                break;
            case properties_i.alignX:
                alignGeomX = default_alignGeom;
                break;
            case properties_i.alignY:
                alignGeomY = default_alignGeom;
                break;
        }
    }

    public static int getDefaultInt(properties_i p)
    {
        return pGetDefaultInt(p);
    }

    private static int pGetDefaultInt(properties_i p)
    {
        int val = 0;
        switch (p)
        {
            case properties_i.enabled:
                val = default_enabled;
                break;
            case properties_i.shapeIndex:
                val = default_shapeIndex;
                break;
            case properties_i.subShapeTipLocIndex:
                val = default_subShapeTipLocIndex;
                break;
            case properties_i.subShape2TipLocIndex:
                val = default_subShape2TipLocIndex;
                break;
            case properties_i.subShape3TipLocIndex:
                val = default_subShape3TipLocIndex;
                break;
            case properties_i.subShapeRefIndex:
                val = default_subShapeRefIndex;
                break;
            case properties_i.posInSubShapeIndex:
                val = default_posInSubShapeIndex;
                break;
            case properties_i.edgeSlide:
                val = default_edgeSlide;
                break;
            case properties_i.proxRays:
                val = default_proximitySideRays;
                break;
            case properties_i.proxSideRaysFallOff:
                val = default_proximitySideRaysFallOff;
                break;
            case properties_i.lwrType:
                val = default_LWRNoiseType;
                break;
            case properties_i.lwr2Type:
                val = default_LWRNoiseType;
                break;
            case properties_i.flipH:
                val = default_flipH;
                break;
            case properties_i.flipV:
                val = default_flipV;
                break;
            case properties_i.alignX:
                val = default_alignGeom;
                break;
            case properties_i.alignY:
                val = default_alignGeom;
                break;
        }

        return val;
    }

    public static readonly decimal default_subShapeHorLength = 0;
    public static readonly decimal default_subShapeHorOffset = 0;
    public static readonly decimal default_subShapeVerLength = 0;
    public static readonly decimal default_subShapeVerOffset = 0;
    private static readonly decimal default_subShape2HorLength = 0;
    private static readonly decimal default_subShape2HorOffset = 0;
    private static readonly decimal default_subShape2VerLength = 0;
    private static readonly decimal default_subShape2VerOffset = 0;
    private static readonly decimal default_subShape3HorLength = 0;
    private static readonly decimal default_subShape3HorOffset = 0;
    private static readonly decimal default_subShape3VerLength = 0;
    private static readonly decimal default_subShape3VerOffset = 0;
    private static readonly decimal default_sideBias = 0;
    private static readonly decimal default_horTipBias = 0;
    private static readonly decimal default_verTipBias = 0;
    private static readonly decimal default_innerCRR = 0;
    private static readonly decimal default_outerCRR = 0;
    public static readonly decimal default_rotation = 0;
    private static readonly decimal default_proximityBias = 0;
    private static readonly decimal default_proximityIsoDistance = 0;
    private static readonly decimal default_proximitySideRaysFallOffMultiplier = 1.0m;
    private static readonly decimal default_edgeSlideTension = 0.35m;
    private static readonly decimal default_horTipBiasNVar = 0;
    private static readonly decimal default_horTipBiasPVar = 0;
    private static readonly decimal default_verTipBiasNVar = 0;
    private static readonly decimal default_verTipBiasPVar = 0;

    private static readonly decimal default_globalHorOffset = 0;
    private static readonly decimal default_globalVerOffset = 0;

    private static readonly decimal default_rayExtension = 1.03m;

    public decimal subShapeHorLength = default_subShapeHorLength;
    public decimal subShapeHorOffset = default_subShapeHorOffset;
    public decimal subShapeVerLength = default_subShapeVerLength;
    public decimal subShapeVerOffset = default_subShapeVerOffset;
    public decimal subShape2HorLength = default_subShape2HorLength;
    public decimal subShape2HorOffset = default_subShape2HorOffset;
    public decimal subShape2VerLength = default_subShape2VerLength;
    public decimal subShape2VerOffset = default_subShape2VerOffset;
    public decimal subShape3HorLength = default_subShape3HorLength;
    public decimal subShape3HorOffset = default_subShape3HorOffset;
    public decimal subShape3VerLength = default_subShape3VerLength;
    public decimal subShape3VerOffset = default_subShape3HorOffset;
    private decimal innerCRR = default_innerCRR;
    private decimal outerCRR = default_outerCRR;
    private decimal globalHorOffset = default_globalHorOffset;
    private decimal globalVerOffset = default_globalVerOffset;
    private decimal sideBias = default_sideBias;
    public decimal horTipBias = default_horTipBias;
    private decimal horTipBiasPVar = default_horTipBias;
    private decimal horTipBiasNVar;
    public decimal verTipBias;
    private decimal verTipBiasPVar;
    private decimal verTipBiasNVar;
    private decimal proximityBias = default_proximityBias;
    private decimal proximityIsoDistance = default_proximityIsoDistance;
    private decimal rotation = default_rotation;
    private decimal edgeSlideTension= default_edgeSlideTension;
    private decimal LWR;
    private decimal LWR2;
    private decimal LWRNoiseFreq = default_LWRNoiseFreq;
    private decimal LWR2NoiseFreq = default_LWRNoiseFreq;
    private decimal proxSideRaysMultiplier = default_proximitySideRaysFallOffMultiplier;
    private decimal rayExtension = default_rayExtension;
    private decimal gcRayExtension;

    public enum properties_decimal
    {
        horLength,
        verLength,
        horOffset,
        verOffset,
        gHorOffset, gVerOffset,
        iCR, oCR,
        sBias, hTBias, hTNVar, hTPVar, vTBias, vTNVar, vTPVar,
        lwr, lwrFreq, lwr2, lwr2Freq,
        eTension, rot,
        pBias, pBiasDist, proxSideRaysMultiplier,
        rayExtension, keyhole_factor
    }

    public decimal getDecimal(properties_decimal p, int _subShapeRef = -1)
    {
        return pGetDecimal(p, _subShapeRef);
    }

    decimal pGetDecimal(properties_decimal p, int _subShapeRef)
    {
        decimal ret = 0;
        switch (p)
        {
            case properties_decimal.horLength:
                switch (_subShapeRef)
                {
                    case 0:
                        ret = subShapeHorLength;
                        break;
                    case 1:
                        ret = subShape2HorLength;
                        break;
                    case 2:
                        ret = subShape3HorLength;
                        break;
                }
                break;
            case properties_decimal.verLength:
                switch (_subShapeRef)
                {
                    case 0:
                        ret = subShapeVerLength;
                        break;
                    case 1:
                        ret = subShape2VerLength;
                        break;
                    case 2:
                        ret = subShape3VerLength;
                        break;
                }
                break;
            case properties_decimal.horOffset:
                switch (_subShapeRef)
                {
                    case 0:
                        ret = subShapeHorOffset;
                        break;
                    case 1:
                        ret = subShape2HorOffset;
                        break;
                    case 2:
                        ret = subShape3HorOffset;
                        break;
                }
                break;
            case properties_decimal.verOffset:
                switch (_subShapeRef)
                {
                    case 0:
                        ret = subShapeVerOffset;
                        break;
                    case 1:
                        ret = subShape2VerOffset;
                        break;
                    case 2:
                        ret = subShape3VerOffset;
                        break;
                }
                break;
            case properties_decimal.gHorOffset:
                ret = globalHorOffset;
                break;
            case properties_decimal.gVerOffset:
                ret = globalVerOffset;
                break;
            case properties_decimal.iCR:
                ret = innerCRR;
                break;
            case properties_decimal.oCR:
                ret = outerCRR;
                break;
            case properties_decimal.sBias:
                ret = sideBias;
                break;
            case properties_decimal.hTBias:
                ret = horTipBias;
                break;
            case properties_decimal.hTNVar:
                ret = horTipBiasNVar;
                break;
            case properties_decimal.hTPVar:
                ret = horTipBiasPVar;
                break;
            case properties_decimal.vTBias:
                ret = verTipBias;
                break;
            case properties_decimal.vTNVar:
                ret = verTipBiasNVar;
                break;
            case properties_decimal.vTPVar:
                ret = verTipBiasPVar;
                break;
            case properties_decimal.lwr:
                ret = LWR;
                break;
            case properties_decimal.lwrFreq:
                ret = LWRNoiseFreq;
                break;
            case properties_decimal.lwr2:
                ret = LWR2;
                break;
            case properties_decimal.lwr2Freq:
                ret = LWR2NoiseFreq;
                break;
            case properties_decimal.eTension:
                ret = edgeSlideTension;
                break;
            case properties_decimal.rot:
                ret = rotation;
                break;
            case properties_decimal.pBias:
                ret = proximityBias;
                break;
            case properties_decimal.pBiasDist:
                ret = proximityIsoDistance;
                break;
            case properties_decimal.proxSideRaysMultiplier:
                ret = proxSideRaysMultiplier;
                break;
            case properties_decimal.rayExtension:
                ret = rayExtension;
                break;
            case properties_decimal.keyhole_factor:
                ret = gcRayExtension;
                break;
        }

        return ret;
    }

    public void setDecimal(properties_decimal p, decimal val, int _subShapeRef = -1)
    {
        pSetDecimal(p, val, _subShapeRef);
    }

    void pSetDecimal(properties_decimal p, decimal val, int _subShapeRef)
    {
        switch (p)
        {
            case properties_decimal.horLength:
                switch (_subShapeRef)
                {
                    case 0:
                        subShapeHorLength = val;
                        break;
                    case 1:
                        subShape2HorLength = val;
                        break;
                    case 2:
                        subShape3HorLength = val;
                        break;
                }
                break;
            case properties_decimal.verLength:
                switch (_subShapeRef)
                {
                    case 0:
                        subShapeVerLength = val;
                        break;
                    case 1:
                        subShape2VerLength = val;
                        break;
                    case 2:
                        subShape3VerLength = val;
                        break;
                }
                break;
            case properties_decimal.horOffset:
                switch (_subShapeRef)
                {
                    case 0:
                        subShapeHorOffset = val;
                        break;
                    case 1:
                        subShape2HorOffset = val;
                        break;
                    case 2:
                        subShape3HorOffset = val;
                        break;
                }
                break;
            case properties_decimal.verOffset:
                switch (_subShapeRef)
                {
                    case 0:
                        subShapeVerOffset = val;
                        break;
                    case 1:
                        subShape2VerOffset = val;
                        break;
                    case 2:
                        subShape3VerOffset = val;
                        break;
                }
                break;
            case properties_decimal.gHorOffset:
                globalHorOffset = val;
                break;
            case properties_decimal.gVerOffset:
                globalVerOffset = val;
                break;
            case properties_decimal.iCR:
                innerCRR = val;
                break;
            case properties_decimal.oCR:
                outerCRR = val;
                break;
            case properties_decimal.sBias:
                sideBias = val;
                break;
            case properties_decimal.hTBias:
                horTipBias = val;
                break;
            case properties_decimal.hTNVar:
                horTipBiasNVar = val;
                break;
            case properties_decimal.hTPVar:
                horTipBiasPVar = val;
                break;
            case properties_decimal.vTBias:
                verTipBias = val;
                break;
            case properties_decimal.vTNVar:
                verTipBiasNVar = val;
                break;
            case properties_decimal.vTPVar:
                verTipBiasPVar = val;
                break;
            case properties_decimal.lwr:
                LWR = val;
                break;
            case properties_decimal.lwrFreq:
                LWRNoiseFreq = val;
                break;
            case properties_decimal.lwr2:
                LWR2 = val;
                break;
            case properties_decimal.lwr2Freq:
                LWR2NoiseFreq = val;
                break;
            case properties_decimal.eTension:
                edgeSlideTension = val;
                break;
            case properties_decimal.rot:
                rotation = val;
                break;
            case properties_decimal.pBias:
                proximityBias = val;
                break;
            case properties_decimal.pBiasDist:
                proximityIsoDistance = val;
                break;
            case properties_decimal.proxSideRaysMultiplier:
                proxSideRaysMultiplier = val;
                break;
            case properties_decimal.rayExtension:
                rayExtension = val;
                break;
            case properties_decimal.keyhole_factor:
                gcRayExtension = val;
                break;
        }
    }

    public void defaultDecimal(properties_decimal p, int _subShapeRef = -1)
    {
        pDefaultDecimal(p, _subShapeRef);
    }

    void pDefaultDecimal(properties_decimal p, int _subShapeRef)
    {
        switch (p)
        {
            case properties_decimal.horLength:
                switch (_subShapeRef)
                {
                    case 0:
                        subShapeHorLength = default_subShapeHorLength;
                        break;
                    case 1:
                        subShape2HorLength = default_subShapeHorLength;
                        break;
                    case 2:
                        subShape3HorLength = default_subShapeHorLength;
                        break;
                }
                break;
            case properties_decimal.verLength:
                switch (_subShapeRef)
                {
                    case 0:
                        subShapeVerLength = default_subShapeVerLength;
                        break;
                    case 1:
                        subShape2VerLength = default_subShapeVerLength;
                        break;
                    case 2:
                        subShape3VerLength = default_subShapeVerLength;
                        break;
                }
                break;
            case properties_decimal.horOffset:
                switch (_subShapeRef)
                {
                    case 0:
                        subShapeHorOffset = default_subShapeHorOffset;
                        break;
                    case 1:
                        subShape2HorOffset = default_subShape2HorOffset;
                        break;
                    case 2:
                        subShape3HorOffset = default_subShape3HorOffset;
                        break;
                }
                break;
            case properties_decimal.verOffset:
                switch (_subShapeRef)
                {
                    case 0:
                        subShapeVerOffset = default_subShapeVerOffset;
                        break;
                    case 1:
                        subShape2VerOffset = default_subShape2VerOffset;
                        break;
                    case 2:
                        subShape3VerOffset = default_subShape3VerOffset;
                        break;
                }
                break;
            case properties_decimal.gHorOffset:
                globalHorOffset = default_globalHorOffset;
                break;
            case properties_decimal.gVerOffset:
                globalVerOffset = default_globalVerOffset;
                break;
            case properties_decimal.iCR:
                innerCRR = default_innerCRR;
                break;
            case properties_decimal.oCR:
                outerCRR = default_outerCRR;
                break;
            case properties_decimal.sBias:
                sideBias = default_sideBias;
                break;
            case properties_decimal.hTBias:
                horTipBias = default_horTipBias;
                break;
            case properties_decimal.hTNVar:
                horTipBiasNVar = default_horTipBiasNVar;
                break;
            case properties_decimal.hTPVar:
                horTipBiasPVar = default_horTipBiasPVar;
                break;
            case properties_decimal.vTBias:
                verTipBias = default_verTipBias;
                break;
            case properties_decimal.vTNVar:
                verTipBiasNVar = default_verTipBiasNVar;
                break;
            case properties_decimal.vTPVar:
                verTipBiasPVar = default_verTipBiasPVar;
                break;
            case properties_decimal.lwr:
                LWR = default_LWR;
                break;
            case properties_decimal.lwrFreq:
                LWRNoiseFreq = default_LWRNoiseFreq;
                break;
            case properties_decimal.lwr2:
                LWR2 = default_LWR;
                break;
            case properties_decimal.lwr2Freq:
                LWR2NoiseFreq = default_LWRNoiseFreq;
                break;
            case properties_decimal.eTension:
                edgeSlideTension = default_edgeSlideTension;
                break;
            case properties_decimal.rot:
                rotation = default_rotation;
                break;
            case properties_decimal.pBias:
                proximityBias = default_proximityBias;
                break;
            case properties_decimal.pBiasDist:
                proximityIsoDistance = default_proximityIsoDistance;
                break;
            case properties_decimal.proxSideRaysMultiplier:
                proxSideRaysMultiplier = default_proximitySideRaysFallOffMultiplier;
                break;
            case properties_decimal.rayExtension:
                rayExtension = default_rayExtension;
                break;
            case properties_decimal.keyhole_factor:
                gcRayExtension = default_rayExtension;
                break;
        }
    }

    public static decimal getDefaultDecimal(properties_decimal p, int _subShapeRef = -1)
    {
        return pGetDefaultDecimal(p, _subShapeRef);
    }

    static decimal pGetDefaultDecimal(properties_decimal p, int _subShapeRef)
    {
        decimal ret = 0;
        switch (p)
        {
            case properties_decimal.horLength:
                switch (_subShapeRef)
                {
                    case 0:
                        ret = default_subShapeHorLength;
                        break;
                    case 1:
                        ret = default_subShapeHorLength;
                        break;
                    case 2:
                        ret = default_subShapeHorLength;
                        break;
                }
                break;
            case properties_decimal.verLength:
                switch (_subShapeRef)
                {
                    case 0:
                        ret = default_subShapeVerLength;
                        break;
                    case 1:
                        ret = default_subShapeVerLength;
                        break;
                    case 2:
                        ret = default_subShapeVerLength;
                        break;
                }
                break;
            case properties_decimal.horOffset:
                switch (_subShapeRef)
                {
                    case 0:
                        ret = default_subShapeHorOffset;
                        break;
                    case 1:
                        ret = default_subShape2HorOffset;
                        break;
                    case 2:
                        ret = default_subShape3HorOffset;
                        break;
                }
                break;
            case properties_decimal.verOffset:
                switch (_subShapeRef)
                {
                    case 0:
                        ret = default_subShapeVerOffset;
                        break;
                    case 1:
                        ret = default_subShape2VerOffset;
                        break;
                    case 2:
                        ret = default_subShape3VerOffset;
                        break;
                }
                break;
            case properties_decimal.gHorOffset:
                ret = default_globalHorOffset;
                break;
            case properties_decimal.gVerOffset:
                ret = default_globalVerOffset;
                break;
            case properties_decimal.iCR:
                ret = default_innerCRR;
                break;
            case properties_decimal.oCR:
                ret = default_outerCRR;
                break;
            case properties_decimal.sBias:
                ret = default_sideBias;
                break;
            case properties_decimal.hTBias:
                ret = default_horTipBias;
                break;
            case properties_decimal.hTNVar:
                ret = default_horTipBiasNVar;
                break;
            case properties_decimal.hTPVar:
                ret = default_horTipBiasPVar;
                break;
            case properties_decimal.vTBias:
                ret = default_verTipBias;
                break;
            case properties_decimal.vTNVar:
                ret = default_verTipBiasNVar;
                break;
            case properties_decimal.vTPVar:
                ret = default_verTipBiasPVar;
                break;
            case properties_decimal.lwr:
                ret = default_LWR;
                break;
            case properties_decimal.lwrFreq:
                ret = default_LWRNoiseFreq;
                break;
            case properties_decimal.lwr2:
                ret = default_LWR;
                break;
            case properties_decimal.lwr2Freq:
                ret = default_LWRNoiseFreq;
                break;
            case properties_decimal.eTension:
                ret = default_edgeSlideTension;
                break;
            case properties_decimal.rot:
                ret = default_rotation;
                break;
            case properties_decimal.pBias:
                ret = default_proximityBias;
                break;
            case properties_decimal.pBiasDist:
                ret = default_proximityIsoDistance;
                break;
            case properties_decimal.proxSideRaysMultiplier:
                ret = default_proximitySideRaysFallOffMultiplier;
                break;
            case properties_decimal.rayExtension:
                ret = default_rayExtension;
                break;
            case properties_decimal.keyhole_factor:
                ret = default_rayExtension;
                break;
        }

        return ret;
    }

    private string layerName = "";

    public enum properties_s
    {
        s_name
    }

    public void setString(properties_s p, string val)
    {
        pSetString(p, val);
    }

    void pSetString(properties_s p, string val)
    {
        switch (p)
        {
            case properties_s.s_name:
                layerName = val;
                break;
            default:
                layerName = "";
                break;
        }
        
    }

    public string getString(properties_s p)
    {
        return pGetString(p);
    }

    string pGetString(properties_s p)
    {
        string ret = "";
        switch (p)
        {
            case properties_s.s_name:
                ret = layerName;
                break;
        }

        return ret;
    }

    public void defaultString(properties_s p)
    {
        pDefaultString(p);
    }

    void pDefaultString(properties_s p)
    {
        switch (p)
        {
            case properties_s.s_name:
                layerName = default_layerName;
                break;
        }
    }
    
    public static string getDefaultString(properties_s p)
    {
        return pGetDefaultString(p);
    }

    static string pGetDefaultString(properties_s p)
    {
        string ret = "";
        switch (p)
        {
            case properties_s.s_name:
                ret = default_layerName;
                break;
        }

        return ret;
    }
    
}
