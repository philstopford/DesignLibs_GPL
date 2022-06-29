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

    public enum subShapeLocations { TL, TR, BL, BR, TS, RS, BS, LS, C }
    
    public static List<string> getPolyFillTypes()
    {
        return pGetPolyFillTypes();
    }

    private static List<string> pGetPolyFillTypes()
    {
        return polyFillTypes;
    }

    private static readonly List<string> availableSubShapePositions = new List<string>
    { "Corner: Top Left", "Corner: Top Right", "Corner: Bottom Left", "Corner: Bottom Right",
        "Middle: Top Side", "Middle: Right Side", "Middle: Bottom Side", "Middle: Left Side",
        "Center"};

    private static readonly List<string> availableTipsLocations = new List<string>
    { "None", "Left", "Right", "Left & Right",
        "Top", "Bottom", "Top & Bottom",
        "Top & Left", "Top & Right", "Top & Left & Right",
        "Bottom & Left", "Bottom & Right", "Bottom & Left & Right",
        "Top & Bottom & Left", "Top & Bottom & Right",
        "All"};

    private static readonly List<string> polyFillTypes = new List<string> { "Even/Odd", "Non-zero", "Positive", "Negative" };

    public enum PolyFill { pftEvenOdd, pftNonZero, pftPositive, pftNegative }

    public enum typeShapes { none, rectangle, L, T, X, U, S, GEOCORE, BOOLEAN }

    private static readonly int default_subShapeTipLocIndex = 0;
    private static readonly int default_subShape2TipLocIndex = 0;

    private static readonly int default_subShape3TipLocIndex = 0;
    private static readonly int default_enabled = 0;
    private static readonly int default_geoCoreShapeEngine = 0;

    private static readonly int default_shapeIndex = (int)typeShapes.none;
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
    private int gCSEngine = default_geoCoreShapeEngine;
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
        gCSEngine,
        subShapeTipLocIndex, subShape2TipLocIndex, subShape3TipLocIndex,
        subShapeRefIndex,posInSubShapeIndex,
        edgeSlide,
        proxRays,proxSideRaysFallOff,
        lwrType, lwr2Type,
        flipH, flipV, alignX, alignY,
    }

    public int getInt(properties_i p)
    {
        return pGetInt(p);
    }

    int pGetInt(properties_i p)
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
            case properties_i.gCSEngine:
                ret = gCSEngine;
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

    public void setInt(properties_i p, int val)
    {
        pSetInt(p, val);
    }

    void pSetInt(properties_i p, int val)
    {
        switch (p)
        {
            case properties_i.enabled:
                enabled = val;
                break;
            case properties_i.shapeIndex:
                shapeIndex = val;
                break;
            case properties_i.gCSEngine:
                gCSEngine = val;
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
            case properties_i.gCSEngine:
                gCSEngine = default_geoCoreShapeEngine;
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
            case properties_i.gCSEngine:
                val = default_geoCoreShapeEngine;
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

    private static readonly decimal default_subShapeHorLength = 0;
    private static readonly decimal default_subShapeHorOffset = 0;
    private static readonly decimal default_subShapeVerLength = 0;
    private static readonly decimal default_subShapeVerOffset = 0;
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
    private static readonly decimal default_rotation = 0;
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

    private decimal subShapeHorLength = default_subShapeHorLength;
    private decimal subShapeHorOffset = default_subShapeHorOffset;
    private decimal subShapeVerLength = default_subShapeVerLength;
    private decimal subShapeVerOffset = default_subShapeVerOffset;
    private decimal subShape2HorLength = default_subShape2HorLength;
    private decimal subShape2HorOffset = default_subShape2HorOffset;
    private decimal subShape2VerLength = default_subShape2VerLength;
    private decimal subShape2VerOffset = default_subShape2VerOffset;
    private decimal subShape3HorLength = default_subShape3HorLength;
    private decimal subShape3HorOffset = default_subShape3HorOffset;
    private decimal subShape3VerLength = default_subShape3VerLength;
    private decimal subShape3VerOffset = default_subShape3HorOffset;
    private decimal innerCRR = default_innerCRR;
    private decimal outerCRR = default_outerCRR;
    private decimal globalHorOffset = default_globalHorOffset;
    private decimal globalVerOffset = default_globalVerOffset;
    private decimal sideBias = default_sideBias;
    private decimal horTipBias = default_horTipBias;
    private decimal horTipBiasPVar = default_horTipBias;
    private decimal horTipBiasNVar;
    private decimal verTipBias;
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
        s0HorLength, s1HorLength, s2HorLength,
        s0VerLength, s1VerLength, s2VerLength,
        s0HorOffset, s1HorOffset, s2HorOffset,
        s0VerOffset, s1VerOffset, s2VerOffset,
        gHorOffset, gVerOffset,
        iCR, oCR,
        sBias, hTBias, hTNVar, hTPVar, vTBias, vTNVar, vTPVar,
        lwr, lwrFreq, lwr2, lwr2Freq,
        eTension, rot,
        pBias, pBiasDist, proxSideRaysMultiplier,
        rayExtension, keyhole_factor
    }

    public decimal getDecimal(properties_decimal p)
    {
        return pGetDecimal(p);
    }

    decimal pGetDecimal(properties_decimal p)
    {
        decimal ret = 0;
        switch (p)
        {
            case properties_decimal.s0HorLength:
                ret = subShapeHorLength;
                break;
            case properties_decimal.s1HorLength:
                ret = subShape2HorLength;
                break;
            case properties_decimal.s2HorLength:
                ret = subShape3HorLength;
                break;
            case properties_decimal.s0VerLength:
                ret = subShapeVerLength;
                break;
            case properties_decimal.s1VerLength:
                ret = subShape2VerLength;
                break;
            case properties_decimal.s2VerLength:
                ret = subShape3VerLength;
                break;
            case properties_decimal.s0HorOffset:
                ret = subShapeHorOffset;
                break;
            case properties_decimal.s1HorOffset:
                ret = subShape2HorOffset;
                break;
            case properties_decimal.s2HorOffset:
                ret = subShape3HorOffset;
                break;
            case properties_decimal.s0VerOffset:
                ret = subShapeVerOffset;
                break;
            case properties_decimal.s1VerOffset:
                ret = subShape2VerOffset;
                break;
            case properties_decimal.s2VerOffset:
                ret = subShape3VerOffset;
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

    public void setDecimal(properties_decimal p, decimal val)
    {
        pSetDecimal(p, val);
    }

    void pSetDecimal(properties_decimal p, decimal val)
    {
        switch (p)
        {
            case properties_decimal.s0HorLength:
                subShapeHorLength = val;
                break;
            case properties_decimal.s1HorLength:
                subShape2HorLength = val;
                break;
            case properties_decimal.s2HorLength:
                subShape3HorLength = val;
                break;
            case properties_decimal.s0VerLength:
                subShapeVerLength = val;
                break;
            case properties_decimal.s1VerLength:
                subShape2VerLength = val;
                break;
            case properties_decimal.s2VerLength:
                subShape3VerLength = val;
                break;
            case properties_decimal.s0HorOffset:
                subShapeHorOffset = val;
                break;
            case properties_decimal.s1HorOffset:
                subShape2HorOffset = val;
                break;
            case properties_decimal.s2HorOffset:
                subShape3HorOffset = val;
                break;
            case properties_decimal.s0VerOffset:
                subShapeVerOffset = val;
                break;
            case properties_decimal.s1VerOffset:
                subShape2VerOffset = val;
                break;
            case properties_decimal.s2VerOffset:
                subShape3VerOffset = val;
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

    public void defaultDecimal(properties_decimal p)
    {
        pDefaultDecimal(p);
    }

    void pDefaultDecimal(properties_decimal p)
    {
        switch (p)
        {
            case properties_decimal.s0HorLength:
                subShapeHorLength = default_subShapeHorLength;
                break;
            case properties_decimal.s1HorLength:
                subShape2HorLength = default_subShape2HorLength;
                break;
            case properties_decimal.s2HorLength:
                subShape3HorLength = default_subShape3HorLength;
                break;
            case properties_decimal.s0VerLength:
                subShapeVerLength = default_subShapeVerLength;
                break;
            case properties_decimal.s1VerLength:
                subShape2VerLength = default_subShape2VerLength;
                break;
            case properties_decimal.s2VerLength:
                subShape3VerLength = default_subShape3VerLength;
                break;
            case properties_decimal.s0HorOffset:
                subShapeHorOffset = default_subShapeHorOffset;
                break;
            case properties_decimal.s1HorOffset:
                subShape2HorOffset = default_subShape2HorOffset;
                break;
            case properties_decimal.s2HorOffset:
                subShape3HorOffset = default_subShape3HorOffset;
                break;
            case properties_decimal.s0VerOffset:
                subShapeVerOffset = default_subShapeVerOffset;
                break;
            case properties_decimal.s1VerOffset:
                subShape2VerOffset = default_subShape2VerOffset;
                break;
            case properties_decimal.s2VerOffset:
                subShape3VerOffset = default_subShape3VerOffset;
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

    public static decimal getDefaultDecimal(properties_decimal p)
    {
        return pGetDefaultDecimal(p);
    }

    static decimal pGetDefaultDecimal(properties_decimal p)
    {
        decimal ret = 0;
        switch (p)
        {
            case properties_decimal.s0HorLength:
                ret = default_subShapeHorLength;
                break;
            case properties_decimal.s1HorLength:
                ret = default_subShape2HorLength;
                break;
            case properties_decimal.s2HorLength:
                ret = default_subShape3HorLength;
                break;
            case properties_decimal.s0VerLength:
                ret = default_subShapeVerLength;
                break;
            case properties_decimal.s1VerLength:
                ret = default_subShape2VerLength;
                break;
            case properties_decimal.s2VerLength:
                ret = default_subShape3VerLength;
                break;
            case properties_decimal.s0HorOffset:
                ret = default_subShapeHorOffset;
                break;
            case properties_decimal.s1HorOffset:
                ret = default_subShape2HorOffset;
                break;
            case properties_decimal.s2HorOffset:
                ret = default_subShape3HorOffset;
                break;
            case properties_decimal.s0VerOffset:
                ret = default_subShapeVerOffset;
                break;
            case properties_decimal.s1VerOffset:
                ret = default_subShape2VerOffset;
                break;
            case properties_decimal.s2VerOffset:
                ret = default_subShape3VerOffset;
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
