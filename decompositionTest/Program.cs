using System.Diagnostics;
using Clipper2Lib;
using geoWrangler;
using geoCoreLib;
using PartitionTestGeometrySource;
using NUnit.Framework;
using utility;

namespace partitionTest;

internal class Program
{
    private static void Main(string[] args)
    {
        problemcase();
        
        Console.WriteLine("Part One");
        partOne();

        Console.WriteLine("Part Two");
        partTwo();

        Console.WriteLine("Part Three (takes a while)");
        partThree();

        Console.WriteLine("Part Four (takes less time than part three)");
        partFour();
        
        Console.WriteLine("Part Five");
        partFive();
    }

    private static void problemcase()
    {
        PathD source = Clipper.MakePath(new double[]
        {
            60, 163,
            60, 173,
            70, 173,
            70, 213,
            60, 213,
            60, 243,
            70, 243,
            70, 304,
            80, 304,
            80, 293,
            90, 293,
            90, 303,
            100, 303,
            100, 273,
            90, 273,
            90, 283,
            80, 283,
            80, 243,
            90, 243,
            90, 233,
            160, 233,
            160, 243,
            150, 243,
            150, 253,
            180, 253,
            180, 243,
            170, 243,
            170, 233,
            210, 233,
            210, 243,
            220, 243,
            220, 315,
            210, 315,
            210, 305,
            200, 305,
            200, 335,
            210, 335,
            210, 325,
            220, 325,
            220, 365,
            210, 365,
            210, 375,
            150, 375,
            150, 385,
            160, 385,
            160, 395,
            150, 395,
            150, 405,
            180, 405,
            180, 395,
            170, 395,
            170, 385,
            210, 385,
            210, 395,
            220, 395,
            220, 466,
            210, 466,
            210, 456,
            200, 456,
            200, 486,
            210, 486,
            210, 476,
            220, 476,
            220, 516,
            210, 516,
            210, 526,
            150, 526,
            150, 536,
            160, 536,
            160, 546,
            150, 546,
            150, 556,
            180, 556,
            180, 546,
            170, 546,
            170, 536,
            210, 536,
            210, 546,
            220, 546,
            220, 618,
            210, 618,
            210, 608,
            200, 608,
            200, 638,
            210, 638,
            210, 628,
            220, 628,
            220, 668,
            210, 668,
            210, 678,
            140, 678,
            140, 668,
            150, 668,
            150, 658,
            120, 658,
            120, 668,
            130, 668,
            130, 678,
            90, 678,
            90, 668,
            80, 668,
            80, 607,
            70, 607,
            70, 618,
            60, 618,
            60, 628,
            70, 628,
            70, 668,
            60, 668,
            60, 698,
            70, 698,
            70, 759,
            80, 759,
            80, 748,
            90, 748,
            90, 758,
            100, 758,
            100, 728,
            90, 728,
            90, 738,
            80, 738,
            80, 698,
            90, 698,
            90, 688,
            160, 688,
            160, 698,
            150, 698,
            150, 708,
            180, 708,
            180, 698,
            170, 698,
            170, 688,
            210, 688,
            210, 698,
            220, 698,
            220, 770,
            210, 770,
            210, 760,
            200, 760,
            200, 790,
            210, 790,
            210, 780,
            220, 780,
            220, 820,
            210, 820,
            210, 830,
            150, 830,
            150, 840,
            160, 840,
            160, 850,
            150, 850,
            150, 860,
            180, 860,
            180, 850,
            170, 850,
            170, 840,
            210, 840,
            210, 850,
            220, 850,
            220, 921,
            210, 921,
            210, 911,
            200, 911,
            200, 941,
            210, 941,
            210, 931,
            220, 931,
            220, 971,
            210, 971,
            210, 981,
            150, 981,
            150, 991,
            160, 991,
            160, 1001,
            150, 1001,
            150, 1011,
            180, 1011,
            180, 1001,
            170, 1001,
            170, 991,
            210, 991,
            210, 1001,
            220, 1001,
            220, 1073,
            210, 1073,
            210, 1063,
            200, 1063,
            200, 1093,
            210, 1093,
            210, 1083,
            220, 1083,
            220, 1123,
            210, 1123,
            210, 1133,
            140, 1133,
            140, 1123,
            150, 1123,
            150, 1113,
            120, 1113,
            120, 1123,
            130, 1123,
            130, 1133,
            90, 1133,
            90, 1123,
            80, 1123,
            80, 1062,
            70, 1062,
            70, 1214,
            80, 1214,
            80, 1203,
            90, 1203,
            90, 1213,
            100, 1213,
            100, 1183,
            90, 1183,
            90, 1193,
            80, 1193,
            80, 1153,
            90, 1153,
            90, 1143,
            160, 1143,
            160, 1153,
            150, 1153,
            150, 1163,
            180, 1163,
            180, 1153,
            170, 1153,
            170, 1143,
            210, 1143,
            210, 1153,
            220, 1153,
            220, 1225,
            210, 1225,
            210, 1215,
            200, 1215,
            200, 1245,
            210, 1245,
            210, 1235,
            220, 1235,
            220, 1275,
            210, 1275,
            210, 1285,
            150, 1285,
            150, 1295,
            160, 1295,
            160, 1305,
            150, 1305,
            150, 1315,
            180, 1315,
            180, 1305,
            170, 1305,
            170, 1295,
            210, 1295,
            210, 1305,
            220, 1305,
            220, 1365,
            230, 1365,
            230, 1355,
            240, 1355,
            240, 1365,
            250, 1365,
            250, 1335,
            240, 1335,
            240, 1345,
            230, 1345,
            230, 1305,
            240, 1305,
            240, 1295,
            300, 1295,
            300, 1285,
            290, 1285,
            290, 1275,
            300, 1275,
            300, 1265,
            270, 1265,
            270, 1275,
            280, 1275,
            280, 1285,
            240, 1285,
            240, 1275,
            230, 1275,
            230, 1203,
            240, 1203,
            240, 1213,
            250, 1213,
            250, 1183,
            240, 1183,
            240, 1193,
            230, 1193,
            230, 1153,
            240, 1153,
            240, 1143,
            310, 1143,
            310, 1153,
            300, 1153,
            300, 1163,
            330, 1163,
            330, 1153,
            320, 1153,
            320, 1143,
            360, 1143,
            360, 1153,
            370, 1153,
            370, 1214,
            380, 1214,
            380, 1203,
            390, 1203,
            390, 1213,
            400, 1213,
            400, 1183,
            390, 1183,
            390, 1193,
            380, 1193,
            380, 1153,
            390, 1153,
            390, 1143,
            450, 1143,
            450, 1133,
            440, 1133,
            440, 1123,
            450, 1123,
            450, 1113,
            420, 1113,
            420, 1123,
            430, 1123,
            430, 1133,
            390, 1133,
            390, 1123,
            380, 1123,
            380, 1062,
            370, 1062,
            370, 1073,
            360, 1073,
            360, 1063,
            350, 1063,
            350, 1093,
            360, 1093,
            360, 1083,
            370, 1083,
            370, 1123,
            360, 1123,
            360, 1133,
            290, 1133,
            290, 1123,
            300, 1123,
            300, 1113,
            270, 1113,
            270, 1123,
            280, 1123,
            280, 1133,
            240, 1133,
            240, 1123,
            230, 1123,
            230, 1051,
            240, 1051,
            240, 1061,
            250, 1061,
            250, 1031,
            240, 1031,
            240, 1041,
            230, 1041,
            230, 1001,
            240, 1001,
            240, 991,
            300, 991,
            300, 981,
            290, 981,
            290, 971,
            300, 971,
            300, 961,
            270, 961,
            270, 971,
            280, 971,
            280, 981,
            240, 981,
            240, 971,
            230, 971,
            230, 900,
            240, 900,
            240, 910,
            250, 910,
            250, 880,
            240, 880,
            240, 890,
            230, 890,
            230, 850,
            240, 850,
            240, 840,
            300, 840,
            300, 830,
            290, 830,
            290, 820,
            300, 820,
            300, 810,
            270, 810,
            270, 820,
            280, 820,
            280, 830,
            240, 830,
            240, 820,
            230, 820,
            230, 748,
            240, 748,
            240, 758,
            250, 758,
            250, 728,
            240, 728,
            240, 738,
            230, 738,
            230, 698,
            240, 698,
            240, 688,
            310, 688,
            310, 698,
            300, 698,
            300, 708,
            330, 708,
            330, 698,
            320, 698,
            320, 688,
            360, 688,
            360, 698,
            370, 698,
            370, 759,
            380, 759,
            380, 748,
            390, 748,
            390, 758,
            400, 758,
            400, 728,
            390, 728,
            390, 738,
            380, 738,
            380, 698,
            390, 698,
            390, 688,
            450, 688,
            450, 678,
            440, 678,
            440, 668,
            450, 668,
            450, 658,
            420, 658,
            420, 668,
            430, 668,
            430, 678,
            390, 678,
            390, 668,
            380, 668,
            380, 607,
            370, 607,
            370, 618,
            360, 618,
            360, 608,
            350, 608,
            350, 638,
            360, 638,
            360, 628,
            370, 628,
            370, 668,
            360, 668,
            360, 678,
            290, 678,
            290, 668,
            300, 668,
            300, 658,
            270, 658,
            270, 668,
            280, 668,
            280, 678,
            240, 678,
            240, 668,
            230, 668,
            230, 596,
            240, 596,
            240, 606,
            250, 606,
            250, 576,
            240, 576,
            240, 586,
            230, 586,
            230, 546,
            240, 546,
            240, 536,
            300, 536,
            300, 526,
            290, 526,
            290, 516,
            300, 516,
            300, 506,
            270, 506,
            270, 516,
            280, 516,
            280, 526,
            240, 526,
            240, 516,
            230, 516,
            230, 445,
            240, 445,
            240, 455,
            250, 455,
            250, 425,
            240, 425,
            240, 435,
            230, 435,
            230, 395,
            240, 395,
            240, 385,
            300, 385,
            300, 375,
            290, 375,
            290, 365,
            300, 365,
            300, 355,
            270, 355,
            270, 365,
            280, 365,
            280, 375,
            240, 375,
            240, 365,
            230, 365,
            230, 293,
            240, 293,
            240, 303,
            250, 303,
            250, 273,
            240, 273,
            240, 283,
            230, 283,
            230, 243,
            240, 243,
            240, 233,
            310, 233,
            310, 243,
            300, 243,
            300, 253,
            330, 253,
            330, 243,
            320, 243,
            320, 233,
            360, 233,
            360, 243,
            370, 243,
            370, 304,
            380, 304,
            380, 293,
            390, 293,
            390, 303,
            400, 303,
            400, 273,
            390, 273,
            390, 283,
            380, 283,
            380, 243,
            390, 243,
            390, 233,
            460, 233,
            460, 243,
            450, 243,
            450, 253,
            480, 253,
            480, 243,
            470, 243,
            470, 233,
            510, 233,
            510, 243,
            520, 243,
            520, 304,
            530, 304,
            530, 293,
            540, 293,
            540, 303,
            550, 303,
            550, 273,
            540, 273,
            540, 283,
            530, 283,
            530, 243,
            540, 243,
            540, 233,
            610, 233,
            610, 243,
            600, 243,
            600, 253,
            630, 253,
            630, 243,
            620, 243,
            620, 233,
            660, 233,
            660, 243,
            670, 243,
            670, 315,
            660, 315,
            660, 305,
            650, 305,
            650, 335,
            660, 335,
            660, 325,
            670, 325,
            670, 365,
            660, 365,
            660, 375,
            600, 375,
            600, 385,
            610, 385,
            610, 395,
            600, 395,
            600, 405,
            630, 405,
            630, 395,
            620, 395,
            620, 385,
            660, 385,
            660, 395,
            670, 395,
            670, 455,
            680, 455,
            680, 445,
            690, 445,
            690, 455,
            700, 455,
            700, 425,
            690, 425,
            690, 435,
            680, 435,
            680, 395,
            690, 395,
            690, 385,
            750, 385,
            750, 375,
            740, 375,
            740, 365,
            750, 365,
            750, 355,
            720, 355,
            720, 365,
            730, 365,
            730, 375,
            690, 375,
            690, 365,
            680, 365,
            680, 293,
            690, 293,
            690, 303,
            700, 303,
            700, 273,
            690, 273,
            690, 283,
            680, 283,
            680, 243,
            690, 243,
            690, 233,
            760, 233,
            760, 243,
            750, 243,
            750, 253,
            780, 253,
            780, 243,
            770, 243,
            770, 233,
            810, 233,
            810, 243,
            820, 243,
            820, 304,
            830, 304,
            830, 293,
            840, 293,
            840, 303,
            850, 303,
            850, 273,
            840, 273,
            840, 283,
            830, 283,
            830, 243,
            840, 243,
            840, 233,
            910, 233,
            910, 243,
            900, 243,
            900, 253,
            930, 253,
            930, 243,
            920, 243,
            920, 233,
            960, 233,
            960, 243,
            970, 243,
            970, 304,
            980, 304,
            980, 293,
            990, 293,
            990, 303,
            1000, 303,
            1000, 273,
            990, 273,
            990, 283,
            980, 283,
            980, 243,
            990, 243,
            990, 233,
            1060, 233,
            1060, 243,
            1050, 243,
            1050, 253,
            1080, 253,
            1080, 243,
            1070, 243,
            1070, 233,
            1110, 233,
            1110, 243,
            1120, 243,
            1120, 315,
            1110, 315,
            1110, 305,
            1100, 305,
            1100, 335,
            1110, 335,
            1110, 325,
            1120, 325,
            1120, 365,
            1110, 365,
            1110, 375,
            1050, 375,
            1050, 385,
            1060, 385,
            1060, 395,
            1050, 395,
            1050, 405,
            1080, 405,
            1080, 395,
            1070, 395,
            1070, 385,
            1110, 385,
            1110, 395,
            1120, 395,
            1120, 455,
            1130, 455,
            1130, 445,
            1140, 445,
            1140, 455,
            1150, 455,
            1150, 425,
            1140, 425,
            1140, 435,
            1130, 435,
            1130, 395,
            1140, 395,
            1140, 385,
            1200, 385,
            1200, 375,
            1190, 375,
            1190, 365,
            1200, 365,
            1200, 355,
            1170, 355,
            1170, 365,
            1180, 365,
            1180, 375,
            1140, 375,
            1140, 365,
            1130, 365,
            1130, 293,
            1140, 293,
            1140, 303,
            1150, 303,
            1150, 273,
            1140, 273,
            1140, 283,
            1130, 283,
            1130, 243,
            1140, 243,
            1140, 233,
            1210, 233,
            1210, 243,
            1200, 243,
            1200, 253,
            1230, 253,
            1230, 243,
            1220, 243,
            1220, 233,
            1260, 233,
            1260, 243,
            1270, 243,
            1270, 304,
            1280, 304,
            1280, 293,
            1290, 293,
            1290, 303,
            1300, 303,
            1300, 273,
            1290, 273,
            1290, 283,
            1280, 283,
            1280, 243,
            1290, 243,
            1290, 233,
            1350, 233,
            1350, 223,
            1340, 223,
            1340, 213,
            1350, 213,
            1350, 203,
            1320, 203,
            1320, 213,
            1330, 213,
            1330, 223,
            1290, 223,
            1290, 213,
            1280, 213,
            1280, 152,
            1270, 152,
            1270, 163,
            1260, 163,
            1260, 153,
            1250, 153,
            1250, 183,
            1260, 183,
            1260, 173,
            1270, 173,
            1270, 213,
            1260, 213,
            1260, 223,
            1190, 223,
            1190, 213,
            1200, 213,
            1200, 203,
            1170, 203,
            1170, 213,
            1180, 213,
            1180, 223,
            1140, 223,
            1140, 213,
            1130, 213,
            1130, 141,
            1140, 141,
            1140, 151,
            1150, 151,
            1150, 121,
            1140, 121,
            1140, 131,
            1130, 131,
            1130, 91,
            1140, 91,
            1140, 81,
            1200, 81,
            1200, 71,
            1190, 71,
            1190, 61,
            1200, 61,
            1200, 51,
            1170, 51,
            1170, 61,
            1180, 61,
            1180, 71,
            1140, 71,
            1140, 61,
            1130, 61,
            1130, 0,
            1120, 0,
            1120, 11,
            1110, 11,
            1110, 1,
            1100, 1,
            1100, 31,
            1110, 31,
            1110, 21,
            1120, 21,
            1120, 61,
            1110, 61,
            1110, 71,
            1050, 71,
            1050, 81,
            1060, 81,
            1060, 91,
            1050, 91,
            1050, 101,
            1080, 101,
            1080, 91,
            1070, 91,
            1070, 81,
            1110, 81,
            1110, 91,
            1120, 91,
            1120, 163,
            1110, 163,
            1110, 153,
            1100, 153,
            1100, 183,
            1110, 183,
            1110, 173,
            1120, 173,
            1120, 213,
            1110, 213,
            1110, 223,
            1040, 223,
            1040, 213,
            1050, 213,
            1050, 203,
            1020, 203,
            1020, 213,
            1030, 213,
            1030, 223,
            990, 223,
            990, 213,
            980, 213,
            980, 152,
            970, 152,
            970, 163,
            960, 163,
            960, 153,
            950, 153,
            950, 183,
            960, 183,
            960, 173,
            970, 173,
            970, 213,
            960, 213,
            960, 223,
            890, 223,
            890, 213,
            900, 213,
            900, 203,
            870, 203,
            870, 213,
            880, 213,
            880, 223,
            840, 223,
            840, 213,
            830, 213,
            830, 152,
            820, 152,
            820, 163,
            810, 163,
            810, 153,
            800, 153,
            800, 183,
            810, 183,
            810, 173,
            820, 173,
            820, 213,
            810, 213,
            810, 223,
            740, 223,
            740, 213,
            750, 213,
            750, 203,
            720, 203,
            720, 213,
            730, 213,
            730, 223,
            690, 223,
            690, 213,
            680, 213,
            680, 141,
            690, 141,
            690, 151,
            700, 151,
            700, 121,
            690, 121,
            690, 131,
            680, 131,
            680, 91,
            690, 91,
            690, 81,
            750, 81,
            750, 71,
            740, 71,
            740, 61,
            750, 61,
            750, 51,
            720, 51,
            720, 61,
            730, 61,
            730, 71,
            690, 71,
            690, 61,
            680, 61,
            680, 0,
            670, 0,
            670, 11,
            660, 11,
            660, 1,
            650, 1,
            650, 31,
            660, 31,
            660, 21,
            670, 21,
            670, 61,
            660, 61,
            660, 71,
            600, 71,
            600, 81,
            610, 81,
            610, 91,
            600, 91,
            600, 101,
            630, 101,
            630, 91,
            620, 91,
            620, 81,
            660, 81,
            660, 91,
            670, 91,
            670, 163,
            660, 163,
            660, 153,
            650, 153,
            650, 183,
            660, 183,
            660, 173,
            670, 173,
            670, 213,
            660, 213,
            660, 223,
            590, 223,
            590, 213,
            600, 213,
            600, 203,
            570, 203,
            570, 213,
            580, 213,
            580, 223,
            540, 223,
            540, 213,
            530, 213,
            530, 152,
            520, 152,
            520, 163,
            510, 163,
            510, 153,
            500, 153,
            500, 183,
            510, 183,
            510, 173,
            520, 173,
            520, 213,
            510, 213,
            510, 223,
            440, 223,
            440, 213,
            450, 213,
            450, 203,
            420, 203,
            420, 213,
            430, 213,
            430, 223,
            390, 223,
            390, 213,
            380, 213,
            380, 152,
            370, 152,
            370, 163,
            360, 163,
            360, 153,
            350, 153,
            350, 183,
            360, 183,
            360, 173,
            370, 173,
            370, 213,
            360, 213,
            360, 223,
            290, 223,
            290, 213,
            300, 213,
            300, 203,
            270, 203,
            270, 213,
            280, 213,
            280, 223,
            240, 223,
            240, 213,
            230, 213,
            230, 141,
            240, 141,
            240, 151,
            250, 151,
            250, 121,
            240, 121,
            240, 131,
            230, 131,
            230, 91,
            240, 91,
            240, 81,
            300, 81,
            300, 71,
            290, 71,
            290, 61,
            300, 61,
            300, 51,
            270, 51,
            270, 61,
            280, 61,
            280, 71,
            240, 71,
            240, 61,
            230, 61,
            230, 0,
            220, 0,
            220, 11,
            210, 11,
            210, 1,
            200, 1,
            200, 31,
            210, 31,
            210, 21,
            220, 21,
            220, 61,
            210, 61,
            210, 71,
            150, 71,
            150, 81,
            160, 81,
            160, 91,
            150, 91,
            150, 101,
            180, 101,
            180, 91,
            170, 91,
            170, 81,
            210, 81,
            210, 91,
            220, 91,
            220, 163,
            210, 163,
            210, 153,
            200, 153,
            200, 183,
            210, 183,
            210, 173,
            220, 173,
            220, 213,
            210, 213,
            210, 223,
            140, 223,
            140, 213,
            150, 213,
            150, 203,
            120, 203,
            120, 213,
            130, 213,
            130, 223,
            90, 223,
            90, 213,
            80, 213,
            80, 152,
            70, 152,
            70, 163,
            60, 163,
        });

        bool abort = false;
        int rayLength = 1000;

        PathsD test = GeoWrangler.rectangular_decomposition(ref abort, source, maxRayLength: rayLength);
        writeToLayout("debug", source, test);

        int x = 2;
    }
    
    private static void partOne()
    {
        // L
        PathD L = TestGeometry.getL();
        PathD rL = TestGeometry.getRL();

        // Reversed orientation.
        PathD L_ccw = new (L);
        L_ccw.Reverse();

        // U
        PathD U = TestGeometry.getU();

        // T
        PathD T = TestGeometry.getT();

        // X
        PathD X = TestGeometry.getX();

        // S
        PathD S = TestGeometry.getS();

        // Negative S
        PathD nS = TestGeometry.getnegS();

        // Complex 1
        PathD C1 = TestGeometry.getComplex1();


        // Complex 2
        PathD C2 = TestGeometry.getComplex2();

        // Complex 3
        PathD C3 = TestGeometry.getComplex3();

        // C3 = GeoWrangler.clockwiseAndReorder(C3);

        bool orth = GeoWrangler.orthogonal(C2, angularTolerance: 0);
        bool orth2 = GeoWrangler.orthogonal(C3, angularTolerance: 0);

        // Complex 10, rot 15
        PathD C10R15 = TestGeometry.getComplex10rot15();

        // Staircase
        PathD S1 = TestGeometry.getStaircase();

        // Rectangular decomposition will return non-orthogonal polygons in the output when encountered.

        int rayLength = 1000;

        bool abort = false;

        PathsD l = GeoWrangler.rectangular_decomposition(ref abort, L, maxRayLength: rayLength);
        writeToLayout("l", L, l);
        // Sign change expected.
        Assert.AreEqual(Clipper.Area(L), -Clipper.Area(l));
        Assert.AreEqual(2, l.Count);
        Assert.AreEqual(0, l[0][0].x);
        Assert.AreEqual(0, l[0][0].y);
        Assert.AreEqual(10, l[0][1].x);
        Assert.AreEqual(0, l[0][1].y);
        Assert.AreEqual(10, l[0][2].x);
        Assert.AreEqual(50, l[0][2].y);
        Assert.AreEqual(0, l[0][3].x);
        Assert.AreEqual(50, l[0][3].y);
        Assert.AreEqual(0, l[0][4].x);
        Assert.AreEqual(0, l[0][4].y);
        Assert.AreEqual(10, l[1][0].x);
        Assert.AreEqual(0, l[1][0].y);
        Assert.AreEqual(60, l[1][1].x);
        Assert.AreEqual(0, l[1][1].y);
        Assert.AreEqual(60, l[1][2].x);
        Assert.AreEqual(20, l[1][2].y);
        Assert.AreEqual(10, l[1][3].x);
        Assert.AreEqual(20, l[1][3].y);
        Assert.AreEqual(10, l[1][4].x);
        Assert.AreEqual(0, l[1][4].y);
        
        PathsD lccw = GeoWrangler.rectangular_decomposition(ref abort, L_ccw, maxRayLength: rayLength);
        writeToLayout("lccw", L_ccw, lccw);
        Assert.AreEqual(Clipper.Area(L_ccw), Clipper.Area(lccw));
        Assert.AreEqual(2, lccw.Count);
        Assert.AreEqual(0, lccw[0][0].x);
        Assert.AreEqual(0, lccw[0][0].y);
        Assert.AreEqual(10, lccw[0][1].x);
        Assert.AreEqual(0, lccw[0][1].y);
        Assert.AreEqual(10, lccw[0][2].x);
        Assert.AreEqual(50, lccw[0][2].y);
        Assert.AreEqual(0, lccw[0][3].x);
        Assert.AreEqual(50, lccw[0][3].y);
        Assert.AreEqual(0, lccw[0][4].x);
        Assert.AreEqual(0, lccw[0][4].y);
        Assert.AreEqual(10, lccw[1][0].x);
        Assert.AreEqual(0, lccw[1][0].y);
        Assert.AreEqual(60, lccw[1][1].x);
        Assert.AreEqual(0, lccw[1][1].y);
        Assert.AreEqual(60, lccw[1][2].x);
        Assert.AreEqual(20, lccw[1][2].y);
        Assert.AreEqual(10, lccw[1][3].x);
        Assert.AreEqual(20, lccw[1][3].y);
        Assert.AreEqual(10, lccw[1][4].x);
        Assert.AreEqual(0, lccw[1][4].y);
        
        PathsD rl = GeoWrangler.rectangular_decomposition(ref abort, rL, maxRayLength: rayLength);
        writeToLayout("rl", rL, rl);
        // Sign change expected.
        Assert.AreEqual(Clipper.Area(rL), -Clipper.Area(rl));
        Assert.AreEqual(2, rl.Count);
        Assert.AreEqual(10, rl[0][0].x);
        Assert.AreEqual(0, rl[0][0].y);
        Assert.AreEqual(60, rl[0][1].x);
        Assert.AreEqual(0, rl[0][1].y);
        Assert.AreEqual(60, rl[0][2].x);
        Assert.AreEqual(50, rl[0][2].y);
        Assert.AreEqual(10, rl[0][3].x);
        Assert.AreEqual(50, rl[0][3].y);
        Assert.AreEqual(10, rl[0][4].x);
        Assert.AreEqual(0, rl[0][4].y);
        Assert.AreEqual(0, rl[1][0].x);
        Assert.AreEqual(0, rl[1][0].y);
        Assert.AreEqual(10, rl[1][1].x);
        Assert.AreEqual(0, rl[1][1].y);
        Assert.AreEqual(10, rl[1][2].x);
        Assert.AreEqual(20, rl[1][2].y);
        Assert.AreEqual(0, rl[1][3].x);
        Assert.AreEqual(20, rl[1][3].y);
        Assert.AreEqual(0, rl[1][4].x);
        Assert.AreEqual(0, rl[1][4].y);

        PathsD u = GeoWrangler.rectangular_decomposition(ref abort, U, maxRayLength: rayLength);
        writeToLayout("u", U, u);
        // Sign change expected.
        Assert.AreEqual(Clipper.Area(U), -Clipper.Area(u));
        Assert.AreEqual(3, u.Count);
        Assert.AreEqual(0, u[0][0].x);
        Assert.AreEqual(0, u[0][0].y);
        Assert.AreEqual(10, u[0][1].x);
        Assert.AreEqual(0, u[0][1].y);
        Assert.AreEqual(10, u[0][2].x);
        Assert.AreEqual(50, u[0][2].y);
        Assert.AreEqual(0, u[0][3].x);
        Assert.AreEqual(50, u[0][3].y);
        Assert.AreEqual(0, u[0][4].x);
        Assert.AreEqual(0, u[0][4].y);
        Assert.AreEqual(60, u[1][0].x);
        Assert.AreEqual(0, u[1][0].y);
        Assert.AreEqual(120, u[1][1].x);
        Assert.AreEqual(0, u[1][1].y);
        Assert.AreEqual(120, u[1][2].x);
        Assert.AreEqual(80, u[1][2].y);
        Assert.AreEqual(60, u[1][3].x);
        Assert.AreEqual(80, u[1][3].y);
        Assert.AreEqual(60, u[1][4].x);
        Assert.AreEqual(0, u[1][4].y);
        Assert.AreEqual(10, u[2][0].x);
        Assert.AreEqual(0, u[2][0].y);
        Assert.AreEqual(60, u[2][1].x);
        Assert.AreEqual(0, u[2][1].y);
        Assert.AreEqual(60, u[2][2].x);
        Assert.AreEqual(20, u[2][2].y);
        Assert.AreEqual(10, u[2][3].x);
        Assert.AreEqual(20, u[2][3].y);
        Assert.AreEqual(10, u[2][4].x);
        Assert.AreEqual(0,u[2][4].y);
        
        PathsD t = GeoWrangler.rectangular_decomposition(ref abort, T, maxRayLength: rayLength);
        writeToLayout("t", T, t);
        // Sign change expected.
        Assert.AreEqual(Clipper.Area(T), -Clipper.Area(t));
        Assert.AreEqual(3, t.Count);
        Assert.AreEqual(60, t[0][0].x);
        Assert.AreEqual(50, t[0][0].y);
        Assert.AreEqual(80, t[0][1].x);
        Assert.AreEqual(50, t[0][1].y);
        Assert.AreEqual(80, t[0][2].x);
        Assert.AreEqual(80, t[0][2].y);
        Assert.AreEqual(60, t[0][3].x);
        Assert.AreEqual(80, t[0][3].y);
        Assert.AreEqual(60, t[0][4].x);
        Assert.AreEqual(50, t[0][4].y);
        Assert.AreEqual(40, t[1][0].x);
        Assert.AreEqual(0, t[1][0].y);
        Assert.AreEqual(60, t[1][1].x);
        Assert.AreEqual(0, t[1][1].y);
        Assert.AreEqual(60, t[1][2].x);
        Assert.AreEqual(80, t[1][2].y);
        Assert.AreEqual(40, t[1][3].x);
        Assert.AreEqual(80, t[1][3].y);
        Assert.AreEqual(40, t[1][4].x);
        Assert.AreEqual(0, t[1][4].y);
        Assert.AreEqual(0, t[2][0].x);
        Assert.AreEqual(50, t[2][0].y);
        Assert.AreEqual(40, t[2][1].x);
        Assert.AreEqual(50, t[2][1].y);
        Assert.AreEqual(40, t[2][2].x);
        Assert.AreEqual(80, t[2][2].y);
        Assert.AreEqual(0, t[2][3].x);
        Assert.AreEqual(80, t[2][3].y);
        Assert.AreEqual(0, t[2][4].x);
        Assert.AreEqual(50, t[2][4].y);
        
        PathsD x = GeoWrangler.rectangular_decomposition(ref abort, X, maxRayLength: rayLength);
        writeToLayout("x", X, x);
        // Sign change expected.
        Assert.AreEqual(Clipper.Area(X), -Clipper.Area(x));
        Assert.AreEqual(3, x.Count);
        Assert.AreEqual(0, x[0][0].x);
        Assert.AreEqual(50, x[0][0].y);
        Assert.AreEqual(60, x[0][1].x);
        Assert.AreEqual(50, x[0][1].y);
        Assert.AreEqual(60, x[0][2].x);
        Assert.AreEqual(80, x[0][2].y);
        Assert.AreEqual(0, x[0][3].x);
        Assert.AreEqual(80, x[0][3].y);
        Assert.AreEqual(0, x[0][4].x);
        Assert.AreEqual(50, x[0][4].y);
        Assert.AreEqual(60, x[1][0].x);
        Assert.AreEqual(20, x[1][0].y);
        Assert.AreEqual(80, x[1][1].x);
        Assert.AreEqual(20, x[1][1].y);
        Assert.AreEqual(80, x[1][2].x);
        Assert.AreEqual(100, x[1][2].y);
        Assert.AreEqual(60, x[1][3].x);
        Assert.AreEqual(100, x[1][3].y);
        Assert.AreEqual(60, x[1][4].x);
        Assert.AreEqual(20, x[1][4].y);
        Assert.AreEqual(80, x[2][0].x);
        Assert.AreEqual(50, x[2][0].y);
        Assert.AreEqual(100, x[2][1].x);
        Assert.AreEqual(50, x[2][1].y);
        Assert.AreEqual(100, x[2][2].x);
        Assert.AreEqual(80, x[2][2].y);
        Assert.AreEqual(80, x[2][3].x);
        Assert.AreEqual(80, x[2][3].y);
        Assert.AreEqual(80, x[2][4].x);
        Assert.AreEqual(50, x[2][4].y);
        
        PathsD s = GeoWrangler.rectangular_decomposition(ref abort, S, maxRayLength: rayLength);
        writeToLayout("s", S, s);
        // Sign change expected.
        Assert.AreEqual(Clipper.Area(S), -Clipper.Area(s));
        Assert.AreEqual(5, s.Count);
        Assert.AreEqual(0, s[0][0].x);
        Assert.AreEqual(50, s[0][0].y);
        Assert.AreEqual(20, s[0][1].x);
        Assert.AreEqual(50, s[0][1].y);
        Assert.AreEqual(20, s[0][2].x);
        Assert.AreEqual(110, s[0][2].y);
        Assert.AreEqual(0, s[0][3].x);
        Assert.AreEqual(110, s[0][3].y);
        Assert.AreEqual(0, s[0][4].x);
        Assert.AreEqual(50, s[0][4].y);
        Assert.AreEqual(0, s[1][0].x);
        Assert.AreEqual(0, s[1][0].y);
        Assert.AreEqual(20, s[1][1].x);
        Assert.AreEqual(0, s[1][1].y);
        Assert.AreEqual(20, s[1][2].x);
        Assert.AreEqual(20, s[1][2].y);
        Assert.AreEqual(0, s[1][3].x);
        Assert.AreEqual(20, s[1][3].y);
        Assert.AreEqual(0, s[1][4].x);
        Assert.AreEqual(0, s[1][4].y);
        Assert.AreEqual(80, s[2][0].x);
        Assert.AreEqual(0, s[2][0].y);
        Assert.AreEqual(100, s[2][1].x);
        Assert.AreEqual(0, s[2][1].y);
        Assert.AreEqual(100, s[2][2].x);
        Assert.AreEqual(60, s[2][2].y);
        Assert.AreEqual(80, s[2][3].x);
        Assert.AreEqual(60, s[2][3].y);
        Assert.AreEqual(80, s[2][4].x);
        Assert.AreEqual(0, s[2][4].y);
        Assert.AreEqual(80, s[3][0].x);
        Assert.AreEqual(80, s[3][0].y);
        Assert.AreEqual(100, s[3][1].x);
        Assert.AreEqual(80, s[3][1].y);
        Assert.AreEqual(100, s[3][2].x);
        Assert.AreEqual(110, s[3][2].y);
        Assert.AreEqual(80, s[3][3].x);
        Assert.AreEqual(110, s[3][3].y);
        Assert.AreEqual(80, s[3][4].x);
        Assert.AreEqual(80, s[3][4].y);
        Assert.AreEqual(20, s[4][0].x);
        Assert.AreEqual(0, s[4][0].y);
        Assert.AreEqual(80, s[4][1].x);
        Assert.AreEqual(0, s[4][1].y);
        Assert.AreEqual(80, s[4][2].x);
        Assert.AreEqual(110, s[4][2].y);
        Assert.AreEqual(20, s[4][3].x);
        Assert.AreEqual(110, s[4][3].y);
        Assert.AreEqual(20, s[4][4].x);
        Assert.AreEqual(0, s[4][4].y);
        
        PathsD ns = GeoWrangler.rectangular_decomposition(ref abort, nS, maxRayLength: rayLength);
        writeToLayout("ns", nS, ns);
        // Sign change expected.
        Assert.AreEqual(Clipper.Area(nS), -Clipper.Area(ns));
        Assert.AreEqual(5, ns.Count);
        Assert.AreEqual(0, ns[0][0].x);
        Assert.AreEqual(-150, ns[0][0].y);
        Assert.AreEqual(20, ns[0][1].x);
        Assert.AreEqual(-150, ns[0][1].y);
        Assert.AreEqual(20, ns[0][2].x);
        Assert.AreEqual(-90, ns[0][2].y);
        Assert.AreEqual(0, ns[0][3].x);
        Assert.AreEqual(-90, ns[0][3].y);
        Assert.AreEqual(0, ns[0][4].x);
        Assert.AreEqual(-150, ns[0][4].y);
        Assert.AreEqual(0, ns[1][0].x);
        Assert.AreEqual(-200, ns[1][0].y);
        Assert.AreEqual(20, ns[1][1].x);
        Assert.AreEqual(-200, ns[1][1].y);
        Assert.AreEqual(20, ns[1][2].x);
        Assert.AreEqual(-180, ns[1][2].y);
        Assert.AreEqual(0, ns[1][3].x);
        Assert.AreEqual(-180, ns[1][3].y);
        Assert.AreEqual(0, ns[1][4].x);
        Assert.AreEqual(-200, ns[1][4].y);
        Assert.AreEqual(80, ns[2][0].x);
        Assert.AreEqual(-200, ns[2][0].y);
        Assert.AreEqual(100, ns[2][1].x);
        Assert.AreEqual(-200, ns[2][1].y);
        Assert.AreEqual(100, ns[2][2].x);
        Assert.AreEqual(-140, ns[2][2].y);
        Assert.AreEqual(80, ns[2][3].x);
        Assert.AreEqual(-140, ns[2][3].y);
        Assert.AreEqual(80, ns[2][4].x);
        Assert.AreEqual(-200, ns[2][4].y);
        Assert.AreEqual(80, ns[3][0].x);
        Assert.AreEqual(-120, ns[3][0].y);
        Assert.AreEqual(100, ns[3][1].x);
        Assert.AreEqual(-120, ns[3][1].y);
        Assert.AreEqual(100, ns[3][2].x);
        Assert.AreEqual(-90, ns[3][2].y);
        Assert.AreEqual(80, ns[3][3].x);
        Assert.AreEqual(-90, ns[3][3].y);
        Assert.AreEqual(80, ns[3][4].x);
        Assert.AreEqual(-120, ns[3][4].y);
        Assert.AreEqual(20, ns[4][0].x);
        Assert.AreEqual(-200, ns[4][0].y);
        Assert.AreEqual(80, ns[4][1].x);
        Assert.AreEqual(-200, ns[4][1].y);
        Assert.AreEqual(80, ns[4][2].x);
        Assert.AreEqual(-90, ns[4][2].y);
        Assert.AreEqual(20, ns[4][3].x);
        Assert.AreEqual(-90, ns[4][3].y);
        Assert.AreEqual(20, ns[4][4].x);
        Assert.AreEqual(-200, ns[4][4].y);
        
        PathsD c1 = GeoWrangler.rectangular_decomposition(ref abort, C1, maxRayLength: rayLength);
        writeToLayout("c1", C1, c1);
        // Sign change expected.
        Assert.AreEqual(Clipper.Area(C1), -Clipper.Area(c1));
        Assert.AreEqual(17, c1.Count);
        // Use area because result is complex and hash not reliable due to floats.
        Assert.AreEqual(4920, Clipper.Area(c1));

        PathsD c2 = GeoWrangler.rectangular_decomposition(ref abort, C2, maxRayLength: rayLength);
        writeToLayout("c2", C2, c2);
        // Sign change expected.
        Assert.AreEqual(Clipper.Area(C2), -Clipper.Area(c2));
        Assert.AreEqual(81, c2.Count);
        // Use area because result is complex and hash not reliable due to floats.
        Assert.AreEqual(24600, Clipper.Area(c2));
        
        PathsD c3 = GeoWrangler.rectangular_decomposition(ref abort, C3, maxRayLength: rayLength);
        writeToLayout("c3", C3, c3);
        // Sign change expected.
        Assert.AreEqual(Clipper.Area(C3), -Clipper.Area(c3));
        Assert.AreEqual(13, c3.Count);
        // Use area because result is complex and hash not reliable due to floats.
        Assert.AreEqual(5424, Clipper.Area(c3));
        
        PathsD c10r15 = GeoWrangler.rectangular_decomposition(ref abort, C10R15, maxRayLength: rayLength);
        writeToLayout("c10r15", C10R15, c10r15);
        // Sign change expected.
        Assert.AreEqual(Clipper.Area(C10R15), Clipper.Area(c10r15));
        Assert.AreEqual(1, c10r15.Count);
        Assert.AreEqual(-14742059.915, Clipper.Area(c10r15), 0.001);

        PathsD s1 = GeoWrangler.rectangular_decomposition(ref abort, S1, maxRayLength: rayLength);
        writeToLayout("s1", S1, s1);
        // Sign change expected.
        Assert.AreEqual(Clipper.Area(S1), -Clipper.Area(s1));
        Assert.AreEqual(4, s1.Count);
        Assert.AreEqual(-50, s1[0][0].x);
        Assert.AreEqual(-50, s1[0][0].y);
        Assert.AreEqual(0, s1[0][1].x);
        Assert.AreEqual(-50, s1[0][1].y);
        Assert.AreEqual(0, s1[0][2].x);
        Assert.AreEqual(0, s1[0][2].y);
        Assert.AreEqual(-50, s1[0][3].x);
        Assert.AreEqual(0, s1[0][3].y);
        Assert.AreEqual(-50, s1[0][4].x);
        Assert.AreEqual(-50, s1[0][4].y);
        Assert.AreEqual(0, s1[1][0].x);
        Assert.AreEqual(-50, s1[1][0].y);
        Assert.AreEqual(100, s1[1][1].x);
        Assert.AreEqual(-50, s1[1][1].y);
        Assert.AreEqual(100, s1[1][2].x);
        Assert.AreEqual(120, s1[1][2].y);
        Assert.AreEqual(0, s1[1][3].x);
        Assert.AreEqual(120, s1[1][3].y);
        Assert.AreEqual(0, s1[1][4].x);
        Assert.AreEqual(-50, s1[1][4].y);
        Assert.AreEqual(150, s1[2][0].x);
        Assert.AreEqual(-50, s1[2][0].y);
        Assert.AreEqual(200, s1[2][1].x);
        Assert.AreEqual(-50, s1[2][1].y);
        Assert.AreEqual(200, s1[2][2].x);
        Assert.AreEqual(300, s1[2][2].y);
        Assert.AreEqual(150, s1[2][3].x);
        Assert.AreEqual(300, s1[2][3].y);
        Assert.AreEqual(150, s1[2][4].x);
        Assert.AreEqual(-50, s1[2][4].y);
        Assert.AreEqual(100, s1[3][0].x);
        Assert.AreEqual(-50, s1[3][0].y);
        Assert.AreEqual(150, s1[3][1].x);
        Assert.AreEqual(-50, s1[3][1].y);
        Assert.AreEqual(150, s1[3][2].x);
        Assert.AreEqual(200, s1[3][2].y);
        Assert.AreEqual(100, s1[3][3].x);
        Assert.AreEqual(200, s1[3][3].y);
        Assert.AreEqual(100, s1[3][4].x);
        Assert.AreEqual(-50, s1[3][4].y);
    }

    private static void partTwo()
    {
        PathsD incoming = new();
        PathD lPieces = new ()
        {
            new(0.00000, 0.00000),
            new(0.00000, 0.05000),
            new(0.01000, 0.05000),
            new(0.01000, 0.00000),
            new(0.00000, 0.00000)
        };


        PathD lPiece2 = new ()
        {
            new(0.01000, 0.00000),
            new(0.01000, 0.02000),
            new(0.06000, 0.02000),
            new(0.06000, 0.00000),
            new(0.01000, 0.00000)
        };

        lPieces.Reverse();
        lPiece2.Reverse();
        incoming.Add(lPieces);
        incoming.Add(lPiece2);
        
        ClipperD c = new(Constants.roundingDecimalPrecision);
        c.AddSubject(incoming);
        PathsD ret = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, ret);
    }

    private static void partThree()
    {
        partThree_moderatelycomplex();
        partThree_morecomplex();
        partThree_evenmorecomplex();
        partThree_gainingacomplex();
        partThree_extremelycomplex();
    }
    
    private static void partThree_moderatelycomplex()
    {
        System.Diagnostics.Stopwatch sw = new();
        int rayLength = 1000;

        Console.WriteLine("  Preparing....");
        sw.Start();

        PathD poly = TestGeometry.moderatelycomplex();

        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Geo clean-up....");
        sw.Restart();
        poly = GeoWrangler.clockwiseAndReorderXY(poly);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Conversion....");
        sw.Restart();
        PathsD done = new() { new(poly) }; // this was originally scaled up by 1000 in the integer pipeline.
        done = Clipper.ScalePaths(done, 1000);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        bool abort = false;

        Console.WriteLine("  Decomposition (vertical)....");
        sw.Restart();
        PathsD ns = GeoWrangler.rectangular_decomposition(ref abort, done, maxRayLength: rayLength);
        
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Writing....");
        writeToLayout("moderatelycomplex", done[0], ns);

        Assert.AreEqual(81, ns.Count);
        // Sign change expected.
        Assert.AreEqual(Clipper.Area(done), -Clipper.Area(ns));
        
        Console.WriteLine("  Decomposition (horizontal)....");
        sw.Restart();
        ns = GeoWrangler.rectangular_decomposition(ref abort, done, maxRayLength: rayLength, vertical: false);

        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Writing....");
        writeToLayout("moderatelycomplex_horizontal", done[0], ns);

        Assert.AreEqual(81, ns.Count);
        // Sign change expected.
        Assert.AreEqual(Clipper.Area(done), -Clipper.Area(ns));

        Console.WriteLine("  Done.");
    }

    private static void partThree_morecomplex()
    {
        System.Diagnostics.Stopwatch sw = new();
        int rayLength = 1000;

        Console.WriteLine("  Preparing....");
        sw.Start();

        PathD poly = TestGeometry.morecomplex();

        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Geo clean-up....");
        sw.Restart();
        poly = GeoWrangler.clockwiseAndReorderXY(poly);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Conversion....");
        sw.Restart();
        PathsD done = new() { new(poly) }; // this was originally scaled up by 1000 in the integer pipeline.
        done = Clipper.ScalePaths(done, 1000);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        bool abort = false;

        Console.WriteLine("  Decomposition (vertical)....");
        sw.Restart();
        PathsD ns = GeoWrangler.rectangular_decomposition(ref abort, done, maxRayLength: rayLength);
        
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Writing....");
        writeToLayout("morecomplex", done[0], ns);

        Assert.AreEqual(161, ns.Count);
        // Sign change expected.
        Assert.AreEqual(Clipper.Area(done), -Clipper.Area(ns));
        
        Console.WriteLine("  Decomposition (horizontal)....");
        sw.Restart();
        ns = GeoWrangler.rectangular_decomposition(ref abort, done, maxRayLength: rayLength, vertical: false);

        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Writing....");
        writeToLayout("morecomplex_horizontal", done[0], ns);

        Assert.AreEqual(161, ns.Count);
        // Sign change expected.
        Assert.AreEqual(Clipper.Area(done), -Clipper.Area(ns));

        Console.WriteLine("  Done.");
    }
    private static void partThree_evenmorecomplex()
    {
        System.Diagnostics.Stopwatch sw = new();
        int rayLength = 1000;

        Console.WriteLine("  Preparing....");
        sw.Start();

        PathD poly = TestGeometry.evenmorecomplex();

        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Geo clean-up....");
        sw.Restart();
        poly = GeoWrangler.clockwiseAndReorderXY(poly);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Conversion....");
        sw.Restart();
        PathsD done = new() { new(poly) }; // this was originally scaled up by 1000 in the integer pipeline.
        done = Clipper.ScalePaths(done, 1000);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        bool abort = false;

        Console.WriteLine("  Decomposition (vertical)....");
        sw.Restart();
        PathsD ns = GeoWrangler.rectangular_decomposition(ref abort, done, maxRayLength: rayLength);
        
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Writing....");
        writeToLayout("evenmorecomplex", done[0], ns);

        Assert.AreEqual(241, ns.Count);
        // Sign change expected.
        Assert.AreEqual(Clipper.Area(done), -Clipper.Area(ns));
        
        Console.WriteLine("  Decomposition (horizontal)....");
        sw.Restart();
        ns = GeoWrangler.rectangular_decomposition(ref abort, done, maxRayLength: rayLength, vertical: false);

        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Writing....");
        writeToLayout("evenmorecomplex_horizontal", done[0], ns);

        Assert.AreEqual(241, ns.Count);
        // Sign change expected.
        Assert.AreEqual(Clipper.Area(done), -Clipper.Area(ns));

        Console.WriteLine("  Done.");
    }

    private static void partThree_gainingacomplex()
    {
        System.Diagnostics.Stopwatch sw = new();
        int rayLength = 1000;

        Console.WriteLine("  Preparing....");
        sw.Start();

        PathD poly = TestGeometry.gainingacomplex();

        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Geo clean-up....");
        sw.Restart();
        poly = GeoWrangler.clockwiseAndReorderXY(poly);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Conversion....");
        sw.Restart();
        PathsD done = new() { new(poly) }; // this was originally scaled up by 1000 in the integer pipeline.
        done = Clipper.ScalePaths(done, 1000);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        bool abort = false;

        Console.WriteLine("  Decomposition (vertical)....");
        sw.Restart();
        PathsD ns = GeoWrangler.rectangular_decomposition(ref abort, done, maxRayLength: rayLength);
        
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Writing....");
        writeToLayout("gainingacomplex", done[0], ns);

        Assert.AreEqual(401, ns.Count);
        // Sign change expected.
        Assert.AreEqual(Clipper.Area(done), -Clipper.Area(ns));
        
        Console.WriteLine("  Decomposition (horizontal)....");
        sw.Restart();
        ns = GeoWrangler.rectangular_decomposition(ref abort, done, maxRayLength: rayLength, vertical: false);

        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Writing....");
        writeToLayout("gainingacomplex_horizontal", done[0], ns);

        Assert.AreEqual(401, ns.Count);
        // Sign change expected.
        Assert.AreEqual(Clipper.Area(done), -Clipper.Area(ns));

        Console.WriteLine("  Done.");
    }
    private static void partThree_extremelycomplex()
    {
        System.Diagnostics.Stopwatch sw = new();
        int rayLength = 1000;

        Console.WriteLine("  Preparing....");
        sw.Start();

        PathD poly = TestGeometry.verycomplex();

        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Geo clean-up....");
        sw.Restart();
        poly = GeoWrangler.clockwiseAndReorderXY(poly);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Conversion....");
        sw.Restart();
        PathsD done = new() { new(poly) }; // this was originally scaled up by 1000 in the integer pipeline.
        done = Clipper.ScalePaths(done, 1000);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        bool abort = false;

        Console.WriteLine("  Decomposition (vertical)....");
        sw.Restart();
        PathsD ns = GeoWrangler.rectangular_decomposition(ref abort, done, maxRayLength: rayLength);
        
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Writing....");
        writeToLayout("extremelycomplex", done[0], ns);

        Assert.AreEqual(200, ns.Count);
        Assert.AreEqual(Clipper.Area(done), Clipper.Area(ns));
        
        Console.WriteLine("  Decomposition (horizontal)....");
        sw.Restart();
        ns = GeoWrangler.rectangular_decomposition(ref abort, done, maxRayLength: rayLength, vertical: false);

        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Writing....");
        writeToLayout("extremelycomplex_horizontal", done[0], ns);

        Assert.AreEqual(607, ns.Count);
        Assert.AreEqual(Clipper.Area(done), Clipper.Area(ns));

        Console.WriteLine("  Done.");
    }

    private static void partFour()
    {
        System.Diagnostics.Stopwatch sw = new();
        Console.WriteLine(" Part 1....");
        Console.WriteLine("  Preparing....");
        sw.Start();
        PathD points_1 = TestGeometry.ortho_fractal_1();
        points_1 = GeoWrangler.close(points_1);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        partFour_do(points_1, "complex_loop", 401, 401, 325380.0);

        Console.WriteLine(" Part 2....");
        Console.WriteLine("  Preparing....");
        sw.Start();
        PathD points_2 = TestGeometry.ortho_fractal_2();
        points_2 = GeoWrangler.close(points_2);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");
        partFour_do(points_2, "complex_loop_rot", 401, 401, 325380);
    }

    private static void partFour_do(PathD points, string baseString, int expectedCountV, int expectedCountH, double expectedArea)
    {
        System.Diagnostics.Stopwatch sw = new();

        bool vertical = true;
        bool abort = false;

        Console.WriteLine("  Keyhole....");
        // Give the keyholder a whirl:
        sw.Restart();
        PathD toDecomp = GeoWrangler.makeKeyHole(points, reverseEval:false, biDirectionalEval:true)[0];
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Query....");
        sw.Restart();
        PathD bounds = GeoWrangler.getBounds(toDecomp);
        PointD dist = GeoWrangler.distanceBetweenPoints_point(bounds[0], bounds[1]);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Decomposition (vertical)....");
        sw.Restart();
        PathsD decompOut = GeoWrangler.rectangular_decomposition(ref abort, toDecomp,
            maxRayLength: (long)Math.Max(Math.Abs(dist.x), Math.Abs(dist.y)) * 1, vertical: vertical);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Writing....");
        writeToLayout(baseString, points, decompOut);
        
        Assert.AreEqual(decompOut.Count, expectedCountV);
        Assert.AreEqual(Clipper.Area(decompOut), expectedArea);

        Console.WriteLine("  Decomposition (horizontal)....");
        sw.Restart();
        decompOut = GeoWrangler.rectangular_decomposition(ref abort, toDecomp,
            maxRayLength: (long)Math.Max(Math.Abs(dist.x), Math.Abs(dist.y)) * 1, vertical: !vertical);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Writing....");
        writeToLayout(baseString + "_horizontal", points, decompOut);

        Assert.AreEqual(decompOut.Count, expectedCountH);
        Assert.AreEqual(Clipper.Area(decompOut), expectedArea);
        Console.WriteLine("  Done.");

        sw.Stop();
    }

    private static void partFive()
    {
        bool vertical = true;
        bool abort = false;

        PathsD polydata = new()
        {
            new()
            {
                new (-40, -30),
                new (-40, 70),
                new (0, 70),
                new (0, 50),
                new (-20, 50),
                new (-20, -10),
                new (0, -10),
                new (0, -30),
                new (-40, -30),
            },
            new()
            {
                new (0, -30),
                new (0, -10),
                new (20, -10),
                new (20, 50),
                new (0, 50),
                new (0, 70),
                new (40, 70),
                new (40, -30),
                new (0, -30),
            },
            new()
            {
                new (-80, -60),
                new (-80, 100),
                new (0, 100),
                new (0, 80),
                new (-60, 80),
                new (-60, -50),
                new (0, -50),
                new (0, -60),
                new (-80, -60),
            },
            new()
            {
                new (0, -60),
                new (0, -50),
                new (60, -50),
                new (60, 80),
                new (0, 80),
                new (0, 100),
                new (80, 100),
                new (80, -60),
                new (0, -60),
            },
            new()
            {
                new (-8, -27),
                new (-8, 40),
                new (9, 40),
                new (9, -27),
                new (-8, -27)
            },
            new()
            {
                new (-14, -4),
                new (-14, 15),
                new (-7, 15),
                new (-7, -4),
                new (-14, -4),
            },
            new()
            {
                new (10, 9),
                new (10, 20),
                new (13, 20),
                new (13, 9),
                new (10, 9),
            },
            new()
            {
                new (48, -1),
                new (48, 31),
                new (55, 31),
                new (55, -1),
                new (48, -1),
            },
            new()
            {
                new (-11, -44),
                new (-11, -39),
                new (16, -39),
                new (16, -44),
                new (-11, -44),
            },
            new()
            {
                new (-51, 3),
                new (-51, 23),
                new (-47, 23),
                new (-47, 3),
                new (-51, 3),
            },
            new()
            {
                new (-16, 76),
                new (-16, 77),
                new (-3, 77),
                new (-3, 76),
                new (-16, 76),
            },
        };
        
        PathsD out_decomp = new();
        for (int i = 0; i < polydata.Count; i++)
        {
         PathD points = new (polydata[i]);
         points = GeoWrangler.removeDuplicates(points);
         points = GeoWrangler.stripCollinear(points);
         points = GeoWrangler.clockwiseAndReorderXY(points);
         
         PathD toDecomp = GeoWrangler.makeKeyHole(GeoWrangler.sliverGapRemoval(points), reverseEval:false, biDirectionalEval:false)[0];
         PathD  bounds = GeoWrangler.getBounds(toDecomp);
         PointD dist = GeoWrangler.distanceBetweenPoints_point(bounds[0], bounds[1]);

         PathsD decompOut = GeoWrangler.rectangular_decomposition(ref abort, toDecomp,
          maxRayLength: (long) Math.Max(Math.Abs(dist.x), Math.Abs(dist.y)), vertical: vertical);
         
         out_decomp.AddRange(decompOut);
        }
        
        Assert.AreEqual(out_decomp.Count, 19);
        Assert.AreEqual(Clipper.Area(out_decomp), 13843.0);
    }


    private static void writeToLayout(string filename, PathD orig, PathsD decomped)
    {
        // Can the system define geometry and write it correctly to Oasis and GDS files.
        GeoCore g = new();
        g.reset();
        GCDrawingfield drawing_ = new("")
        {
            accyear = 2018,
            accmonth = 12,
            accday = 5,
            acchour = 2,
            accmin = 10,
            accsec = 10,
            modyear = 2018,
            modmonth = 12,
            modday = 5,
            modhour = 2,
            modmin = 10,
            modsec = 10,
            databaseunits = 1000,
            userunits = 1E-3, // 0.001 / 1E-6;
            libname = "noname"
        };

        GCCell gcell = drawing_.addCell();
        gcell.accyear = 2018;
        gcell.accmonth = 12;
        gcell.accday = 5;
        gcell.acchour = 2;
        gcell.accmin = 10;
        gcell.accsec = 10;
        gcell.modyear = 2018;
        gcell.modmonth = 12;
        gcell.modday = 5;
        gcell.modhour = 2;
        gcell.modmin = 10;
        gcell.modsec = 10;

        gcell.cellName = "test";

        gcell.addPolygon(GeoWrangler.path64FromPathD(orig), 1, 0);

        for (int i = 0; i < decomped.Count; i++)
        {
            gcell.addBox(GeoWrangler.path64FromPathD(decomped[i]), 1, 1);
        }

        g.setDrawing(drawing_);
        g.setValid(true);

        gds.gdsWriter gw = new(g, "../../../../../decomp_out/" + filename + "_partitiontest.gds");
        gw.save();

        oasis.oasWriter ow = new(g, "../../../../../decomp_out/" + filename + "_partitiontest.oas");
        ow.save();
    }
}