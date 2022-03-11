// See https://aka.ms/new-console-template for more information

using ClipperLib1Test;
using ClipperLib2Test;
using ClipperLibTest;

Console.WriteLine("Comparison tests for Clipper 1 and Clipper 2");

Clipper2Test.colinearTest();

// Compare open line clipping
ComparisonTest.lineClipTest1();
ComparisonTest.lineClipTest2();

Clipper1Test.unionTest();
Clipper2Test.unionTest();

Clipper1Test.notTest();
Clipper2Test.notTest();

Clipper1Test.edgeOffsetTest();
Clipper2Test.edgeOffsetTest();

Clipper1Test.rightChordTest();
Clipper2Test.rightChordTest();

Clipper1Test.leftChordTest();
Clipper2Test.leftChordTest();

Clipper1Test.zFillCallbackTest();
Clipper2Test.zFillCallbackTest();

Clipper1Test.coincident_openPathTest();
Clipper2Test.coincident_openPathTest();

Clipper1Test.keyHole_test2();
Clipper2Test.keyHole_test2();

Clipper1Test.openPath_clipTest1();
Clipper2Test.openPath_clipTest1();

Clipper1Test.openPath_clipTest2();
Clipper2Test.openPath_clipTest2();

Clipper1Test.offsetTest();
Clipper2Test.offsetTest();

Clipper1Test.keyHole_test1();
Clipper2Test.keyHole_test1();
