// See https://aka.ms/new-console-template for more information

using ClipperLib1Test;
using ClipperLib2Test;

Console.WriteLine("Comparison tests for Clipper 1 and Clipper 2");

Clipper1Test.test5();
Clipper2Test.test5();

Clipper1Test.test4();
Clipper2Test.test4();

Clipper1Test.test3();
Clipper2Test.test3();

Clipper1Test.test2();
Clipper2Test.test2();

Clipper1Test.offsetTest();
Clipper2Test.offsetTest();

Clipper1Test.test1();
Clipper2Test.test1();
