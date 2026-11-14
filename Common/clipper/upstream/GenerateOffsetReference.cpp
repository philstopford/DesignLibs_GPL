// Build like GenerateReference.cpp (no ClipFileLoad needed), then run without arguments.
#include "clipper2/clipper.h"
#include <iostream>
using namespace Clipper2Lib;
static void paths(const Paths64& pp) {
  std::cout << '['; bool first = true;
  for (const auto& p : pp) {
    if (!first) std::cout << ','; first = false;
    std::cout << '['; bool firstPt = true;
    for (const auto& pt : p) {
      if (!firstPt) std::cout << ','; firstPt = false;
      std::cout << '[' << pt.x << ',' << pt.y << ']';
    }
    std::cout << ']';
  }
  std::cout << ']';
}
int main() {
  Paths64 sets[] = {
    {MakePath({0,0, 0,100, 100,100, 100,0})},
    {MakePath({0,0, 1000,0, 1000,1000, 0,1000}), MakePath({100,100, 100,200, 200,200, 200,100})},
    {MakePath({0,0, 100,0}), MakePath({0,1000, 500,1000, 500,1500})},
    {Path64{}, MakePath({500,500}), MakePath({0,0, 100,0, 100,100, 0,100})},
    {MakePath({0,0, 100,100, 0,100, 100,0})},
    {MakePath({-360,-100, 120,260, 300,-140, -340,160})},
    {Ellipse(Point64(0,0), 1000,800,80), MakePath({-100,-100, -100,100, 100,100, 100,-100})},
    {MakePath({33621,18600, 33889,18931, 33622,18601, 33118,18033}),
     MakePath({33058,3177, 33032,3337, 33017,3570, 33040,17609, 33118,18033, 32723,17588, 32698,3608, 32710,3376, 32737,3177})}
  };
  std::cout << '['; int n = 0;
  for (const auto& input : sets)
    for (int jt = 0; jt < 4; ++jt)
      for (int et = 0; et < 5; ++et)
        for (double delta : { -10.0, 0.0, 10.0 }) {
          ClipperOffset co(2.0, 0.25);
          co.AddPaths(input, JoinType(jt), EndType(et));
          Paths64 result; co.Execute(delta, result);
          if (n++) std::cout << ',';
          std::cout << "{\"number\":" << n << ",\"joinType\":" << jt << ",\"endType\":" << et << ",\"delta\":" << delta;
          std::cout << ",\"subject\":"; paths(input);
          std::cout << ",\"result\":"; paths(result);
          std::cout << '}';
        }
  std::cout << "]\n";
}
