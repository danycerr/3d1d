Point(1) = {0.15, 0.05, 0., 0. };
For i In {1:5:1}
Rotate {{0, 0, 1}, {0.05, 0.05, 0}, Pi/3} {
  Duplicata { Point{i}; }
}
EndFor
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {4, 1, 2, 3};
Line(6) = {1,6};
Line(7) = {6,5};
Line(8) = {5,4};

Line Loop(9) = {4, 6, 7, 8};

Transfinite Line {3, 4, 1, 2} = 4 Using Bump 1;
Transfinite Line {4,6,7,8} = 4 Using Bump 1;
Plane Surface(6) = {5};
Plane Surface(7) = {9};
Transfinite Surface {6};
Recombine Surface {6};
        surfaceVector[] = Extrude {0, 0, 1} {
         Surface{6};
         Layers{6};
         Recombine;
        };
Recombine Surface {7};
        surfaceVector[] = Extrude {0, 0, 1} {
         Surface{7};
         Layers{6};
         Recombine;
        };
        
Physical Volume(54) = {2, 1};
