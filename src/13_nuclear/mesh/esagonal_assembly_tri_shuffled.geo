Point(1) = {0.15, 0.05, 0., 0. };
For i In {1:5:1}
Rotate {{0, 0, 1}, {0.05, 0.05, 0}, Pi/3} {
  Duplicata { Point{i}; }
}
EndFor
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};
Line Loop(7) = {4, 5, 6, 1, 2, 3};
Plane Surface(8) = {7};
Extrude {0, 0, 1} {
  Surface{8};
}
Transfinite Line {15, 14, 13, 12, 11, 10, 4, 5, 6, 1, 2, 3} = 3 Using Progression 1;
Transfinite Line {13, 1} = 4 Using Progression 1;
Physical Volume(41) = {1};
