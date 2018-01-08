Kn = 1;

r = .5/Kn;
R = 15*r;
Rx = R/2;
dx = R/4;

refine_r = 6; // 6--8
refine_c = 4;
N_R = 16; // 12--24
N_theta = N_R;

q_r = Exp(Log(refine_r*R/r)/N_R);
p_r = Exp(Log(refine_r*R/r/4)/N_R);
q_c = Exp(Log(refine_c)/N_theta);
Printf("q_r = %f, q_c = %f", q_r, q_c);
Printf("total cells = %.0f", N_R*N_theta);

Point(1) = {0, 0, 0};
Point(2) = {r, 0, 0};
Point(3) = {0, r, 0};
Point(4) = {R, 0, 0};
Point(5) = {dx, R, 0};
Point(6) = {-r, 0, 0};
Point(7) = {-Rx, 0, 0};
Point(10) = {dx, 0, 0};

Ellipse(1) = {2, 1, 3, 3};
Ellipse(2) = {3, 1, 3, 6};
Ellipse(3) = {4, 10, 5, 5};
Ellipse(4) = {5, 10, 5, 7};
Line(5) = {2, 4};
Line(6) = {6, 7};
Line(7) = {3, 5};

Line Loop(1) = {-1, 5, 3, -7};
Line Loop(2) = {7, 4, -6, -2};
Plane Surface(1) = {1};
Plane Surface(2) = {2};

Transfinite Line {1, -2} = N_theta+1;
Transfinite Line {3, -4} = N_theta+1 Using Progression q_c;
Transfinite Line {5} = N_R+1 Using Progression q_r;
Transfinite Line {6} = N_R+1 Using Progression p_r;
Transfinite Line {7} = N_R+1 Using Progression q_r;
Transfinite Surface {1} = {2, 4, 5, 3};
Transfinite Surface {2} = {5, 3, 6, 7};
Recombine Surface {1, 2};

Extrude {0, 0, 1/Kn} {
  Surface{1, 2};
  Layers{1};
  Recombine;
}

Physical Volume("internal") = {1, 2};
//Physical Surface("front") = {29, 51};
//Physical Surface("back") = {1, 2};
Physical Surface("body") = {16, 50};
Physical Surface("free") = {24, 42};
Physical Surface("symmetry") = {20, 46};
