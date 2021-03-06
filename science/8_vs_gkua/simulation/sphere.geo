//Kn = 1;

r = .5/Kn;
R = 10*r;

refine_r = 10; // 6 -- 8
refine_c = 4;
N_R = 24; // 12--24
N_theta = N_R;
N_phi = N_R/4;

q_r = Exp(Log(refine_r*R/r)/N_R);
q_c = Exp(Log(refine_c)/N_theta);
Printf("q_r = %f, q_c = %f", q_r, q_c);
Printf("total cells = %.0f", N_R*N_theta*N_phi);

Point(1) = {0, 0, 0};
Point(2) = {r, 0, 0};
Point(3) = {0, r, 0};
Point(4) = {R, 0, 0};
Point(5) = {0, R, 0};
Point(6) = {-r, 0, 0};
Point(7) = {-R, 0, 0};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 6};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 7};
Line(5) = {2, 4};
Line(6) = {6, 7};

Line Loop(5) = {-1, 5, 3, 4, -6, -2};
Plane Surface(6) = {5};

Transfinite Line {1, -2, 3, -4} = N_theta+1 Using Progression q_c;
Transfinite Line {5, 6} = N_R+1 Using Progression q_r;
Transfinite Surface {6} = {2,4,6,7};
Recombine Surface {6};

Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {
    Surface{6};
    Layers{N_phi};
    Recombine;
}

Physical Volume("internal") = {1};
Physical Surface("bottom") = {28};
Physical Surface("back") = {6};
Physical Surface("body") = {16,27};
Physical Surface("free") = {20,23};
