//Kn = 1;

r = .5/Kn;
R = 15*r;
Rx = R/3;
dx = (R-Rx)/2;
sw = 2.15*r;

refine_r = 8; // 6--8
refine_c = 2;
N_R = 48; // 24--48
N_theta = Round(N_R/3);

outer = .5*N_R;
inner = .8*N_R;
delta = Round(N_theta/8);

q_r = Exp(Log(refine_r*R/r*12)/N_R);
q_R = Exp(Log(refine_c)/N_theta);
q_sw = Exp(Log(refine_c)/N_theta);
Printf("q_r = %f, q_R = %f, delta = %.0f", q_r, q_R, delta);
Printf("total cells = %.0f", (inner+outer)*N_theta);

Point(1) = {0, 0, 0};
Point(2) = {r, 0, 0};
Point(3) = {0, r, 0};
Point(4) = {R, 0, 0};
Point(5) = {dx, R, 0};
Point(6) = {-r, 0, 0};
Point(7) = {-Rx, 0, 0};
Point(10) = {dx, 0, 0};
Point(11) = {.9*sw, 0, 0};
Point(12) = {-sw, 0, 0};
Point(13) = {.9*sw, 2.2*sw, 0};
Point(14) = {2*sw, 0, 0};

Ellipse(1) = {2, 1, 3, 3};
Ellipse(2) = {3, 1, 3, 6};
Ellipse(3) = {4, 10, 5, 5};
Ellipse(4) = {5, 10, 5, 7};
Ellipse(5) = {12, 11, 13, 13};
Ellipse(6) = {13, 11, 13, 14};
Line(10) = {14, 2};
Line(11) = {14, 4};
Line(12) = {6, 12};
Line(13) = {7, 12};

Line Loop(1) = {1, 10, 6, 5, 12, 2};
Line Loop(2) = {6, 11, 3, 4, 13, 5};
Plane Surface(1) = {1};
Plane Surface(2) = {2};

Transfinite Line {1, -2} = N_theta+1;                                               // inner back+front
Transfinite Line {3} = N_theta-delta*0+1 Using Progression q_R;                     // outer back
Transfinite Line {-4} = N_theta+delta*0+1 Using Progression q_R;                    // outer front
Transfinite Line {-6} = N_theta-delta+1 Using Progression q_sw;                     // sw back
Transfinite Line {5} = N_theta+delta+1 Using Progression q_sw;                      // sw front
Transfinite Line {12} = inner+1 Using Bump 0.2;                                     // sw--inner front
Transfinite Line {-13} = outer+1 Using Progression q_r;                             // sw--outer front
Transfinite Line {11} = outer+1 Using Progression Exp(Log(refine_r*R/r/10)/N_R);    // sw--outer back
Transfinite Line {-10} = inner+1 Using Progression Exp(Log(refine_r*R/r/5)/N_R);    // sw--inner back
Transfinite Surface {1} = {2, 14, 12, 6};
Transfinite Surface {2} = {14, 4, 7, 12};
Recombine Surface {1, 2};

Extrude {0, 0, 1/Kn} {
  Surface{1, 2};
  Layers{1};
  Recombine;
}

Physical Volume("internal") = {1, 2};
//Physical Surface("front") = {45, 77};
//Physical Surface("back") = {1, 2};
Physical Surface("body") = {24, 28};
Physical Surface("free") = {64, 68};
Physical Surface("symmetry") = {32, 44, 60, 72};

plot_image = 0;
If (plot_image == 1)
    Merge "Ma-70km.png";
    c = 1.15;
    scale = R/(View[0].MaxY)/c;
    View[0].TransformYY = scale;
    View[0].TransformXX = scale;
    View[0].OffsetX = -R/c/2.848;
EndIf
