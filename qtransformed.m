function qd = qtransform(x,Q)
R = [1 0 0;
     0 1 0;
     0 0 2];
 T = [cos(x)^2 sin(x)^2 2*sin(x)*cos(x);
      sin(x)^2 cos(x)^2 -2*sin(x)*cos(x);
      -sin(x)*cos(x) sin(x)*cos(x) cos(x)^2 - sin(x)^2];
  qd = T\Q*R*T/R;
end

  