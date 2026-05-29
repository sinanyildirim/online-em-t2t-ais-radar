function med = point_mediator(X1,X2,dist_ratio)
    vect_dir = (X2-X1)./norm(X2-X1);
    vect_dir_med = [-vect_dir(2); vect_dir(1)];
    midpoint = (X1+X2)/2;
    dist = dist_ratio*norm(X2-X1);
    med = midpoint + vect_dir_med*dist;
end