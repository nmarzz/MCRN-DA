function dx = FLor95(T, x)
% Evaluates the right hand side of the Lorenz '96 system
% \frac{dx_i}{dt} = -x_{i-2}x_{i-1} + x_{i-1}x_{i+1} - x_i + F
% F = 8, DIM=40

    F=8;
    DIM = 40;
    dx = zeros(DIM,1);

    for j=3:DIM-1
    dx(j) = -x(j)+x(j-1)*(x(j+1)-x(j-2)) + F;
    end
    dx(1) = -x(1)+x(DIM)*(x(2)-x(DIM-1)) + F;
    dx(2) = -x(2)+x(1)*(x(3)-x(DIM)) + F;
    dx(DIM) = -x(DIM)+x(DIM-1)*(x(1)-x(DIM-2)) + F;

    return
end
