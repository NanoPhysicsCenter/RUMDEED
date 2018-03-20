function [ x_0, y_0, z_0 ] = Find_closest_point( x_a, ~, z_a, eta_1, a )
x_0 = 0.0d0;
y_0 = 0.0d0;

s = 0.25;
for i = 1:100
    z_0 = eta_1/sqrt(1-eta_1^2)*sqrt(x_0^2+a^2*(1-eta_1^2));
    b = 1/(z_0 - z_a) * sqrt(1-eta_1^2)/eta_1 * sqrt(x_0^2 + a^2*(1-eta_1^2));
    x_new = x_a * b/(1+b);
    
    x_old = x_0;
    x_0 = s*x_new + (1-s)*x_old;
    
    if (abs((x_0 - x_new)/x_new) < 1E-4)
        break
    end
end

%z_0 = -z_0;
end

