function I = gaussian3d(x, y, z, center, sigma)
% I = gaussian3d(x, y, z, [x0, y0, z0], sigma)
% I = gaussian3d(x, y, z, [x0, y0, z0], [sigmax, sigmay, sigmaz])

% https://math.stackexchange.com/questions/434629/3-d-generalization-of-the-gaussian-point-spread-function
x = x-center(1);
y = y-center(2);
z = z-center(3);

if numel(sigma) == 1
    N = 1/(sigma^3*(2*pi)^(3/2));
    I = N*exp(-(x.^2+y.^2+z.^2)/(2*sigma^2));
elseif numel(sigma) == 3
    N = 1/(sigma(1)*sigma(2)*sigma(3)*(2*pi)^(3/2));
    I = exp(-x.^2/(2*sigma(1)^2));
    I = I.*exp(-y.^2/(2*sigma(2)^2));
    I = I.*exp(-z.^2/(2*sigma(3)^2));
    I = N*I;
end