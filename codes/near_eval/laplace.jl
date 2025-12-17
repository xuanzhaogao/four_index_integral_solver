function laplace2d_DT(x, y, nx, ny)
    return (nx * x + ny * y) / (x^2 + y^2)
end

function laplace3d_DT(x, y, z, nx, ny, nz)
    return (nx * x + ny * y + nz * z) / (sqrt(x^2 + y^2 + z^2))^3
end