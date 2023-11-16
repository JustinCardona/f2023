function z = Newton(f, Dinvf, z0, tol)
    z = z0;
    e = f(z);
    while norm(e) > tol
        z = z - Dinvf(z)*e;
        e = f(z);
    end
end

