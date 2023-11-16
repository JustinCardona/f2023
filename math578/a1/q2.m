function q2()
    fprintf("%d\n", BinarySearch(@(d) abs(d) * (1 + d^2) / (1 - d^2), 10, 0, 1, 1e-6));
end

