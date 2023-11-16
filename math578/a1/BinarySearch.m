function x  = BinarySearch(f, goal, min, max, tol)
    x = (min+max)/2;
    y = f(x);
    while abs(y-goal) > tol
        if y < goal
            min = x;
        else
            max = x;
        end
        x = (min+max)/2;
        y = f(x);
    end
end