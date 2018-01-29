
for alpha = 0.25
    i = 1;
    for t = 0:0.05:1
        true_solution(i) = t^8 - 3*t^(4+alpha/2)+9/4*t^alpha;
        i = i +1;
    end
    time = 0:0.05:1;
    plot(time, true_solution);
    hold on
end