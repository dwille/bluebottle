load part_data.cgns

% plot force
figure
plot(time, FZ)
hold on
plot(time, mean(FZ), 'k', 'LineWidth', 2)
plot(time, mean(FZ) + std(FZ), 'r', 'LineWidth', 2)
plot(time, mean(FZ) - std(FZ), 'r', 'LineWidth', 2)
title('F_z')
xlabel('Time')
ylabel('F_z')

% plot velocity
figure
plot(time, Wp)
hold on
plot(time, mean(Wp), 'k', 'LineWidth', 2)
plot(time, mean(Wp) + std(Wp), 'r', 'LineWidth', 2)
plot(time, mean(Wp) - std(Wp), 'r', 'LineWidth', 2)
title('W')
xlabel('Time')
ylabel('W')

% plot position
figure
plot(time, Zp)
hold on
plot(time, mean(Zp), 'k', 'LineWidth', 2)
plot(time, mean(Zp) + std(Zp), 'r', 'LineWidth', 2)
plot(time, mean(Zp) - std(Zp), 'r', 'LineWidth', 2)
title('Z')
xlabel('Time')
ylabel('Z')
