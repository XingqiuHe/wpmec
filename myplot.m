clear all;

mode = [3];
marker_list = ['^', 'o', 's', 'x', '+', 'v'];

if ismember(1, mode)
  load('data/MAP_5_30_1.mat');

  figure('Position', [200 200 450 350]);
  hold on; grid on;

  scatter(locAP(:,1), locAP(:,2), 80, '^', 'filled', 'DisplayName', 'AP');

  scatter(locWD(:,1), locWD(:,2), 40, 'o', 'filled', 'DisplayName', 'WD');

  xlabel('X Coordinate (m)');
  ylabel('Y Coordinate (m)');
  legend('Location', 'east');
  axis equal;
  axis square;

  xlim([0 10]);
  ylim([0 10]);
  xticks(0:2:10);
  yticks(0:2:10);

  set(gca, 'FontName', 'Times New Roman');  
  set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

  filename = "fig/location";
  exportgraphics(gcf, filename + ".png", 'Resolution', 600);
  exportgraphics(gcf, filename + ".eps");
end

if ismember(2, mode)
  load('data/queue_dynamics.mat');

  % plot Q dynamics
  figure('Position', [200 200 485 300]);
  hold on; grid on;

  y = mean(Q_PH, 2);
  x = 1:length(y);
  lw = 0.8;
  plot(x, y, 'linewidth', lw, 'DisplayName', '$Q(t)$ w/ PH');
  plot(x, mean(Q_act_PH, 2), 'linewidth', lw, 'DisplayName', '$Q^{act}(t)$ w/ PH');
  plot(x, mean(q_PH, 2), 'linewidth', lw, 'DisplayName', '$q(t)$ w/ PH');
  plot(x, mean(Q, 2), 'linewidth', lw, 'DisplayName', '$Q(t)$ w/o PH');

  xlabel('Time Slot');
  ylabel('Queue Length (bit)');
  lgd = legend('Location', 'northeast');
  set(lgd, 'Interpreter', 'latex');

  set(gca, 'FontName', 'Times New Roman');  
  set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

  filename = "fig/Q_dynamics";
  exportgraphics(gcf, filename + ".png", 'Resolution', 600);
  exportgraphics(gcf, filename + ".eps");
  
  % % plot B dynamics
  % figure('Position', [200 200 485 300]);
  % hold on; grid on;
  %
  % plot(x, mean(B_remain_PH, 2), 'DisplayName', '$B(t)$ w/ PH');
  % plot(x, mean(B_remain, 2), 'DisplayName', '$B(t)$ w/o PH');
  %
  % xlabel('Time Slot');
  % ylabel('Remaining Energy (J)');
  % lgd = legend('Location', 'northeast');
  % set(lgd, 'Interpreter', 'latex');
  %
  % % set(gca, 'FontName', 'Times New Roman');  
  % % set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');
  %
  % filename = "fig/B_dynamics";
  % exportgraphics(gcf, filename + ".png", 'Resolution', 600);
  % exportgraphics(gcf, filename + ".eps");

  % plot heatmap
  figure('Position', [200 200 485 300]);
  hold on; grid on;

  avg_B_remain = mean(B_remain_PH, 1);
  [~, idx] = sort(avg_B_remain);
  B_remain_sorted = B_remain_PH(:, idx);

  [M, N] = size(B_remain_sorted);
  x_edges = (1:M)-0.5;
  y_edges = (1:N)-0.5;

  imagesc(x_edges, y_edges, B_remain_sorted');
  colormap('viridis');
  colorbar;

  xlabel('Time Slot');
  ylabel('WD Index');

  set(gca, 'FontName', 'Times New Roman');  
  set(gca, 'YDir', 'normal');
  set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

  % yticks(1:N+1);
  % xticks(1:M+1);

  filename = "fig/B_heatmap";
  exportgraphics(gcf, filename + ".png", 'Resolution', 600);
  exportgraphics(gcf, filename + ".eps");
end

if ismember(3, mode)
  load('data/different_V.mat');
  alg_list = ["Prop", "LCO", "FO", "Greedy"];
  name_list = ["Prop", "LCO", "FO", "Myopic"];
  len = length(seed_list);

  figure('Position', [200 200 485 300]);
  hold on; grid on;
 
  for i = 1:4
    y = L_data.(alg_list(i));
    y = mean(reshape(y, len, []), 1);
    plot(V_list, y, 'DisplayName', name_list(i), 'marker', marker_list(i), 'linewidth', 1.0);
  end

  lgd = legend('Location', 'north');
  lgd.NumColumns = 4;

  set(gca, 'YScale', 'log');
  ylim([50, 20000]);
  xticks(5000:2000:15000);
  xlabel('V');
  ylabel('Latency (ms)');

  set(gca, 'FontName', 'Times New Roman');  
  set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

  filename = "fig/V_latency";
  exportgraphics(gcf, filename + ".png", 'Resolution', 600);
  exportgraphics(gcf, filename + ".eps");

  figure('Position', [200 200 485 300]);
  hold on; grid on;

  for i = 1:4
    y = E_data.(alg_list(i));
    y = mean(reshape(y, len, []), 1)*1000;
    plot(V_list, y, 'DisplayName', name_list(i), 'marker', marker_list(i), 'linewidth', 1.0);
  end
  lgd = legend('Location', 'north');
  lgd.NumColumns = 4;

  xticks(5000:2000:15000);
  xlabel('V');
  ylabel('Energy Consumption (mJ)');

  set(gca, 'FontName', 'Times New Roman');  
  set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');
  
  filename = "fig/V_energy";
  exportgraphics(gcf, filename + ".png", 'Resolution', 600);
  exportgraphics(gcf, filename + ".eps");
end

if ismember(4, mode)
  load('data/arrival_rate.mat');
  alg_list = ["Prop", "LCO", "FO", "Greedy"];
  name_list = ["Prop", "LCO", "FO", "Myopic"];
  len = length(seed_list);

  figure('Position', [200 200 485 300]);
  hold on; grid on;

  for i = 1:4
    y = L_data.(alg_list(i));
    y = reshape(y, len, []);
    y = sort(y, 1);
    y = mean(y(2:end-1, :), 1);
    plot(discount_list*1500, y, 'DisplayName', name_list(i), 'marker', marker_list(i), 'linewidth', 1.0);
  end

  lgd = legend('Location', 'northwest');

  xlim([150, 1500]);
  xticks(150:150:1500);
  set(gca, 'YScale', 'log');
  xlabel('Average Arrival Rate (bits/slot)');
  ylabel('Latency (ms)');

  set(gca, 'FontName', 'Times New Roman');  
  set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

  filename = "fig/arrival_rate_latency";
  exportgraphics(gcf, filename + ".png", 'Resolution', 600);
  exportgraphics(gcf, filename + ".eps");

  figure('Position', [200 200 485 300]);
  hold on; grid on;

  for i = 1:4
    y = E_data.(alg_list(i));
    y = reshape(y, len, []);
    y = sort(y, 1);
    y = mean(y(2:end-1, :), 1)*1000;
    plot(discount_list*1500, y, 'DisplayName', name_list(i), 'marker', marker_list(i), 'linewidth', 1.0);
  end
  lgd = legend('Location', 'northwest');

  xlim([150, 1500]);
  xticks(150:150:1500);
  xlabel('Average Arrival Rate (bits/slot)');
  ylabel('Energy Consumption (mJ)');

  set(gca, 'FontName', 'Times New Roman');  
  set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');
  
  filename = "fig/arrival_rate_energy";
  exportgraphics(gcf, filename + ".png", 'Resolution', 600);
  exportgraphics(gcf, filename + ".eps");
end

if ismember(5, mode)
  load('data/different_M.mat');
  alg_list = ["Prop", "LCO", "FO", "Greedy"];
  name_list = ["Prop", "LCO", "FO", "Myopic"];
  len = length(seed_list);

  figure('Position', [200 200 485 300]);
  hold on; grid on;

  for i = 1:4
    y = L_data.(alg_list(i));
    y = reshape(y, len, []);
    y = sort(y, 1);
    y = mean(y(2:end-1, :), 1);
    plot(M_list, y, 'DisplayName', name_list(i), 'marker', marker_list(i), 'linewidth', 1.0);
  end

  lgd = legend('Location', 'best');

  xlim([1, 10]);
  set(gca, 'YScale', 'log');
  xlabel('Number of APs');
  ylabel('Latency (ms)');

  set(gca, 'FontName', 'Times New Roman');  
  set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

  filename = "fig/M_latency";
  exportgraphics(gcf, filename + ".png", 'Resolution', 600);
  exportgraphics(gcf, filename + ".eps");

  figure('Position', [200 200 485 300]);
  hold on; grid on;

  for i = 1:4
    y = E_data.(alg_list(i));
    y = reshape(y, len, []);
    y = sort(y, 1);
    y = mean(y(2:end-1, :), 1)*1000;
    plot(M_list, y, 'DisplayName', name_list(i), 'marker', marker_list(i), 'linewidth', 1.0);
  end
  lgd = legend('Location', 'best');

  xlim([1, 10]);
  xlabel('Number of APs');
  ylabel('Energy Consumption (mJ)');

  set(gca, 'FontName', 'Times New Roman');  
  set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');
  
  filename = "fig/M_energy";
  exportgraphics(gcf, filename + ".png", 'Resolution', 600);
  exportgraphics(gcf, filename + ".eps");
end

if ismember(6, mode)
  load('data/different_N.mat');
  alg_list = ["Prop", "LCO", "FO", "Greedy"];
  name_list = ["Prop", "LCO", "FO", "Myopic"];
  len = length(seed_list);

  figure('Position', [200 200 485 300]);
  hold on; grid on;

  for i = 1:4
    y = L_data.(alg_list(i));
    y = reshape(y, len, []);
    y = sort(y, 1);
    y = mean(y(2:end-1, :), 1);
    plot(N_list, y, 'DisplayName', name_list(i), 'marker', marker_list(i), 'linewidth', 1.0);
  end

  lgd = legend('Location', 'southeast');

  lgd.Units = 'normalized';
  pos = lgd.Position;
  pos(2) = pos(2) + 0.1;
  lgd.Position = pos;

  ylim([100, 20000]);

  set(gca, 'YScale', 'log');
  xlabel('Number of WDs');
  ylabel('Latency (ms)');

  set(gca, 'FontName', 'Times New Roman');  
  set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

  filename = "fig/N_latency";
  exportgraphics(gcf, filename + ".png", 'Resolution', 600);
  exportgraphics(gcf, filename + ".eps");

  figure('Position', [200 200 485 300]);
  hold on; grid on;

  for i = 1:4
    y = E_data.(alg_list(i));
    y = reshape(y, len, []);
    y = sort(y, 1);
    y = mean(y(2:end-1, :), 1)*1000;
    plot(N_list, y, 'DisplayName', name_list(i), 'marker', marker_list(i), 'linewidth', 1.0);
  end
  lgd = legend('Location', 'best');

  xlabel('Number of WDs');
  ylabel('Energy Consumption (mJ)');

  set(gca, 'FontName', 'Times New Roman');  
  set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');
  
  filename = "fig/N_energy";
  exportgraphics(gcf, filename + ".png", 'Resolution', 600);
  exportgraphics(gcf, filename + ".eps");
end
