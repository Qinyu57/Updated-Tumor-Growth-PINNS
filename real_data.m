data = load('demo4.mat');
x_coords = [-2, -1.5, -1, 0, 1, 1.5, 2];
y_coords = [-2, -1.5, -1, 0, 1, 1.5, 2];
rho_values = zeros(7,7);

for i = 1:7
    for j = 1:7
        % Find index where x ≈ x_coords(i) and y ≈ y_coords(j)
        idx = find(abs(data.X - x_coords(i)) < 0.01 & abs(data.Y - y_coords(j)) < 0.01);
        if ~isempty(idx)
            rho_values(i,j) = data.rho(idx(1));
        end
    end
end

% Write to Excel
xlswrite('rho_data.xlsx', rho_values);