% validation

nodeloc = importdata('B737_data/B737.inp');
nodeloc_x = nodeloc.data(:,2);
nodeloc_y = nodeloc.data(:,3);
nodeloc_z = nodeloc.data(:,4);
scatter3(nodeloc_x,nodeloc_y,nodeloc_z);

% fid = fopen('B737_data/B737.rpt');
% % data = importdata('B737_data/B737.rpt','%f');
% % fclose('all');
% % data_region1 = data(12:5792,:);
