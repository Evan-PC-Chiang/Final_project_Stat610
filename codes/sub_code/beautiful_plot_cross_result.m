close all
data_type = data_matrix(:,1);

figure()
t = tiledlayout(7, 1, 'TileSpacing', 'none', 'Padding', 'compact');


ax4 = nexttile(4,[4 1]);
hold on
q = quiver(GF_all.Xobs,-GF_all.Zobs,model_st_horz_mean,model_lt_vert_mean);
q.Color = '#58B2DC';
fault = GF.Fault;
count = sum(data.Fault(:,5) ~= max(data.Fault(:,5)));
seg_data = zeros(count,4);

for k=1:count
   a=[fault(k,2) fault(k,2)+fault(k,3)*cos(fault(k,4)*pi/180)]; 
   b=-[fault(k,1) fault(k,1)+fault(k,3)*sin(fault(k,4)*pi/180)];
   seg_data(k,:) = [a(1) b(1) a(2) b(2)];
   plot(a,b,'r','LineWidth',1.5)
   text(mean(a),mean(b),sprintf('%.2f',locking_mean(k)),'Fontsize',12, ...
       'HorizontalAlignment', 'center','VerticalAlignment', 'middle','FontName','Helvetica')
end

[unique_vals, ~, group_indices] = unique(round(seg_slip_mean(1:count),6));
for j=1:length(unique_vals)
    segments = seg_data(group_indices == j,:);
    points = round([segments(:, 1:2); segments(:, 3:4)],0); % [x1 y1; x2 y2]
    [a,b] = ismember(0,points(:,2));
    if ismember(0,points(:,2))
        text_y = 2;
        if height(points) >4
            text_x = mean([points(b,1) points(b+2,1)]);
        else
            text_x = points(b,1);
        end
    else
        text_y = -22;
        unique_points = unique(points, 'rows'); % Find unique points
        % Count occurrences of each unique point
        counts = zeros(size(unique_points, 1), 1);
            for i = 1:size(unique_points, 1)
                counts(i) = sum(ismember(points, unique_points(i, :), 'rows'));
            end
        start_end_points = unique_points(counts == 1, :);
        text_x = mean(start_end_points(:,1));
    end
        % Size and padding for the rounded rectangle
    padding = 1;
    widthh = 3.8;
    heightt = 1.8;

    %text_y = mean(start_end_points(:,2));
    %Add a rounded rectangle as background
    rectangle('Position', [text_x - widthh / 2, text_y - heightt / 2, widthh, heightt], ...
              'Curvature', [0.08, 0.08], ... % Rounded corners (1 = full circle)
              'FaceColor', [0.9 0.9 0.9], ... % Light gray background
              'EdgeColor', 'black');
    text(text_x+0.05,text_y,sprintf('%.2f',unique_vals(j)),'FontSize',16, ...
        'HorizontalAlignment', 'center','VerticalAlignment', 'middle','FontName','Helvetica')
end

text(-67,-22,'Invicid Mantle','FontName','Helvetica','FontSize',22)
text(-67,-18,'Elastic Crust','FontName','Helvetica','FontSize',22)
text(-68,-8.5,'Eurasia Plate','FontName','Helvetica','FontSize',20, ...
    'Rotation',90,'HorizontalAlignment', 'center','VerticalAlignment', 'top')
text(-66.5,-8.5,'(Fixed)','FontName','Helvetica','FontSize',20, ...
    'Rotation',90,'HorizontalAlignment', 'center','VerticalAlignment', 'top')
text(69,-10,'Philippine Sea','FontName','Helvetica','FontSize',20, ...
    'Rotation',270,'HorizontalAlignment', 'center','VerticalAlignment', 'top')
text(67.5,-10,'Plate (Moving)','FontName','Helvetica','FontSize',20, ...
    'Rotation',270,'HorizontalAlignment', 'center','VerticalAlignment', 'top')

if operator == 1
    GF.Build_gemo.plot_Shear
else
end
x1 = min(-68); x2 = max(70);
plot([x1 x2], [0 0],'k','LineWidth',2)
plot([x1 x2], [-max(fault(:,1)) -max(fault(:,1))],'k','LineWidth',2)

xlim([x1,x2])
xlabel('Distance along profile (km)')
ylabel('Depth (km)','FontSize',16)
pbaspect([x2-x1 28 1]);
ylim([data.H-4 4]);

%
% Customize ticks
grid on
ax = gca; % Get current axes
ax.LineWidth = 2;
ax.Box = 'on'; % Remove top and right border
ax.GridLineStyle = '--';
ax.XAxis.TickDirection = 'in';
ax.YAxis.TickDirection = 'out';
ax.TickLength = [0.002 0.005];
ax.FontSize = 16;
ax.FontName = 'Helvetica';

ax1 = nexttile(1);
plot(GF_all.xobs,model_st_horz_mean(1,:),'blue','LineWidth',1.5)
hold on
errorbar(GF.xobs(data_type==1),data_matrix(data_type==1,3),2*data_matrix(data_type==1,6),'.','MarkerSize',12)
xlim([x1,x2])
ylim([-65 5])
ylabel(' ','FontSize',16)
text(-67,-62,'Short-term East-West','FontName','Helvetica','FontSize',22, ...
    'Rotation',0,'HorizontalAlignment', 'left','VerticalAlignment', 'bottom')
grid on
ax = gca; % Get current axes
ax.LineWidth = 2;
ax.Box = 'on'; % Remove top and right border
ax.GridLineStyle = '--';
ax.XAxis.TickDirection = 'in';
ax.YAxis.TickDirection = 'out';
ax.TickLength = [0.002 0.005];
ax.FontSize = 16;
ax.FontName = 'Helvetica';
ax.YTick = ([-60 -40 -20 0]);
pbaspect([(x2-x1)*2 14 1]);

ax2 = nexttile(2);
plot(GF_all.xobs,model_st_vert_mean(1,:),'blue','LineWidth',1.5)
hold on
errorbar(GF.xobs(data_type==1),data_matrix(data_type==1,5),2*data_matrix(data_type==1,8),'.','MarkerSize',12)
xlim([x1,x2])
ylim([-20 30])
ylabel('Velocity (mm/yr)','FontSize',16)
text(-67,-19,'Short-term Vertical','FontName','Helvetica','FontSize',22, ...
    'Rotation',0,'HorizontalAlignment', 'left','VerticalAlignment', 'bottom')
%title('Short-term Vertical')
grid on
ax = gca; % Get current axes
ax.LineWidth = 2;
ax.Box = 'on'; % Remove top and right border
ax.GridLineStyle = '--';
ax.XAxis.TickDirection = 'in';
ax.YAxis.TickDirection = 'out';
ax.TickLength = [0.002 0.005];
ax.FontSize = 16;
ax.FontName = 'Helvetica';
ax.YTick = ([-10 0 10 20]);
pbaspect([(x2-x1)*2 14 1]);

ax3 = nexttile(3);
plot(GF_all.xobs,model_lt_vert_mean(1,:),'blue','LineWidth',1.5)
hold on
errorbar(GF.xobs(data_type==2),data_matrix(data_type==2,5),2*data_matrix(data_type==2,8),'.','MarkerSize',12)
xlim([x1,x2])
ylim([-20 20])
ylabel(' ','FontSize',16)
text(-67,-19,'Long-term Vertical','FontName','Helvetica','FontSize',22, ...
    'Rotation',0,'HorizontalAlignment', 'left','VerticalAlignment', 'bottom')
%title('Long-term Vertical')
grid on
ax = gca; % Get current axes
ax.LineWidth = 2;
ax.Box = 'on'; % Remove top and right border
ax.GridLineStyle = '--';
ax.XAxis.TickDirection = 'in';
ax.YAxis.TickDirection = 'out';
ax.TickLength = [0.002 0.005];
ax.FontSize = 16;
ax.FontName = 'Helvetica';
ax.YTick = ([-16 -8 0 8 16]);
ax.XTickLabel = [];
pbaspect([(x2-x1)*2 14 1]);

set(gcf, 'Position', [100, 1100, 2399, 845]);


myfigure = strcat('./figure/', GF.Build_gemo.extractedString, '_', model,'_', kind, '.jpg');
print(gcf, myfigure, '-djpeg', '-r300'); % Save as jpeg
