function [fitresult, gof] = createFit1(x, y, z)

%% Fit: '球面拟合'.
[xData, yData, zData] = prepareSurfaceData( x, y, z );

% 设置求解初始
ft = fittype( '-sqrt(r^2-x^2-y^2)', 'independent', {'x', 'y'}, 'dependent', 'z' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = 300;

% 拟合数据
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );

% 绘图
figure( 'Name', '曲面拟合' );
h = plot( fitresult, [xData, yData], zData );
legend( h, '曲面拟合', '(x,y,z)', 'Location', 'NorthEast', 'Interpreter', 'none' );
% 标签
xlabel( 'x', 'Interpreter', 'none' );
ylabel( 'y', 'Interpreter', 'none' );
zlabel( 'z', 'Interpreter', 'none' );
grid on;


