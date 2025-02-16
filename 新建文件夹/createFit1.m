function [fitresult, gof] = createFit1(x, y, z)

%% Fit: '�������'.
[xData, yData, zData] = prepareSurfaceData( x, y, z );

% ��������ʼ
ft = fittype( '-sqrt(r^2-x^2-y^2)', 'independent', {'x', 'y'}, 'dependent', 'z' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = 300;

% �������
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );

% ��ͼ
figure( 'Name', '�������' );
h = plot( fitresult, [xData, yData], zData );
legend( h, '�������', '(x,y,z)', 'Location', 'NorthEast', 'Interpreter', 'none' );
% ��ǩ
xlabel( 'x', 'Interpreter', 'none' );
ylabel( 'y', 'Interpreter', 'none' );
zlabel( 'z', 'Interpreter', 'none' );
grid on;


