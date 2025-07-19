%% 读入显示器特性化参数
clear;

%% 载入数据
RGB = load('D:\颜色实验用程序\Projector_files\CalibrationData\15-Jan-2025\T1\RGBcal.txt');
XYZcal = load('D:\颜色实验用程序\Projector_files\CalibrationData\15-Jan-2025\T1\XYZcal.txt');
xyz = XYZcal;

%% 黑点与三原色参考点
hd = xyz(1,:);            % 黑点
r_xyz = xyz(26,:) - hd;   % R 基准块（xyz第26行）
g_xyz = xyz(40,:) - hd;   % G 基准块（xyz第40行）
b_xyz = xyz(53,:) - hd;   % B 基准块（xyz第53行）

mjuzheng = [r_xyz; g_xyz; b_xyz]';   % 正向矩阵
mni = inv(mjuzheng);                % 逆矩阵

%% 提取原色块数据，转换为归一化RGB值
Rz = zeros(13,1); Gz = zeros(13,1); Bz = zeros(13,1);

for i = 15:27
    delta = xyz(i,:) - hd;
    Rz(i-14) = mni(1,:) * delta';
end

for i = 28:40
    delta = xyz(i,:) - hd;
    Gz(i-27) = mni(2,:) * delta';
end

for i = 41:53
    delta = xyz(i,:) - hd;
    Bz(i-40) = mni(3,:) * delta';
end

%% GOG拟合每个通道
br = RGB(15:27,1)/255;
bg = RGB(28:40,2)/255;
bb = RGB(41:53,3)/255;


fr = fittype('(kgr*t + (1-kgr))^gamr', ...
    'independent','t', ...
    'coefficients',{'kgr','gamr'});

fg = fittype('(kgg*t + (1-kgg))^gamg', ...
    'independent','t', ...
    'coefficients',{'kgg','gamg'});

fb = fittype('(kgb*t + (1-kgb))^gamb', ...
    'independent','t', ...
    'coefficients',{'kgb','gamb'});

% 设置约束范围（kg在0~1之间，gamma在0.5~3之间为常见区间）
opts = fitoptions('Method','NonlinearLeastSquares',...
    'Lower', [0, 0.5], ...
    'Upper', [1, 3], ...
    'StartPoint', [0.5, 2]);

cfunr = fit(br, Rz, fr, opts);
cfung = fit(bg, Gz, fg, opts);
cfunb = fit(bb, Bz, fb, opts);


% 可视化拟合结果
xi = linspace(0,1,100);
figure; plot(br, Rz, 'ro', xi, cfunr(xi), 'b-'); title('R通道GOG拟合');
figure; plot(bg, Gz, 'go', xi, cfung(xi), 'b-'); title('G通道GOG拟合');
figure; plot(bb, Bz, 'bo', xi, cfunb(xi), 'b-'); title('B通道GOG拟合');

%% 应用 GOG 拟合将 RGB → XYZ
cfitR = coeffvalues(cfunr);
cfitG = coeffvalues(cfung);
cfitB = coeffvalues(cfunb);

delt = [];
for i = 1:84
    Rin = RGB(i,1)/255;
    Gin = RGB(i,2)/255;
    Bin = RGB(i,3)/255;

    Rzz = (cfitR(1)*Rin + (1 - cfitR(1)))^cfitR(2);
    Gzz = (cfitG(1)*Gin + (1 - cfitG(1)))^cfitG(2);
    Bzz = (cfitB(1)*Bin + (1 - cfitB(1)))^cfitB(2);

    % 保证暗通道为 0
    if RGB(i,1)==0, Rzz = 0; end
    if RGB(i,2)==0, Gzz = 0; end
    if RGB(i,3)==0, Bzz = 0; end

    Xan = r_xyz(1)*Rzz + g_xyz(1)*Gzz + b_xyz(1)*Bzz + hd(1);
    Yan = r_xyz(2)*Rzz + g_xyz(2)*Gzz + b_xyz(2)*Bzz + hd(2);
    Zan = r_xyz(3)*Rzz + g_xyz(3)*Gzz + b_xyz(3)*Bzz + hd(3);

    xyz1 = [Xan, Yan, Zan];
    xyz2 = XYZcal(i,:);
    
    uv1 = xyz2uv(xyz1);
    uv2 = xyz2uv(xyz2);
    
    delt(i) = deltaE(uv1, uv2, 12);
end

% 显示 ΔE 色差图
figure; plot(delt, 'k.-'); title('ΔE 色差'); xlabel('Patch Index'); ylabel('ΔE');
fprintf('平均 ΔE: %.4f\n最大 ΔE: %.4f\n', mean(delt), max(delt));



