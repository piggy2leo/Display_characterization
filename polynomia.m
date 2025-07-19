%% 初始化
clear;
basePath = 'D:\颜色实验用程序\Projector_files\CalibrationData\15-Jan-2025\T1\';

RGB = load(fullfile(basePath, 'RGBcal.txt'));    % 84行 RGB
xyz = load(fullfile(basePath, 'XYZcal.txt'));    % 84行 XYZ
XYZcal = xyz;                                     % ground truth

%% 读入黑白点及三原色块
hd = xyz(1,:);     % 黑点（暗场）
bd = xyz(14,:);    % 白点
red_xyz   = xyz(15:27, :);
green_xyz = xyz(28:40, :);
blue_xyz  = xyz(41:53, :);

%% 构建 XYZ→RGB 的逆矩阵 mni
% 使用 R26、G40、B53（在 xyz 中分别为第26,40,53行）
r_xyz = xyz(26,:) - hd;
g_xyz = xyz(40,:) - hd;
b_xyz = xyz(53,:) - hd;

matrix = [r_xyz; g_xyz; b_xyz]';   % 正向矩阵：RGB → XYZ
mni = inv(matrix);                % 逆矩阵：XYZ → RGB

%% 计算归一化 RGB 值（每个通道分开处理）
Rz = zeros(13,1); Gz = zeros(13,1); Bz = zeros(13,1);

for i = 15:27
    delta_xyz = xyz(i,:) - hd;
    Rz(i-14) = mni(1,:) * delta_xyz';
end

for i = 28:40
    delta_xyz = xyz(i,:) - hd;
    Gz(i-27) = mni(2,:) * delta_xyz';
end

for i = 41:53
    delta_xyz = xyz(i,:) - hd;
    Bz(i-40) = mni(3,:) * delta_xyz';
end

%% 多项式拟合每个通道（8阶）
br = RGB(15:27,1)/255;
bg = RGB(28:40,2)/255;
bb = RGB(41:53,3)/255;

pr = polyfit(br, Rz, 8);
pg = polyfit(bg, Gz, 8);
pb = polyfit(bb, Bz, 8);

%% 可视化拟合（可选）
xi = linspace(0,1,100);
figure, plot(br, Rz, 'r*', xi, polyval(pr, xi), 'b-'); title('R通道拟合');
figure, plot(bg, Gz, 'g*', xi, polyval(pg, xi), 'b-'); title('G通道拟合');
figure, plot(bb, Bz, 'b*', xi, polyval(pb, xi), 'b-'); title('B通道拟合');

%% 用拟合模型将 RGBcal 拟合 → 得到 XYZ 重建值
Rin = RGB(1:84,1)/255;
Gin = RGB(1:84,2)/255;
Bin = RGB(1:84,3)/255;

Rzz = polyval(pr, Rin);
Gzz = polyval(pg, Gin);
Bzz = polyval(pb, Bin);

% 保证暗通道值为 0
Rzz(RGB(:,1)==0) = 0;
Gzz(RGB(:,2)==0) = 0;
Bzz(RGB(:,3)==0) = 0;

% 用 mjuzheng 正向矩阵还原 XYZ
Xan = r_xyz(1)*Rzz + g_xyz(1)*Gzz + b_xyz(1)*Bzz + hd(1);
Yan = r_xyz(2)*Rzz + g_xyz(2)*Gzz + b_xyz(2)*Bzz + hd(2);
Zan = r_xyz(3)*Rzz + g_xyz(3)*Gzz + b_xyz(3)*Bzz + hd(3);

XYZ_reconstructed = [Xan, Yan, Zan];

%% ΔE（uv）色差计算
uv1 = xyz2uv(XYZ_reconstructed);
uv2 = xyz2uv(XYZcal(1:84,:));
delt = deltaE(uv1, uv2, 12);

figure, plot(delt, 'k.-'); title('ΔE 色差'); xlabel('patch'); ylabel('ΔE');

