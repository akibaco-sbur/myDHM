%% 反射型オフアキシスDHM 模擬干渉縞生成スクリプト
clear; clc; close all;

% 1. 今回の光学系パラメータ
lambda = 532e-9;          % 波長: 532 nm
NA     = 0.30;            % 対物レンズ開口数 (渋谷光学 LDJ10X)
M      = 10;              % 倍率 (チューブレンズ f=200mm 組合せ時)
px_cam = 3.45e-6;         % カメラ画素ピッチ仮定値 (3.45 um)
dx     = px_cam / M;      % 試料面上の実効ピクセルサイズ (0.345 um/px)
N      = 512;             % 画像サイズ (512 x 512 px)

[X, Y] = meshgrid((-N/2 : N/2-1) * dx, (-N/2 : N/2-1) * dx);

% 2. 正解サンプル形状 h_true (単位: m)
% 例: CDピットを模した「深さ -100 nm の凹み」と、比較用の「高さ +80 nm の凸部」
h_true = zeros(N, N);
pit_mask  = (abs(X + 15e-6) < 3e-6) & (abs(Y) < 12e-6);   % 凹ピット (-100 nm)
bump_mask = (X - 15e-6).^2 + Y.^2 < (8e-6)^2;             % 凸バンプ (+80 nm)
h_true(pit_mask)  = -100e-9;
h_true(bump_mask) = +80e-9;

% 3. 反射型における物体光の複素振幅 (往復光路差 2*h)
% 凸(h > 0)は光路が短くなり位相が進む(+)、凹(h < 0)は位相が遅れる(-)
phi_obj = (4 * pi / lambda) * h_true;
U_obj_ideal = exp(1i * phi_obj);

% 4. 対物レンズの開口数(NA)による空間周波数制限 (横分解能 R ≈ 1.08 um の再現)
df = 1 / (N * dx);
[Fx, Fy] = meshgrid((-N/2 : N/2-1) * df, (-N/2 : N/2-1) * df);
F_rad = sqrt(Fx.^2 + Fy.^2);
f_cutoff = NA / lambda;                    % カットオフ周波数
CTF = double(F_rad <= f_cutoff);           % コヒーレント伝達関数
U_obj = ifft2(ifftshift(fftshift(fft2(U_obj_ideal)) .* CTF));

% 5. 対物レンズ由来の球面収差（お椀型の位相歪み）を意図的に付加
phi_aberration = 1.5e9 * (X.^2 + Y.^2);
U_obj = U_obj .* exp(1i * phi_aberration);

% 6. オフアキシス参照光 (キャリア周波数: 1周期あたり約 5 px になるよう傾ける)
fx_ref = 1 / (3.2 * dx);
fy_ref = 1 / (3.2 * dx);
U_ref = 1.0 * exp(1i * 2 * pi * (fx_ref * X + fy_ref * Y));

% 7. 干渉縞強度の計算 + コヒーレンス低下(V=0.8) + センサノイズ + 8-bit量子化
V = 0.8; % コヒーレンス長(1~2mm)によるビジビリティ低下
I_raw = abs(U_obj).^2 + abs(U_ref).^2 + 2 * V * real(U_obj .* conj(U_ref));
I_noisy = I_raw + 0.05 * randn(N, N);      % ガウシアンノイズ付加

% 0〜255 の 8-bit 画像へスケーリングして保存
I_norm = (I_noisy - min(I_noisy(:))) / (max(I_noisy(:)) - min(I_noisy(:)));
holo_img = uint8(I_norm * 255);

imwrite(holo_img, 'sim_hologram.png');
save('ground_truth_nm.mat', 'h_true', 'dx', 'lambda');

figure;
subplot(1,2,1); imshow(holo_img); title('模擬オフアキシス反射干渉縞');
subplot(1,2,2); imagesc(h_true*1e9); axis image; colorbar; title('正解高さ [nm]');