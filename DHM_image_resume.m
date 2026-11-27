%% オフアキシスDHM 像再生スクリプト（デジタル参照波補償・適応フィルタ・IRLS収差除去対応版）
clear; clc; close all;

%% 1. 実行時ダイアログによる Mode 選択
mode_list = {
    'Mode 1 : [機能検証1] MATLAB模擬画像（高さ換算・符号・定量精度の検証）'
    'Mode 2 : [機能検証2] 外部実測画像（ノイズ耐性・収差除去・Z伝搬の検証）'
    'Mode 3 : [本番計測]  実機反射型DHM（CDピット等の実測データ再生）'
};

[mode, tf] = listdlg('PromptString', '実行するモードを選択してください:', ...
                     'SelectionMode', 'single', ...
                     'ListString', mode_list, ...
                     'ListSize', [460, 110], ...
                     'Name', 'DHM 像再生モード選択');
if ~tf
    disp('モード選択がキャンセルされました。処理を終了します。');
    return;
end

%% 2. モード別の初期設定
lambda     = 532e-9;      % 光源波長 [m]
NA         = 0.30;        % 対物レンズ開口数
M          = 10;          % 実効倍率
px_cam     = 3.45e-6;     % カメラ画素ピッチ [m]
poly_order = 2;           % 背景収差フィッティング次数 (1:傾き, 2:お椀型球面収差)

switch mode
    case 1  % 【機能検証1】MATLAB模擬反射画像
        default_label    = 'Verify1_Sim_HeightCheck';
        holo_path        = 'sim_hologram.png';
        ask_bg           = false;
        eval_unit        = 'nm_reflection';
        phase_sign       = -1;
        filter_ratio     = 1.0;
        dc_mask_ratio    = 1.2;
        default_z_um     = 0;
        interactive_sign = false;

    case 2  % 【機能検証2】外部実測画像（未知光学系対応・デジタル参照波補償）
        default_label    = 'Verify2_RealNoise_Refocus';
        holo_path        = '';
        ask_bg           = true;
        eval_unit        = 'rad_algorithm';
        lambda           = 632.8e-9;
        phase_sign       = 1;
        filter_ratio     = 0.42;
        dc_mask_ratio    = 0.25;
        poly_order       = 4;     % 実レンズの高次収差・DCT境界歪みまで除去
        default_z_um     = 0;
        interactive_sign = false;
        
    case 3  % 【本番計測】実機反射型DHM
        default_label    = 'Real_Reflection_Sample01';
        holo_path        = '';
        ask_bg           = true;
        eval_unit        = 'nm_reflection';
        phase_sign       = 1;
        filter_ratio     = 0.85;
        dc_mask_ratio    = 1.1;
        default_z_um     = 0;
        interactive_sign = true;
end

%% 3. 入力画像・背景ホログラム・ラベル・保存先の指定
if isempty(holo_path) || ~isfile(holo_path)
    [f_in, p_in] = uigetfile({'*.png;*.tif;*.tiff;*.bmp;*.jpg', '画像ファイル'}, ...
                             sprintf('[Mode %d] 1/2: 試料ありの干渉縞画像を選択してください', mode));
    if isequal(f_in, 0), disp('キャンセルされました。'); return; end
    holo_path = fullfile(p_in, f_in);
else
    p_in = pwd;
end

bg_path = '';
if ask_bg
    [f_bg, p_bg] = uigetfile({'*.png;*.tif;*.tiff;*.bmp;*.jpg', '画像ファイル'}, ...
                             sprintf('[Mode %d] 2/2: 背景ホログラムを選択（無い場合は「キャンセル」でスキップ）', mode), ...
                             p_in);
    if ~isequal(f_bg, 0)
        bg_path = fullfile(p_bg, f_bg);
    end
end

ans_dlg = inputdlg({'出力データ・画像に付けるラベル名:', ...
                    '数値伝搬距離 z_prop [um] (リフォーカス検証用 / 通常は 0):'}, ...
                   sprintf('Mode %d パラメータ確認', mode), [1 55], ...
                   {default_label, num2str(default_z_um)});
if isempty(ans_dlg), disp('キャンセルされました。'); return; end

sample_label = strtrim(ans_dlg{1});
z_prop       = str2double(ans_dlg{2}) * 1e-6;
if isnan(z_prop), z_prop = 0; end

save_dir = uigetdir(p_in, sprintf('[Mode %d] 結果データと画像の保存先フォルダを選択', mode));
if isequal(save_dir, 0), disp('保存先選択がキャンセルされました。'); return; end

dx = px_cam / M;

%% 4. 画像読み込みと2次元フーリエ変換
img_raw = double(imread(holo_path));
if size(img_raw, 3) == 3, img_raw = mean(img_raw, 3); end
[Ny, Nx] = size(img_raw);

dfx = 1 / (Nx * dx);
dfy = 1 / (Ny * dx);
[Fx, Fy] = meshgrid((-floor(Nx/2) : ceil(Nx/2)-1) * dfx, ...
                    (-floor(Ny/2) : ceil(Ny/2)-1) * dfy);

% 画像中央92%を平坦に保つTukey窓（外周4%ずつを滑らかにゼロへ減衰）
win = tukeywin_custom(Ny, 0.08) * tukeywin_custom(Nx, 0.08)';
img_ac = img_raw - mean(img_raw(:));
F_holo = fftshift(fft2(img_ac .* win));

%% 5. 1次光ピークの自動検出と適応的空間フィルタリング（補正①）
cx0 = floor(Nx/2) + 1;
cy0 = floor(Ny/2) + 1;
[Px, Py] = meshgrid(1:Nx, 1:Ny);
R_px = sqrt((Px - cx0).^2 + (Py - cy0).^2);

if mode == 2
    % 外部画像では光学系パラメータが未知のため、画像サイズ基準でDC除外半径を設定
    dc_mask_radius = round(min(Nx, Ny) * 0.04);
else
    f_cutoff = NA / lambda;
    dc_mask_radius = round((f_cutoff / dfx) * dc_mask_ratio);
end

spec_search = abs(F_holo);
spec_search(R_px <= dc_mask_radius) = 0;
spec_search(:, 1:cx0-1) = 0;
spec_search(cy0:end, cx0) = 0;

[~, max_idx] = max(spec_search(:));
[peak_y, peak_x] = ind2sub([Ny, Nx], max_idx);
d_carrier = hypot(peak_x - cx0, peak_y - cy0);

if mode == 2
    r_filter_px = max(12, round(d_carrier * filter_ratio));
else
    f_cutoff = NA / lambda;
    % dc_notchで0次光を除去しているため、上限を 0.48 -> 0.65 に緩和して高周波を通す
    r_filter_px = min(round((f_cutoff / dfx) * filter_ratio), round(d_carrier * 0.65));
end

% シフト後の1次光中心フィルタ（Raised-Cosine型アポダイゼーションでギブス現象を抑制）
F_shift = circshift(F_holo, [cy0 - peak_y, cx0 - peak_x]);
r_inner = r_filter_px * 0.65;
r_outer = r_filter_px * 1.15;
filter_mask = zeros(Ny, Nx);
filter_mask(R_px <= r_inner) = 1.0;
taper_idx = (R_px > r_inner) & (R_px < r_outer);
filter_mask(taper_idx) = 0.5 * (1 + cos(pi * (R_px(taper_idx) - r_inner) / (r_outer - r_inner)));

% 元の0次光（DC）位置がフィルタ端に掛かる場合のノッチ抑圧
dc_in_shift_x = cx0 + (cx0 - peak_x);
dc_in_shift_y = cy0 + (cy0 - peak_y);
R_dc_shift = sqrt((Px - dc_in_shift_x).^2 + (Py - dc_in_shift_y).^2);
dc_notch = 1 - exp(-(R_dc_shift.^2) / (2 * (0.22 * d_carrier)^2));

F_filtered = F_shift .* filter_mask .* dc_notch;

%% 6. 角スペクトル法（ASM）によるデジタルリフォーカス
if z_prop ~= 0
    kz = sqrt(max(0, (1/lambda)^2 - Fx.^2 - Fy.^2));
    F_filtered = F_filtered .* exp(-1i * 2 * pi * z_prop * kz);
end

U_obj = ifft2(ifftshift(F_filtered));

%% 7. 背景ホログラム減算（指定された場合のみ実行）
if ~isempty(bg_path) && isfile(bg_path)
    bg_raw = double(imread(bg_path));
    if size(bg_raw, 3) == 3, bg_raw = mean(bg_raw, 3); end
    if ~isequal(size(bg_raw), [Ny, Nx])
        error('背景ホログラムの画像サイズ (%dx%d) が本画像 (%dx%d) と一致しません。', ...
              size(bg_raw,2), size(bg_raw,1), Nx, Ny);
    end
    bg_ac = bg_raw - mean(bg_raw(:));
    F_bg = fftshift(fft2(bg_ac .* win));
    F_bg_shift = circshift(F_bg, [cy0 - peak_y, cx0 - peak_x]) .* filter_mask .* dc_notch;
    U_bg = ifft2(ifftshift(F_bg_shift));
    U_obj = U_obj .* conj(U_bg);
    fprintf('[情報] 背景ホログラム減算を適用しました: %s\n', bg_path);
else
    fprintf('[情報] 背景ホログラムなし（単一ショット数値位相補償を実行）\n');
end

%% 7.5 複素空間でのデジタル参照波補償（補正②：アンラップ前の球面波・傾き除去）
[U_obj, phi_comp_ref] = compensate_complex_wavefront(U_obj, poly_order);

% 窓関数で振幅が減衰した外周ガードバンド（5%）を除去し、ゼロ振幅位相ノイズを遮断
margin_y = max(4, round(Ny * 0.05));
margin_x = max(4, round(Nx * 0.05));
U_obj        = U_obj(margin_y+1 : end-margin_y, margin_x+1 : end-margin_x);
phi_comp_ref = phi_comp_ref(margin_y+1 : end-margin_y, margin_x+1 : end-margin_x);
[Ny, Nx]     = size(U_obj);

%% 8. 位相抽出・アンラップ・全画面ロバスト収差除去（補正③）
[phi_wrapped, phi_clean] = reconstruct_phase(U_obj, phase_sign, poly_order);
x_um = (0:Nx-1) * dx * 1e6;
y_um = (0:Ny-1) * dx * 1e6;

if interactive_sign
    h_preview = (lambda / (4 * pi)) * phi_clean * 1e9;
    fig_pre = create_light_figure('符号確認プレビュー (仮符号: +1)', [300, 300, 900, 380]);
    subplot(1, 2, 1);
    imagesc(x_um, y_um, h_preview); axis image; colormap(gca, 'turbo');
    style_light_axes(gca, '仮復元高さマップ [nm]', 'X [\mum]', 'Y [\mum]', '高さ [nm]');
    subplot(1, 2, 2);
    plot(x_um, h_preview(round(Ny/2), :), 'Color', [0.00, 0.36, 0.65], 'LineWidth', 1.6); grid on;
    style_light_axes(gca, '中央断面プロファイル', 'X [\mum]', '高さ [nm]', '');
    
    ans_btn = questdlg('凹凸の向きは正しいですか？（CDピットが凹になっていれば「そのまま」を選択）', ...
                       '位相符号の確定', 'そのまま (+1)', '反転する (-1)', 'そのまま (+1)');
    close(fig_pre);
    if strcmp(ans_btn, '反転する (-1)')
        phase_sign = -1;
        [phi_wrapped, phi_clean] = reconstruct_phase(U_obj, phase_sign, poly_order);
    end
end

%% 9. 物理量換算 & 振幅フラットフィールド正規化（補正④） & データ保存
height_nm = (lambda / (4 * pi)) * phi_clean * 1e9;
if strcmp(eval_unit, 'nm_reflection')
    z_data  = height_nm;
    z_label = '高さ [nm]';
else
    z_data  = phi_clean;
    z_label = '復元位相 [rad]';
end

amp_raw = abs(U_obj);
amp_env = imgaussfilt_simple(amp_raw, max(15, round(min(Ny, Nx) * 0.08)));
amplitude = amp_raw ./ max(amp_env, median(amp_env(:)) * 0.1);

recon.mode        = mode;
recon.eval_unit   = eval_unit;
recon.label       = sample_label;
recon.timestamp   = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
recon.source_file = holo_path;
recon.bg_file     = bg_path;
recon.params      = struct('lambda', lambda, 'NA', NA, 'M', M, 'px_cam', px_cam, ...
                           'dx', dx, 'z_prop', z_prop, 'phase_sign', phase_sign, ...
                           'filter_ratio', filter_ratio, 'r_filter_px', r_filter_px, ...
                           'poly_order', poly_order);
recon.height_nm   = height_nm;
recon.phi_clean   = phi_clean;
recon.phi_comp_ref = phi_comp_ref;
recon.amplitude   = amplitude;
recon.amplitude_raw = amp_raw;
recon.x_um        = x_um;
recon.y_um        = y_um;

out_mat = fullfile(save_dir, sprintf('recon_%s.mat', sample_label));
save(out_mat, 'recon');

%% 10. 個別ウィンドウ出力（ライトモード） & 画像一括保存
fig1 = create_light_figure(['1. FFTスペクトル - ' sample_label], [60, 520, 480, 420]);
imagesc(log(1 + abs(F_holo))); axis image; colormap(gca, 'gray'); hold on;
viscircles_simple(peak_x, peak_y, r_filter_px, [0.85, 0.15, 0.15]);
style_light_axes(gca, ['2D FFTスペクトル [' sample_label ']'], 'X [px]', 'Y [px]', '');

fig2 = create_light_figure(['2. 再生振幅像 - ' sample_label], [560, 520, 480, 420]);
amp_disp_max = quantile_simple(amplitude(:), 0.995);
imagesc(x_um, y_um, amplitude, [0, max(amp_disp_max, 1.2)]); axis image; colormap(gca, 'gray');
style_light_axes(gca, ['再生振幅像 (フラット補正済) [' sample_label ']'], 'X [\mum]', 'Y [\mum]', '正規化振幅');

fig3 = create_light_figure(['3. ラップ位相 - ' sample_label], [1060, 520, 480, 420]);
imagesc(x_um, y_um, phi_wrapped, [-pi, pi]); axis image; colormap(gca, 'parula');
style_light_axes(gca, ['ラップ位相 (-π ～ +π) [' sample_label ']'], 'X [\mum]', 'Y [\mum]', '位相 [rad]');

fig4 = create_light_figure(['4. 復元マップ - ' sample_label], [60, 40, 480, 420]);
imagesc(x_um, y_um, z_data); axis image; colormap(gca, 'turbo');
style_light_axes(gca, ['復元マップ [' sample_label ']'], 'X [\mum]', 'Y [\mum]', z_label);

fig5 = create_light_figure(['5. 断面プロファイル - ' sample_label], [560, 40, 480, 420]);
mid_y = round(Ny / 2);
plot(x_um, z_data(mid_y, :), 'Color', [0.00, 0.36, 0.65], 'LineWidth', 1.8);
grid on; xlim([x_um(1), x_um(end)]);
style_light_axes(gca, sprintf('断面プロファイル (Y=%.1f um) [%s]', y_um(mid_y), sample_label), ...
                 'X [\mum]', z_label, '');

fig6 = create_light_figure(['6. 3D表面形状 - ' sample_label], [1060, 40, 680, 540]);
% 高解像度メッシュ化（最大512×512グリッド）とアンチエイリアス平滑化
step = max(1, round(min(Nx, Ny) / 512));
if step > 1
    z_3d_disp = imgaussfilt_simple(z_data, step * 0.45);
else
    z_3d_disp = z_data;
end
surf(x_um(1:step:end), y_um(1:step:end), z_3d_disp(1:step:end, 1:step:end), ...
     'EdgeColor', 'none', 'FaceColor', 'interp');
shading interp;
colormap(gca, 'turbo');
view(-35, 45);
camlight('headlight');
lighting gouraud;
material dull;
axis tight;
style_light_axes(gca, ['3D表面形状 [' sample_label ']'], 'X [\mum]', 'Y [\mum]', z_label);
zlabel(z_label, 'Color', 'k', 'FontSize', 11, 'FontWeight', 'bold');

exportgraphics(fig1, fullfile(save_dir, sprintf('%s_01_fft.png', sample_label)),     'Resolution', 300);
exportgraphics(fig2, fullfile(save_dir, sprintf('%s_02_amp.png', sample_label)),     'Resolution', 300);
exportgraphics(fig3, fullfile(save_dir, sprintf('%s_03_wrap.png', sample_label)),    'Resolution', 300);
exportgraphics(fig4, fullfile(save_dir, sprintf('%s_04_map.png', sample_label)),     'Resolution', 300);
exportgraphics(fig5, fullfile(save_dir, sprintf('%s_05_profile.png', sample_label)), 'Resolution', 300);
exportgraphics(fig6, fullfile(save_dir, sprintf('%s_06_3d.png', sample_label)),      'Resolution', 600);
fprintf('[完了] Mode %d (%s) の再生データと画像を保存しました:\n  %s\n', mode, sample_label, save_dir);

%% 11. AI診断用定量ログのCSV自動追記
log_csv_path = fullfile(save_dir, 'dhm_ai_diagnostic_log.csv');
diag_summary = export_ai_diagnostic_log(log_csv_path, recon, F_holo, peak_x, peak_y, ...
    r_filter_px, d_carrier, phi_wrapped, z_data, z_label);
fprintf('[AIログ出力] %s\n  自動診断ステータス: %s\n', log_csv_path, diag_summary.auto_flags);

%% --- 補助関数群 ---
function [U_comp, phi_ref] = compensate_complex_wavefront(U_obj, poly_order)
    % 複素振幅の隣接位相差から傾き・2次球面波曲率を直接推定し、unwrap前に相殺する
    [Ny, Nx] = size(U_obj);
    [X, Y] = meshgrid(linspace(-1, 1, Nx), linspace(-1, 1, Ny));
    dx_n = 2 / max(1, Nx - 1);
    dy_n = 2 / max(1, Ny - 1);

    % 隣接画素間の位相差（ラップの影響を受けない位相勾配）
    gx = angle(U_obj(:, 2:end) .* conj(U_obj(:, 1:end-1))) / dx_n;
    gy = angle(U_obj(2:end, :) .* conj(U_obj(1:end-1, :))) / dy_n;
    wx0 = abs(U_obj(:, 2:end) .* U_obj(:, 1:end-1));
    wy0 = abs(U_obj(2:end, :) .* U_obj(1:end-1, :));

    Xm = 0.5 * (X(:, 2:end) + X(:, 1:end-1));
    Ym = Y(:, 2:end);
    Xn = X(2:end, :);
    Yn = 0.5 * (Y(2:end, :) + Y(1:end-1, :));

    % 計算高速化と境界ハン窓減衰の回避のため中央90%領域を間引きサンプリング
    st = max(1, round(min(Ny, Nx) / 256));
    rx = round(Ny*0.05)+1 : st : round(Ny*0.95);
    cx = round(Nx*0.05)+1 : st : round(Nx*0.95)-1;
    ry = round(Ny*0.05)+1 : st : round(Ny*0.95)-1;
    cy = round(Nx*0.05)+1 : st : round(Nx*0.95);

    xm_s = Xm(rx, cx); ym_s = Ym(rx, cx); gx_s = gx(rx, cx); wx_s = wx0(rx, cx);
    xn_s = Xn(ry, cy); yn_s = Yn(ry, cy); gy_s = gy(ry, cy); wy_s = wy0(ry, cy);

    xm_v = xm_s(:); ym_v = ym_s(:); gx_v = gx_s(:); wx_v = wx_s(:) / (median(wx_s(:)) + eps);
    xn_v = xn_s(:); yn_v = yn_s(:); gy_v = gy_s(:); wy_v = wy_s(:) / (median(wy_s(:)) + eps);

    if poly_order == 1
        Mx = [ones(size(xm_v)), zeros(size(xm_v))];
        My = [zeros(size(xn_v)), ones(size(xn_v))];
    else
        Mx = [ones(size(xm_v)), zeros(size(xm_v)), 2*xm_v, ym_v, zeros(size(xm_v))];
        My = [zeros(size(xn_v)), ones(size(xn_v)), zeros(size(xn_v)), xn_v, 2*yn_v];
    end
    M_all = [Mx; My];
    b_all = [gx_v; gy_v];
    w_base = [wx_v; wy_v];

    % エッジの急勾配を除外するため、各方向ごとに勾配が中央値以下の平坦領域に初期重みを限定
    flat_x = double(abs(gx_v - median(gx_v)) <= prctile_simple(abs(gx_v - median(gx_v)), 60));
    flat_y = double(abs(gy_v - median(gy_v)) <= prctile_simple(abs(gy_v - median(gy_v)), 60));
    w_base = w_base .* [flat_x; flat_y];

    % 反復重み付き最小二乗法（IRLS / Tukey Bisquare）で物体エッジの影響を排除
    w_rob = ones(size(b_all));
    for iter = 1:6
        W = sqrt(max(0, w_base .* w_rob));
        theta = (M_all .* W) \ (b_all .* W);
        res = b_all - M_all * theta;
        mad_val = median(abs(res - median(res))) * 1.4826 + 1e-6;
        u = res / (4.685 * mad_val);
        w_rob = (abs(u) < 1) .* ((1 - u.^2).^2);
    end

    if poly_order == 1
        phi_ref = theta(1)*X + theta(2)*Y;
    else
        phi_ref = theta(1)*X + theta(2)*Y + theta(3)*(X.^2) + theta(4)*(X.*Y) + theta(5)*(Y.^2);
    end

    U_comp = U_obj .* exp(-1i * phi_ref);
    % 背景の基準位相（ピストン成分）を 0 rad に合わせ、±pi 境界での折り返しノイズを防止
    piston_phase = angle(mean(U_comp(:)));
    U_comp = U_comp .* exp(-1i * piston_phase);
    
end

function [phi_wrapped, phi_clean] = reconstruct_phase(U_obj, phase_sign, poly_order)
    [Ny, Nx] = size(U_obj);
    phi_wrapped = phase_sign * angle(U_obj);
    phi_unwrapped = unwrap2D_dct(phi_wrapped);
    
    % 全画素対象のIRLS多項式フィッティングで微細な残留うねりを除去
    [X_norm, Y_norm] = meshgrid(linspace(-1, 1, Nx), linspace(-1, 1, Ny));
    st = max(1, round(min(Ny, Nx) / 256));
    Xs = X_norm(1:st:end, 1:st:end);
    Ys = Y_norm(1:st:end, 1:st:end);
    Zs = phi_unwrapped(1:st:end, 1:st:end);

    xm = Xs(:); ym = Ys(:); zm = Zs(:);
    if poly_order == 1
        A_fit = [ones(size(xm)), xm, ym];
        A_all = [ones(numel(X_norm),1), X_norm(:), Y_norm(:)];
    elseif poly_order == 2
        A_fit = [ones(size(xm)), xm, ym, xm.^2, xm.*ym, ym.^2];
        A_all = [ones(numel(X_norm),1), X_norm(:), Y_norm(:), X_norm(:).^2, X_norm(:).*Y_norm(:), Y_norm(:).^2];
    else
        % 4次多項式基底（球面収差・非点収差・コマ収差・DCT境界ベル型歪みを除去）
        A_fit = [ones(size(xm)), xm, ym, ...
                 xm.^2, xm.*ym, ym.^2, ...
                 xm.^3, (xm.^2).*ym, xm.*(ym.^2), ym.^3, ...
                 xm.^4, (xm.^3).*ym, (xm.^2).*(ym.^2), xm.*(ym.^3), ym.^4];
        X1 = X_norm(:); Y1 = Y_norm(:);
        A_all = [ones(size(X1)), X1, Y1, ...
                 X1.^2, X1.*Y1, Y1.^2, ...
                 X1.^3, (X1.^2).*Y1, X1.*(Y1.^2), Y1.^3, ...
                 X1.^4, (X1.^3).*Y1, (X1.^2).*(Y1.^2), X1.*(Y1.^3), Y1.^4];
    end

    w = ones(size(zm));
    for iter = 1:6
        W = sqrt(max(0, w));
        c = (A_fit .* W) \ (zm .* W);
        res = zm - A_fit * c;
        mad_val = median(abs(res - median(res))) * 1.4826 + 1e-6;
        u = res / (4.685 * mad_val);
        w = (abs(u) < 1) .* ((1 - u.^2).^2);
    end

    phi_bg_fit = reshape(A_all * c, Ny, Nx);
    phi_clean  = phi_unwrapped - phi_bg_fit;
    phi_clean  = phi_clean - median(phi_clean(:));
end

function q = quantile_simple(v, p)
    vs = sort(v(:));
    idx = max(1, min(numel(vs), round(p * numel(vs))));
    q = vs(idx);
end

function fig = create_light_figure(fig_title, pos)
    fig = figure('Name', fig_title, 'NumberTitle', 'off', 'Color', 'w', 'Position', pos);
    try theme(fig, 'light'); catch, end
end

function style_light_axes(ax, title_str, x_str, y_str, cb_str)
    set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k', ...
            'GridColor', [0.75, 0.75, 0.75], 'GridAlpha', 0.8, 'FontSize', 10.5, 'LineWidth', 1.0);
    title(ax, title_str, 'Color', 'k', 'FontSize', 11.5, 'FontWeight', 'bold', 'Interpreter', 'none');
    xlabel(ax, x_str, 'Color', 'k', 'FontSize', 11);
    ylabel(ax, y_str, 'Color', 'k', 'FontSize', 11);
    if ~isempty(cb_str)
        cb = colorbar(ax);
        cb.Color = 'k'; cb.Label.String = cb_str; cb.Label.Color = 'k'; cb.Label.FontSize = 10.5;
    end
end

function phi_unwrap = unwrap2D_dct(phi_wrap)
    [M, N] = size(phi_wrap);
    dx = [zeros(M,1), angle(exp(1i * diff(phi_wrap, 1, 2)))];
    dy = [zeros(1,N); angle(exp(1i * diff(phi_wrap, 1, 1)))];
    rho = diff([dx, zeros(M,1)], 1, 2) + diff([dy; zeros(1,N)], 1, 1);
    dct_rho = dct2_custom(rho);
    [kx, ky] = meshgrid(0:N-1, 0:M-1);
    denom = 2 * (cos(pi * kx / N) + cos(pi * ky / M) - 2);
    denom(1, 1) = 1;
    dct_phi = dct_rho ./ denom;
    dct_phi(1, 1) = 0;
    phi_unwrap = idct2_custom(dct_phi);
end

function B = dct2_custom(A), B = dct1_custom(dct1_custom(A).').'; end
function B = idct2_custom(A), B = idct1_custom(idct1_custom(A).').'; end

function B = dct1_custom(A)
    M = size(A, 1);
    v = [A(1:2:M, :); A(M - mod(M,2) : -2 : 2, :)];
    V = fft(v);
    w = exp(-1i * pi * (0:M-1).' / (2 * M)) * sqrt(2 / M);
    w(1) = w(1) / sqrt(2);
    B = real(bsxfun(@times, V, w));
end

function B = idct1_custom(A)
    M = size(A, 1);
    V = zeros(size(A), 'like', 1i);
    V(1, :) = sqrt(M) * A(1, :);
    if M > 1
        k = (1:M-1).';
        w = sqrt(M / 2) * exp(1i * pi * k / (2 * M));
        V(2:M, :) = bsxfun(@times, A(2:M, :) - 1i * A(M:-1:2, :), w);
    end
    v = real(ifft(V));
    B = zeros(M, size(A, 2));
    B(1:2:M, :) = v(1:ceil(M/2), :);
    B(2:2:M, :) = v(M:-1:ceil(M/2)+1, :);
end

function out = imgaussfilt_simple(img, sigma)
    r = ceil(3 * sigma);
    g = exp(-(-r:r).^2 / (2 * sigma^2));
    g = g / sum(g);
    out = conv2(conv2(img, g, 'same'), g', 'same');
end

function viscircles_simple(cx, cy, r, col)
    th = linspace(0, 2*pi, 100);
    plot(cx + r*cos(th), cy + r*sin(th), 'Color', col, 'LineWidth', 1.6);
end

function w = tukeywin_custom(N, alpha)
if nargin < 2, alpha = 0.08; end
if alpha <= 0
    w = ones(N, 1);
    return;
elseif alpha >= 1
    w = hann(N);
    return;
end
n = (0:N-1)';
per = alpha * (N - 1) / 2;
tl = floor(per) + 1;
th = N - tl + 1;
w = ones(N, 1);
w(1:tl) = 0.5 * (1 + cos(pi * (n(1:tl) / per - 1)));
w(th:N) = 0.5 * (1 + cos(pi * ((n(th:N) - (N - 1)) / per + 1)));
end

function diag = export_ai_diagnostic_log(csv_path, recon, F_holo, peak_x, peak_y, ...
                                         r_filter_px, d_carrier, phi_wrapped, z_data, z_label)
    [Ny, Nx] = size(z_data);
    cx0 = floor(Nx/2) + 1;
    cy0 = floor(Ny/2) + 1;

    % 1. FFT・フィルタ関連指標
    filter_to_carrier = r_filter_px / max(d_carrier, 1);
    dc_power = abs(F_holo(cy0, cx0));
    peak_power = abs(F_holo(peak_y, peak_x));
    carrier_snr_db = 20 * log10(max(peak_power, eps) / max(median(abs(F_holo(:))), eps));

    % 2. 収差補償・ダイナミックレンジ指標
    pre_comp_span_rad = max(recon.phi_comp_ref(:)) - min(recon.phi_comp_ref(:));
    final_phi_span_rad = max(recon.phi_clean(:)) - min(recon.phi_clean(:));
    
    % 試料本体（外れ値）を除外した背景画素のみで、中央領域 vs 外周領域の残留うねりを評価
    phi_med = median(recon.phi_clean(:));
    phi_mad = median(abs(recon.phi_clean(:) - phi_med)) * 1.4826 + 1e-9;
    is_bg_phi = abs(recon.phi_clean - phi_med) < (3.0 * phi_mad);

    [Xn, Yn] = meshgrid(linspace(-1, 1, Nx), linspace(-1, 1, Ny));
    r_norm = hypot(Xn, Yn);
    center_bg = recon.phi_clean(is_bg_phi & (r_norm < 0.45));
    outer_bg  = recon.phi_clean(is_bg_phi & (r_norm > 0.70 & r_norm < 0.95));
    if ~isempty(center_bg) && ~isempty(outer_bg)
        bowl_residual_rad = abs(median(center_bg) - median(outer_bg));
    else
        bowl_residual_rad = 0;
    end
    
    % お椀型残留収差のチェック（中央領域 vs 外周領域の中央値差）
    [Xn, Yn] = meshgrid(linspace(-1, 1, Nx), linspace(-1, 1, Ny));
    r_norm = hypot(Xn, Yn);
    center_med = median(recon.phi_clean(r_norm < 0.3));
    outer_med  = median(recon.phi_clean(r_norm > 0.75 & r_norm < 0.95));
    bowl_residual_rad = abs(center_med - outer_med);

    % 3. アンラップ整合性（位相残差・渦の数）と境界スパイク指標
    w_dx = angle(exp(1i * diff(phi_wrapped, 1, 2)));
    w_dy = angle(exp(1i * diff(phi_wrapped, 1, 1)));
    curl_loop = w_dx(1:end-1, :) + w_dy(:, 2:end) - w_dx(2:end, :) - w_dy(:, 1:end-1);
    residue_count = nnz(abs(curl_loop) > pi);

    inner_data = z_data(3:end-2, 3:end-2);
    border_vals = [z_data(1,:), z_data(end,:), z_data(:,1)', z_data(:,end)'];
    inner_iqr = max(prctile_simple(inner_data(:), 95) - prctile_simple(inner_data(:), 5), 1e-6);
    border_spike_ratio = max(abs(border_vals - median(inner_data(:)))) / inner_iqr;

    % 4. 背景平坦度（リンギング・十字ノイズ・照明ムラ指標）
    % 外れ値（段差本体）を除いた背景領域の標準偏差
    med_z = median(inner_data(:));
    mad_z = median(abs(inner_data(:) - med_z)) * 1.4826 + 1e-9;
    bg_mask = abs(inner_data - med_z) < (3.0 * mad_z);
    bg_rms = std(inner_data(bg_mask));
    signal_span = prctile_simple(inner_data(:), 99.5) - prctile_simple(inner_data(:), 0.5);
    signal_to_bg_rms = signal_span / max(bg_rms, 1e-9);
    amp_cv = std(recon.amplitude(:)) / max(mean(recon.amplitude(:)), 1e-9);

    % 5. AI向け自動診断フラグ生成
    flags = {};
    if filter_to_carrier > 0.68
        flags{end+1} = 'WARN_DC_LEAK_RISK';
    end
    if bowl_residual_rad > 0.35
        flags{end+1} = 'WARN_BOWL_ABERRATION_RESIDUAL';
    end
    if residue_count > 10
        flags{end+1} = sprintf('WARN_UNWRAP_RESIDUES(%d)', residue_count);
    end
    if border_spike_ratio > 1.5
        flags{end+1} = 'WARN_BORDER_SPIKE';
    end
    if signal_to_bg_rms < 5.0
        flags{end+1} = 'WARN_LOW_SB_OR_RINGING';
    end
    if isempty(flags)
        auto_flags = 'OK';
    else
        auto_flags = strjoin(flags, ';');
    end

    % 1行テーブル作成とCSV追記
    T_row = table( ...
        string(recon.timestamp), ...
        recon.mode, ...
        string(recon.label), ...
        string(z_label), ...
        Nx, Ny, ...
        round(d_carrier, 2), ...
        r_filter_px, ...
        round(filter_to_carrier, 3), ...
        round(carrier_snr_db, 1), ...
        round(pre_comp_span_rad, 3), ...
        round(final_phi_span_rad, 3), ...
        round(bowl_residual_rad, 4), ...
        residue_count, ...
        round(border_spike_ratio, 3), ...
        round(signal_span, 3), ...
        round(bg_rms, 4), ...
        round(signal_to_bg_rms, 2), ...
        round(amp_cv, 3), ...
        string(auto_flags), ...
        'VariableNames', { ...
            'timestamp', 'mode', 'label', 'unit', 'Nx', 'Ny', ...
            'd_carrier_px', 'r_filter_px', 'filter_to_carrier', 'carrier_snr_db', ...
            'pre_comp_span_rad', 'final_phi_span_rad', 'bowl_residual_rad', ...
            'unwrap_residues', 'border_spike_ratio', 'signal_span', ...
            'bg_rms', 'signal_to_bg_rms', 'amp_cv', 'auto_flags' ...
        });

    if isfile(csv_path)
        writetable(T_row, csv_path, 'WriteMode', 'append', 'WriteVariableNames', false);
    else
        writetable(T_row, csv_path, 'WriteMode', 'overwrite', 'WriteVariableNames', true);
    end

    diag.auto_flags = auto_flags;
end

function q = prctile_simple(v, p)
    vs = sort(v(:));
    idx = max(1, min(numel(vs), round((p / 100) * numel(vs))));
    q = vs(idx);
end