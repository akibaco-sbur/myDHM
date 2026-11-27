%% DHM 再生データ 自動モード判定・機能検証/計測評価スクリプト
clear; clc; close all;

%% 1. 評価対象ファイルの選択
[f_mat, p_mat] = uigetfile({'*.mat', 'DHM再生データ (recon_*.mat)'}, ...
                           '評価する再生データ (recon_*.mat) を選択してください');
if isequal(f_mat, 0), disp('キャンセルされました。'); return; end
recon_file = fullfile(p_mat, f_mat);

eval_save_dir = uigetdir(p_mat, '【評価】比較画像とレポートの保存先フォルダを選択してください');
if isequal(eval_save_dir, 0), disp('保存先選択がキャンセルされました。'); return; end

loaded = load(recon_file);
recon = loaded.recon;
mode  = recon.mode;
[Ny, Nx] = size(recon.height_nm);
dx = recon.params.dx;

% 外周10%を除外した有効領域マスク
valid = false(Ny, Nx);
my = round(Ny * 0.1); mx = round(Nx * 0.1);
valid(my:end-my, mx:end-mx) = true;

report_lines = {
    '======================================================'
    sprintf('  DHM 評価レポート [Mode %d]', mode)
    '======================================================'
    sprintf('対象ラベル   : %s', recon.label)
    sprintf('再生日時     : %s', recon.timestamp)
    sprintf('評価日時     : %s', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')))
    sprintf('元画像       : %s', recon.source_file)
    sprintf('設定パラメータ: z_prop = %.1f um, phase_sign = %+d, filter_ratio = %.2f', ...
            recon.params.z_prop*1e6, recon.params.phase_sign, recon.params.filter_ratio)
    '------------------------------------------------------'
};

%% 2. モード別評価の実行
switch mode
    case 1  % 【機能検証1】MATLAB模擬反射画像の真値比較（高さ換算・符号・定量精度）
        gt_file = 'ground_truth_nm.mat';
        if ~isfile(gt_file)
            [f_gt, p_gt] = uigetfile('*.mat', '真値ファイル (ground_truth_nm.mat) を選択');
            if isequal(f_gt, 0), return; end
            gt_file = fullfile(p_gt, f_gt);
        end
        gt = load(gt_file);
        h_gt = gt.h_true * 1e9;
        h_rec = recon.height_nm;

        % ガードバンド除去後の再生サイズ [Ny, Nx] に合わせて真値マップを中央トリミング
        [Ny_gt, Nx_gt] = size(h_gt);
        if ~isequal([Ny_gt, Nx_gt], [Ny, Nx])
            if Ny_gt >= Ny && Nx_gt >= Nx
                sy = floor((Ny_gt - Ny) / 2);
                sx = floor((Nx_gt - Nx) / 2);
                h_gt = h_gt(sy+1 : sy+Ny, sx+1 : sx+Nx);
            else
                error('真値データ (%dx%d) が再生データ (%dx%d) より小さいため比較できません。', ...
                      Nx_gt, Ny_gt, Nx, Ny);
            end
        end
        
        [Xm, Ym] = meshgrid((-floor(Nx/2):ceil(Nx/2)-1)*dx, (-floor(Ny/2):ceil(Ny/2)-1)*dx);
        pit_core  = (abs(Xm + 15e-6) < 1.5e-6) & (abs(Ym) < 10e-6);
        bump_core = ((Xm - 15e-6).^2 + Ym.^2) < (6.0e-6)^2;
        bg_core   = valid & (h_gt == 0) & (abs(Xm + 15e-6) > 5e-6) & (((Xm - 15e-6).^2 + Ym.^2) > (11e-6)^2);
        
        pit_meas  = mean(h_rec(pit_core));
        bump_meas = mean(h_rec(bump_core));
        bg_mean   = mean(h_rec(bg_core));
        bg_std    = std(h_rec(bg_core));
        rmse_all  = sqrt(mean((h_rec(valid) - h_gt(valid)).^2));
        
        mode_lines = {
            '[機能検証1: 高さ換算係数・符号・定量精度の検証]'
            sprintf('  1. ピット深さ (真値 -100.0 nm) : %+.2f nm  (誤差: %+.2f nm)', pit_meas, pit_meas + 100)
            sprintf('  2. バンプ高さ (真値  +80.0 nm) : %+.2f nm  (誤差: %+.2f nm)', bump_meas, bump_meas - 80)
            sprintf('  3. 背景平坦度 (真値    0.0 nm) : %+.2f nm ± %.2f nm (SD)', bg_mean, bg_std)
            sprintf('  4. 有効領域全体 RMSE           : %.2f nm', rmse_all)
        };
        
        mid_y = round(Ny / 2);
        fig1 = create_light_fig(['検証1: 真値比較 - ' recon.label], [120, 300, 620, 420]);
        plot(recon.x_um, h_gt(mid_y, :), '--', 'Color', [0.2, 0.2, 0.2], 'LineWidth', 1.6, 'DisplayName', '真値');
        hold on;
        plot(recon.x_um, h_rec(mid_y, :), '-', 'Color', [0.00, 0.40, 0.75], 'LineWidth', 2.0, 'DisplayName', 'DHM復元値');
        grid on; xlim([recon.x_um(1), recon.x_um(end)]);
        legend('Location', 'northwest', 'TextColor', 'k', 'Color', 'w', 'EdgeColor', 'k');
        style_eval_ax(gca, ['真値 vs 復元高さ [' recon.label ']'], 'X [\mum]', '高さ [nm]', '');
        
        err_map = h_rec - h_gt; err_map(~valid) = NaN;
        fig2 = create_light_fig(['検証1: 誤差マップ - ' recon.label], [760, 300, 520, 420]);
        imagesc(recon.x_um, recon.y_um, err_map, 'AlphaData', ~isnan(err_map));
        axis image; colormap(gca, 'turbo');
        style_eval_ax(gca, ['誤差分布マップ [' recon.label ']'], 'X [\mum]', 'Y [\mum]', '誤差 [nm]');
        set(gca, 'Color', [0.92, 0.92, 0.92]);
        
        exportgraphics(fig1, fullfile(eval_save_dir, sprintf('%s_eval_01_profile.png', recon.label)), 'Resolution', 300);
        exportgraphics(fig2, fullfile(eval_save_dir, sprintf('%s_eval_02_errmap.png', recon.label)),  'Resolution', 300);

    case 2  % 【機能検証2】Mode 1未検証機能の成立確認（ノイズ耐性・収差除去・リフォーカス鮮鋭度）
        phi_v = recon.phi_clean(valid);
        
        % (A) 実レンズ収差除去の成立判定: 外周フレーム背景の平坦度(SD)
        bg_frame = valid & ((abs(linspace(-1,1,Nx)) > 0.65) | (abs(linspace(-1,1,Ny)') > 0.65));
        bg_sd_rad = std(recon.phi_clean(bg_frame));
        
        % (B) 実測ノイズ下での位相アンラップ連続性: 隣接画素の最大位相ジャンプ
        max_jump = max([max(abs(diff(recon.phi_clean(valid), 1, 1)), [], 'all'), ...
                        max(abs(diff(recon.phi_clean(valid), 1, 2)), [], 'all')]);
        
        % (C) ASMリフォーカス合焦評価指標: 勾配エネルギー (Tenengrad) と 振幅コントラスト
        [gx, gy] = gradient(recon.phi_clean(valid));
        sharpness_grad = mean(gx.^2 + gy.^2, 'all');
        amp_v = recon.amplitude(valid);
        amp_contrast = std(amp_v) / (mean(amp_v) + eps);
        
        mode_lines = {
            '[機能検証2: Mode 1未検証機能（ノイズ耐性・収差除去・Z伝搬）の成立判定]'
            sprintf('  1. [アンラップ連続性] 最大隣接位相差 : %.3f rad (pi=3.14未満で破綻なし)', max_jump)
            sprintf('  2. [収差平坦化精度]   背景残差 SD    : %.3f rad (お椀型歪み除去の指標)', bg_sd_rad)
            sprintf('  3. [Zリフォーカス指標] 位相勾配鮮鋭度 : %.5f (デフォーカス比較で最大値探索)', sharpness_grad)
            sprintf('  4. [Zリフォーカス指標] 振幅コントラスト: %.4f (合焦面付近で変化)', amp_contrast)
            sprintf('  5. [位相レンジ]       99%%-1%% 変位幅  : %.3f rad', prctile(phi_v, 99) - prctile(phi_v, 1))
        };
        
        [~, max_idx] = max(recon.phi_clean(:) .* valid(:));
        [peak_y, ~] = ind2sub([Ny, Nx], max_idx);
        
        fig1 = create_light_fig(['検証2: 復元位相断面 - ' recon.label], [120, 300, 620, 420]);
        plot(recon.x_um, recon.phi_clean(peak_y, :), 'Color', [0.00, 0.45, 0.35], 'LineWidth', 1.8);
        grid on; xlim([recon.x_um(1), recon.x_um(end)]);
        style_eval_ax(gca, sprintf('位相断面プロファイル (Y=%.1f um, z=%.0f um) [%s]', ...
                      recon.y_um(peak_y), recon.params.z_prop*1e6, recon.label), ...
                      'X [\mum]', '復元位相 [rad]', '');
        
        fig2 = create_light_fig(['検証2: 背景平坦度・位相分布 - ' recon.label], [760, 300, 520, 420]);
        histogram(phi_v, 60, 'FaceColor', [0.20, 0.45, 0.70], 'EdgeColor', 'none');
        grid on;
        style_eval_ax(gca, ['有効領域 位相分布ヒストグラム [' recon.label ']'], '復元位相 [rad]', '画素数', '');
        
        exportgraphics(fig1, fullfile(eval_save_dir, sprintf('%s_eval_01_phase_profile.png', recon.label)), 'Resolution', 300);
        exportgraphics(fig2, fullfile(eval_save_dir, sprintf('%s_eval_02_histogram.png', recon.label)),     'Resolution', 300);

    case 3  % 【本番計測】実機反射型DHMの形状計測レポート
        h_v = recon.height_nm(valid);
        bg_frame = valid & ((abs(linspace(-1,1,Nx)) > 0.65) | (abs(linspace(-1,1,Ny)') > 0.65));
        bg_roughness = std(recon.height_nm(bg_frame));
        
        mode_lines = {
            '[本番計測: 実機反射型DHM 表面形状統計]'
            sprintf('  1. 最大高さ (99%%tile)     : %+.2f nm', prctile(h_v, 99))
            sprintf('  2. 最小深さ (1%%tile)      : %+.2f nm', prctile(h_v, 1))
            sprintf('  3. 段差・高低差 (PV値)    : %.2f nm', prctile(h_v, 99) - prctile(h_v, 1))
            sprintf('  4. 背景面のノイズ床 (Rq)  : %.2f nm', bg_roughness)
        };
        
        mid_y = round(Ny / 2);
        fig1 = create_light_fig(['本番計測: 断面プロファイル - ' recon.label], [120, 300, 620, 420]);
        plot(recon.x_um, recon.height_nm(mid_y, :), 'Color', [0.00, 0.36, 0.65], 'LineWidth', 1.8);
        grid on; xlim([recon.x_um(1), recon.x_um(end)]);
        style_eval_ax(gca, sprintf('計測断面プロファイル (Y=%.1f um) [%s]', recon.y_um(mid_y), recon.label), ...
                      'X [\mum]', '高さ [nm]', '');
        exportgraphics(fig1, fullfile(eval_save_dir, sprintf('%s_eval_01_profile.png', recon.label)), 'Resolution', 300);
end

report_lines = [report_lines; mode_lines; {'======================================================'}];

%% 3. レポート出力 & 保存
report_txt_path = fullfile(eval_save_dir, sprintf('%s_eval_report.txt', recon.label));
fid = fopen(report_txt_path, 'w');
for i = 1:numel(report_lines)
    fprintf('%s\n', report_lines{i});
    if fid ~= -1, fprintf(fid, '%s\n', report_lines{i}); end
end
if fid ~= -1, fclose(fid); end
fprintf('\n[完了] 評価結果を保存しました: %s\n', eval_save_dir);

%% --- 描画補助関数 ---
function fig = create_light_fig(fig_title, pos)
    fig = figure('Name', fig_title, 'NumberTitle', 'off', 'Color', 'w', 'Position', pos);
    try, theme(fig, 'light'); catch, end
end

function style_eval_ax(ax, title_str, x_str, y_str, cb_str)
    set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'FontSize', 10.5);
    title(ax, title_str, 'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'none');
    xlabel(ax, x_str, 'Color', 'k', 'FontSize', 11);
    ylabel(ax, y_str, 'Color', 'k', 'FontSize', 11);
    if ~isempty(cb_str)
        cb = colorbar(ax); cb.Color = 'k'; cb.Label.String = cb_str; cb.Label.Color = 'k';
    end
end