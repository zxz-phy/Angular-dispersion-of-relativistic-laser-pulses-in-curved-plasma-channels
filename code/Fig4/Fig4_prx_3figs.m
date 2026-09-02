clc; clear; close all;

script_dir = fileparts(mfilename('fullpath'));
figs_root = fileparts(fileparts(script_dir));
data_dir = fullfile(figs_root, 'data', 'Fig4');
output_dir = fullfile(figs_root, 'Figs', 'Fig4');
if ~isfolder(output_dir), mkdir(output_dir); end

red = load(fullfile(data_dir, ...
    'Prism_Thick_dispersion_relation_L3mm_S2000R20000_1e3lines.mat'));
blue = load(fullfile(data_dir, ...
    'Prism_Thick_dispersion_relation_L3mm_S6000R20000_1e3lines.mat'));
lam_red = red.lam_prism; theta_red = -red.angle_prism;
lam_blue = blue.lam_prism; theta_blue = -blue.angle_prism;

% Preserve the smoothing and finite-difference processing of the selected source.
theta_blue_s = smooth(theta_blue, 100);
lam_blue_s = smooth(lam_blue, 100);
dtheta_blue = diff(theta_blue_s) ./ (lam_blue_s(2) - lam_blue_s(1));
dtheta_blue(1:67) = smooth(dtheta_blue(1:67), 5);
dtheta_blue(1:67) = smooth(dtheta_blue(1:67), 105);
dtheta_red = diff(theta_red) ./ (lam_red(2) - lam_red(1));
dlam_blue = (lam_blue_s(1:end-1) + lam_blue_s(2:end)) / 2;
dlam_red = (lam_red(1:end-1) + lam_red(2:end)) / 2;

fs = 18;
fg = figure(4);
set(fg, 'Color', 'w', 'Units', 'pixels', 'Position', [300, 80, 655, 860]);
layout = tiledlayout(fg, 3, 1, 'TileSpacing', 'none', 'Padding', 'compact');

ax_a = nexttile(layout); plot(ax_a, lam_red, theta_red, 'r-', 'LineWidth', 3);
xlim(ax_a, [0.8, 10]); ylim(ax_a, [1.5, 3.5]); yticks(ax_a, [2, 2.5, 3]);
ylabel(ax_a, '$\theta~(^{\circ})$', 'Interpreter', 'latex');
legend(ax_a, '$S=2$ mm, $R=20$ mm', 'Interpreter', 'latex', ...
    'Location', 'southeast', 'Box', 'off');
add_panel_label(ax_a, '(a)', fs);

ax_b = nexttile(layout); plot(ax_b, lam_blue, theta_blue, 'b-', 'LineWidth', 3);
xlim(ax_b, [0.8, 10]); ylim(ax_b, [18.5, 20.5]); yticks(ax_b, [19, 19.5, 20]);
ylabel(ax_b, '$\theta~(^{\circ})$', 'Interpreter', 'latex');
legend(ax_b, '$S=6$ mm, $R=20$ mm', 'Interpreter', 'latex', ...
    'Location', 'northeast', 'Box', 'off');
add_panel_label(ax_b, '(b)', fs);

ax_c = nexttile(layout); hold(ax_c, 'on');
plot(ax_c, dlam_red, dtheta_red, 'r-', 'LineWidth', 3);
plot(ax_c, dlam_blue, dtheta_blue, 'b-', 'LineWidth', 3);
xlim(ax_c, [0.8, 10]); ylim(ax_c, [-0.6, 0.6]); yticks(ax_c, [-0.4, 0, 0.4]);
xlabel(ax_c, '$\lambda~(\mu\mathrm{m})$', 'Interpreter', 'latex');
ylabel(ax_c, '$d\theta/d\lambda$', 'Interpreter', 'latex');
add_panel_label(ax_c, '(c)', fs);

set([ax_a, ax_b, ax_c], 'FontName', 'Times New Roman', 'FontSize', fs, ...
    'Box', 'on', 'TickDir', 'out', 'XTick', [0.8, 5, 10]);
set([ax_a, ax_b], 'XTickLabel', []);

exportgraphics(fg, fullfile(output_dir, 'Fig4.png'), 'Resolution', 600);
print(fg, fullfile(output_dir, 'Fig4.eps'), '-depsc', '-painters');
exportgraphics(fg, fullfile(output_dir, 'Fig4.pdf'), 'ContentType', 'vector');
fprintf('Saved Fig. 4 to %s\n', output_dir);

function add_panel_label(ax, label, fs)
text(ax, 0.015, 0.12, label, 'Units', 'normalized', ...
    'FontName', 'Times New Roman', 'FontSize', fs + 4, ...
    'FontWeight', 'bold', 'Clipping', 'off');
end
