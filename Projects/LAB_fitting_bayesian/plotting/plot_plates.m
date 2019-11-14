function plot_plates(posterior_L, q_method)


figure('position', [300, 500, 1000, 400], 'color', 'w'); 
a1 = axes('position', [0.1, 0.2, 0.25, 0.7]);
plot(posterior_L.zPlate, posterior_L.p_zPlate); xlabel('zPlate (km)');
set(a1, 'ycolor', 'none', 'box', 'off')
a2 = axes('position', [0.4, 0.2, 0.25, 0.7]);
plot(posterior_L.Tp, posterior_L.p_Tp); xlabel('Tp (C)')
set(a2, 'ycolor', 'none', 'box', 'off')
a3 = axes('position', [0.7, 0.2, 0.25, 0.7]);
plot(posterior_L.zLAB, posterior_L.p_zLAB); xlabel('zLAB (km)')
set(a3, 'ycolor', 'none', 'box', 'off')
axes('position', [0, 0, 1, 1], 'visible', 'off')
text(0.39, 0.95, 'Probabilities of zPlate, Tp, and zLAB', 'fontsize', 14)


%%
Tp_bins = posterior_L.Tp + [0, diff(posterior_L.Tp)/2];
Tp_bins(1) = 0; Tp_bins(end) = max([posterior_L.MC_Tp + 1, Tp_bins]);
[n_Tp, b_Tp] = histc(posterior_L.MC_Tp, Tp_bins);

zPlate_bins = posterior_L.zPlate + [0, diff(posterior_L.zPlate)/2];
zPlate_bins(1) = 0; zPlate_bins(end) = ...
    max([posterior_L.MC_zPlate + 1, zPlate_bins]);
[n_zPlate, b_zPlate] = histc(posterior_L.MC_zPlate, zPlate_bins);

max_n_bin = max(n_Tp) + max(n_zPlate);
min_n_bin = min(n_Tp) + min(n_zPlate);

figure('position', [200, 100, 1000, 800], 'color', 'w');
a1 = axes('position', [0.1, 0.1, 0.25, 0.7]);  hold on
a2 = axes('position', [0.4, 0.1, 0.25, 0.7]);  hold on
a3 = axes('position', [0.7, 0.1, 0.25, 0.7]);  hold on
max_zplate = 0;

Box = load(posterior_L.Files.SV_Box);
VBR = load(posterior_L.Files.VBR_Box);

for i_Tp = 1:size(Box.Box, 1)
    for i_zPlate = 1:size(Box.Box, 2)
        b = Box.Box(i_Tp, i_zPlate);
        z = b.run_info.Z_km;
        t = b.Frames(end).T;
        % Plot Vs and Q at the frequency closest to T = 80 s.
        [~, i_f] = min(abs(1 ./ VBR.VBR(i_Tp, i_zPlate).in.SV.f - 80));
        vs = VBR.VBR(i_Tp, i_zPlate).out.anelastic.(q_method).V(:, i_f) ./1e3;
        q = VBR.VBR(i_Tp, i_zPlate).out.anelastic.(q_method).Q(:, i_f);
        n_this_combo = sum(b_Tp == i_Tp & b_zPlate == i_zPlate);
        transp_val = min(1, (n_this_combo - min_n_bin) ...
            / (max_n_bin - min_n_bin) * 4);
        [~, i_LAB] = min(abs(posterior_L.zLAB_grid(i_Tp, i_zPlate) - z));
        if transp_val < 0.05
            continue
        end
        
        axes(a1);
        p = patch([t; flipud(t)+10], [z; flipud(z)], 'red', ...
             'linestyle', 'none');
        p.FaceAlpha = transp_val;
        l = patch(t(i_LAB) + 20 * [-1, 1, 1, -1],  ...
            posterior_L.zLAB_grid(i_Tp, i_zPlate) + 1.5 * [-1, -1, 1, 1], ...
            'k', 'linestyle', 'none');
        l.FaceAlpha = transp_val;
        axis ij
        ylabel('Depth (km)')
        xlabel('Temperature (C)')
        
        
        
        axes(a2);
        p = patch([vs; flipud(vs)+0.01], [z; flipud(z)], 'red', ...
             'linestyle', 'none');
        p.FaceAlpha = transp_val;
        l = patch(vs(i_LAB) + 0.02 * [-1, 1, 1, -1],  ...
            posterior_L.zLAB_grid(i_Tp, i_zPlate) + 1.5 * [-1, -1, 1, 1], ...
            'k', 'linestyle', 'none');
        l.FaceAlpha = transp_val;
        axis ij
        ylabel('Depth (km)')
        xlabel('Vs (km/s)')
%         
%         
%         
        axes(a3);
        p = patch([q; flipud(q)+20], [z; flipud(z)], 'red', ...
             'linestyle', 'none');
        p.FaceAlpha = transp_val;
        l = patch(q(i_LAB) + 20 * [-1, 1, 1, -1],  ...
            posterior_L.zLAB_grid(i_Tp, i_zPlate) + 1.5 * [-1, -1, 1, 1], ...
            'k', 'linestyle', 'none');
        l.FaceAlpha = transp_val;
        axis ij
        ylabel('Depth (km)')
        xlabel('Q')
        xlim([0, 2000])
        
        
        max_zplate = max(max_zplate, posterior_L.zPlate(i_zPlate));
        
    end
end

set(a1, 'ylim', [0, min(max(z), max_zplate + 100)], 'xaxislocation', 'top')
set(a2, 'ylim', [0, min(max(z), max_zplate + 100)], 'xaxislocation', 'top')
set(a3, 'ylim', [0, min(max(z), max_zplate + 100)], 'xaxislocation', 'top')


axes('position', [0, 0, 1, 1], 'visible', 'off')
text(0.25, 0.965, 'Suite of Likely Geotherms with corresponding Vs and Q profiles', ...
    'fontsize', 14)
zLAB = posterior_L.zLAB(posterior_L.p_zLAB == max(posterior_L.p_zLAB));
Tp = posterior_L.Tp(posterior_L.p_Tp == max(posterior_L.p_Tp));
text(0.32, 0.93, ...
    sprintf('given likeliest zLAB = %.0f km, Tp = %.0f %cC', zLAB, Tp, 176), ...
    'fontsize', 14)
text(0.4, 0.895, ...
    sprintf('using %s', strrep(q_method, '_', ' ')), 'fontsize', 14)


end