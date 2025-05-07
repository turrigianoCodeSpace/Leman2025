function DLC_huntplot_mf(m_mat, f_mat, type, exp, hs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Type inputs are either: 1 (pool by day) or 2 (pool by hunting session)
%Exp inputs are either: 1 (time to capture), 2 (latency to attack) or 3 (pursuit
%duration)
%HS inputs are either: 0 (no sessions plotted), 1,2,3 for individual
%hunting, or 4 (all sessions)
%session
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compare Conditions by hunting session
if type == 2 && hs <= 4
    if exp == 1 && hs == 4
        %% Behavior Curve

        if size(m_mat,2) ~= size(f_mat,2)
            disp('Conditions are not the same size!')
            keyboard
        end
        sess2plot = size(m_mat,2);

        figure;
        plot(1:sess2plot, mean(m_mat(:,1:sess2plot)),'Color', [0.07,0.62,1.00], 'LineWidth', 3)
        hold on
        plot(1:sess2plot, mean(f_mat(:,1:sess2plot)),'Color', [1.00,0.07,0.65], 'LineWidth', 3)
        errorbar(mean(m_mat(:,1:sess2plot)),std(m_mat(:,1:sess2plot))./sqrt(numel(m_mat(:,1:sess2plot)/sess2plot)) ,'k', 'LineWidth', 2.5,'LineStyle','none')
        errorbar(mean(f_mat(:,1:sess2plot)),std(f_mat(:,1:sess2plot))./sqrt(numel(f_mat(:,1:sess2plot)/sess2plot)) ,'k', 'LineWidth', 2.5,'LineStyle','none')

        xlabel('Hunting Session')
        ylabel('Time to Capture (s)')
        xlim([0.5 sess2plot+0.5])
        xticks([1 2 3 4])
        xticklabels({'1','2','3'})
        m_max = max(m_mat(:,1:sess2plot));
        f_max = max(f_mat(:,1:sess2plot));
        max_mf = max([m_max,f_max]);

        ylim([0 80])
        legend('MALES (n=11)','FEMALES (n=9)')
        title('Average Time to Capture on First Day of Hunting')
        legend boxoff
        set(gca,'box','off')
        set(findall(gcf,'-property','FontSize'),'FontSize',18)
        set(gca,'TickDir','out')
        hAx=gca;
        hAx.LineWidth=3;
        set(gca,'FontSize',18,'Fontweight','bold','LineWidth',6,'box','off','tickdir','out')

    end
end
end
