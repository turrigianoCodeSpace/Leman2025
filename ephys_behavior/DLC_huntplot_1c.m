function DLC_huntplot_1c(anim_mat, type, exp, hs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Type inputs are either: 1 (pool by day) or 2 (pool by hunting session)
%Exp inputs are either: 1 (time to capture), 2 (latency to attack) or 3 (pursuit
%duration)
%HS inputs are either: 0 (no sessions plotted), 1,2,3 for individual
%hunting, or 4 (all sessions)
%session
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if exp == 1 && hs == 4
            %% Behavior Curve

            if size(anim_mat,2) ~= size(anim_mat,2)
                disp('Conditions are not the same size!')
                keyboard
            end
            sess2plot = size(anim_mat,2);

            figure;
            plot(1:sess2plot, mean(anim_mat(:,1:sess2plot)),'Color', [1.00,0.00,1.00], 'LineWidth', 3)
            hold on
            errorbar(mean(anim_mat(:,1:sess2plot)),std(anim_mat(:,1:sess2plot))./sqrt(numel(anim_mat(:,1:sess2plot)/sess2plot)) ,'k', 'LineWidth', 2.5,'LineStyle','none')
            
            xlabel('Hunting Session')
            ylabel('Time to Capture (s)')
            xlim([0.5 sess2plot+0.5])
            xticks([1 2 3 4])
            xticklabels({'1','2','3'})
            
            title('Average Time to Capture on First Day of Hunting')
            set(gca,'box','off')
            set(findall(gcf,'-property','FontSize'),'FontSize',18)
            set(gca,'TickDir','out')
            hAx=gca;
            hAx.LineWidth=3;
            set(gca,'FontSize',18,'Fontweight','bold','LineWidth',6,'box','off','tickdir','out')
        end
    
    
    if exp == 1
        disp('Stats for Time to Capture')
    elseif exp == 2
        disp('Stats for Latency to Attack')
    elseif exp == 3
        disp('Stats for Pursuit Duration')
    end
    

    disp('Comparing Session 1 to Session 3')
    wsr_1v3 = signrank(anim_mat(:,1),anim_mat(:,3))
    
end