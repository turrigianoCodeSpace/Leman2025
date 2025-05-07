function DLC_huntplot_1c(anim_mat, type, exp, hs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Type inputs are either: 1 (pool by day) or 2 (pool by hunting session)
%Exp inputs are either: 1 (time to capture), 2 (latency to attack) or 3 (pursuit
%duration)
%HS inputs are either: 0 (no sessions plotted), 1,2,3 for individual
%hunting, or 4 (all sessions)
%session
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(anim_mat,2) ~= size(anim_mat,2)
    disp('Conditions are not the same size!')
    keyboard
end
sess2plot = size(anim_mat,2);

        if hs == 4 && exp > 1
            %% Bar Plots
            hs1 = anim_mat(:,1);
            hs2 = anim_mat(:,2);
            hs3 = anim_mat(:,3);

            sem_hs1 = std(hs1)/sqrt(numel(hs1));
            sem_hs2 = std(hs2)/sqrt(numel(hs2));
            sem_hs3 = std(hs3)/sqrt(numel(hs3));

            figure;
            bar(0.5, mean(hs1),'k','FaceAlpha',0.3,'LineWidth',4,'Barwidth',0.4,'FaceColor', 'k')
            hold on,
            bar(1.25, mean(hs2),'k','FaceAlpha',0.3,'LineWidth',4,'Barwidth',0.4,'FaceColor', 'k')
            bar(2, mean(hs3),'k','FaceAlpha',0.3,'LineWidth',4,'Barwidth',0.4,'FaceColor', 'k')
            plot(0.35*ones(1,length(hs1)), hs1,'ko','LineWidth',1.5);
            plot(1.15*ones(1,length(hs2)), hs2,'ko','LineWidth',1.5);
            plot(1.95*ones(1,length(hs3)), hs3,'ko','LineWidth',1.5);
            errorbar(0.5,mean(hs1), sem_hs1, 'k', 'LineWidth',2)
            errorbar(1.25,mean(hs2), sem_hs2, 'k', 'LineWidth',2)
            errorbar(2,mean(hs3), sem_hs3, 'k', 'LineWidth',2)
            xlabel('Hunting Session')
            xlim([0 2.25])
            xticks([0.5 1.25 2])
            xticklabels({'1'; '2';'3'})
            set(findall(gcf,'-property','FontSize'),'FontSize',18)
            %title('HS1 vs HS3')
            set(gca,'FontSize',18,'Fontweight','bold','LineWidth',6,'box','off','tickdir','out')


    if exp == 2
        ylabel('Latency to Attack (s)')
        ylim([0 40])
    elseif exp == 3
        ylabel('Pursuit Duration (s)')
        ylim([0 150])
    end


        end

        if exp == 1 && hs == 4
            %% Behavior Curve

            if size(anim_mat,2) ~= size(anim_mat,2)
                disp('Conditions are not the same size!')
                keyboard
            end
            sess2plot = size(anim_mat,2);

            figure;
            plot(1:sess2plot, mean(anim_mat(:,1:sess2plot)),'Color', 'm', 'LineWidth', 3)
            hold on
            errorbar(mean(anim_mat(:,1:sess2plot)),std(anim_mat(:,1:sess2plot))./sqrt(numel(anim_mat(:,1:sess2plot)/sess2plot)) ,'k', 'LineWidth', 2.5,'LineStyle','none')
            
            xlabel('Hunting Session')
            ylabel('Time to Capture (s)')
            xlim([0.5 sess2plot+0.5])
            xticks([1 2 3])
            xticklabels({'1','2','3'})
            ylim([0 70])
            set(gca,'box','off')
            set(findall(gcf,'-property','FontSize'),'FontSize',18)
            set(gca,'TickDir','out')
            hAx=gca;
            hAx.LineWidth=3;
            set(gca,'FontSize',18,'Fontweight','bold','LineWidth',6,'box','off','tickdir','out')
        end
    
    
end
