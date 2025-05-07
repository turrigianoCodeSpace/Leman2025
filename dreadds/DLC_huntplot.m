function DLC_huntplot(ctl_mat, dr_mat, type, exp, hs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Type inputs are either: 1 (pool by day) or 2 (pool by hunting session)
%Exp inputs are either: 1 (time to capture), 2 (latency to attack) or 3 (pursuit
%duration)
%HS inputs are either: 0 (no sessions plotted), 1,2,3 for individual
%hunting, or 4 (all sessions)
%session
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compare conditions across day

if type == 1
    ctld1 = mean(ctl_mat(:,1:3),'all');
    drd1 = mean(dr_mat(:,1:3),'all');


    figure
    bar(0.25,mean(ctld1),0.4,'FaceColor', [0.7 0.7 0.7])
    hold on,
    bar(0.75,mean(drd1),0.4,'FaceColor', '#FF6347')
    plot(0.15*ones(1,length(ctld1)), ctl_mat(1:3),'kx','LineWidth',4);
    plot(0.65*ones(1,length(drd1)), dr_mat(1:3),'kx','LineWidth',4);
    errorbar(0.25,mean(ctld1), std(ctld1)/sqrt(numel(ctld1)), 'k', 'LineWidth',4)
    errorbar(0.75,mean(drd1), std(drd1)/sqrt(numel(drd1)), 'k', 'LineWidth',4)

    if exp == 1
        ylabel('Time to Capture (s)')
    elseif exp == 2
        ylabel('Latency to Attack (s)')
    elseif exp == 3
        ylabel('Pursuit Duration (s)')
    end

    xlim([0 1])
    xticks([0.25 0.75])
    xticklabels({'Cre^-'; 'Cre^+'});
    set(findall(gcf,'-property','FontSize'),'FontSize',18,'FontWeight','bold')
    set(gca,'box','off')
    title('Day 1')
    hAx=gca;
    hAx.LineWidth=3;
    set(gca,'TickDir','out')


end

%% Compare Conditions by hunting session
if type == 2 && hs <= 4

    if hs >0 && hs <=4
        if hs == 3
            i = hs;

            ctlhs = ctl_mat(:,i);
            drhs = dr_mat(:,i);

            figure
            bar(0.25,mean(ctlhs),0.4,'FaceColor', [0.7 0.7 0.7])
            hold on,
            bar(0.75,mean(drhs),0.4,'FaceColor', '#FF6347')
            plot(0.15*ones(1,length(ctlhs)), ctlhs,'kx','LineWidth',4);
            plot(0.65*ones(1,length(drhs)), drhs,'kx','LineWidth',4);
            errorbar(0.25,mean(ctlhs), std(ctlhs)/sqrt(numel(ctlhs)), 'k', 'LineWidth',4)
            errorbar(0.75,mean(drhs), std(drhs)/sqrt(numel(drhs)), 'k', 'LineWidth',4)

            if exp == 1
                ylabel('Time to Capture (s)')
            elseif exp == 2
                ylabel('Latency to Attack (s)')
            elseif exp == 3
                ylabel('Pursuit Duration (s)')
            end

            xlim([0 1])
            xticks([0.25 0.75])
            xticklabels({'Cre^-'; 'Cre^+'});
            set(findall(gcf,'-property','FontSize'),'FontSize',18,'FontWeight','bold')
            set(gca,'box','off')
            title(['Hunting Session ' num2str(i) ''])
            hAx=gca;
            hAx.LineWidth=3;
            set(gca,'TickDir','out')

        elseif hs == 4
            hs2plot = 1:3;
            for i = 1:max(hs2plot)
                ctlhs = ctl_mat(:,i);
                drhs = dr_mat(:,i);
                figure
                bar(0.25,mean(ctlhs),0.4,'FaceColor', [0.7 0.7 0.7])
                hold on,
                bar(0.75,mean(drhs),0.4,'FaceColor', '#FF6347')
                plot(0.15*ones(1,length(ctlhs)), ctlhs,'kx','LineWidth',4);
                plot(0.65*ones(1,length(drhs)), drhs,'kx','LineWidth',4);
                errorbar(0.25,mean(ctlhs), std(ctlhs)/sqrt(numel(ctlhs)), 'k', 'LineWidth',4)
                errorbar(0.75,mean(drhs), std(drhs)/sqrt(numel(drhs)), 'k', 'LineWidth',4)

                if exp == 1
                    ylabel('Time to Capture (s)')
                elseif exp == 2
                    ylabel('Latency to Attack (s)')
                elseif exp == 3
                    ylabel('Pursuit Duration (s)')
                end

                xlim([0 1])
                xticks([0.25 0.75])
                xticklabels({'Cre^-'; 'Cre^+'});
                set(findall(gcf,'-property','FontSize'),'FontSize',18,'FontWeight','bold')
                set(gca,'box','off')
                title(['Hunting Session ' num2str(i) ''])
                hAx=gca;
                hAx.LineWidth=3;
                set(gca,'TickDir','out')

            end

        end

        if exp == 1 && hs == 4
            %% Behavior Curve

            if size(ctl_mat,2) ~= size(dr_mat,2)
                disp('Conditions are not the same size!')
                keyboard
            end
            sess2plot = size(ctl_mat,2);

            figure;
            plot(1:sess2plot, mean(dr_mat(:,1:sess2plot)),'Color', '#FF6347', 'LineWidth', 3)
            hold on
            plot(1:sess2plot, mean(ctl_mat(:,1:sess2plot)),'Color', [0.7 0.7 0.7], 'LineWidth', 3)
            errorbar(mean(dr_mat(:,1:sess2plot)),std(dr_mat(:,1:sess2plot))./sqrt(numel(dr_mat(:,1:sess2plot)/sess2plot)) ,'k', 'LineWidth', 2.5,'LineStyle','none')
            errorbar(mean(ctl_mat(:,1:sess2plot)),std(ctl_mat(:,1:sess2plot))./sqrt(numel(ctl_mat(:,1:sess2plot)/sess2plot)) ,'k', 'LineWidth', 2.5,'LineStyle','none')
            
            xlabel('Hunting Session')
            ylabel('Time to Capture (s)')
            xlim([0.5 sess2plot+0.5])
            xticks([1 2 3 4])
            xticklabels({'1','2','3'})
            legend('Cre^+', 'Cre^-')
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