function DLC_huntplot_xp(ctl_mat, xpro_mat, type, exp, hs)
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
    ctld1 = mean(ctl_mat(:,1:3),2);
    xprod1 = mean(xpro_mat(:,1:3),2);

    ctld2 = mean(ctl_mat(:,4:5),2);
    xprod2 = mean(xpro_mat(:,4:5),2);

    %DAY 1
    figure
    bar(0.25,mean(ctld1),0.4,'FaceColor', [0.7 0.7 0.7])
    hold on,
    bar(0.75,mean(xprod1),0.4,'FaceColor', [0.7 0 1])
    plot(0.15*ones(1,length(ctld1)), ctld1,'kx','LineWidth',4);
    plot(0.65*ones(1,length(xprod1)), xprod1,'kx','LineWidth',4);
    errorbar(0.25,mean(ctld1), std(ctld1)/sqrt(numel(ctld1)), 'k', 'LineWidth',4)
    errorbar(0.75,mean(xprod1), std(xprod1)/sqrt(numel(xprod1)), 'k', 'LineWidth',4)

    if exp == 1
        ylabel('Time to Capture (s)')
    elseif exp == 2
        ylabel('Latency to Attack (s)')
    elseif exp == 3
        ylabel('Pursuit Duration (s)')
    end

    xlim([0 1])
    xticks([0.25 0.75])
    xticklabels({'CTL'; 'XPRO'});
    set(findall(gcf,'-property','FontSize'),'FontSize',18,'FontWeight','bold')
    set(gca,'box','off')
    title('Day 1')
    hAx=gca;
    hAx.LineWidth=3;
    set(gca,'TickDir','out')

    %PROBE
    figure
    bar(0.25,mean(ctld2),0.4,'FaceColor', [0.7 0.7 0.7])
    hold on,
    bar(0.75,mean(xprod2),0.4,'FaceColor', [0.7 0 1])
    plot(0.15*ones(1,length(ctld2)), ctld2,'kx','LineWidth',4);
    plot(0.65*ones(1,length(xprod2)), xprod2,'kx','LineWidth',4);
    errorbar(0.25,mean(ctld2), std(ctld2)/sqrt(numel(ctld2)), 'k', 'LineWidth',4)
    errorbar(0.75,mean(xprod2), std(xprod2)/sqrt(numel(xprod2)), 'k', 'LineWidth',4)

    if exp == 1
        ylabel('Time to Capture (s)')
    elseif exp == 2
        ylabel('Latency to Attack (s)')
    elseif exp == 3
        ylabel('Pursuit Duration (s)')
    end

    xlim([0 1])
    xticks([0.25 0.75])
    xticklabels({'CTL'; 'XPRO'});
    set(findall(gcf,'-property','FontSize'),'FontSize',18,'FontWeight','bold')
    set(gca,'box','off')
    title('PROBE')
    hAx=gca;
    hAx.LineWidth=3;
    set(gca,'TickDir','out')


end


    
end