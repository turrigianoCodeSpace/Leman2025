function DLC_huntplot(m_mat, f_mat, type, exp, hs)
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
    md1 = mean(m_mat(:,1:3),'all');
    fd1 = mean(f_mat(:,1:3),'all');


    figure
    bar(0.25,mean(md1),0.4,'FaceColor', [0.07,0.62,1.00])
    hold on,
    bar(0.75,mean(fd1),0.4,'FaceColor', '#FF6347')
    plot(0.15*ones(1,length(md1)), m_mat(1:3),'kx','LineWidth',4);
    plot(0.65*ones(1,length(fd1)), f_mat(1:3),'kx','LineWidth',4);
    errorbar(0.25,mean(md1), std(md1)/sqrt(numel(md1)), 'k', 'LineWidth',4)
    errorbar(0.75,mean(fd1), std(fd1)/sqrt(numel(fd1)), 'k', 'LineWidth',4)

    if exp == 1
        ylabel('Time to Capture (s)')
    elseif exp == 2
        ylabel('Latency to Attack (s)')
    elseif exp == 3
        ylabel('Pursuit Duration (s)')
    end

    xlim([0 1])
    xticks([0.25 0.75])
    xticklabels({'Males'; 'Females'});
    set(findall(gcf,'-property','FontSize'),'FontSize',18,'FontWeight','bold')
    set(gca,'box','off')
    title('Day 1')
    hAx=gca;
    hAx.LineWidth=3;
    set(gca,'TickDir','out')

    %Calculate Stats
    wrs_d1 = ranksum(reshape(m_mat(:,1:3),[1 numel(m_mat(:,1:3))]),reshape(f_mat(:,1:3),[1 numel(f_mat(:,1:3))]))

end

%% Compare Conditions by hunting session
if type == 2 && hs <= 4

    %% Compare individual sessions
    if hs >0 && hs <=4
        if hs < 4
            i = hs;

            mhs = m_mat(:,i);
            fhs = f_mat(:,i);

            figure
            bar(0.25,mean(mhs),0.4,'FaceColor', [0.07,0.62,1.00])
            hold on,
            bar(0.75,mean(fhs),0.4,'FaceColor', [1.00,0.07,0.65])
            plot(0.15*ones(1,length(mhs)), mhs,'kx','LineWidth',4);
            plot(0.65*ones(1,length(fhs)), fhs,'kx','LineWidth',4);
            errorbar(0.25,mean(mhs), std(mhs)/sqrt(numel(mhs)), 'k', 'LineWidth',4)
            errorbar(0.75,mean(fhs), std(fhs)/sqrt(numel(fhs)), 'k', 'LineWidth',4)

            if exp == 1
                ylabel('Time to Capture (s)')
            elseif exp == 2
                ylabel('Latency to Attack (s)')
            elseif exp == 3
                ylabel('Pursuit Duration (s)')
            end

            xlim([0 1])
            xticks([0.25 0.75])
            xticklabels({'Male'; 'Females'});
            set(findall(gcf,'-property','FontSize'),'FontSize',18,'FontWeight','bold')
            set(gca,'box','off')
            title(['Hunting Session ' num2str(i) ''])
            hAx=gca;
            hAx.LineWidth=3;
            set(gca,'TickDir','out')

        elseif hs == 4
            hs2plot = 1:3;
            for i = 1:max(hs2plot)
                mhs = m_mat(:,i);
                fhs = f_mat(:,i);
                figure
                bar(0.25,mean(mhs),0.4,'FaceColor', [0.07,0.62,1.00])
                hold on,
                bar(0.75,mean(fhs),0.4,'FaceColor', [1.00,0.07,0.65])
                plot(0.15*ones(1,length(mhs)), mhs,'kx','LineWidth',4);
                plot(0.65*ones(1,length(fhs)), fhs,'kx','LineWidth',4);
                errorbar(0.25,mean(mhs), std(mhs)/sqrt(numel(mhs)), 'k', 'LineWidth',4)
                errorbar(0.75,mean(fhs), std(fhs)/sqrt(numel(fhs)), 'k', 'LineWidth',4)

                if exp == 1
                    ylabel('Time to Capture (s)')
                elseif exp == 2
                    ylabel('Latency to Attack (s)')
                elseif exp == 3
                    ylabel('Pursuit Duration (s)')
                end

                xlim([0 1])
                xticks([0.25 0.75])
                xticklabels({'Male'; 'Females'});
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

            if size(m_mat,2) ~= size(f_mat,2)
                disp('Conditions are not the same size!')
                keyboard
            end
            sess2plot = size(m_mat,2);

            figure;
            plot(1:sess2plot, mean(f_mat(:,1:sess2plot)),'Color', [1.00,0.07,0.65], 'LineWidth', 3)
            hold on
            plot(1:sess2plot, mean(m_mat(:,1:sess2plot)),'Color', [0.07,0.62,1.00], 'LineWidth', 3)
            errorbar(mean(f_mat(:,1:sess2plot)),std(f_mat(:,1:sess2plot))./sqrt(numel(f_mat(:,1:sess2plot)/sess2plot)) ,'k', 'LineWidth', 2.5,'LineStyle','none')
            errorbar(mean(m_mat(:,1:sess2plot)),std(m_mat(:,1:sess2plot))./sqrt(numel(m_mat(:,1:sess2plot)/sess2plot)) ,'k', 'LineWidth', 2.5,'LineStyle','none')
            
            xlabel('Hunting Session')
            ylabel('Time to Capture (s)')
            xlim([0.5 sess2plot+0.5])
            xticks([1 2 3 4])
            xticklabels({'1','2','3'})
            m_max = max(m_mat(:,1:sess2plot));
            f_max = max(f_mat(:,1:sess2plot));
            max_mf = max([m_max,f_max]);

            ylim([0 80])
            legend('FEMALES', 'MALES')
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
    
    if exp == 1
        disp('Stats for Time to Capture')
    elseif exp == 2
        disp('Stats for Latency to Attack')
    elseif exp == 3
        disp('Stats for Pursuit Duration')
    end
    
    disp('Comparing between Conditions')
    wrs_hs1 = ranksum(m_mat(:,1),f_mat(:,1))
    wrs_hs2 = ranksum(m_mat(:,2),f_mat(:,2))
    wrs_hs3 = ranksum(m_mat(:,3),f_mat(:,3))

    disp('Comparing Session 1 to Session 3')
    wsr_m1v3 = signrank(m_mat(:,1),m_mat(:,3))
    wsr_f1v3 = signrank(f_mat(:,1),f_mat(:,3))
    
end
