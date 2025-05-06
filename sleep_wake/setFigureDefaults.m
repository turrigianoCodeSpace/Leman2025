function setFigureDefaults

if ismac
    set(0,'defaultAxesLineWidth',2);
    set(0,'defaultAxesFontSize',16);
    set(0,'defaultAxesXColor','k');
    set(0,'defaultAxesYColor','k');
    set(0,'defaultFigureColor','w');
    set(0,'defaultFigureUnit','normalized');
    set(0,'defaultFigurePosition',[.1 .2 .6 .7]);
elseif ispc
    set(0,'defaultAxesLineWidth',2);
    set(0,'defaultAxesFontSize',14);
    set(0,'defaultAxesXColor','k');
    set(0,'defaultAxesYColor','k');
    set(0,'defaultFigureColor','w');
    set(0,'defaultFigureUnit','normalized');
    set(0,'defaultFigurePosition',[.2 .3 .5 .5]);
end

end