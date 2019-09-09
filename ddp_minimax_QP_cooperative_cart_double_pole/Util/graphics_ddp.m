function  stop = graphics_ddp(figures,x,u,cost,L,Vx,Vxx,fx,fxx,fu,fuu,trace,init)
stop = 0;

if figures == 0
    return;
end

n  = size(x,1);
N  = size(x,2);
nL = size(L,2);
m  = size(u,1);

cost  = sum(cost,1);

% === first figure
if figures ~= 0  && ( mod(trace(end,1)-1,figures) == 0 || init == 2 )
    
    fig1 = findobj(0,'name','iLQG');
    if  isempty(fig1)
        fig1 = figure();
        set(fig1,'NumberTitle','off','Name','iLQG','KeyPressFcn',@Kpress,'user',0,'toolbar','none');
        fprintf('Type ESC in the graphics window to terminate early.\n')
    end
    
    if size(trace,1) == 1
        set(fig1,'user',0);
    end
    
    set(0,'currentfigure',fig1);
    clf(fig1);
    
    ax1   = subplot(2,2,1);
    set(ax1,'XAxisL','top','YAxisL','right','xlim',[1 N],'xtick',[])
    line(1:N,cost,'linewidth',4,'color',.5*[1 1 1]);
    ax2 = axes('Position',get(ax1,'Position'));
    plot((1:N),x','linewidth',2);
    set(ax2,'xlim',[1 N],'Ygrid','on','YMinorGrid','off','color','none');
    set(ax1,'Position',get(ax2,'Position'));
    double_title(ax1,ax2,'state','running cost')
    
    axL = subplot(2,2,3);
    CO = get(axL,'colororder');
    set(axL,'nextplot','replacechildren','colororder',CO(1:min(n,7),:))
    Lp = reshape(permute(L,[2 1 3]), [nL*m N-1])';
    plot(axL,1:N-1,Lp,'linewidth',1,'color',0.7*[1 1 1]);
    ylim  = get(axL,'Ylim');
    ylim  = [-1 1]*max(abs(ylim));
    set(axL,'XAxisL','top','YAxisL','right','xlim',[1 N],'xtick',[],'ylim',ylim)
    axu = axes('Position',get(axL,'Position'));
    plot(axu,(1:N-1),u(:,1:N-1)','linewidth',2);
    ylim  = get(axu,'Ylim');
    ylim  = [-1 1]*max(abs(ylim));
    set(axu,'xlim',[1 N],'Ygrid','on','YMinorGrid','off','color','none','ylim',ylim);
    set(axL,'Position',get(axu,'Position'));
    double_title(axu,axL,'controls','gains')
    xlabel 'timesteps'
    
    T        = trace(:,1);
    mT       = max(T);
    
    ax1      = subplot(2,2,2);
    set(ax1,'XAxisL','top','YAxisL','right','xlim',[1 mT+eps],'xtick',[])
    hV = line(T,trace(:,7),'linewidth',4,'color',.5*[1 1 1]);
    ax2 = axes('Position',get(ax1,'Position'));
    hT = semilogy(T,max(0, trace(:,2:5)),'.-','linewidth',2,'markersize',10);
    set(ax2,'xlim',[1 mT+eps],'Ygrid','on','YMinorGrid','off','color','none');
    set(ax1,'Position',get(ax2,'Position'));
    double_title(ax1,ax2,'convergence trace','total cost')
    
    subplot(2,2,4);
    plot(T,trace(:,6),'.-','linewidth',2);
    title 'actual/expected reduction ratio'
    set(gca,'xlim',[0 mT+1],'ylim',[0 2],'Ygrid','on');
    xlabel 'iterations'
    
    set(findobj(fig1,'-property','FontSize'),'FontSize',8)
    stop = get(fig1,'user');
end

if figures < 0  &&  (mod(abs(trace(end,1))-1,figures) == 0 || init == 2) && ~isempty(Vx)
    
    fig2 = findobj(0,'name','iLQG - derivatives');
    if  isempty(fig2)
        fig2 = figure();
        set(fig2,'NumberTitle','off','Name','iLQG - derivatives','KeyPressFcn',@Kpress,'user',0);
    end
    
    if size(trace,1) == 1
        set(fig2,'user',0);
    end
    
    set(0,'currentfigure',fig2);
    clf(fig2);
    
    subplot(2,3,1);
    plot(1:N,Vx','linewidth',2);
    set(gca,'xlim',[1 N]);
    title 'V_x'
    grid on;
    
    subplot(2,3,4);
    z = reshape(Vxx,nL^2,N)';
    zd = (1:nL+1:nL^2);
    plot(1:N,z(:,setdiff(1:nL^2,zd)),'color',.5*[1 1 1]);
    hold on;
    plot(1:N,z(:,zd),'linewidth',2);
    hold off
    grid on;
    set(gca,'xlim',[1 N]);
    title 'V_{xx}'
    xlabel 'timesteps'
    
    subplot(2,3,2);
    Nfx     = size(fx,3);
    z = reshape(fx,nL^2,Nfx)';
    zd = (1:n+1:n^2);
    plot(1:Nfx,z(:,setdiff(1:n^2,zd)),'color',.5*[1 1 1]);
    hold on;
    plot(1:Nfx,z,'linewidth',2);
    set(gca,'xlim',[1 Nfx+eps]);
    hold off
    grid on;
    title 'f_{x}'
    
    if numel(fxx) > 0
        fxx = fxx(:,:,:,1:N-1);
        subplot(2,3,5);
        z  = reshape(fxx,[numel(fxx)/(N-1) N-1])';
        plot(1:N-1,z);
        title 'f_{xx}'
        grid on;
        set(gca,'xlim',[1 N-1+eps]);
    end
    
    subplot(2,3,3);
    Nfu     = size(fu,3);
    z = reshape(fu,nL*m,Nfu)';
    plot(1:Nfu,z','linewidth',2);
    set(gca,'xlim',[1 Nfu]);
    title 'f_u'
    grid on;
    
    if numel(fuu) > 0
        subplot(2,3,6);
        fuu = fuu(:,:,:,1:N-1);
        z  = reshape(fuu,[numel(fuu)/(N-1) N-1])';
        plot(1:N-1,z);
        title 'f_{uu}'
        grid on;
        set(gca,'xlim',[1 N-1+eps]);
    end
    
    set(findobj(fig2,'-property','FontSize'),'FontSize',8)
    stop = stop + get(fig2,'user');
end

if init == 1
    figure(fig1);
elseif init == 2
    strings  = {'V','\lambda','\alpha','\partial_uV','\Delta{V}'};
    legend([hV; hT],strings,'Location','Best');
end

drawnow;

function Kpress(src,evnt)
if strcmp(evnt.Key,'escape')
    set(src,'user',1)
end

function double_title(ax1, ax2, title1, title2)

t1 = title(ax1, title1);
set(t1,'units','normalized')
pos1 = get(t1,'position');
t2 = title(ax2, title2);
set(t2,'units','normalized')
pos2 = get(t2,'position');
[pos1(2),pos2(2)] = deal(min(pos1(2),pos2(2)));
pos1(1)  = 0.05;
set(t1,'pos',pos1,'HorizontalAlignment','left')
pos2(1)  = 1-0.05;
set(t2,'pos',pos2,'HorizontalAlignment','right')