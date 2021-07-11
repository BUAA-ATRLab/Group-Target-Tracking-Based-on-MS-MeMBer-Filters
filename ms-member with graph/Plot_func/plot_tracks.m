function fHand = plot_tracks(model, truth)

fHand = figure;
subplot(2,1,1) 

%% First plot the birth densities 

for k =  1:truth.K
        if ~isempty(truth.X{k})
            htaracks = plot(truth.X{k}(1,:), truth.X{k}(2,:), '.k','LineWidth',3, 'MarkerSize', 14);              % plot target tracks
        end
        hold on
end
    axis([model.obs.xrange model.obs.yrange])
    ylabel('Y coordinate [m]')
    xlabel('X coordinate [m]')
    grid on; box on; set(gca,'GridLineStyle','--');
    title('True target tracks')
    if ~isempty(hsensors)
        legend([hbirth hsensors htaracks], {'Birth pdf conf. 90%','Sensor positions', 'Target tracks'})
    else
        legend([hbirth htaracks], {'Birth pdf conf. 90%', 'Target tracks'})
    end
    set(gca,'FontSize',25,'fontWeight','bold');
    %axis equal
subplot(2,1,2)   
    plot(1:truth.K, truth.N,'-k','LineWidth',3, 'MarkerSize', 8)
    axis([0 truth.K 1 9])
    ylabel('True cardinality')
    xlabel('Time sample [s]')
    grid on; box on; set(gca,'GridLineStyle','--');
    set(gca,'FontSize',25,'fontWeight','bold');
set(findall(fHand,'type','text'),'fontSize',25,'fontWeight','bold')


end

