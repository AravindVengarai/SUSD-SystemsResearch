% Credit to Sean Maxon from Georgia Tech Systems Research Laboratory %

function graph2Dsusd(obj_fun, xx, axis_limits)
    % Generating Function Values for Contour %
    x_min = axis_limits(1);
    x_max = axis_limits(2);
    y_min = axis_limits(3);
    y_max = axis_limits(4);
    
    step = (x_max - x_min) / 200;

    z = obj_fun;

    xp = x_min:step:x_max;
    [X,Y] = meshgrid(xp,xp);
    Z = z([X(:),Y(:)]');
    Z = reshape(Z, size(X));
    
    % Figure window setup %
    title_text = 'Testing';
    disp(['Running ' title_text])
    f = figure('name',title_text,'pos',[200 100 700 600]);
    num_swarms = size(xx,3);   
    if num_swarms <= 5
        colors = zeros(10, 3);
        colors(1,:) = [1 0 0];
        colors(2,:) = [0 1 0];
        colors(3,:) = [0 0 1];
        colors(4,:) = [1 1 0];
        colors(5,:) = [0 1 1];

    else
        colors = jet(num_swarms);
    end

    % Number of iterations %
    MAX_ITER = size(xx,4);
    axlist = {1,MAX_ITER};
    panlist = {1,MAX_ITER};
    
    % Plotting Loop %
    for i = 1:MAX_ITER
        % Setup axis for plotting %
        panlist{i} = uipanel(f,'Position',[0 0.1 1 .9],...
                             'Tag','Panel with SUSD Plot');
        axlist{i} = axes('Parent',panlist{i},'Position',[.03 .05 .9 .9]);
        axis equal
        grid off
        ax = axlist{i};%gca;
        ax.XAxis.TickValues = x_min:10:x_max;
        ax.YAxis.TickValues = y_min:10:y_max;
        box = [x_min x_max y_min y_max];
        axis(box)
        hold on
        
        % Title %
        title(['Iteration ' num2str(i-1)])
        
        % Colorbar
        colorbar()
        
        % show the agent positions %
        for swarm=1:num_swarms
            plot(xx(1,:,swarm,i), xx(2,:,swarm,i), 'o', 'MarkerFaceColor',colors(swarm,:),'MarkerEdgeColor','k'); hold on;
        end
        % show the objective function %
        contour(X,Y,Z);

        % show graph in real time %
        drawnow;
    end
    
    if MAX_ITER < 2
       return 
    end
    
    % Construct GUI slider to view iterations %
    for i = 2:MAX_ITER
        set(panlist{i},'Visible','off')
    end
    ui = uicontrol('Parent',f,...
            'Style', 'slider',...
            'Min',0, ...
            'Max',MAX_ITER-1, ...
            'Value',[0],...
            'SliderStep',[1/(MAX_ITER-1) 1/(MAX_ITER-1)], ...
            'Position', [50 20 610 20],...
            'Callback', @src,...
            'UserData',struct('plots',[panlist{:}],'total',MAX_ITER));
    addlistener(ui, 'Value', 'PostSet', @drag);
end

function src(source,event)
    panlist = source.UserData.plots;
    MAX_ITER = source.UserData.total;
    for i = 1:MAX_ITER
        set(panlist(i),'Visible','off')
    end
    val = 1+floor(get(source, 'Value'));
    set(panlist(val),'Visible','on')
end

function drag(object,event)
    source = event.AffectedObject;
    panlist = source.UserData.plots;
    MAX_ITER = source.UserData.total;
    for i = 1:MAX_ITER
        set(panlist(i),'Visible','off')
    end
    val = 1+floor(get(source, 'Value'));
    set(panlist(val),'Visible','on')
end
    