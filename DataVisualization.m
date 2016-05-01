function [ A ] = DataVisualization(params, lr1, lr2) %rope1, rope2,

%lr2 = @(t) params.lr20;
%lr2 = params.lr20*ones(length(lr1),1);

figure(1);
numframes=length(lr1);
A=moviein(numframes); % create the movie matrix
% set(gca,'NextPlot','replacechildren')
% axis equal % fix the axes
for i=1:numframes
    
    % plot of rope length
    subplot(1,2,1)
    
    hold on;
    %fplot(lr1,[1,10]);
    %fplot(lr2,[1,10]);
    plot([1:length(lr1)],lr1);
    plot([1:length(lr1)],lr2);
    plot([i,i],[0,30],'r');
    
    % Visualization of actual excavator configuration
    [alpha,l,pos_load] = ExcavatorConfig(lr1(i),lr2(i),params );
    subplot(1,2,2);
    hold on;
    axis equal;
    axis([0,15,0,30]);
    % rope 1, base to top pulley to load position
    plot([0,params.l3*cos(params.ang_base),pos_load(1)],[0,params.l3*sin(params.ang_base),pos_load(2)],'linewidth',1.5);
    % sidearm, stops at base arm
    %plot([params.l1*cos(params.ang_base),pos_load(1)],[params.l1*sin(params.ang_base),pos_load(2)],'linewidth',1.5);

    % sidearm
    BC_vect = [params.l1*cos(params.ang_base) - pos_load(1);
               params.l1*sin(params.ang_base) - pos_load(2)];
    Z = params.l4*BC_vect/norm(BC_vect) + [pos_load(1); pos_load(2)];   % upper endpoint of sidearm
    plot([Z(1),pos_load(1)],[Z(2),pos_load(2)],'linewidth',1.5);
    
    
    A(:,i)=getframe(gcf);
    clf;
end

axes('pos',[0 0 1 1],'visible','off');
movie(A,1,2) % Play the MATLAB movie
% save movie.mat A % save the MATLAB movie to a file
% mpgwrite(A,jet,'movie.mpg'); % Convert the movie to MPEG format
% % Notice the MPEG file is about a quarter of the size of the MATLAB movie file
% unix('mpeg_play movie.mpg') % Play the MPEG movie 

end

