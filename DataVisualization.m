function [ ] = DataVisualization(params) %rope1, rope2,

lr1 = @(t) 19 + 0.5*t;
lr2 = @(t) 12;

figure(1);
numframes=10;
A=moviein(numframes); % create the movie matrix
% set(gca,'NextPlot','replacechildren')
% axis equal % fix the axes
for i=1:numframes
    
    % plot of rope length
    subplot(1,2,1)
    
    hold on;
    fplot(lr1,[1,10]);
    fplot(lr2,[1,10]);
    plot([i,i],[0,30],'r');
    
    % Visualization of actual excavator configuration
    [alpha,l,pos_load] = ExcavatorConfig(lr1(i),lr2(i),params );
    subplot(1,2,2);
    hold on;
    plot([0,params.l3*cos(params.ang_base),pos_load(1)],[0,params.l3*sin(params.ang_base),pos_load(2)],'linewidth',1.5);
    plot([params.l1*cos(params.ang_base),pos_load(1)],[params.l1*sin(params.ang_base),pos_load(2)],'linewidth',1.5);
    
    A(:,i)=getframe(gcf);
    clf;
end
movie(A,1,2) % Play the MATLAB movie
% save movie.mat A % save the MATLAB movie to a file
% mpgwrite(A,jet,'movie.mpg'); % Convert the movie to MPEG format
% % Notice the MPEG file is about a quarter of the size of the MATLAB movie file
% unix('mpeg_play movie.mpg') % Play the MPEG movie 

end

