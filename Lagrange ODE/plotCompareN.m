function plotCompareN()
% plot results related to compareN

load('data/compareN.mat');

% retrieve data
NxIter = zeros(length(N),length(rkOrder));
NxFval = zeros(length(N),length(rkOrder));
NxPerr = zeros(length(N),length(rkOrder));

for k=1:length(rkOrder)
    for n=1:length(N)
        if res_sqp{k}{n}{3} == 1 & res_sqp{k}{n}{2} < 1
            % converged
            NxIter(n,k) = res_sqp{k}{n}{4}.iterations;
            NxFval(n,k) = res_sqp{k}{n}{2};
            NxPerr(n,k) = norm(res_sqp{k}{n}{1}-p_opt)/norm(p_opt);
        end
        % else leave it 0
    end
end


% plot iter over N
figure;
hold on;

for k=1:length(rkOrder)
    plot(N,NxIter(:,k));
end

% plot fval over N
figure;
hold on;

for k=1:length(rkOrder)
    plot(N,NxFval(:,k));
end

% relative p error over N
figure;
hold on;

for k=1:length(rkOrder)
    plot(N,NxPerr(:,k));
end

end
