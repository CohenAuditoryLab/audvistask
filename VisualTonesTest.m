f_0 = [];

for k = 1 : 400
    [~, ~, ~, ~, bursts] = VisualTones(1000, 2000, 0, 'None');
    f = bursts(2:2:end);
    f_0(end+1:end+length(f)) = f;
end 

% figure()
% histogram(f_0)
% title('Coherence Test, 0%')
figure()
histogram(f_0)
title('Interval Normalcy Test')

%%%%%%%%%%%%%%%%%%

% f_50 = zeros(1, 1000);
% i_50 = zeros(1, 1000);
% 
% for k = 1 : 1000
%     [f, int] = VisualTones(1000, 10000, .25);
%     f_50(k) = f;
%     i_50(k) = int;
% end 
% 
% figure()
% histogram(f_50)
% title('Coherence Test, 25%')
% figure()
% histogram(i_50)
% title('Interval Normalcy Test, 25%')
% 
% %%%%%%%%%%%%%%%%%%%%
% 
% f_50 = zeros(1, 1000);
% i_50 = zeros(1, 1000);
% 
% for k = 1 : 1000
%     [f, int] = VisualTones(1000, 10000, .5);
%     f_50(k) = f;
%     i_50(k) = int;
% end 
% 
% figure()
% histogram(f_50)
% title('Coherence Test, 50%')
% figure()
% histogram(i_50)
% title('Interval Normalcy Test, 50%')
% 
% %%%%%%%%%%%%%%%%%%%%
% 
% f_100 = zeros(1, 1000);
% i_100 = zeros(1, 1000);
% 
% for k = 1 : 1000
%     [f, int] = VisualTones(1000, 10000, .75);
%     f_100(k) = f;
%     i_100(k) = int;
% end 
% 
% figure()
% histogram(f_100)
% title('Coherence Test, 75%')
% figure()
% histogram(i_100)
% title('Interval Normalcy Test, 75%')
% 
% %%%%%%%%%%%%%%%%%
% f_100 = zeros(1, 1000);
% i_100 = zeros(1, 1000);
% 
% for k = 1 : 1000
%     [f, int] = VisualTones(1000, 10000, 1);
%     f_100(k) = f;
%     i_100(k) = int;
% end 
% 
% figure()
% histogram(f_100)
% title('Coherence Test, 100%')
% figure()
% histogram(i_100)
% title('Interval Normalcy Test, 100%')
