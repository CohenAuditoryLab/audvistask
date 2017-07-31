delete(instrfindall);
s = serial('/dev/tty.usbserial', 'BaudRate', 9600, 'DataBits', 8,...
    'Parity', 'None', 'StopBits', 1, 'FlowControl', 'None');
fopen(s);
pause(2)

t = 'hello';

for i = 1:100
    fprintf(s, t); %'%s\n', 'hello')
    disp('sent!');
    pause(.2)
end
   
%d = fread(s, 6);
fclose(s);

%-------------------------

% data = 'hello';
% st = whos('data');
% 
% delete(instrfindall);
% s = tcpip('0.0.0.0',61900,'NetworkRole','server'); %130.91.169.182 %169.254.74.85
% set(s,'OutputBufferSize',st.bytes);
% fopen(s);
% disp('opened');
% 
% fwrite(s, data(:), 'char');
% fclose(s);

%--------------------------
% delete(instrfindall);
% % Configuration and connection
% t = udp('130.91.168.1', 6666, 'LocalPort', 6665);
% %tcpip('130.91.169.182', 61900, 'NetworkRole', 'Client');
% 
% % Open socket and wait before sending data
% fopen(t);
% pause(0.2);
% 
% 
% % Send data every 500ms
% for i=0:100    
%     DataToSend = 'hello';
%     fwrite(t, DataToSend);
% %     DataToRead = char(fread(t, 6))';
% %     disp(['entry ' num2str(i) ': ' DataToRead]);
%      pause (1);
% end
% 
% % Close and delete connection
% fclose(t);
% %delete(t);