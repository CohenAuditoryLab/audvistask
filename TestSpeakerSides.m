delete(instrfindall);

s = serial('/dev/tty.usbserial', 'BaudRate', 9600, 'DataBits', 8, ...
    'Parity', 'None', 'StopBits', 1, 'FlowControl', 'None', 'Terminator',...
    'LF');

fopen(s);

pause(2);

%then try again with TWO
fprintf(s, '%s\n', 'START.1000.2000.50.All.TWO.STOP');

pause(2);

fprintf(s, '%s\n', 'go');

%issue a print 'no' to stop once you check and a done to quit 