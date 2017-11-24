function UnlinkedSampleTones

delete(instrfindall);

s = serial('/dev/tty.usbserial', 'BaudRate', 9600, 'DataBits', 8, ...
    'Parity', 'None', 'StopBits', 1, 'FlowControl', 'None', 'Terminator',...
    'LF');

fopen(s);

%high, left speaker, left blinks
pause(15);
fprintf(s, '%s\n', 'START.1000.2000.100.TARGET.ONE.STOP');
pause(2);
fprintf(s, '%s\n', 'go');
pause(5.6);
fprintf(s, '%s\n', 'no');

%high, right speaker, right blinks
pause(7);
fprintf(s, '%s\n', 'START.1000.2000.100.TARGET.TWO.STOP');
pause(2);
fprintf(s, '%s\n', 'go');
pause(5.6);
fprintf(s, '%s\n', 'no');

%low, left speaker, right blinks
pause(7);
fprintf(s, '%s\n', 'START.1000.2000.0.MASKER.ONE.STOP');
pause(2);
fprintf(s, '%s\n', 'go');
pause(5.6);
fprintf(s, '%s\n', 'no');

%low, right speaker, left blinks
pause(7);
fprintf(s, '%s\n', 'START.1000.2000.0.MASKER.TWO.STOP');
pause(2);
fprintf(s, '%s\n', 'go');
pause(5.6);
fprintf(s, '%s\n', 'no');

pause(7);
fprintf(s, '%s\n', 'done');

end

