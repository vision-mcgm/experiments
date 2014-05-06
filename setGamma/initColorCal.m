%vsgInit; % try this if init does not work
colorCALinit;
[ErrorCode]= colorCALautocalibrate;
CIEcolour = colorCALread;
