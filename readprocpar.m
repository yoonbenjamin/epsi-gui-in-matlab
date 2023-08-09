function a = readprocpar(f, param)
% syntax a = readprocpar(filepath, param)
% Before the end of file, read a line in propar file, if the first 
% character in that line is parameter read and return the value on the
% next line!    
% g = fopen([f '.fid\procpar']);      % use for pc
    g = fopen(strcat(f,'.fid/procpar'));    % use for mac
    while(~feof(g))
        s = fgets(g);
        if (strcmp(strtok(s),param))
            a = strread(fgets(g));
        end
    end
