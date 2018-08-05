function [Matrix,Matrix_full,count,NR]=read_Curry_file3_AC(filename,type,needContinue,Continued)
% This function to read the curry file.
% Input: the curry file name and the information type you need to read from that file.
%        the Type can be: LOCATION, NORMALS, CONTRIB, STRENGTHS, ERRORS, DEVIATIONS, and so on. 
% Output: A is the matrix of data read from the file
%         count is the number of the data read from the file
if ~Continued
    begin=0;
    matched=0;
    found=0;
    located=0;
    A=[];
    count=0;
    line_num=0;
    got_info=0;
    already_jump=0;
    already_judged=0;
    fid=fopen(filename,'rt');
    now_read=1;
    start=1;
    timePoint=1;
    end_header=0;
else
    load workspace;
    found=0;
    located=1;
    now_read=now_read-1;
    timePoint=timePoint+1;
    start=timenow+1;
end

keyNum=1;
NRnum=1;
fid=fopen(filename,'rt');
while (~feof(fid) & ~found)
    TLine=fgets(fid);
    %disp(TLine);
    line_num=line_num+1;
    
    
    if (line_num<250 & ~begin)
        if ~isempty( findstr(['POINT_KEYWORDS START'],TLine) )  
            begin=1; 
        end
    elseif (line_num<250 & ~isempty( findstr(['POINT_KEY'],TLine) ) )
        key{keyNum}=TLine(strfind(TLine,'=')+2:end-1);
        keyNum=keyNum+1;
    elseif (line_num<250 & ~isempty( findstr(['POINT_NR'],TLine) ) )
        NR(NRnum)=str2num(TLine(strfind(TLine,'=')+1:end));
        NRnum=NRnum+1;
    elseif (~isempty( findstr(['POINT_TRAFO END_LIST'],TLine)) )
        %end_header=1;
    elseif (line_num<250 & ~isempty( findstr(['LIST_NR_ROWS'],TLine) ) )
        NR_locations(1)=str2num(TLine(strfind(TLine,'=')+1:end));
        end_header=1;
    elseif (~end_header & line_num>250)
        disp('May be a wrong file name')
        break
    end
    if (end_header)
        %Now judge whether the type is included
        if (~already_judged)
            for i=1:15
                if ( ~isempty(findstr([type,'_LIST'],key{i}))) 
                    matched=1; match_line=i+1;break
                end
            end
            already_judged=1;
        end
        if (~matched)
            disp('May be a wrong type')

            break
        else
            %if (~already_jump)
            %    fseek(fid,46*37,-1);
                already_jump=1;
                %end

            key_now=key{now_read};
            if ( ~((key_now(2)=='O') * (key_now(3)=='_')) )
                if (~isempty(findstr([key_now,' START_LIST'],TLine)) & length(TLine)>4 )
                    located=1; 
                end
                if (located)
                    row=NR_locations(1);

                    if (now_read==1 | now_read==2)
                        col=3;
                    elseif (now_read==3 | now_read==4)
                        col=1;
                    end
                    if (now_read==2 | now_read==3 | now_read==4)
                        time_rep=NR(2);
                    else
                        time_rep=1;    
                    end
                    for timenow=start:time_rep
                        if (strcmp(key_now,[type,'_LIST']) & timenow==timePoint)
                            [Matrix,count]=fscanf(fid,'%f',[col,row]);
                            Matrix=Matrix';
                            found=1;
                            Matrix_full = Matrix;
                            Matrix = Matrix(1:NR(1),:);
                            %disp([fgets(fid),num2str(timenow)]);
                            
                            break
                        else
                            [temp,count_tmp]=fscanf(fid,'%f',[row,col]);
                            %disp([fgets(fid),num2str(timenow)]);
                            fgets(fid);
                        end
                    end
                    line_num=line_num+time_rep*(row+2);
                    now_read=now_read+1;
                    located=0;
                end
            end
        end
    end
end
if ~needContinue 
    fclose(fid);
    if (Continued) 
        delete workspace
    end
else
    save workspace
end