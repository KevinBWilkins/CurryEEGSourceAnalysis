function [Matrix,count,NR]=read_Curry_file4(filename,type,needContinue,Continued)
% This function to read the curry file.
% Input: the curry file name and the information type you need to read from that file.
%        the Type can be: LOCATION, NORMALS, CONTRIB, STRENGTHS, ERRORS, DEVIATIONS, and so on. 
% Output: A is the matrix of data read from the file
%         count is the number of the data read from the file

if ~Continued
    fpos=0;
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
%     fpostmp=fpos;
    load workspace
%     fpos=fpostmp;
%     clear fpostmp
    found=0;
    located=1;
    now_read=now_read-1;
    timePoint=timePoint+1;
    start=timenow+1;
end

keyNum=1;
NRnum=1;
fid=fopen(filename,'rt');

%ATTENTION------------------------------------------
%------------ATTENTION------------------------------
%------------------------ATTENTION
% fseek(fid,fpos,'bof');  %may not work for old curry files
fseek(fid,fpos-500,'bof');  %may not work for old curry files
%ATTENTION------------------------------------------
%------------ATTENTION------------------------------
%------------------------ATTENTION------------------


TLine=[];
while (~feof(fid) & end_header==0)
    if fpos==0 
        TLine=fgets(fid); 
    end
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
        end_header=1;
        break;
    elseif (~end_header & line_num>250)
        disp('May be a wrong file name')
        break
    end
end



for i_key = now_read:length(key)
    key_now=key{i_key};
    if (strcmp(key_now,[type,'_LIST']))
        now_read=i_key;
        
        if (now_read==1 | now_read==2)
            col=3;
        elseif (now_read==3 | now_read==4 | now_read==5)
            col=1;
        end
        
        if (now_read==2 | now_read==3 | now_read==4 | now_read==5)
            time_rep=NR(2);
        else
            time_rep=1;    
        end
        
        break;
    end
end

while (~feof(fid) & ~found)
    TLine=fgets(fid) ;
    
    if (~isempty(findstr([type,'_LIST START_LIST'],TLine)) | ~isempty(findstr([type,'_LIST NEW_TIMEPOINT'],TLine)))
        row=NR(1);
        
        for timenow=start:time_rep
            if (timenow==timePoint)
                
                [Matrix,count]=fscanf(fid,'%f',[col,row]);
                Matrix=Matrix';
                
                Matrixfirstvalues(timePoint) = Matrix(1);
                
                found=1;
                
                
                %---------------------------------------------------
                %ATTENTION------------------------------------------
                %------------ATTENTION------------------------------
                %------------------------ATTENTION------------------
                %TLine=fgets(fid);   %only needed for curry 5.0 files
                
                %ATTENTION------------------------------------------
                %------------ATTENTION------------------------------
                %------------------------ATTENTION------------------
                
                
                
                %TLine=fgets(fid)
                break
            else
                [temp,count_tmp]=fscanf(fid,'%f',[col,row]);
                %                             disp([fgets(fid),num2str(timenow)])
                TLine=fgets(fid);
                %TLine=fgets(fid)
            end
        end
        now_read=now_read+1;
    end


%     key_now=key{now_read};
%     
%     if ( key_now(2:3)=='O_')
%         now_read=now_read+1;
%     else
%         if (~isempty(findstr([type,'_LIST'],TLine)) & length(TLine)>4 )
%             row=NR(1);
%             
%             if (now_read==1 | now_read==2)
%                 col=3;
%             elseif (now_read==3 | now_read==4 | now_read==5)
%                 col=1;
%             end
%             if (now_read==2 | now_read==3 | now_read==4)
%                 time_rep=NR(2);
%             else
%                 time_rep=1;    
%             end
%             for timenow=start:time_rep
%                 if (strcmp(key_now,[type,'_LIST']) & timenow==timePoint)
%                     [Matrix,count]=fscanf(fid,'%f',[col,row]);
%                     Matrix=Matrix';
%                     found=1;
%                     fgets(fid);
%                     fgets(fid);
%                     break
%                 end
%             end
% %             line_num=line_num+time_rep*(row+2);
%             now_read=now_read+1;
%         end
%     end
end
fpos = ftell(fid);
if ~needContinue 
%     fclose(fid)
    if (Continued) 
        delete workspace
    end
else
    save workspace
end
fclose(fid);
