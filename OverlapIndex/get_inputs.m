function output = get_inputs (matrix)

output =[];
for t=7:11
    for task =1:2
        tmp = squeeze(matrix(:,task,:,t));
        output = cat(2, output, tmp);
    end
   
end
save ('E:\Documents\paper\overlap-time\tmp_output.txt','output', '-ascii');

     