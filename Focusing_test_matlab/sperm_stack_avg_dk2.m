function [average] = sperm_stack_avg_dk2(If,length_frame_stack)
%% Average Selected Bin Files


fprintf('Reading file: ...');
fprintf('Done \n');

A=If{1,1};

%read in all the other images
for i = 2:length_frame_stack
     
    fprintf('Reading file: %i...',i);
    B = If{1,i};
    
    fprintf('Done \n');
    A = A + B;		
end


C = A/length_frame_stack;

average = C;
end


