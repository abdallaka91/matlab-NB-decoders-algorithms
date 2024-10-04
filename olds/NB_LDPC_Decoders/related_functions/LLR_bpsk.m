function [data] = LLR_bpsk(Q, cdata, sigma)
   p = log2(Q);
   data_to_compare = (-1).^de2bi(0:Q-1, p,2);
   N = length(cdata)/p;
   data = zeros(Q,N);
   for i = 1:N
      input = cdata(p*(i-1)+1:p*i);
      temp = data_to_compare - repmat(input, Q, 1);
      temp = temp.^2;
      square_distances = sum(temp, 2);
      llr_one_symbol = -square_distances/(2*sigma^2);
      llr_one_symbol = llr_one_symbol - llr_one_symbol(1);
      data(:,i) = llr_one_symbol;
   end

% function [data] = LLR_bpsk(Q, cdata, sigma)
% p = log2(Q);
% data_to_compare = (-1).^de2bi(0:Q-1, p,2);
% N = length(cdata)/p;
% data = zeros(Q,N);
% for i = 1:N
%     input = cdata(p*(i-1)+1:p*i);
%     temp = data_to_compare - repmat(input, Q, 1);
%     temp = temp.^2;
%     square_distances = sum(temp, 2);
%     llr_one_symbol = -square_distances/(2*sigma^2);
%     data(:,i) = llr_one_symbol;
% end
% 
% for i = 1:N
%     data(:,i) = data(:,i) - max(data(:,i));
% end