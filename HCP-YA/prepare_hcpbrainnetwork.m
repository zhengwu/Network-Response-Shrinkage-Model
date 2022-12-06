
clear all;
close all;

load HCP_networkdata.mat
subN = length(all_id);

mean_network = mean(all_network,3);
figure, imagesc(log(mean_network+1));

mask = mean_network>0.2;
figure, imagesc(mask);

%get non-zeros in the mask
idx = 1;
for i=1:size(mask,1)-1
    for j=i+1:size(mask,1)
        if(mask(i,j)>0)
            no_zeros_pairs(idx,1) = i;
            no_zeros_pairs(idx,2) = j;
            idx = idx +1;
        end
    end
end


Npair = idx - 1;
%vectorized network data
for i=1:subN
    for j=1:Npair
        all_network_vector(i,j) = all_network(no_zeros_pairs(j,1),no_zeros_pairs(j,2),i);
    end
end


%get the prior connectivity
thrd = 0;
shared_network = zeros(size(mask));
for i=1:subN
   binary_network(:,:,i) =  all_network(:,:,i)>thrd;
   shared_network = shared_network + binary_network(:,:,i);
end

binary_shared_network = shared_network>(subN-0.5);
figure;imagesc(binary_shared_network);


%get non-zeros in the binary_shared_network
idx = 1;
for i=1:size(binary_shared_network,1)
    for j=1:size(binary_shared_network,1)
        if(binary_shared_network(i,j)>0)
            E_no_zeros_pairs(idx,1) = i;
            E_no_zeros_pairs(idx,2) = j;
            idx = idx +1;
        end
    end
end


save prepared_fa_networkdata X_measure all_network subN all_id mask no_zeros_pairs all_network_vector E_no_zeros_pairs;



