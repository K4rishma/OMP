kgrid = data{1,4};
Scan = data{1,1};
Trans = data{1,2};
medium = data{1,3};
pdi = data{1,8};
CI = data{1,7};
BF = data{1,6};
% Comparing for zero angle
ka = find(Scan.Angles == 0);
BF_angle(:,:,1:Scan.Ne) = BF(:,:,ka,1:Scan.Ne);

for j = 1:Scan.Ne
    BF_grid(:,:,j) = imresize(BF_angle(:,:,j),[kgrid.Nx kgrid.Ny]);
end


avg_no = Opts.n_ave;
skip = Opts.frameSkip;
beam_lines = BF_grid(:,kgrid.Ny/2 - 2 : kgrid.Ny/2 + 2, 1:avg_no + skip);

avg_abs_1 = mean(abs(beam_lines(:,:,1:avg_no)),3); % desired
avg_abs_2 = mean(abs(beam_lines(:,:,2:end)),3);

avg_absReal_1 = mean(abs(real(beam_lines(:,:,1:avg_no))),3); 
avg_absReal_2 = mean(abs(real(beam_lines(:,:,2:end))),3);

avg_Real_1 = mean(real(beam_lines(:,:,1:avg_no)),3); 
avg_Real_2 = mean(real(beam_lines(:,:,2:end)),3);

%observation points - hard coded 
mid = zeros(kgrid.Ny,1);

mid(61:69,1) = 1; % based on circular disc of width 10 pixels

figure; 
imagesc([avg_abs_1 mid avg_abs_2]); title('avg_abs');

figure; 
imagesc(avg_absReal_1 mid avg_absReal_2]); title('avg_absReal');

figure; 
imagesc([avg_Real_1 mid avg_Real_2]); title('avg_Real');


subplot(3,1,1);
    plot(data_lines(1:end - kgrid.Nt/2,1));
    xlabel('Depth');
    ylabel('Amplitute');
    title(['Channel lines at angle ' num2str(Scan.Angles(1)) ])
    
subplot(3,1,2);
    plot(noise_filter_lines(1:end - kgrid.Nt/2,1));
    xlabel('Depth');
    ylabel('Amplitute');
    title(['Noise filter lines at angle ' num2str(Scan.Angles(1)) ])
    
subplot(3,1,3);
    plot(hilber_env_lines(1:end - kgrid.Nt/2,1));
    xlabel('Depth');
    ylabel('Amplitute');
    title(['Hilbert envelope lines at angle ' num2str(Scan.Angles(1)) ])