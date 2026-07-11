ccc

fname = '/Users/lac/Desktop/temp/MandaHomo/examples/Manda/sourceFile.txt'
fid = fopen(fname,'r');
N = fscanf(fid, '%d', 1);
source = fscanf(fid, '%f', [3, N])';
fclose(fid);

fname = '/Users/lac/Desktop/temp/MandaHomo/examples/Manda/source_corrector_file'
fid = fopen(fname,'r');
Nsrccorr = fscanf(fid, '%d', 1);
srcGn = fscanf(fid, '%f', [36, Inf])';
fclose(fid);


src = zeros(Nsrccorr,6,6);
for is = 1:Nsrccorr
    src(is,:,:) = reshape(srcGn(is,:),6,6);
end

m = [1 1 1 sqrt(2)*0 sqrt(2)*0 sqrt(2)*0];
M = zeros(Nsrccorr,3,3);
for is = 1:Nsrccorr
    corsrc = squeeze(src(is,:,:));
    for i = 1:6
        mcor(i) = sum(corsrc(:,i).*m(:));
    end
    M(is,:,:) = [mcor(1) mcor(6)/sqrt(2) mcor(5)/sqrt(2); ...
        mcor(6)/sqrt(2) mcor(2) mcor(4)/sqrt(2);...
        mcor(5)/sqrt(2) mcor(4)/sqrt(2) mcor(3)];
end


%%

fname = '/Users/lac/Desktop/temp/MandaHomo/examples/Manda/receiver_corrector_file'
fid = fopen(fname,'r');
Nrcvcorr = fscanf(fid,'%d',1);
rcvGn = fscanf(fid,'%f',[18 Inf])';
fclose(fid);

fname = '/Users/lac/Desktop/temp/MandaHomo/examples/Manda/monitorFile.txt'
fid = fopen(fname,'r');
N = fscanf(fid, '%d', 1);
Monitor = fscanf(fid, '%f', [3, N])';
fclose(fid);

fname = '/Users/lac/Desktop/temp/MandaHomo/examples/Manda/receiver_corrector_file'
fid = fopen(fname,'r');
Nrcvcorr = fscanf(fid,'%d',1);
rcvGn = fscanf(fid,'%f',[18 Inf])';
fclose(fid);

rcv = zeros(Nrcvcorr,3,6);
for ir = 1:Nsrccorr
    rcv(ir,:,:) = reshape(rcvGn(ir,:),3,6);
end

for ir = 1:Nsrccorr
    Cr   = squeeze(rcv(ir,:,:));                       % 3 x 6
    s2   = sqrt(2);
    eK   = epsilon .* [1 1 1 s2 s2 s2];                  % Nt x 6 (Kelvin)
    %u_corr = desloc + eK * Cr.';                         % Nt x 3
    % (u_corr(:,i) = u_hom(:,i) + Σ_J rcvGn(i,J)·eK(:,J))
end