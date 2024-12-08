function plotWeibullDistribution(sphereRadii,R_factor)

data = readmatrix("PH102 PSD.csv");

size = R_factor*data(:,1);
vol_per = data(:,2);
vol_frac = vol_per/100;
count_unity = vol_frac./(4/3*pi*(size/2).^3);
pdf_exp = count_unity/sum(count_unity);

linewidth = 1;
gcaFontsize = 20;
labelFontsize = 24;
legendFontsize = 22;

% 1. Generate particle sizes from a Weibull distribution
particles = sphereRadii*1e6;  % 1000 particle sizes

% 2. Define bin edges (adjust these based on your range)
binEdges = vertcat(0,size);  % Example bin edges in microns

% 3. Calculate particle volumes (assuming spherical particles)
volumes = (4/3) * pi * (particles / 2).^3;

% 4. Bin the particles by size
[N, edges, binIdx] = histcounts(particles, binEdges);
N = N';

% 5. Sum the volumes in each bin
binVolumes = zeros(length(N),1);  % Initialize array to store volume per bin

for i = 1:length(N)
    binVolumes(i) = sum(volumes(binIdx == i));
end

% 6. Calculate total volume
totalVolume = sum(volumes);

% 7. Calculate volume percentage per bin
volumePercent = (binVolumes / totalVolume) * 100;
volumeFraction = volumePercent/100;
countUnity = volumeFraction./(4/3*pi*(size/2).^3);
pdf_fit = countUnity/sum(countUnity);

figure()
tiledlayout(2,2)

nexttile()
plot(size,vol_per,'-o','LineWidth',linewidth)
hold on
plot(size,volumePercent,'-s','LineWidth',linewidth)
set(gcf,'color','w');
set(gca, 'FontSize', gcaFontsize)
xlabel('size [$\mu$m]','Interpreter','latex','FontSize', labelFontsize);
ylabel('volume distribution [\%]','Interpreter','latex','FontSize', labelFontsize);
set(gca, 'TickLabelInterpreter','latex','XMinorTick','on','YMinorTick','on','Fontsize',gcaFontsize)
xlim([0 max(size)/6])
box on
hl = legend('experiment','simulation');
set(hl,'FontSize',legendFontsize,'Interpreter','latex','Location','NorthEast')

nexttile()
plot(size,pdf_exp,'-o','LineWidth',linewidth)
hold on
plot(size,pdf_fit,'-s','LineWidth',linewidth)
set(gcf,'color','w');
set(gca, 'FontSize', gcaFontsize)
xlabel('size [$\mu$m]','Interpreter','latex','FontSize', labelFontsize);
ylabel('probability density function','Interpreter','latex','FontSize', labelFontsize);
set(gca, 'TickLabelInterpreter','latex','XMinorTick','on','YMinorTick','on','Fontsize',gcaFontsize)
xlim([0 max(size)/6])
ylim([0 0.15])
box on
hl = legend('experiment','simulation');
set(hl,'FontSize',legendFontsize,'Interpreter','latex','Location','NorthEast')

end