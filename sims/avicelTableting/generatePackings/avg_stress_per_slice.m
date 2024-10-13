function [havg,SZZavg,SRRavg,tablet_height] = avg_stress_per_slice(data,stroke_factor,slices)

SXX = -data(:,1);
SYY = -data(:,2);
SZZ = -data(:,3);
V = data(:,4);
Z = data(:,5);

die_radius = 4e-3;
die_height = 1e-2;
die_area = pi*die_radius^2;
tablet_height = die_height*(1 - stroke_factor);

h = linspace(0,tablet_height,slices+1);

Nsum = zeros(length(h)-1,1);
SXXsum = zeros(length(h)-1,1);
SYYsum = zeros(length(h)-1,1);
SZZsum = zeros(length(h)-1,1);
Vsum = zeros(length(h)-1,1);

for i = 1 : length(Z)
    z = Z(i);
    v = V(i);
    sxx = SXX(i);
    syy = SYY(i);
    szz = SZZ(i);
    for j = 1 : length(h)-1
        if z >= h(j) && z <= h(j+1)
            Nsum(j) = Nsum(j) + 1;
            SXXsum(j) = SXXsum(j) + sxx;
            SYYsum(j) = SYYsum(j) + syy;
            SZZsum(j) = SZZsum(j) + szz;
            Vsum(j) = Vsum(j) + v;
            break
        end
    end
end

SXXavg = zeros(length(SXXsum),1);
SYYavg = zeros(length(SYYsum),1);
SZZavg = zeros(length(SZZsum),1);
havg = zeros(length(h)-1,1);

for i = 1 : length(h)-1
    Vslice = (h(i+1)-h(i))*die_area;
    havg(i) = (h(i+1)-h(i))/2 + h(i);
    SXXavg(i) = (SXXsum(i)/Nsum(i)) * (Vsum(i)/Vslice) * 1e-6;
    SYYavg(i) = (SYYsum(i)/Nsum(i)) * (Vsum(i)/Vslice) * 1e-6;
    SZZavg(i) = (SZZsum(i)/Nsum(i)) * (Vsum(i)/Vslice) * 1e-6;
end

SRRavg = (SXXavg + SYYavg)./2;

end