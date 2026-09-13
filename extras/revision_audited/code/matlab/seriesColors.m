function colors=seriesColors(n)
S=zstyle();p=[S.col.Mn;S.col.Co;S.col.Ni;S.col.Cu;S.col.Zn;S.col.Na;S.ink2];colors=p(1:n,:);
end
