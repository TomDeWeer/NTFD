function fTangent = calc_fTangent(nF,fNormal)
    % via rechterhandregel: t cross n = ez
    fTangent = [];
    for i = 1:nF
        fNormali = fNormal(2*i-1:2*i);
        fTangenti = [fNormali(2); -fNormali(1)];
        % teken bepalen adhv rechterhand
        n = [fNormali; 0];
        t = [fTangenti; 0];
        ez = cross(t,n);
        fTangenti = fTangenti * sign(ez(3));
        fTangent = [fTangent fTangenti];
    end
end



