function Vark = mix_Cf(Vark,Cs0,BCs,dz)
    Vark.Cs= Cs0./((1./Vark.Kd).*Vark.phi + (1-Vark.phi));
    Vark.Cs = BC_setghosts(Vark.Cs,BCs.val_Cs,BCs.type_Cs,dz);
%   calculate Cf and Bulk
    Vark.Cf = Vark.Cs./Vark.Kd; 
    Vark.Cbulk= Vark.Cf.*Vark.phi + Vark.Cs .* (1-Vark.phi);
end