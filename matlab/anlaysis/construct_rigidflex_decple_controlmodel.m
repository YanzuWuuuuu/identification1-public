function [G,G11,G12,G21,G22] = construct_rigidflex_decple_controlmodel(nrigidcontrol,nflexcontrol,w2_ss,zeta_w_2_ss,...
    Br,Cr,Bf,Cf,Tu,Ty,Pcnull,Pbnull)

Ar1 = [zeros(nrigidcontrol,nrigidcontrol),eye(nrigidcontrol,nrigidcontrol);
      zeros(nrigidcontrol,nrigidcontrol),zeros(nrigidcontrol,nrigidcontrol)];

Br1 = [zeros(nrigidcontrol,nrigidcontrol);eye(nrigidcontrol,nrigidcontrol)];
Cr1 = [eye(nrigidcontrol,nrigidcontrol),zeros(nrigidcontrol,nrigidcontrol)];
Dr1 = zeros(nrigidcontrol,nrigidcontrol);

Af1 = [zeros(nflexcontrol,nflexcontrol),eye(nflexcontrol,nflexcontrol);
      -w2_ss(4:3+nflexcontrol,4:3+nflexcontrol),-zeta_w_2_ss(4:3+nflexcontrol,4:3+nflexcontrol)];

Bf1 = [zeros(nflexcontrol,nflexcontrol);eye(nflexcontrol,nflexcontrol)];
Cf1 = [eye(nflexcontrol,nflexcontrol),zeros(nflexcontrol,nflexcontrol)];
Df1 = zeros(nflexcontrol,nflexcontrol);

Pr = ss(Ar1,Br1,Cr1,Dr1);
Pf = ss(Af1,Bf1,Cf1,Df1);

G11 = Ty*(Cr*Pr*Br + Cf*Pf*Bf)*Tu;
G12 = Ty*Cf*Pf*Bf*Pbnull;
G21 = Pcnull*Cf*Pf*Bf*Tu;
G22 = Pcnull*Cf*Pf*Bf*Pbnull;

G = [G11,G12;
     G21,G22];
end