function [K,M,F] = Beam_element_assembly(E,I,Ar,ro,L,Nel,q)
%  BEAM_ELEMENT_ASSEMBLY This function is used to assemble the mass and stiffness matrix.
% Input :
% E= Modulus of elasticity
% I= Area momemt of intertia
% Ar= c/s area of the beam
% ro= density of the beam
% L= length of the beam
% Nel= number of elements
% q= external load on beam (uniform)

% The output of this function is Stiffness, mass matrices and nodal force vector.

m=ro*Ar;                % m= mass per unit length
l=L/Nel;
Mel=m*l/420*[156 22*l 54 -13*l;
    22*l 4*l^2 13*l -3*l^2;     % Mel= Beam Element mass matrix
    54 13*l 156 -22*l;
    -13*l -3*l^2 -22*l 4*l^2];
Kel=(E*I/(l^3))*[12 6*l -12 6*l;
    6*l 4*l^2 -6*l 2*l^2;
    -12 -6*l 12 -6*l;         % Kel=stiffness matrix of beam element
    6*l 2*l^2 -6*l 4*l^2];
% Generation of nodal force vector
% beam is uniformly loaded
% Gauss quardature is used for the integration
zeta=[-1/sqrt(3) 1/sqrt(3)];
Fel(1)=q*l/2*(1-3/4*(1+zeta(1))^2+1/4*(1+zeta(1))^3+1-3/4*(1+zeta(2))^2+1/4*(1+zeta(2))^3);
Fel(2)=q*l*l/2*(1/2*(1+zeta(1))-1/2*(1+zeta(1))^2+1/8*(1+zeta(1))^3+1/2*(1+zeta(2))-1/2*(1+zeta(2))^2+1/8*(1+zeta(2))^3);
Fel(3)=q*l/2*(3/4*(1+zeta(1))^2-1/4*(1+zeta(1))^3+3/4*(1+zeta(2))^2-1/4*(1+zeta(2))^3);
Fel(4)=q*l*l/2*(-1/4*(1+zeta(1))^2+1/8*(1+zeta(1))^3-1/4*(1+zeta(2))^2+1/8*(1+zeta(2))^3);
Fel=Fel';
Ndof=2*(Nel+1);
K=zeros(Ndof,Ndof);
M=zeros(Ndof,Ndof);
F=zeros(Ndof,1);
for i=1:Nel
    K(2*i-1:4+2*(i-1),2*i-1:4+2*(i-1))=K(2*i-1:4+2*(i-1),2*i-1:4+2*(i-1))+Kel;
    M(2*i-1:4+2*(i-1),2*i-1:4+2*(i-1))=M(2*i-1:4+2*(i-1),2*i-1:4+2*(i-1))+Mel;    % efficeint approach
    F(2*i-1:4+2*(i-1),1)=F(2*i-1:4+2*(i-1),1)+Fel;
end
end