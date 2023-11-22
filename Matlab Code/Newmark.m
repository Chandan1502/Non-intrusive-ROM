function [U ,Udot,Uddot] = Newmark(M,C,K,F,dt,N,Uo,Udoto )
%% Newmark Beta Method
% Size of M, C and K is Ndof x Ndof.
% Size of F is Ndof x N
% Uo is intial displacement
% Udoto is initial velocity

[row,~]=size(F);
dof=row;
f=F;     

U=zeros(dof,N); % Displacement
Udot=zeros(dof,N); % Velocity
Uddot=zeros(dof,N);

U(:,1)=Uo; 
Udot(:,1)=Udoto;
Uddot(:,1)=(M)\(f(:,1)-C*Udot(:,1)-K*U(:,1));

delta=.5;
alpha=0.25;

a0=1/(alpha*dt^2);
a1=delta/(alpha*dt);
a2=1/(alpha*dt);
a3=(1/(2*alpha)-1);
a4=delta/alpha-1;
a5=dt/2*(delta/alpha-2);
a6=dt*(1-delta);
a7=delta*dt;
MM=a0*M+a1*C+K;
MMM=inv(MM);

for i=1:N-1
    Mt=M*(a0*U(:,i)+a2*Udot(:,i)+a3*Uddot(:,i));
    Ct=C*(a1*U(:,i)+a4*Udot(:,i)+a5*Uddot(:,i));
    U(:,i+1)=(MMM)*(f(:,i+1)+Mt+Ct);
    Udot(:,i+1)=Udot(:,i)+dt*((1-delta)*Uddot(:,i)+delta*a0*(U(:,i+1)-U(:,i)-Udot(:,i)*dt)-delta/alpha*(1/2-alpha)*Uddot(:,i));
    Uddot(:,i+1)=a0*(U(:,i+1)-U(:,i)-Udot(:,i)*dt)-(1/alpha)*(1/2-alpha)*Uddot(:,i);
end

