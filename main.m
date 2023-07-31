clc
clear all
close all

xMin = 0;
xMax = 1;

rMin = 0;%15e-3;
rMax = 1;%19e-3;
ncellx =30;
ncellr = 30;

uLid = [1,0];

%% initialization
% generate mesh
% 1 left - 2 down - 3 right - 4 up
dx = (xMax - xMin) / ncellx;
dr = (rMax - rMin) / ncellr;
[x,r] = meshgrid(linspace(xMin,xMax,ncellx+1),linspace(rMin,rMax,ncellr+1));

for i = 1:ncellr
    for j = 1:ncellx
        volume(i,j) = dr * dx; % 2-dim volume
        faceArea(i,j,1) = dr;
        faceArea(i,j,3) = dr;
        faceArea(i,j,4) = dx;
        faceArea(i,j,2) = dx;

    end
end
% determine the B.C. type of each face:
bType = -1*ones(ncellr,ncellx,4);
bType(1,:,4) = 0;	% 4 means up boundary
bType(:,1,1) = 0;	% 0 means left wall
bType(:,ncellx,3) = 0;	% 0 means right wall
bType(ncellr,:,2) = 4;	% 0 means down wall

% normal vecter, colume vector
normal = [[-1;0],[0;1],[1;0],[0;-1]];

% neighbour index offset and distance
ioffset = [0,1,0,-1];
joffset = [-1,0,1,0];
d = [dx,dr,dx,dr];

disp('mesh has been generated*********')

nu = 1e-2; % kinematic viscosity = dynamic viscosity / rho

disp('properties are assigned *******')

u = zeros(ncellr,ncellx,2);
p = zeros(ncellr,ncellx);

disp('initialization is done********')
relaxP = 0.6; % Jasak suggest 0.2 on thesis page 149
relaxU = 0.2;

%% calculation
% variables initialization
F = zeros(ncellr,ncellx,4);
VbyA = zeros(ncellr,ncellx);
VbyA_f = zeros(ncellr,ncellx,4);
HbyA_f = zeros(ncellr,ncellx,4);

for iter = 1:10000

    gradP = grad(ncellx,ncellr,dx,dr,p);

    A = sparse(zeros(ncellx*ncellr,ncellx*ncellr));
    b = zeros(ncellx*ncellr,2);
   
    for i = 1:ncellr
        for j = 1:ncellx
            imat = (i-1)*ncellx + j;
            for k = 1:4
                jmat = (i+ioffset(k)-1)*ncellx + j+joffset(k); % determine which neighbour delt with in this cycle
                if bType(i,j,k) == -1
                    uf(1) = 0.5* (u(i,j,1) + u(i+ioffset(k),j+joffset(k),1));
                    uf(2) = 0.5* (u(i,j,2) + u(i+ioffset(k),j+joffset(k),2));
                    F(i,j,k) = uf*normal(:,k)*faceArea(i,j,k);

                    A(imat,jmat) = F(i,j,k)/2 - abs(F(i,j,k)/2) -nu*faceArea(i,j,k)/d(k);
                    A(imat,imat) = A(imat,imat) + (F(i,j,k)/2 +abs(F(i,j,k)/2) + nu*faceArea(i,j,k)/d(k));

                elseif bType(i,j,k) == 0
                    F(i,j,k) = 0;
                    A(imat,imat) = A(imat,imat) + nu*faceArea(i,j,k)/(d(k)/2);

                elseif bType(i,j,k) == 4
                    F(i,j,k) = uLid*normal(:,k)*faceArea(i,j,k);
                    A(imat,imat) = A(imat,imat) + nu*faceArea(i,j,k)/(d(k)/2);
                    b(imat,1) = b(imat,1) - F(i,j,k)*uLid(1) +nu*faceArea(i,j,k)*uLid(1)/(d(k)/2);
                    b(imat,2) = b(imat,2) - F(i,j,k)*uLid(2) +nu*faceArea(i,j,k)*uLid(2)/(d(k)/2);

                else
                    error('Error occurred.\n Undefined boundary!');
                end
            end

            b(imat,1) = b(imat,1) - gradP(i,j,1)*volume(i,j);
            b(imat,2) = b(imat,2) - gradP(i,j,2)*volume(i,j);

            %relax

            b(imat,1) =b(imat,1) + u(i,j,1)*A(imat,imat) * (1-relaxU)/relaxU;
            b(imat,2) =b(imat,2) + u(i,j,2)*A(imat,imat) * (1-relaxU)/relaxU;
            A(imat,imat) = A(imat,imat) /relaxU;
        end
    end

    %     aP = reshape(diag(A),ncellx,ncellr)'; % diagnol element of coeff matrix A
    for i = 1:ncellr
        for j = 1:ncellx
            imat = (i-1)*ncellx + j;
            aP(i,j) = A(imat,imat);
        end
    end

    VbyA = volume./(aP); % pesdo conductive coeff of equ of pressure = VP/aP = volvin = volume velocity inverse = V* rAU from Openfoam
%     ux = (reshape(A\b(:,1),ncellx,ncellr))';% cause the inner cycle index is j, so the reshape shuould be in order of ncellx(for j) and then ncellr.
%     ur = (reshape(A\b(:,2),ncellx,ncellr))';
    [u1Sol, u1Flag] = bicgstab(A,b(:,1));
    [u2Sol, u2Flag] = bicgstab(A,b(:,2));
    ux = (reshape(u1Sol,ncellx,ncellr))';
    ur = (reshape(u2Sol,ncellx,ncellr))';

    u(:,:,1) = ux;
    u(:,:,2) = ur;

   
%     disp('momentum equation is solved********')



    Ap = sparse(zeros(ncellx*ncellr,ncellx*ncellr));
    bp = zeros(ncellx*ncellr,1);
    for i = 1:ncellr
        for j = 1:ncellx
            imat = (i-1)*ncellx+j;
            for k = 1:4
                jmat = (i-1+ioffset(k))*ncellx + j+joffset(k);
                if bType(i,j,k) == -1
                    VbyA_f(i,j,k) = (VbyA(i,j) + VbyA(i+ioffset(k),j+joffset(k)))/2.0;% pesdo condc coeff of equ pressure interpolated to face

                    uP = [u(i,j,1),u(i,j,2)]; % predict value
                    uN = [u(i + ioffset(k) ,j + joffset(k),1),u(i + ioffset(k) ,j + joffset(k),2)];
                    gradPP = [gradP(i,j,1),gradP(i,j,2)]; % pressure on cell centroid
                    gradPN = [gradP(i + ioffset(k) ,j + joffset(k),1),gradP(i + ioffset(k) ,j + joffset(k),2)];
                    VbyAP = VbyA(i,j);
                    VbyAN = VbyA(i + ioffset(k) ,j + joffset(k));
                    sn = normal(:,k);  % surface normal
% 
                    HbyA_f(i,j,k) = ((uP + VbyAP*gradPP + uN + VbyAN*gradPN)/2.0)*sn;
%                     HbyA_f(i,j,k) = (HbyA(i,j) + HbyA(i+ioffset(k),j+joffset(k))/2.0;

                    Ap(imat,imat) = Ap(imat,imat) - VbyA_f(i,j,k)*faceArea(i,j,k)/d(k);
                    Ap(imat,jmat) = VbyA_f(i,j,k)*faceArea(i,j,k)/d(k);

                    bp(imat) = bp(imat) + HbyA_f(i,j,k)*faceArea(i,j,k);

                elseif bType(i,j,k) == 0
                    Ap(imat,imat) = Ap(imat,imat) - 0;
                    bp(imat) = bp(imat)+0;
                    % bp(imat) = bp(imat) + VbyA(i,j)*([gradP(i,j,1),gradP(i,j,2)]*normal(:,k))*faceArea(i,j,k);
                    %  VbyA(i,j)*([gradP(i,j,1),gradP(i,j,2)]*normal(:,k))*faceArea(i,j,k)

                elseif bType(i,j,k) == 4
                    Ap(imat,imat) = Ap(imat,imat) - 0;
                    bp(imat) = bp(imat) + 0;
                    % bp(imat) = bp(imat) + VbyA(i,j)*([gradP(i,j,1),gradP(i,j,2)]*normal(:,k))*faceArea(i,j,k);
                    %  VbyA(i,j)*([gradP(i,j,1),gradP(i,j,2)]*normal(:,k))*faceArea(i,j,k)


                else
                    error('Error occurred. \n Undefined boundary!');
                end
            end
        end
    end

    pOld = p;
%     p = (reshape(Ap\bp,ncellx,ncellr))';
    [pSol, pFlag] = bicgstab(Ap,bp,1e-5,100);
    p = (reshape(pSol,ncellx,ncellr))';
    pUnCor = p;
    p = pOld + relaxP *(p-pOld); % Jasak thesis equ 3.145

%     disp('continuity equation is solved********')



    gradPOld = gradP;
    gradP = grad(ncellx,ncellr,dx,dr,p);


    for i=1:ncellr
        for j=1:ncellx
            u(i,j,:)= u(i,j,:)-VbyA(i,j)*(gradP(i,j,:)-gradPOld(i,j,:));
        end
    end

    maxU = max(max((u(:,:,1))));
    if iter >= 2
        res=abs(maxU - maxUOld);
        if mod(iter,10) == 0
            disp([iter,res])

        end
        if res<= 1e-5
            break
        elseif res >=100
            break
        end
    end
    maxUOld = maxU;

end

paint(u(:,:,1))
figure
paint(u(:,:,2))
figure
paint(p)