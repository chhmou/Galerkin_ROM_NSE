% Generate ROM mass matrix
% this code generate trilinear terms in R^m 
% this code consider the efficient way to generate the tril
% i.e. TriD is of size mxmxm, r+1<=m<= d

clear temp


tic
% MassROM = zeros(N,N);
% StiffROM = zeros(N,N);
% GradDivROM = zeros(N,N);
% % create mass, graddiv and stiffness matrices
% for i=1:N
%     for j=1:N
%         MassROM(i,j) = PhiR(:,i)' * (MassMatrix * PhiR(:,j) );
%         StiffROM(i,j) = PhiR(:,i)' * (StiffnessMatrix * PhiR(:,j) );
%         GradDivROM(i,j) = PhiR(:,i)' * (GradDivMatrix * PhiR(:,j) );
%     end
% end
DdcModes=Dim+1; 

% MassROM = PodU2d'*MassMatrix*PodU2d;
% StiffROM = PodU2d'*StiffnessMatrix*PodU2d;
% GradDivROM = PodU2d'*GradDivMatrix*PodU2d;




PodU2All = zeros(length(Uave),DdcModes);
PodU2All(:,2:DdcModes) = PhiR(:,1:DdcModes-1);
PodU2All(:,1) = Uave;


%
% create 3D matrix for the b(u,v,w) term
%== Comment: TriLinROM = TriLinROM2
TriLinRomFull = zeros(DdcModes,DdcModes,DdcModes);
TriLinRomFull2 = zeros(DdcModes,DdcModes,DdcModes);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:r
    %i;
    for k=1:DdcModes
        %  k;
        tic
        for j=1:DdcModes
            
            
            for itrg=1:size(elnode,1)
                % Description of triangle.
                triag_no = itrg ;
                cotri(1:3,1) = nodeco(elnode(triag_no, 1:3), 1) ;
                cotri(1:3,2) = nodeco(elnode(triag_no, 1:3), 2) ;
                Jmat = [(cotri(2,1) - cotri(1,1)), (cotri(3,1) - cotri(1,1)) ; ...
                    (cotri(2,2) - cotri(1,2)) , (cotri(3,2) - cotri(1,2)) ] ;
                detJ = abs(Jmat(1,1)*Jmat(2,2) - Jmat(1,2)*Jmat(2,1));
     
               % JInv = inv(Jmat) ;
             
                
                [quad_pts, quad_wghts] = feval('quad_75') ;
                
                nqpts = size(quad_pts,1) ;
                quad_wghts = detJ * quad_wghts ;
                
                basvals = basisValues(:,:,itrg);
                Gradtrue1 = gradBasisValues(:,:,:,itrg);
                
                Gradx(:,:) = Gradtrue1(:,1,:) ;
                Grady(:,:) = Gradtrue1(:,2,:) ;
                
                Vstart = [2*(elnode(triag_no,1:3) - 1) + 1 , 2*(elnode(triag_no,4:6) + nVert - 1) + 1 ] ;
                Vel1i = PodU2All([Vstart], i) ;
                Vel1k = PodU2All([Vstart], k) ;
                Vel1j = PodU2All([Vstart], j) ;
                
                Vstart = Vstart + 1 ;
                Vel2i = PodU2All([Vstart], i) ;
                Vel2k = PodU2All([Vstart], k) ;
                Vel2j = PodU2All([Vstart], j) ;
                
                vi1 = Vel1i.' * basvals ;
                vi2 = Vel2i.' * basvals ;
                vj1 = Vel1j.' * basvals ;
                vj2 = Vel2j.' * basvals ;
                vk1 = Vel1k.' * basvals ;
                vk2 = Vel2k.' * basvals ;
                
                vj1x = Vel1j.' * Gradx;
                vj1y = Vel1j.' * Grady;
                vj2x = Vel2j.' * Gradx;
                vj2y = Vel2j.' * Grady;
                
                vi1x = Vel1i.' * Gradx;
                vi1y = Vel1i.' * Grady;
                vi2x = Vel2i.' * Gradx;
                vi2y = Vel2i.' * Grady;
                
                temp = vi1 .* ( vk1 .* vj1x + vk2 .* vj1y) + vi2 .* (vk1 .* vj2x + vk2 .* vj2y);
                TriLinRomFull(i,j,k)=TriLinRomFull(i,j,k)+ quad_wghts * temp.';
                temp2 = vj1 .* ( vk1 .* vi1x + vk2 .* vi1y) + vj2 .* (vk1 .* vi2x + vk2 .* vi2y);
                TriLinRomFull2(i,j,k)=TriLinRomFull2(i,j,k)+ quad_wghts * (0.5 * temp - 0.5 * temp2 ).';
                
                
            end
        end
        toc
    end
    
end

TrD =TriLinRomFull(2:r+1,2:DdcModes,2:DdcModes);

% save the data to file
%save('../ROM_DATA/Re500/Mesh25KTriDim12.mat','TrD','Dim')

toc
