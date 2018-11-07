%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo of "Fuzzy region-based active contour driven by global and local
%  fitting energy for image segmentation" submitting to Applied Soft
%  Computing 
% Jiangxiong Fang(fangchj2002@163.com)
% East China University of Technology
% 6th, August, 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u,deltaF]= GLFIF(Img,LImg,u0,Ksigma,lambda1,lambda2,alpha1,alpha2,g)

  u1 = u0.^2;
  u2 = (1-u0).^2;
         
  Iu1 = Img.*u1;
  Iu2 = Img.*u2;  
  
  c1 = sum(sum(Iu1))/sum(sum(u1));
  c2 = sum(sum(Iu2))/sum(sum(u2));  
  
  k1 = Iu1./u1;
  k2 = Iu2./u2;  
  
   Ku1 = imfilter(u1,Ksigma,'replicate'); 
   Ku2 = imfilter(u2,Ksigma,'replicate'); 
         
   KI1 = imfilter(Iu1,Ksigma,'replicate');
   KI2 = imfilter(Iu2,Ksigma,'replicate');  
         
   s1 = KI1./Ku1;
   s2 = KI2./Ku2; 
   
   kim = c1.*u1+c2.*u2;
   DcH = (LImg-kim).*LImg;
   F3_old = DcH;
   
   sim = s1.*u1+s2.*u2;
   DsH = (LImg-sim).*LImg;
   F4_old = DsH;
   
   un= 1./(1+(lambda1*((Img-c1).^2)+(alpha1*c2+alpha2*s2))./(lambda2*((Img-c2).^2)+(alpha1*c1+alpha2*s1)));
   
   
   un1 = un.^2;
   un2 = (1-un).^2;
   
   delta_u1 = un1-u1;
   delta_u2 = un2-u2;
   
   delta_F1 = un1.*delta_u1.*((Img-c1).^2)./(un1+delta_u1);
   delta_F2 = un2.*delta_u2.*((Img-c2).^2)./(un2+delta_u2);
   
   NIu1 = LImg.*un1;
   NIu2 = LImg.*un2;
   
   Nc1 = sum(sum(NIu1))/sum(sum(un1));
   Nc2 = sum(sum(NIu1))/sum(sum(un1));  
   
   Nk1 = NIu1./un1;
   Nk2 = NIu2./un2; 
   
   Nkim =  un1.*Nc1+un2.*Nc2; 
   F3_new = (LImg-Nkim).*LImg;

   
   
    NKu1 = imfilter(un1,Ksigma,'replicate'); 
    NKu2 = imfilter(un2,Ksigma,'replicate'); 
         
    NKI1 = imfilter(NIu1,Ksigma,'replicate');
    NKI2 = imfilter(NIu2,Ksigma,'replicate');  
         
    Ns1 = NKI1./NKu1;
    Ns2 = NKI2./NKu2; 
    
    
    Nsim =  un1.*Ns1+un2.*Ns2; 
    F4_new = (LImg-Nsim).*LImg;
    
    
    deltaF = lambda1*delta_F1.*g+lambda2*delta_F2.*g+alpha1*(F3_new-F3_old).*g+alpha2*(F4_new-F4_old).*g;
    idx = find(deltaF<0);
    u0(idx)=un(idx); 
    deltaF = sum(sum(deltaF));  
    u = u0;   
    u = imfilter(u,Ksigma,'replicate');  
 
end



