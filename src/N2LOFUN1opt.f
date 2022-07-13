c     real*8 cont(2) 
c     cont1(1) means 1s0 cont(2) means 3s1
       subroutine phasetheo 
             use phase
             implicit real*8 (a-h,o-z)
c     we use common conta(2)to pass on this para
           

c     we use phaseshifts to pass out the useful phaseshifts
c     to be specific,it's  
             real*8 elab(41)
             parameter (n=64)
             parameter (n6=6*n,na=2*(n+1)*2*(n+1),naa=6*n/2*(n+1))
             dimension vv(n6),s(n),u(n),a(na),b(na),aa(naa),qq(n),eq(n)
             external n2lo500
             common /alpha/ melab
             common /einject/ elab 
             common /crdwrt/ kread,kwrite,kpunch,kda(9)          
             integer,save:: i=0
             call phases (n2lo500,vv,s,u,a,b,aa,qq,eq)
             i=i+1
             write(*,*) i
            end
