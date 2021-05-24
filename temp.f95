!! calculane_membrane
        !$OMP PARALLEL DO PRIVATE(jj, kk, dist, A, Av, s, su, sd, Area, ddif, adif, step)
        DO ii=1,Np
    
            dist(1:N-1,3,1:2)=pn(2:N,1:2, ii)-pn(1:N-1,1:2, ii) ! (i+1) - i
            dist(N,3,1:2)    =pn(1,1:2, ii)  -pn(N,1:2, ii)
    
    
            dist(:,3,3)=dist(:,3,1)*dist(:,3,1)+dist(:,3,2)*dist(:,3,2)
    
            dist(1:N-1,4,1:3) = dist(2:N,3,1:3) ! (i+2) - (i+1)
            dist(N,4,1:3)     = dist(1,3,1:3)
    
            dist(2:N,2,1:3) = dist(1:N-1,3,1:3) ! (i) - (i-1)
            dist(1,2,1:3)   = dist(N,3,1:3)
    
            dist(3:N,1,1:3) = dist(1:N-2,3,1:3) ! (i-1) - (i-2)
            dist(2,1,1:3)   = dist(N,3,1:3)
            dist(1,1,1:3)   = dist(N-1,3,1:3)
    
    
            dist(:,:,4)=SQRT(dist(:,:,3))
    

            
            DO jj=1,N
                A(jj,1)= +1./dist(jj,2,4)
                
                A(jj,2)= -(dist(jj,2,1)*dist(jj,1,1)+dist(jj,2,2)*dist(jj,1,2))/(dist(jj,2,3)*dist(jj,1,4)) &
                     -(dist(jj,3,1)*dist(jj,2,1)+dist(jj,3,2)*dist(jj,2,2))/(dist(jj,3,4)*dist(jj,2,3)) &
                     -1./dist(jj,3,4)
            
                A(jj,3)= +(dist(jj,3,1)*dist(jj,2,1)+dist(jj,3,2)*dist(jj,2,2))/(dist(jj,3,3)*dist(jj,2,4)) &
                     +(dist(jj,4,1)*dist(jj,3,1)+dist(jj,4,2)*dist(jj,3,2))/(dist(jj,4,4)*dist(jj,3,3)) &
                     +1./dist(jj,2,4)
                     
                A(jj,4)= -1./dist(jj,3,4)
        
            END DO
        
            
            Av(:,1,1)=A(:,1)*dist(:,1,1)/dist(:,1,4)
            Av(:,1,2)=A(:,1)*dist(:,1,2)/dist(:,1,4)
    
            Av(:,2,1)=A(:,2)*dist(:,2,1)/dist(:,2,4)
            Av(:,2,2)=A(:,2)*dist(:,2,2)/dist(:,2,4)
    
            Av(:,3,1)=A(:,3)*dist(:,3,1)/dist(:,3,4)
            Av(:,3,2)=A(:,3)*dist(:,3,2)/dist(:,3,4)
    
            Av(:,4,1)=A(:,4)*dist(:,4,1)/dist(:,4,4)
            Av(:,4,2)=A(:,4)*dist(:,4,2)/dist(:,4,4)
            
            s=(1.0-m_r0/dist(:,3,4))
    
            su(:,1)=s*dist(:,3,1)
            su(:,2)=s*dist(:,3,2)
    
            sd(2:N,:)=-su(1:N-1,:)
            sd(1,:)  =-su(N,:)
            
            Area=0.5*DABS(SUM(    +(pn(2:N,2, ii)+pn(1:N-1,2, ii))*(pn(2:N,1, ii)-pn(1:N-1,1, ii)) ) &
                        +(pn(1,2, ii)+pn(N,2, ii))*(pn(1,1, ii)-pn(N,1, ii)) )
                
            

            ddif(2:N-1,1)=pn(3:N,2, ii)-pn(1:N-2,2, ii)
            ddif(1,1)=pn(2,2, ii)-pn(N,2, ii)
            ddif(N,1)=pn(1,2, ii)-pn(N-1,2, ii)
            
            ddif(2:N-1,2)=pn(1:N-2,1, ii)-pn(3:N,1, ii)
            ddif(1,2)=pn(N,1, ii)-pn(2,1, ii)
            ddif(N,2)=pn(N-1,1, ii)-pn(1,1, ii)
            
            
            adif(:,1:2)=ddif(:,:)
            adif(:,3)=SQRT(adif(:,1)**2+adif(:,2)**2)
            

            step=INT(1./(INT(f_check(:,ii)) +1 ))
                
            adif(:,1)=adif(:,1)/adif(:,3)*step
            adif(:,2)=adif(:,2)/adif(:,3)*step
            
            f_mb(:,:,ii)=0.5*m_kb*(Av(:,1,:)+Av(:,2,:)+Av(:,3,:)+Av(:,4,:))
            f_mA(:,:,ii)=-0.5*m_kA*(Area-m_A0)*ddif ! Esse 0.5 vem da derivada da area
            f_ms(:,:,ii)=m_ks*(su+sd)
            
            f_sAf(:,:,ii)=-Af(ii)*adif(:,1:2)
            
            dist(:,3,1)=pn(1:N,1,ii)-pn(Nv,1,ii)
            dist(:,3,2)=pn(1:N,2,ii)-pn(Nv,2,ii)
            dist(:,3,4)=DSQRT(dist(:,3,1)*dist(:,3,1)+dist(:,3,2)*dist(:,3,2))
    
            s=(1.0-n_r0/dist(:,3,4))

            f_ns(:,1,ii)=-n_ks*s*dist(:,3,1)
            f_ns(:,2,ii)=-n_ks*s*dist(:,3,2)
            
            f_nF(:,1,ii)=-((-DSIGN(0.5d0,dist(:,3,4)-n_re)+0.5)*n_Fr*(dist(:,3,4)-n_re)/n_re)*dist(:,3,1)/dist(:,3,4)
            f_nF(:,2,ii)=-((-DSIGN(0.5d0,dist(:,3,4)-n_re)+0.5)*n_Fr*(dist(:,3,4)-n_re)/n_re)*dist(:,3,2)/dist(:,3,4)
    
            f_NNs(1,ii)=-SUM(f_ns(:,1,ii))
            f_NNs(2,ii)=-SUM(f_ns(:,2,ii))
            
            f_NNF(1,ii)=-SUM(f_nF(:,1,ii))
            f_NNF(2,ii)=-SUM(f_nF(:,2,ii))
            
        END DO
        !$OMP END PARALLEL DO
        
        !$OMP PARALLEL DO
		DO ii=1, Np
            f(1:N,:,ii)=f_mb(:,:,ii)+f_mA(:,:,ii)+f_ms(:,:,ii)+f_mF(:,:,ii)+f_ns(:,:,ii) &
                        +f_nF(:,:,ii)+f_cF(:,:,ii)+f_sAf(:,:,ii)+f_wF(:,:,ii)
            f(Nv,:,ii)=f_NNs(:,ii)+f_NNF(:,ii)
        END DO
        !$OMP END PARALLEL DO
        !! execute_movement
		!$OMP PARALLEL DO
		DO ii=1, Np
			pn(:,4,ii)=mu0*f(:,2,ii)
			pn(:,3,ii)=mu0*f(:,1,ii)
			pn(:,2,ii)=pn(:,2,ii)+pn(:,4,ii)*dt
			pn(:,1,ii)=pn(:,1,ii)+pn(:,3,ii)*dt
		END DO
		!$OMP END PARALLEL DO
		
		wn(:,1)=wn(:,1)+wn(:,3)*dt
        wn(:,2)=wn(:,2)+wn(:,4)*dt
        
        step_var=step_var+1
        IF (step_var == rsp) THEN
            step_var=0
            step_pos=step_pos+1
            pr(:,:,:,step_pos)=pn
            wr(:,:,step_pos)=wn
            fr_mb(:,:,:,step_pos)=f_mb
            fr_mA(:,:,:,step_pos)=f_mA
            fr_ms(:,:,:,step_pos)=f_ms
            fr_mF(:,:,:,step_pos)=f_mF
            fr_ns(:,:,:,step_pos)=f_ns
            fr_nF(:,:,:,step_pos)=f_nF
            fr_NNs(:,:,step_pos)=f_NNs
            fr_NNF(:,:,step_pos)=f_NNF
            fr_cF(:,:,:,step_pos)=f_cF
            fr_sAf(:,:,:,step_pos)=f_sAf
            fr_wF(:,:,:,step_pos)=f_wF
            
            IF(step_pos == Nret) THEN
                EXIT
            END IF
        END IF
