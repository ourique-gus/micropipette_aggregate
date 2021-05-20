SUBROUTINE Integrate(p, w, N, Np, Nw, Lx, Ly, dx, dy, &
                m_kb, m_kA, m_ks, m_A0, m_r0, m_Fr, m_re, &
                n_ks, n_r0, n_Fr, n_re, n_ms, n_dp, &
                c_Fr, c_Fa, c_re, c_rl, &
                w_Fr,w_re, w_xr, w_xl, w_yb, w_yt, &
                v0, ddt, ddr, tau, mu0, Af, dt, Ndt, Nret, Ncores, pr, wr, &
                fr_mb, fr_mA, fr_ms, fr_mF, fr_ns, fr_nF, &
                fr_NNs, fr_NNF, fr_cF, fr_sAf, fr_wF, &
                nanvalue)
                

USE omp_lib
IMPLICIT NONE
    REAL*8, INTENT(IN) :: p(:, :, :), w(:,:), Lx, Ly, dx, dy
    REAL*8, INTENT(IN) :: m_kb, m_kA, m_ks, m_A0, m_r0, m_Fr, m_re
    REAL*8, INTENT(IN) :: n_ks, n_r0, n_Fr, n_re, n_ms, n_dp
    REAL*8, INTENT(IN) :: c_Fr, c_Fa, c_re, c_rl
    REAL*8, INTENT(IN) :: w_Fr, w_re, w_xr, w_xl, w_yb, w_yt
    REAL*8, INTENT(IN) :: v0, ddt, ddr, tau, mu0, Af(:), dt
    INTEGER, INTENT(IN) :: N, Np, Nw, Ndt, Nret, Ncores
    REAL*8, INTENT(OUT) :: pr(N+1, 5,Np, Nret), wr(Nw,4, Nret)
    REAl*8, INTENT(OUT) :: fr_mb(N, 2, Np, Nret), fr_mA(N, 2,Np, Nret), fr_ms(N, 2, Np, Nret), fr_mF(N,2,Np, Nret)
    REAl*8, INTENT(OUT) :: fr_ns(N, 2, Np, Nret), fr_nF(N, 2,Np, Nret), fr_NNs(2, Np, Nret), fr_NNF(2, Np, Nret)
    REAl*8, INTENT(OUT) :: fr_cF(N,2, Np, Nret)
    REAl*8, INTENT(OUT) :: fr_sAf(N, 2, Np, Nret), fr_wF(N, 2, Np, Nret)
    LOGICAL, INTENT(OUT) :: nanvalue
    REAL*8 :: f(N+1,2,Np), pn(N+1,5,Np), wn(Nw,4)
    REAl*8 :: f_mb(N,2,Np), f_mA(N,2,Np), f_ms(N,2,Np), f_mF(N,2,Np), f_ns(N,2,Np), f_nF(N,2,Np)
    REAl*8 :: f_NNs(2, Np), f_NNF(2, Np)
    REAl*8 :: f_cF(N, 2, Np)
    REAl*8 :: f_sAf(N,2, Np), f_wF(N, 2, Np)
    REAL*8 :: dist(N,4,4), A(N,4), Av(N,4,2), s(N), su(N,2), sd(N,2), Req, Frep, Area, ddif(N,2), adif(N,3), step(N)
    INTEGER :: i, j, ii, jj, kk, ll, Nv, rsp, step_var, step_pos,  f_check(N, Np)
    INTEGER :: Nc(INT(Lx/dx)*INT(Ly/dy)), Ncpos(N,Np), pos(2, N, INT(Lx/dx)*INT(Ly/dy)), ct(9,INT(Lx/dx)*INT(Ly/dy))
    INTEGER :: uu, tt, posv, Nx, Ny, id1, id2, idp1, idp2, nei(2,N, N, Np), nn(N,Np)
    INTEGER :: wNc(INT(Lx/dx)*INT(Ly/dy)), wpos(Nw, INT(Lx/dx)*INT(Ly/dy)), wnei(N,N, Np), wnn(N,Np)
    REAL*8 :: rv(2), r, faux, pi
    REAL*8 :: sqdt, sq3, Acte, otau, rdn(N)
    INTEGER*8 :: b_time, s_time, e_time, r_time, Ntries
    
    CALL SYSTEM_CLOCK(b_time, r_time)
    Ntries=10000
    
    
    pi=2*DACOS(0.d0)
    
    nanvalue=.FALSE.
    
    CALL omp_set_num_threads( Ncores )
    
    sqdt=DSQRT(dt)
    sq3=DSQRT(3.0d0)
    Acte=DSQRT(2*ddr*dt)
    otau=1./tau
    
    pn=p
    wn=w
    
    Nv=N+1

    Nx=INT(Lx/dx)
    Ny=INT(Ly/dy)
    

    ct=-1
    DO ii=1, Nx
        DO jj=1, Ny
            kk=ii+(jj-1)*Ny
            IF ((ii-1 > 0) .AND. (jj-1 > 0)) THEN
                ct(1,kk)=(ii-1)+((jj-1)-1)*Ny
            END IF
            IF (jj-1 > 0) THEN
                ct(2,kk)=(ii)+((jj-1)-1)*Ny
            END IF
            IF (ii+1 <= Nx) THEN
                ct(3,kk)=(ii+1)+((jj-1)-1)*Ny
                ct(6,kk)=(ii+1)+((jj)-1)*Ny
            END IF
            IF (ii-1 > 0) THEN
                ct(4,kk)=(ii-1)+((jj)-1)*Ny
            END IF
            ct(5,kk)=kk
            IF ((ii-1 > 0) .AND. (jj+1 <= Ny)) THEN
                ct(7,kk)=(ii-1)+((jj+1)-1)*Ny
            END IF
            IF (jj+1 <= Ny) THEN
                ct(8,kk)=(ii)+((jj+1)-1)*Ny
            END IF
            IF ((ii+1 <= Nx) .AND. (jj+1 <= Ny)) THEN
                ct(9,kk)=(ii+1)+((jj+1)-1)*Ny
            END IF

        END DO
    END DO
    
    step_pos=0
    step_var=0
    rsp=Ndt/Nret
    
    !! TIME = 0.19
    
    DO tt=1, Ndt
    
        f=0
        f_mb=0
        f_mA=0
        f_ms=0
        f_mF=0
        f_ns=0
        f_nF=0
        f_NNs=0
        f_NNF=0
        f_cF=0
        f_sAf=0
        f_wF=0
        f_check=0

        !$OMP PARALLEL SECTIONS PRIVATE(ii,jj)

        !$OMP SECTION
        !! put_cell_in_box
        Nc=0
        DO ii=1, Np
            DO jj=1,N
                posv=INT(pn(jj, 1,ii)/dx)+1+INT(pn(jj, 2,ii)/dy)*Ny
                Nc(posv)=Nc(posv)+1
                pos(1, Nc(posv), posv)=jj
                pos(2, Nc(posv), posv)=ii
                Ncpos(jj,ii)=Nc(posv)
            END DO
        END DO
        
        !$OMP SECTION
        !! put_wall_in_box
        wNc=0
        DO ii=1,Nw
            posv=INT(wn(ii,1)/dx)+1+INT(wn(ii,2)/dy)*Ny
            wNc(posv)=wNc(posv)+1
            wpos(wNc(posv),posv)=ii
        END DO
        
        !$OMP END PARALLEL SECTIONS
    
    
        !! find_neighbours_to_interact
        nn=0
        wnn=0
        
        !$OMP PARALLEL DO PRIVATE(i,ii,jj,kk, ll, posv) 
        DO ii=1,Np
            DO jj=1,N
                posv=INT(pn(jj,1,ii)/dx)+1+INT(pn(jj,2,ii)/dy)*Ny
                DO kk=1,4
                    DO ll=1, Nc(ct(kk,posv))
                        nn(jj,ii)=nn(jj,ii)+1
                        nei( :,nn(jj,ii), jj, ii)=pos( :, ll, ct(kk,posv) )
                    END DO
                END DO
                DO ll=Ncpos(jj,ii)+1, Nc(posv)
                    nn(jj,ii)=nn(jj,ii)+1
                    nei( :,nn(jj, ii), jj, ii)=pos( :, ll, posv )
                END DO
                DO kk=1,9
                    DO ll=1, wNc(ct(kk,posv))
                            wnn(jj,ii)=wnn(jj,ii)+1
                            wnei(wnn(jj,ii), jj, ii)=wpos( ll, ct(kk,posv) )
                    END DO
                END DO
            END DO
        END DO
        !$OMP END PARALLEL DO
        
        !! calculate_interactions
        !$OMP PARALLEL DO PRIVATE(ii,jj,kk, rv, r, faux) REDUCTION ( + : f_mF, f_cF )
        DO ii=1, Np
            DO jj=1,N
                DO kk=1,nn(jj,ii)
                    id2=nei(1,kk,jj,ii)
                    idp2=nei(2,kk,jj,ii)
                    rv=pn(id2,1:2,idp2)-pn(jj,1:2, ii)
                    r=DSQRT(rv(1)*rv(1)+rv(2)*rv(2))
                    IF (ii==nei(2,kk,jj,ii)) THEN
                        IF (r <= m_re) THEN
                            faux=m_Fr*(r-m_re)/m_re
                            f_mF(jj,:,ii)=f_mF(jj,:,ii)+faux*rv/r
                            f_mF(id2,:,idp2)=f_mF(id2,:,idp2)-faux*rv/r
                        END IF
                    ELSE
                        IF (r < c_re) THEN
                            faux=c_Fr*(r-c_re)/c_re
                            f_cF(jj,:,ii)=f_cF(jj,:,ii)+faux*rv/r
                            f_cF(id2,:,idp2)=f_cF(id2,:,idp2)-faux*rv/r  
                            f_check(jj,ii)=1
                            f_check(id2,idp2)=1  
                        ELSEIF (r < c_rl) THEN
                            faux=c_Fa*(r-c_re)/(c_rl-c_re)  
                            f_cF(jj,:,ii)=f_cF(jj,:,ii)+faux*rv/r
                            f_cF(id2,:,idp2)=f_cF(id2,:,idp2)-faux*rv/r  
                            f_check(jj,ii)=1
                            f_check(id2,idp2)=1
                        END IF
                    END IF 
                END DO
                DO kk=1,wnn(jj,ii)
                    rv=pn(jj,1:2,ii)-wn(wnei(kk,jj,ii),1:2)
                    r=DSQRT(rv(1)*rv(1)+rv(2)*rv(2))
                    IF (r <= w_re) THEN
                        faux=w_Fr*(r-w_re)/w_re
                        f_wF(jj,:,ii)=f_wF(jj,:,ii)-faux*rv/(r+1E-6)
                        f_check(jj,ii)=1
                    END IF
                END DO
            END DO
        END DO
        !$OMP END PARALLEL DO

        !! calculane_membrane
        !$OMP PARALLEL DO PRIVATE(ii, jj, kk, dist, A, Av, s, su, sd, Area, ddif, adif, step)
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
       
        f(1:N,:,:)=f_mb+f_mA+f_ms+f_mF+f_ns+f_nF+f_cF+f_sAf+f_wF
        f(Nv,:,:)=f_NNs+f_NNF
        
        !! execute_movement
		!$OMP PARALLEL DO PRIVATE(kk,ii, rdn) SHARED(otau, v0)
		DO ii=1, Np
			pn(:,3,ii)=mu0*f(:,1,ii)
			pn(:,4,ii)=mu0*f(:,2,ii)
			pn(:,1,ii)=pn(:,1,ii)+pn(:,3,ii)*dt
			pn(:,2,ii)=pn(:,2,ii)+pn(:,4,ii)*dt

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
        
                
        
        IF (ANY(ISNAN(pn))) THEN
            nanvalue=.TRUE.
            EXIT
        END IF

    END DO
    
    !WRITE(*,*) f_mF
    
    
RETURN

END SUBROUTINE Integrate








