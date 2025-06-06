	program ising


	real*8 DE, p, dran_u,xi
	real*8 T,h                   !Paso entre temperaturas: h
        real*8 E, Ec, E_m, Ec_m      !Suma E(S)/N^2  (promediado en spines)
        real*8 M, Mc, M_m, Mc_m      !Suma spines/N^2  (promediado en spines)

        real*8 cesp,chi
        
	integer i,j,N,k,q
        integer i_dran,l,paso,p_medida,i_temp
        
	parameter(N=64)     !DEBE SER PAR PARA CALCULAR BIEN funci¢n de correlaci¢n
	real*8 S(0:N+1,0:N+1)
	real*8 f(N)  !Funci¢n de correlaci¢n
	
	h=0.3d0  !Por defecto paso entre temperatura h=0.3 pero cerca de Tc veremos que lo reducimos

	open(12,file='energia_udvol_cada100pMC.dat')
	open(17,file='magnetizacion_udvol_cada100pMC.dat')
	
	open(13,file='Em_T.dat')
	open(18,file='Mm_T.dat')
	
	open(14,file='caloresp16.dat')
	open(15,file='susceptibilidad16.dat')
	open(16,file='f.dat')


	open(19,file='lc.dat')

	!C lculo de exponentes criticos
	open(20,file='exp_m.dat')
	open(21,file='exp_c.dat')
	open(22,file='exp_chi.dat')

	call dran_ini(761389)     !1897263


	!Condiciones iniciales
	do i=1,N
	  do j=1,N
	    S(i,j)=1.0d0
	  end do
	end do

	!Condiciones de contorno periodicas
	do i=1,N
	   s(0,i)=s(N,i)
	   s(N+1,i)=s(1,i)
	   s(i,0)=s(i,N)
	   s(i,N+1)=s(i,1)
	end do

	!Tomamos datos cada 100 pasos montecarlo:
	p_medida=10000    !N£mero de veces que se hacen los 100pMC  (10000 veces 100pMC)
	paso=N*N*100  !Pasos totales en 100pMC
	
!!!!!!!Bucle para c lculos para distintas Temperaturas:
      T=1.5d0
           !Caso 1: 1 sola Temperatura
      !do i_temp=1,1
           !Caso 2: varias (11) Temperaturas
      do i_temp=1,11
         write(17,*) 'T=', T
         write(17,*)
         write(12,*) 'T=', T
         write(12,*)
         write(16,*) 'T=', T
         write(16,*)

	!Inicializo:
	M_m=0.0d0      !Magnetizaci¢n media
	Mc_m=0.0d0  !M^2 media

	E_m=0.0d0      !Energ¡a total promediada: <E>
	Ec_m=0.0d0     !Energ¡a total al cuadrado promediada: <E^2>
	f=0.0d0        !Funci¢n de correlaci¢n

!!!!!!ALGORITMO:
      do l=1,p_medida
	do k=1,paso
	  !Spin aleatorio:
		i=i_dran(N)
		j=i_dran(N)

          !Condiciones peri¢dicas:
	     if (i.eq.1) then
	       S(i-1,j)=S(N,j)
	     else if (i.eq.N) then
	       S(i+1,j)=S(1,j)
	     end if

	     if (j.eq.1) then
	       S(i,j-1)=S(i,N)
	     else if (j.eq.N) then
               S(i,j+1)=S(i,1)
	     end if

          !Calculo DeltaE y hago p=min(1,DE)
	     DE=2.0*S(i,j)*(S(i+1,j)+S(i-1,j)+S(i,j-1)+S(i,j+1))

	     p=min(exp(-DE/T),1.0)
          !Genero n£mero aleatorio uniforme en [0,1] y cambio el signo del spin seg£n si es mayor o menor que p
           xi=dran_u()
           if (xi.lt.p) then
             S(i,j)=-S(i,j)
           endif
        enddo                    !Fin paso k
        !_____________________________________________________________

	M=0.0
	E=0.0

!!!!!!!!!!!Calculo M, E, f  cada 100pMC
	do i=1,N
	  do j=1,N
	    M=M+S(i,j)
	    E=E-0.5d0*S(i,j)*(S(i,j+1)+S(i,j-1)+S(i+1,j)+S(i-1,j))
          end do
	end do
              !! M: Mag total
        if (mod(l,100).eq.0) then
	      write(17,*) l,abs(M)/(dfloat(N)**2.0d0)  !Magnetizaci¢n por unidad de volumen: m
	endif
	
	M_m=M_m+abs(M) !<M> (luego hay que dividir, promediar en cada 100pMC)
	
	Mc=M**2.0d0  !M^2
	Mc_m=Mc_m+Mc  !<M^2>

              !! E : Energ¡a total (no es por unidad de volumen)
        if (mod(l,100).eq.0) then
	   write(12,*) l,E/(dfloat(N)**2.0d0) !Energ¡a por ud volumen
	endif
	
	E_m=E_m+E   !<E>
	
	Ec=E**2.0
	Ec_m=Ec_m+Ec
	
	!Calculo funci¢n correlaci¢n. NOTA: N debe ser PAR
	do q=1,N/2   !Debido a las cond peri¢dicas, dos espines en l¡nea no pueden estar a una distancia mayor que N/2 (si ponemos N estar¡amos sobrecontando).
                             ! M x distancia: (N/2,N/2), pero aqu¡ solo la calculamos en 1 dimensi¢n
	  do i=1,N
	    do j=1,N
               if (i+q.gt.N) then
	         f(q)=f(q)+(S(i,j)*S(i+q-N,j))
	       else
	         f(q)=f(q)+(S(i,j)*S(i+q,j))
	       endif
            end do
	  end do
	  
	end do

      end do                     !Fin algoritmo
      !_________________________________________________________________
      
         write(17,*)
         write(12,*)

      !Funci¢n de correlaci¢n media (en cada 100pMC), por ud volumen
         do q=1,N/2 !NOTA: N debe ser PAR
	    write(16,*) dfloat(q), f(q)/(p_medida*dfloat(N)**2.d0)
         end do
         write(16,*)

      !E y M medias (en cada 100mPC), por ud volumen
	E_m=E_m/p_medida   !Energ¡a media total
	Ec_m=Ec_m/p_medida
	
	M_m=M_m/p_medida     !Magnetizaci¢n media total
	Mc_m=Mc_m/p_medida
	
	write(13,*) T, E_m/(dfloat(N)**2.0d0)
	write(18,*) T, M_m/(dfloat(N)**2.0d0)

       	!Calor espec¡fico
	cesp=(Ec_m-E_m**2.0d0)/(T*dfloat(N))**2.0d0
	write(14,*) T, cesp

        !Susceptibilidad magn‚tica (chi)
	chi=(Mc_m-M_m**2.0d0)/(T*dfloat(N)**2.0d0)
	write(15,*) T, chi


        !Calculo exp cr¡ticos: tomo datos s¢lo de T<Tc
        !Esto es porque puede que el exp cr¡t sea no entero y no se puede elevar un n£mero negativo a este exponente.
	if (T.lt.2.269d0) then
		write(20,*) T,M_m
		write(21,*) T,cesp
		write(22,*) T,chi
	end if
	
	
       !Cerca de temperatura cr¡tica (Tc=2.269) hay m s cambio de la funci¢n: necesitamos m s datos --> h=0.1  paso temperatura
	if ((T.ge.(2.0d0)).and.(T.lt.(2.5d0))) then
		h=0.1
	else
		h=0.3
	end if

	T=T+h

      end do

      close(12)
      close(13)
      close(14)
      close(15)
      close(16)
      close(17)
      close(18)
      close(19)
      close(20)
      close(21)
      close(22)

      stop
      end

      include 'dranxor2_new.f'
