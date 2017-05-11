! Programa que calcula perfiles de densidad promedio e histogramas.
! Calcula el parámetro de orden usando una celda de tamaño fijo y 
! desplazándola un deltap
!
! Programador: Jorge M
! date: 17-03-17

module variables
 implicit none
 integer:: nat,nbloq,nblancos=9,lines=0
 integer:: nbinx,nbiny,nbinz,nrbin,nbinp
 double precision :: dd=0.02,deltap=0.3
 double precision :: lxi,lyi,lzi
 double precision :: lxf,lyf,lzf
 double precision :: lx,ly,lz
 double precision,dimension(:,:),allocatable:: rz,ry,rx
 double precision,dimension(:,:),allocatable:: ez,ey,ex
 double precision,dimension(:),allocatable :: histox,histoy,histoz
 integer,dimension(:),allocatable :: frecx,frecy,frecz
 integer,dimension(:),allocatable :: nqz
 double precision,dimension(:),allocatable :: p2s,p2d 
 double precision :: dirx,diry,dirz
 double precision,dimension(:),allocatable :: qxxp,qxyp,qxzp,qyyp,qyzp,qzzp
end
!***********************
!*********MAIN**********
program densi
use variables
 implicit none

 call leedata
write(*,*) 'salí leedata'
 call posiciones
write(*,*) 'salí posi'
 call p2_nematico
write(*,*) 'salí p2'
! call densidad
! call histo

end
!*********************
!*********************
subroutine leedata
use variables
 implicit none

 integer :: i,ios

 open(1,file='HISTORY.lammpstrj',status='old',action='read')
  do
   read(1,*,iostat=ios)
   if(ios/=0) exit     
   lines=lines+1
  enddo
 close(1)

 open(1,file='HISTORY.lammpstrj',status='old',action='read')
  do i=1,3
   read(1,*)
  enddo
  read(1,*)nat
  read(1,*)
  read(1,*)lxi,lxf
  read(1,*)lyi,lyf
  read(1,*)lzi,lzf
  lx = lxf-lxi
  ly = lyf-lyi
  lz = lzf-lzi
 close(1)

 nbloq = lines/(nat+nblancos)
 allocate(rz(nbloq,nat),ry(nbloq,nat),rx(nbloq,nat))
 allocate(ez(nbloq,nat),ey(nbloq,nat),ex(nbloq,nat))
 
 nbinp = int(lz/deltap)
 allocate(qxxp(0:nbinp),qxyp(0:nbinp),qxzp(0:nbinp),qyyp(0:nbinp),qyzp(0:nbinp),qzzp(0:nbinp))
 allocate(nqz(0:nbinp),p2s(0:nbinp),p2d(0:nbinp))

end
!***********************
subroutine posiciones
use variables
 implicit none

 integer::i,j,id,tipo

 open(1,file='HISTORY.lammpstrj',status='old',action='read')
  do j=1,nbloq
   do i=1,nblancos
    read(1,*)
   enddo
   do i=1,nat
    read(1,*)id,tipo,rx(j,id),ry(j,id),rz(j,id),ex(j,id),ey(j,id),ez(j,id)
   enddo
  enddo
 close(1)
end
!************************
subroutine p2_nematico
use variables
 implicit none

 integer :: i,j,bin,n
 double precision :: qq(3,3),eigval(3),eigvec(3,3)
 double precision :: eigval1,eigval2,eigval3
 double precision :: s,nraiz,igamma
 integer :: ini,fin,celda
 double precision :: qxxb,qxyb,qxzb,qyyb,qyzb,qzzb

n=0
p2d = 0.0d0

do j=1,nbloq
  nqz = 0
  qxxp = 0.0d0
  qxyp = 0.0d0
  qxzp = 0.0d0
  qyyp = 0.0d0
  qyzp = 0.0d0
  qzzp = 0.0d0
 do i=1,nat
  bin = rz(j,i)/deltap
  !Sabemos que algunas partículas están ligeramente fuera de la caja.
  !dado el deltap pueden caer en otro bin. Aqui las consideramos dentro
  !del bin extremo.
  if(bin > nbinp)  bin = nbinp
  nqz(bin)  = nqz(bin) + 1
  qxxp(bin) = qxxp(bin) + ex(j,i)*ex(j,i)
  qxyp(bin) = qxyp(bin) + ex(j,i)*ey(j,i)
  qxzp(bin) = qxzp(bin) + ex(j,i)*ez(j,i)
  qyyp(bin) = qyyp(bin) + ey(j,i)*ey(j,i)
  qyzp(bin) = qyzp(bin) + ey(j,i)*ez(j,i)
  qzzp(bin) = qzzp(bin) + ez(j,i)*ez(j,i)
 enddo
 
 ini   = 0 
 fin   = 10 !numero de bins que conforman una celda
 do i=ini, fin
  n = n + nqz(i)
  qxxb = qxxb + qxxp(i)
  qxyb = qxyb + qxyp(i)
  qxzb = qxzb + qxzp(i)
  qyyb = qyyb + qyyp(i)
  qyzb = qyzb + qyzp(i)
  qzzb = qzzb + qzzp(i)
 enddo

 celda = 0
 do while(ini<nbinp)
  if(n>0) then
   qq(1,1) = 1.50*qxxb/dble(n)-0.50
   qq(1,2) = qxyb/dble(n)
   qq(1,3) = qxzb/dble(n)
   qq(2,1) = qq(1,2)
   qq(2,2) = 1.50*qyyb/dble(n)-0.50
   qq(2,3) = qyzb/dble(n)
   qq(3,1) = qq(1,3)
   qq(3,2) = qq(2,3)
   qq(3,3) = 1.50*qzzb/dble(n)-0.50
   call jacobi(3,qq,eigval,eigvec)
   eigval1 = eigval(1)
   eigval2 = eigval(2)
   eigval3 = eigval(3)

   p2d(celda) = p2d(celda) + eigval3

   ini   = ini + 1
   fin   = fin + 1
   celda = celda + 1
   if(fin>nbinp)  fin  = 0
   n     = n + nqz(fin) - nqz(ini-1)
   qxxb  = qxxb + qxxp(fin) - qxxp(ini-1)
   qxyb  = qxyb + qxyp(fin) - qxyp(ini-1)
   qxzb  = qxzb + qxzp(fin) - qxzp(ini-1)
   qyyb  = qyyb + qyyp(fin) - qyyp(ini-1)
   qyzb  = qyzb + qyzp(fin) - qyzp(ini-1)
   qzzb  = qzzb + qzzp(fin) - qzzp(ini-1)
  endif
 enddo

enddo

p2d = p2d/dble(nbloq)
open(21,file='p2d.dat',status='replace') 
 do i=0,celda-1
  s = dble(i)*deltap
  write(21,*) s,p2d(i)
 enddo
close(21)

end
!************************
subroutine densidad
use variables
 implicit none

 integer:: i,j,binx,biny,binz

 nbinx = int(lx/dd)
 nbiny = int(ly/dd)
 nbinz = int(lz/dd)

 allocate(histox(0:nbinx))
 allocate(histoy(0:nbiny))
 allocate(histoz(0:nbinz))

 histox=0.0
 histoy=0.0
 histoz=0.0

 do j=1,nbloq
  do i=1,nat

   binx = int(rx(j,i)/dd)
   biny = int(ry(j,i)/dd)
   binz = int(rz(j,i)/dd)
   
   if(binx>0 .and. binx<nbinx) histox(binx)=histox(binx)+1.0
   if(biny>0 .and. biny<nbiny) histoy(biny)=histoy(biny)+1.0
   if(binz>0 .and. binz<nbinz) histoz(binz)=histoz(binz)+1.0

  enddo
 enddo

 histox = histox/(ly*lz*dd*nbloq)
 histoy = histoy/(lx*lz*dd*nbloq)
 histoz = histoz/(lx*ly*dd*nbloq)

 open(3,file='dens.dat',status='replace',action='write')
  do i=0,nbinz
   write(3,001)real(i*dd),histoz(i)
  enddo
  write(3,*)'&'
  do i=0,nbiny
   write(3,001)real(i*dd),histoy(i)
  enddo
  write(3,*)'&'
  do i=0,nbinx
   write(3,001)real(i*dd),histox(i)
  enddo
 close(3)

001 format (f5.2,1X,f7.4)

end
!***********************
subroutine histo
use variables
 implicit none

 double precision :: drho=0.001
 integer :: i
 integer,dimension(:),allocatable:: hbinz,hbinx,hbiny 

 allocate(hbinz(nbinz),hbinx(nbinx),hbiny(nbiny))

 hbinz = int(histoz/drho)
 hbinx = int(histox/drho)
 hbiny = int(histoy/drho)
 
 do i=0,nbinz 
  write(77,*) hbinz(i)
 enddo 
 do i=0,nbinx
  write(78,*) hbinx(i)
 enddo 
 do i=0,nbiny
  write(79,*) hbiny(i)
 enddo

 call system('sort -n -o fort.77 fort.77 && uniq -c fort.77 > outz')
 call system('sort -n -o fort.78 fort.78 && uniq -c fort.78 > outx')
 call system('sort -n -o fort.79 fort.79 && uniq -c fort.79 > outy')

 call system("awk '{print $2*0.001,$1}' < outz > densz.dat")
 call system("awk '{print $2*0.001,$1}' < outx > densx.dat")
 call system("awk '{print $2*0.001,$1}' < outy > densy.dat")

 call system('rm fort* out*')

end 
!**************************
subroutine jacobi(np,a,d,v)
implicit none
integer          :: i,j,k,ip,iq,n,np,nrot,maxrot
double precision :: sm,tresh,ss,c,t,theta,tau,h,gg,p
double precision ::  a(np,np),d(np),v(np,np),b(np),z(np)
!
!     setup and initialization
!
      maxrot = 100
      nrot = 0
      n = np
      do ip = 1, n
         do iq = 1, n
            v(ip,iq) = 0.0d0
         end do
         v(ip,ip) = 1.0d0
      end do
      do ip = 1, n
         b(ip) = a(ip,ip)
         d(ip) = b(ip)
         z(ip) = 0.0d0
      end do
!
!     perform the jacobi rotations
!
      do i = 1, maxrot
         sm = 0.0d0
         do ip = 1, n-1
            do iq = ip+1, n
               sm = sm + abs(a(ip,iq))
            end do
         end do
         if (sm .eq. 0.0d0)  goto 10
         if (i .lt. 4) then
            tresh = 0.2d0*sm / n**2
         else
            tresh = 0.0d0
         end if
         do ip = 1, n-1
            do iq = ip+1, n
               gg = 100.0d0 * abs(a(ip,iq))
               if (i.gt.4 .and. abs(d(ip))+gg.eq.abs(d(ip)).and.abs(d(iq))+gg.eq.abs(d(iq))) then
                  a(ip,iq) = 0.0d0
               else if (abs(a(ip,iq)) .gt. tresh) then
                  h = d(iq) - d(ip)
                  if (abs(h)+gg .eq. abs(h)) then
                     t = a(ip,iq) / h
                  else
                     theta = 0.5d0*h / a(ip,iq)
                     t = 1.0d0 / (abs(theta)+sqrt(1.0d0+theta**2))
                     if (theta .lt. 0.0d0)  t = -t
                  end if
                  c  = 1.0d0 / sqrt(1.0d0+t**2)
                  ss = t*c
                  tau = ss/(1.0d0+c)
                  h = t * a(ip,iq)
                  z(ip) = z(ip) - h
                  z(iq) = z(iq) + h
                  d(ip) = d(ip) - h
                  d(iq) = d(iq) + h
                  a(ip,iq) = 0.0d0
                  do j = 1, ip-1
                     gg = a(j,ip)
                     h = a(j,iq)
                     a(j,ip) = gg - ss*(h+gg*tau)
                     a(j,iq) = h + ss*(gg-h*tau)
                  end do
                  do j = ip+1, iq-1
                     gg = a(ip,j)
                     h = a(j,iq)
                     a(ip,j) = gg - ss*(h+gg*tau)
                     a(j,iq) = h + ss*(gg-h*tau)
                  end do
                  do j = 1, n
                     gg = v(j,ip)
                     h = v(j,iq)
                     v(j,ip) = gg - ss*(h+gg*tau)
                     v(j,iq) = h + ss*(gg-h*tau)
                  end do
                  nrot = nrot + 1
               end if
            end do
         end do
         do ip = 1, n
            b(ip) = b(ip) + z(ip)
            d(ip) = b(ip)
            z(ip) = 0.0d0
         end do
      end do
!
!     print warning if not converged
!
   10 continue
      if (nrot .eq. maxrot) then
         write (*,20)
   20    format (/,' JACOBI  --  Matrix Diagonalization not Converged')
      end if
!
!     sort the eigenvalues and vectors
!
      do i = 1, n-1
         k = i
         p = d(i)
         do j = i+1, n
            if (d(j) .lt. p) then
               k = j
               p = d(j)
            end if
         end do
         if (k .ne. i) then
            d(k) = d(i)
            d(i) = p
            do j = 1, n
               p = v(j,i)
               v(j,i) = v(j,k)
               v(j,k) = p
            end do
         end if
      end do

end
