c
c   subroutine wig2(l1,is,l2,a,id1)
c
c   calculates:
c
c
c   a(m1+l1+1,m2+l2+1)= (-1)**m1  ( l1    is    l2 )
c                                 (-m1   m1-m2  m2 )
c
c   for all m1 (-l1.le.m1.le.l1)
c   and all m2 (-l2.le.m2.le.l2)
c
c   program and algorithm by j. h. woodhouse
c
      subroutine wig2(l1,is,l2,a,id1)
      implicit double precision (a-h,o-z)
      double precision a
      dimension a(id1,*) !second index was 1, APV change
      dfloat(n)=n
c
c        zero left hand side of matrix
c
      l12p1=2*l1+1
      l2p1=l2+1
      do 510 im1=1,l12p1
      do 510 im2=1,l2p1
510   a(im1,im2)=0.d0
      if((is.lt.iabs(l1-l2)).or.(is.gt.(l1+l2))) goto 610
c
c        compute the corner value
c
      r1=1.d0/dsqrt(dfloat(2*max0(l1,l2)+1))
      ldel=iabs(l1-l2)
      num=is-ldel
      if(num.eq.0) goto 100
      do 40 n=1,num
      isc=ldel+n
40    r1=r1*dsqrt(dfloat(l1-isc+l2+1)/dfloat(l1+l2+isc+1))
100   a(1,1)=r1
c
c        compute first row
c
      num1=min0(is-l1+l2,l2)
      if(num1.eq.0) goto 120
      do 110 n=1,num1
      m2=-l2+n
110   a(1,n+1)=-a(1,n)*dsqrt(dfloat((l1+is+m2)*(is-l1-m2+1))
     1   /dfloat((l2-m2+1)*(l2+m2)))
120   continue
c
c       compute first column
c
      num2=is-l2+l1
      if(num2.eq.0) goto 130
      do 140 n=1,num2
      m1=-l1+n
140   a(n+1,1)=-a(n,1)*dsqrt(dfloat((l2+is+m1)*(-l2+is-m1+1))
     1   /dfloat((l1+m1)*(l1-m1+1)))
130   continue
      iss=is*(is+1)-l1*(l1+1)-l2*(l2+1)
c
c       iterate south-east
c
      numd=min0(2*is+1,is+l1)
      if(numd.eq.0) goto 610
      do 300 nd=1,numd
      m1b=max0(-l1,is-l2-nd+1)
      m2b=max0(-l2,-is-l1+nd-1)
      numit=min0(-m2b,l1-m1b)
      if(numit.eq.0) goto 300
      do 350 nit=1,numit
      m1=m1b+nit
      m2=m2b+nit
      im1=m1+l1+1
      im2=m2+l2+1
      if((im1.lt.3).or.(im2.lt.3)) goto 320
      r1=a(im1-2,im2-2)
      goto 310
320   r1=0.d0
310   a(im1,im2)=-r1*dsqrt(dfloat((l1+m1-1)*(l1-m1+2)
     1    *(l2-m2+2)*(l2+m2-1)))
      a(im1,im2)=a(im1,im2)+dfloat(iss+2*(m1-1)*(m2-1))
     1    *a(im1-1,im2-1)
      a(im1,im2)=a(im1,im2)/dsqrt(dfloat((l1+m1)*(l1-m1+1)
     1    *(l2-m2+1)*(l2+m2)))
350   continue
300   continue
c
c       fill the rest of array
c
610   sgn=1.d0
      it=l1+l2+is
      if((it/2)*2.ne.it) sgn=-1.d0
      do 420 im1=1,l12p1
      if(l2.eq.0) goto 420
      do 410 m2=1,l2
      m1=im1-l1-1
      im1m=l1+1-m1
      im2m=l2+1-m2
      im2=l2+1+m2
410   a(im1,im2)=sgn*a(im1m,im2m)
420   continue
      l22p1=2*l2+1
      ind=mod(l2+is+1,2)+1
      ind2=l12p1-ind+1
      if(ind.gt.ind2) goto 720
      do 710 i=ind,ind2,2
      do 710 j=1,l22p1
      a(i,j)=-a(i,j)
710   continue
  720 return
      end
