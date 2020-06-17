!==================================
  MODULE simpraise_param
!==================================
  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER :: dp=SELECTED_REAL_KIND(14)	! precision up to 14 decimal digits
  INTEGER, PARAMETER :: tmax=1
  INTEGER, PARAMETER :: n_rep=1

  INTEGER :: nloci_X,nloci_Y,nloci_XY
  INTEGER :: nloci_sel,nloci_tot
  INTEGER :: nchrom
  LOGICAL :: file_exists

  TYPE candidate
    INTEGER :: self,dad,mum,gener,rank,sel
	REAL(KIND=dp) :: AddVal_X,AddVal_Y
	REAL(KIND=dp) :: DomDev_X,DomDev_Y
    REAL(KIND=dp) :: GenotVal_X,GenotVal_Y
	REAL(KIND=dp) :: fiped,fimol
    INTEGER :: genet_contrib ! contribution in number of gametes to next generation
    CHARACTER(LEN=1) :: sex
  END TYPE candidate
  
  TYPE mapping
    INTEGER :: chrom
    REAL(KIND=dp) :: dista
  END TYPE mapping

  TYPE(candidate), ALLOCATABLE, DIMENSION(:) :: ped
  
  TYPE(mapping), ALLOCATABLE, DIMENSION(:) :: gen_d_map

  ! whether individuals are dioecious [=0] or monoecious [=1]
  INTEGER(KIND=1) :: mono_dioecious
  ! whether simulation starts from founders [=0] OR resumes from advance generation [=1]
  INTEGER(KIND=1) :: status_gener

  INTEGER :: jrep,n_traits,counter
  INTEGER :: mnum,fnum,tot_founder
  INTEGER :: num_cand
  INTEGER :: jlast
  INTEGER :: current_ind, current_ind2
  INTEGER :: time, generation
  INTEGER :: dim_coa_mat
  ! jini was originally a vector of tmax elements, now tmax=1
  INTEGER :: jini

  INTEGER, DIMENSION(2) :: iseed

  INTEGER(KIND=1), ALLOCATABLE, DIMENSION(:) :: sel_Xs_map
  
  INTEGER, ALLOCATABLE, DIMENSION(:) :: gene_type
  INTEGER, ALLOCATABLE, DIMENSION(:) :: sel_sires,sel_dams
  INTEGER, ALLOCATABLE, DIMENSION(:) :: tot_Xs_map
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: par_alcopies,des_alcopies
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: mating_design
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: par_genome,des_genome
  
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: number_Xcross

  REAL(KIND=dp) :: xrand
  REAL(KIND=dp), DIMENSION(2) :: xdum
  REAL, PARAMETER :: zero=0.0_dp
    
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: freqini_vector
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: newa,olda
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: MGenot_X,SMGenot_X,SqMGenot_X
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: MGenot_Y,SMGenot_Y,SqMGenot_Y
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: MGenot_xx,SMGenot_xx,SqMGenot_xx
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: MGenot_yy,SMGenot_yy,SqMGenot_yy
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: MGenot_xy,SMGenot_xy,SqMGenot_xy
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: MGenot_yx,SMGenot_yx,SqMGenot_yx
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: VGenot_X,SVGenot_X,SqVGenot_X
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: VGenot_Y,SVGenot_Y,SqVGenot_Y
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: CovGenot_XY,SCovGenot_XY,SqCovGenot_XY
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: VaGenic_X,SVaGenic_X,SqVaGenic_X
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: VdGenic_X,SVdGenic_X,SqVdGenic_X
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: VaGenic_Y,SVaGenic_Y,SqVaGenic_Y
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: VdGenic_Y,SVdGenic_Y,SqVdGenic_Y
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Vadd_X,SVadd_X,SqVadd_X
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Vadd_Y,SVadd_Y,SqVadd_Y
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: CovAdd_XY,SCovAdd_XY,SqCovAdd_XY
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Vdom_X,SVdom_X,SqVdom_X
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Vdom_Y,SVdom_Y,SqVdom_Y
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: CovDesq_X11,SCovDesq_X11,SqCovDesq_X11
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: CovDesq_X12,SCovDesq_X12,SqCovDesq_X12
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: CovDesq_X22,SCovDesq_X22,SqCovDesq_X22
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: CovDesq_Y33,SCovDesq_Y33,SqCovDesq_Y33
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: CovDesq_Y23,SCovDesq_Y23,SqCovDesq_Y23
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: CovDesq_Y22,SCovDesq_Y22,SqCovDesq_Y22
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: CovDesq_XY13,SCovDesq_XY13,SqCovDesq_XY13
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: CovDesq_XY22,SCovDesq_XY22,SqCovDesq_XY22
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: CovDesq_XY12,SCovDesq_XY12,SqCovDesq_XY12
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: CovDesq_XY23,SCovDesq_XY23,SqCovDesq_XY23
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: CovPlei_XY22,SCovPlei_XY22,SqCovPlei_XY22
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Fimol_tot,SFimol_tot,SqFimol_tot
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Fiped_tot,SFiped_tot,SqFiped_tot
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Coaped,SCoaped,SqCoaped
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Coamol,SCoamol,SqCoamol
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Coa_newg, SCoa_newg, SqCoa_newg
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: gen_r_map

  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: newg
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: allele_freq,S_allele_freq,Sq_allele_freq
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: gene_effect
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: Fbylocus,SFbylocus,SqFbylocus
  
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:,:) :: par_Genot,des_Genot

  CHARACTER(LEN=50) :: rname
  CHARACTER(LEN=30) :: fname1
  CHARACTER(LEN=30) :: fname2
  CHARACTER(LEN=30) :: fname3
  CHARACTER(LEN=10) :: ctime
  CHARACTER(LEN=50) :: heading
!==================================
  END MODULE simpraise_param
!==================================
!==================================
  MODULE simpraise_stat
!==================================
  SAVE
  INTEGER, PARAMETER, PRIVATE :: dp=SELECTED_REAL_KIND(14)

  CONTAINS
!==================================
  FUNCTION real_mean(values)
!==================================
  IMPLICIT NONE
  REAL(KIND=dp) :: real_mean
  INTEGER :: n
  REAL(KIND=dp), DIMENSION(:) :: values
  n=SIZE(values)
  IF (n==0) THEN
    PRINT *, 'ERROR!!, attempt to calculate mean of zero elements in',&
    &' FUNCTION real_mean'
  ENDIF
  real_mean=SUM(values)/n
!==================================
  END FUNCTION real_mean
!==================================
!==================================
  FUNCTION real_var(values)
! calculates variance of a population
!==================================
  IMPLICIT NONE
  REAL(KIND=dp) :: real_var
  INTEGER :: n
  REAL(KIND=dp), DIMENSION(:) :: values
  n=SIZE(values)
  IF (n<2) THEN
    PRINT *, 'ERROR!!, attempt to calculate variance of only one element in',&
    &' FUNCTION real_var'
  ENDIF
  IF(SUM(values**2)>((SUM(values)**2)/n))THEN
    real_var=(SUM(values**2)-(SUM(values)**2)/n)/(n)
  ELSE
    real_var=0.0_dp
  END IF
!==================================
  END FUNCTION real_var
!==================================
!==================================
  FUNCTION real_cov(data1,data2)
!==================================
  IMPLICIT NONE
  REAL(KIND=dp) :: real_cov
  INTEGER :: n
  REAL(KIND=dp), DIMENSION(:) :: data1,data2
  n=SIZE(data1)
  IF (n<2) THEN
    PRINT *,'ERROR!!, attempt to calculate covariance of only one pair in',&
    &' FUNCTION real_cov'
  ELSE IF (SIZE(data1)/=SIZE(data2)) THEN
    PRINT *,'ERROR!!, unequal number of elements for covariance in',&
    &' FUNCTION real_cov'
  ENDIF
  IF(SUM(data1*data2)==(SUM(data1)*SUM(data2))/n)THEN
    real_cov=0.0_dp
  ELSE
    real_cov=(SUM(data1*data2)-(SUM(data1)*SUM(data2))/n)/n
  END IF
!==================================
  END FUNCTION real_cov
!==================================
!==================================
  FUNCTION real_cor(data1,data2)
!==================================
  IMPLICIT NONE
  INTEGER :: n
  REAL(KIND=dp) :: real_cor
  REAL(KIND=dp) :: var_1,var_2,covar_12
  REAL(KIND=dp), DIMENSION(:) :: data1,data2
  n=SIZE(data1)
  IF (n<2) THEN
    PRINT *, 'ERROR!!, attempt to calculate correlation of only one pair in',&
    &' FUNCTION real_cor'
  ELSE IF (SIZE(data1)/=SIZE(data2)) THEN
    PRINT *, 'ERROR!!, unequal number of elements for correlation in',&
    &' FUNCTION real_cor'
  ENDIF
  var_1=real_var(data1)
  var_2=real_var(data2)
  covar_12=(SUM(data1*data2)-(SUM(data1)*SUM(data2))/n)/n
  IF(var_1*var_2>0.0_dp)THEN
    real_cor=covar_12/DSQRT(var_1*var_2)
  ELSE
    real_cor=0.0_dp
  END IF
!==================================
  END FUNCTION real_cor
!==================================
!==================================
  FUNCTION real_reg(data1,data2)
! regression of data1 on data2
! data1 = b*data2 + e
!==================================
  IMPLICIT NONE
  INTEGER :: n
  REAL(KIND=dp) :: real_reg
  REAL(KIND=dp), DIMENSION(:) :: data1,data2
  n=SIZE(data1)
  IF (n<2) THEN
    PRINT *, 'ERROR!!, attempt to calculate regression of only one pair in',&
      &' FUNCTION real_reg'
    STOP
  ELSE IF (SIZE(data1)/=SIZE(data2)) THEN
    PRINT *,'ERROR!!, unequal number of elements for regression in',&
      &' FUNCTION real_reg'
    STOP
  ENDIF
  real_reg=real_cov(data1,data2)/real_var(data2)
!==================================
  END FUNCTION real_reg
!==================================
!==================================
  FUNCTION n_rec(row,col,max_dim)
!==================================
  IMPLICIT NONE
  INTEGER :: n_rec
  INTEGER :: max_dim
  INTEGER :: row,col
  n_rec=(((row-1)*max_dim)-(((row*row-row)/2)))+col
!==================================
  END FUNCTION n_rec
!==================================
!==================================
  FUNCTION mean_rep(data1,elements)
!==================================
  IMPLICIT NONE
  INTEGER :: elements
  REAL(KIND=dp) :: mean_rep
  REAL(KIND=dp) :: data1
  mean_rep=data1/elements
!==================================
  END FUNCTION mean_rep
!==================================
!==================================
  FUNCTION sem_rep(data1,data2,elements)
!==================================
  IMPLICIT NONE
  INTEGER :: elements
  REAL(KIND=dp) :: sem_rep
  REAL(KIND=dp) :: data1,data2
  IF(elements<=1)THEN
    sem_rep=0.0_dp
  ELSE
    sem_rep=(data1-((data2**2)/elements))/(elements-1)
  END IF
!==================================
  END FUNCTION sem_rep
!==================================
!==================================
  FUNCTION kosambi_rec(m_dist)
! estimates recombination fraction
! from distance in centi-Morgans
!==================================
  IMPLICIT NONE
  REAL(KIND=dp) :: m_dist
  REAL(KIND=dp) :: kosambi_rec
  
  kosambi_rec=0.5_dp*((exp(0.04_dp*m_dist)-1)/(exp(0.04_dp*m_dist)+1))
!==================================
  END FUNCTION kosambi_rec
!==================================
!==================================
  FUNCTION which_branch(wb_min,wb_max)
! randomly picks one
! of the two branches
!==================================
  IMPLICIT NONE
  INTEGER :: which_branch
  INTEGER :: wb_min,wb_max
  REAL(KIND=dp) :: xrand_branch
  
  CALL RANDOM_NUMBER(xrand_branch)
  which_branch=FLOOR(xrand_branch*wb_max)+wb_min
!==================================
  END FUNCTION which_branch
!==================================
!==================================
  END MODULE simpraise_stat
!==================================
!==================================
  MODULE simpraise_normal_table
!==================================
  SAVE
  INTEGER, PARAMETER, PRIVATE :: dp=SELECTED_REAL_KIND(14)
  INTEGER, PARAMETER, PRIVATE :: rr=selected_real_kind(p=14)

  CONTAINS
!==================================
  REAL(KIND=dp) FUNCTION n_dev(x)
!==================================
  IMPLICIT NONE
  REAL(KIND=dp), INTENT(INOUT) :: x
  INTEGER :: ix
!
  ix=0
  CALL random_number(x)
  n_dev=gcef(x,ix)
  IF(ix/=0) then
    print *, 'gcef fail'
    stop
  END if
!==================================
  END FUNCTION n_dev
!==================================
!----------------------------------
  REAL(KIND=dp) FUNCTION gcef(q,ix)
! calculates normal deviate x from lower tail proportion p
! by using order statistics
! LSR 20/07/04 (adapted from JAW)
!----------------------------------
  IMPLICIT NONE
  INTEGER, OPTIONAL :: ix
  REAL(KIND=dp) :: q
  REAL(KIND=dp) :: p
  REAL(KIND=dp) :: pp
  REAL(KIND=dp) :: u
  REAL(KIND=dp) :: t
  REAL(KIND=dp) :: x
  REAL(KIND=dp) :: zero=0.0_dp
  REAL(KIND=dp) :: one=1.0_dp
  REAL(KIND=dp) :: half=0.5_dp
  REAL(KIND=dp) :: a1=2.515517_dp
  REAL(KIND=dp) :: a2=0.802853_dp
  REAL(KIND=dp) :: a3=0.010328_dp
  REAL(KIND=dp) :: b1=1.432788_dp
  REAL(KIND=dp) :: b2=0.189269_dp
  REAL(KIND=dp) :: b3=0.001308_dp

  p=REAL(q,kind=dp)
  SELECT CASE((p>=zero).and.(p<=one))
  CASE(.true.)
  IF(PRESENT(ix)) ix=0
  IF(p>=one) THEN
    gcef=7.0_dp
  ELSE IF(p<1.0e-10_dp) THEN
    gcef=-7.0_dp
  ELSE
! calculates for probabilities gt half by using remainder
    IF(p>half) THEN
      pp=one-p
    ELSE
      pp=p
    END IF
    u=LOG(one/(pp*pp))
    t=DSQRT(u)
    x=(a1+(a2*t)+(a3*u))
    x=x/(one+(b1*t)+(b2*u)+b3*EXP(3.0_dp*LOG(t)))
    IF(p>half) THEN
      gcef=t-x
    ELSE
      gcef=x-t
    END IF
  END IF

  CASE DEFAULT
  IF(PRESENT(ix)) ix=1
    PRINT *, 'GCEF: probability out of bounds; ignore results'
    gcef=zero
  END SELECT
!----------------------------------
  END FUNCTION gcef
!----------------------------------
!==================================
  END MODULE simpraise_normal_table
!==================================
!==================================
  MODULE simpraise_routines
!==================================
  USE simpraise_param
  USE simpraise_stat
  USE simpraise_normal_table
  SAVE

  CONTAINS
!==================================
  SUBROUTINE setup_base 
!==================================
  IMPLICIT NONE
  INTEGER :: i,j,k,m

  DO i=1,tot_founder
    ped(i)%self=0
    ped(i)%dad=0
    ped(i)%mum=0
    ped(i)%gener=0

!   IMPORTANT for function of %rank
!   %rank indicates the position of the candidate in its
!   generation, in between 1 and num_cand. The same position
!   is allocated for this candidate in the matrices par_genome,
!   newa, olda, par_alcopies, par_Genot, mating_design

    ped(i)%rank=0
    ped(i)%sel=0
    ped(i)%AddVal_X=0.0_dp
    ped(i)%AddVal_Y=0.0_dp
    ped(i)%DomDev_X=0.0_dp
    ped(i)%DomDev_Y=0.0_dp
    ped(i)%GenotVal_X=0.0_dp
    ped(i)%GenotVal_Y=0.0_dp
    ped(i)%fiped=0.0_dp
    ped(i)%fimol=0.0_dp
    ped(i)%genet_contrib=0
  END DO

  generation=0
  
  FORALL(j=1:num_cand,k=1:nloci_tot,m=1:2) par_genome(j,k,m)=0
  FORALL(j=1:num_cand,k=1:nloci_sel) par_alcopies(j,k)=0
  FORALL(j=1:num_cand,k=1:nloci_sel,m=1:2) par_Genot(j,k,m)=0.0_dp
  FORALL(j=1:SIZE(sel_sires)) sel_sires(j)=0
  FORALL(j=1:SIZE(sel_dams)) sel_dams(j)=0

! FOUNDERS ARE SORTED BY SEX

  jlast=0
  IF(mono_dioecious==0)THEN
      DO i=1,tot_founder
        jlast=jlast+1
        ped(jlast)%self=jlast
        ped(jlast)%dad=0
        ped(jlast)%mum=0
        ped(jlast)%gener=generation
        ped(jlast)%sex='x'
        counter=0
        selected_loci_unisexuals: DO j=1,nloci_tot,2 ! uneven numbers
          counter=counter+1
          selected_loci_allele_unisexuals: DO k=1,2
            CALL random_number(xrand)
            IF(xrand<freqini_vector(counter))THEN
              par_genome(i,j,k)=1
            ELSE IF(xrand>=freqini_vector(counter))THEN
              par_genome(i,j,k)=0
            END IF
          END DO selected_loci_allele_unisexuals
        END DO selected_loci_unisexuals
        neutral_loci_unisexuals: DO j=2,nloci_tot,2 ! even numbers
          par_genome(i,j,1)=i
          par_genome(i,j,2)=i+mnum+fnum
        END DO neutral_loci_unisexuals
        sel_sires(i)=jlast
        sel_dams(i)=jlast ! equal vector to sel_sires
        ped(jlast)%rank=i
        ped(jlast)%sel=1
      END DO ! number of founders by population
  ELSE IF(mono_dioecious==1)THEN
!     loop sires
      DO i=1,mnum ! number of sires by population
        jlast=jlast+1
        ped(jlast)%self=jlast
        ped(jlast)%dad=0
        ped(jlast)%mum=0
        ped(jlast)%gener=generation
        ped(jlast)%sex='m'
        counter=0
        selected_loci_sires: DO j=1,nloci_tot,2 ! uneven numbers
          counter=counter+1
          selected_loci_allele_sires: DO k=1,2
            CALL random_number(xrand)
            IF(xrand<freqini_vector(counter))THEN
              par_genome(i,j,k)=1
            ELSE IF(xrand>=freqini_vector(counter))THEN
              par_genome(i,j,k)=0
            END IF
          END DO selected_loci_allele_sires
        END DO selected_loci_sires
        neutral_loci_sires: DO j=2,nloci_tot,2 ! even numbers
          par_genome(i,j,1)=i
          par_genome(i,j,2)=i+mnum+fnum
        END DO neutral_loci_sires
        sel_sires(i)=jlast
        ped(jlast)%rank=i
        ped(jlast)%sel=1
      END DO ! number of sires by population
!     loop dams
      DO i=1,fnum ! number of dams by population
        jlast=jlast+1
        ped(jlast)%self=jlast
        ped(jlast)%dad=0
        ped(jlast)%mum=0
        ped(jlast)%gener=generation
        ped(jlast)%sex='f'
        counter=0
        selected_loci_dams: DO j=1,nloci_tot,2 ! uneven numbers
          counter=counter+1
          selected_loci_allele_dams: DO k=1,2
            CALL random_number(xrand)
            IF(xrand<freqini_vector(counter))THEN
              par_genome(i+mnum,j,k)=1
            ELSE IF(xrand>=freqini_vector(counter))THEN
              par_genome(i+mnum,j,k)=0
            END IF
          END DO selected_loci_allele_dams
        END DO selected_loci_dams
        neutral_loci_dams: DO j=2,nloci_tot,2 ! even numbers
          par_genome(i+mnum,j,1)=i+mnum
          par_genome(i+mnum,j,2)=i+mnum+mnum+fnum
        END DO neutral_loci_dams
        sel_dams(i)=jlast
        ped(jlast)%rank=i+mnum
        ped(jlast)%sel=1
      END DO ! number of dams by population
  END IF ! mono_dioecious, whether 0 or 1
!==================================
  END SUBROUTINE setup_base
!==================================
!==================================
  SUBROUTINE sumup_founderGenot
!==================================
  IMPLICIT NONE
  INTEGER :: i,j
  INTEGER :: self_rank

    DO i=1,tot_founder
      counter=0
      self_rank=ped(i)%rank
      DO j=1,nloci_tot,2
        counter=counter+1
        par_alcopies(self_rank,counter)=SUM(par_genome(self_rank,j,:))
        IF(gene_type(counter)==1)THEN
          IF(par_alcopies(self_rank,counter)==0)THEN
            par_Genot(self_rank,counter,1)=-gene_effect(counter,1)
          ELSE IF(par_alcopies(self_rank,counter)==1)THEN
            par_Genot(self_rank,counter,1)=gene_effect(counter,2)
          ELSE IF(par_alcopies(self_rank,counter)==2)THEN
            par_Genot(self_rank,counter,1)=gene_effect(counter,1)
          END IF
        ELSE IF(gene_type(counter)==2)THEN
          IF(par_alcopies(self_rank,counter)==0)THEN
            par_Genot(self_rank,counter,1)=-gene_effect(counter,1)
            par_Genot(self_rank,counter,2)=-gene_effect(counter,3)
          ELSE IF(par_alcopies(self_rank,counter)==1)THEN
            par_Genot(self_rank,counter,1)=gene_effect(counter,2)
            par_Genot(self_rank,counter,2)=gene_effect(counter,4)
          ELSE IF(par_alcopies(self_rank,counter)==2)THEN
            par_Genot(self_rank,counter,1)=gene_effect(counter,1)
            par_Genot(self_rank,counter,2)=gene_effect(counter,3)
          END IF
        ELSE IF(gene_type(counter)==3)THEN
          IF(par_alcopies(self_rank,counter)==0)THEN
            par_Genot(self_rank,counter,2)=-gene_effect(counter,3)
          ELSE IF(par_alcopies(self_rank,counter)==1)THEN
            par_Genot(self_rank,counter,2)=gene_effect(counter,4)
          ELSE IF(par_alcopies(self_rank,counter)==2)THEN
            par_Genot(self_rank,counter,2)=gene_effect(counter,3)
          END IF
        END IF
      END DO
    END DO
!==================================
  END SUBROUTINE sumup_founderGenot
!==================================
!==================================
  SUBROUTINE create_offspring
! modified from METAGENE for simpraise
! LSR LSR 24/07/09
!==================================
  IMPLICIT NONE
  INTEGER :: i,j,k,m
  INTEGER :: max_fam_size,act_fam_size
  INTEGER :: sex_sample
  INTEGER :: branch,sire,dam,num_offs_recom,num_off
  INTEGER, ALLOCATABLE, DIMENSION(:) :: self_males,self_fmales
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: sire_mirror,dam_mirror
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: sire_gametes,dam_gametes

  jini=jlast+1
  
  max_fam_size=0
!  DO i=1, SIZE(mating_design,DIM=1) ! rows, female or maternal contributions
!    DO j=1, SIZE(mating_design,DIM=2) ! columns, male or paternal contributions
!      IF(mating_design(i,j)>max_fam_size)THEN
!        max_fam_size=mating_design(i,j)
!      END IF
!    END DO
!  END DO
  max_fam_size=MAXVAL(mating_design)

  IF(mono_dioecious==0)THEN ! rows: mnum+fnum x cols: mnum+fnum

    ALLOCATE(self_males(num_cand),self_fmales(num_cand))
    ALLOCATE(sire_gametes(max_fam_size,nloci_tot),dam_gametes(max_fam_size,nloci_tot))
    ALLOCATE(sire_mirror(2,nloci_tot),dam_mirror(2,nloci_tot))
    FORALL(i=1:max_fam_size,j=1:nloci_tot) sire_gametes(i,j)=0
    FORALL(i=1:max_fam_size,j=1:nloci_tot) dam_gametes(i,j)=0
    FORALL(i=1:2,j=1:nloci_tot) sire_mirror(i,j)=0
    FORALL(i=1:2,j=1:nloci_tot) dam_mirror(i,j)=0
    FORALL(j=1:num_cand,k=1:nloci_tot,m=1:2) des_genome(j,k,m)=0  

      FORALL(i=1:num_cand) self_males(i)=sel_sires(i)
      FORALL(i=1:num_cand) self_fmales(i)=sel_dams(i)
      num_off=0
      sires_monoecious: DO sire=1,mnum+fnum ! adding +fnum allows for selfing and reciprocals
        dams_monoecious: DO dam=1,mnum+fnum ! adding +fnum allows for selfing and reciprocals
          ! TODOS do not allow for selfs and change dim 1 for females and dim 2 for males
          IF(mating_design(dam,sire)>0)THEN
            DO j=1,nloci_tot
              DO k=1,2
                sire_mirror(k,j)=par_genome(ped(self_males(sire))%rank,j,k)
                dam_mirror(k,j)=par_genome(ped(self_fmales(dam))%rank,j,k)
              END DO
            END DO
            act_fam_size=mating_design(dam,sire)
            gametes_monoecious: DO i=1,act_fam_size
              ! do sire side
              CALL mask_rec(gen_r_map, sel_Xs_map, nloci_sel)
              CALL verify_Xcross(sire, dam, sel_Xs_map)
              branch=which_branch(1,2)
              CALL branch_rec(sel_Xs_map, tot_Xs_map, nloci_sel, branch)
              FORALL(j=1:nloci_tot) sire_gametes(i,j)=sire_mirror(tot_Xs_map(j),j)
              ! do dam side
              CALL mask_rec(gen_r_map, sel_Xs_map, nloci_sel)
              CALL verify_Xcross(sire, dam, sel_Xs_map)
              branch=which_branch(1,2)
              CALL branch_rec(sel_Xs_map, tot_Xs_map, nloci_sel, branch)
              FORALL(j=1:nloci_tot) dam_gametes(i,j)=dam_mirror(tot_Xs_map(j),j)
            END DO gametes_monoecious
            offspring_monoecious: DO i=1,act_fam_size
              jlast=jlast+1
              num_off=num_off+1
              ped(jlast)%self=jlast
              ped(jlast)%dad=self_males(sire)
              ped(jlast)%mum=self_fmales(dam)
              ped(jlast)%gener=generation
              ped(jlast)%sex='x'
              ped(jlast)%rank=num_off
              ped(jlast)%sel=0
              DO j=1,nloci_tot
                des_genome(num_off,j,1)=sire_gametes(i,j)
                des_genome(num_off,j,2)=dam_gametes(i,j)
              END DO
            END DO offspring_monoecious
          END IF ! IF(mating_design(dam,sire)>0)THEN
        END DO dams_monoecious
      END DO sires_monoecious 
      
  ELSE IF(mono_dioecious==1)THEN

    ALLOCATE(self_males(num_cand/2),self_fmales(num_cand/2))
    ALLOCATE(sire_gametes(max_fam_size,nloci_tot),dam_gametes(max_fam_size,nloci_tot))
    ALLOCATE(sire_mirror(2,nloci_tot),dam_mirror(2,nloci_tot))
    FORALL(i=1:max_fam_size,j=1:nloci_tot) sire_gametes(i,j)=0
    FORALL(i=1:max_fam_size,j=1:nloci_tot) dam_gametes(i,j)=0
    FORALL(i=1:2,j=1:nloci_tot) sire_mirror(i,j)=0
    FORALL(i=1:2,j=1:nloci_tot) dam_mirror(i,j)=0
    FORALL(j=1:num_cand,k=1:nloci_tot,m=1:2) des_genome(j,k,m)=0

      FORALL(i=1:num_cand/2) self_males(i)=sel_sires(i)
      FORALL(i=1:num_cand/2) self_fmales(i)=sel_dams(i)
      counter=0
      num_off=0
      sires: DO sire=1,mnum
        dams: DO dam=1,fnum
          ! TODOS change dim 1 for females and dim 2 for males
          IF(mating_design(dam,sire)>0)THEN
            DO j=1,nloci_tot
              DO k=1,2
                sire_mirror(k,j)=par_genome(ped(self_males(sire))%rank,j,k)
                dam_mirror(k,j)=par_genome(ped(self_fmales(dam))%rank,j,k)
              END DO
            END DO
            act_fam_size=mating_design(dam,sire)
            gametes_dioecious: DO i=1,act_fam_size
              ! do sire side
              CALL mask_rec(gen_r_map, sel_Xs_map, nloci_sel)
              CALL verify_Xcross(sire, dam, sel_Xs_map)
              branch=which_branch(1,2)
              CALL branch_rec(sel_Xs_map, tot_Xs_map, nloci_sel, branch)
              FORALL(j=1:nloci_tot) sire_gametes(i,j)=sire_mirror(tot_Xs_map(j),j)
              ! do dam side
              CALL mask_rec(gen_r_map, sel_Xs_map, nloci_sel)
              CALL verify_Xcross(sire, dam, sel_Xs_map)
              branch=which_branch(1,2)
              CALL branch_rec(sel_Xs_map, tot_Xs_map, nloci_sel, branch)
              FORALL(j=1:nloci_tot) dam_gametes(i,j)=dam_mirror(tot_Xs_map(j),j)
            END DO gametes_dioecious
            offspring_dioecious_males: DO i=1,act_fam_size-1,2
              jlast=jlast+1
              num_off=num_off+1
              ped(jlast)%self=jlast
              ped(jlast)%dad=self_males(sire)
              ped(jlast)%mum=self_fmales(dam)
              ped(jlast)%gener=generation
              ped(jlast)%sex='m'
              ped(jlast)%rank=num_off
              ped(jlast)%sel=0
              DO j=1,nloci_tot
                des_genome(num_off,j,1)=sire_gametes(i,j)
                des_genome(num_off,j,2)=dam_gametes(i,j)
              END DO
            END DO offspring_dioecious_males
            offspring_dioecious: DO i=2,act_fam_size,2
              jlast=jlast+1
              num_off=num_off+1
              ped(jlast)%self=jlast
              ped(jlast)%dad=self_males(sire)
              ped(jlast)%mum=self_fmales(dam)
              ped(jlast)%gener=generation
              ped(jlast)%sex='f'
              ped(jlast)%rank=num_off
              ped(jlast)%sel=0
              DO j=1,nloci_tot
                des_genome(num_off,j,1)=sire_gametes(i,j)
                des_genome(num_off,j,2)=dam_gametes(i,j)
              END DO
            END DO offspring_dioecious
          END IF ! IF(mating_design(dam,sire)==1)THEN
        END DO dams
      END DO sires

  END IF ! mono_dioecious

  DEALLOCATE(dam_gametes,sire_gametes)
  DEALLOCATE(sire_mirror,dam_mirror)
  DEALLOCATE(self_males,self_fmales)
!==================================
  END SUBROUTINE create_offspring
!==================================
!==================================
  SUBROUTINE sumup_offspringGenot
! modified from METAGENE for simpraise
!==================================
  IMPLICIT NONE
  INTEGER :: i,j,k,m

  FORALL(j=1:num_cand,k=1:nloci_sel,m=1:2) des_Genot(j,k,m)=0.0_dp
  FORALL(j=1:num_cand,k=1:nloci_sel) des_alcopies(j,k)=0
    DO j=1,num_cand
      counter=0
      DO k=1,nloci_tot,2
        counter=counter+1
        des_alcopies(j,counter)=SUM(des_genome(j,k,:))
        IF(gene_type(counter)==1)THEN
          IF(des_alcopies(j,counter)==0)THEN
            des_Genot(j,counter,1)=-gene_effect(counter,1)
          ELSE IF(des_alcopies(j,counter)==1)THEN
            des_Genot(j,counter,1)=gene_effect(counter,2)
          ELSE IF(des_alcopies(j,counter)==2)THEN
            des_Genot(j,counter,1)=gene_effect(counter,1)
          END IF
        ELSE IF(gene_type(counter)==2)THEN
          IF(des_alcopies(j,counter)==0)THEN
            des_Genot(j,counter,1)=-gene_effect(counter,1)
            des_Genot(j,counter,2)=-gene_effect(counter,3)
          ELSE IF(des_alcopies(j,counter)==1)THEN
            des_Genot(j,counter,1)=gene_effect(counter,2)
            des_Genot(j,counter,2)=gene_effect(counter,4)
          ELSE IF(des_alcopies(j,counter)==2)THEN
            des_Genot(j,counter,1)=gene_effect(counter,1)
            des_Genot(j,counter,2)=gene_effect(counter,3)
          END IF
        ELSE IF(gene_type(counter)==3)THEN
          IF(des_alcopies(j,counter)==0)THEN
            des_Genot(j,counter,2)=-gene_effect(counter,3)
          ELSE IF(des_alcopies(j,counter)==1)THEN
            des_Genot(j,counter,2)=gene_effect(counter,4)
          ELSE IF(des_alcopies(j,counter)==2)THEN
            des_Genot(j,counter,2)=gene_effect(counter,3)
          END IF
        END IF ! (gene_type(counter)==1)
      END DO ! k=1,nloci_tot,2
    END DO ! j=1,num_cand
!==================================
  END SUBROUTINE sumup_offspringGenot
!==================================
!==================================
  SUBROUTINE allelic_freq
! calculates allelic frequencies
! per locus across genome
!==================================
  IMPLICIT NONE
  INTEGER :: i,j

  IF(time==0)THEN
    DO i=1,nloci_sel
      allele_freq(1,i)=REAL(SUM(par_alcopies(:,i)),dp)/(2*(tot_founder))
    END DO
  ELSE IF(time>=1)THEN
    DO i=1,nloci_sel
      allele_freq(2,i)=REAL(SUM(des_alcopies(:,i)),dp)/(2*num_cand)
    END DO
  END IF
!==================================
  END SUBROUTINE allelic_freq
!==================================
!==================================
  SUBROUTINE base_coa_reset
!==================================
  IMPLICIT NONE
  INTEGER :: i,j
  INTEGER :: rec_tria1

  newa=0.0_dp
  DO j=1,mnum+fnum
    rec_tria1=n_rec(j,j,num_cand)
!     these are additive relationships
!     therefore double as much as kinships
    newa(rec_tria1)=1.0_dp
  END DO

  Fiped_tot(1)=0.0_dp
  Coaped(1)=0.5_dp/(mnum+fnum)
  Coamol(1)=0.5_dp/(mnum+fnum)
!==================================
  END SUBROUTINE base_coa_reset
!==================================
!==================================
  SUBROUTINE update_var_and_means
!==================================
  IMPLICIT NONE
  INTEGER :: i,j,k,m,tot_ind,init
  INTEGER :: count_X,count_Y,count_XY
  INTEGER :: self_rank
  REAL(KIND=dp) :: rvar1,rvar2,rvar3,rvar4,rvar5,rvar6,rvar7,rvar8,rvar9,rvar0
  REAL(KIND=dp) :: alpha
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Genot_data,Copy_data,data_1,data_2
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Genot_xx,Genot_yy,Genot_xy,Genot_yx
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: VgenX_by_LocX,VgenY_by_LocY
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: VgenX_by_LocXY,VgenY_by_LocXY
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: VdomX_by_LocX,VdomY_by_LocY
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: VdomX_by_LocXY,VdomY_by_LocXY
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:,:) :: AddVal,DomDev

  IF(time==0)THEN
    tot_ind=mnum+fnum
  ELSE IF(time>=1)THEN
    tot_ind=num_cand
  END IF
  init=jini

  ALLOCATE(Genot_data(tot_ind),Copy_data(tot_ind),data_1(tot_ind),data_2(tot_ind))
  ALLOCATE(AddVal(nloci_sel,tot_ind,2))
  ALLOCATE(DomDev(nloci_sel,tot_ind,2))
  FORALL(i=1:tot_ind) Genot_data(i)=0.0_dp
  FORALL(i=1:tot_ind) Copy_data(i)=0.0_dp
  FORALL(i=1:tot_ind) data_1(i)=0.0_dp
  FORALL(i=1:tot_ind) data_2(i)=0.0_dp
  FORALL(i=1:nloci_sel,j=1:tot_ind,k=1:2) AddVal(i,j,k)=0.0_dp
  FORALL(i=1:nloci_sel,j=1:tot_ind,k=1:2) DomDev(i,j,k)=0.0_dp      
  IF(nloci_X>0)THEN
    ALLOCATE(Genot_xx(nloci_X))
    ALLOCATE(VgenX_by_LocX(nloci_X))
    ALLOCATE(VdomX_by_LocX(nloci_X))
    FORALL(i=1:nloci_X) Genot_xx(i)=0.0_dp
    FORALL(i=1:nloci_X) VgenX_by_LocX(i)=0.0_dp
    FORALL(i=1:nloci_X) VdomX_by_LocX(i)=0.0_dp
  END IF
  IF(nloci_Y>0)THEN
    ALLOCATE(Genot_yy(nloci_Y))
    ALLOCATE(VgenY_by_LocY(nloci_Y))
    ALLOCATE(VdomY_by_LocY(nloci_Y))
    FORALL(i=1:nloci_Y) Genot_yy(i)=0.0_dp
    FORALL(i=1:nloci_Y) VgenY_by_LocY(i)=0.0_dp
    FORALL(i=1:nloci_Y) VdomY_by_LocY(i)=0.0_dp
  END IF
  IF(nloci_XY>0)THEN
    ALLOCATE(Genot_xy(nloci_XY),Genot_yx(nloci_XY))
    ALLOCATE(VgenX_by_LocXY(nloci_XY))
    ALLOCATE(VdomX_by_LocXY(nloci_XY))
    ALLOCATE(VgenY_by_LocXY(nloci_XY))
    ALLOCATE(VdomY_by_LocXY(nloci_XY))
    FORALL(i=1:nloci_XY) Genot_xy(i)=0.0_dp
    FORALL(i=1:nloci_XY) Genot_yx(i)=0.0_dp
    FORALL(i=1:nloci_XY) VgenX_by_LocXY(i)=0.0_dp
    FORALL(i=1:nloci_XY) VdomX_by_LocXY(i)=0.0_dp
    FORALL(i=1:nloci_XY) VgenY_by_LocXY(i)=0.0_dp
    FORALL(i=1:nloci_XY) VdomY_by_LocXY(i)=0.0_dp
  END IF
  
! Genic variance & dominance variance

  IF(time==0)THEN
    count_X=0
    count_Y=0
    count_XY=0
    DO i=1,nloci_sel
      IF(gene_type(i)==1)THEN
        count_X=count_X+1
        FORALL(j=1:tot_ind) Genot_data(j)=par_Genot(ped(j)%rank,i,1)
        FORALL(j=1:tot_ind) Copy_data(j)=REAL(par_alcopies(ped(j)%rank,i),dp)
        IF(real_var(Copy_data)==0.0_dp) alpha=0.0_dp
        IF(real_var(Copy_data)>0.0_dp) alpha=real_reg(Genot_data,Copy_data)
        FORALL(j=1:tot_ind) data_2(j)=alpha*Copy_data(j)
        rvar1=real_mean(Genot_data)
        rvar2=real_mean(data_2)
        IF(rvar1==rvar2)THEN
          FORALL(j=1:tot_ind) data_1(j)=data_2(j)
        ELSE
          FORALL(j=1:tot_ind) data_1(j)=data_2(j)+rvar1-rvar2
        END IF
        FORALL(j=1:tot_ind) AddVal(i,j,1)=data_1(j)
        VgenX_by_LocX(count_X)=real_var(data_1)
        IF(gene_effect(i,2)==0.0_dp)THEN
          VdomX_by_LocX(count_X)=0.0_dp
        ELSE
          FORALL(j=1:tot_ind) data_2(j)=Genot_data(j)-data_1(j)
          rvar3=real_mean(data_2)
          FORALL(j=1:tot_ind) DomDev(i,j,1)=data_2(j)-rvar3
          FORALL(j=1:tot_ind) data_2(j)=DomDev(i,j,1)
          VdomX_by_LocX(count_X)=real_var(data_2)
        END IF
        Genot_xx(count_X)=SUM(par_Genot(:,i,1))
      ELSE IF(gene_type(i)==2)THEN
        count_XY=count_XY+1
        FORALL(j=1:tot_ind) Genot_data(j)=par_Genot(ped(j)%rank,i,1)
        FORALL(j=1:tot_ind) Copy_data(j)=REAL(par_alcopies(ped(j)%rank,i),dp)
        IF(real_var(Copy_data)==0.0_dp) alpha=0.0_dp
        IF(real_var(Copy_data)>0.0_dp) alpha=real_reg(Genot_data,Copy_data)
        FORALL(j=1:tot_ind) data_2(j)=alpha*Copy_data(j)
        rvar1=real_mean(Genot_data)
        rvar2=real_mean(data_2)
        IF(rvar1==rvar2)THEN
          FORALL(j=1:tot_ind) data_1(j)=data_2(j)
        ELSE
          FORALL(j=1:tot_ind) data_1(j)=data_2(j)+rvar1-rvar2
        END IF
        FORALL(j=1:tot_ind) AddVal(i,j,1)=data_1(j)
        VgenX_by_LocXY(count_XY)=real_var(data_1)
        IF(gene_effect(i,2)==0.0_dp)THEN
          VdomX_by_LocXY(count_XY)=0.0_dp
        ELSE
          FORALL(j=1:tot_ind) data_2(j)=Genot_data(j)-data_1(j)
          rvar3=real_mean(data_2)
          FORALL(j=1:tot_ind) DomDev(i,j,1)=data_2(j)-rvar3
          FORALL(j=1:tot_ind) data_2(j)=DomDev(i,j,1)
          VdomX_by_LocXY(count_XY)=real_var(data_2)
        END IF
        FORALL(j=1:tot_ind) Genot_data(j)=par_Genot(ped(j)%rank,i,2)
        IF(real_var(Copy_data)==0.0_dp) alpha=0.0_dp
        IF(real_var(Copy_data)>0.0_dp) alpha=real_reg(Genot_data,Copy_data)
        FORALL(j=1:tot_ind) data_2(j)=alpha*Copy_data(j)
        rvar1=real_mean(Genot_data)
        rvar2=real_mean(data_2)
        IF(rvar1==rvar2)THEN
          FORALL(j=1:tot_ind) data_1(j)=data_2(j)
        ELSE
          FORALL(j=1:tot_ind) data_1(j)=data_2(j)+rvar1-rvar2
        END IF
        FORALL(j=1:tot_ind) AddVal(i,j,2)=data_1(j)
        VgenY_by_LocXY(count_XY)=real_var(data_1)
        IF(gene_effect(i,4)==0.0_dp)THEN
          VdomY_by_LocXY(count_XY)=0.0_dp
        ELSE
          FORALL(j=1:tot_ind) data_2(j)=Genot_data(j)-data_1(j)
          rvar3=real_mean(data_2)
          FORALL(j=1:tot_ind) DomDev(i,j,2)=data_2(j)-rvar3
          FORALL(j=1:tot_ind) data_2(j)=DomDev(i,j,2)
          VdomY_by_LocXY(count_XY)=real_var(data_2)
        END IF
        Genot_xy(count_XY)=SUM(par_Genot(:,i,1))
        Genot_yx(count_XY)=SUM(par_Genot(:,i,2))
      ELSE IF(gene_type(i)==3)THEN
        count_Y=count_Y+1
        FORALL(j=1:tot_ind) Genot_data(j)=par_Genot(ped(j)%rank,i,2)
        FORALL(j=1:tot_ind) Copy_data(j)=REAL(par_alcopies(ped(j)%rank,i),dp)
        IF(real_var(Copy_data)==0.0_dp) alpha=0.0_dp
        IF(real_var(Copy_data)>0.0_dp) alpha=real_reg(Genot_data,Copy_data)
        FORALL(j=1:tot_ind) data_2(j)=alpha*Copy_data(j)
        rvar1=real_mean(Genot_data)
        rvar2=real_mean(data_2)
        IF(rvar1==rvar2)THEN
          FORALL(j=1:tot_ind) data_1(j)=data_2(j)
        ELSE
          FORALL(j=1:tot_ind) data_1(j)=data_2(j)+rvar1-rvar2
        END IF
        FORALL(j=1:tot_ind) AddVal(i,j,2)=data_1(j)
        VgenY_by_LocY(count_Y)=real_var(data_1)
        IF(gene_effect(i,4)==0.0_dp)THEN
          VdomY_by_LocY(count_Y)=0.0_dp
        ELSE
          FORALL(j=1:tot_ind) data_2(j)=Genot_data(j)-data_1(j)
          rvar3=real_mean(data_2)
          FORALL(j=1:tot_ind) DomDev(i,j,2)=data_2(j)-rvar3
          FORALL(j=1:tot_ind) data_2(j)=DomDev(i,j,2)
          VdomY_by_LocY(count_Y)=real_var(data_2)
        END IF
        Genot_yy(count_Y)=SUM(par_Genot(:,i,2))
      END IF
    END DO
    current_ind=0
    DO i=1,jlast
        current_ind=current_ind+1
        self_rank=ped(i)%rank
        ped(i)%AddVal_X=SUM(AddVal(:,current_ind,1))
        ped(i)%AddVal_Y=SUM(AddVal(:,current_ind,2))
        ped(i)%DomDev_X=SUM(DomDev(:,current_ind,1))
        ped(i)%DomDev_Y=SUM(DomDev(:,current_ind,2))
        ped(i)%GenotVal_X=SUM(par_Genot(self_rank,:,1))
        ped(i)%GenotVal_Y=SUM(par_Genot(self_rank,:,2))
    END DO
! genic variances before covariance due to linkage
    IF((nloci_X>0).AND.(nloci_XY==0))THEN
      VaGenic_X(time+1)=SUM(VgenX_by_LocX,MASK=VgenX_by_LocX.NE.zero)! Va_genic_X
      VdGenic_X(time+1)=SUM(VdomX_by_LocX,MASK=VdomX_by_LocX.NE.zero)! Vd_genic_X
    ELSE IF((nloci_X>0).AND.(nloci_XY>0))THEN
      VaGenic_X(time+1)=SUM(VgenX_by_LocX,MASK=VgenX_by_LocX.NE.zero)+SUM(VgenX_by_LocXY,MASK=VgenX_by_LocXY.NE.zero) ! Va_genic_X
      VdGenic_X(time+1)=SUM(VdomX_by_LocX,MASK=VdomX_by_LocX.NE.zero)+SUM(VdomX_by_LocXY,MASK=VdomX_by_LocXY.NE.zero) ! Vd_genic_X
    END IF
    IF((nloci_Y>0).AND.(nloci_XY==0))THEN
      VaGenic_Y(time+1)=SUM(VgenY_by_LocY,MASK=VgenY_by_LocY.NE.zero)! Va_genic_Y
      VdGenic_Y(time+1)=SUM(VdomY_by_LocY,MASK=VdomY_by_LocY.NE.zero)! Vd_genic_Y
    ELSE IF((nloci_Y>0).AND.(nloci_XY>0))THEN
      VaGenic_Y(time+1)=SUM(VgenY_by_LocY,MASK=VgenY_by_LocY.NE.zero)+SUM(VgenY_by_LocXY,MASK=VgenY_by_LocXY.NE.zero) ! Va_genic_Y
      VdGenic_Y(time+1)=SUM(VdomY_by_LocY,MASK=VdomY_by_LocY.NE.zero)+SUM(VdomY_by_LocXY,MASK=VdomY_by_LocXY.NE.zero) ! Vd_genic_Y
    END IF
    ! calculates genetic mean from each type of loci
    IF(nloci_X>0) MGenot_xx(time+1)=SUM(Genot_xx)/tot_ind
    IF(nloci_Y>0) MGenot_yy(time+1)=SUM(Genot_yy)/tot_ind
    IF(count_XY>0) MGenot_xy(time+1)=SUM(Genot_xy)/tot_ind
    IF(count_XY>0) MGenot_yx(time+1)=SUM(Genot_yx)/tot_ind
  ELSE IF(time>=1)THEN ! IF(time==0)THEN
    count_X=0
    count_Y=0
    count_XY=0
    DO i=1,nloci_sel
      IF(gene_type(i)==1)THEN
        count_X=count_X+1
        FORALL(j=1:tot_ind) Genot_data(j)=des_Genot(j,i,1)
        FORALL(j=1:tot_ind) Copy_data(j)=REAL(des_alcopies(j,i),dp)
        IF(real_var(Copy_data)==0.0_dp) alpha=0.0_dp
        IF(real_var(Copy_data)>0.0_dp) alpha=real_reg(Genot_data,Copy_data)
        FORALL(j=1:tot_ind) data_2(j)=alpha*Copy_data(j)
        rvar1=real_mean(Genot_data)
        rvar2=real_mean(data_2)
        IF(rvar1==rvar2)THEN
          FORALL(j=1:tot_ind) data_1(j)=data_2(j)
        ELSE
          FORALL(j=1:tot_ind) data_1(j)=data_2(j)+rvar1-rvar2
        END IF
        FORALL(j=1:tot_ind) AddVal(i,j,1)=data_1(j)
        VgenX_by_LocX(count_X)=real_var(data_1)
        IF(gene_effect(i,2)==0.0_dp)THEN
          VdomX_by_LocX(count_X)=0.0_dp
        ELSE
          FORALL(j=1:tot_ind) data_2(j)=Genot_data(j)-data_1(j)
          rvar3=real_mean(data_2)
          FORALL(j=1:tot_ind) DomDev(i,j,1)=data_2(j)-rvar3
          FORALL(j=1:tot_ind) data_2(j)=DomDev(i,j,1)
          VdomX_by_LocX(count_X)=real_var(data_2)
        END IF
        Genot_xx(count_X)=SUM(des_Genot(:,i,1))
      ELSE IF(gene_type(i)==2)THEN
        count_XY=count_XY+1
        FORALL(j=1:tot_ind) Genot_data(j)=des_Genot(j,i,1)
        FORALL(j=1:tot_ind) Copy_data(j)=REAL(des_alcopies(j,i),dp)
        IF(real_var(Copy_data)==0.0_dp) alpha=0.0_dp
        IF(real_var(Copy_data)>0.0_dp) alpha=real_reg(Genot_data,Copy_data)
        FORALL(j=1:tot_ind) data_2(j)=alpha*Copy_data(j)
        rvar1=real_mean(Genot_data)
        rvar2=real_mean(data_2)
        IF(rvar1==rvar2)THEN
          FORALL(j=1:tot_ind) data_1(j)=data_2(j)
        ELSE
          FORALL(j=1:tot_ind) data_1(j)=data_2(j)+rvar1-rvar2
        END IF
        FORALL(j=1:tot_ind) AddVal(i,j,1)=data_1(j)
        VgenX_by_LocXY(count_XY)=real_var(data_1)
        IF(gene_effect(i,2)==0.0_dp)THEN
          VdomX_by_LocXY(count_XY)=0.0_dp
        ELSE
          FORALL(j=1:tot_ind) data_2(j)=Genot_data(j)-data_1(j)
          rvar3=real_mean(data_2)
          FORALL(j=1:tot_ind) DomDev(i,j,1)=data_2(j)-rvar3
          FORALL(j=1:tot_ind) data_2(j)=DomDev(i,j,1)
          VdomX_by_LocXY(count_XY)=real_var(data_2)
        END IF
        FORALL(j=1:tot_ind) Genot_data(j)=des_Genot(j,i,2)
        IF(real_var(Copy_data)==0.0_dp) alpha=0.0_dp
        IF(real_var(Copy_data)>0.0_dp) alpha=real_reg(Genot_data,Copy_data)
        FORALL(j=1:tot_ind) data_2(j)=alpha*Copy_data(j)
        rvar1=real_mean(Genot_data)
        rvar2=real_mean(data_2)
        IF(rvar1==rvar2)THEN
          FORALL(j=1:tot_ind) data_1(j)=data_2(j)
        ELSE
          FORALL(j=1:tot_ind) data_1(j)=data_2(j)+rvar1-rvar2
        END IF
        FORALL(j=1:tot_ind) AddVal(i,j,2)=data_1(j)
        VgenY_by_LocXY(count_XY)=real_var(data_1)
        IF(gene_effect(i,4)==0.0_dp)THEN
          VdomY_by_LocXY(count_XY)=0.0_dp
        ELSE
          FORALL(j=1:tot_ind) data_2(j)=Genot_data(j)-data_1(j)
          rvar3=real_mean(data_2)
          FORALL(j=1:tot_ind) DomDev(i,j,2)=data_2(j)-rvar3
          FORALL(j=1:tot_ind) data_2(j)=DomDev(i,j,2)
          VdomY_by_LocXY(count_XY)=real_var(data_2)
        END IF
        Genot_xy(count_XY)=SUM(des_Genot(:,i,1))
        Genot_yx(count_XY)=SUM(des_Genot(:,i,2))
      ELSE IF(gene_type(i)==3)THEN
        count_Y=count_Y+1
        FORALL(j=1:tot_ind) Genot_data(j)=des_Genot(j,i,2)
        FORALL(j=1:tot_ind) Copy_data(j)=REAL(des_alcopies(j,i),dp)
        IF(real_var(Copy_data)==0.0_dp) alpha=0.0_dp
        IF(real_var(Copy_data)>0.0_dp) alpha=real_reg(Genot_data,Copy_data)
        FORALL(j=1:tot_ind) data_2(j)=alpha*Copy_data(j)
        rvar1=real_mean(Genot_data)
        rvar2=real_mean(data_2)
        IF(rvar1==rvar2)THEN
          FORALL(j=1:tot_ind) data_1(j)=data_2(j)
        ELSE
          FORALL(j=1:tot_ind) data_1(j)=data_2(j)+rvar1-rvar2
        END IF
        FORALL(j=1:tot_ind) AddVal(i,j,2)=data_1(j)
        VgenY_by_LocY(count_Y)=real_var(data_1)
        IF(gene_effect(i,4)==0.0_dp)THEN
          VdomY_by_LocY(count_Y)=0.0_dp
        ELSE
          FORALL(j=1:tot_ind) data_2(j)=Genot_data(j)-data_1(j)
          rvar3=real_mean(data_2)
          FORALL(j=1:tot_ind) DomDev(i,j,2)=data_2(j)-rvar3
          FORALL(j=1:tot_ind) data_2(j)=DomDev(i,j,2)
          VdomY_by_LocY(count_Y)=real_var(data_2)
        END IF
        Genot_yy(count_Y)=SUM(des_Genot(:,i,2))
      END IF
    END DO
    counter=0
    DO i=init,jlast
        counter=counter+1
        ped(i)%AddVal_X=SUM(AddVal(:,counter,1))
        ped(i)%AddVal_Y=SUM(AddVal(:,counter,2))
        ped(i)%DomDev_X=SUM(DomDev(:,counter,1))
        ped(i)%DomDev_Y=SUM(DomDev(:,counter,2))
        ped(i)%GenotVal_X=SUM(des_Genot(counter,:,1))
        ped(i)%GenotVal_Y=SUM(des_Genot(counter,:,2))
    END DO
! genic variances before covariance due to linkage
    IF((nloci_X>0).AND.(nloci_XY==0))THEN
      VaGenic_X(time+1)=SUM(VgenX_by_LocX,MASK=VgenX_by_LocX.NE.zero)! Va_genic_X
      VdGenic_X(time+1)=SUM(VdomX_by_LocX,MASK=VdomX_by_LocX.NE.zero)! Vd_genic_X
    ELSE IF((nloci_X>0).AND.(nloci_XY>0))THEN
      VaGenic_X(time+1)=SUM(VgenX_by_LocX,MASK=VgenX_by_LocX.NE.zero)+SUM(VgenX_by_LocXY,MASK=VgenX_by_LocXY.NE.zero) ! Va_genic_X
      VdGenic_X(time+1)=SUM(VdomX_by_LocX,MASK=VdomX_by_LocX.NE.zero)+SUM(VdomX_by_LocXY,MASK=VdomX_by_LocXY.NE.zero) ! Vd_genic_X
    END IF
    IF((nloci_Y>0).AND.(nloci_XY==0))THEN
      VaGenic_Y(time+1)=SUM(VgenY_by_LocY,MASK=VgenY_by_LocY.NE.zero)! Va_genic_Y
      VdGenic_Y(time+1)=SUM(VdomY_by_LocY,MASK=VdomY_by_LocY.NE.zero)! Vd_genic_Y
    ELSE IF((nloci_Y>0).AND.(nloci_XY>0))THEN
      VaGenic_Y(time+1)=SUM(VgenY_by_LocY,MASK=VgenY_by_LocY.NE.zero)+SUM(VgenY_by_LocXY,MASK=VgenY_by_LocXY.NE.zero) ! Va_genic_Y
      VdGenic_Y(time+1)=SUM(VdomY_by_LocY,MASK=VdomY_by_LocY.NE.zero)+SUM(VdomY_by_LocXY,MASK=VdomY_by_LocXY.NE.zero) ! Vd_genic_Y
    END IF
    ! calculates genetic mean from each type of loci
    IF(nloci_X>0) MGenot_xx(time+1)=SUM(Genot_xx)/tot_ind
    IF(nloci_Y>0) MGenot_yy(time+1)=SUM(Genot_yy)/tot_ind
    IF(count_XY>0) MGenot_xy(time+1)=SUM(Genot_xy)/tot_ind
    IF(count_XY>0) MGenot_yx(time+1)=SUM(Genot_yx)/tot_ind
  END IF ! IF(time==0)THEN

! genic variances before covariance due to linkage
  IF(nloci_X>0) DEALLOCATE(VgenX_by_LocX)
  IF(nloci_Y>0) DEALLOCATE(VgenY_by_LocY)
  IF(nloci_XY>0) DEALLOCATE(VgenX_by_LocXY,VgenY_by_LocXY)
  IF(nloci_X>0) DEALLOCATE(VdomX_by_LocX)
  IF(nloci_Y>0) DEALLOCATE(VdomY_by_LocY)
  IF(nloci_XY>0) DEALLOCATE(VdomX_by_LocXY,VdomY_by_LocXY)
  DEALLOCATE(Genot_data,Copy_data)
  DEALLOCATE(AddVal,DomDev)

! the sum Vadd_+Vdom_ DOES NOT have to be equal to Vgenot_
! there is a covariance component between additive values (a) and
! dominance deviations (d) when considering finite populations and
! deviations from Hardy-Weinberg equilibrium, more precisely
! there are four components of desequilibrium covariance:
! cov(a,a), 2*cov(a,d) and cov(d,d)

! additive and dominance variances (with covariances due to linkage)

!   additive means, variances and covariances
    counter=0
    DO i=init,jlast
        counter=counter+1
        data_1(counter)=ped(i)%AddVal_X
        data_2(counter)=ped(i)%AddVal_Y
    END DO
    Vadd_X(time+1)=real_var(data_1)
    IF(n_traits==2)THEN
      Vadd_Y(time+1)=real_var(data_2)
      CovAdd_XY(time+1)=real_cov(data_1,data_2)
    END IF
!   dominance variances
    counter=0
    DO i=init,jlast
        counter=counter+1
        data_1(counter)=ped(i)%DomDev_X
        data_2(counter)=ped(i)%DomDev_Y
    END DO
    Vdom_X(time+1)=real_var(data_1)
    Vdom_Y(time+1)=real_var(data_2)
!   genotypic variances and covariances
    counter=0
    DO i=init,jlast
        counter=counter+1
        data_1(counter)=ped(i)%GenotVal_X
        data_2(counter)=ped(i)%GenotVal_Y
    END DO
    VGenot_X(time+1)=real_var(data_1)
    MGenot_X(time+1)=real_mean(data_1)
    IF(n_traits==2)THEN
      VGenot_Y(time+1)=real_var(data_2)
      MGenot_Y(time+1)=real_mean(data_2)
      CovGenot_XY(time+1)=real_cov(data_1,data_2)
    END IF

  IF(nloci_X>0) DEALLOCATE(Genot_xx)
  IF(nloci_Y>0) DEALLOCATE(Genot_yy)
  IF(nloci_XY>0) DEALLOCATE(Genot_xy,Genot_yx)

! calculates linkage desequilibrium covariances across genomes & candidates
! considering a single covariance between genotypic values (due to non-random
! associations (a,a), (d,d) and (a,d)

  IF(time==0)THEN
    rvar1=0.0_dp
    rvar2=0.0_dp
    rvar3=0.0_dp
    rvar4=0.0_dp
    rvar5=0.0_dp
    rvar6=0.0_dp
    rvar7=0.0_dp
    rvar8=0.0_dp
    rvar9=0.0_dp
    rvar0=0.0_dp
    DO i=1,nloci_sel-1
      DO j=i+1,nloci_sel
        IF((gene_type(i)==1).AND.(gene_type(j)==1))THEN
          FORALL(k=1:tot_ind) data_1(k)=par_Genot(ped(k)%rank,i,1)
          FORALL(k=1:tot_ind) data_2(k)=par_Genot(ped(k)%rank,j,1)
          rvar1=rvar1+real_cov(data_1,data_2)
        ELSE IF((gene_type(i)==1).AND.(gene_type(j)==2))THEN
          FORALL(k=1:tot_ind) data_1(k)=par_Genot(ped(k)%rank,i,1)
          FORALL(k=1:tot_ind) data_2(k)=par_Genot(ped(k)%rank,j,1)
          rvar2=rvar2+real_cov(data_1,data_2)
          FORALL(k=1:tot_ind) data_1(k)=par_Genot(ped(k)%rank,i,1)
          FORALL(k=1:tot_ind) data_2(k)=par_Genot(ped(k)%rank,j,2)
          rvar9=rvar9+real_cov(data_1,data_2)
        ELSE IF((gene_type(i)==2).AND.(gene_type(j)==1))THEN
          FORALL(k=1:tot_ind) data_1(k)=par_Genot(ped(k)%rank,i,1)
          FORALL(k=1:tot_ind) data_2(k)=par_Genot(ped(k)%rank,j,1)
          rvar2=rvar2+real_cov(data_1,data_2)
          FORALL(k=1:tot_ind) data_1(k)=par_Genot(ped(k)%rank,i,2)
          FORALL(k=1:tot_ind) data_2(k)=par_Genot(ped(k)%rank,j,1)
          rvar9=rvar9+real_cov(data_1,data_2)
        ELSE IF((gene_type(i)==2).AND.(gene_type(j)==2))THEN
          FORALL(k=1:tot_ind) data_1(k)=par_Genot(ped(k)%rank,i,1)
          FORALL(k=1:tot_ind) data_2(k)=par_Genot(ped(k)%rank,j,1)
          rvar5=rvar5+real_cov(data_1,data_2)
          FORALL(k=1:tot_ind) data_1(k)=par_Genot(ped(k)%rank,i,2)
          FORALL(k=1:tot_ind) data_2(k)=par_Genot(ped(k)%rank,j,2)
          rvar6=rvar6+real_cov(data_1,data_2)
          FORALL(k=1:tot_ind) data_1(k)=par_Genot(ped(k)%rank,i,1)
          FORALL(k=1:tot_ind) data_2(k)=par_Genot(ped(k)%rank,j,2)
          rvar8=rvar8+real_cov(data_1,data_2)
          FORALL(k=1:tot_ind) data_1(k)=par_Genot(ped(k)%rank,i,2)
          FORALL(k=1:tot_ind) data_2(k)=par_Genot(ped(k)%rank,j,1)
          rvar8=rvar8+real_cov(data_1,data_2)
        ELSE IF((gene_type(i)==3).AND.(gene_type(j)==3))THEN
          FORALL(k=1:tot_ind) data_1(k)=par_Genot(ped(k)%rank,i,2)
          FORALL(k=1:tot_ind) data_2(k)=par_Genot(ped(k)%rank,j,2)
          rvar3=rvar3+real_cov(data_1,data_2)
        ELSE IF((gene_type(i)==2).AND.(gene_type(j)==3))THEN
          FORALL(k=1:tot_ind) data_1(k)=par_Genot(ped(k)%rank,i,2)
          FORALL(k=1:tot_ind) data_2(k)=par_Genot(ped(k)%rank,j,2)
          rvar4=rvar4+real_cov(data_1,data_2)
          FORALL(k=1:tot_ind) data_1(k)=par_Genot(ped(k)%rank,i,1)
          FORALL(k=1:tot_ind) data_2(k)=par_Genot(ped(k)%rank,j,2)
          rvar0=rvar0+real_cov(data_1,data_2)
        ELSE IF((gene_type(i)==3).AND.(gene_type(j)==2))THEN
          FORALL(k=1:tot_ind) data_1(k)=par_Genot(ped(k)%rank,i,2)
          FORALL(k=1:tot_ind) data_2(k)=par_Genot(ped(k)%rank,j,2)
          rvar4=rvar4+real_cov(data_1,data_2)
          FORALL(k=1:tot_ind) data_1(k)=par_Genot(ped(k)%rank,i,2)
          FORALL(k=1:tot_ind) data_2(k)=par_Genot(ped(k)%rank,j,1)
          rvar0=rvar0+real_cov(data_1,data_2)
        ELSE IF((gene_type(i)==1).AND.(gene_type(j)==3))THEN
          FORALL(k=1:tot_ind) data_1(k)=par_Genot(ped(k)%rank,i,1)
          FORALL(k=1:tot_ind) data_2(k)=par_Genot(ped(k)%rank,j,2)
          rvar7=rvar7+real_cov(data_1,data_2)
        ELSE IF((gene_type(i)==3).AND.(gene_type(j)==1))THEN
          FORALL(k=1:tot_ind) data_1(k)=par_Genot(ped(k)%rank,i,2)
          FORALL(k=1:tot_ind) data_2(k)=par_Genot(ped(k)%rank,j,1)
          rvar7=rvar7+real_cov(data_1,data_2)
        END IF
      END DO ! j=i+1,nloci_sel
    END DO ! i=1,nloci_sel-1
    CovDesq_X11(time+1)=rvar1
    CovDesq_X12(time+1)=rvar2
    CovDesq_Y33(time+1)=rvar3
    CovDesq_Y23(time+1)=rvar4
    CovDesq_X22(time+1)=rvar5
    CovDesq_Y22(time+1)=rvar6
    CovDesq_XY13(time+1)=rvar7
    CovDesq_XY22(time+1)=rvar8
    CovDesq_XY12(time+1)=rvar9
    CovDesq_XY23(time+1)=rvar0
    rvar1=0.0_dp
    DO i=1,nloci_sel
      IF(gene_type(i)==2)THEN
        FORALL(k=1:tot_ind) data_1(k)=par_Genot(ped(k)%rank,i,1)
        FORALL(k=1:tot_ind) data_2(k)=par_Genot(ped(k)%rank,i,2)
        rvar1=rvar1+real_cov(data_1,data_2)
      END IF
    END DO ! i=1,nloci_sel
    CovPlei_XY22(time+1)=rvar1
  ELSE IF(time>=1)THEN
    rvar1=0.0_dp
    rvar2=0.0_dp
    rvar3=0.0_dp
    rvar4=0.0_dp
    rvar5=0.0_dp
    rvar6=0.0_dp
    rvar7=0.0_dp
    rvar8=0.0_dp
    rvar9=0.0_dp
    rvar0=0.0_dp
    DO i=1,nloci_sel-1
      DO j=i+1,nloci_sel
        IF((gene_type(i)==1).AND.(gene_type(j)==1))THEN
          FORALL(k=1:tot_ind) data_1(k)=des_Genot(k,i,1)
          FORALL(k=1:tot_ind) data_2(k)=des_Genot(k,j,1)
          rvar1=rvar1+real_cov(data_1,data_2)
        ELSE IF((gene_type(i)==1).AND.(gene_type(j)==2))THEN
          FORALL(k=1:tot_ind) data_1(k)=des_Genot(k,i,1)
          FORALL(k=1:tot_ind) data_2(k)=des_Genot(k,j,1)
          rvar2=rvar2+real_cov(data_1,data_2)
          FORALL(k=1:tot_ind) data_1(k)=des_Genot(k,i,1)
          FORALL(k=1:tot_ind) data_2(k)=des_Genot(k,j,2)
          rvar9=rvar9+real_cov(data_1,data_2)
        ELSE IF((gene_type(i)==2).AND.(gene_type(j)==1))THEN
          FORALL(k=1:tot_ind) data_1(k)=des_Genot(k,i,1)
          FORALL(k=1:tot_ind) data_2(k)=des_Genot(k,j,1)
          rvar2=rvar2+real_cov(data_1,data_2)
          FORALL(k=1:tot_ind) data_1(k)=des_Genot(k,i,2)
          FORALL(k=1:tot_ind) data_2(k)=des_Genot(k,j,1)
          rvar9=rvar9+real_cov(data_1,data_2)
        ELSE IF((gene_type(i)==2).AND.(gene_type(j)==2))THEN
          FORALL(k=1:tot_ind) data_1(k)=des_Genot(k,i,1)
          FORALL(k=1:tot_ind) data_2(k)=des_Genot(k,j,1)
          rvar5=rvar5+real_cov(data_1,data_2)
          FORALL(k=1:tot_ind) data_1(k)=des_Genot(k,i,2)
          FORALL(k=1:tot_ind) data_2(k)=des_Genot(k,j,2)
          rvar6=rvar6+real_cov(data_1,data_2)
          FORALL(k=1:tot_ind) data_1(k)=des_Genot(k,i,1)
          FORALL(k=1:tot_ind) data_2(k)=des_Genot(k,j,2)
          rvar8=rvar8+real_cov(data_1,data_2)
          FORALL(k=1:tot_ind) data_1(k)=des_Genot(k,i,2)
          FORALL(k=1:tot_ind) data_2(k)=des_Genot(k,j,1)
          rvar8=rvar8+real_cov(data_1,data_2)
        ELSE IF((gene_type(i)==3).AND.(gene_type(j)==3))THEN
          FORALL(k=1:tot_ind) data_1(k)=des_Genot(k,i,2)
          FORALL(k=1:tot_ind) data_2(k)=des_Genot(k,j,2)
          rvar3=rvar3+real_cov(data_1,data_2)
        ELSE IF((gene_type(i)==2).AND.(gene_type(j)==3))THEN
          FORALL(k=1:tot_ind) data_1(k)=des_Genot(k,i,2)
          FORALL(k=1:tot_ind) data_2(k)=des_Genot(k,j,2)
          rvar4=rvar4+real_cov(data_1,data_2)
          FORALL(k=1:tot_ind) data_1(k)=des_Genot(k,i,1)
          FORALL(k=1:tot_ind) data_2(k)=des_Genot(k,j,2)
          rvar0=rvar0+real_cov(data_1,data_2)
        ELSE IF((gene_type(i)==3).AND.(gene_type(j)==2))THEN
          FORALL(k=1:tot_ind) data_1(k)=des_Genot(k,i,2)
          FORALL(k=1:tot_ind) data_2(k)=des_Genot(k,j,2)
          rvar4=rvar4+real_cov(data_1,data_2)
          FORALL(k=1:tot_ind) data_1(k)=des_Genot(k,i,2)
          FORALL(k=1:tot_ind) data_2(k)=des_Genot(k,j,1)
          rvar0=rvar0+real_cov(data_1,data_2)
        ELSE IF((gene_type(i)==1).AND.(gene_type(j)==3))THEN
          FORALL(k=1:tot_ind) data_1(k)=des_Genot(k,i,1)
          FORALL(k=1:tot_ind) data_2(k)=des_Genot(k,j,2)
          rvar7=rvar7+real_cov(data_1,data_2)
        ELSE IF((gene_type(i)==3).AND.(gene_type(j)==1))THEN
          FORALL(k=1:tot_ind) data_1(k)=des_Genot(k,i,2)
          FORALL(k=1:tot_ind) data_2(k)=des_Genot(k,j,1)
          rvar7=rvar7+real_cov(data_1,data_2)
        END IF
      END DO ! j=i+1,nloci_sel
    END DO ! i=1,nloci_sel-1
    CovDesq_X11(time+1)=rvar1
    CovDesq_X12(time+1)=rvar2
    CovDesq_Y33(time+1)=rvar3
    CovDesq_Y23(time+1)=rvar4
    CovDesq_X22(time+1)=rvar5
    CovDesq_Y22(time+1)=rvar6
    CovDesq_XY13(time+1)=rvar7
    CovDesq_XY22(time+1)=rvar8
    CovDesq_XY12(time+1)=rvar9
    CovDesq_XY23(time+1)=rvar0
    rvar1=0.0_dp
    DO i=1,nloci_sel
      IF(gene_type(i)==2)THEN
        FORALL(k=1:tot_ind) data_1(k)=des_Genot(k,i,1)
        FORALL(k=1:tot_ind) data_2(k)=des_Genot(k,i,2)
        rvar1=rvar1+real_cov(data_1,data_2)
      END IF
    END DO ! i=1,nloci_sel
    CovPlei_XY22(time+1)=rvar1

  END IF ! (time==0) or (time>=1)

! Desequilibrium Covariance for trait X:   rvar1+rvar2+rvar5
! Desequilibrium & Pleiotropy Covariance for traits XY: rvar7+rvar8
! Desequilibrium Covariance for trait Y:   rvar3+rvar4+rvar6

  DEALLOCATE(data_1,data_2)
!==================================
  END SUBROUTINE update_var_and_means
!==================================
!==================================
!==================================
  SUBROUTINE rank_candidates
! ranks candidates on index values
! based on heapsort algorithm taking
!           Nlog(2)N steps
!
! from NUMERICAL RECIPES IN FORTRAN 77
! modified to make index a random variable
!==================================
  IMPLICIT NONE
  INTEGER :: i,j,m,position,count_males,count_fmales
  INTEGER :: dad_of_sel,mum_of_sel
  INTEGER :: rank_dad_of_sel,rank_mum_of_sel
  INTEGER, ALLOCATABLE, DIMENSION(:) :: self_males,self_fmales
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: rank_males,rank_fmales

! variables for heapsort ranking algorithm
  INTEGER :: nsort,lsort,irsort,rrself
  REAL(KIND=dp) :: rrindx

  IF(mono_dioecious==0)THEN ! monoecious is TRUE
    ALLOCATE(self_males(num_cand),self_fmales(num_cand))
    ALLOCATE(rank_males(num_cand),rank_fmales(num_cand))
    FORALL(i=1:SIZE(self_males)) self_males(i)=0
    FORALL(i=1:SIZE(self_fmales)) self_fmales(i)=0
    FORALL(i=1:SIZE(rank_males)) rank_males(i)=0.0_dp
    FORALL(i=1:SIZE(rank_fmales)) rank_fmales(i)=0.0_dp
    FORALL(j=1:num_cand) sel_sires(j)=0
    FORALL(j=1:num_cand) sel_dams(j)=0
      count_males=0
      count_fmales=0
      DO i=jini,jlast
          count_males=count_males+1
          self_males(count_males)=ped(i)%self
          CALL RANDOM_NUMBER(xrand)
          rank_males(count_males)=xrand
      END DO
!     ranking according to xrand
!     ----------ascending numerical order---------
      nsort=num_cand
      lsort=nsort/2+1
      irsort=nsort
104   CONTINUE
        IF(lsort.GT.1)THEN
          lsort=lsort-1
          rrself=self_males(lsort)
          rrindx=rank_males(lsort)
        ELSE
          rrindx=rank_males(irsort)
          rrself=self_males(irsort)
          rank_males(irsort)=rank_males(1)
          self_males(irsort)=self_males(1)
          irsort=irsort-1
          IF(irsort.EQ.1)THEN
            rank_males(1)=rrindx
            self_males(1)=rrself
            GOTO 304
          END IF
        END IF
        i=lsort
        j=lsort+lsort
204     IF(j.LE.irsort)THEN
          IF(j.LT.irsort)THEN
            IF(rank_males(j).LT.rank_males(j+1)) j=j+1
          END IF
          IF(rrindx.LT.rank_males(j))THEN
            rank_males(i)=rank_males(j)
            self_males(i)=self_males(j)
            i=j
            j=j+j
          ELSE
            j=irsort+1
          END IF
          GOTO 204
        END IF
        rank_males(i)=rrindx
        self_males(i)=rrself
        GOTO 104
304   CONTINUE
!     ascending numerical order, thus loop goes backwards
!     goes from num_cand to 'mnum+fnum' elements further down
      nsort=num_cand-mnum-fnum+1
      position=0
      DO i=num_cand,nsort,-1
        position=position+1
        sel_sires(position)=self_males(i)
        ! stores individuals who will be parents
        ped(sel_sires(position))%sel=1
      END DO
      FORALL(j=1:num_cand) sel_dams(j)=sel_sires(j)
  ELSE IF(mono_dioecious==1)THEN ! dioecious is TRUE
    ALLOCATE(self_males(num_cand/2),self_fmales(num_cand/2))
    ALLOCATE(rank_males(num_cand/2),rank_fmales(num_cand/2))
    FORALL(i=1:SIZE(self_males)) self_males(i)=0
    FORALL(i=1:SIZE(self_fmales)) self_fmales(i)=0
    FORALL(i=1:SIZE(rank_males)) rank_males(i)=0.0_dp
    FORALL(i=1:SIZE(rank_fmales)) rank_fmales(i)=0.0_dp
    FORALL(j=1:num_cand/2) sel_sires(j)=0
    FORALL(j=1:num_cand/2) sel_dams(j)=0
      count_males=0
      count_fmales=0
      DO i=jini,jlast
        IF(ped(i)%sex=='m')THEN
          count_males=count_males+1
          self_males(count_males)=ped(i)%self
          CALL RANDOM_NUMBER(xrand)
          rank_males(count_males)=xrand
        ELSE IF(ped(i)%sex=='f')THEN
          count_fmales=count_fmales+1
          self_fmales(count_fmales)=ped(i)%self
          CALL RANDOM_NUMBER(xrand)
          rank_fmales(count_fmales)=xrand
        END IF
      END DO
!     ranking according to xrand for males
!     ----------ascending numerical order---------
      nsort=count_males
      lsort=nsort/2+1
      irsort=nsort
102   CONTINUE
        IF(lsort.GT.1)THEN
          lsort=lsort-1
          rrself=self_males(lsort)
          rrindx=rank_males(lsort)
        ELSE
          rrindx=rank_males(irsort)
          rrself=self_males(irsort)
          rank_males(irsort)=rank_males(1)
          self_males(irsort)=self_males(1)
          irsort=irsort-1
          IF(irsort.EQ.1)THEN
            rank_males(1)=rrindx
            self_males(1)=rrself
            GOTO 302
          END IF
        END IF
        i=lsort
        j=lsort+lsort
202     IF(j.LE.irsort)THEN
          IF(j.LT.irsort)THEN
            IF(rank_males(j).LT.rank_males(j+1)) j=j+1
          END IF
          IF(rrindx.LT.rank_males(j))THEN
            rank_males(i)=rank_males(j)
            self_males(i)=self_males(j)
            i=j
            j=j+j
          ELSE
            j=irsort+1
          END IF
          GOTO 202
        END IF
        rank_males(i)=rrindx
        self_males(i)=rrself
        GOTO 102
302   CONTINUE
!     ascending numerical order, thus loop goes backwards
      nsort=count_males-mnum+1
      position=0
      DO i=count_males,nsort,-1
        position=position+1
        sel_sires(position)=self_males(i)
        ! stores individuals who will be parents
        ped(sel_sires(position))%sel=1
      END DO
!     ranking according to xrand for females
!     ----------ascending numerical order---------
      nsort=count_fmales
      lsort=nsort/2+1
      irsort=nsort
103   CONTINUE
        IF(lsort.GT.1)THEN
          lsort=lsort-1
          rrself=self_fmales(lsort)
          rrindx=rank_fmales(lsort)
        ELSE
          rrindx=rank_fmales(irsort)
          rrself=self_fmales(irsort)
          rank_fmales(irsort)=rank_fmales(1)
          self_fmales(irsort)=self_fmales(1)
          irsort=irsort-1
          IF(irsort.EQ.1)THEN
            rank_fmales(1)=rrindx
            self_fmales(1)=rrself
            GOTO 303
          END IF
        END IF
        i=lsort
        j=lsort+lsort
203     IF(j.LE.irsort)THEN
          IF(j.LT.irsort)THEN
            IF(rank_fmales(j).LT.rank_fmales(j+1)) j=j+1
          END IF
          IF(rrindx.LT.rank_fmales(j))THEN
            rank_fmales(i)=rank_fmales(j)
            self_fmales(i)=self_fmales(j)
            i=j
            j=j+j
          ELSE
          j=irsort+1
          END IF
          GOTO 203
        END IF
        rank_fmales(i)=rrindx
        self_fmales(i)=rrself
        GOTO 103
303   CONTINUE
!     ascending numerical order, thus loop goes backwards
      nsort=count_fmales-fnum+1
      position=0
      DO i=count_fmales,nsort,-1
        position=position+1
        sel_dams(position)=self_fmales(i)
        ! stores individuals who will be parents
        ped(sel_dams(position))%sel=1
      END DO
  END IF ! mono_dioecious 
  
  DEALLOCATE(self_males,self_fmales,rank_males,rank_fmales)
!==================================
  END SUBROUTINE rank_candidates
!==================================
!==================================
  SUBROUTINE update_ARM
! updated LSR 20/07/09 LSR
! calculates the upper triangle and diagonal of
! the Additive Relationship Matrix
! Result takes the form of a vector
! row index (i) MUST be <= col index (j)
! if we look for any pair(i,j) where i>j
! then take instead its symmetrical pair(j,i)
! Note that ARM elements are twice the kinship coefficient
! Inbreeding is given as diagonal value minus 1
! Coancestries are given as half the additive relationships
!==================================
  IMPLICIT NONE

  INTEGER :: i,j,m
  INTEGER :: ii,si,di,jj,sj,dj ! refer to %selfs
  INTEGER :: ii_rank,si_rank,di_rank ! refer to %ranks between 1 and num_cand
  INTEGER :: jj_rank,sj_rank,dj_rank ! refer to %ranks between 1 and num_cand
  INTEGER :: rec_tria1
  INTEGER :: rec_tria2
  INTEGER :: rec_tria3
  INTEGER :: rec_tria4
  INTEGER :: rec_tria5
  REAL(KIND=dp) :: Fiped_sum
  REAL(KIND=dp) :: Coaped_sum

  FORALL(i=1:SIZE(newa)) olda(i)=newa(i)
  FORALL(i=1:SIZE(newa)) newa(i)=0.0_dp

    Fiped_sum=0.0_dp
    Coaped_sum=0.0_dp
    ii_rank=0
    across_cohort: DO i=jini,jlast
        ii_rank=ii_rank+1
        rec_tria1=n_rec(ii_rank,ii_rank,num_cand) ! self diagonal position
        si=ped(i)%dad
        di=ped(i)%mum
        si_rank=ped(si)%rank
        di_rank=ped(di)%rank
        IF(si_rank<=di_rank) rec_tria2=n_rec(si_rank,di_rank,num_cand)    ! parental relationship matrix position
        IF(di_rank<si_rank) rec_tria2=n_rec(di_rank,si_rank,num_cand)     ! parental relationship matrix position
        ped(i)%fiped=0.5_dp*olda(rec_tria2)
        newa(rec_tria1)=1.0_dp+0.5_dp*olda(rec_tria2)
        Fiped_sum=Fiped_sum+ped(i)%fiped
!       it accumulates through the diagonal elements of the ARM
!       half of the additive relationship, which is kinship or coancestry
        Coaped_sum=Coaped_sum+newa(rec_tria1)/2
        jj_rank=ii_rank
        down_cohort: DO j=i+1,jlast
            jj_rank=jj_rank+1
            sj=ped(j)%dad
            dj=ped(j)%mum
            rec_tria1=n_rec(ii_rank,jj_rank,num_cand) ! position of i - j pair in the matrix
            sj_rank=ped(sj)%rank
            dj_rank=ped(dj)%rank
            IF(si_rank<=sj_rank) rec_tria2=n_rec(si_rank,sj_rank,num_cand)
            IF(sj_rank<si_rank) rec_tria2=n_rec(sj_rank,si_rank,num_cand)
            IF(si_rank<=dj_rank) rec_tria3=n_rec(si_rank,dj_rank,num_cand)
            IF(dj_rank<si_rank) rec_tria3=n_rec(dj_rank,si_rank,num_cand)
            IF(di_rank<=sj_rank) rec_tria4=n_rec(di_rank,sj_rank,num_cand)
            IF(sj_rank<di_rank) rec_tria4=n_rec(sj_rank,di_rank,num_cand)
            IF(di_rank<=dj_rank) rec_tria5=n_rec(di_rank,dj_rank,num_cand)
            IF(dj_rank<di_rank) rec_tria5=n_rec(dj_rank,di_rank,num_cand)
            newa(rec_tria1)=(olda(rec_tria2)+olda(rec_tria3)+olda(rec_tria4)+olda(rec_tria5))/4
!           it accumulates through the upper triangle,
!           but coancestry average has to take into account
!           two equivalent triangles: upper and lower triangles
            Coaped_sum=Coaped_sum+newa(rec_tria1)
        END DO down_cohort
    END DO across_cohort
    Fiped_tot(time+1)=Fiped_sum/num_cand
    Coaped(time+1)=Coaped_sum/(num_cand*num_cand)
    
  INQUIRE(FILE="arm.r", EXIST=file_exists)
  IF(file_exists)THEN ! it's replaced
    OPEN(UNIT=1,FILE='arm.r', status="replace", action="write")
  ELSE ! it's created anew
    OPEN(UNIT=1,FILE='arm.r', status="new", action="write")
  END IF
  DO i=1,SIZE(newa)
    WRITE(1,'(F13.10)')newa(i)
  END DO
  CLOSE(1)
!==================================
  END SUBROUTINE update_ARM
!==================================
!==================================
  SUBROUTINE update_molmat
! updates the IBDs by neutral markers
!==================================
  IMPLICIT NONE
  INTEGER :: i,j,m,ii,current_locus
  INTEGER :: founder_1,founder_2
  INTEGER :: founder_3,founder_4
  INTEGER :: inbred_events,inbred_events_ind
  INTEGER :: rec_tria1
  INTEGER, ALLOCATABLE, DIMENSION (:) :: InbreedEventsbyLocus
  INTEGER, ALLOCATABLE, DIMENSION (:,:) :: Count_neutral
  INTEGER, ALLOCATABLE, DIMENSION (:,:) :: Count_IBD
  REAL(KIND=dp) :: SCoa_dig, SCoa_off
  REAL(KIND=dp), ALLOCATABLE, DIMENSION (:,:) :: IBD_pop

  ALLOCATE(InbreedEventsbyLocus(nloci_sel))
  ALLOCATE(Count_neutral(nloci_sel,(mnum+fnum)*2))
  ALLOCATE(IBD_pop(nloci_sel,(mnum+fnum)*2))
  ALLOCATE(Count_IBD(nloci_sel,dim_coa_mat))
  newg=0.0_dp

    FORALL(i=1:nloci_sel) InbreedEventsbyLocus(i)=0
    FORALL(i=1:nloci_sel,j=1:(mnum+fnum)*2) Count_neutral(i,j)=0
    FORALL(i=1:nloci_sel,j=1:(mnum+fnum)*2) IBD_pop(i,j)=0.0_dp
    FORALL(i=1:nloci_sel,j=1:dim_coa_mat) Count_IBD(i,j)=0
    current_ind=0
    inbred_events=0
    DO i=jini,jlast
        current_ind=current_ind+1
        inbred_events_ind=0
        current_locus=0
        DO j=2,nloci_tot,2
          current_locus=current_locus+1
          founder_1=des_genome(current_ind,j,1)
          founder_2=des_genome(current_ind,j,2)
          IF(founder_1==founder_2)THEN
            inbred_events=inbred_events+1
            inbred_events_ind=inbred_events_ind+1
            InbreedEventsbyLocus(current_locus)=InbreedEventsbyLocus(current_locus)+1
          END IF
          Count_neutral(current_locus,founder_1)=Count_neutral(current_locus,founder_1)+1
          Count_neutral(current_locus,founder_2)=Count_neutral(current_locus,founder_2)+1
          current_ind2=0
          DO ii=i,jlast
              current_ind2=current_ind2+1
              rec_tria1=n_rec(current_ind,current_ind2,num_cand) ! position of i - j pair in the matrix
              founder_3=des_genome(current_ind2,j,1)
              founder_4=des_genome(current_ind2,j,2)
              IF(founder_1==founder_3) Count_IBD(current_locus,rec_tria1)=Count_IBD(current_locus,rec_tria1)+1
              IF(founder_1==founder_4) Count_IBD(current_locus,rec_tria1)=Count_IBD(current_locus,rec_tria1)+1
              IF(founder_2==founder_3) Count_IBD(current_locus,rec_tria1)=Count_IBD(current_locus,rec_tria1)+1
              IF(founder_2==founder_4) Count_IBD(current_locus,rec_tria1)=Count_IBD(current_locus,rec_tria1)+1
          END DO ! ii=i,jlast
        END DO ! j=2,nloci_tot,2
        ped(i)%fimol=REAL(inbred_events_ind,dp)/(nloci_sel)
    END DO ! i=jini,jlast
    FORALL(i=1:nloci_sel,j=1:(mnum+fnum)*2) IBD_pop(i,j)=(REAL(Count_neutral(i,j),dp)/(2*num_cand))**2
    FORALL(i=1:nloci_sel) Fbylocus(time,i)=REAL(InbreedEventsbyLocus(i),dp)/num_cand
    Fimol_tot(time+1)=REAL(inbred_events,dp)/(num_cand*nloci_sel)
    Coamol(time+1)=SUM(IBD_pop)/(nloci_sel)
    FORALL(i=1:nloci_sel,j=1:dim_coa_mat) newg(i,j)=REAL(Count_IBD(i,j),dp)/4

    current_ind=0
    SCoa_dig=0.0_dp
    SCoa_off=0.0_dp
    DO i=jini,jlast
        current_ind=current_ind+1
        current_ind2=0
        DO ii=i,jlast
            current_ind2=current_ind2+1
            rec_tria1=n_rec(current_ind,current_ind2,num_cand)
            IF(i==ii) SCoa_dig=SCoa_dig+SUM(newg(:,rec_tria1))
            IF(i/=ii) SCoa_off=SCoa_off+SUM(newg(:,rec_tria1))
        END DO ! ii=i,jlast
    END DO ! i=jini,jlast
    Coa_newg(time+1)=(2*SCoa_off+SCoa_dig)/(nloci_sel*num_cand*num_cand)
 
  DEALLOCATE(InbreedEventsbyLocus)
  DEALLOCATE(Count_neutral)
  DEALLOCATE(IBD_pop)
  DEALLOCATE(Count_IBD)
!==================================
  END SUBROUTINE update_molmat
!==================================
!==================================
  SUBROUTINE new_genome
! fill parental chromosomes with
! selected offspring following
! ped(.)%rank order
!==================================
  IMPLICIT NONE
  INTEGER :: i,j,k

  FORALL(i=1:num_cand,j=1:nloci_tot,k=1:2) par_genome(i,j,k)=0

  IF(mono_dioecious==0)THEN
      counter=0
      DO i=jini,jlast
          counter=counter+1
          DO j=1,nloci_tot
            DO k=1,2
              par_genome(counter,j,k)=des_genome(counter,j,k)
            END DO
          END DO
      END DO
  ELSE IF(mono_dioecious==1)THEN
      counter=0
      DO i=jini,jlast
          counter=counter+1
          IF(ped(i)%sex=='m')THEN
            DO j=1,nloci_tot
              DO k=1,2
                par_genome(counter,j,k)=des_genome(counter,j,k)
              END DO
            END DO
          ELSE IF(ped(i)%sex=='f')THEN
            DO j=1,nloci_tot
              DO k=1,2
                par_genome(counter,j,k)=des_genome(counter,j,k)
              END DO
            END DO
          END IF
      END DO
  END IF
!==================================
  END SUBROUTINE new_genome
!==================================
!==================================
  SUBROUTINE backup_ped
! pedigree can grow over generations
! so size of file is not fixed
!==================================
  IMPLICIT NONE
  INTEGER :: i
  
  INQUIRE(FILE="ped.r", EXIST=file_exists)

  IF(status_gener==1)THEN ! ped file must exists from previous run
    IF(file_exists)THEN
      OPEN(UNIT=1,FILE='ped.r', status="replace", action="write")
    ELSE
      WRITE(*,*)'ERROR! ped file does not exists'
      STOP
    END IF
  ELSE IF(status_gener==0)THEN ! ped file if exists must be replaced
    IF(file_exists)THEN ! it's replaced
      OPEN(UNIT=1,FILE='ped.r', status="replace", action="write")
    ELSE ! it's created anew
      OPEN(UNIT=1,FILE='ped.r', status="new", action="write")
    END IF
  END IF

  INQUIRE(FILE="ped_size.r", EXIST=file_exists)

  IF(file_exists)THEN
    OPEN(UNIT=2,FILE='ped_size.r', status="replace", action="write")
  ELSE
    OPEN(UNIT=2,FILE='ped_size.r', status="new", action="write")
  END IF

  WRITE(2,'(I8)')jlast
    
! backup pedigree info into sequential file
  DO i=1,jlast
    WRITE(1,'(6I8,8F12.5,A2)') &
&     ped(i)%self,&
&     ped(i)%dad,&
&     ped(i)%mum,&
&     ped(i)%gener,&
&     ped(i)%rank,&
&     ped(i)%sel,&
&     ped(i)%AddVal_X,&
&     ped(i)%AddVal_Y,&
&     ped(i)%DomDev_X,&
&     ped(i)%DomDev_Y,&
&     ped(i)%GenotVal_X,&
&     ped(i)%GenotVal_Y,&
&     ped(i)%fiped,&
&     ped(i)%fimol,&
&     ped(i)%sex
  END DO

  CLOSE(1)  
!==================================
  END SUBROUTINE backup_ped
!==================================
!==================================
  SUBROUTINE backup_genomes
!==================================
  IMPLICIT NONE
  INTEGER :: i,j,k,rec_length

  rec_length=KIND(par_genome(1,1,1))
  counter=0
      
! backup genome info into direct access file
  IF(time==0)THEN
    OPEN(UNIT=1,FILE='gen0.r',FORM='unformatted',ACCESS='direct',RECL=rec_length)
    DO i=1,tot_founder
      DO j=1,nloci_tot
        DO k=1,2
          counter=counter+1
          WRITE(1,REC=counter) par_genome(i,j,k)
        END DO
      END DO
    END DO
  ELSE IF(time>=1)THEN
    OPEN(UNIT=1,FILE='gen1.r',FORM='unformatted',ACCESS='direct',RECL=rec_length)
    DO i=1,num_cand
      DO j=1,nloci_tot
        DO k=1,2
          counter=counter+1
          WRITE(1,REC=counter) des_genome(i,j,k)
        END DO
      END DO
    END DO
  END IF
  
  CLOSE(1)  
  
  ! backup genome info into sequential file, just for checking
  IF(time==0)THEN
    INQUIRE(FILE="gen0_sec.r", EXIST=file_exists)
    IF(file_exists)THEN
      OPEN(UNIT=1,FILE='gen0_sec.r', status="REPLACE", action="write")
    ELSE
      OPEN(UNIT=1,FILE='gen0_sec.r', status="new", action="write")
    END IF
    DO i=1,tot_founder
      DO j=1,nloci_tot
        DO k=1,2
          WRITE(1,'(5I8)') ped(i)%self, ped(i)%rank, j, k, par_genome(ped(i)%rank,j,k)
        END DO
      END DO
    END DO
  ELSE IF(time>=1)THEN
    INQUIRE(FILE="gen0_sec.r", EXIST=file_exists)
    IF(file_exists)THEN
      OPEN(UNIT=1,FILE='gen1_sec.r', status="REPLACE", action="write")
    ELSE
      OPEN(UNIT=1,FILE='gen1_sec.r', status="new", action="write")
    END IF
    DO i=1,num_cand
      DO j=1,nloci_tot
        DO k=1,2
          WRITE(1,'(5I8)') ped(i)%self, ped(i)%rank, j, k, des_genome(ped(i)%rank,j,k)
        END DO
      END DO
    END DO
  END IF
  
  CLOSE(1)  
!==================================
  END SUBROUTINE backup_genomes
!==================================  
!==================================
  SUBROUTINE recover_ped
!==================================
  IMPLICIT NONE
  INTEGER :: i
  INTEGER :: status_ped ! I/O status
  INTEGER :: count_sires, count_dams
  INTEGER :: max_gener
  INTEGER, ALLOCATABLE, DIMENSION(:) :: gener_vector

  ALLOCATE(gener_vector(jlast))
  gener_vector=0
  
  INQUIRE(FILE="ped.r", EXIST=file_exists)
  
  IF(file_exists)THEN
    OPEN(UNIT=1,FILE='ped.r', status="old", action="read")
  ELSE
    WRITE(*,*)'ERROR! file ped.r does not exist: ', file_exists
    STOP
  END IF

  ! recover pedigree info into sequential file
  DO i=1,jlast
    READ(1,'(6I8,8F12.5,A2)') &
&     ped(i)%self,&
&     ped(i)%dad,&
&     ped(i)%mum,&
&     ped(i)%gener,&
&     ped(i)%rank,&
&     ped(i)%sel,&
&     ped(i)%AddVal_X,&
&     ped(i)%AddVal_Y,&
&     ped(i)%DomDev_X,&
&     ped(i)%DomDev_Y,&
&     ped(i)%GenotVal_X,&
&     ped(i)%GenotVal_Y,&
&     ped(i)%fiped,&
&     ped(i)%fimol,&
&     ped(i)%sex
    gener_vector(i)=ped(i)%gener
  END DO

  max_gener=gener_vector(jlast)
  generation=max_gener+1
  
  FORALL(i=1:SIZE(sel_sires)) sel_sires(i)=0
  FORALL(i=1:SIZE(sel_dams)) sel_dams(i)=0

  ! recover from pedigree sel_sires & sel_dams
  IF(mono_dioecious==0)THEN
    count_sires=0
    count_dams=0
    DO i=jlast,jlast-num_cand+1,-1 ! takes the last num_cand individuals
      IF(ped(i)%sel==1)THEN
        count_sires=count_sires+1
        count_dams=count_dams+1
        sel_sires(count_sires)=ped(i)%self
        sel_dams(count_dams)=ped(i)%self
      END IF
    END DO
  ELSE IF(mono_dioecious==1)THEN
    count_sires=0
    count_dams=0
    DO i=jlast,jlast-num_cand+1,-1 ! takes the last num_cand individuals
      IF((ped(i)%sel==1).AND.(ped(i)%sex=='m'))THEN
        count_sires=count_sires+1
        sel_sires(count_sires)=ped(i)%self
      ELSE IF((ped(i)%sel==1).AND.(ped(i)%sex=='f'))THEN
        count_dams=count_dams+1
        sel_dams(count_dams)=ped(i)%self
      END IF
    END DO
  END IF
  
  CLOSE(1)  
  DEALLOCATE(gener_vector)
!==================================
  END SUBROUTINE recover_ped
!==================================
!==================================
  SUBROUTINE recover_genomes
!==================================
  IMPLICIT NONE
  INTEGER :: i,j,k,rec_length

  rec_length=KIND(par_genome(1,1,1)) ! par_genome(1,1,1) is just an example
      
! read genome info from sequential file
! and into parental matrice, as this is going
! to be used to generate offspring
  
  counter=0
  OPEN(UNIT=1,FILE='gen1.r',FORM='unformatted',ACCESS='direct',RECL=rec_length)
  DO i=1,num_cand
    DO j=1,nloci_tot
      DO k=1,2
        counter=counter+1
        READ(1,REC=counter) par_genome(i,j,k)
      END DO
    END DO
  END DO
  
  CLOSE(1)  
!==================================
  END SUBROUTINE recover_genomes
!==================================  
!==================================
  SUBROUTINE recover_ARM
!==================================
  IMPLICIT NONE
  INTEGER :: i,j,k,rec_length
        
  INQUIRE(FILE="arm.r", EXIST=file_exists)
  IF(file_exists)THEN
    OPEN(UNIT=1,FILE='arm.r', status="old", action="read")
  ELSE
    WRITE(*,*)'ERROR, arm.r does not exist ', file_exists
    STOP
  END IF
  DO i=1,SIZE(newa)
    READ(1,'(F13.10)')newa(i)
  END DO
  CLOSE(1)
!==================================
  END SUBROUTINE recover_ARM
!==================================
!==================================
  SUBROUTINE map_recomb
!==================================
  IMPLICIT NONE
  INTEGER :: i
  
  DO i=1,nloci_sel
    IF(gen_d_map(i)%dista>0.0_dp)THEN ! it's a locus in-between two other loci
      gen_r_map(i)=kosambi_rec(gen_d_map(i)%dista)
    ELSE! it's at the start of a chromosome
      gen_r_map(i)=0.5_dp
    END IF
  END DO
!==================================
  END SUBROUTINE map_recomb
!==================================
!==================================
  SUBROUTINE mask_rec(r_vector, mask_vector, r_size)
! creates a mask of 0 and 1 indicating
! where recombination occurs in terms
! of nloci_sel intervals: the last
! interval is included because there
! is one last neutral locus
!==================================
  IMPLICIT NONE
  INTEGER :: i_mask
  INTEGER, INTENT(IN) :: r_size
  REAL(KIND=dp) :: xrand_mask
  REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: r_vector
  INTEGER(KIND=1), DIMENSION(:), INTENT(OUT) :: mask_vector
    
  mask_vector=0
  DO i_mask=1,r_size
    CALL RANDOM_NUMBER(xrand_mask)
    IF(xrand_mask<=r_vector(i_mask))THEN
      mask_vector(i_mask)=1
    ELSEIF(xrand_mask>r_vector(i_mask))THEN
      mask_vector(i_mask)=0
    END IF
  END DO
!==================================
  END SUBROUTINE mask_rec
!==================================  
!==================================
  SUBROUTINE branch_rec(in_xcross, out_branches, branch_size, branch_ini)
! transforms mask vector of mask_rec
! into a vector of branches indicating
! the parental haplotype according to
! recombination, and expands the action
! to intermediate neutral loci located
! in-between each pair of consecutive
! selected loci
!==================================
  IMPLICIT NONE
  INTEGER :: i_brec, branch_temp
  INTEGER, INTENT(IN) :: branch_ini
  INTEGER, INTENT(IN) :: branch_size
  INTEGER(KIND=1), DIMENSION(:), INTENT(IN) :: in_xcross
  INTEGER, DIMENSION(branch_size) :: temp_branches
  INTEGER, DIMENSION(:), INTENT(OUT) :: out_branches

  temp_branches=0
  temp_branches(1)=branch_ini
  branch_temp=branch_ini
  IF(branch_size>1)THEN
    DO i_brec=1,branch_size-1
      IF(in_xcross(i_brec)==0)THEN
        temp_branches(i_brec+1)=branch_temp
      ELSE IF(in_xcross(i_brec)==1)THEN
        IF(branch_temp==1)THEN
          branch_temp=2
        ELSE IF(branch_temp==2)THEN
          branch_temp=1
        END IF
        temp_branches(i_brec+1)=branch_temp
      END IF
    END DO
  END IF
! expand to neutrals
  counter=0
  out_branches=0
! store selected loci
  DO i_brec=1,branch_size*2,2
    counter=counter+1
    out_branches(i_brec)=temp_branches(counter)
  END DO
  counter=0
! interpolation of neutrals
  DO i_brec=2,branch_size*2,2
    counter=counter+1
    IF(i_brec<(branch_size*2))THEN
      IF(out_branches(i_brec-1)==out_branches(i_brec+1))THEN
        out_branches(i_brec)=out_branches(i_brec-1)
      ELSE IF(out_branches(i_brec-1)/=out_branches(i_brec+1))THEN
        branch_temp=which_branch(1,2)
        IF(branch_temp==1)THEN
          out_branches(i_brec)=out_branches(i_brec-1)
        ELSE IF(branch_temp==2)THEN
          out_branches(i_brec)=out_branches(i_brec+1)
        END IF
      END IF
    ELSE IF(i_brec==(branch_size*2))THEN
    ! this needs further rearrangement to account for recombination
      IF(in_xcross(branch_size)==1)THEN ! there is recombination between last selected and neutral
        branch_temp=which_branch(1,2)
        out_branches(i_brec)=branch_temp
      ELSE IF(in_xcross(branch_size)==0)THEN ! there is NO recombination between last selected and neutral
        out_branches(i_brec)=out_branches(i_brec-1)
      END IF
    END IF
  END DO
!==================================
  END SUBROUTINE branch_rec
!==================================
!==================================
  SUBROUTINE verify_Xcross(verx_sir, verx_dam, verx_vector)
!==================================
  IMPLICIT NONE
  INTEGER :: i_ver
  INTEGER, INTENT(IN) :: verx_sir,verx_dam
  INTEGER(KIND=1), DIMENSION(:), INTENT(IN) :: verx_vector
    
  DO i_ver=1,nloci_sel
    number_Xcross(verx_sir,verx_dam,i_ver)=&
    &number_Xcross(verx_sir,verx_dam,i_ver)+verx_vector(i_ver)
  END DO
!==================================
  END SUBROUTINE verify_Xcross
!==================================
!==================================
  END MODULE simpraise_routines
!==================================
!==================================
  PROGRAM simpraise
!==================================
  USE simpraise_param
  USE simpraise_stat
  USE simpraise_normal_table
  USE simpraise_routines
  IMPLICIT NONE
  SAVE

  INTEGER :: i,j
  OPEN(UNIT=31,FILE='insim.txt',STATUS='old')
  READ(31,'(A)') heading
  READ(31,'(A)') rname
  READ(31,'(A)') heading
  READ(31,'(A)') fname1
  READ(31,'(A)') heading
  READ(31,'(A)') fname2
!----------------------------------------------------
! mono_dioecious sets whether the population follows
! a monoic mating regime (mono_dioecious=0,
! both sexes within each individual and therefore
! selfin possible) or a dioecious mating regime
! (mono_dioecious=1, only one sexe per individual and
! no selfing is allowed)
!----------------------------------------------------
  READ(31,'(A)') heading
  READ(31,*) mono_dioecious
  READ(31,'(A)') heading
  READ(31,*) status_gener
  READ(31,'(A)') heading
  READ(31,*) mnum,fnum
  READ(31,'(A)') heading
  READ(31,*) num_cand
  READ(31,'(A)') heading
  READ(31,*) n_traits
  READ(31,'(A)') heading
  READ(31,*) nloci_X,nloci_XY,nloci_Y

  nloci_sel=nloci_X+nloci_Y+nloci_XY
  nloci_tot=nloci_sel+nloci_sel ! given that there are as many neutrals as selected

  IF(mono_dioecious==0)THEN
    ! adding +fnum et +mnum allows potentially for selfing
    ! it is at the elaboration of mating_design that selfing is to be avoided
    ALLOCATE(mating_design(mnum+fnum,mnum+fnum))
  ELSE IF(mono_dioecious==1)THEN
    ALLOCATE(mating_design(mnum,fnum))
  END IF ! mono_dioecious
  ALLOCATE(gene_type(nloci_sel),gene_effect(nloci_sel,4))
  ALLOCATE(freqini_vector(nloci_sel))
  freqini_vector=0.0_dp

  READ(31,'(A)') heading
  DO counter=1,nloci_sel
    READ(31,*) gene_type(counter),gene_effect(counter,1),&
    &gene_effect(counter,2),gene_effect(counter,3),gene_effect(counter,4)
  END DO
  READ(31,'(A)') heading
  DO counter=1,nloci_sel
    READ(31,*) freqini_vector(counter)
  END DO
  
  ALLOCATE(gen_d_map(nloci_sel))
  gen_d_map(:)%chrom=0
  gen_d_map(:)%dista=0.0_dp
  ALLOCATE(gen_r_map(nloci_sel))
  gen_r_map=0.0_dp
  ALLOCATE(sel_Xs_map(nloci_sel))
  sel_Xs_map=0
  ALLOCATE(tot_Xs_map(nloci_tot))
  tot_Xs_map=0
  
  IF(mono_dioecious==0)THEN
    ALLOCATE(number_Xcross(mnum+fnum,mnum+fnum,nloci_sel))
  ELSE IF(mono_dioecious==1)THEN
    ALLOCATE(number_Xcross(mnum,fnum,nloci_sel))
  END IF ! mono_dioecious
  number_Xcross=0
      
  READ(31,'(A)') heading
! first variable is chromosome, second is map distance within chromosome
  nchrom=0
  DO counter=1,nloci_sel
    READ(31,*) gen_d_map(counter)%chrom,gen_d_map(counter)%dista
    IF(gen_d_map(counter)%chrom>nchrom)THEN
      nchrom=gen_d_map(counter)%chrom
    END IF
  END DO
  
  CALL map_recomb
  
  tot_founder=mnum+fnum

  IF(num_cand<(mnum+fnum))THEN
    PRINT *, 'ERROR!!, less scored candidates than needed for parental replacement'
    STOP
  END IF
  IF((n_traits<1).OR.(n_traits>2))THEN
    PRINT *, 'ERROR!!, n_traits out of range'
    STOP
  END IF
  IF(MOD(num_cand,2)/=0)THEN
    PRINT *,'ERROR!!, non-even num_cand'
    STOP
  END IF
  IF(SUM(gene_type)/=(nloci_X+2*nloci_XY+3*nloci_Y))THEN
    PRINT *,'ERROR!!, number & type of gene do not match with READ'
    STOP
  END IF
  IF((n_traits==2).AND.(nloci_XY+nloci_Y==0))THEN
    PRINT *,'ERROR!!, THERE MUST BE genes for two traits'
    STOP
  END IF

  IF(mono_dioecious==1)THEN
    IF(MOD(fnum,mnum)/=0)THEN
      PRINT *,'ERROR!!, non-integer mating ratio'
      STOP
    ENDIF
    IF((MOD(fnum,2)/=0).AND.(MOD(mnum,2)/=0))THEN
      PRINT *,'ERROR!!, non-even initial census'
      STOP
    END IF
    IF(MOD(num_cand,fnum)/=0)THEN
      PRINT *,'ERROR!!, non-even num_cand/fnum'
      STOP
    END IF
    IF(MOD(num_cand,mnum)/=0)THEN
      PRINT *,'ERROR!!, non-even num_cand/mnum'
      STOP
    END IF
    IF(mnum+fnum<4)THEN
      PRINT *,'ERROR!!, census too small'
      STOP
    END IF
  END IF

  READ(31,'(A)') heading
  ! reads and stores genetic mating_design
  ! which contains also contributions per cross
  DO i=1,SIZE(mating_design,DIM=1) ! rows or females or maternal contributions
    DO j=1,SIZE(mating_design,DIM=2) ! columns or males or paternal contributions
      READ(31,*) mating_design(i,j)
    END DO
  END DO
  IF(SUM(mating_design)/=num_cand)THEN
    PRINT *,'ERROR!!, SUM(mating_design) not equal to num_cand'
    STOP
  END IF

!-----------------------------------------------
! be good and set and record a random number seed
! alternative version for compatibility GNU gfortran
  CALL system_clock(iseed(1))
  iseed(2)=-0.5*iseed(1)
  i=SIZE(iseed)
  CALL random_seed
  CALL random_seed(SIZE=i)
  CALL random_seed(PUT=iseed(1:i))
  CALL random_number(xdum)
  iseed(:)=INT(1.0e8*xdum(:))
  IF(1.0e8*xdum(1)-REAL(iseed(1),KIND(xdum))<0.5) iseed(1)=-iseed(1)
  IF(1.0e8*xdum(2)-REAL(iseed(2),KIND(xdum))>0.5) iseed(2)=-iseed(2)
  CALL random_seed(PUT=iseed(1:i))
!-----------------------------------------------
  
  OPEN(UNIT=21,FILE=fname1,STATUS='unknown')
  WRITE(21,'(A)') 'Random number seeds ... '
  WRITE(21,'(2I16)') iseed

  CALL write_input

  IF(status_gener==0)THEN
    ! ped contains founders + 1st num_cand descendants
    ALLOCATE(ped(tmax*num_cand+tot_founder))
  ELSE IF(status_gener==1)THEN
    ! ped must be dimensioned for previously accumulated records
    ! + those to be generated in current run
    INQUIRE(FILE="ped_size.r", EXIST=file_exists)
    IF(file_exists)THEN
      OPEN(UNIT=2,FILE='ped_size.r', status="old", action="read")
      READ(2,*)counter
      jlast=counter
      ALLOCATE(ped(jlast+tmax*num_cand))
      CLOSE(2)
    ELSE
      WRITE(*,*)'ERROR! file ped_size.r does not exist: ', file_exists
      STOP
    END IF
  END IF ! status_gener
  
  ALLOCATE(par_genome(num_cand,nloci_tot,2))
  ALLOCATE(des_genome(num_cand,nloci_tot,2))
  ALLOCATE(par_Genot(num_cand,nloci_sel,2))
  ALLOCATE(des_Genot(num_cand,nloci_sel,2))
  ALLOCATE(par_alcopies(num_cand,nloci_sel))
  ALLOCATE(des_alcopies(num_cand,nloci_sel))

  jini=0

  IF(MOD(num_cand*num_cand-num_cand,2)/=0)THEN
    PRINT *,'ERROR!!, non-even (num_cand*num_cand-num_cand)/2'
    STOP
  END IF

  dim_coa_mat=(num_cand*num_cand-num_cand)/2+num_cand
  ALLOCATE(newa(dim_coa_mat))
  ALLOCATE(olda(dim_coa_mat))
  ALLOCATE(newg(nloci_sel,dim_coa_mat))
  newa=0.0_dp
  olda=0.0_dp
  newg=0.0_dp

  ALLOCATE(allele_freq(tmax+1,nloci_sel),&
  & S_allele_freq(tmax+1,nloci_sel),&
  & Sq_allele_freq(tmax+1,nloci_sel))

  allele_freq=0.0_dp
  S_allele_freq=0.0_dp
  Sq_allele_freq=0.0_dp

  IF(mono_dioecious==0)THEN
    ALLOCATE(sel_sires(num_cand))
    ALLOCATE(sel_dams(num_cand))
  ELSE IF(mono_dioecious==1)THEN
    ALLOCATE(sel_sires(num_cand/2))
    ALLOCATE(sel_dams(num_cand/2))
  END IF

  ALLOCATE(Coamol(tmax+1),&
  & SCoamol(tmax+1),&
  & SqCoamol(tmax+1))
  Coamol=0.0_dp
  SCoamol=0.0_dp
  SqCoamol=0.0_dp

  ALLOCATE(Coa_newg(tmax+1),&
  & SCoa_newg(tmax+1),&
  & SqCoa_newg(tmax+1))
  Coa_newg=0.0_dp
  SCoa_newg=0.0_dp
  SqCoa_newg=0.0_dp

  ALLOCATE(Coaped(tmax+1),&
  & SCoaped(tmax+1),&
  & SqCoaped(tmax+1))
  Coaped=0.0_dp
  SCoaped=0.0_dp
  SqCoaped=0.0_dp

  ALLOCATE(Fbylocus(tmax,nloci_sel),&
  & SFbylocus(tmax,nloci_sel),&
  & SqFbylocus(tmax,nloci_sel))
  Fbylocus=0.0_dp
  SFbylocus=0.0_dp
  SqFbylocus=0.0_dp

  ALLOCATE(Fimol_tot(tmax+1),&
  & SFimol_tot(tmax+1),&
  & SqFimol_tot(tmax+1))
  Fimol_tot=0.0_dp
  SFimol_tot=0.0_dp
  SqFimol_tot=0.0_dp

  ALLOCATE(Fiped_tot(tmax+1),&
  & SFiped_tot(tmax+1),&
  & SqFiped_tot(tmax+1))
  Fiped_tot=0.0_dp
  SFiped_tot=0.0_dp
  SqFiped_tot=0.0_dp

  ALLOCATE(CovDesq_X11(tmax+1),&
  & SCovDesq_X11(tmax+1),&
  & SqCovDesq_X11(tmax+1))
  CovDesq_X11=0.0_dp
  SCovDesq_X11=0.0_dp
  SqCovDesq_X11=0.0_dp

  ALLOCATE(CovDesq_X12(tmax+1),&
  & SCovDesq_X12(tmax+1),&
  & SqCovDesq_X12(tmax+1))
  CovDesq_X12=0.0_dp
  SCovDesq_X12=0.0_dp
  SqCovDesq_X12=0.0_dp

  ALLOCATE(CovDesq_X22(tmax+1),&
  & SCovDesq_X22(tmax+1),&
  & SqCovDesq_X22(tmax+1))
  CovDesq_X22=0.0_dp
  SCovDesq_X22=0.0_dp
  SqCovDesq_X22=0.0_dp

  ALLOCATE(CovDesq_XY13(tmax+1),&
  & SCovDesq_XY13(tmax+1),&
  & SqCovDesq_XY13(tmax+1))
  CovDesq_XY13=0.0_dp
  SCovDesq_XY13=0.0_dp
  SqCovDesq_XY13=0.0_dp

  ALLOCATE(CovDesq_XY22(tmax+1),&
  & SCovDesq_XY22(tmax+1),&
  & SqCovDesq_XY22(tmax+1))
  CovDesq_XY22=0.0_dp
  SCovDesq_XY22=0.0_dp
  SqCovDesq_XY22=0.0_dp

  ALLOCATE(CovDesq_XY12(tmax+1),&
  & SCovDesq_XY12(tmax+1),&
  & SqCovDesq_XY12(tmax+1))
  CovDesq_XY12=0.0_dp
  SCovDesq_XY12=0.0_dp
  SqCovDesq_XY12=0.0_dp

  ALLOCATE(CovDesq_XY23(tmax+1),&
  & SCovDesq_XY23(tmax+1),&
  & SqCovDesq_XY23(tmax+1))
  CovDesq_XY23=0.0_dp
  SCovDesq_XY23=0.0_dp
  SqCovDesq_XY23=0.0_dp

  ALLOCATE(CovPlei_XY22(tmax+1),&
  & SCovPlei_XY22(tmax+1),&
  & SqCovPlei_XY22(tmax+1))
  CovPlei_XY22=0.0_dp
  SCovPlei_XY22=0.0_dp
  SqCovPlei_XY22=0.0_dp

  ALLOCATE(CovDesq_Y33(tmax+1),&
  & SCovDesq_Y33(tmax+1),&
  & SqCovDesq_Y33(tmax+1))
  CovDesq_Y33=0.0_dp
  SCovDesq_Y33=0.0_dp
  SqCovDesq_Y33=0.0_dp

  ALLOCATE(CovDesq_Y23(tmax+1),&
  & SCovDesq_Y23(tmax+1),&
  & SqCovDesq_Y23(tmax+1))
  CovDesq_Y23=0.0_dp
  SCovDesq_Y23=0.0_dp
  SqCovDesq_Y23=0.0_dp

  ALLOCATE(CovDesq_Y22(tmax+1),&
  & SCovDesq_Y22(tmax+1),&
  & SqCovDesq_Y22(tmax+1))
  CovDesq_Y22=0.0_dp
  SCovDesq_Y22=0.0_dp
  SqCovDesq_Y22=0.0_dp

  ALLOCATE(MGenot_X(tmax+1),&
  & SMGenot_X(tmax+1),&
  & SqMGenot_X(tmax+1))
  MGenot_X=0.0_dp
  SMGenot_X=0.0_dp
  SqMGenot_X=0.0_dp

  ALLOCATE(MGenot_Y(tmax+1),&
  & SMGenot_Y(tmax+1),&
  & SqMGenot_Y(tmax+1))
  MGenot_Y=0.0_dp
  SMGenot_Y=0.0_dp
  SqMGenot_Y=0.0_dp

  ALLOCATE(MGenot_xx(tmax+1),&
  & SMGenot_xx(tmax+1),&
  & SqMGenot_xx(tmax+1))
  MGenot_xx=0.0_dp
  SMGenot_xx=0.0_dp
  SqMGenot_xx=0.0_dp

  ALLOCATE(MGenot_yy(tmax+1),&
  & SMGenot_yy(tmax+1),&
  & SqMGenot_yy(tmax+1))
  MGenot_yy=0.0_dp
  SMGenot_yy=0.0_dp
  SqMGenot_yy=0.0_dp

  ALLOCATE(MGenot_xy(tmax+1),&
  & SMGenot_xy(tmax+1),&
  & SqMGenot_xy(tmax+1))
  MGenot_xy=0.0_dp
  SMGenot_xy=0.0_dp
  SqMGenot_xy=0.0_dp

  ALLOCATE(MGenot_yx(tmax+1),&
  & SMGenot_yx(tmax+1),&
  & SqMGenot_yx(tmax+1))
  MGenot_yx=0.0_dp
  SMGenot_yx=0.0_dp
  SqMGenot_yx=0.0_dp

  ALLOCATE(VaGenic_X(tmax+1),&
  & SVaGenic_X(tmax+1),&
  & SqVaGenic_X(tmax+1))
  VaGenic_X=0.0_dp
  SVaGenic_X=0.0_dp
  SqVaGenic_X=0.0_dp

  ALLOCATE(VaGenic_Y(tmax+1),&
  & SVaGenic_Y(tmax+1),&
  & SqVaGenic_Y(tmax+1))
  VaGenic_Y=0.0_dp
  SVaGenic_Y=0.0_dp
  SqVaGenic_Y=0.0_dp

  ALLOCATE(VdGenic_X(tmax+1),&
  & SVdGenic_X(tmax+1),&
  & SqVdGenic_X(tmax+1))
  VdGenic_X=0.0_dp
  SVdGenic_X=0.0_dp
  SqVdGenic_X=0.0_dp

  ALLOCATE(VdGenic_Y(tmax+1),&
  & SVdGenic_Y(tmax+1),&
  & SqVdGenic_Y(tmax+1))
  VdGenic_Y=0.0_dp
  SVdGenic_Y=0.0_dp
  SqVdGenic_Y=0.0_dp

  ALLOCATE(Vadd_X(tmax+1),&
  & SVadd_X(tmax+1),&
  & SqVadd_X(tmax+1))
  Vadd_X=0.0_dp
  SVadd_X=0.0_dp
  SqVadd_X=0.0_dp

  ALLOCATE(Vadd_Y(tmax+1),&
  & SVadd_Y(tmax+1),&
  & SqVadd_Y(tmax+1))
  Vadd_Y=0.0_dp
  SVadd_Y=0.0_dp
  SqVadd_Y=0.0_dp

  ALLOCATE(CovAdd_XY(tmax+1),&
  & SCovAdd_XY(tmax+1),&
  & SqCovAdd_XY(tmax+1))
  CovAdd_XY=0.0_dp
  SCovAdd_XY=0.0_dp
  SqCovAdd_XY=0.0_dp

  ALLOCATE(Vdom_X(tmax+1),&
  & SVdom_X(tmax+1),&
  & SqVdom_X(tmax+1))
  Vdom_X=0.0_dp
  SVdom_X=0.0_dp
  SqVdom_X=0.0_dp

  ALLOCATE(Vdom_Y(tmax+1),&
  & SVdom_Y(tmax+1),&
  & SqVdom_Y(tmax+1))
  Vdom_Y=0.0_dp
  SVdom_Y=0.0_dp
  SqVdom_Y=0.0_dp

  ALLOCATE(VGenot_X(tmax+1),&
  & SVGenot_X(tmax+1),&
  & SqVGenot_X(tmax+1))
  VGenot_X=0.0_dp
  SVGenot_X=0.0_dp
  SqVGenot_X=0.0_dp

  ALLOCATE(VGenot_Y(tmax+1),&
  & SVGenot_Y(tmax+1),&
  & SqVGenot_Y(tmax+1))
  VGenot_Y=0.0_dp
  SVGenot_Y=0.0_dp
  SqVGenot_Y=0.0_dp
  
  ALLOCATE(CovGenot_XY(tmax+1),&
  & SCovGenot_XY(tmax+1),&
  & SqCovGenot_XY(tmax+1))
  CovGenot_XY=0.0_dp
  SCovGenot_XY=0.0_dp
  SqCovGenot_XY=0.0_dp

  CALL DATE_AND_TIME(time=ctime)
  PRINT *, ctime,' : starting simulation ... '
    
! n_rep is set as parameter=1
  DO jrep=1,n_rep
    CALL do_a_rep
    CALL save_rep
    CALL DATE_AND_TIME(time=ctime)
    PRINT *, ctime,' : finished replicate ... ',jrep
  END DO

  CALL write_output1
  CLOSE(21)
  OPEN(UNIT=21,FILE=fname2,STATUS='unknown')
  CALL write_output2
  CLOSE(21) ! output files
  CLOSE(31) ! input files

  DEALLOCATE(ped)
  DEALLOCATE(newa,olda,newg)
  DEALLOCATE(par_genome)
  DEALLOCATE(des_genome)
  DEALLOCATE(par_Genot,des_Genot)
  DEALLOCATE(par_alcopies,des_alcopies)
  DEALLOCATE(freqini_vector)
  DEALLOCATE(sel_sires,sel_dams)
  DEALLOCATE(allele_freq,S_allele_freq,Sq_allele_freq)
  DEALLOCATE(gene_effect,gene_type)
  DEALLOCATE(MGenot_X,SMGenot_X,SqMGenot_X)
  DEALLOCATE(MGenot_Y,SMGenot_Y,SqMGenot_Y)
  DEALLOCATE(MGenot_xx,SMGenot_xx,SqMGenot_xx)
  DEALLOCATE(MGenot_yy,SMGenot_yy,SqMGenot_yy)
  DEALLOCATE(MGenot_xy,SMGenot_xy,SqMGenot_xy)
  DEALLOCATE(MGenot_yx,SMGenot_yx,SqMGenot_yx)
  DEALLOCATE(VGenot_X,SVGenot_X,SqVGenot_X)
  DEALLOCATE(VGenot_Y,SVGenot_Y,SqVGenot_Y)
  DEALLOCATE(CovGenot_XY,SCovGenot_XY,SqCovGenot_XY)
  DEALLOCATE(Vadd_X,SVadd_X,SqVadd_X)
  DEALLOCATE(Vadd_Y,SVadd_Y,SqVadd_Y)
  DEALLOCATE(Vdom_X,SVdom_X,SqVdom_X)
  DEALLOCATE(Vdom_Y,SVdom_Y,SqVdom_Y)
  DEALLOCATE(CovAdd_XY,SCovAdd_XY,SqCovAdd_XY)
  DEALLOCATE(VaGenic_X,SVaGenic_X,SqVaGenic_X)
  DEALLOCATE(VaGenic_Y,SVaGenic_Y,SqVaGenic_Y)
  DEALLOCATE(VdGenic_X,SVdGenic_X,SqVdGenic_X)
  DEALLOCATE(VdGenic_Y,SVdGenic_Y,SqVdGenic_Y)
  DEALLOCATE(CovDesq_X11,SCovDesq_X11,SqCovDesq_X11)
  DEALLOCATE(CovDesq_X12,SCovDesq_X12,SqCovDesq_X12)
  DEALLOCATE(CovDesq_X22,SCovDesq_X22,SqCovDesq_X22)
  DEALLOCATE(CovDesq_Y33,SCovDesq_Y33,SqCovDesq_Y33)
  DEALLOCATE(CovDesq_Y23,SCovDesq_Y23,SqCovDesq_Y23)
  DEALLOCATE(CovDesq_Y22,SCovDesq_Y22,SqCovDesq_Y22)
  DEALLOCATE(CovDesq_XY13,SCovDesq_XY13,SqCovDesq_XY13)
  DEALLOCATE(CovDesq_XY22,SCovDesq_XY22,SqCovDesq_XY22)
  DEALLOCATE(CovDesq_XY12,SCovDesq_XY12,SqCovDesq_XY12)
  DEALLOCATE(CovDesq_XY23,SCovDesq_XY23,SqCovDesq_XY23)
  DEALLOCATE(CovPlei_XY22,SCovPlei_XY22,SqCovPlei_XY22)
  DEALLOCATE(Coamol,SCoamol,SqCoamol)
  DEALLOCATE(Coa_newg,SCoa_newg,SqCoa_newg)
  DEALLOCATE(Coaped,SCoaped,SqCoaped)
  DEALLOCATE(Fbylocus,SFbylocus,SqFbylocus)
  DEALLOCATE(Fimol_tot,SFimol_tot,SqFimol_tot)
  DEALLOCATE(Fiped_tot,SFiped_tot,SqFiped_tot)
  DEALLOCATE(mating_design)
  DEALLOCATE(gen_d_map,gen_r_map,sel_Xs_map,tot_Xs_map)
  DEALLOCATE(number_Xcross)

  CONTAINS
!==================================
  SUBROUTINE do_a_rep
! LSR 15/01/2016
!==================================
  IMPLICIT NONE  
  
  IF(status_gener==0)THEN ! starts from founders
    time=0
    generation=0
    jini=1
    CALL setup_base
    CALL backup_genomes ! of founders
    CALL sumup_founderGenot
    CALL allelic_freq
    CALL update_var_and_means
    CALL base_coa_reset
    time=1
    generation=1
    CALL create_offspring
    CALL sumup_offspringGenot
    CALL allelic_freq
    CALL update_var_and_means
    CALL update_ARM
    CALL update_molmat
    CALL rank_candidates
    CALL new_genome
    CALL backup_genomes ! of descendants
    CALL backup_ped
  ELSE IF(status_gener==1)THEN ! starts from previous descendants
    time=1
    CALL recover_genomes
    CALL recover_ped
    CALL create_offspring
    CALL sumup_offspringGenot
    CALL allelic_freq
    CALL update_var_and_means
    CALL update_ARM
    CALL update_molmat
    CALL rank_candidates
    CALL new_genome
    CALL backup_genomes
    CALL backup_ped
  END IF
!==================================
  END SUBROUTINE do_a_rep
!==================================
!==================================
  SUBROUTINE save_rep
!==================================
  IMPLICIT NONE
  INTEGER :: i,j,k

    DO i=1,tmax+1
      SMGenot_X(i)=SMGenot_X(i)+MGenot_X(i)
      SqMGenot_X(i)=SqMGenot_X(i)+(MGenot_X(i)*MGenot_X(i))
      SMGenot_Y(i)=SMGenot_Y(i)+MGenot_Y(i)
      SqMGenot_Y(i)=SqMGenot_Y(i)+(MGenot_Y(i)*MGenot_Y(i))
      SMGenot_xx(i)=SMGenot_xx(i)+MGenot_xx(i)
      SqMGenot_xx(i)=SqMGenot_xx(i)+(MGenot_xx(i)*MGenot_xx(i))
      SMGenot_yy(i)=SMGenot_yy(i)+MGenot_yy(i)
      SqMGenot_yy(i)=SqMGenot_yy(i)+(MGenot_yy(i)*MGenot_yy(i))
      SMGenot_xy(i)=SMGenot_xy(i)+MGenot_xy(i)
      SqMGenot_xy(i)=SqMGenot_xy(i)+(MGenot_xy(i)*MGenot_xy(i))
      SMGenot_yx(i)=SMGenot_yx(i)+MGenot_yx(i)
      SqMGenot_yx(i)=SqMGenot_yx(i)+(MGenot_yx(i)*MGenot_yx(i))
      SVGenot_X(i)=SVGenot_X(i)+VGenot_X(i)
      SqVGenot_X(i)=SqVGenot_X(i)+(VGenot_X(i)*VGenot_X(i))
      SVGenot_Y(i)=SVGenot_Y(i)+VGenot_Y(i)
      SqVGenot_Y(i)=SqVGenot_Y(i)+(VGenot_Y(i)*VGenot_Y(i))
      SCovGenot_XY(i)=SCovGenot_XY(i)+CovGenot_XY(i)
      SqCovGenot_XY(i)=SqCovGenot_XY(i)+(CovGenot_XY(i)*CovGenot_XY(i))
      SVadd_X(i)=SVadd_X(i)+Vadd_X(i)
      SqVadd_X(i)=SqVadd_X(i)+(Vadd_X(i)*Vadd_X(i))
      SVadd_Y(i)=SVadd_Y(i)+Vadd_Y(i)
      SqVadd_Y(i)=SqVadd_Y(i)+(Vadd_Y(i)*Vadd_Y(i))
      SCovAdd_XY(i)=SCovAdd_XY(i)+CovAdd_XY(i)
      SqCovAdd_XY(i)=SqCovAdd_XY(i)+(CovAdd_XY(i)*CovAdd_XY(i))
      SVdom_X(i)=SVdom_X(i)+Vdom_X(i)
      SqVdom_X(i)=SqVdom_X(i)+(Vdom_X(i)*Vdom_X(i))
      SVdom_Y(i)=SVdom_Y(i)+Vdom_Y(i)
      SqVdom_Y(i)=SqVdom_Y(i)+(Vdom_Y(i)*Vdom_Y(i))
      SVaGenic_X(i)=SVaGenic_X(i)+VaGenic_X(i)
      SqVaGenic_X(i)=SqVaGenic_X(i)+(VaGenic_X(i)*VaGenic_X(i))
      SVaGenic_Y(i)=SVaGenic_Y(i)+VaGenic_Y(i)
      SqVaGenic_Y(i)=SqVaGenic_Y(i)+(VaGenic_Y(i)*VaGenic_Y(i))
      SVdGenic_X(i)=SVdGenic_X(i)+VdGenic_X(i)
      SqVdGenic_X(i)=SqVdGenic_X(i)+(VdGenic_X(i)*VdGenic_X(i))
      SVdGenic_Y(i)=SVdGenic_Y(i)+VdGenic_Y(i)
      SqVdGenic_Y(i)=SqVdGenic_Y(i)+(VdGenic_Y(i)*VdGenic_Y(i))
      SCovDesq_X11(i)=SCovDesq_X11(i)+CovDesq_X11(i)
      SqCovDesq_X11(i)=SqCovDesq_X11(i)+(CovDesq_X11(i)*CovDesq_X11(i))
      SCovDesq_X12(i)=SCovDesq_X12(i)+CovDesq_X12(i)
      SqCovDesq_X12(i)=SqCovDesq_X12(i)+(CovDesq_X12(i)*CovDesq_X12(i))
      SCovDesq_X22(i)=SCovDesq_X22(i)+CovDesq_X22(i)
      SqCovDesq_X22(i)=SqCovDesq_X22(i)+(CovDesq_X22(i)*CovDesq_X22(i))
      SCovDesq_Y33(i)=SCovDesq_Y33(i)+CovDesq_Y33(i)
      SqCovDesq_Y33(i)=SqCovDesq_Y33(i)+(CovDesq_Y33(i)*CovDesq_Y33(i))
      SCovDesq_Y23(i)=SCovDesq_Y23(i)+CovDesq_Y23(i)
      SqCovDesq_Y23(i)=SqCovDesq_Y23(i)+(CovDesq_Y23(i)*CovDesq_Y23(i))
      SCovDesq_Y22(i)=SCovDesq_Y22(i)+CovDesq_Y22(i)
      SqCovDesq_Y22(i)=SqCovDesq_Y22(i)+(CovDesq_Y22(i)*CovDesq_Y22(i))
      SCovDesq_XY13(i)=SCovDesq_XY13(i)+CovDesq_XY13(i)
      SqCovDesq_XY13(i)=SqCovDesq_XY13(i)+(CovDesq_XY13(i)*CovDesq_XY13(i))
      SCovDesq_XY22(i)=SCovDesq_XY22(i)+CovDesq_XY22(i)
      SqCovDesq_XY22(i)=SqCovDesq_XY22(i)+(CovDesq_XY22(i)*CovDesq_XY22(i))
      SCovDesq_XY12(i)=SCovDesq_XY12(i)+CovDesq_XY12(i)
      SqCovDesq_XY12(i)=SqCovDesq_XY12(i)+(CovDesq_XY12(i)*CovDesq_XY12(i))
      SCovDesq_XY23(i)=SCovDesq_XY23(i)+CovDesq_XY23(i)
      SqCovDesq_XY23(i)=SqCovDesq_XY23(i)+(CovDesq_XY23(i)*CovDesq_XY23(i))
      SCovPlei_XY22(i)=SCovPlei_XY22(i)+CovPlei_XY22(i)
      SqCovPlei_XY22(i)=SqCovPlei_XY22(i)+(CovPlei_XY22(i)*CovPlei_XY22(i))
      SFiped_tot(i)=SFiped_tot(i)+Fiped_tot(i)
      SqFiped_tot(i)=SqFiped_tot(i)+(Fiped_tot(i)*Fiped_tot(i))
      SFimol_tot(i)=SFimol_tot(i)+Fimol_tot(i)
      SqFimol_tot(i)=SqFimol_tot(i)+(Fimol_tot(i)*Fimol_tot(i))
      SCoaped(i)=SCoaped(i)+Coaped(i)
      SqCoaped(i)=SqCoaped(i)+(Coaped(i)*Coaped(i))
      SCoamol(i)=SCoamol(i)+Coamol(i)
      SqCoamol(i)=SqCoamol(i)+(Coamol(i)*Coamol(i))
      SCoa_newg(i)=SCoa_newg(i)+Coa_newg(i)
      SqCoa_newg(i)=SqCoa_newg(i)+(Coa_newg(i)*Coa_newg(i))
    END DO

    DO i=1,tmax
      DO j=1,nloci_sel
        IF(Fbylocus(i,j)>0.0_dp)THEN
          SFbylocus(i,j)=SFbylocus(i,j)+Fbylocus(i,j)
          SqFbylocus(i,j)=SqFbylocus(i,j)+(Fbylocus(i,j)*Fbylocus(i,j))
        END IF
      END DO
    END DO

  DO i=1,tmax+1
    DO j=1,nloci_sel
      S_allele_freq(i,j)=S_allele_freq(i,j)+allele_freq(i,j)
      Sq_allele_freq(i,j)=Sq_allele_freq(i,j)+(allele_freq(i,j)*allele_freq(i,j))
    END DO
  END DO
!==================================
  END SUBROUTINE save_rep
!==================================
!==================================
  SUBROUTINE write_input
!==================================
  IMPLICIT NONE

  WRITE (21,'(A,T50)') 'simpraise'
  WRITE (21,'(A,T50)') 'leopoldo.sanchez-rodriguez@inra.fr'
  WRITE (21,*)' '
  WRITE (21,'(A,A)') 'FILE ',fname1
  WRITE (21,'(A,A)') 'NAME_OF_THIS_RUN ',rname
  WRITE (21,*) ' '
  WRITE (21,'(A)') 'INPUT_PARAMETERS'
  WRITE (21,'(A)') '************************************************'
  IF(mono_dioecious==0)THEN
    WRITE (21,'(A)') '* monoecious'
  ELSE IF(mono_dioecious==1)THEN
    WRITE (21,'(A)') '* dioecious'
  END IF
  WRITE (21,'(A,T50,I6)') 'Total_#_traits: ',n_traits
  WRITE (21,'(A,T50,I6)') 'Total_#_selected_loci: ',nloci_sel
  WRITE (21,'(A,T50,I6)') 'Loci_trait_X: ',nloci_X
  WRITE (21,'(A,T50,I6)') 'Pleiotropic_loci: ',nloci_XY
  WRITE (21,'(A,T50,I6)') 'Loci_trait_Y: ',nloci_Y
  WRITE (21,'(A,T50,I6)') 'Tmax: ',tmax
  WRITE (21,'(A,T53,10(I4,1X))') '#_male_parents: ',mnum
  WRITE (21,'(A,T53,10(I4,1X))') '#_female_parents: ',fnum
  WRITE (21,'(A,T53,10(I4,1X))') '#_scored_candidates: ',num_cand
  WRITE (21,'(A,T50,I6)') '#_replicates: ',n_rep
  WRITE (21,'(A)') '************************************************'
  WRITE (21,*) ' '
!==================================
  END SUBROUTINE write_input
!==================================
!==================================
  SUBROUTINE write_output1
!==================================
  IMPLICIT NONE
  INTEGER :: i,j,k

    DO i=1,tmax+1
      MGenot_X(i)=sem_rep(SqMGenot_X(i),SMGenot_X(i),n_rep)
      SqMGenot_X(i)=mean_rep(SMGenot_X(i),n_rep)
      MGenot_Y(i)=sem_rep(SqMGenot_Y(i),SMGenot_Y(i),n_rep)
      SqMGenot_Y(i)=mean_rep(SMGenot_Y(i),n_rep)
      MGenot_xx(i)=sem_rep(SqMGenot_xx(i),SMGenot_xx(i),n_rep)
      SqMGenot_xx(i)=mean_rep(SMGenot_xx(i),n_rep)
      MGenot_yy(i)=sem_rep(SqMGenot_yy(i),SMGenot_yy(i),n_rep)
      SqMGenot_yy(i)=mean_rep(SMGenot_yy(i),n_rep)
      MGenot_xy(i)=sem_rep(SqMGenot_xy(i),SMGenot_xy(i),n_rep)
      SqMGenot_xy(i)=mean_rep(SMGenot_xy(i),n_rep)
      MGenot_yx(i)=sem_rep(SqMGenot_yx(i),SMGenot_yx(i),n_rep)
      SqMGenot_yx(i)=mean_rep(SMGenot_yx(i),n_rep)
      VGenot_X(i)=sem_rep(SqVGenot_X(i),SVGenot_X(i),n_rep)
      SqVGenot_X(i)=mean_rep(SVGenot_X(i),n_rep)
      VGenot_Y(i)=sem_rep(SqVGenot_Y(i),SVGenot_Y(i),n_rep)
      SqVGenot_Y(i)=mean_rep(SVGenot_Y(i),n_rep)
      CovGenot_XY(i)=sem_rep(SqCovGenot_XY(i),SCovGenot_XY(i),n_rep)
      SqCovGenot_XY(i)=mean_rep(SCovGenot_XY(i),n_rep)
      CovAdd_XY(i)=sem_rep(SqCovAdd_XY(i),SCovAdd_XY(i),n_rep)
      SqCovAdd_XY(i)=mean_rep(SCovAdd_XY(i),n_rep)
      VaGenic_X(i)=sem_rep(SqVaGenic_X(i),SVaGenic_X(i),n_rep)
      SqVaGenic_X(i)=mean_rep(SVaGenic_X(i),n_rep)
      VaGenic_Y(i)=sem_rep(SqVaGenic_Y(i),SVaGenic_Y(i),n_rep)
      SqVaGenic_Y(i)=mean_rep(SVaGenic_Y(i),n_rep)
      VdGenic_X(i)=sem_rep(SqVdGenic_X(i),SVdGenic_X(i),n_rep)
      SqVdGenic_X(i)=mean_rep(SVdGenic_X(i),n_rep)
      VdGenic_Y(i)=sem_rep(SqVdGenic_Y(i),SVdGenic_Y(i),n_rep)
      SqVdGenic_Y(i)=mean_rep(SVdGenic_Y(i),n_rep)
      Vadd_X(i)=sem_rep(SqVadd_X(i),SVadd_X(i),n_rep)
      SqVadd_X(i)=mean_rep(SVadd_X(i),n_rep)
      Vadd_Y(i)=sem_rep(SqVadd_Y(i),SVadd_Y(i),n_rep)
      SqVadd_Y(i)=mean_rep(SVadd_Y(i),n_rep)
      Vdom_X(i)=sem_rep(SqVdom_X(i),SVdom_X(i),n_rep)
      SqVdom_X(i)=mean_rep(SVdom_X(i),n_rep)
      Vdom_Y(i)=sem_rep(SqVdom_Y(i),SVdom_Y(i),n_rep)
      SqVdom_Y(i)=mean_rep(SVdom_Y(i),n_rep)
      CovDesq_X11(i)=sem_rep(SqCovDesq_X11(i),SCovDesq_X11(i),n_rep)
      SqCovDesq_X11(i)=mean_rep(SCovDesq_X11(i),n_rep)
      CovDesq_X12(i)=sem_rep(SqCovDesq_X12(i),SCovDesq_X12(i),n_rep)
      SqCovDesq_X12(i)=mean_rep(SCovDesq_X12(i),n_rep)
      CovDesq_X22(i)=sem_rep(SqCovDesq_X22(i),SCovDesq_X22(i),n_rep)
      SqCovDesq_X22(i)=mean_rep(SCovDesq_X22(i),n_rep)
      CovDesq_Y33(i)=sem_rep(SqCovDesq_Y33(i),SCovDesq_Y33(i),n_rep)
      SqCovDesq_Y33(i)=mean_rep(SCovDesq_Y33(i),n_rep)
      CovDesq_Y23(i)=sem_rep(SqCovDesq_Y23(i),SCovDesq_Y23(i),n_rep)
      SqCovDesq_Y23(i)=mean_rep(SCovDesq_Y23(i),n_rep)
      CovDesq_Y22(i)=sem_rep(SqCovDesq_Y22(i),SCovDesq_Y22(i),n_rep)
      SqCovDesq_Y22(i)=mean_rep(SCovDesq_Y22(i),n_rep)
      CovDesq_XY13(i)=sem_rep(SqCovDesq_XY13(i),SCovDesq_XY13(i),n_rep)
      SqCovDesq_XY13(i)=mean_rep(SCovDesq_XY13(i),n_rep)
      CovDesq_XY22(i)=sem_rep(SqCovDesq_XY22(i),SCovDesq_XY22(i),n_rep)
      SqCovDesq_XY22(i)=mean_rep(SCovDesq_XY22(i),n_rep)
      CovDesq_XY12(i)=sem_rep(SqCovDesq_XY12(i),SCovDesq_XY12(i),n_rep)
      SqCovDesq_XY12(i)=mean_rep(SCovDesq_XY12(i),n_rep)
      CovDesq_XY23(i)=sem_rep(SqCovDesq_XY23(i),SCovDesq_XY23(i),n_rep)
      SqCovDesq_XY23(i)=mean_rep(SCovDesq_XY23(i),n_rep)
      CovPlei_XY22(i)=sem_rep(SqCovPlei_XY22(i),SCovPlei_XY22(i),n_rep)
      SqCovPlei_XY22(i)=mean_rep(SCovPlei_XY22(i),n_rep)
      Fiped_tot(i)=sem_rep(SqFiped_tot(i),SFiped_tot(i),n_rep)
      SqFiped_tot(i)=mean_rep(SFiped_tot(i),n_rep)
      Fimol_tot(i)=sem_rep(SqFimol_tot(i),SFimol_tot(i),n_rep)
      SqFimol_tot(i)=mean_rep(SFimol_tot(i),n_rep)
      Coaped(i)=sem_rep(SqCoaped(i),SCoaped(i),n_rep)
      SqCoaped(i)=mean_rep(SCoaped(i),n_rep)
      Coamol(i)=sem_rep(SqCoamol(i),SCoamol(i),n_rep)
      SqCoamol(i)=mean_rep(SCoamol(i),n_rep)
      Coa_newg(i)=sem_rep(SqCoa_newg(i),SCoa_newg(i),n_rep)
      SqCoa_newg(i)=mean_rep(SCoa_newg(i),n_rep)
    END DO
    DO i=1,tmax
      DO j=1,nloci_sel
        Fbylocus(i,j)=sem_rep(SqFbylocus(i,j),SFbylocus(i,j),n_rep)
        SqFbylocus(i,j)=mean_rep(SFbylocus(i,j),n_rep)
      END DO
    END DO
  DO i=1,tmax+1
    DO j=1,nloci_sel
      allele_freq(i,j)=sem_rep(Sq_allele_freq(i,j),S_allele_freq(i,j),n_rep)
      Sq_allele_freq(i,j)=mean_rep(S_allele_freq(i,j),n_rep)
    END DO
  END DO

  WRITE (21,'(A)') 'OUTPUT VARIABLES'
  WRITE (21,'(A)') '************************************************'
  WRITE (21,'(A)') ' '

    WRITE (21,'(A)') 'Genotypic_mean_X'
    DO i=1,tmax+1
      WRITE(21,*)SqMGenot_X(i),MGenot_X(i)
    END DO
    WRITE (21,'(A)') 'Genotypic_mean_Y'
    DO i=1,tmax+1
      WRITE(21,*)SqMGenot_Y(i),MGenot_Y(i)
    END DO
    WRITE (21,'(A)') 'Genotypic_mean_xx'
    DO i=1,tmax+1
      WRITE(21,*)SqMGenot_xx(i),MGenot_xx(i)
    END DO
    WRITE (21,'(A)') 'Genotypic_mean_yy'
    DO i=1,tmax+1
      WRITE(21,*)SqMGenot_yy(i),MGenot_yy(i)
    END DO
    WRITE (21,'(A)') 'Genotypic_mean_xy'
    DO i=1,tmax+1
      WRITE(21,*)SqMGenot_xy(i),MGenot_xy(i)
    END DO
    WRITE (21,'(A)') 'Genotypic_mean_yx'
    DO i=1,tmax+1
      WRITE(21,*)SqMGenot_yx(i),MGenot_yx(i)
    END DO
    WRITE (21,'(A)') 'Genotypic_variance_X'
    DO i=1,tmax+1
      WRITE(21,*)SqVGenot_X(i),VGenot_X(i)
    END DO
    WRITE (21,'(A)') 'Genotypic_variance_Y'
    DO i=1,tmax+1
      WRITE(21,*)SqVGenot_Y(i),VGenot_Y(i)
    END DO
    WRITE (21,'(A)') 'Genotypic_covariance_XY'
    DO i=1,tmax+1
      WRITE(21,*)SqCovGenot_XY(i),CovGenot_XY(i)
    END DO
    WRITE (21,'(A)') 'Additive_variance_X'
    DO i=1,tmax+1
      WRITE(21,*)SqVadd_X(i),Vadd_X(i)
    END DO
    WRITE (21,'(A)') 'Additive_variance_Y'
    DO i=1,tmax+1
      WRITE(21,*)SqVadd_Y(i),Vadd_Y(i)
    END DO
    WRITE (21,'(A)') 'Additive_covariance_XY'
    DO i=1,tmax+1
      WRITE(21,*)SqCovAdd_XY(i),CovAdd_XY(i)
    END DO
    WRITE (21,'(A)') 'Dominance_variance_X'
    DO i=1,tmax+1
      WRITE(21,*)SqVdom_X(i),Vdom_X(i)
    END DO
    WRITE (21,'(A)') 'Dominance_variance_Y'
    DO i=1,tmax+1
      WRITE(21,*)SqVdom_Y(i),Vdom_Y(i)
    END DO
    WRITE (21,'(A)') 'Additive_Genic_X'
    DO i=1,tmax+1
      WRITE(21,*)SqVaGenic_X(i),VaGenic_X(i)
    END DO
    WRITE (21,'(A)') 'Additive_Genic_Y'
    DO i=1,tmax+1
      WRITE(21,*)SqVaGenic_Y(i),VaGenic_Y(i)
    END DO
    WRITE (21,'(A)') 'Dominance_Genic_X'
    DO i=1,tmax+1
      WRITE(21,*)SqVdGenic_X(i),VdGenic_X(i)
    END DO
    WRITE (21,'(A)') 'Dominance_Genic_Y'
    DO i=1,tmax+1
      WRITE(21,*)SqVdGenic_Y(i),VdGenic_Y(i)
    END DO
    WRITE (21,'(A)') 'LD_covariance_X11'
    DO i=1,tmax+1
      WRITE(21,*)SqCovDesq_X11(i),CovDesq_X11(i)
    END DO
    WRITE (21,'(A)') 'LD_covariance_X12'
    DO i=1,tmax+1
      WRITE(21,*)SqCovDesq_X12(i),CovDesq_X12(i)
    END DO
    WRITE (21,'(A)') 'LD_covariance_X22'
    DO i=1,tmax+1
      WRITE(21,*)SqCovDesq_X22(i),CovDesq_X22(i)
    END DO
    WRITE (21,'(A)') 'LD_covariance_Y33'
    DO i=1,tmax+1
      WRITE(21,*)SqCovDesq_Y33(i),CovDesq_Y33(i)
    END DO
    WRITE (21,'(A)') 'LD_covariance_Y23'
    DO i=1,tmax+1
      WRITE(21,*)SqCovDesq_Y23(i),CovDesq_Y23(i)
    END DO
    WRITE (21,'(A)') 'LD_covariance_Y22'
    DO i=1,tmax+1
      WRITE(21,*)SqCovDesq_Y22(i),CovDesq_Y22(i)
    END DO
    WRITE (21,'(A)') 'LD_covariance_XY13'
    DO i=1,tmax+1
      WRITE(21,*)SqCovDesq_XY13(i),CovDesq_XY13(i)
    END DO
    WRITE (21,'(A)') 'LD_covariance_XY22'
    DO i=1,tmax+1
      WRITE(21,*)SqCovDesq_XY22(i),CovDesq_XY22(i)
    END DO
    WRITE (21,'(A)') 'LD_covariance_XY12'
    DO i=1,tmax+1
      WRITE(21,*)SqCovDesq_XY12(i),CovDesq_XY12(i)
    END DO
    WRITE (21,'(A)') 'LD_covariance_XY23'
    DO i=1,tmax+1
      WRITE(21,*)SqCovDesq_XY23(i),CovDesq_XY23(i)
    END DO
    WRITE (21,'(A)') 'Pleiotropic_covariance_XY22'
    DO i=1,tmax+1
      WRITE(21,*)SqCovPlei_XY22(i),CovPlei_XY22(i)
    END DO
    WRITE (21,'(A)') 'Fiped_tot'
    DO i=1,tmax+1
      WRITE(21,*)SqFiped_tot(i),Fiped_tot(i)
    END DO
    WRITE (21,'(A)') 'Fimol_tot'
    DO i=1,tmax+1
      WRITE(21,*)SqFimol_tot(i),Fimol_tot(i)
    END DO
    WRITE (21,'(A)') 'Coaped'
    DO i=1,tmax+1
      WRITE(21,*)SqCoaped(i),Coaped(i)
    END DO
    WRITE (21,'(A)') 'Coamol'
    DO i=1,tmax+1
      WRITE(21,*)SqCoamol(i),Coamol(i)
    END DO
    WRITE (21,'(A)') 'Coa_newg'
    DO i=1,tmax+1
      WRITE(21,*)SqCoa_newg(i),Coa_newg(i)
    END DO
  !==================================
  END SUBROUTINE write_output1
!==================================
!==================================
  SUBROUTINE write_output2
!==================================
  IMPLICIT NONE
  INTEGER :: i,j,k
  REAL(KIND=dp) :: sum_across
  REAL(KIND=dp), DIMENSION(nchrom) :: s_xcross_chrom

  WRITE (21,'(A,T50)') 'simpraise'
  WRITE (21,'(A,T50)') 'leopoldo.sanchez-rodriguez@inra.fr'
  WRITE (21,*)' '
  WRITE (21,'(A,A)') 'FILE ',fname2
  WRITE (21,*)' '
  WRITE (21,'(A,A)') 'NAME_OF_THIS_RUN ',rname
  WRITE (21,*) ' '
  WRITE (21,'(A,T50,I6)') 'Tmax: ',tmax
  WRITE (21,'(A,T50,I6)') 'Nloci_neu: ',nloci_sel
  WRITE (21,'(A,T50,I6)') 'Nloci_sel: ',nloci_sel

  WRITE (21,'(A)') 'inbreeding/pop/gener/locus'
    DO i=1,tmax
      DO j=1,nloci_sel
        WRITE(21,*)SqFbylocus(i,j),Fbylocus(i,j)
      END DO
    END DO

  WRITE (21,'(A)') 'allele_freqs/gener/locus'
  DO i=1,tmax+1
    DO j=1,nloci_sel
      WRITE(21,*)Sq_allele_freq(i,j),allele_freq(i,j)
    END DO
  END DO

  WRITE (21,'(A)') 'cross-overs'
  DO i=1,nloci_sel
    WRITE(21,*) i, gen_d_map(i)%dista
  END DO
  sum_across=0.0_dp
  s_xcross_chrom=0.0_dp
  DO i=1,SIZE(mating_design,DIM=1) ! rows
    DO j=1,SIZE(mating_design,DIM=2) ! columns
      IF(mating_design(i,j)>0)THEN
        DO k=1, nloci_sel
          sum_across=REAL(number_Xcross(i,j,k),dp)/REAL((2*mating_design(i,j)),dp)
          WRITE(21,*) k, sum_across
          IF(gen_d_map(k)%dista>0.0_dp)THEN
            s_xcross_chrom(gen_d_map(k)%chrom)=s_xcross_chrom(gen_d_map(k)%chrom)+REAL(number_Xcross(i,j,k),dp)
          END IF
        END DO
      END IF
    END DO
  END DO
  WRITE (21,'(A)') 'cross-overs/chromosome'
  sum_across=0.0_dp
  DO i=1,nchrom
    sum_across=s_xcross_chrom(i)/(2*num_cand)
    WRITE(21,*) i,sum_across
  END DO
!==================================
  END SUBROUTINE write_output2
!==================================
!==================================
  END PROGRAM simpraise
!==================================