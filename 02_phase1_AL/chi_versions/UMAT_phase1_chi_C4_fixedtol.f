!DEC$ FREEFORM
!
MODULE Abaqus_Interface
  !INCLUDE 'ABA_PARAM.INC'
  PRIVATE
  INTEGER,PARAMETER,PUBLIC::rk=KIND(0.D0),ik=KIND(1)
END MODULE Abaqus_Interface
!     
! Include this module for explicit declaration of variables and interface
! definition for DEBUG purpose 
! dummy interface for debugging without ABAQUS -> DEBUG
! MODULE ABAQUS_INTERFACE
!   PRIVATE
!   INTEGER,PARAMETER,PUBLIC::rk=KIND(0.D0),ik=KIND(1)
! END MODULE
!     
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$                  
MODULE Abaqus_Interface_extended
! subroutine interface to compute eigenvalues of stress
! (Abaqus utility routine)
  INTERFACE
     SUBROUTINE SPRINC(S,PS,LSTR,NDI,NSHR)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk)::S(*),PS(*)
       INTEGER(KIND=ik)::LSTR,NDI,NSHR
     END SUBROUTINE SPRINC
  END INTERFACE
! subroutine interface to compute stress invariants
! (Abaqus utility routine)
  INTERFACE 
     SUBROUTINE SINV(S,SINV1,SINV2,NDI,NSHR)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk)::S(*),SINV1,SINV2
       INTEGER(KIND=ik)::NDI,NSHR
     END SUBROUTINE SINV
  END INTERFACE
! subroutine interface to rotate stress and strain tensor
! (Abaqus utility routine)
  INTERFACE
     SUBROUTINE ROTSIG(S,R,SPRIME,LSTR,NDI,NSHR)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk)::S(*),R(*),SPRIME(*)
       INTEGER(KIND=ik)::LSTR,NDI,NSHR
     END SUBROUTINE ROTSIG
  END INTERFACE
END MODULE Abaqus_Interface_extended
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
MODULE LAPACK_utility_functions
  INTERFACE GESV 
    SUBROUTINE SGESV( N, NRHS, A, LDA, PIV, B, LDB, INFO )
      INTEGER, INTENT(IN) :: LDA, LDB, NRHS, N
      INTEGER, INTENT(OUT) :: INFO
      INTEGER, INTENT(OUT) :: PIV(*)
      REAL(KIND(1.0E0)), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
    END SUBROUTINE
    SUBROUTINE DGESV( N, NRHS, A, LDA, PIV, B, LDB, INFO )
      INTEGER, INTENT(IN) :: LDA, LDB, NRHS, N
      INTEGER, INTENT(OUT) :: INFO
      INTEGER, INTENT(OUT) :: PIV(*)
      REAL(KIND(1.0D0)), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
    END SUBROUTINE
 END INTERFACE GESV
 INTERFACE GELSS
    SUBROUTINE SGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, &
         & WORK, LWORK, INFO)
      INTEGER, INTENT(IN) :: M, N, NRHS, LDA, LDB, LWORK
      INTEGER, INTENT(OUT) :: RANK, INFO
      REAL(KIND(1.0E0)), INTENT(IN) :: RCOND
      REAL(KIND(1.0E0)), INTENT(OUT) :: WORK(*), S(*)
      REAL(KIND(1.0E0)), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
    END SUBROUTINE SGELSS
    SUBROUTINE DGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, &
         & WORK, LWORK, INFO)
      INTEGER, INTENT(IN) :: M, N, NRHS, LDA, LDB, LWORK
      INTEGER, INTENT(OUT) :: RANK, INFO
      REAL(KIND(1.0D0)), INTENT(IN) :: RCOND
      REAL(KIND(1.0D0)), INTENT(OUT) :: WORK(*), S(*)
      REAL(KIND(1.0D0)), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
    END SUBROUTINE DGELSS
 END INTERFACE GELSS
END MODULE LAPACK_utility_functions
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
MODULE Constants
  USE Abaqus_Interface
  INTEGER(KIND=ik),PARAMETER::ONE_int=1_ik, TWO_int=2_ik, THREE_int=3_ik
  INTEGER(KIND=ik),PARAMETER::FOUR_int=4_ik, FIVE_int=5_ik, SIX_int=6_ik
  INTEGER(KIND=ik),PARAMETER::EIGHT_int=8_ik, NINE_int=9_ik, ZERO_int=0_ik
  INTEGER(KIND=ik),PARAMETER::TEN_int=10_ik,TWELVE_int=12_ik
  INTEGER(KIND=ik),PARAMETER::TWENTY_THREE_int=23_ik,TWENTY_FOUR_int=24_ik
  !
  REAL(KIND=rk),PARAMETER::ZERO_r=0.0_rk,ONE_r=1.0_rk,TWO_r=2.0_rk
  REAL(KIND=rk),PARAMETER::THREE_r=3.0_rk,FOUR_r=4.0_rk,SIX_r=6.0_rk
  REAL(KIND=rk),PARAMETER::FIVE_r=5.0_rk,NINE_r=9.0_rk,SIXTEEN_r=16.0_rk
  REAL(KIND=rk),PARAMETER::ONE_by_FIVE_r=ONE_r/FIVE_r
  REAL(KIND=rk),PARAMETER::NINE_by_SIXTEEN_r=NINE_r/SIXTEEN_r
  REAL(KIND=rk),PARAMETER::THREE_by_TWO_r=THREE_r/TWO_r
  REAL(KIND=rk),PARAMETER::ONE_HALF_r=ONE_r/TWO_r,ONE_TENTH_r=0.1_rk
  REAL(KIND=rk),PARAMETER::ONE_by_THREE_r=ONE_r/THREE_r
  REAL(KIND=rk),PARAMETER::TWO_by_THREE_r=TWO_r/THREE_r
  REAL(KIND=rk),PARAMETER::ONE_FOURTH_r=ONE_r/FOUR_r
  REAL(KIND=rk),PARAMETER::ONE_HUNDRED_EIGHTY_r=180.0_rk
  REAL(KIND=rk),PARAMETER::ONE_HUNDREDTH_r=0.01_rk
  REAL(KIND=rk),PARAMETER::sqrt_two_r=sqrt(TWO_r),ZERO_p_SIX_r=0.6_rk
  REAL(KIND=rk),PARAMETER::sqrt_three_r=sqrt(THREE_r)
  REAL(KIND=rk),PARAMETER::ONE_by_sqrt_two_r=ONE_r/sqrt_two_r
  REAL(KIND=rk),PARAMETER::ONE_by_sqrt_three_r=ONE_r/sqrt_three_r
  REAL(KIND=rk),PARAMETER::PI=3.14159265358979323846264338327950288_rk
  !  
  INTEGER(KIND=rk),DIMENSION(9,2),PARAMETER::INDEX_NINE_FEAP=&
       &RESHAPE((/1,2,3,1,2,2,3,1,3,1,2,3,2,1,3,2,3,1/),(/9,2/))
  INTEGER(KIND=rk),DIMENSION(9,2),PARAMETER::INDEX_NINE_ABAQ=&
       &RESHAPE((/1,2,3,1,2,1,3,2,3,1,2,3,2,1,3,1,3,2/),(/9,2/))
  INTEGER(KIND=rk),DIMENSION(6,2),PARAMETER::INDEX_SIX_FEAP=&
       &RESHAPE((/1,2,3,1,2,1,1,2,3,2,3,3/),(/6,2/))
  INTEGER(KIND=rk),DIMENSION(6,2),PARAMETER::INDEX_SIX_ABAQ=&
       &RESHAPE((/1,2,3,1,1,2,1,2,3,2,3,3/),(/6,2/))
  INTEGER(KIND=rk),DIMENSION(3,3),PARAMETER::INDEX_TENS_VEC_FEAP=&
       &RESHAPE((/1,5,9,4,2,7,8,6,3/),(/3,3/))
  INTEGER(KIND=rk),DIMENSION(3,3),PARAMETER::INDEX_TENS_VEC_FEAP_S=&
       &RESHAPE((/1,4,6,4,2,5,6,5,3/),(/3,3/))
  INTEGER(KIND=rk),DIMENSION(3,3),PARAMETER::INDEX_TENS_VEC_ABAQ=&
       &RESHAPE((/1,5,7,4,2,9,6,8,3/),(/3,3/))
  INTEGER(KIND=rk),DIMENSION(3,3),PARAMETER::INDEX_TENS_VEC_ABAQ_S=&
       &RESHAPE((/1,4,5,4,2,6,5,6,3/),(/3,3/))
  ! real-valued parameters defining tolerances in iterative loop
  REAL(KIND=rk),PARAMETER::eps_iter=(EPSILON(ZERO_r))**(ONE_HALF_r)
  REAL(KIND=rk),PARAMETER::zero_distinct=eps_iter**(THREE_by_TWO_r)
  REAL(KIND=rk),PARAMETER::chi=1.0E-5_rk            ! Eq.40 χ
  REAL(KIND=rk),PARAMETER::eps_min_nw=zero_distinct  ! Eq.39/40 ε_min
  REAL(KIND=rk),PARAMETER::tol_kkt=1.0E-10_rk        ! ε_c, Alg.2 line 8
  REAL(KIND=rk),PARAMETER::beta_eta=10.0_rk          ! Eq.36 β
  REAL(KIND=rk),PARAMETER::gamma_eta=0.25_rk         ! Eq.36 γ
  REAL(KIND=rk),PARAMETER::tiny_zero=eps_iter**THREE_int
  REAL(KIND=rk),PARAMETER::tol_nw_red=(EPSILON(ZERO_r))**(NINE_by_SIXTEEN_r)
  REAL(KIND=rk),PARAMETER::rel_step_size=TWO_r**(-TWENTY_THREE_int)
  REAL(KIND=rk),PARAMETER::up_bound_relax_time=ONE_r/eps_iter
  REAL(KIND=rk),PARAMETER::lb_invert=zero_distinct**(ONE_HALF_r)
  REAL(KIND=rk),PARAMETER::alpha_bt=eps_iter**(ONE_HALF_r)
  REAL(KIND=rk),PARAMETER::eps_one=EPSILON(ZERO_r)
  REAL(KIND=rk),PARAMETER::tol_angle=eps_iter/ONE_HUNDREDTH_r
  REAL(KIND=rk),DIMENSION(3,3),PARAMETER::IDENT=RESHAPE((/ONE_r,&
       &ZERO_r,ZERO_r,ZERO_r,ONE_r,ZERO_r,ZERO_r,ZERO_r,ONE_r/),(/3,3/))
END MODULE Constants
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
MODULE auxiliary_routines_interface
  INTERFACE arrange_condense_4_order_tensor
     SUBROUTINE arrange_condense_4_order_tensor(A_4ot,A_4ot_9b9,A_4ot_6b6,&
          &output_type)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(3,3,3,3),INTENT(IN)::A_4ot
       REAL(KIND=rk),DIMENSION(9,9),INTENT(OUT)::A_4ot_9b9
       REAL(KIND=rk),DIMENSION(6,6),INTENT(OUT)::A_4ot_6b6
       CHARACTER(3),INTENT(IN)::output_type
     END SUBROUTINE arrange_condense_4_order_tensor
  END INTERFACE arrange_condense_4_order_tensor
  !
  INTERFACE back_tracking_step_length_scaling_parameter
     SUBROUTINE back_tracking_step_length_scaling_parameter(back_tracking,&
          &back_track_count,lambda,lambda_two,x,x_save,dx,obj_f,obj_f_save,&
          &obj_f_two,grad_f_save,ndim)
       USE Abaqus_Interface
       IMPLICIT NONE
       INTEGER(KIND=ik),INTENT(IN)::ndim
       INTEGER(KIND=ik),INTENT(INOUT)::back_track_count
       REAL(KIND=rk),INTENT(INOUT)::lambda,lambda_two,obj_f_two
       REAL(KIND=rk),INTENT(IN)::obj_f,obj_f_save
       REAL(KIND=rk),DIMENSION(ndim),INTENT(INOUT)::x
       REAL(KIND=rk),DIMENSION(ndim),INTENT(IN)::x_save,dx,grad_f_save
       LOGICAL,INTENT(INOUT)::back_tracking
     END SUBROUTINE back_tracking_step_length_scaling_parameter
  END INTERFACE back_tracking_step_length_scaling_parameter
  !
  INTERFACE M33INV
     SUBROUTINE M33INV (A, AINV, OK_FLAG)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(3,3), INTENT(IN)  :: A
       REAL(KIND=rk),DIMENSION(3,3), INTENT(OUT) :: AINV
       LOGICAL, INTENT(OUT) :: OK_FLAG
     END SUBROUTINE M33INV
  END INTERFACE M33INV
  !
  INTERFACE pushrc
     SUBROUTINE pushrc(S,AA,S0,AA0,DFGRD,iaa,env)
       USE Abaqus_Interface
       IMPLICIT NONE
       INTEGER(KIND=ik),INTENT(IN)::iaa
       REAL(KIND=rk),DIMENSION(6),INTENT(IN)::S0
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::DFGRD
       REAL(KIND=rk),DIMENSION(6,6),INTENT(IN)::AA0
       REAL(KIND=rk),DIMENSION(6),INTENT(OUT)::S
       REAL(KIND=rk),DIMENSION(6,6),INTENT(OUT)::AA
       CHARACTER(4),INTENT(IN)::env
     END SUBROUTINE pushrc
  END INTERFACE pushrc
  !
  INTERFACE transform_symm_tens_vec_form
     SUBROUTINE transform_symm_tens_vec_form(A_tens,A_tens_vec,transform_id,&
          &env)
       USE Abaqus_Interface
       IMPLICIT NONE
       INTEGER(KIND=ik),INTENT(IN)::transform_id
       REAL(KIND=rk),DIMENSION(3,3)::A_tens
       REAL(KIND=rk),DIMENSION(6)::A_tens_vec
       CHARACTER(4),INTENT(IN)::env
     END SUBROUTINE transform_symm_tens_vec_form
  END INTERFACE transform_symm_tens_vec_form
  !
  INTERFACE transform_tens_vec_form
     SUBROUTINE transform_tens_vec_form(A_tens,A_tens_vec,transform_id,&
          &env)
       USE Abaqus_Interface
       IMPLICIT NONE
       INTEGER(KIND=ik),INTENT(IN)::transform_id
       REAL(KIND=rk),DIMENSION(3,3)::A_tens
       REAL(KIND=rk),DIMENSION(9)::A_tens_vec
       CHARACTER(4),INTENT(IN)::env
     END SUBROUTINE transform_tens_vec_form
  END INTERFACE transform_tens_vec_form
! ----------------------------------------------------------------------
  INTERFACE approp_index
     FUNCTION approp_index(ndim,env,funct_name)
       USE Abaqus_Interface
       IMPLICIT NONE
       INTEGER(KIND=ik),INTENT(IN)::ndim
       CHARACTER(4),INTENT(IN)::env
       CHARACTER(*),INTENT(IN)::funct_name
       INTEGER(KIND=ik),DIMENSION(ndim,2)::approp_index
     END FUNCTION approp_index
  END INTERFACE approp_index
  !
  INTERFACE comp_det_Matrix
     FUNCTION comp_det_Matrix(Matrix)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::Matrix
       REAL(KIND=rk)::comp_det_Matrix
     END FUNCTION comp_det_Matrix
  END INTERFACE comp_det_Matrix
  !
  INTERFACE comp_sec_ord_tens_rank_one
     FUNCTION comp_sec_ord_tens_rank_one(a_vec,b_vec)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(3),INTENT(IN)::a_vec,b_vec
       REAL(KIND=rk),DIMENSION(3,3)::comp_sec_ord_tens_rank_one
     END FUNCTION comp_sec_ord_tens_rank_one
  END INTERFACE comp_sec_ord_tens_rank_one
  !
  INTERFACE compute_trace_second_order_tensor
     FUNCTION compute_trace_second_order_tensor(A_tens)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::A_tens
       REAL(KIND=rk)::compute_trace_second_order_tensor
     END FUNCTION compute_trace_second_order_tensor
  END INTERFACE compute_trace_second_order_tensor
  !
  INTERFACE double_contraction_fourth_order_second_order_tensor_matmul
     FUNCTION double_contraction_fourth_order_second_order_tensor_matmul &
          &(A_fot_6x6,b_sot_6x1)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(6,6),INTENT(IN)::A_fot_6x6
       REAL(KIND=rk),DIMENSION(6),INTENT(IN)::b_sot_6x1
       REAL(KIND=rk),DIMENSION(6)::double_contraction_fourth_order_second_order_tensor_matmul
     END FUNCTION double_contraction_fourth_order_second_order_tensor_matmul
  END INTERFACE double_contraction_fourth_order_second_order_tensor_matmul
  !
  INTERFACE double_contraction_fourth_order_tensor_matmul
     FUNCTION double_contraction_fourth_order_tensor_matmul(A_fot,B_fot)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(6,6),INTENT(IN)::A_fot,B_fot
       REAL(KIND=rk),DIMENSION(6,6)::double_contraction_fourth_order_tensor_matmul
     END FUNCTION double_contraction_fourth_order_tensor_matmul
  END INTERFACE double_contraction_fourth_order_tensor_matmul
  !
  INTERFACE double_contraction_fourth_order_tensor_second_order_tensor
     FUNCTION double_contraction_fourth_order_tensor_second_order_tensor(&
          &A4_tensor,B2_tensor)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::B2_tensor
       REAL(KIND=rk),DIMENSION(3,3,3,3),INTENT(IN)::A4_tensor
       REAL(KIND=rk),DIMENSION(3,3)::double_contraction_fourth_order_tensor_second_order_tensor
     END FUNCTION double_contraction_fourth_order_tensor_second_order_tensor
  END INTERFACE double_contraction_fourth_order_tensor_second_order_tensor
  !
  INTERFACE double_contraction_fourth_order_tensor
     FUNCTION double_contraction_fourth_order_tensor(A4_tensor,B4_tensor)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(3,3,3,3),INTENT(IN)::A4_tensor,B4_tensor
       REAL(KIND=rk),DIMENSION(3,3,3,3)::double_contraction_fourth_order_tensor
     END FUNCTION double_contraction_fourth_order_tensor
  END INTERFACE double_contraction_fourth_order_tensor
  !
  INTERFACE double_contraction_second_order_tensors
     FUNCTION double_contraction_second_order_tensors(A_tens,B_tens)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::A_tens,B_tens
       REAL(KIND=rk)::double_contraction_second_order_tensors
     END FUNCTION double_contraction_second_order_tensors
  END INTERFACE double_contraction_second_order_tensors
  !
  INTERFACE DYADIC_PRODUCT_SECOND_ORDER_TENSORS
     FUNCTION DYADIC_PRODUCT_SECOND_ORDER_TENSORS(A_tens,B_tens)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::A_tens,B_tens
       REAL(KIND=rk),DIMENSION(3,3,3,3)::DYADIC_PRODUCT_SECOND_ORDER_TENSORS
     END FUNCTION DYADIC_PRODUCT_SECOND_ORDER_TENSORS
  END INTERFACE DYADIC_PRODUCT_SECOND_ORDER_TENSORS
  !
  INTERFACE EulAngFromRotMat
     FUNCTION EulAngFromRotMat(RotMat)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::RotMat
       REAL(KIND=rk),DIMENSION(3)::EulAngFromRotMat
     END FUNCTION EulAngFromRotMat
  END INTERFACE EulAngFromRotMat
  !
  INTERFACE expand_fourth_order_tensor_matrix2tensor
     FUNCTION expand_fourth_order_tensor_matrix2tensor(A_tens_mat,env)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(9,9),INTENT(IN)::A_tens_mat
       REAL(KIND=rk),DIMENSION(3,3,3,3)::expand_fourth_order_tensor_matrix2tensor
       CHARACTER(4),INTENT(IN)::env
     END FUNCTION expand_fourth_order_tensor_matrix2tensor
  END INTERFACE expand_fourth_order_tensor_matrix2tensor
  !
  INTERFACE expand_minor_symm_fourth_order_tensor_matrix2tensor
     FUNCTION expand_minor_symm_fourth_order_tensor_matrix2tensor(A6b6,&
          &env)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(6,6),INTENT(IN)::A6b6
       REAL(KIND=rk),DIMENSION(3,3,3,3)::expand_minor_symm_fourth_order_tensor_matrix2tensor
       CHARACTER(4),INTENT(IN)::env
     END FUNCTION expand_minor_symm_fourth_order_tensor_matrix2tensor
  END INTERFACE expand_minor_symm_fourth_order_tensor_matrix2tensor
  !
  INTERFACE extract_symm_fourth_ord_tens_minor_symm
     FUNCTION extract_symm_fourth_ord_tens_minor_symm(A_tens_9b9)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(9,9),INTENT(IN)::A_tens_9b9
       REAL(KIND=rk),DIMENSION(6,6)::extract_symm_fourth_ord_tens_minor_symm
     END FUNCTION extract_symm_fourth_ord_tens_minor_symm
  END INTERFACE extract_symm_fourth_ord_tens_minor_symm
  !
  INTERFACE left_double_contraction_fourth_order_tensor_second_order_tensor
     FUNCTION left_double_contraction_fourth_order_tensor_second_order_tensor(&
          &A4_tensor,B2_tensor)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(3,3,3,3),INTENT(IN)::A4_tensor
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::B2_tensor
       REAL(KIND=rk),DIMENSION(3,3)::left_double_contraction_fourth_order_tensor_second_order_tensor
     END FUNCTION left_double_contraction_fourth_order_tensor_second_order_tensor
  END INTERFACE left_double_contraction_fourth_order_tensor_second_order_tensor
  !
  INTERFACE map_fourth_order_tensor_2_matrix_format
     FUNCTION map_fourth_order_tensor_2_matrix_format(A_tens,env)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(3,3,3,3),INTENT(IN)::A_tens
       REAL(KIND=rk),DIMENSION(9,9)::map_fourth_order_tensor_2_matrix_format
       CHARACTER(4),INTENT(IN)::env
     END FUNCTION map_fourth_order_tensor_2_matrix_format
  END INTERFACE map_fourth_order_tensor_2_matrix_format
  !
  INTERFACE pertub_F_numtang_material_tang
     FUNCTION pertub_F_numtang_material_tang(F_inv,pertub_ind)
       USE Abaqus_Interface
       IMPLICIT NONE
       INTEGER(KIND=ik),DIMENSION(2),INTENT(IN)::pertub_ind
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::F_inv
       REAL(KIND=rk),DIMENSION(3,3)::pertub_F_numtang_material_tang
     END FUNCTION pertub_F_numtang_material_tang
  END INTERFACE pertub_F_numtang_material_tang
  !
  INTERFACE polar_decompostion_def_grad
     FUNCTION polar_decompostion_def_grad(DefGrad,NDI,NSHR,env)
       USE Abaqus_Interface
       IMPLICIT NONE
       INTEGER(KIND=ik),INTENT(IN)::NDI,NSHR
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::DefGrad
       REAL(KIND=rk),DIMENSION(3,3)::polar_decompostion_def_grad
       CHARACTER(4),INTENT(IN)::env
     END FUNCTION polar_decompostion_def_grad
  END INTERFACE polar_decompostion_def_grad
  !
  INTERFACE ROTATION_MATRIX_BUNGE
     FUNCTION ROTATION_MATRIX_BUNGE(angles)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(3),INTENT(IN)::angles
       REAL(KIND=rk),DIMENSION(3,3)::ROTATION_MATRIX_BUNGE
     END FUNCTION ROTATION_MATRIX_BUNGE
  END INTERFACE ROTATION_MATRIX_BUNGE
  !
  INTERFACE SOLVE_SYSTEM_OF_EQUATION
     FUNCTION SOLVE_SYSTEM_OF_EQUATION(Mat,rhs,nrhs,prob_dim)
       USE Abaqus_Interface
       IMPLICIT NONE
       INTEGER(KIND=ik),INTENT(IN)::nrhs,prob_dim
       REAL(KIND=rk),DIMENSION(prob_dim,prob_dim),INTENT(IN)::Mat! prob_dim,prob_dim
       REAL(KIND=rk),DIMENSION(prob_dim,nrhs),INTENT(IN)::rhs ! prob_dim,nrhs
       REAL(KIND=rk),DIMENSION(prob_dim,nrhs)::SOLVE_SYSTEM_OF_EQUATION ! prob_dim,nrhs
     END FUNCTION SOLVE_SYSTEM_OF_EQUATION
  END INTERFACE SOLVE_SYSTEM_OF_EQUATION
  !
  INTERFACE SOLVE_ILL_COND_SYSTEM_OF_EQUATION
     FUNCTION SOLVE_ILL_COND_SYSTEM_OF_EQUATION(Mat,rhs,nrhs,prob_dim,tol)
       USE Abaqus_Interface
       IMPLICIT NONE
       INTEGER(KIND=ik),INTENT(IN)::nrhs,prob_dim
       REAL(KIND=rk),INTENT(IN)::tol
       REAL(KIND=rk),DIMENSION(prob_dim,prob_dim),INTENT(IN)::Mat
       REAL(KIND=rk),DIMENSION(prob_dim,nrhs),INTENT(IN)::rhs
       REAL(KIND=rk),DIMENSION(prob_dim,nrhs)::SOLVE_ILL_COND_SYSTEM_OF_EQUATION
     END FUNCTION SOLVE_ILL_COND_SYSTEM_OF_EQUATION
  END INTERFACE SOLVE_ILL_COND_SYSTEM_OF_EQUATION
  !
  INTERFACE spec_tensor_product_second_order_tensors
     FUNCTION spec_tensor_product_second_order_tensors(A_tens,B_tens,&
          &env)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::A_tens,B_tens
       REAL(KIND=rk),DIMENSION(6,6)::spec_tensor_product_second_order_tensors
       CHARACTER(4),INTENT(IN)::env
     END FUNCTION spec_tensor_product_second_order_tensors
  END INTERFACE spec_tensor_product_second_order_tensors
  !
  INTERFACE spec_tensor_product_unsym_second_order_tensors
     FUNCTION spec_tensor_product_unsym_second_order_tensors(A_tens,B_tens,&
          &env)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::A_tens,B_tens
       REAL(KIND=rk),DIMENSION(9,9)::spec_tensor_product_unsym_second_order_tensors
       CHARACTER(4),INTENT(IN)::env
     END FUNCTION spec_tensor_product_unsym_second_order_tensors
  END INTERFACE spec_tensor_product_unsym_second_order_tensors
  !
  INTERFACE sq_prod_unsym_sec_ord_tens_dble_contr_fourth_order_sym_Ident
     FUNCTION sq_prod_unsym_sec_ord_tens_dble_contr_fourth_order_sym_Ident(A,&
          &B,env)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::A,B
       REAL(KIND=rk),DIMENSION(9,9)::sq_prod_unsym_sec_ord_tens_dble_contr_fourth_order_sym_Ident
       CHARACTER(4),INTENT(IN)::env
     END FUNCTION sq_prod_unsym_sec_ord_tens_dble_contr_fourth_order_sym_Ident
  END INTERFACE sq_prod_unsym_sec_ord_tens_dble_contr_fourth_order_sym_Ident
  !
  INTERFACE square_product_second_order_tensors
     FUNCTION square_product_second_order_tensors(A_tens,B_tens,env)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::A_tens,B_tens
       REAL(KIND=rk),DIMENSION(6,6)::square_product_second_order_tensors
       CHARACTER(4),INTENT(IN)::env
     END FUNCTION square_product_second_order_tensors
  END INTERFACE square_product_second_order_tensors
  !
  INTERFACE square_product_unsym_second_order_tensors
     FUNCTION square_product_unsym_second_order_tensors(A_tens,B_tens,&
          &env)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::A_tens,B_tens
       REAL(KIND=rk),DIMENSION(9,9)::square_product_unsym_second_order_tensors
       CHARACTER(4),INTENT(IN)::env
     END FUNCTION square_product_unsym_second_order_tensors
  END INTERFACE square_product_unsym_second_order_tensors
  !
  INTERFACE symm_tensor_product_second_order_tensors
     FUNCTION symm_tensor_product_second_order_tensors(A_tens,B_tens,&
          &env)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::A_tens,B_tens
       REAL(KIND=rk),DIMENSION(6,6)::symm_tensor_product_second_order_tensors
       CHARACTER(4),INTENT(IN)::env
     END FUNCTION symm_tensor_product_second_order_tensors
  END INTERFACE symm_tensor_product_second_order_tensors
END MODULE auxiliary_routines_interface
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
MODULE constitutive_routines_interface
  INTERFACE adjust_visco_parameter_miehe_vp_log
     SUBROUTINE adjust_visco_parameter_miehe_vp_log (p_exp_vp_t,&
          &tau_relax_vp_t,mode,p_exp_vp_target,tau_relax_vp_target)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),INTENT(IN)::p_exp_vp_target,tau_relax_vp_target
       REAL(KIND=rk),INTENT(INOUT)::p_exp_vp_t,tau_relax_vp_t
       CHARACTER(17),INTENT(IN)::mode
     END SUBROUTINE adjust_visco_parameter_miehe_vp_log
  END INTERFACE adjust_visco_parameter_miehe_vp_log
  !
  INTERFACE assembly_residual_jacobian_newton
     SUBROUTINE assembly_residual_jacobian_newton (residual_nw,&
          & jacobian_nw,active_slip_system,deriv_visco_plast_potent,&
          & Phi_alpha,dPhi_alpha_d_delta_lambda_ip1_beta,delta_lambda_i,&
          & delta_lambda_ip1,plastic_slip,algo,dt,eta,alpha_slip,&
          & mat_param_inelast,theta_KK_t,hard_law)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(24),INTENT(INOUT)::residual_nw
       REAL(KIND=rk),DIMENSION(24),INTENT(INOUT)::deriv_visco_plast_potent
       REAL(KIND=rk),DIMENSION(24,24),INTENT(INOUT)::jacobian_nw
       REAL(KIND=rk),INTENT(IN)::Phi_alpha,dt,eta,theta_KK_t
       REAL(KIND=rk),DIMENSION(13),INTENT(IN)::mat_param_inelast
       REAL(KIND=rk),DIMENSION(24),INTENT(IN)::dPhi_alpha_d_delta_lambda_ip1_beta
       REAL(KIND=rk),DIMENSION(24),INTENT(IN)::delta_lambda_i,delta_lambda_ip1
       REAL(KIND=rk),DIMENSION(24),INTENT(IN)::plastic_slip
       LOGICAL,DIMENSION(24),INTENT(INOUT)::active_slip_system
       INTEGER(KIND=ik),INTENT(IN)::algo,alpha_slip,hard_law
     END SUBROUTINE assembly_residual_jacobian_newton
  END INTERFACE assembly_residual_jacobian_newton
  !
  INTERFACE check_eigval_distinct
     SUBROUTINE check_eigval_distinct(n_distinct,index_distinct,lambda,&
          &tol)
       USE Abaqus_Interface
       IMPLICIT NONE
       INTEGER(KIND=ik),INTENT(OUT)::n_distinct
       INTEGER(KIND=ik),DIMENSION(3),INTENT(OUT)::index_distinct
       REAL(KIND=rk),DIMENSION(3),INTENT(IN)::lambda
       REAL(KIND=rk),INTENT(IN)::tol
     END SUBROUTINE check_eigval_distinct
  END INTERFACE check_eigval_distinct
  !
  INTERFACE coefficients_for_projection
     SUBROUTINE coefficients_for_projection(theta,xi,eta,e_vect,d_vect,&
          &lambda,index_distinct,n_distinct)
       USE Abaqus_Interface
       IMPLICIT NONE
       INTEGER(KIND=ik),INTENT(IN)::n_distinct
       INTEGER(KIND=ik),DIMENSION(3),INTENT(IN)::index_distinct
       REAL(KIND=rk),DIMENSION(3),INTENT(IN)::lambda,e_vect,d_vect
       REAL(KIND=rk),DIMENSION(3,3),INTENT(OUT)::theta,xi
       REAL(KIND=rk),INTENT(OUT)::eta
     END SUBROUTINE coefficients_for_projection
  END INTERFACE coefficients_for_projection
  !
  INTERFACE crystal_plasticity_augm_lagr
     SUBROUTINE crystal_plasticity_augm_lagr(S2PK,C_tang,SDV,F,dt,&
          &elast_const,mat_param_inelast,m_SH,symm_class,orth_basis,&
          &orth_proj,algo,numtang,NDI,NSHR,NSTATV,env,hard_law,theta_KK,&
          &conv_issue_dil)
       USE Abaqus_Interface
       IMPLICIT NONE
       INTEGER(KIND=ik),INTENT(IN)::algo,symm_class,NDI,NSHR,numtang
       INTEGER(KIND=ik),INTENT(IN)::hard_law,NSTATV
       REAL(KIND=rk),DIMENSION(3,3),INTENT(OUT)::S2PK
       REAL(KIND=rk),DIMENSION(6,6),INTENT(OUT)::C_tang
       REAL(KIND=rk),DIMENSION(NSTATV),INTENT(INOUT)::SDV
       REAL(KIND=rk),INTENT(IN)::dt,m_SH,theta_KK
       REAL(KIND=rk),DIMENSION(3),INTENT(IN)::elast_const
       REAL(KIND=rk),DIMENSION(13),INTENT(IN)::mat_param_inelast
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::F
       REAL(KIND=rk),DIMENSION(6,3,3),INTENT(IN)::orth_basis
       REAL(KIND=rk),DIMENSION(6,6,6),INTENT(IN)::orth_proj
       LOGICAL,INTENT(OUT)::conv_issue_dil
       CHARACTER(4),INTENT(IN)::env
     END SUBROUTINE crystal_plasticity_augm_lagr
  END INTERFACE crystal_plasticity_augm_lagr
  !
  INTERFACE determine_incremental_lagrange_multiplier_augmented_lagrange
     SUBROUTINE determine_incremental_lagrange_multiplier_augmented_lagrange(&
          &F_p_inv_iter,plastic_slip_iter,d_delta_lambda_ip1_d_C,&
          &active_slip_system,unimod_C,F_p_inv_n,C_inv,Je,plastic_slip_n,&
          &algo,dt,elast_const,mat_param_inelast,m_SH,symm_class,&
          &orth_basis,orth_proj,NDI,NSHR,env,hard_law,theta_KK,&
          &conv_issue_max_ati)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(24),INTENT(OUT)::plastic_slip_iter
       REAL(KIND=rk),DIMENSION(3,3),INTENT(OUT)::F_p_inv_iter
       REAL(KIND=rk),DIMENSION(24,6),INTENT(OUT)::d_delta_lambda_ip1_d_C
       REAL(KIND=rk),INTENT(IN)::Je,dt,m_SH,theta_KK
       REAL(KIND=rk),DIMENSION(3),INTENT(IN)::elast_const
       REAL(KIND=rk),DIMENSION(13),INTENT(IN)::mat_param_inelast
       REAL(KIND=rk),DIMENSION(24),INTENT(IN)::plastic_slip_n
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::unimod_C,F_p_inv_n,C_inv
       REAL(KIND=rk),DIMENSION(6,3,3),INTENT(IN)::orth_basis
       REAL(KIND=rk),DIMENSION(6,6,6),INTENT(IN)::orth_proj
       INTEGER(KIND=ik),INTENT(IN)::algo,symm_class,NDI,NSHR,hard_law
       LOGICAL,DIMENSION(24),INTENT(OUT)::active_slip_system
       LOGICAL,INTENT(OUT)::conv_issue_max_ati
       CHARACTER(4),INTENT(IN)::env
     END SUBROUTINE determine_incremental_lagrange_multiplier_augmented_lagrange
  END INTERFACE determine_incremental_lagrange_multiplier_augmented_lagrange
  !
  INTERFACE determine_incremental_lagrange_multiplier_miehe_viscoplastic
     SUBROUTINE determine_incremental_lagrange_multiplier_miehe_viscoplastic(&
          &F_p_inv_iter,plastic_slip_iter,d_delta_lambda_ip1_d_C,&
          &active_slip_system,unimod_C,F_p_inv_n,C_inv,Je,plastic_slip_n,&
          &algo,dt,elast_const,mat_param_inelast,m_SH,symm_class,&
          &orth_basis,orth_proj,NDI,NSHR,env,hard_law,theta_KK,&
          &aver_adjust_per_step,conv_issue)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(24),INTENT(OUT)::plastic_slip_iter
       REAL(KIND=rk),DIMENSION(3,3),INTENT(OUT)::F_p_inv_iter
       REAL(KIND=rk),DIMENSION(24,6),INTENT(OUT)::d_delta_lambda_ip1_d_C
       REAL(KIND=rk),INTENT(IN)::Je,dt,m_SH,theta_KK
       REAL(KIND=rk),INTENT(INOUT)::aver_adjust_per_step
       REAL(KIND=rk),DIMENSION(3),INTENT(IN)::elast_const
       REAL(KIND=rk),DIMENSION(13),INTENT(IN)::mat_param_inelast
       REAL(KIND=rk),DIMENSION(24),INTENT(IN)::plastic_slip_n
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::unimod_C,F_p_inv_n,C_inv
       REAL(KIND=rk),DIMENSION(6,3,3),INTENT(IN)::orth_basis
       REAL(KIND=rk),DIMENSION(6,6,6),INTENT(IN)::orth_proj
       INTEGER(KIND=ik),INTENT(IN)::algo,symm_class,NDI,NSHR,hard_law
       LOGICAL,INTENT(OUT)::conv_issue
       LOGICAL,DIMENSION(24),INTENT(OUT)::active_slip_system
       CHARACTER(4),INTENT(IN)::env
     END SUBROUTINE determine_incremental_lagrange_multiplier_miehe_viscoplastic
  END INTERFACE determine_incremental_lagrange_multiplier_miehe_viscoplastic
  !
  INTERFACE diagonal_funct_deriv_projectors
     SUBROUTINE diagonal_funct_deriv_projectors(e_vect,d_vect,f_vect,Proj_tens,&
          &lambda,index_distinct,n_distinct,Ce,m_SH,tol)
       USE Abaqus_Interface
       IMPLICIT NONE
       INTEGER(KIND=ik),INTENT(IN)::n_distinct
       INTEGER(KIND=ik),DIMENSION(3),INTENT(IN)::index_distinct
       REAL(KIND=rk),DIMENSION(3),INTENT(IN)::lambda
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::Ce
       REAL(KIND=rk),INTENT(IN)::tol,m_SH
       REAL(KIND=rk),DIMENSION(3),INTENT(OUT)::e_vect,d_vect,f_vect
       REAL(KIND=rk),DIMENSION(3,3,3),INTENT(OUT)::Proj_tens
     END SUBROUTINE diagonal_funct_deriv_projectors
  END INTERFACE diagonal_funct_deriv_projectors
  !
  INTERFACE eval_yield_function_linearization
     SUBROUTINE eval_yield_function_linearization(Phi_alpha,&
          &dPhi_alpha_d_delta_lambda_ip1_beta,tau_rss,&
          &d_tau_alpha_d_delta_lambda_ip1_beta,plastic_slip,i_slip,algo,&
          &mat_param_inelast,hard_law)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),INTENT(OUT)::Phi_alpha
       REAL(KIND=rk),DIMENSION(24),INTENT(OUT)::dPhi_alpha_d_delta_lambda_ip1_beta
       REAL(KIND=rk),INTENT(IN)::tau_rss
       REAL(KIND=rk),DIMENSION(13),INTENT(IN)::mat_param_inelast
       REAL(KIND=rk),DIMENSION(24),INTENT(IN)::d_tau_alpha_d_delta_lambda_ip1_beta
       REAL(KIND=rk),DIMENSION(24),INTENT(IN)::plastic_slip
       INTEGER(KIND=ik),INTENT(IN)::i_slip,algo,hard_law
     END SUBROUTINE eval_yield_function_linearization
  END INTERFACE eval_yield_function_linearization
  !
  INTERFACE hardening_function_cailletaud_forest
     SUBROUTINE hardening_function_cailletaud_forest(Y_alpha,&
          &dY_alpha_d_delta_lambda_beta,plastic_slip,alpha_slip,&
          &hard_param_cf)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),INTENT(OUT)::Y_alpha
       REAL(KIND=rk),DIMENSION(24),INTENT(OUT)::dY_alpha_d_delta_lambda_beta
       REAL(KIND=rk),DIMENSION(10),INTENT(IN)::hard_param_cf
       REAL(KIND=rk),DIMENSION(24),INTENT(IN)::plastic_slip
       INTEGER(KIND=ik),INTENT(IN)::alpha_slip
     END SUBROUTINE hardening_function_cailletaud_forest
  END INTERFACE hardening_function_cailletaud_forest
  !
  INTERFACE hardening_function_schmidt_baldassari
     SUBROUTINE hardening_function_schmidt_baldassari(Y_alpha,&
          &dY_alpha_d_delta_lambda_beta,plastic_slip,hard_param_sb)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),INTENT(OUT)::Y_alpha
       REAL(KIND=rk),DIMENSION(24),INTENT(OUT)::dY_alpha_d_delta_lambda_beta
       REAL(KIND=rk),DIMENSION(3),INTENT(IN)::hard_param_sb
       REAL(KIND=rk),DIMENSION(24),INTENT(IN)::plastic_slip
     END SUBROUTINE hardening_function_schmidt_baldassari
  END INTERFACE hardening_function_schmidt_baldassari
  !
  INTERFACE isotr_elasticity_mandel_stress_volumetric_isochoric_decoupled
     SUBROUTINE isotr_elasticity_mandel_stress_volumetric_isochoric_decoupled(&
          &M_hyd,M_dev,Je,unimod_Ce,kappa,mu)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),INTENT(OUT)::M_hyd
       REAL(KIND=rk),DIMENSION(3,3),INTENT(OUT)::M_dev
       REAL(KIND=rk),INTENT(IN)::Je,kappa,mu
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::unimod_Ce
     END SUBROUTINE isotr_elasticity_mandel_stress_volumetric_isochoric_decoupled
  END INTERFACE isotr_elasticity_mandel_stress_volumetric_isochoric_decoupled
  !
  INTERFACE linear_elasticity_subspace_projection
     SUBROUTINE linear_elasticity_subspace_projection(stress,tangent,strain,&
          &elast_const,symm_class,orth_basis,orth_proj)
       USE Abaqus_Interface
       IMPLICIT NONE
       INTEGER(KIND=ik),INTENT(IN)::symm_class
       REAL(KIND=rk),DIMENSION(3),INTENT(IN)::elast_const
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::strain
       REAL(KIND=rk),DIMENSION(6,3,3),INTENT(IN)::orth_basis
       REAL(KIND=rk),DIMENSION(6,6,6),INTENT(IN)::orth_proj
       REAL(KIND=rk),DIMENSION(3,3),INTENT(OUT)::stress
       REAL(KIND=rk),DIMENSION(6,6),INTENT(OUT)::tangent
     END SUBROUTINE linear_elasticity_subspace_projection
  END INTERFACE linear_elasticity_subspace_projection
  !
  INTERFACE newton_back_track_solver_crystal_visco_plasticity
     SUBROUTINE newton_back_track_solver_crystal_visco_plasticity(&
          &delta_lambda_ip1,conv_issue,unimod_C,F_p_inv_n,Je,&
          &plastic_slip_n,delta_lambda_i,algo,dt,elast_const,&
          &mat_param_inelast,m_SH,symm_class,orth_basis,orth_proj,NDI,&
          &NSHR,env,hard_law,theta_KK,eta,tol_nw,n_iter_max_nw,n_max_bt)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),INTENT(IN)::Je,dt,m_SH,theta_KK,eta,tol_nw
       REAL(KIND=rk),DIMENSION(3),INTENT(IN)::elast_const
       REAL(KIND=rk),DIMENSION(13),INTENT(IN)::mat_param_inelast
       REAL(KIND=rk),DIMENSION(24),INTENT(IN)::delta_lambda_i
       REAL(KIND=rk),DIMENSION(24),INTENT(IN)::plastic_slip_n
       REAL(KIND=rk),DIMENSION(24),INTENT(INOUT)::delta_lambda_ip1
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::unimod_C,F_p_inv_n
       REAL(KIND=rk),DIMENSION(6,3,3),INTENT(IN)::orth_basis
       REAL(KIND=rk),DIMENSION(6,6,6),INTENT(IN)::orth_proj
       INTEGER(KIND=ik),INTENT(IN)::algo,symm_class,NDI,NSHR,hard_law
       INTEGER(KIND=ik),INTENT(IN)::n_iter_max_nw,n_max_bt
       LOGICAL,INTENT(OUT)::conv_issue
       CHARACTER(4),INTENT(IN)::env
     END SUBROUTINE newton_back_track_solver_crystal_visco_plasticity
  END INTERFACE newton_back_track_solver_crystal_visco_plasticity
  !
  INTERFACE orthonormal_basis_subspace
     SUBROUTINE orthonormal_basis_subspace(orth_basis,orth_proj,symm_class,&
          &Rot_tens)
       USE Abaqus_Interface
       IMPLICIT NONE
       INTEGER(KIND=ik),INTENT(IN)::symm_class
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::Rot_tens
       REAL(KIND=rk),DIMENSION(6,3,3),INTENT(OUT)::orth_basis
       REAL(KIND=rk),DIMENSION(6,6,6),INTENT(OUT)::orth_proj
     END SUBROUTINE orthonormal_basis_subspace
  END INTERFACE orthonormal_basis_subspace
  !
  INTERFACE projection_tensors_stress_tangent_transformation
     SUBROUTINE projection_tensors_stress_tangent_transformation(P_tens,&
          &T_dc_L,Proj_tens,d_vect,f_vect,eta,xi,theta,n_distinct,T_stress,&
          &env)
       USE Abaqus_Interface
       IMPLICIT NONE
       INTEGER(KIND=ik),INTENT(IN)::n_distinct
       REAL(KIND=rk),DIMENSION(6,6),INTENT(OUT)::P_tens,T_dc_L
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::T_stress
       REAL(KIND=rk),DIMENSION(3,3,3),INTENT(IN)::Proj_tens
       REAL(KIND=rk),DIMENSION(3),INTENT(IN)::d_vect,f_vect
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::xi,theta
       REAL(KIND=rk),INTENT(IN)::eta
       CHARACTER(4),INTENT(IN)::env
     END SUBROUTINE projection_tensors_stress_tangent_transformation
  END INTERFACE projection_tensors_stress_tangent_transformation
  !
  INTERFACE residual_tangent_computation_nw_crystal_plasticity
     SUBROUTINE residual_tangent_computation_nw_crystal_plasticity(&
          &residual_nw,jacobian_nw,active_slip_system,&
          &deriv_visco_plast_potent,unimod_C,F_p_inv_iter,F_p_inv_n,Je,&
          &plastic_slip_iter,delta_lambda_i,delta_lambda_ip1,algo,dt,&
          &m_SH,elast_const,symm_class,orth_basis,orth_proj,&
          &mat_param_inelast,NDI,NSHR,env,eta,hard_law,theta_KK,Phi_alpha_all)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(24),INTENT(OUT)::residual_nw
       REAL(KIND=rk),DIMENSION(24),INTENT(OUT)::Phi_alpha_all
       REAL(KIND=rk),DIMENSION(24),INTENT(OUT)::deriv_visco_plast_potent
       REAL(KIND=rk),DIMENSION(24,24),INTENT(OUT)::jacobian_nw
       REAL(KIND=rk),INTENT(IN)::Je,dt,eta,m_SH,theta_KK
       REAL(KIND=rk),DIMENSION(3),INTENT(IN)::elast_const
       REAL(KIND=rk),DIMENSION(13),INTENT(IN)::mat_param_inelast
       REAL(KIND=rk),DIMENSION(24),INTENT(IN)::plastic_slip_iter
       REAL(KIND=rk),DIMENSION(24),INTENT(IN)::delta_lambda_i
       REAL(KIND=rk),DIMENSION(24),INTENT(IN)::delta_lambda_ip1
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::unimod_C,F_p_inv_iter
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::F_p_inv_n
       REAL(KIND=rk),DIMENSION(6,3,3),INTENT(IN)::orth_basis
       REAL(KIND=rk),DIMENSION(6,6,6),INTENT(IN)::orth_proj
       INTEGER(KIND=ik),INTENT(IN)::algo,symm_class,NDI,NSHR,hard_law
       LOGICAL,DIMENSION(24),INTENT(OUT)::active_slip_system
       CHARACTER(4),INTENT(IN)::env
     END SUBROUTINE residual_tangent_computation_nw_crystal_plasticity
  END INTERFACE residual_tangent_computation_nw_crystal_plasticity
  !
  INTERFACE residual_tangent_nw_crystal_plasticity_generic_solver
     SUBROUTINE residual_tangent_nw_crystal_plasticity_generic_solver(&
          &residual_nw,jacobian_nw,unimod_C,F_p_inv_n,Je,plastic_slip_n,&
          &delta_lambda_i,delta_lambda_ip1,algo,dt,m_SH,elast_const,&
          &symm_class,orth_basis,orth_proj,mat_param_inelast,NDI,NSHR,&
          &env,eta,hard_law,theta_KK)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(24),INTENT(OUT)::residual_nw
       REAL(KIND=rk),DIMENSION(24,24),INTENT(OUT)::jacobian_nw
       REAL(KIND=rk),INTENT(IN)::Je,dt,eta,m_SH,theta_KK
       REAL(KIND=rk),DIMENSION(3),INTENT(IN)::elast_const
       REAL(KIND=rk),DIMENSION(13),INTENT(IN)::mat_param_inelast
       REAL(KIND=rk),DIMENSION(24),INTENT(IN)::plastic_slip_n
       REAL(KIND=rk),DIMENSION(24),INTENT(IN)::delta_lambda_i
       REAL(KIND=rk),DIMENSION(24),INTENT(IN)::delta_lambda_ip1
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::unimod_C
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::F_p_inv_n
       REAL(KIND=rk),DIMENSION(6,3,3),INTENT(IN)::orth_basis
       REAL(KIND=rk),DIMENSION(6,6,6),INTENT(IN)::orth_proj
       INTEGER(KIND=ik),INTENT(IN)::algo,symm_class,NDI,NSHR,hard_law
       CHARACTER(4),INTENT(IN)::env
     END SUBROUTINE residual_tangent_nw_crystal_plasticity_generic_solver
  END INTERFACE residual_tangent_nw_crystal_plasticity_generic_solver
  !
  INTERFACE resolved_shear_stress_linearization
     SUBROUTINE resolved_shear_stress_linearization(tau_alpha_rss,&
          &d_tau_alpha_d_delta_lambda_ip1_beta,M_dev,N_alpha,&
          &d_Mandel_stress_d_delta_lambda_ip1_beta,env)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),INTENT(OUT)::tau_alpha_rss
       REAL(KIND=rk),DIMENSION(24),INTENT(OUT)::d_tau_alpha_d_delta_lambda_ip1_beta
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::M_dev,N_alpha
       REAL(KIND=rk),DIMENSION(9,24),INTENT(IN)::d_Mandel_stress_d_delta_lambda_ip1_beta
       CHARACTER(4),INTENT(IN)::env
     END SUBROUTINE resolved_shear_stress_linearization
  END INTERFACE resolved_shear_stress_linearization
  !
  INTERFACE slip_system_fcc
     SUBROUTINE slip_system_fcc(normal_dir,slip_dir,index)
       USE Abaqus_Interface
       IMPLICIT NONE
       INTEGER(KIND=ik),INTENT(IN)::index
       REAL(KIND=rk),DIMENSION(3),INTENT(OUT)::normal_dir,slip_dir
     END SUBROUTINE slip_system_fcc
  END INTERFACE slip_system_fcc
  !
  INTERFACE stress_computation_single_crystal_plasticity
     SUBROUTINE stress_computation_single_crystal_plasticity(S2PK,&
          &S2PK_IC,F_p_inv_n1,plastic_slip_iter,d_delta_lambda_ip1_d_C,&
          &active_slip_system,Ce_n1,Ce_n1_inv,Je_n1,C,det_Fp_inv,&
          &F_p_inv_iter,F,dt,F_p_inv_n,plastic_slip_n,elast_const,&
          &mat_param_inelast,m_SH,symm_class,orth_basis,orth_proj,NDI,&
          &NSHR,algo,env,hard_law,theta_KK,aver_adj_stp,conv_issue_dil)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),INTENT(OUT)::Je_n1,det_Fp_inv
       REAL(KIND=rk),DIMENSION(24),INTENT(OUT)::plastic_slip_iter
       REAL(KIND=rk),DIMENSION(3,3),INTENT(OUT)::S2PK,S2PK_IC,F_p_inv_n1
       REAL(KIND=rk),DIMENSION(3,3),INTENT(OUT)::Ce_n1,Ce_n1_inv,C
       REAL(KIND=rk),DIMENSION(3,3),INTENT(OUT)::F_p_inv_iter
       REAL(KIND=rk),DIMENSION(24,6),INTENT(OUT)::d_delta_lambda_ip1_d_C
       REAL(KIND=rk),INTENT(IN)::dt,m_SH,theta_KK
       REAL(KIND=rk),INTENT(INOUT)::aver_adj_stp
       REAL(KIND=rk),DIMENSION(3),INTENT(IN)::elast_const
       REAL(KIND=rk),DIMENSION(13),INTENT(IN)::mat_param_inelast
       REAL(KIND=rk),DIMENSION(24),INTENT(IN)::plastic_slip_n
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::F,F_p_inv_n
       REAL(KIND=rk),DIMENSION(6,3,3),INTENT(IN)::orth_basis
       REAL(KIND=rk),DIMENSION(6,6,6),INTENT(IN)::orth_proj
       INTEGER(KIND=ik),INTENT(IN)::algo,symm_class,NDI,NSHR,hard_law
       LOGICAL,INTENT(OUT)::conv_issue_dil
       LOGICAL,DIMENSION(24),INTENT(OUT)::active_slip_system
       CHARACTER(4),INTENT(IN)::env
     END SUBROUTINE stress_computation_single_crystal_plasticity
  END INTERFACE stress_computation_single_crystal_plasticity
  !
  INTERFACE S2PK_IC_cubic_elasticity_seth_hill
     SUBROUTINE S2PK_IC_cubic_elasticity_seth_hill(S_2PK,C_mat_tang,Ce,m_SH,&
          &elast_const,symm_class,orth_basis,orth_proj,NDI,NSHR,env)
       USE Abaqus_Interface
       IMPLICIT NONE
       INTEGER(KIND=ik),INTENT(IN)::symm_class,NDI,NSHR
       REAL(KIND=rk),INTENT(IN)::m_SH
       REAL(KIND=rk),DIMENSION(3),INTENT(IN)::elast_const
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::Ce
       REAL(KIND=rk),DIMENSION(6,3,3),INTENT(IN)::orth_basis
       REAL(KIND=rk),DIMENSION(6,6,6),INTENT(IN)::orth_proj
       REAL(KIND=rk),DIMENSION(6),INTENT(OUT)::S_2PK
       REAL(KIND=rk),DIMENSION(6,6),INTENT(OUT)::C_mat_tang
       CHARACTER(4),INTENT(IN)::env
     END SUBROUTINE S2PK_IC_cubic_elasticity_seth_hill
  END INTERFACE S2PK_IC_cubic_elasticity_seth_hill
  !
  INTERFACE update_internal_variable_nw_iter
     SUBROUTINE update_internal_variable_nw_iter(delta_lambda_ip1,&
          &plastic_slip_iter,F_p_inv_iter,dd_delta_lambda_ip1,&
          &plastic_slip_n,F_p_inv_n)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(24),INTENT(INOUT)::delta_lambda_ip1
       REAL(KIND=rk),DIMENSION(24),INTENT(OUT)::plastic_slip_iter
       REAL(KIND=rk),DIMENSION(3,3),INTENT(OUT)::F_p_inv_iter
       REAL(KIND=rk),DIMENSION(24),INTENT(IN)::dd_delta_lambda_ip1
       REAL(KIND=rk),DIMENSION(24),INTENT(IN)::plastic_slip_n
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::F_p_inv_n
     END SUBROUTINE update_internal_variable_nw_iter
  END INTERFACE update_internal_variable_nw_iter
! ----------------------------------------------------------------------
  INTERFACE aux_prod_projector
     FUNCTION aux_prod_projector(A_tens,lambda,index_distinct,p_app,&
          &q_app)
       USE Abaqus_Interface
       IMPLICIT NONE
       INTEGER(KIND=ik),INTENT(IN)::p_app,q_app
       INTEGER(KIND=ik),DIMENSION(3),INTENT(IN)::index_distinct
       REAL(KIND=rk),DIMENSION(3),INTENT(IN)::lambda
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::A_tens
       REAL(KIND=rk),DIMENSION(3,3)::aux_prod_projector
     END FUNCTION aux_prod_projector
  END INTERFACE aux_prod_projector
  !
  INTERFACE comp_add_stress_contrib_tang
     FUNCTION comp_add_stress_contrib_tang(STRESS)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(6),INTENT(IN)::STRESS
       REAL(KIND=rk),DIMENSION(6,6)::comp_add_stress_contrib_tang
     END FUNCTION comp_add_stress_contrib_tang
  END INTERFACE comp_add_stress_contrib_tang
  !
  INTERFACE comp_d_delta_lambda_ip1_d_C
     FUNCTION comp_d_delta_lambda_ip1_d_C(active_slip_system,&
          &jacobian_nw,d_delta_lambda_i_d_C,F_p_inv_iter,unimod_C,Je,&
          &C_inv,elast_const,m_SH,symm_class,orth_basis,orth_proj,eta,&
          &dt,algo,deriv_visco_plast_potent,env,NDI,NSHR)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),INTENT(IN)::Je,m_SH,eta,dt
       REAL(KIND=rk),DIMENSION(3),INTENT(IN)::elast_const
       REAL(KIND=rk),DIMENSION(24),INTENT(IN)::deriv_visco_plast_potent
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::F_p_inv_iter,unimod_C
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::C_inv
       REAL(KIND=rk),DIMENSION(24,6),INTENT(IN)::d_delta_lambda_i_d_C
       REAL(KIND=rk),DIMENSION(24,24),INTENT(IN)::jacobian_nw
       REAL(KIND=rk),DIMENSION(6,3,3),INTENT(IN)::orth_basis
       REAL(KIND=rk),DIMENSION(6,6,6),INTENT(IN)::orth_proj
       REAL(KIND=rk),DIMENSION(24,6)::comp_d_delta_lambda_ip1_d_C
       INTEGER(KIND=ik),INTENT(IN)::algo,symm_class,NDI,NSHR
       CHARACTER(4),INTENT(IN)::env
       LOGICAL,DIMENSION(24),INTENT(IN)::active_slip_system
     END FUNCTION comp_d_delta_lambda_ip1_d_C
  END INTERFACE comp_d_delta_lambda_ip1_d_C
  !
  INTERFACE comp_deriv_d_C_e_d_F_p_inv
     FUNCTION comp_deriv_d_C_e_d_F_p_inv(F_p_inv,C,env)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::F_p_inv,C
       REAL(KIND=rk),DIMENSION(3,3,3,3)::comp_deriv_d_C_e_d_F_p_inv
       CHARACTER(4),INTENT(IN)::env
     END FUNCTION comp_deriv_d_C_e_d_F_p_inv
  END INTERFACE comp_deriv_d_C_e_d_F_p_inv
  !
  INTERFACE comp_derivative_d_S2PK_IC_d_Ce_IC
     FUNCTION comp_derivative_d_S2PK_IC_d_Ce_IC (Je,Ce,Ce_inv,kappa,mu,&
          &env)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),INTENT(IN)::Je,kappa,mu
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::Ce,Ce_inv
       REAL(KIND=rk),DIMENSION(3,3,3,3)::comp_derivative_d_S2PK_IC_d_Ce_IC
       CHARACTER(4),INTENT(IN)::env
     END FUNCTION comp_derivative_d_S2PK_IC_d_Ce_IC
  END INTERFACE comp_derivative_d_S2PK_IC_d_Ce_IC
  !
  INTERFACE consistent_algorithmic_tangent_augmented_lagrangian_penalty
     FUNCTION consistent_algorithmic_tangent_augmented_lagrangian_penalty(&
          &C,F_p_inv_n,active_slip_system,d_delta_lambda_ip1_d_C,&
          &det_Fp_inv,F_p_inv_iter,Je_n1,Ce_n1,Ce_n1_inv,F_p_inv_n1,&
          &S2PK_IC,elast_const,m_SH,symm_class,orth_basis,orth_proj,NDI,&
          &NSHR,env)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),INTENT(IN)::det_Fp_inv,Je_n1,m_SH
       REAL(KIND=rk),DIMENSION(3),INTENT(IN)::elast_const
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::C,F_p_inv_n,F_p_inv_iter
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::Ce_n1,Ce_n1_inv
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::S2PK_IC,F_p_inv_n1
       REAL(KIND=rk),DIMENSION(24,6),INTENT(IN)::d_delta_lambda_ip1_d_C
       REAL(KIND=rk),DIMENSION(6,3,3),INTENT(IN)::orth_basis
       REAL(KIND=rk),DIMENSION(6,6,6),INTENT(IN)::orth_proj
       REAL(KIND=rk),DIMENSION(6,6)::consistent_algorithmic_tangent_augmented_lagrangian_penalty
       INTEGER(KIND=ik),INTENT(IN)::symm_class,NDI,NSHR
       LOGICAL,DIMENSION(24),INTENT(IN)::active_slip_system
       CHARACTER(4),INTENT(IN)::env
     END FUNCTION consistent_algorithmic_tangent_augmented_lagrangian_penalty
  END INTERFACE consistent_algorithmic_tangent_augmented_lagrangian_penalty
  !
  INTERFACE CSDA_material_tangent_single_crystal_plast
     FUNCTION CSDA_material_tangent_single_crystal_plast(S2PK,F,dt,&
          &F_p_inv_n,plastic_slip_n,elast_const,mat_param_inelast,m_SH,&
          &symm_class,orth_basis,orth_proj,numtang,algo,NDI,NSHR,env,&
          &hard_law,theta_KK,aver_adj_stp)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),INTENT(IN)::dt,m_SH,theta_KK,aver_adj_stp
       REAL(KIND=rk),DIMENSION(3),INTENT(IN)::elast_const
       REAL(KIND=rk),DIMENSION(13),INTENT(IN)::mat_param_inelast
       REAL(KIND=rk),DIMENSION(24),INTENT(IN)::plastic_slip_n
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::S2PK,F,F_p_inv_n
       REAL(KIND=rk),DIMENSION(6,3,3),INTENT(IN)::orth_basis
       REAL(KIND=rk),DIMENSION(6,6,6),INTENT(IN)::orth_proj
       REAL(KIND=rk),DIMENSION(6,6)::CSDA_material_tangent_single_crystal_plast
       INTEGER(KIND=ik),INTENT(IN)::symm_class,algo,numtang,NDI,NSHR
       INTEGER(KIND=ik),INTENT(IN)::hard_law
       CHARACTER(4),INTENT(IN)::env
     END FUNCTION CSDA_material_tangent_single_crystal_plast
  END INTERFACE CSDA_material_tangent_single_crystal_plast
  !
  INTERFACE deriv_d_S_d_Fp_inv
     FUNCTION deriv_d_S_d_Fp_inv(S2PK_IC,F_p_inv,env)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::S2PK_IC,F_p_inv
       REAL(KIND=rk),DIMENSION(3,3,3,3)::deriv_d_S_d_Fp_inv
       CHARACTER(4),INTENT(IN)::env
     END FUNCTION deriv_d_S_d_Fp_inv
  END INTERFACE deriv_d_S_d_Fp_inv
  !
  INTERFACE derivative_inverse_plastic_def_grad_right_cauchy_green_tensor_c
     FUNCTION derivative_inverse_plastic_def_grad_right_cauchy_green_tensor_c(&
          &F_p_inv_n,active_slip_system,d_delta_lambda_ip1_d_C,env)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::F_p_inv_n
       REAL(KIND=rk),DIMENSION(24,6),INTENT(IN)::d_delta_lambda_ip1_d_C
       REAL(KIND=rk),DIMENSION(3,3,3,3)::derivative_inverse_plastic_def_grad_right_cauchy_green_tensor_c
       LOGICAL,DIMENSION(24),INTENT(IN)::active_slip_system
       CHARACTER(4),INTENT(IN)::env
     END FUNCTION derivative_inverse_plastic_def_grad_right_cauchy_green_tensor_c
  END INTERFACE derivative_inverse_plastic_def_grad_right_cauchy_green_tensor_c
  !
  INTERFACE derivative_mandel_stress_ce_sh_elasticity
     FUNCTION derivative_mandel_stress_ce_sh_elasticity(S_2PK_IC,&
          &C_mat_elast,Ce,env)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::S_2PK_IC,Ce
       REAL(KIND=rk),DIMENSION(6,6),INTENT(IN)::C_mat_elast
       REAL(KIND=rk),DIMENSION(3,3,3,3)::derivative_mandel_stress_ce_sh_elasticity
       CHARACTER(4),INTENT(IN)::env
     END FUNCTION derivative_mandel_stress_ce_sh_elasticity
  END INTERFACE derivative_mandel_stress_ce_sh_elasticity
  !
  INTERFACE extract_hardening_parameters_crystal_plast
     FUNCTION extract_hardening_parameters_crystal_plast(&
          &mat_param_inelast,hard_law,dim_hard)
       USE Abaqus_Interface
       IMPLICIT NONE
       INTEGER(KIND=ik),INTENT(IN)::hard_law,dim_hard
       REAL(KIND=rk),DIMENSION(13),INTENT(IN)::mat_param_inelast
       REAL(KIND=rk),DIMENSION(dim_hard)::extract_hardening_parameters_crystal_plast
     END FUNCTION extract_hardening_parameters_crystal_plast
  END INTERFACE extract_hardening_parameters_crystal_plast
  !
  INTERFACE extract_viscoplastic_parameters_crystal_plast
     FUNCTION extract_viscoplastic_parameters_crystal_plast(&
          &mat_param_inelast)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(13),INTENT(IN)::mat_param_inelast
       REAL(KIND=rk),DIMENSION(3)::extract_viscoplastic_parameters_crystal_plast
     END FUNCTION extract_viscoplastic_parameters_crystal_plast
  END INTERFACE extract_viscoplastic_parameters_crystal_plast
  !
  INTERFACE interaction_matrix
     FUNCTION interaction_matrix(hard_param_im)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(6),INTENT(IN)::hard_param_im
       REAL(KIND=rk),DIMENSION(24,24)::interaction_matrix
     END FUNCTION interaction_matrix
  END INTERFACE interaction_matrix
  !
  INTERFACE linearization_mandel_stress
     FUNCTION linearization_mandel_stress(F_p_inv_n,&
          &unimod_C_F_p_inv_iter,mu,env)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),INTENT(IN)::mu
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::F_p_inv_n
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::unimod_C_F_p_inv_iter
       REAL(KIND=rk),DIMENSION(9,24)::linearization_mandel_stress
       CHARACTER(4),INTENT(IN)::env
     END FUNCTION linearization_mandel_stress
  END INTERFACE linearization_mandel_stress
  !
  INTERFACE linearization_mandel_stress_elasticity_SH_strain
     FUNCTION linearization_mandel_stress_elasticity_SH_strain(&
          &F_p_inv_n,unimod_C_F_p_inv_iter,Ce,S_2PK_IC,C_mat_elast,env)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::F_p_inv_n,Ce,S_2PK_IC
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::unimod_C_F_p_inv_iter
       REAL(KIND=rk),DIMENSION(6,6),INTENT(IN)::C_mat_elast
       REAL(KIND=rk),DIMENSION(9,24)::linearization_mandel_stress_elasticity_SH_strain
       CHARACTER(4),INTENT(IN)::env
     END FUNCTION linearization_mandel_stress_elasticity_SH_strain
  END INTERFACE linearization_mandel_stress_elasticity_SH_strain
  !
  INTERFACE mod_d_F_p_inv_d_C_isochoric_proj_constraint
     FUNCTION mod_d_F_p_inv_d_C_isochoric_proj_constraint(&
          &d_F_p_inv_i1_d_C,det_Fp_inv,F_p_inv_iter)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),INTENT(IN)::det_Fp_inv
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::F_p_inv_iter
       REAL(KIND=rk),DIMENSION(3,3,3,3),INTENT(IN)::d_F_p_inv_i1_d_C
       REAL(KIND=rk),DIMENSION(3,3,3,3)::mod_d_F_p_inv_d_C_isochoric_proj_constraint
     END FUNCTION mod_d_F_p_inv_d_C_isochoric_proj_constraint
  END INTERFACE mod_d_F_p_inv_d_C_isochoric_proj_constraint
  !
  INTERFACE S2PK_IC_isotr_elasticity
     FUNCTION S2PK_IC_isotr_elasticity(Ce,Ce_inv,Je,kappa,mu)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),INTENT(IN)::Je,kappa,mu
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::Ce,Ce_inv
       REAL(KIND=rk),DIMENSION(3,3)::S2PK_IC_isotr_elasticity
     END FUNCTION S2PK_IC_isotr_elasticity
  END INTERFACE S2PK_IC_isotr_elasticity
  !
  INTERFACE update_plastic_deformation_gradient
     FUNCTION update_plastic_deformation_gradient(F_p_inv_n,&
          &delta_lambda_ip1)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(24),INTENT(IN)::delta_lambda_ip1
       REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::F_p_inv_n
       REAL(KIND=rk),DIMENSION(3,3)::update_plastic_deformation_gradient
     END FUNCTION update_plastic_deformation_gradient
  END INTERFACE update_plastic_deformation_gradient
  !
  INTERFACE update_visco_in_inelast_param_crystal_plast
     FUNCTION update_visco_in_inelast_param_crystal_plast(&
          &mat_param_visco_t,mat_param_inelast)
       USE Abaqus_Interface
       IMPLICIT NONE
       REAL(KIND=rk),DIMENSION(3),INTENT(IN)::mat_param_visco_t
       REAL(KIND=rk),DIMENSION(13),INTENT(IN)::mat_param_inelast
       REAL(KIND=rk),DIMENSION(13)::update_visco_in_inelast_param_crystal_plast
     END FUNCTION update_visco_in_inelast_param_crystal_plast
  END INTERFACE update_visco_in_inelast_param_crystal_plast
END MODULE constitutive_routines_interface
!
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!                        Auxiliary_routines
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! **********************************************************************
! Rearrange fourth order tensor into 9x9 matrix format or 6x6 matrix, if 
! tensor possesses minor symmetries
! **********************************************************************
SUBROUTINE arrange_condense_4_order_tensor(A_4ot,A_4ot_9b9,A_4ot_6b6,&
     &output_type)
  !
  USE Abaqus_Interface
  USE Constants
  !
  ! declaration of variables
  !
  ! A_4ot.........fourth order tensor
  ! A_4ot_9b9.....representation of tensor coefficients in 9x9 matrix
  !               format
  ! A_4ot_6b6.....representation of tensor coefficients in 6x6 matrix
  !               format due to minor symmetries in fourth order tensor
  ! output_type...identifier for output quantity of interest
  !
  IMPLICIT NONE
  !
  INTEGER(KIND=ik)::m_ac4ot,n_ac4ot,i_ac4ot,j_ac4ot,k_ac4ot,l_ac4ot
  !
  REAL(KIND=rk),DIMENSION(3,3,3,3),INTENT(IN)::A_4ot
  REAL(KIND=rk),DIMENSION(9,9),INTENT(OUT)::A_4ot_9b9
  REAL(KIND=rk),DIMENSION(6,6),INTENT(OUT)::A_4ot_6b6
  !
  CHARACTER(3),INTENT(IN)::output_type
  !
  SELECT CASE (output_type)
  CASE ('9b9')
     DO m_ac4ot=ONE_int,NINE_int
        DO n_ac4ot=ONE_int,NINE_int
           i_ac4ot=INDEX_NINE_ABAQ(m_ac4ot,ONE_int)
           j_ac4ot=INDEX_NINE_ABAQ(m_ac4ot,TWO_int)
           k_ac4ot=INDEX_NINE_ABAQ(n_ac4ot,ONE_int)
           l_ac4ot=INDEX_NINE_ABAQ(n_ac4ot,TWO_int)
           A_4ot_9b9(m_ac4ot,n_ac4ot)=A_4ot(i_ac4ot,j_ac4ot,k_ac4ot,&
                &l_ac4ot)
        END DO
     END DO
     A_4ot_6b6=ZERO_r
  CASE ('6b6')
     DO m_ac4ot=ONE_int,SIX_int
        DO n_ac4ot=ONE_int,SIX_int
           i_ac4ot=INDEX_SIX_ABAQ(m_ac4ot,ONE_int)
           j_ac4ot=INDEX_SIX_ABAQ(m_ac4ot,TWO_int)
           k_ac4ot=INDEX_SIX_ABAQ(n_ac4ot,ONE_int)
           l_ac4ot=INDEX_SIX_ABAQ(n_ac4ot,TWO_int)
           A_4ot_6b6(m_ac4ot,n_ac4ot)=A_4ot(i_ac4ot,j_ac4ot,k_ac4ot,&
                &l_ac4ot)
        END DO
     END DO
     A_4ot_9b9=ZERO_r
  CASE DEFAULT
     WRITE(16,'(A)')'Unknown identifier in arrange_condense_4_order_tensor'
     WRITE(16,'(A)')'STOPPING COMPUTATION'
     STOP
  END SELECT
  !
END SUBROUTINE arrange_condense_4_order_tensor
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! contains logic for line search with back-tracking based on quadratic/
! cubic interpolation
! adopted from: Press, et al., "Numerical recipes in Fortran",Vol.1, 
! 1992, Cambridge University Press
! **********************************************************************
SUBROUTINE back_tracking_step_length_scaling_parameter(back_tracking,&
     &back_track_count,lambda,lambda_two,x,x_save,dx,obj_f,obj_f_save,&
     &obj_f_two,grad_f_save,ndim)
  !
  USE Abaqus_Interface
  USE Constants
  !
  ! variable declaration
  ! back_tracking..........switch variable to activate back-tracking
  ! back_track_count.......number of back-track iterations
  ! lambda (_two)..........scaling parameter for iterative improvement
  ! x (_save)..............independent variables
  ! dx.....................iterative improvement (Newton-step)
  ! obj_f (_save,two)......value of objective function to be reduced
  ! grad_f_save............gradient of objective function
  !
  IMPLICIT NONE
  INTEGER(KIND=ik),INTENT(IN)::ndim
  INTEGER(KIND=ik),INTENT(INOUT)::back_track_count
  REAL(KIND=rk),INTENT(INOUT)::lambda,lambda_two,obj_f_two
  REAL(KIND=rk),INTENT(IN)::obj_f,obj_f_save
  REAL(KIND=rk),DIMENSION(ndim),INTENT(INOUT)::x
  REAL(KIND=rk),DIMENSION(ndim),INTENT(IN)::x_save,dx,grad_f_save
  REAL(KIND=rk)::grad_f_dp_dx,obj_f_one,lambda_one
  REAL(KIND=rk),DIMENSION(2)::rhs,param_ab
  REAL(KIND=rk),DIMENSION(2,2)::Mat
  LOGICAL,INTENT(INOUT)::back_tracking
  !
  grad_f_dp_dx=DOT_PRODUCT(grad_f_save,dx)
  !
  IF (obj_f<obj_f_save+lambda*alpha_bt*grad_f_dp_dx) THEN
     back_tracking=.FALSE.
     back_track_count=ZERO_int
  ELSE
     IF (.NOT.(back_tracking)) THEN
        ! quadratic interpolation
        back_tracking=.TRUE.
        lambda=MAX(-(grad_f_dp_dx)/TWO_r/(obj_f-obj_f_save-grad_f_dp_dx),&
             &ONE_TENTH_r)
        lambda=MIN(lambda,ONE_HALF_r)
        x=x_save+lambda*dx
        lambda_two=ONE_r
        obj_f_two=obj_f
        back_track_count=back_track_count+ONE_int
        !
     ELSE
        ! cubic interpolation
        obj_f_one=obj_f 
        lambda_one=lambda
        !
        ! solve for the parameters of the cubic interpolation
        Mat=RESHAPE((/ONE_r/lambda_one**TWO_int,&
             &-lambda_two/lambda_one**TWO_int,-ONE_r/lambda_two**TWO_int,&
             &lambda_one/lambda_two**TWO_int/),(/2,2/))
        rhs=(/obj_f_one-grad_f_dp_dx*lambda_one-obj_f_save,obj_f_two-&
             &grad_f_dp_dx*lambda_two-obj_f_save/)
        param_ab=MATMUL(Mat,rhs); param_ab=param_ab/(lambda_one-lambda_two)
        !
        ! compute lambda and ensure bounds on step length parameter
        ! 0.1*lambda_one <= lambda <= 0.5*lambda_one
        lambda=(-param_ab(2)+SQRT(param_ab(2)**TWO_int-THREE_r*param_ab(1)*&
             &grad_f_dp_dx))/(THREE_r*param_ab(1))
        lambda=MIN(lambda,ONE_HALF_r*lambda_one)
        lambda=MAX(lambda,ONE_TENTH_r*lambda_one)
        !
        x=x_save+lambda*dx;
        !
        ! shift values of objective function and step length parameter lambda
        obj_f_two=obj_f_one
        lambda_two=lambda_one
        back_track_count=back_track_count+ONE_int
        !
     END IF
  END IF
  !
END SUBROUTINE back_tracking_step_length_scaling_parameter
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! ......................................................................
!  M33INV  -  Compute the inverse of a 3x3 matrix.
!
!  Programmer:   David G. Simpson
!                NASA Goddard Space Flight Center
!                Greenbelt, Maryland  20771
!
! **********************************************************************
SUBROUTINE M33INV (A, AINV, OK_FLAG)
  !
  USE Abaqus_Interface
  USE Constants
  !
  !
  !  A........input 3x3 matrix to be inverted
  !  AINV.....output 3x3 inverse of matrix A
  !  OK_FLAG..(output) .TRUE. if the input matrix could be inverted, and
  !           .FALSE. if the input matrix is singular.
  IMPLICIT NONE
  REAL(KIND=rk),DIMENSION(3,3), INTENT(IN)  :: A
  REAL(KIND=rk),DIMENSION(3,3), INTENT(OUT) :: AINV
  LOGICAL, INTENT(OUT) :: OK_FLAG
  REAL(KIND=rk):: DET
  REAL(KIND=rk),DIMENSION(3,3) :: COFACTOR
  !
  DET =  A(1,1)*A(2,2)*A(3,3)  &
       & - A(1,1)*A(2,3)*A(3,2) &
       & - A(1,2)*A(2,1)*A(3,3) &
       & + A(1,2)*A(2,3)*A(3,1) &
       & + A(1,3)*A(2,1)*A(3,2) &
       & - A(1,3)*A(2,2)*A(3,1)

  IF (ABS(DET) .LE. zero_distinct) THEN 
     AINV = ZERO_r
     OK_FLAG = .FALSE.
     RETURN
  END IF
  !
  COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
  COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
  COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
  COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
  COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
  COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
  COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
  COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
  COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))
  !
  AINV = TRANSPOSE(COFACTOR) / DET
  !
  OK_FLAG = .TRUE.
  !
END SUBROUTINE M33INV
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! perform push forward from 2PK to Kirchhoff stress and associated
! moduli (adopted from FEAP js_mathtools),(Abaqus/FEAP storage order)
! **********************************************************************
SUBROUTINE pushrc(S,AA,S0,AA0,DFGRD,iaa,env)
  !
  USE Abaqus_Interface
  USE Constants
  !
  ! NOTE: storage order in Abaqus: 11,22,33,12,13,23
  !       storage order in FEAP:   11,22,33,12,23,13
  !  S^ab    <= So^AB F^a_A F^b_B
  !  AA^abcd <= AAo^ABCD F^a_A F^b_B F^c_C F^d_D
  !  iaa = 0: transform stresses only
  !  iaa = 1: transform moduli only
  !  iaa = 2: transform stresses and moduli
  !
  ! declaration of variables
  ! S............Kirchhoff stress tensor (tau) (6x1 vector)
  ! AA...........spatial tangent (c) (6x6 matrix)
  ! S0...........2nd Piola-Kirchhoff stress (S) (6x1 vector)
  ! AA0..........material tangent (C) (6x6 matrix)
  ! DFGRD........deformation gradient (3x3 matrix)
  ! iaa..........task switch (see above)
  ! env..........storage order based on the environment in which the
  !              routine is called ('ABAQ','FEAP')
  IMPLICIT NONE
  !
  INTEGER(KIND=ik),INTENT(IN)::iaa
  INTEGER(KIND=ik)::n_prc,i_prc,j_prc,l1,l2,i1,i3,i4,m_prc,k_prc
  INTEGER(KIND=ik)::l_prc,l3,i5,l4,i2,i6
  INTEGER(KIND=ik),DIMENSION(3,3)::kk,kks
  INTEGER(KIND=ik),DIMENSION(6,2)::ind_lc
  REAL(KIND=rk),DIMENSION(6),INTENT(IN)::S0
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::DFGRD
  REAL(KIND=rk),DIMENSION(6,6),INTENT(IN)::AA0
  REAL(KIND=rk),DIMENSION(9)::DFGRD_vec
  REAL(KIND=rk),DIMENSION(6),INTENT(OUT)::S
  REAL(KIND=rk),DIMENSION(6,6),INTENT(OUT)::AA
  !
  CHARACTER(4),INTENT(IN)::env
  !
  ! set variables based on calling environment
  SELECT CASE (env)
  CASE ('FEAP')
     ind_lc=INDEX_SIX_FEAP
     kk=INDEX_TENS_VEC_FEAP; kks=INDEX_TENS_VEC_FEAP_S
     DO n_prc=ONE_int,NINE_int
        DFGRD_vec(n_prc)=&
             &DFGRD(INDEX_NINE_FEAP(n_prc,1),INDEX_NINE_FEAP(n_prc,2))
     END DO
  CASE ('ABAQ')
     ind_lc=INDEX_SIX_ABAQ
     kk=INDEX_TENS_VEC_ABAQ; kks=INDEX_TENS_VEC_ABAQ_S
     DO n_prc=ONE_int,NINE_int
        DFGRD_vec(n_prc)=&
             &DFGRD(INDEX_NINE_ABAQ(n_prc,1),INDEX_NINE_ABAQ(n_prc,2))
     END DO
  CASE DEFAULT
     WRITE(16,'(A)')'Unknown env var. in pushrc'
     WRITE(16,'(A)')'STOPPING COMPUTATION'
     STOP
  END SELECT
  !
  ! initialize the stress and the tangent
  S=ZERO_r; AA=ZERO_r
  ! transformation based on task id
  SELECT CASE (iaa)
  CASE (0) ! transform stress only
     DO n_prc=ONE_int,SIX_int
        i_prc=ind_lc(n_prc,ONE_int)
        j_prc=ind_lc(n_prc,TWO_int)
        DO l1=ONE_int,THREE_int
           i3=kk(i_prc,l1)
           DO l2=ONE_int,THREE_int
              i1=kks(l1,l2)
              i4=kk(j_prc,l2)
              S(n_prc)=S(n_prc)+S0(i1)*DFGRD_vec(i3)*DFGRD_vec(i4)
           END DO ! l2
        END DO ! l1
     END DO ! n_prc
     !
  CASE (1) ! transform moduli only
     DO n_prc=ONE_int,SIX_int
        i_prc=ind_lc(n_prc,ONE_int)
        j_prc=ind_lc(n_prc,TWO_int)
        DO l1=ONE_int,THREE_int
           i3=kk(i_prc,l1)
           DO l2=ONE_int,THREE_int
              i1=kks(l1,l2)
              i4=kk(j_prc,l2)
              DO m_prc=ONE_int,SIX_int
                 k_prc=ind_lc(m_prc,ONE_int)
                 l_prc=ind_lc(m_prc,TWO_int)
                 DO l3=ONE_int,THREE_int
                    i5=kk(k_prc,l3)
                    DO l4=ONE_int,THREE_int
                       i2=kks(l3,l4)
                       i6=kk(l_prc,l4)
                       AA(n_prc,m_prc) = AA(n_prc,m_prc) + AA0(i1,i2)*&
                            &DFGRD_vec(i3)*DFGRD_vec(i4)*DFGRD_vec(i5)*&
                            &DFGRD_vec(i6)
                    END DO ! l4
                 END DO ! l3
              END DO ! m_prc
           END DO ! l2
        END DO ! l1
     END DO ! n_prc
     !
  CASE (2) ! transform stress and moduli
     DO n_prc=ONE_int,SIX_int
        i_prc=ind_lc(n_prc,ONE_int)
        j_prc=ind_lc(n_prc,TWO_int)
        DO l1=ONE_int,THREE_int
           i3=kk(i_prc,l1)
           DO l2=ONE_int,THREE_int
              i1=kks(l1,l2)
              i4=kk(j_prc,l2)
              S(n_prc)=S(n_prc)+S0(i1)*DFGRD_vec(i3)*DFGRD_vec(i4)
              DO m_prc=ONE_int,SIX_int
                 k_prc=ind_lc(m_prc,ONE_int)
                 l_prc=ind_lc(m_prc,TWO_int)
                 DO l3=ONE_int,THREE_int
                    i5=kk(k_prc,l3)
                    DO l4=ONE_int,THREE_int
                       i2=kks(l3,l4)
                       i6=kk(l_prc,l4)
                       AA(n_prc,m_prc) = AA(n_prc,m_prc) + AA0(i1,i2)*&
                            &DFGRD_vec(i3)*DFGRD_vec(i4)*DFGRD_vec(i5)*&
                            &DFGRD_vec(i6)
                    END DO ! l4
                 END DO ! l3
              END DO ! m_prc
           END DO ! l2
        END DO ! l1
     END DO ! n_prc
     !
  CASE DEFAULT
     WRITE(16,'(A)')'Unknown iaa var. in pushrc'
     WRITE(16,'(A)')'STOPPING COMPUTATION'
     STOP
  END SELECT
END SUBROUTINE pushrc
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! transform symmetric second order tensor to vector format and vice
! versa (Abaqus and FEAP storage order)
! **********************************************************************
SUBROUTINE transform_symm_tens_vec_form(A_tens,A_tens_vec,transform_id,&
     &env)
  !
  USE Abaqus_Interface
  !
  ! NOTE: storage order in Abaqus: 11,22,33,12,13,23
  !       storage order in FEAP:   11,22,33,12,23,13
  ! declaration of variables
  !
  ! A_tens.........symmetric second order tensor
  ! A_tens_vec.....vector format of symmetric second order tensor
  ! transform_id...0=tensor to vector, 1=vector to tensor
  ! env............storage order based on the environment in which the
  !                routine is called ('ABAQ','FEAP')
  !
  IMPLICIT NONE
  !
  INTEGER(KIND=ik),INTENT(IN)::transform_id
  !
  REAL(KIND=rk),DIMENSION(3,3)::A_tens
  REAL(KIND=rk),DIMENSION(6)::A_tens_vec
  !
  CHARACTER(4),INTENT(IN)::env
  !
  SELECT CASE (transform_id)
  CASE (0) ! tensor to vector transformation
     A_tens_vec(1)=A_tens(1,1)
     A_tens_vec(2)=A_tens(2,2)
     A_tens_vec(3)=A_tens(3,3)
     A_tens_vec(4)=A_tens(1,2)
     SELECT CASE (env)
     CASE ('FEAP')
        A_tens_vec(5)=A_tens(2,3)
        A_tens_vec(6)=A_tens(1,3)
     CASE ('ABAQ')
        A_tens_vec(5)=A_tens(1,3)
        A_tens_vec(6)=A_tens(2,3)
     CASE DEFAULT
        WRITE(16,'(A)')'Unknown env var. in transform_symm_tens_vec_form'
        WRITE(16,'(A)')'STOPPING COMPUTATION'
        STOP
     END SELECT
     !
  CASE (1) ! vector to tensor transformation
     A_tens(1,1)=A_tens_vec(1)
     A_tens(2,2)=A_tens_vec(2)
     A_tens(3,3)=A_tens_vec(3)
     A_tens(1,2)=A_tens_vec(4);       A_tens(2,1)=A_tens(1,2)
     SELECT CASE (env)
     CASE ('FEAP')
        A_tens(2,3)=A_tens_vec(5);       A_tens(3,2)=A_tens(2,3)
        A_tens(1,3)=A_tens_vec(6);       A_tens(3,1)=A_tens(1,3)
     CASE ('ABAQ')
        A_tens(1,3)=A_tens_vec(5);       A_tens(3,1)=A_tens(1,3)
        A_tens(2,3)=A_tens_vec(6);       A_tens(3,2)=A_tens(2,3)
     CASE DEFAULT
        WRITE(16,'(A)')'Unknown env var. in transform_symm_tens_vec_form'
        WRITE(16,'(A)')'STOPPING COMPUTATION'
        STOP
     END SELECT
     !                    
  CASE DEFAULT
     WRITE(16,'(A)')'Unknown transformation id in transform_symm_tens_vec_form'
     WRITE(16,'(A)')'STOPPING COMPUTATION'
     STOP
  END SELECT
  !
END SUBROUTINE transform_symm_tens_vec_form
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! transform symmetric second order tensor to vector format and vice
! versa (Abaqus and FEAP storage order)
! **********************************************************************
SUBROUTINE transform_tens_vec_form(A_tens,A_tens_vec,transform_id,&
     &env)
  !
  USE Abaqus_Interface
  !
  ! NOTE: storage order in Abaqus: 11,22,33,12,21,13,31,23,32
  !       storage order in FEAP:   11,22,33,12,21,23,32,13,31
  ! declaration of variables
  !
  ! A_tens.........second order tensor
  ! A_tens_vec.....vector format of second order tensor
  ! transform_id...0=tensor to vector, 1=vector to tensor
  ! env............storage order based on the environment in which the
  !                routine is called ('ABAQ','FEAP')
  !
  IMPLICIT NONE
  !
  INTEGER(KIND=ik),INTENT(IN)::transform_id
  !
  REAL(KIND=rk),DIMENSION(3,3)::A_tens
  REAL(KIND=rk),DIMENSION(9)::A_tens_vec
  !
  CHARACTER(4),INTENT(IN)::env
  !
  SELECT CASE (transform_id)
  CASE (0) ! tensor to vector transformation
     A_tens_vec(1)=A_tens(1,1)
     A_tens_vec(2)=A_tens(2,2)
     A_tens_vec(3)=A_tens(3,3)
     A_tens_vec(4)=A_tens(1,2)
     A_tens_vec(5)=A_tens(2,1)
     SELECT CASE (env)
     CASE ('FEAP')
        A_tens_vec(6)=A_tens(2,3);
        A_tens_vec(7)=A_tens(3,2);
        A_tens_vec(8)=A_tens(1,3);
        A_tens_vec(9)=A_tens(3,1);
     CASE ('ABAQ')
        A_tens_vec(6)=A_tens(1,3);
        A_tens_vec(7)=A_tens(3,1);
        A_tens_vec(8)=A_tens(2,3);
        A_tens_vec(9)=A_tens(3,2);
     CASE DEFAULT
        WRITE(16,'(A)')'Unknown env var. in transform_tens_vec_form'
        WRITE(16,'(A)')'STOPPING COMPUTATION'
        STOP
     END SELECT
     !
  CASE (1) ! vector to tensor transformation
     A_tens(1,1)=A_tens_vec(1)
     A_tens(2,2)=A_tens_vec(2)
     A_tens(3,3)=A_tens_vec(3)
     A_tens(1,2)=A_tens_vec(4)
     A_tens(2,1)=A_tens_vec(5)
     SELECT CASE (env)
     CASE ('FEAP')
        A_tens(2,3)=A_tens_vec(6)
        A_tens(3,2)=A_tens_vec(7)
        A_tens(1,3)=A_tens_vec(8)
        A_tens(3,1)=A_tens_vec(9)
     CASE ('ABAQ')
        A_tens(1,3)=A_tens_vec(6)
        A_tens(3,1)=A_tens_vec(7)
        A_tens(2,3)=A_tens_vec(8)
        A_tens(3,2)=A_tens_vec(9)
     CASE DEFAULT
        WRITE(16,'(A)')'Unknown env var. in transform_tens_vec_form'
        WRITE(16,'(A)')'STOPPING COMPUTATION'
        STOP
     END SELECT
     !                    
  CASE DEFAULT
     WRITE(16,'(A)')'Unknown transformation id in transform_tens_vec_form'
     WRITE(16,'(A)')'STOPPING COMPUTATION'
     STOP
  END SELECT
  !
END SUBROUTINE transform_tens_vec_form
! ......................................................................
! **********************************************************************
!
! ######################################################################
!
! **********************************************************************
! assign the appropriate matrix of indices to convert tensor to vector
! **********************************************************************
FUNCTION approp_index(ndim,env,funct_name)
  !
  USE Abaqus_Interface
  USE Constants
  !
  IMPLICIT NONE
  INTEGER(KIND=ik),INTENT(IN)::ndim
  INTEGER(KIND=ik),DIMENSION(ndim,2)::approp_index,index_t
  CHARACTER(4),INTENT(IN)::env
  CHARACTER(*),INTENT(IN)::funct_name
  !
  SELECT CASE (ndim)
  CASE (SIX_int)
     SELECT CASE (env)
     CASE ('FEAP')
        index_t=INDEX_SIX_FEAP
     CASE ('ABAQ')
        index_t=INDEX_SIX_ABAQ
     CASE DEFAULT
        WRITE(16,'(2A)')'Unknown env var. in: ',funct_name
        WRITE(16,'(A)')'STOPPING COMPUTATION'
        STOP
     END SELECT
     !
  CASE (NINE_int)
     SELECT CASE (env)
     CASE ('FEAP')
        index_t=INDEX_NINE_FEAP
     CASE ('ABAQ')
        index_t=INDEX_NINE_ABAQ
     CASE DEFAULT
        WRITE(16,'(2A)')'Unknown env var. in: ',funct_name
        WRITE(16,'(A)')'STOPPING COMPUTATION'
        STOP
     END SELECT
     !
  CASE DEFAULT
  END SELECT
  approp_index=index_t
  !
END FUNCTION approp_index
! ...................................................................
! *******************************************************************
!
! *******************************************************************
! computation of determinant of a 3x3 matrix (Sarrus-rule)
! *******************************************************************
FUNCTION comp_det_Matrix(Matrix)
  !
  USE Abaqus_Interface
  !
  ! declaration of variables
  !
  ! comp_det_Matrix....determinant of the matrix
  ! Matrix.............matrix
  !
  IMPLICIT NONE
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::Matrix
  REAL(KIND=rk)::comp_det_Matrix
  !
  comp_det_Matrix = Matrix(1,1)*Matrix(2,2)*Matrix(3,3) + Matrix(2,1)*&
       & Matrix(3,2)*Matrix(1,3) + Matrix(3,1)*Matrix(1,2)*Matrix(2,3)&
       & - Matrix(1,3)*Matrix(2,2)*Matrix(3,1) - Matrix(2,3)*Matrix(3,2)*&
       & Matrix(1,1) - Matrix(3,3)*Matrix(1,2)*Matrix(2,1)

END FUNCTION comp_det_Matrix
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! compute second order tensor of rank one
! **********************************************************************
FUNCTION comp_sec_ord_tens_rank_one(a_vec,b_vec)
  !
  USE Abaqus_Interface
  USE Constants
  !
  IMPLICIT NONE
  REAL(KIND=rk),DIMENSION(3),INTENT(IN)::a_vec,b_vec
  REAL(KIND=rk),DIMENSION(3,3)::comp_sec_ord_tens_rank_one,A_temp
  INTEGER(KIND=ik)::m_csot,n_csot
  !
  A_temp=ZERO_r
  DO m_csot=ONE_int,THREE_int
     DO n_csot=ONE_int,THREE_int
        A_temp(m_csot,n_csot)=a_vec(m_csot)*b_vec(n_csot)
     END DO
  END DO
  comp_sec_ord_tens_rank_one=A_temp
  !
END FUNCTION comp_sec_ord_tens_rank_one
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! compute trace of a second order tensor
! **********************************************************************
FUNCTION compute_trace_second_order_tensor(A_tens)
  !
  USE Abaqus_Interface
  USE Constants
  !
  IMPLICIT NONE
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::A_tens
  REAL(KIND=rk)::tr_A,compute_trace_second_order_tensor
  INTEGER(KIND=ik)::i_cts
  !
  tr_A=ZERO_r
  DO i_cts=ONE_int,THREE_int
     tr_A=tr_A+A_tens(i_cts,i_cts)
  END DO
  compute_trace_second_order_tensor=tr_A
  !
END FUNCTION compute_trace_second_order_tensor
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! computes the double contraction of a fourth order tensor (given in 6x6
! format) and a symmetric second order tensor (6x1 format) by matrix
! multiplication
! **********************************************************************
FUNCTION double_contraction_fourth_order_second_order_tensor_matmul &
     &(A_fot_6x6,b_sot_6x1)
  !
  USE Abaqus_Interface
  USE Constants
  !
  ! declaration of variables
  ! A_fot_6x6.....symmetric fourth order tensor (stored in 6x6 matrix)
  ! b_sot_6x1.....symmetric second order tensor (stored in 6x1 vector)
  !
  IMPLICIT NONE
  !
  INTEGER(KIND=ik)::i_count
  REAL(KIND=rk),DIMENSION(6,6),INTENT(IN)::A_fot_6x6
  REAL(KIND=rk),DIMENSION(6),INTENT(IN)::b_sot_6x1
  REAL(KIND=rk),DIMENSION(6)::b_temp
  REAL(KIND=rk),DIMENSION(6)::double_contraction_fourth_order_second_order_tensor_matmul
  INTRINSIC MATMUL
  !
  b_temp=b_sot_6x1
  DO i_count=FOUR_int,SIX_int
     b_temp(i_count)=TWO_r*b_temp(i_count)
  END DO
  !
  double_contraction_fourth_order_second_order_tensor_matmul=MATMUL(&
       & A_fot_6x6,b_temp)
  !
END FUNCTION double_contraction_fourth_order_second_order_tensor_matmul
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! computes the double contraction of two fourth order tensors with minor
! symmetries as matrix product, input: 6x6 matrices 
! **********************************************************************
FUNCTION double_contraction_fourth_order_tensor_matmul(A_fot,B_fot)
  !
  USE Abaqus_Interface
  USE Constants
  !
  ! declaration of variables
  ! ()_fot_6x6.....symmetric fourth order tensors (stored in 6x6 matrix)
  !
  IMPLICIT NONE
  !
  INTEGER(KIND=ik)::i_count
  REAL(KIND=rk),DIMENSION(6,6),INTENT(IN)::A_fot,B_fot
  REAL(KIND=rk),DIMENSION(6,6)::A_fot_mod,B_fot_mod
  REAL(KIND=rk),DIMENSION(6)::temp
  REAL(KIND=rk),DIMENSION(6,6)::double_contraction_fourth_order_tensor_matmul
  !
  INTRINSIC MATMUL
  !
  A_fot_mod=A_fot; B_fot_mod=B_fot;
  DO i_count=FOUR_int,SIX_int
     temp=A_fot_mod(:,i_count)
     A_fot_mod(:,i_count)=temp*sqrt_two_r
     temp=B_fot_mod(i_count,:)
     B_fot_mod(i_count,:)=temp*sqrt_two_r
  END DO
  !
  double_contraction_fourth_order_tensor_matmul=MATMUL(A_fot_mod,&
       & B_fot_mod)
  !
END FUNCTION double_contraction_fourth_order_tensor_matmul
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! compute the double contraction of a fourth and a second order tensor
! **********************************************************************
FUNCTION double_contraction_fourth_order_tensor_second_order_tensor(&
     &A4_tensor,B2_tensor)
  !
  USE Abaqus_Interface
  USE Constants
  !
  IMPLICIT NONE
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::B2_tensor
  REAL(KIND=rk),DIMENSION(3,3,3,3),INTENT(IN)::A4_tensor
  REAL(KIND=rk),DIMENSION(3,3)::double_contraction_fourth_order_tensor_second_order_tensor
  REAL(KIND=rk),DIMENSION(3,3)::C2_tensor
  INTEGER(KIND=ik)::i_dcf,j_dcf,k_dcf,l_dcf
  !
  C2_tensor=ZERO_r
  DO i_dcf=ONE_int,THREE_int
     DO j_dcf=ONE_int,THREE_int
        DO k_dcf=ONE_int,THREE_int
           DO l_dcf=ONE_int,THREE_int
              C2_tensor(i_dcf,j_dcf)=C2_tensor(i_dcf,j_dcf)+&
                   &A4_tensor(i_dcf,j_dcf,k_dcf,l_dcf)*B2_tensor(k_dcf,l_dcf)
           END DO
        END DO
     END DO
  END DO
  !
  double_contraction_fourth_order_tensor_second_order_tensor=C2_tensor
  !
END FUNCTION double_contraction_fourth_order_tensor_second_order_tensor
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! compute the contraction of two fourth order tensors
! **********************************************************************
FUNCTION double_contraction_fourth_order_tensor(A4_tensor,B4_tensor)
  !
  USE Abaqus_Interface
  USE Constants
  !
  IMPLICIT NONE
  REAL(KIND=rk),DIMENSION(3,3,3,3),INTENT(IN)::A4_tensor,B4_tensor
  REAL(KIND=rk),DIMENSION(3,3,3,3)::C4_tensor
  REAL(KIND=rk),DIMENSION(3,3,3,3)::double_contraction_fourth_order_tensor
  INTEGER(KIND=ik)::m_dcfo,n_dcfo,i_dcfo,j_dcfo,k_dcfo,l_dcfo
  !
  C4_tensor=ZERO_r
  DO i_dcfo=ONE_int,THREE_int
     DO j_dcfo=ONE_int,THREE_int
        DO k_dcfo=ONE_int,THREE_int
           DO l_dcfo=ONE_int,THREE_int
              DO m_dcfo=ONE_int,THREE_int
                 DO n_dcfo=ONE_int,THREE_int
                    C4_tensor(i_dcfo,j_dcfo,k_dcfo,l_dcfo)=&
                         &C4_tensor(i_dcfo,j_dcfo,k_dcfo,l_dcfo)+&
                         &A4_tensor(i_dcfo,j_dcfo,m_dcfo,n_dcfo)*&
                         &B4_tensor(m_dcfo,n_dcfo,k_dcfo,l_dcfo)
                 END DO ! n
              END DO ! m
           END DO ! l
        END DO ! k
     END DO ! j
  END DO ! i
  double_contraction_fourth_order_tensor=C4_tensor
  !
END FUNCTION double_contraction_fourth_order_tensor
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! Compute the double contraction of two second order tensors
! **********************************************************************
FUNCTION double_contraction_second_order_tensors(A_tens,B_tens)
  !
  USE Abaqus_Interface
  USE Constants
  !
  ! declaration of variables
  !
  ! ()_tens......second order tensors
  !
  IMPLICIT NONE
  !
  INTEGER(KIND=ik)::i_count,j_count
  !
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::A_tens,B_tens
  REAL(KIND=rk)::contr_res
  REAL(KIND=rk)::double_contraction_second_order_tensors
  !
  contr_res=ZERO_r
  DO i_count=ONE_int,THREE_int
     DO j_count=ONE_int,THREE_int
        contr_res=contr_res+A_tens(i_count,j_count)*&
             & B_tens(i_count,j_count)
     END DO
  END DO
  !
  double_contraction_second_order_tensors=contr_res
  !
END FUNCTION double_contraction_second_order_tensors
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! Compute the dyadic product of two second order tensors
! **********************************************************************
FUNCTION DYADIC_PRODUCT_SECOND_ORDER_TENSORS(A_tens,B_tens)
  !
  USE Abaqus_Interface
  USE Constants
  !
  ! declaration of variables
  !
  ! (..)_tens.....second order tensors (input)
  !
  IMPLICIT NONE
  !
  INTEGER(KIND=ik)::i_DPSOT,j_DPSOT,k_DPSOT,l_DPSOT
  !
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::A_tens,B_tens
  REAL(KIND=rk),DIMENSION(3,3,3,3)::DYADIC_PRODUCT_SECOND_ORDER_TENSORS
  REAL(KIND=rk),DIMENSION(3,3,3,3)::A_dyad_B
  !
  A_dyad_B=ZERO_r
  !
  DO i_DPSOT=ONE_int,THREE_int
     DO j_DPSOT=ONE_int,THREE_int
        DO k_DPSOT=ONE_int,THREE_int
           DO l_DPSOT=ONE_int,THREE_int
              A_dyad_B(i_DPSOT,j_DPSOT,k_DPSOT,l_DPSOT)=&
                   &A_tens(i_DPSOT,j_DPSOT)*B_tens(k_DPSOT,l_DPSOT)
           END DO ! l
        END DO ! k
     END DO ! j
  END DO ! i
  !
  DYADIC_PRODUCT_SECOND_ORDER_TENSORS=A_dyad_B
  !
END FUNCTION DYADIC_PRODUCT_SECOND_ORDER_TENSORS
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! compute the Euler angles (Bungle convention) from a given rotation
! matrix
! based on: Bunge, "Texture analysis in materials science", Cuvillier
!           Verlag, Göttingen, 1993, p.21f.
! **********************************************************************
FUNCTION EulAngFromRotMat(RotMat)
  !
  USE Abaqus_Interface
  USE Constants
  !
  IMPLICIT NONE
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::RotMat
  REAL(KIND=rk),DIMENSION(3)::EulAngFromRotMat
  REAL(KIND=rk)::Phi,phi1,phi2,sphi
  !
  Phi=ACOS(RotMat(3,3))
  IF (ABS(ONE_r-RotMat(3,3)).LT.tol_angle) THEN
     Phi=ZERO_r
     ! note that this choice of phi1 and phi2 is arbitrary, because due
     ! to Ptolemy's identities (sum and difference formulas for sine and
     ! cosine) only the sum of phi1 and phi2 enters the computation
     ! zxz convention, thus if rotation about x is absent, only phi1+phi2
     ! can be determined
     phi1=ATAN2(RotMat(1,2),RotMat(1,1))
     phi2=ZERO_r
  ELSE
     sphi=SIN(Phi)
     phi2=ATAN2(RotMat(1,3)/sphi,RotMat(2,3)/sphi)
     phi1=ATAN2(RotMat(3,1)/sphi,-RotMat(3,2)/sphi)
  END IF
  !
  ! return angles in degree
  EulAngFromRotMat=(/phi1,Phi,phi2/)*(ONE_HUNDRED_EIGHTY_r/PI)
  !
END FUNCTION EulAngFromRotMat
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! expand a fourth order tensor (9x9 matrix notation) to tensor notation
! **********************************************************************
FUNCTION expand_fourth_order_tensor_matrix2tensor(A_tens_mat,env)
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:approp_index
  !
  IMPLICIT NONE
  REAL(KIND=rk),DIMENSION(9,9),INTENT(IN)::A_tens_mat
  REAL(KIND=rk),DIMENSION(3,3,3,3)::expand_fourth_order_tensor_matrix2tensor
  REAL(KIND=rk),DIMENSION(3,3,3,3)::A_tens
  INTEGER(KIND=ik)::i_efo,j_efo
  INTEGER(KIND=ik),DIMENSION(9,2)::index_t
  CHARACTER(4),INTENT(IN)::env
  !
  index_t=approp_index(NINE_int,env,'expand_fourth_order_tensor_matrix2tensor')
  !
  DO i_efo=ONE_int,NINE_int
     DO j_efo=ONE_int,NINE_int
        A_tens(index_t(i_efo,ONE_int),index_t(i_efo,TWO_int),index_t(j_efo,ONE_int),&
             &index_t(j_efo,TWO_int))=A_tens_mat(i_efo,j_efo)
     END DO
  END DO
  expand_fourth_order_tensor_matrix2tensor=A_tens
  !
END FUNCTION expand_fourth_order_tensor_matrix2tensor
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! expand a fourth order tensor with minor symmetries (6x6 matrix notation)
! to tensor notation
! **********************************************************************
FUNCTION expand_minor_symm_fourth_order_tensor_matrix2tensor(A6b6,env)
  !
  USE Abaqus_Interface
  USE Constants
  !
  IMPLICIT NONE
  REAL(KIND=rk),DIMENSION(6,6),INTENT(IN)::A6b6
  REAL(KIND=rk),DIMENSION(3,3,3,3)::A_tens
  REAL(KIND=rk),DIMENSION(3,3,3,3)::expand_minor_symm_fourth_order_tensor_matrix2tensor
  CHARACTER(4),INTENT(IN)::env
  !
  SELECT CASE (env)
  CASE ('FEAP') ! ordering 11,22,33,12,21,23,32,13,31
     A_tens(1,1,1,1)=A6b6(1,1)
     A_tens(1,1,2,2)=A6b6(1,2)
     A_tens(1,1,3,3)=A6b6(1,3)
     A_tens(1,1,1,2)=A6b6(1,4)
     A_tens(1,1,2,1)=A6b6(1,4)
     A_tens(1,1,2,3)=A6b6(1,5)
     A_tens(1,1,3,2)=A6b6(1,5)
     A_tens(1,1,1,3)=A6b6(1,6)
     A_tens(1,1,3,1)=A6b6(1,6)
     !
     A_tens(2,2,1,1)=A6b6(2,1)
     A_tens(2,2,2,2)=A6b6(2,2)
     A_tens(2,2,3,3)=A6b6(2,3)
     A_tens(2,2,1,2)=A6b6(2,4)
     A_tens(2,2,2,1)=A6b6(2,4)
     A_tens(2,2,2,3)=A6b6(2,5)
     A_tens(2,2,3,2)=A6b6(2,5)
     A_tens(2,2,1,3)=A6b6(2,6)
     A_tens(2,2,3,1)=A6b6(2,6)
     !
     A_tens(3,3,1,1)=A6b6(3,1)
     A_tens(3,3,2,2)=A6b6(3,2)
     A_tens(3,3,3,3)=A6b6(3,3)
     A_tens(3,3,1,2)=A6b6(3,4)
     A_tens(3,3,2,1)=A6b6(3,4)
     A_tens(3,3,2,3)=A6b6(3,5)
     A_tens(3,3,3,2)=A6b6(3,5)
     A_tens(3,3,1,3)=A6b6(3,6)
     A_tens(3,3,3,1)=A6b6(3,6)
     !
     A_tens(1,2,1,1)=A6b6(4,1)
     A_tens(1,2,2,2)=A6b6(4,2)
     A_tens(1,2,3,3)=A6b6(4,3)
     A_tens(1,2,1,2)=A6b6(4,4)
     A_tens(1,2,2,1)=A6b6(4,4)
     A_tens(1,2,2,3)=A6b6(4,5)
     A_tens(1,2,3,2)=A6b6(4,5)
     A_tens(1,2,1,3)=A6b6(4,6)
     A_tens(1,2,3,1)=A6b6(4,6)
     !
     A_tens(2,1,1,1)=A6b6(4,1)
     A_tens(2,1,2,2)=A6b6(4,2)
     A_tens(2,1,3,3)=A6b6(4,3)
     A_tens(2,1,1,2)=A6b6(4,4)
     A_tens(2,1,2,1)=A6b6(4,4)
     A_tens(2,1,2,3)=A6b6(4,5)
     A_tens(2,1,3,2)=A6b6(4,5)
     A_tens(2,1,1,3)=A6b6(4,6)
     A_tens(2,1,3,1)=A6b6(4,6)
     !
     A_tens(2,3,1,1)=A6b6(5,1)
     A_tens(2,3,2,2)=A6b6(5,2)
     A_tens(2,3,3,3)=A6b6(5,3)
     A_tens(2,3,1,2)=A6b6(5,4)
     A_tens(2,3,2,1)=A6b6(5,4)
     A_tens(2,3,2,3)=A6b6(5,5)
     A_tens(2,3,3,2)=A6b6(5,5)
     A_tens(2,3,1,3)=A6b6(5,6)
     A_tens(2,3,3,1)=A6b6(5,6)
     !
     A_tens(3,2,1,1)=A6b6(5,1)
     A_tens(3,2,2,2)=A6b6(5,2)
     A_tens(3,2,3,3)=A6b6(5,3)
     A_tens(3,2,1,2)=A6b6(5,4)
     A_tens(3,2,2,1)=A6b6(5,4)
     A_tens(3,2,2,3)=A6b6(5,5)
     A_tens(3,2,3,2)=A6b6(5,5)
     A_tens(3,2,1,3)=A6b6(5,6)
     A_tens(3,2,3,1)=A6b6(5,6)
     !
     A_tens(1,3,1,1)=A6b6(6,1)
     A_tens(1,3,2,2)=A6b6(6,2)
     A_tens(1,3,3,3)=A6b6(6,3)
     A_tens(1,3,1,2)=A6b6(6,4)
     A_tens(1,3,2,1)=A6b6(6,4)
     A_tens(1,3,2,3)=A6b6(6,5)
     A_tens(1,3,3,2)=A6b6(6,5)
     A_tens(1,3,1,3)=A6b6(6,6)
     A_tens(1,3,3,1)=A6b6(6,6)
     !
     A_tens(3,1,1,1)=A6b6(6,1)
     A_tens(3,1,2,2)=A6b6(6,2)
     A_tens(3,1,3,3)=A6b6(6,3)
     A_tens(3,1,1,2)=A6b6(6,4)
     A_tens(3,1,2,1)=A6b6(6,4)
     A_tens(3,1,2,3)=A6b6(6,5)
     A_tens(3,1,3,2)=A6b6(6,5)
     A_tens(3,1,1,3)=A6b6(6,6)
     A_tens(3,1,3,1)=A6b6(6,6)
     !
  CASE ('ABAQ') ! ordering 11,22,33,12,21,13,31,23,32
     A_tens(1,1,1,1)=A6b6(1,1)
     A_tens(1,1,2,2)=A6b6(1,2)
     A_tens(1,1,3,3)=A6b6(1,3)
     A_tens(1,1,1,2)=A6b6(1,4)
     A_tens(1,1,2,1)=A6b6(1,4)
     A_tens(1,1,1,3)=A6b6(1,5)
     A_tens(1,1,3,1)=A6b6(1,5)
     A_tens(1,1,2,3)=A6b6(1,6)
     A_tens(1,1,3,2)=A6b6(1,6)
     !
     A_tens(2,2,1,1)=A6b6(2,1)
     A_tens(2,2,2,2)=A6b6(2,2)
     A_tens(2,2,3,3)=A6b6(2,3)
     A_tens(2,2,1,2)=A6b6(2,4)
     A_tens(2,2,2,1)=A6b6(2,4)
     A_tens(2,2,1,3)=A6b6(2,5)
     A_tens(2,2,3,1)=A6b6(2,5)
     A_tens(2,2,2,3)=A6b6(2,6)
     A_tens(2,2,3,2)=A6b6(2,6)
     !
     A_tens(3,3,1,1)=A6b6(3,1)
     A_tens(3,3,2,2)=A6b6(3,2)
     A_tens(3,3,3,3)=A6b6(3,3)
     A_tens(3,3,1,2)=A6b6(3,4)
     A_tens(3,3,2,1)=A6b6(3,4)
     A_tens(3,3,1,3)=A6b6(3,5)
     A_tens(3,3,3,1)=A6b6(3,5)
     A_tens(3,3,2,3)=A6b6(3,6)
     A_tens(3,3,3,2)=A6b6(3,6)
     !
     A_tens(1,2,1,1)=A6b6(4,1)
     A_tens(1,2,2,2)=A6b6(4,2)
     A_tens(1,2,3,3)=A6b6(4,3)
     A_tens(1,2,1,2)=A6b6(4,4)
     A_tens(1,2,2,1)=A6b6(4,4)
     A_tens(1,2,1,3)=A6b6(4,5)
     A_tens(1,2,3,1)=A6b6(4,5)
     A_tens(1,2,2,3)=A6b6(4,6)
     A_tens(1,2,3,2)=A6b6(4,6)
     !
     A_tens(2,1,1,1)=A6b6(4,1)
     A_tens(2,1,2,2)=A6b6(4,2)
     A_tens(2,1,3,3)=A6b6(4,3)
     A_tens(2,1,1,2)=A6b6(4,4)
     A_tens(2,1,2,1)=A6b6(4,4)
     A_tens(2,1,1,3)=A6b6(4,5)
     A_tens(2,1,3,1)=A6b6(4,5)
     A_tens(2,1,2,3)=A6b6(4,6)
     A_tens(2,1,3,2)=A6b6(4,6)
     !
     A_tens(1,3,1,1)=A6b6(5,1)
     A_tens(1,3,2,2)=A6b6(5,2)
     A_tens(1,3,3,3)=A6b6(5,3)
     A_tens(1,3,1,2)=A6b6(5,4)
     A_tens(1,3,2,1)=A6b6(5,4)
     A_tens(1,3,1,3)=A6b6(5,5)
     A_tens(1,3,3,1)=A6b6(5,5)
     A_tens(1,3,2,3)=A6b6(5,6)
     A_tens(1,3,3,2)=A6b6(5,6)
     !
     A_tens(3,1,1,1)=A6b6(5,1)
     A_tens(3,1,2,2)=A6b6(5,2)
     A_tens(3,1,3,3)=A6b6(5,3)
     A_tens(3,1,1,2)=A6b6(5,4)
     A_tens(3,1,2,1)=A6b6(5,4)
     A_tens(3,1,1,3)=A6b6(5,5)
     A_tens(3,1,3,1)=A6b6(5,5)
     A_tens(3,1,2,3)=A6b6(5,6)
     A_tens(3,1,3,2)=A6b6(5,6)
     !
     A_tens(2,3,1,1)=A6b6(6,1)
     A_tens(2,3,2,2)=A6b6(6,2)
     A_tens(2,3,3,3)=A6b6(6,3)
     A_tens(2,3,1,2)=A6b6(6,4)
     A_tens(2,3,2,1)=A6b6(6,4)
     A_tens(2,3,1,3)=A6b6(6,5)
     A_tens(2,3,3,1)=A6b6(6,5)
     A_tens(2,3,2,3)=A6b6(6,6)
     A_tens(2,3,3,2)=A6b6(6,6)
     !
     A_tens(3,2,1,1)=A6b6(6,1)
     A_tens(3,2,2,2)=A6b6(6,2)
     A_tens(3,2,3,3)=A6b6(6,3)
     A_tens(3,2,1,2)=A6b6(6,4)
     A_tens(3,2,2,1)=A6b6(6,4)
     A_tens(3,2,1,3)=A6b6(6,5)
     A_tens(3,2,3,1)=A6b6(6,5)
     A_tens(3,2,2,3)=A6b6(6,6)
     A_tens(3,2,3,2)=A6b6(6,6)
     !
  CASE DEFAULT
     WRITE(16,'(A)')'Unknown env var. in: "expand_minor_symm_fourth_order_tensor_matrix2tensor"'
     WRITE(16,'(A)')'STOPPING COMPUTATION'
     STOP
  END SELECT
  expand_minor_symm_fourth_order_tensor_matrix2tensor=A_tens
  !
END FUNCTION expand_minor_symm_fourth_order_tensor_matrix2tensor
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! extract the symmtric part of A due to minor symmetry of A 
! input:  9 by 9 matrix format
! output: 6 by 6 matrix format
! **********************************************************************
FUNCTION extract_symm_fourth_ord_tens_minor_symm(A_tens_9b9)
  !
  USE Abaqus_Interface
  USE Constants
  !
  IMPLICIT NONE
  REAL(KIND=rk),DIMENSION(9,9),INTENT(IN)::A_tens_9b9
  REAL(KIND=rk),DIMENSION(6,6)::A_tens_6b6
  REAL(KIND=rk),DIMENSION(6,6)::extract_symm_fourth_ord_tens_minor_symm
  !
  ! upper left block
  A_tens_6b6(1:3,1:3)=A_tens_9b9(1:3,1:3)
  !
  ! diagonal terms (shear)
  A_tens_6b6(4,4)=A_tens_9b9(4,4)
  A_tens_6b6(5,5)=A_tens_9b9(6,6)
  A_tens_6b6(6,6)=A_tens_9b9(8,8)
  !
  ! off-diagonal terms (shear-normal)
  A_tens_6b6(4,1:3)=A_tens_9b9(4,1:3)
  A_tens_6b6(1:3,4)=A_tens_9b9(1:3,4)
  A_tens_6b6(5,1:3)=A_tens_9b9(6,1:3)
  A_tens_6b6(1:3,5)=A_tens_9b9(1:3,6)
  A_tens_6b6(6,1:3)=A_tens_9b9(8,1:3)
  A_tens_6b6(1:3,6)=A_tens_9b9(1:3,8)
  !
  ! off-diagonal terms (shear-shear)
  A_tens_6b6(4,5:6)=(/A_tens_9b9(4,6),A_tens_9b9(4,8)/)
  A_tens_6b6(5:6,4)=(/A_tens_9b9(6,4),A_tens_9b9(8,4)/)
  A_tens_6b6(5,6)=A_tens_9b9(6,8)
  A_tens_6b6(6,5)=A_tens_9b9(8,6)
  !
  extract_symm_fourth_ord_tens_minor_symm=A_tens_6b6
  !
END FUNCTION extract_symm_fourth_ord_tens_minor_symm
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! compute the left double contraction of a fourth order and a second
! order tensor C=B:A
! **********************************************************************
FUNCTION left_double_contraction_fourth_order_tensor_second_order_tensor(&
     &A4_tensor,B2_tensor)
  !
  USE Abaqus_Interface
  USE Constants
  !
  IMPLICIT NONE
  REAL(KIND=rk),DIMENSION(3,3,3,3),INTENT(IN)::A4_tensor
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::B2_tensor
  REAL(KIND=rk),DIMENSION(3,3)::left_double_contraction_fourth_order_tensor_second_order_tensor
  REAL(KIND=rk),DIMENSION(3,3)::C2_tensor
  INTEGER(KIND=ik)::i_ldcfo,j_ldcfo,k_ldcfo,l_ldcfo
  !
  C2_tensor=ZERO_r
  DO i_ldcfo=ONE_int,THREE_int
     DO j_ldcfo=ONE_int,THREE_int
        DO k_ldcfo=ONE_int,THREE_int
           DO l_ldcfo=ONE_int,THREE_int
              C2_tensor(i_ldcfo,j_ldcfo)=C2_tensor(i_ldcfo,j_ldcfo)+&
                   &B2_tensor(k_ldcfo,l_ldcfo)*&
                   &A4_tensor(k_ldcfo,l_ldcfo,i_ldcfo,j_ldcfo)
           END DO
        END DO
     END DO
  END DO
  left_double_contraction_fourth_order_tensor_second_order_tensor=C2_tensor
  !
END FUNCTION left_double_contraction_fourth_order_tensor_second_order_tensor
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! maps the components of a fourth order tensor to a 9x9 matrix format
! **********************************************************************
FUNCTION map_fourth_order_tensor_2_matrix_format(A_tens,env)
  !
  USE Abaqus_Interface
  USE Constants
  !
  IMPLICIT NONE
  REAL(KIND=rk),DIMENSION(3,3,3,3),INTENT(IN)::A_tens
  REAL(KIND=rk),DIMENSION(9,9)::A9b9
  REAL(KIND=rk),DIMENSION(9,9)::map_fourth_order_tensor_2_matrix_format
  CHARACTER(4),INTENT(IN)::env
  !
  SELECT CASE (env)
  CASE ('FEAP') ! ordering 11,22,33,12,21,23,32,13,31
     A9b9(1,1)=A_tens(1,1,1,1)
     A9b9(1,2)=A_tens(1,1,2,2)
     A9b9(1,3)=A_tens(1,1,3,3)
     A9b9(1,4)=A_tens(1,1,1,2)
     A9b9(1,5)=A_tens(1,1,2,1)
     A9b9(1,6)=A_tens(1,1,2,3)
     A9b9(1,7)=A_tens(1,1,3,2)
     A9b9(1,8)=A_tens(1,1,1,3)
     A9b9(1,9)=A_tens(1,1,3,1)
     !
     A9b9(2,1)=A_tens(2,2,1,1)
     A9b9(2,2)=A_tens(2,2,2,2)
     A9b9(2,3)=A_tens(2,2,3,3)
     A9b9(2,4)=A_tens(2,2,1,2)
     A9b9(2,5)=A_tens(2,2,2,1)
     A9b9(2,6)=A_tens(2,2,2,3)
     A9b9(2,7)=A_tens(2,2,3,2)
     A9b9(2,8)=A_tens(2,2,1,3)
     A9b9(2,9)=A_tens(2,2,3,1)
     !
     A9b9(3,1)=A_tens(3,3,1,1)
     A9b9(3,2)=A_tens(3,3,2,2)
     A9b9(3,3)=A_tens(3,3,3,3)
     A9b9(3,4)=A_tens(3,3,1,2)
     A9b9(3,5)=A_tens(3,3,2,1)
     A9b9(3,6)=A_tens(3,3,2,3)
     A9b9(3,7)=A_tens(3,3,3,2)
     A9b9(3,8)=A_tens(3,3,1,3)
     A9b9(3,9)=A_tens(3,3,3,1)
     !
     A9b9(4,1)=A_tens(1,2,1,1)
     A9b9(4,2)=A_tens(1,2,2,2)
     A9b9(4,3)=A_tens(1,2,3,3)
     A9b9(4,4)=A_tens(1,2,1,2)
     A9b9(4,5)=A_tens(1,2,2,1)
     A9b9(4,6)=A_tens(1,2,2,3)
     A9b9(4,7)=A_tens(1,2,3,2)
     A9b9(4,8)=A_tens(1,2,1,3)
     A9b9(4,9)=A_tens(1,2,3,1)
     !
     A9b9(5,1)=A_tens(2,1,1,1)
     A9b9(5,2)=A_tens(2,1,2,2)
     A9b9(5,3)=A_tens(2,1,3,3)
     A9b9(5,4)=A_tens(2,1,1,2)
     A9b9(5,5)=A_tens(2,1,2,1)
     A9b9(5,6)=A_tens(2,1,2,3)
     A9b9(5,7)=A_tens(2,1,3,2)
     A9b9(5,8)=A_tens(2,1,1,3)
     A9b9(5,9)=A_tens(2,1,3,1)
     !
     A9b9(6,1)=A_tens(2,3,1,1)
     A9b9(6,2)=A_tens(2,3,2,2)
     A9b9(6,3)=A_tens(2,3,3,3)
     A9b9(6,4)=A_tens(2,3,1,2)
     A9b9(6,5)=A_tens(2,3,2,1)
     A9b9(6,6)=A_tens(2,3,2,3)
     A9b9(6,7)=A_tens(2,3,3,2)
     A9b9(6,8)=A_tens(2,3,1,3)
     A9b9(6,9)=A_tens(2,3,3,1)
     !
     A9b9(7,1)=A_tens(3,2,1,1)
     A9b9(7,2)=A_tens(3,2,2,2)
     A9b9(7,3)=A_tens(3,2,3,3)
     A9b9(7,4)=A_tens(3,2,1,2)
     A9b9(7,5)=A_tens(3,2,2,1)
     A9b9(7,6)=A_tens(3,2,2,3)
     A9b9(7,7)=A_tens(3,2,3,2)
     A9b9(7,8)=A_tens(3,2,1,3)
     A9b9(7,9)=A_tens(3,2,3,1)
     !
     A9b9(8,1)=A_tens(1,3,1,1)
     A9b9(8,2)=A_tens(1,3,2,2)
     A9b9(8,3)=A_tens(1,3,3,3)
     A9b9(8,4)=A_tens(1,3,1,2)
     A9b9(8,5)=A_tens(1,3,2,1)
     A9b9(8,6)=A_tens(1,3,2,3)
     A9b9(8,7)=A_tens(1,3,3,2)
     A9b9(8,8)=A_tens(1,3,1,3)
     A9b9(8,9)=A_tens(1,3,3,1)
     !
     A9b9(9,1)=A_tens(3,1,1,1)
     A9b9(9,2)=A_tens(3,1,2,2)
     A9b9(9,3)=A_tens(3,1,3,3)
     A9b9(9,4)=A_tens(3,1,1,2)
     A9b9(9,5)=A_tens(3,1,2,1)
     A9b9(9,6)=A_tens(3,1,2,3)
     A9b9(9,7)=A_tens(3,1,3,2)
     A9b9(9,8)=A_tens(3,1,1,3)
     A9b9(9,9)=A_tens(3,1,3,1)
     !
  CASE ('ABAQ') ! ordering 11,22,33,12,21,13,31,23,32
     A9b9(1,1)=A_tens(1,1,1,1)
     A9b9(1,2)=A_tens(1,1,2,2)
     A9b9(1,3)=A_tens(1,1,3,3)
     A9b9(1,4)=A_tens(1,1,1,2)
     A9b9(1,5)=A_tens(1,1,2,1)
     A9b9(1,6)=A_tens(1,1,1,3)
     A9b9(1,7)=A_tens(1,1,3,1)
     A9b9(1,8)=A_tens(1,1,2,3)
     A9b9(1,9)=A_tens(1,1,3,2)
     !
     A9b9(2,1)=A_tens(2,2,1,1)
     A9b9(2,2)=A_tens(2,2,2,2)
     A9b9(2,3)=A_tens(2,2,3,3)
     A9b9(2,4)=A_tens(2,2,1,2)
     A9b9(2,5)=A_tens(2,2,2,1)
     A9b9(2,6)=A_tens(2,2,1,3)
     A9b9(2,7)=A_tens(2,2,3,1)
     A9b9(2,8)=A_tens(2,2,2,3)
     A9b9(2,9)=A_tens(2,2,3,2)
     !
     A9b9(3,1)=A_tens(3,3,1,1)
     A9b9(3,2)=A_tens(3,3,2,2)
     A9b9(3,3)=A_tens(3,3,3,3)
     A9b9(3,4)=A_tens(3,3,1,2)
     A9b9(3,5)=A_tens(3,3,2,1)
     A9b9(3,6)=A_tens(3,3,1,3)
     A9b9(3,7)=A_tens(3,3,3,1)
     A9b9(3,8)=A_tens(3,3,2,3)
     A9b9(3,9)=A_tens(3,3,3,2)
     !
     A9b9(4,1)=A_tens(1,2,1,1)
     A9b9(4,2)=A_tens(1,2,2,2)
     A9b9(4,3)=A_tens(1,2,3,3)
     A9b9(4,4)=A_tens(1,2,1,2)
     A9b9(4,5)=A_tens(1,2,2,1)
     A9b9(4,6)=A_tens(1,2,1,3)
     A9b9(4,7)=A_tens(1,2,3,1)
     A9b9(4,8)=A_tens(1,2,2,3)
     A9b9(4,9)=A_tens(1,2,3,2)
     !
     A9b9(5,1)=A_tens(2,1,1,1)
     A9b9(5,2)=A_tens(2,1,2,2)
     A9b9(5,3)=A_tens(2,1,3,3)
     A9b9(5,4)=A_tens(2,1,1,2)
     A9b9(5,5)=A_tens(2,1,2,1)
     A9b9(5,6)=A_tens(2,1,1,3)
     A9b9(5,7)=A_tens(2,1,3,1)
     A9b9(5,8)=A_tens(2,1,2,3)
     A9b9(5,9)=A_tens(2,1,3,2)
     !
     A9b9(6,1)=A_tens(1,3,1,1)
     A9b9(6,2)=A_tens(1,3,2,2)
     A9b9(6,3)=A_tens(1,3,3,3)
     A9b9(6,4)=A_tens(1,3,1,2)
     A9b9(6,5)=A_tens(1,3,2,1)
     A9b9(6,6)=A_tens(1,3,1,3)
     A9b9(6,7)=A_tens(1,3,3,1)
     A9b9(6,8)=A_tens(1,3,2,3)
     A9b9(6,9)=A_tens(1,3,3,2)
     !
     A9b9(7,1)=A_tens(3,1,1,1)
     A9b9(7,2)=A_tens(3,1,2,2)
     A9b9(7,3)=A_tens(3,1,3,3)
     A9b9(7,4)=A_tens(3,1,1,2)
     A9b9(7,5)=A_tens(3,1,2,1)
     A9b9(7,6)=A_tens(3,1,1,3)
     A9b9(7,7)=A_tens(3,1,3,1)
     A9b9(7,8)=A_tens(3,1,2,3)
     A9b9(7,9)=A_tens(3,1,3,2)
     !
     A9b9(8,1)=A_tens(2,3,1,1)
     A9b9(8,2)=A_tens(2,3,2,2)
     A9b9(8,3)=A_tens(2,3,3,3)
     A9b9(8,4)=A_tens(2,3,1,2)
     A9b9(8,5)=A_tens(2,3,2,1)
     A9b9(8,6)=A_tens(2,3,1,3)
     A9b9(8,7)=A_tens(2,3,3,1)
     A9b9(8,8)=A_tens(2,3,2,3)
     A9b9(8,9)=A_tens(2,3,3,2)
     !
     A9b9(9,1)=A_tens(3,2,1,1)
     A9b9(9,2)=A_tens(3,2,2,2)
     A9b9(9,3)=A_tens(3,2,3,3)
     A9b9(9,4)=A_tens(3,2,1,2)
     A9b9(9,5)=A_tens(3,2,2,1)
     A9b9(9,6)=A_tens(3,2,1,3)
     A9b9(9,7)=A_tens(3,2,3,1)
     A9b9(9,8)=A_tens(3,2,2,3)
     A9b9(9,9)=A_tens(3,2,3,2)
     !
  CASE DEFAULT
     WRITE(16,'(A)')'Unknown env var. in: "map_fourth_order_tensor_2_matrix_format"'
     WRITE(16,'(A)')'STOPPING COMPUTATION'
     STOP
  END SELECT
  map_fourth_order_tensor_2_matrix_format=A9b9
  !
END FUNCTION map_fourth_order_tensor_2_matrix_format
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! compute the symmetric pertubation of the deformation gradient to
! compute the material tangent (derivative of second Piola-Kirchhoff
! stress w.r.t. right Cauchy-Green tensor C, connecting material time
! derivative of S2PK and material time derivative of Green-Lagrange
! strain tensor or 1/2 \dot{C}) by means of numerical differentiation 
! proposed by: Miehe, Numerical computation of algorithmic (consistent) 
!     tangent moduli in large-strain computational inelasticity, CMAME, 
!     134, 1996, p.223-240
! **********************************************************************
FUNCTION pertub_F_numtang_material_tang(F_inv,pertub_ind)
  !
  USE Abaqus_Interface
  USE Constants
  !
  ! variable declaration
  ! F_inv........inverse of deformation gradient
  ! pertub_ind...indices to be pertubed
  IMPLICIT NONE
  INTEGER(KIND=ik),DIMENSION(2),INTENT(IN)::pertub_ind
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::F_inv
  REAL(KIND=rk),DIMENSION(3,3)::pertub
  REAL(KIND=rk),DIMENSION(3,3)::pertub_F_numtang_material_tang
  INTEGER(KIND=ik)::i_pfn,j_pfn
  !
  DO i_pfn=ONE_int,THREE_int
     DO j_pfn=ONE_int,THREE_int
        pertub(i_pfn,j_pfn)=(F_inv(pertub_ind(1),i_pfn)*&
             &IDENT(pertub_ind(2),j_pfn)+F_inv(pertub_ind(2),i_pfn)*&
             &IDENT(pertub_ind(1),j_pfn))
     END DO
  END DO
  !
  pertub_F_numtang_material_tang=pertub
  !
END FUNCTION pertub_F_numtang_material_tang
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! compute polar decomposition of deformation gradient
! based on: Simo, Hughes, Computational Inelasticity, Springer, 2000,
!           p. 243
! derivation relies on multiple use of Cayley-Hamilton theorem
! **********************************************************************
FUNCTION polar_decompostion_def_grad(DefGrad,NDI,NSHR,env)
  !
  USE Abaqus_Interface
  USE Abaqus_Interface_extended
  USE Constants
  USE auxiliary_routines_interface, only:transform_symm_tens_vec_form
  !
  IMPLICIT NONE
  INTEGER(KIND=ik),INTENT(IN)::NDI,NSHR
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::DefGrad
  REAL(KIND=rk),DIMENSION(3,3)::polar_decompostion_def_grad
  REAL(KIND=rk)::i1,i2,i3,D
  REAL(KIND=rk),DIMENSION(3)::C_princ,U_princ
  REAL(KIND=rk),DIMENSION(6)::C_vec
  REAL(KIND=rk),DIMENSION(3,3)::C,C_sq,U,U_inv
  !
  CHARACTER(4),INTENT(IN)::env
  !
  C=MATMUL(TRANSPOSE(DefGrad),DefGrad)
  C_sq=MATMUL(C,C)
  !
  ! compute principle values of right Cauchy-Green and right stretch
  CALL transform_symm_tens_vec_form(C,C_vec,ZERO_int,env)
  CALL SPRINC(C_vec,C_princ,ONE_int,NDI,NSHR)
  U_princ=SQRT(C_princ)
  !
  ! invariants of right stretch tensor
  i1=SUM(U_princ); i3=U_princ(1)*U_princ(2)*U_princ(3)
  i2=U_princ(1)*U_princ(2)+U_princ(1)*U_princ(3)+U_princ(2)*U_princ(3)
  !
  D=i1*i2-i3
  !
  ! compute the right stretch tensor and its inverse
  U=(-C_sq+(i1**(TWO_int)-i2)*C+i1*i3*IDENT)/D
  U_inv=(C-i1*U+i2*IDENT)/i3
  !
  polar_decompostion_def_grad=MATMUL(DefGrad,U_inv)
  !
END FUNCTION polar_decompostion_def_grad
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! Function to solve linear system of equations
! **********************************************************************
FUNCTION SOLVE_SYSTEM_OF_EQUATION(Mat,rhs,nrhs,prob_dim)
  !
  USE Abaqus_Interface
  USE LAPACK_utility_functions
  USE Constants
  !
  ! variable declaration
  ! --------------------
  ! Mat...........square matrix (prob_dim,prob_dim)
  ! rhs...........vector/matrix of right hand side (prob_dim,nrhs)
  ! nrhs..........number of right hand sides
  ! prob_dim......number of unknowns
  ! 
  IMPLICIT NONE
  !
  INTEGER(KIND=ik),INTENT(IN)::nrhs,prob_dim
  INTEGER(KIND=ik)::info,rank,info_GELSS,i_ssoe
  INTEGER(KIND=ik),PARAMETER::lwork=50
  INTEGER(KIND=ik),DIMENSION(prob_dim)::ipiv_SSOE
  !
  REAL(KIND=rk),DIMENSION(prob_dim,prob_dim),INTENT(IN)::Mat
  REAL(KIND=rk),DIMENSION(prob_dim,prob_dim)::Mat_work
  REAL(KIND=rk),DIMENSION(prob_dim,nrhs),INTENT(IN)::rhs
  REAL(KIND=rk),DIMENSION(prob_dim,nrhs)::rhs_work,SOLVE_SYSTEM_OF_EQUATION
  REAL(KIND=rk),DIMENSION(prob_dim)::sing_val
  REAL(KIND=rk),DIMENSION(lwork)::work
  !
  Mat_work=Mat
  rhs_work=rhs
  info=ONE_int
  !
  ! solve linear system of equations
  CALL GESV(prob_dim,nrhs,Mat_work,prob_dim,ipiv_SSOE,rhs_work,prob_dim,info)
  !
  IF (info .EQ. ZERO_int) THEN
     SOLVE_SYSTEM_OF_EQUATION=rhs_work
  ELSE
     !
     ! output for error handling
     WRITE(16,'(A)')'ill-conditioned system of equations,'
     WRITE(16,'(A)')'switching to SVD (SOLVE_SYSTEM_OF_EQUATION)'
     !
     Mat_work=Mat
     rhs_work=rhs
     info_GELSS=ONE_int
     !
     ! use SVD to solve ill-conditioned system of equation (possibly overdetermined)
     CALL GELSS(prob_dim,prob_dim,nrhs,Mat_work,prob_dim,rhs_work, &
          & prob_dim,sing_val,-eps_iter,rank,work,lwork,info_GELSS)
     !
     IF (info_GELSS .EQ. ZERO_int) THEN
        SOLVE_SYSTEM_OF_EQUATION=rhs_work
     ELSE
        WRITE(16,'(A)')'NOT ABLE TO OBTAIN A SOLUTION in: SOLVE_SYSTEM_OF_EQUATION,'
        WRITE(16,'(A)')'THINK ABOUT ERROR HANDLING'
        WRITE(16,'(A)')'STOP COMPUTATION'
        STOP
     END IF
     !     
  END IF
  !
END FUNCTION SOLVE_SYSTEM_OF_EQUATION
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! Function to solve ill-conditioned linear system of equations
! **********************************************************************
FUNCTION SOLVE_ILL_COND_SYSTEM_OF_EQUATION(Mat,rhs,nrhs,prob_dim,tol)
  !
  USE Abaqus_Interface
  USE LAPACK_utility_functions
  USE Constants
  !
  ! ill-conditioned system (overdetermined systems)
  ! www.sciencedirect.com/science/article/pii/S0893965900000860
  ! https://people.sc.fsu.edu/~jburkardt/f_src/qr_solve/qr_solve.html
  ! http://www.netlib.org/napack/
  !
  ! variable declaration
  ! --------------------
  ! Mat...........square matrix (prob_dim,prob_dim)
  ! rhs...........vector/matrix of right hand side (prob_dim,nrhs)
  ! nrhs..........number of right hand sides
  ! prob_dim......number of unknowns
  ! 
  IMPLICIT NONE
  !
  INTEGER(KIND=ik),INTENT(IN)::nrhs,prob_dim
  INTEGER(KIND=ik)::info,rank,info_GELSS
  INTEGER(KIND=ik),PARAMETER::lwork=150
  INTEGER(KIND=ik),DIMENSION(prob_dim)::ipiv_SSOE
  !
  REAL(KIND=rk),INTENT(IN)::tol
  REAL(KIND=rk),DIMENSION(prob_dim,prob_dim),INTENT(IN)::Mat
  REAL(KIND=rk),DIMENSION(prob_dim,prob_dim)::Mat_work
  REAL(KIND=rk),DIMENSION(prob_dim,nrhs),INTENT(IN)::rhs
  REAL(KIND=rk),DIMENSION(prob_dim,nrhs)::rhs_work
  REAL(KIND=rk),DIMENSION(prob_dim,nrhs)::SOLVE_ILL_COND_SYSTEM_OF_EQUATION
  REAL(KIND=rk),DIMENSION(prob_dim)::sing_val
  REAL(KIND=rk),DIMENSION(lwork)::work
  !
  Mat_work=Mat
  rhs_work=rhs
  info_GELSS=ONE_int
  !
  ! use SVD to solve ill-conditioned system of equation (possibly overdetermined)
  CALL GELSS(prob_dim,prob_dim,nrhs,Mat_work,prob_dim,rhs_work, &
       & prob_dim,sing_val,-eps_iter,rank,work,lwork,info_GELSS)
  !
  IF (info_GELSS .EQ. ZERO_int) THEN
     SOLVE_ILL_COND_SYSTEM_OF_EQUATION=rhs_work
  ELSE
     WRITE(16,'(A)')'NOT ABLE TO OBTAIN A SOLUTION in: SOLVE_ILL_COND_SYSTEM_OF_EQUATION,'
     WRITE(16,'(A)')'THINK ABOUT ERROR HANDLING'
     WRITE(16,'(A)')'STOP COMPUTATION'
     STOP
  END IF
  !
END FUNCTION SOLVE_ILL_COND_SYSTEM_OF_EQUATION

! **********************************************************************
!
! **********************************************************************
! rotation matrix based on Euler angles (Bunge convention)
! **********************************************************************
FUNCTION ROTATION_MATRIX_BUNGE(angles)
  !
  ! IMPORTANT NOTE: the definition of individual rotation matrices and
  ! the order R3*R2*R1 (fixed lab axis) corresponds to definition of
  ! Euler angles according to Bunge, "Texture analysis in materials
  ! science", Cuvillier Verlag, Göttingen, 1993
  !
  USE Abaqus_Interface
  USE Constants
  !
  ! declaration of variables
  !
  ! angles.........Euler angles defining the orientation
  !                [psi_deg,theta_deg,phi_deg]
  !
  IMPLICIT NONE
  !
  REAL(KIND=rk),DIMENSION(3),INTENT(IN)::angles
  REAL(KIND=rk),DIMENSION(3)::angles_t
  REAL(KIND=rk),DIMENSION(3,3)::ROTATION_MATRIX_BUNGE!,R1,R2,R3,
  REAL(KIND=rk)::psi_deg,theta_deg,phi_deg
  !
  INTRINSIC MATMUL
  angles_t=angles/ONE_HUNDRED_EIGHTY_r*PI
  psi_deg=angles_t(ONE_int)
  theta_deg=angles_t(TWO_int)
  phi_deg=angles_t(THREE_int)
  !
  ! R1=RESHAPE((/cos(psi_deg),-sin(psi_deg),ZERO_r,&
  !      &sin(psi_deg),cos(psi_deg),ZERO_r,&
  !      &ZERO_r,ZERO_r,ONE_r/),(/3,3/))
  ! !
  ! R2=RESHAPE((/ONE_r,ZERO_r,ZERO_r,&
  !      &ZERO_r,cos(theta_deg),-sin(theta_deg),&
  !      &ZERO_r,sin(theta_deg),cos(theta_deg)/),(/3,3/))
  ! !
  ! R3=RESHAPE((/cos(phi_deg),-sin(phi_deg),ZERO_r,&
  !      &sin(phi_deg),cos(phi_deg),ZERO_r,&
  !      &ZERO_r,ZERO_r,ONE_r/),(/3,3/))
  ! !
  ! ROTATION_MATRIX_BUNGE=MATMUL(R3,MATMUL(R2,R1))
  !
  ROTATION_MATRIX_BUNGE=RESHAPE((/&
       & -sin(phi_deg)*cos(theta_deg)*sin(psi_deg)+cos(phi_deg)*cos(psi_deg),& ! R11
       & -cos(phi_deg)*cos(theta_deg)*sin(psi_deg)-sin(phi_deg)*cos(psi_deg),& ! R21
       & sin(theta_deg)*sin(psi_deg),& ! R31
       & sin(phi_deg)*cos(theta_deg)*cos(psi_deg)+cos(phi_deg)*sin(psi_deg),& ! R12
       & cos(phi_deg)*cos(theta_deg)*cos(psi_deg)-sin(phi_deg)*sin(psi_deg),& ! R22
       & -sin(theta_deg)*cos(psi_deg),& ! R32
       & sin(phi_deg)*sin(theta_deg),& ! R13
       & cos(phi_deg)*sin(theta_deg),& ! R23
       & cos(theta_deg)/),(/3,3/)) ! R33
  !
END FUNCTION ROTATION_MATRIX_BUNGE
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! computes the tensor product of two second order tensors A, B such
! that (A dyad_underline B):C = A*C^(T)*B^(T); symbol: dyad underline
! **********************************************************************
FUNCTION spec_tensor_product_second_order_tensors (A_tens,B_tens,env)
  !
  ! NOTE: storage order in Abaqus: 11,22,33,12,13,23
  !       storage order in FEAP:   11,22,33,12,23,13
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:approp_index
  !
  ! declaration of variables
  !
  ! ()_tens......symmetric second order tensor
  !
  IMPLICIT NONE
  !
  INTEGER(KIND=ik)::m_stp,n_stp,i_stp,j_stp,k_stp,l_stp
  INTEGER(KIND=ik),DIMENSION(6,2)::index_t
  !
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::A_tens,B_tens
  REAL(KIND=rk),DIMENSION(6,6)::tens_prod
  REAL(KIND=rk),DIMENSION(6,6)::spec_tensor_product_second_order_tensors
  !
  CHARACTER(4),INTENT(IN)::env
  !
  index_t=approp_index(SIX_int,env,'spec_product_second_order_tensors')
  !
  DO m_stp=ONE_int,SIX_int
     DO n_stp=ONE_int,SIX_int
        i_stp=index_t(m_stp,ONE_int); j_stp=index_t(m_stp,TWO_int)
        k_stp=index_t(n_stp,ONE_int); l_stp=index_t(n_stp,TWO_int)
        tens_prod(m_stp,n_stp)=A_tens(i_stp,l_stp)*B_tens(j_stp,k_stp)
     END DO
  END DO
  !
  spec_tensor_product_second_order_tensors=tens_prod
  !
END FUNCTION spec_tensor_product_second_order_tensors
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! computes the tensor product of two second order tensors A, B such
! that (A dyad_underline B):C = A*C^(T)*B^(T); symbol: dyad underline
! **********************************************************************
FUNCTION spec_tensor_product_unsym_second_order_tensors(A_tens,B_tens,&
     &env)
  !
  ! NOTE: storage order in Abaqus: 11,22,33,12,21,13,31,23,32
  !       storage order in FEAP:   11,22,33,12,21,23,32,13,31
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:approp_index
  !
  ! declaration of variables
  !
  ! ()_tens......unsymmetric second order tensor
  !
  IMPLICIT NONE
  !
  INTEGER(KIND=ik)::m_stp,n_stp,i_stp,j_stp,k_stp,l_stp
  INTEGER(KIND=ik),DIMENSION(9,2)::index_t
  !
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::A_tens,B_tens
  REAL(KIND=rk),DIMENSION(9,9)::tens_prod
  REAL(KIND=rk),DIMENSION(9,9)::spec_tensor_product_unsym_second_order_tensors
  !
  CHARACTER(4),INTENT(IN)::env
  !
  index_t=approp_index(NINE_int,env,'spec_product_unsym_second_order_tensors')
  !
  DO m_stp=ONE_int,NINE_int
     DO n_stp=ONE_int,NINE_int
        i_stp=index_t(m_stp,ONE_int); j_stp=index_t(m_stp,TWO_int)
        k_stp=index_t(n_stp,ONE_int); l_stp=index_t(n_stp,TWO_int)
        tens_prod(m_stp,n_stp)=A_tens(i_stp,l_stp)*B_tens(j_stp,k_stp)
     END DO
  END DO
  !
  spec_tensor_product_unsym_second_order_tensors=tens_prod
  !
END FUNCTION spec_tensor_product_unsym_second_order_tensors
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! computes the square product of two unsymmetric second order tensors 
! A, B such that (A sp B):I4^(sym) = d/d_C {(A sp B):C} = d/d_C {A*C*B^T},
! where C is a symmetric second order tens
! output: 4-order tensor, coefficients stored in a 9x9 matrix format
! **********************************************************************
FUNCTION sq_prod_unsym_sec_ord_tens_dble_contr_fourth_order_sym_Ident(A,&
     &B,env)
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:approp_index
  !
  IMPLICIT NONE
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::A,B
  REAL(KIND=rk),DIMENSION(9,9)::sq_prod_unsym_sec_ord_tens_dble_contr_fourth_order_sym_Ident
  REAL(KIND=rk),DIMENSION(9,9)::sq_prod_dc_4o_sym_Ident
  INTEGER(KIND=ik)::m_spuso,n_spuso,i_spuso,j_spuso,k_spuso,l_spuso
  INTEGER(KIND=ik),DIMENSION(9,2)::index_t
  CHARACTER(4),INTENT(IN)::env
  !
  ! initialization
  sq_prod_dc_4o_sym_Ident=ZERO_r
  index_t=approp_index(NINE_int,env,'sq_prod_unsym_sec_ord_tens_dble_contr_fourth_order_sym_Ident')
  !
  DO m_spuso=ONE_r,NINE_int
     DO n_spuso=ONE_r,NINE_int
        i_spuso=index_t(m_spuso,ONE_int)
        j_spuso=index_t(m_spuso,TWO_int)
        k_spuso=index_t(n_spuso,ONE_int)
        l_spuso=index_t(n_spuso,TWO_int)
        sq_prod_dc_4o_sym_Ident(m_spuso,n_spuso)=ONE_HALF_r*(&
             &A(i_spuso,k_spuso)*B(j_spuso,l_spuso)+A(i_spuso,l_spuso)*&
             &B(j_spuso,k_spuso))
     END DO
  END DO
  sq_prod_unsym_sec_ord_tens_dble_contr_fourth_order_sym_Ident=sq_prod_dc_4o_sym_Ident
END FUNCTION sq_prod_unsym_sec_ord_tens_dble_contr_fourth_order_sym_Ident
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! computes the square product of two second order tensors A, B such
! that (A sp B):C = A*C*B^(T)
! **********************************************************************
FUNCTION square_product_second_order_tensors(A_tens,B_tens,env)
  !
  ! NOTE: storage order in Abaqus: 11,22,33,12,13,23
  !       storage order in FEAP:   11,22,33,12,23,13
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:approp_index
  !
  ! declaration of variables
  !
  ! ()_tens......symmetric second order tensor
  !
  IMPLICIT NONE
  !
  INTEGER(KIND=ik)::m_stp,n_stp,i_stp,j_stp,k_stp,l_stp
  INTEGER(KIND=ik),DIMENSION(6,2)::index_t
  !
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::A_tens,B_tens
  REAL(KIND=rk),DIMENSION(6,6)::tens_prod
  REAL(KIND=rk),DIMENSION(6,6)::square_product_second_order_tensors
  !
  CHARACTER(4),INTENT(IN)::env
  !
  index_t=approp_index(SIX_int,env,'square_product_second_order_tensors')
  !
  DO m_stp=ONE_int,SIX_int
     DO n_stp=ONE_int,SIX_int
        i_stp=index_t(m_stp,ONE_int); j_stp=index_t(m_stp,TWO_int)
        k_stp=index_t(n_stp,ONE_int); l_stp=index_t(n_stp,TWO_int)
        tens_prod(m_stp,n_stp)=A_tens(i_stp,k_stp)*B_tens(j_stp,l_stp)
     END DO
  END DO
  !
  square_product_second_order_tensors=tens_prod
  !
END FUNCTION square_product_second_order_tensors
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! computes the square product of two second order tensors A, B such that
! (A sp B):C = A*C*B^(T)
! **********************************************************************
FUNCTION square_product_unsym_second_order_tensors(A_tens,B_tens,env)
  !
  ! NOTE: storage order in Abaqus: 11,22,33,12,21,13,31,23,32
  !       storage order in FEAP:   11,22,33,12,21,23,32,13,31
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:approp_index
  !
  ! declaration of variables
  !
  ! ()_tens......unsymmetric second order tensor
  !
  IMPLICIT NONE
  !
  INTEGER(KIND=ik)::m_stp,n_stp,i_stp,j_stp,k_stp,l_stp
  INTEGER(KIND=ik),DIMENSION(9,2)::index_t
  !
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::A_tens,B_tens
  REAL(KIND=rk),DIMENSION(9,9)::tens_prod
  REAL(KIND=rk),DIMENSION(9,9)::square_product_unsym_second_order_tensors
  !
  CHARACTER(4),INTENT(IN)::env
  !
  index_t=approp_index(NINE_int,env,'square_product_unsym_second_order_tensors')
  !
  DO m_stp=ONE_int,NINE_int
     DO n_stp=ONE_int,NINE_int
        i_stp=index_t(m_stp,ONE_int); j_stp=index_t(m_stp,TWO_int)
        k_stp=index_t(n_stp,ONE_int); l_stp=index_t(n_stp,TWO_int)
        tens_prod(m_stp,n_stp)=A_tens(i_stp,k_stp)*B_tens(j_stp,l_stp)
     END DO
  END DO
  !
  square_product_unsym_second_order_tensors=tens_prod
  !
END FUNCTION square_product_unsym_second_order_tensors
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! symmetric tensor product of two second order tensors A and B, such that
! ((A sp B) + (A dyad_underline B)):C = A*C*B^(T) + A*C^(T)*B^(T)
! **********************************************************************
FUNCTION symm_tensor_product_second_order_tensors (A_tens,B_tens,env)
  !
  ! NOTE: storage order in Abaqus: 11,22,33,12,13,23
  !       storage order in FEAP:   11,22,33,12,23,13
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:spec_tensor_product_second_order_tensors,&
       & square_product_second_order_tensors
  !
  ! declaration of variables
  !
  ! ()_tens......symmetric second order tensor
  !
  IMPLICIT NONE
  !
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::A_tens,B_tens
  REAL(KIND=rk),DIMENSION(6,6)::symm_tensor_product_second_order_tensors
  !
  CHARACTER(4),INTENT(IN)::env
  !
  symm_tensor_product_second_order_tensors=&
       & square_product_second_order_tensors(A_tens,B_tens,env) + &
       & spec_tensor_product_second_order_tensors(A_tens,B_tens,env)
  !
END FUNCTION symm_tensor_product_second_order_tensors
! ......................................................................
! **********************************************************************
! END                    Auxiliary_routines
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!                      Constitutive_routines
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! **********************************************************************
! adjust parameters controlling the increase and decrease of the visco-
! plastic parameters
! **********************************************************************
SUBROUTINE adjust_visco_parameter_miehe_vp_log (p_exp_vp_t,tau_relax_vp_t,&
     & mode,p_exp_vp_target,tau_relax_vp_target)
  !
  USE Abaqus_Interface
  USE Constants
  !
  ! declaration of variables
  ! p_exp_vp_t......rate sensitivity exponent temp. val. (input/update)
  ! tau_relax_vp_t..relaxation time temporary value (input and update)
  ! _target.........prescribed values of rate sensitivity exponent and 
  !                 relaxation time
  ! mode............mode of operation
  !
  IMPLICIT NONE
  REAL(KIND=rk),INTENT(IN)::p_exp_vp_target,tau_relax_vp_target
  REAL(KIND=rk),INTENT(INOUT)::p_exp_vp_t,tau_relax_vp_t
  REAL(KIND=rk)::alpha,beta,alpha_bar,beta_bar
  !
  CHARACTER(17),INTENT(IN)::mode
  !
  ! decrease and increase of rate-sensitivity exponent
  alpha=ONE_by_FIVE_r; alpha_bar=TWO_r
  ! increase and decrease of characteristic time tau_relax_vp 
  beta=FIVE_r; beta_bar=ONE_HALF_r
  !
  SELECT CASE (mode)
  CASE ('incr_exp_decr_tau')
     p_exp_vp_t=min(alpha_bar*p_exp_vp_t,p_exp_vp_target)
     tau_relax_vp_t=max(beta_bar*tau_relax_vp_t,tau_relax_vp_target)
  CASE ('decr_exp_incr_tau')
     p_exp_vp_t=max(alpha*p_exp_vp_t,TWO_r)
     tau_relax_vp_t=min(beta*tau_relax_vp_t,up_bound_relax_time)
  CASE DEFAULT
     WRITE(16,'(A)')'Unknown mode of operation used in adjust_visco_parameter_miehe_vp_log'
  END SELECT
  !
END SUBROUTINE adjust_visco_parameter_miehe_vp_log
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! assemble the residual vector and the jacobian for newton method to
! solve the discrete (augmented) consistency condition or viscoplastic
! residual equations
! **********************************************************************
SUBROUTINE assembly_residual_jacobian_newton (residual_nw,jacobian_nw,&
     &active_slip_system,deriv_visco_plast_potent,Phi_alpha,&
     &dPhi_alpha_d_delta_lambda_ip1_beta,delta_lambda_i,delta_lambda_ip1,&
     &plastic_slip,algo,dt,eta,alpha_slip,mat_param_inelast,theta_KK_t,&
     &hard_law)
  !
  USE Abaqus_Interface
  USE Constants
  USE constitutive_routines_interface, only:extract_hardening_parameters_crystal_plast,&
       & extract_viscoplastic_parameters_crystal_plast, &
       & hardening_function_cailletaud_forest,&
       & hardening_function_schmidt_baldassari
  !
  ! declaration of variables
  ! residual_nw.................residual vector
  ! jacobian_nw.................jacobian matrix
  ! deriv_visco_plast_potent....derivative of the viscoplast potential
  !                             aux. var in rate-independent case
  ! active_slip_system..........list of active slip systems
  ! Phi_alpha...................value yield funct. of alpha slip system
  ! dPhi_alpha_d_delta_lambda_ip1_beta...derivative of alpha yield funct.
  !                                      w.r.t. beta incr. lagr. multipl
  ! delta_lambda_...............incremental lagrange multipliers for
  !                             each system
  ! plastic_slip................accumulated plastic slip on each system
  ! algo........................id for algorithm employed
  !                             0...augmented-lagrangian
  !                             1...NCP-function (Kanzow-Kleinmichel)
  !                             2...Perzyna viscoplastic formulation
  !                             3...Ortiz/Miehe viscoplastic formulation
  ! dt..........................time increment
  ! eta.........................pseudo-viscosity in augmented lagrangian
  ! hard_law....................id of hardening law
  !                    0....nonlin. Taylor-type (Schmidt-Baldassari)
  !                    1....nonlin. interaction matrix (Cailletaud,Forest)
  !
  IMPLICIT NONE
  REAL(KIND=rk),DIMENSION(24),INTENT(INOUT)::residual_nw
  REAL(KIND=rk),DIMENSION(24),INTENT(INOUT)::deriv_visco_plast_potent
  REAL(KIND=rk),DIMENSION(24,24),INTENT(INOUT)::jacobian_nw
  REAL(KIND=rk),INTENT(IN)::Phi_alpha,dt,eta,theta_KK_t
  REAL(KIND=rk),DIMENSION(13),INTENT(IN)::mat_param_inelast
  REAL(KIND=rk),DIMENSION(24),INTENT(IN)::dPhi_alpha_d_delta_lambda_ip1_beta
  REAL(KIND=rk),DIMENSION(24),INTENT(IN)::delta_lambda_i,delta_lambda_ip1
  REAL(KIND=rk),DIMENSION(24),INTENT(IN)::plastic_slip
  REAL(KIND=rk)::val_max_funct,aux_funct,p_exp_vp
  REAL(KIND=rk)::tau_relax_vp,K_norm_stress_vp,tau_crss_alpha
  REAL(KIND=rk)::Y_alpha,aux_aux_funct,Y
  REAL(KIND=rk),DIMENSION(3)::mat_param_visco,hard_param_sb
  REAL(KIND=rk),DIMENSION(10)::hard_param_cf
  REAL(KIND=rk),DIMENSION(24)::dY_alpha_d_delta_lambda_beta
  !
  LOGICAL,DIMENSION(24),INTENT(INOUT)::active_slip_system
  INTEGER(KIND=ik),INTENT(IN)::algo,alpha_slip,hard_law
  !
  SELECT CASE (algo)
  CASE (ZERO_int) ! augmented lagrangian approach
     ! residual
     val_max_funct=max(ZERO_r,delta_lambda_i(alpha_slip)+dt*eta*Phi_alpha)
     residual_nw(alpha_slip)=delta_lambda_ip1(alpha_slip)-val_max_funct
     !
     ! jacobian
     IF (val_max_funct > eps_one) THEN
        active_slip_system(alpha_slip)=.TRUE.
        jacobian_nw(alpha_slip,:)=jacobian_nw(alpha_slip,:)-dt*eta*&
             &dPhi_alpha_d_delta_lambda_ip1_beta
     END IF
     !
  CASE (ONE_int)  ! NCP-function (Kanzow-Kleinmichel)
     ! residual
     aux_funct=sqrt(Phi_alpha**TWO_int+(TWO_r-theta_KK_t)*Phi_alpha*&
          & delta_lambda_ip1(alpha_slip)+delta_lambda_ip1(alpha_slip)&
          & **TWO_int)
     residual_nw(alpha_slip)=aux_funct+Phi_alpha-delta_lambda_ip1(alpha_slip)
     deriv_visco_plast_potent(alpha_slip)=((Phi_alpha+ &
          & delta_lambda_ip1(alpha_slip))/aux_funct+ONE_r)
     ! jacobian
     jacobian_nw(alpha_slip,:)=jacobian_nw(alpha_slip,:)+((Phi_alpha+&
          &(ONE_r-ONE_HALF_r*theta_KK_t)*delta_lambda_ip1(alpha_slip))&
          &/aux_funct+ONE_r)*dPhi_alpha_d_delta_lambda_ip1_beta
     jacobian_nw(alpha_slip,alpha_slip)=jacobian_nw(alpha_slip,alpha_slip)&
          &+((delta_lambda_ip1(alpha_slip)+(ONE_r-ONE_HALF_r*theta_KK_t)*&
          &Phi_alpha)/aux_funct-ONE_r)
     IF (delta_lambda_ip1(alpha_slip) > eps_one) THEN
        active_slip_system(alpha_slip)=.TRUE.
     END IF
     !
  CASE (TWO_int)  ! Perzyna viscoplastic formulation
     ! extract material parameters for rate-dependent behavior
     mat_param_visco=extract_viscoplastic_parameters_crystal_plast(&
          &mat_param_inelast)
     p_exp_vp=mat_param_visco(1); tau_relax_vp=mat_param_visco(2)
     K_norm_stress_vp=mat_param_visco(3)
     !
     ! residual
     val_max_funct=max(ZERO_r,Phi_alpha/K_norm_stress_vp)
     residual_nw(alpha_slip)=delta_lambda_ip1(alpha_slip)-&
          &dt/tau_relax_vp*val_max_funct**p_exp_vp
     !
     ! jacobian
     IF (val_max_funct > eps_one) THEN
        active_slip_system(alpha_slip)=.TRUE.
        aux_funct=dt/tau_relax_vp*p_exp_vp/K_norm_stress_vp*&
             &((val_max_funct)**(p_exp_vp-ONE_r))
        deriv_visco_plast_potent(alpha_slip)=aux_funct
        jacobian_nw(alpha_slip,:)=jacobian_nw(alpha_slip,:)-aux_funct*&
             &dPhi_alpha_d_delta_lambda_ip1_beta
     END IF
     !
  CASE (THREE_int)! Ortiz/Miehe viscoplastic formulation
     ! extract material parameters for rate-dependent behavior
     mat_param_visco=extract_viscoplastic_parameters_crystal_plast(&
          &mat_param_inelast)
     p_exp_vp=mat_param_visco(1); tau_relax_vp=mat_param_visco(2)
     !
     ! current hardening (change in yield stress) and its derivative 
     ! w.r.t. increments in lagrange multipliers
     SELECT CASE (hard_law)
     CASE (ZERO_int) ! nonlinear Taylor-type hardening
        hard_param_sb=extract_hardening_parameters_crystal_plast(&
             &mat_param_inelast,hard_law,THREE_int)
        CALL hardening_function_schmidt_baldassari(Y_alpha,&
             &dY_alpha_d_delta_lambda_beta,plastic_slip,hard_param_sb)
        Y=hard_param_sb(1)
        !
     CASE (ONE_int)  ! nonlinear hardening based on interaction matrix
        hard_param_cf=extract_hardening_parameters_crystal_plast(&
             &mat_param_inelast,hard_law,TEN_int)
        CALL hardening_function_cailletaud_forest(Y_alpha,&
             &dY_alpha_d_delta_lambda_beta,plastic_slip,alpha_slip,&
             &hard_param_cf)
        Y=hard_param_cf(1)
        !
     CASE DEFAULT
        WRITE(16,'(A)')'Unknown value of hard_law in: "assembly_residual_jacobian_newton"'
        WRITE(16,'(A)')'STOPPING COMPUTATION'
        STOP
     END SELECT
     !
     ! residual
     tau_crss_alpha=Y+Y_alpha
     val_max_funct=max(ZERO_r,Phi_alpha)
     aux_funct=max(lb_invert,ONE_r+delta_lambda_ip1(alpha_slip)*&
          &tau_relax_vp/dt)
     aux_aux_funct=val_max_funct/tau_crss_alpha+ONE_r
     residual_nw(alpha_slip)=log(aux_funct)-p_exp_vp*log(aux_aux_funct)
     !
     ! jacobian
     jacobian_nw(alpha_slip,alpha_slip)=jacobian_nw(alpha_slip,alpha_slip)&
          &+ONE_r/aux_funct*tau_relax_vp/dt
     IF (val_max_funct > eps_one) THEN
        active_slip_system(alpha_slip)=.TRUE.
        deriv_visco_plast_potent(alpha_slip)=p_exp_vp/aux_aux_funct
        jacobian_nw(alpha_slip,:)=jacobian_nw(alpha_slip,:)-&
             &deriv_visco_plast_potent(alpha_slip)*(&
             &dPhi_alpha_d_delta_lambda_ip1_beta*tau_crss_alpha-&
             &val_max_funct*dY_alpha_d_delta_lambda_beta)/&
             &tau_crss_alpha**TWO_int
     END IF
  CASE DEFAULT
     WRITE(16,'(A)')'Unknown algorithm id used in "assembly_residual_jacobian_newton"'
  END SELECT
  !
END SUBROUTINE assembly_residual_jacobian_newton
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! Determine the number of distinct eigenvalues
! **********************************************************************
SUBROUTINE check_eigval_distinct(n_distinct,index_distinct,lambda,tol)
  !
  USE Abaqus_Interface
  USE Constants
  !
  ! declaration of variables
  !
  ! lambda...........vector of eigenvalues
  ! tol..............tolerance to equality of eigenvalues
  ! index_distinct...index of the distinct eigenvalues in lambda
  ! n_distinct.......number of distinct eigenvalues
  !
  IMPLICIT NONE
  !
  INTEGER(KIND=ik),INTENT(OUT)::n_distinct
  INTEGER(KIND=ik),DIMENSION(3),INTENT(OUT)::index_distinct
  !
  REAL(KIND=rk),DIMENSION(3),INTENT(IN)::lambda
  REAL(KIND=rk),INTENT(IN)::tol
  !
  INTRINSIC ABS
  !
  ! three distinct eigenvalues
  IF (ABS(lambda(ONE_int)-lambda(TWO_int)) > tol .AND. &
       & ABS(lambda(ONE_int)-lambda(THREE_int)) > tol .AND. &
       & ABS(lambda(THREE_int)-lambda(TWO_int)) > tol) THEN
     index_distinct=(/ONE_int,TWO_int,THREE_int/)
     n_distinct=THREE_int
     ! two distinct eigenvalues
  ELSE IF (ABS(lambda(ONE_int)-lambda(TWO_int)) <= tol .AND. &
       & ABS(lambda(ONE_int)-lambda(THREE_int)) > tol) THEN
     index_distinct=(/ONE_int,THREE_int,THREE_int/)
     n_distinct=TWO_int
  ELSE IF (ABS(lambda(TWO_int)-lambda(THREE_int)) <= tol .AND. &
       & ABS(lambda(ONE_int)-lambda(TWO_int)) > tol) THEN
     index_distinct=(/ONE_int,TWO_int,TWO_int/)
     n_distinct=TWO_int
  ELSE IF (ABS(lambda(ONE_int)-lambda(THREE_int)) <= tol .AND. &
       & ABS(lambda(TWO_int)-lambda(THREE_int)) > tol) THEN
     index_distinct=(/TWO_int,THREE_int,THREE_int/)
     n_distinct=TWO_int
     ! one eigenvalue
  ELSE IF (ABS(lambda(ONE_int)-lambda(TWO_int)) <= tol .AND. &
       & ABS(lambda(ONE_int)-lambda(THREE_int)) <= tol) THEN
     index_distinct=(/ONE_int,ONE_int,ONE_int/)
     n_distinct=ONE_int
  END IF
  !
END SUBROUTINE check_eigval_distinct
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! compute coefficients for projection operators
! **********************************************************************
SUBROUTINE coefficients_for_projection(theta,xi,eta,e_vect,d_vect,&
     &lambda,index_distinct,n_distinct)
  !
  USE Abaqus_Interface
  USE Constants
  !
  ! declaration of variables
  !
  ! lambda...........vector of eigenvalues
  ! index_distinct...index of the distinct eigenvalues in lambda
  ! n_distinct.......number of distinct eigenvalues
  ! ()_vect..........diagonal functions and their 1st derivative
  ! theta,xi,eta.....coefficients
  !
  IMPLICIT NONE
  !
  INTEGER(KIND=ik),INTENT(IN)::n_distinct
  INTEGER(KIND=ik),DIMENSION(3),INTENT(IN)::index_distinct
  INTEGER(KIND=ik)::s_cfp,t_cfp,u_cfp
  !
  REAL(KIND=rk),DIMENSION(3),INTENT(IN)::lambda,e_vect,d_vect
  REAL(KIND=rk),DIMENSION(3,3),INTENT(OUT)::theta,xi
  REAL(KIND=rk),INTENT(OUT)::eta
  !
  eta=ZERO_r
  IF (n_distinct == THREE_int) THEN
     DO s_cfp=ONE_int,THREE_int
        DO t_cfp=ONE_int,THREE_int
           IF (s_cfp /= t_cfp) THEN
              DO u_cfp=ONE_int,THREE_int
                 IF ((u_cfp /= t_cfp).AND.(u_cfp /= s_cfp)) THEN
                    eta=eta+e_vect(s_cfp)/(TWO_r*(lambda(s_cfp)-&
                         &lambda(t_cfp))*(lambda(s_cfp)-lambda(u_cfp)))
                 END IF
              END DO ! u
           END IF
        END DO ! t
     END DO ! s
  END IF
  !
  theta=ZERO_r; xi=ZERO_r
  DO s_cfp=ONE_int,n_distinct
     DO t_cfp=ONE_int,n_distinct
        IF (s_cfp /= t_cfp) THEN
           theta(s_cfp,t_cfp)=(e_vect(s_cfp)-e_vect(t_cfp))/&
                &(lambda(index_distinct(s_cfp))-&
                &lambda(index_distinct(t_cfp)))
           xi(s_cfp,t_cfp)=(theta(s_cfp,t_cfp)-ONE_HALF_r*&
                &d_vect(t_cfp))/(lambda(index_distinct(s_cfp))-&
                &lambda(index_distinct(t_cfp)))
        END IF
     END DO ! t
  END DO ! s
  !
END SUBROUTINE coefficients_for_projection
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! extension of crystal plasticity model of Schmidt-Baldassari, CMAME
! (192), 2003
! including anisotropic elasticity and nonlinear latent and cross
! hardening based on interaction matrix
! identical to Casal, Forest, CMS (45), 2009
! **********************************************************************
SUBROUTINE crystal_plasticity_augm_lagr(S2PK,C_tang,SDV,F,dt,&
     &elast_const,mat_param_inelast,m_SH,symm_class,orth_basis,&
     &orth_proj,algo,numtang,NDI,NSHR,NSTATV,env,hard_law,theta_KK,&
     &conv_issue_dil)
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:transform_tens_vec_form, &
       & polar_decompostion_def_grad, EulAngFromRotMat
  USE constitutive_routines_interface, only: CSDA_material_tangent_single_crystal_plast,&
       & consistent_algorithmic_tangent_augmented_lagrangian_penalty,&
       & stress_computation_single_crystal_plasticity
  !
  ! variable declaration
  ! S2PK............second Piola-Kirchhoff stress (ref. conf.)
  ! C_tang..........algorithmic tangent in material setting
  ! SDV.............vector of state dependent variables (internal var.)
  ! F...............deformation gradient (t_{n+1})
  ! dt..............time increment
  ! elast_const.....material parameter for elastic behavior
  ! mat_param_inelast..material parameter assoc. with inelastic behavior
  ! m_SH............Seth-Hill parameter
  ! symm_class......symmetry included in elastic law
  ! orth_basis......all second order, orthogonal basis N(i), rotated in
  !                 case of cubic elasticity
  ! orth_proj.......symmetric dyadic products of orthogonal basis i.e.
  !                 N(i) dyadic N(i)
  ! algo............id for algorithm employed
  !                 0...augmented-lagrangian
  !                 1...NCP-function (Kanzow-Kleinmichel)
  !                 2...Perzyna viscoplastic formulation
  !                 3...Ortiz/Miehe viscoplastic formulation
  ! numtang.........id for computation of consistent tangent
  !                 0...analytical
  !                 1...CSDA (currently not implemented, definition of
  !                     complex variables required)
  !                 2...FD
  ! hard_law........id of hardening law
  !                 0....nonlin. Taylor-type (Schmidt-Baldassari)
  !                 1....nonlin. interaction matrix (Cailletaud,Forest)
  ! theta_KK........algorithmic parameter for Kanzow-Kleinmichel
  !                 NCP-function
  ! conv_issue_dil..error flag indicating non-convergence in
  !                 determination of incremental lagr. multiplier
  IMPLICIT NONE
  INTEGER(KIND=ik),INTENT(IN)::algo,symm_class,NDI,NSHR,numtang
  INTEGER(KIND=ik),INTENT(IN)::hard_law,NSTATV
  REAL(KIND=rk),DIMENSION(3,3),INTENT(OUT)::S2PK
  REAL(KIND=rk),DIMENSION(6,6),INTENT(OUT)::C_tang
  REAL(KIND=rk),DIMENSION(NSTATV),INTENT(INOUT)::SDV
  REAL(KIND=rk),INTENT(IN)::dt,m_SH,theta_KK
  REAL(KIND=rk),DIMENSION(3),INTENT(IN)::elast_const
  REAL(KIND=rk),DIMENSION(13),INTENT(IN)::mat_param_inelast
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::F
  REAL(KIND=rk),DIMENSION(6,3,3),INTENT(IN)::orth_basis
  REAL(KIND=rk),DIMENSION(6,6,6),INTENT(IN)::orth_proj
  REAL(KIND=rk)::Je_n1,det_Fp_inv,aver_adj_stp,aver_adj_stp_save
  REAL(KIND=rk),DIMENSION(3,3)::F_p_inv_n,F_p_inv_n1,S2PK_IC,Ce_n1
  REAL(KIND=rk),DIMENSION(3,3)::Ce_n1_inv,C,F_p_inv_iter,Fe_n1,Re_n1
  REAL(KIND=rk),DIMENSION(24)::plastic_slip_n,plastic_slip_iter
  REAL(KIND=rk),DIMENSION(24,6)::d_delta_lambda_ip1_d_C
  LOGICAL,INTENT(OUT)::conv_issue_dil
  LOGICAL,DIMENSION(24)::active_slip_system
  !
  CHARACTER(4),INTENT(IN)::env
  !    
  ! internal variables at t_n
  CALL transform_tens_vec_form(F_p_inv_n,SDV(1:9),ONE_int,env)
  plastic_slip_n=SDV(10:33)
  aver_adj_stp=MAX(ONE_r,SDV(NSTATV))
  !
  ! compute stress and internal variables
  CALL stress_computation_single_crystal_plasticity(S2PK,S2PK_IC,&
       &F_p_inv_n1,plastic_slip_iter,d_delta_lambda_ip1_d_C,&
       &active_slip_system,Ce_n1,Ce_n1_inv,Je_n1,C,det_Fp_inv,&
       &F_p_inv_iter,F,dt,F_p_inv_n,plastic_slip_n,elast_const,&
       &mat_param_inelast,m_SH,symm_class,orth_basis,orth_proj,NDI,NSHR,&
       &algo,env,hard_law,theta_KK,aver_adj_stp,conv_issue_dil)
  !
  ! determine polar decomposition of elastic deformation gradient
  Fe_n1=MATMUL(F,F_p_inv_n1)
  Re_n1=polar_decompostion_def_grad(Fe_n1,NDI,NSHR,env)
  !
  ! assign internal variable to history field
  CALL transform_tens_vec_form(F_p_inv_n1,SDV(1:9),ZERO_int,env)
  SDV(10:33)=plastic_slip_iter
  SDV(34)=SUM(plastic_slip_iter)
  SDV(35:37)=EulAngFromRotMat(Transpose(Re_n1))
  SDV(NSTATV)=aver_adj_stp
  !
  ! compute tangent
  SELECT CASE (numtang)
  CASE (ZERO_int) ! analytical tangent
     C_tang=consistent_algorithmic_tangent_augmented_lagrangian_penalty(&
          &C,F_p_inv_n,active_slip_system,d_delta_lambda_ip1_d_C,&
          &det_Fp_inv,F_p_inv_iter,Je_n1,Ce_n1,Ce_n1_inv,F_p_inv_n1,&
          &S2PK_IC,elast_const,m_SH,symm_class,orth_basis,orth_proj,NDI,&
          &NSHR,env)
  CASE (ONE_int)  ! CSDA tangent
     WRITE(16,'(A)')'CSDA tangent computation is not implemented in "crystal_plasticity_augm_lagr"'
     WRITE(16,'(A)')'declaration of complex variables required'
     WRITE(16,'(A)')'STOPPING COMPUTATION'
     STOP
  CASE (TWO_int)  ! FD tangent
     aver_adj_stp_save=aver_adj_stp
     !
     C_tang=CSDA_material_tangent_single_crystal_plast(S2PK,F,dt,&
          &F_p_inv_n,plastic_slip_n,elast_const,mat_param_inelast,m_SH,&
          &symm_class,orth_basis,orth_proj,numtang,algo,NDI,NSHR,env,&
          &hard_law,theta_KK,aver_adj_stp_save)
  CASE DEFAULT
     WRITE(16,'(A)')'Unknown numtang ID in "crystal_plasticity_augm_lagr"'
     WRITE(16,'(A)')'STOPPING COMPUTATION'
     STOP
  END SELECT
  !
END SUBROUTINE crystal_plasticity_augm_lagr
! ......................................................................
! **********************************************************************
!
! ********************************************************************
! determine the incremental lagrange multipliers based on augmented
! lagrangian formulation (fixed-point iteration with nested newton loop)
! proposed by Schmidt-Baldassari, CMAME (192), 2003
! ********************************************************************
SUBROUTINE determine_incremental_lagrange_multiplier_augmented_lagrange(&
     &F_p_inv_iter,plastic_slip_iter,d_delta_lambda_ip1_d_C,&
     &active_slip_system,unimod_C,F_p_inv_n,C_inv,Je,plastic_slip_n,&
     &algo,dt,elast_const,mat_param_inelast,m_SH,symm_class,orth_basis,&
     &orth_proj,NDI,NSHR,env,hard_law,theta_KK,conv_issue_max_ati)
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:SOLVE_SYSTEM_OF_EQUATION,&
       & back_tracking_step_length_scaling_parameter
  USE constitutive_routines_interface, only:comp_d_delta_lambda_ip1_d_C,&
       & residual_tangent_computation_nw_crystal_plasticity, &
       & residual_tangent_nw_crystal_plasticity_generic_solver, &
       & update_internal_variable_nw_iter
  !
  ! variable declaration
  ! F_p_inv_()..........inverse of plastic part of deformation gradient
  ! plastic_slip_().....accumulated plastic slip on each system
  ! d_delta_lambda_ip1_d_C..deriv. of incr. lagrange multiplier w.r.t.
  !                         right Cauchy-Green tensor
  ! active_slip_system..list of active slip systems
  ! unimod_C............right Cauchy-Green (unimodular in case of vol/
  !                     isochoric split elasticity)
  ! C_inv...............inverse of right Cauchy-Green tensor
  ! Je..................determinant of elastic part of deformation
  !                     gradient
  ! algo................id for algorithm employed
  !                     0...augmented-lagrangian
  !                     1...NCP-function (Kanzow-Kleinmichel)
  !                     2...Perzyna viscoplastic formulation
  !                     3...Ortiz/Miehe viscoplastic formulation
  ! dt..................time increment
  ! elast_const.........material parameter for elastic behavior
  ! mat_param_inelast...material parameter assoc. with inelastic behavior
  ! m_SH................Seth-Hill parameter
  ! symm_class..........symmetry included in elastic law
  ! orth_basis..........all second order, orthogonal basis N(i), rotated
  !                     in case of cubic elasticity
  ! orth_proj...........symmetric dyadic products of orthogonal basis
  !                     i.e. N(i) dyadic N(i)
  ! env.................environment from which material is called
  ! hard_law............id of hardening law
  !                     0....nonlin. Taylor-type (Schmidt-Baldassari)
  !                     1....nonlin. interaction matrix (Cailletaud,Forest)
  ! theta_KK............algorithmic parameter for Kanzow-Kleinmichel
  !                     NCP-function
  ! conv_issue_max_ati..error flag indicating non-convergence in
  !                     determination of incremental lagr. multiplier,
  !                     number of attempts exceeds maximum value
  !
  IMPLICIT NONE
  REAL(KIND=rk),DIMENSION(24),INTENT(OUT)::plastic_slip_iter
  REAL(KIND=rk),DIMENSION(3,3),INTENT(OUT)::F_p_inv_iter
  REAL(KIND=rk),DIMENSION(24,6),INTENT(OUT)::d_delta_lambda_ip1_d_C
  REAL(KIND=rk),INTENT(IN)::Je,dt,m_SH,theta_KK
  REAL(KIND=rk),DIMENSION(3),INTENT(IN)::elast_const
  REAL(KIND=rk),DIMENSION(13),INTENT(IN)::mat_param_inelast
  REAL(KIND=rk),DIMENSION(24),INTENT(IN)::plastic_slip_n
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::unimod_C,F_p_inv_n,C_inv
  REAL(KIND=rk),DIMENSION(6,3,3),INTENT(IN)::orth_basis
  REAL(KIND=rk),DIMENSION(6,6,6),INTENT(IN)::orth_proj
  REAL(KIND=rk)::eta,eta_save,mult_fact,lambda,lambda_two,obj_f_two
  REAL(KIND=rk)::obj_f,obj_f_save,tol_nw,tol_fp,mult_fact_init
  REAL(KIND=rk)::norm_diff_delta_lambda_i_ip1,aver_red_fp
  ! --- Phase 1 (Algorithm 2) additions ---
  REAL(KIND=rk),DIMENSION(24)::Phi_alpha_all,phi_tilde_plus
  REAL(KIND=rk)::norm_R,norm_Phi_tilde_inf,norm_Phi_tilde_prev,Phi_ref
  LOGICAL::kkt_ok
  ! ----------------------------------------
  REAL(KIND=rk),DIMENSION(24)::delta_lambda_ip1,delta_lambda_i
  REAL(KIND=rk),DIMENSION(24)::residual_nw,grad_f,grad_f_save
  REAL(KIND=rk),DIMENSION(24)::delta_lambda_ip1_save
  REAL(KIND=rk),DIMENSION(24)::deriv_visco_plast_potent
  REAL(KIND=rk),DIMENSION(24)::dd_delta_lambda_ip1,diff_delta_lambda
  REAL(KIND=rk),DIMENSION(24,1)::dummy_var
  REAL(KIND=rk),DIMENSION(24,6)::d_delta_lambda_i_d_C
  REAL(KIND=rk),DIMENSION(24,24)::jacobian_nw
  INTEGER(KIND=ik),INTENT(IN)::algo,symm_class,NDI,NSHR,hard_law
  INTEGER(KIND=ik)::l_iter_fp,l_attempt,l_iter_nw,back_track_count
  INTEGER(KIND=ik),PARAMETER::n_iter_max_nw=20,n_iter_max_fp=40
  INTEGER(KIND=ik),PARAMETER::n_attempt_max=20,n_max_bt=10,n_memory=3
  REAL(KIND=rk),DIMENSION(n_memory)::mem_mult_fact,mem_red_fp
  LOGICAL,DIMENSION(24),INTENT(OUT)::active_slip_system
  LOGICAL,INTENT(OUT)::conv_issue_max_ati
  LOGICAL::n_conv_fp,n_conv_nw,back_tracking,back_track_off,conv_issue
  !
  CHARACTER(4),INTENT(IN)::env
  !
  ! initialize variables for iterative procedure
  F_p_inv_iter=F_p_inv_n
  plastic_slip_iter=plastic_slip_n
  delta_lambda_ip1=ZERO_r
  delta_lambda_i=ZERO_r
  !
  ! variable for tangent computation
  d_delta_lambda_i_d_C=ZERO_r
  !
  ! parameter for augmented lagrangian approach
  ! Alg.2 line 1: Dt*eta(0) ~ 10/G  (paper-faithful, consistent MPa units)
  ! isotropic: G = elast_const(2) = mu ; cubic: G = elast_const(3) = C44
  IF (symm_class==ONE_int) THEN
     eta=(TWO_r*FIVE_r)/elast_const(TWO_int)
  ELSE
     eta=(TWO_r*FIVE_r)/elast_const(THREE_int)
  END IF
  eta_save=eta
  mult_fact=TWO_r
  mult_fact_init=mult_fact
  !
  ! parameters for newton and fixed-point iteration
  ! Step 0 / Alg.2 line 2: initial adaptive tolerance, Eq.(40)
  CALL residual_tangent_nw_crystal_plasticity_generic_solver(&
       &residual_nw,jacobian_nw,unimod_C,F_p_inv_n,Je,plastic_slip_n,&
       &delta_lambda_i,delta_lambda_ip1,algo,dt,m_SH,elast_const,&
       &symm_class,orth_basis,orth_proj,mat_param_inelast,NDI,NSHR,&
       &env,eta,hard_law,theta_KK)
  norm_R=SQRT(DOT_PRODUCT(residual_nw,residual_nw))
  tol_nw=eps_min_nw   ! C4: adaptive tolerance DISABLED, fixed tight tol
  tol_fp=zero_distinct
  ! Eq.(36) needs previous violation; initialise large so first step can grow eta
  norm_Phi_tilde_prev=HUGE(ONE_r)
  !
  ! initialization fixed-point iteration
  l_iter_fp=ZERO_int
  n_conv_fp=.TRUE.
  l_attempt=ZERO_int
  conv_issue_max_ati=.FALSE.
  !
  ! memory functions
  mem_mult_fact=mult_fact
  mem_red_fp=ZERO_r
  !
  ! fixed-point iteration (l_iter_fp)
  DO WHILE (n_conv_fp)
     l_iter_nw=ZERO_int
     n_conv_nw=.TRUE.
     ! initialization related to back-tracking
     lambda=ONE_r
     lambda_two=lambda
     back_tracking=.FALSE.
     back_track_count=ZERO_int
     obj_f_two=ZERO_r
     back_track_off=.FALSE.
     ! parameter for backtracking
     conv_issue=.FALSE.
     !
     ! nested newton iteration (l_iter_nw)
     DO WHILE (n_conv_nw)
        !
        ! evaluate residual and jacobian
        CALL residual_tangent_nw_crystal_plasticity_generic_solver(&
             &residual_nw,jacobian_nw,unimod_C,F_p_inv_n,Je,plastic_slip_n,&
             &delta_lambda_i,delta_lambda_ip1,algo,dt,m_SH,elast_const,&
             &symm_class,orth_basis,orth_proj,mat_param_inelast,NDI,NSHR,&
             &env,eta,hard_law,theta_KK)
        !
        ! compute quantities for backtracking
        obj_f=ONE_HALF_r*dot_product(residual_nw,residual_nw)
        grad_f=MATMUL(TRANSPOSE(jacobian_nw),residual_nw)
        !
        ! perform back-tracking if necessary
        IF ((l_iter_nw .GT. ZERO_int).AND.(.NOT.(back_track_off))) THEN
           CALL back_tracking_step_length_scaling_parameter(back_tracking,&
                &back_track_count,lambda,lambda_two,delta_lambda_ip1,&
                &delta_lambda_ip1_save,dd_delta_lambda_ip1,obj_f,&
                &obj_f_save,obj_f_two,grad_f_save,TWENTY_FOUR_int)
        END IF
        !
        ! safe guard to avoid infinite loop
        IF (back_track_count>n_max_bt) THEN
           n_conv_nw=.FALSE.
           conv_issue=.TRUE.
           l_attempt=l_attempt+ONE_int
        END IF
        !
        ! regular Newton step (switched off while back-tracking)
        IF (.NOT.(back_tracking)) THEN
           obj_f_save=obj_f
           grad_f_save=grad_f
           delta_lambda_ip1_save=delta_lambda_ip1
           !
           ! compute increment
           dummy_var(:,1)=residual_nw
           dummy_var=SOLVE_SYSTEM_OF_EQUATION(jacobian_nw,&
                &-dummy_var,ONE_int,TWENTY_FOUR_int)
           dd_delta_lambda_ip1=dummy_var(:,1)
           !
           ! Alg.2 line 4: residual-based exit with positivity guard
           ! exit when ||R|| <= eps^(i) AND all Dlambda >= 0
           delta_lambda_ip1=delta_lambda_ip1+dd_delta_lambda_ip1
           norm_R=SQRT(DOT_PRODUCT(residual_nw,residual_nw))
           ! --- DIAG: residual at every Newton iteration ---
           WRITE(16,'(A,1X,I4,1X,I4,1X,ES16.8,1X,ES16.8,1X,I3)') 'NEWTON',&
                &l_iter_fp,l_iter_nw,norm_R,tol_nw,COUNT(active_slip_system)
           ! --- END DIAG ---
           IF ((norm_R<tol_nw).AND.(ALL(delta_lambda_ip1>=-tol_kkt))) THEN
              n_conv_nw=.FALSE.
              l_attempt=ZERO_int
           ELSE
              l_iter_nw=l_iter_nw+ONE_int
           END IF
           !
           ! catch non-convergence or back-track failure
           IF ((l_iter_nw>n_iter_max_nw).OR.(conv_issue)) THEN
              n_conv_nw=.FALSE.
              l_attempt=l_attempt+ONE_int
           END IF
        END IF
        !
     END DO ! l_iter_nw
     !
     ! update internal variables and tangent and active slip system comp.
     CALL update_internal_variable_nw_iter(delta_lambda_ip1,&
          &plastic_slip_iter,F_p_inv_iter,dd_delta_lambda_ip1,&
          &plastic_slip_n,F_p_inv_n)
     CALL residual_tangent_computation_nw_crystal_plasticity(residual_nw,&
          &jacobian_nw,active_slip_system,deriv_visco_plast_potent,&
          &unimod_C,F_p_inv_iter,F_p_inv_n,Je,plastic_slip_iter,delta_lambda_i,&
          &delta_lambda_ip1,algo,dt,m_SH,elast_const,symm_class,orth_basis,&
          &orth_proj,mat_param_inelast,NDI,NSHR,env,eta,hard_law,theta_KK,Phi_alpha_all)
     !
     
     ! --- DIAG: residual per fixed-point iteration (FIXPT) ---
     WRITE(16,'(A,1X,I4,1X,ES16.8,1X,I3)') 'FIXPT',l_iter_fp,&
          &SQRT(DOT_PRODUCT(residual_nw,residual_nw)),COUNT(active_slip_system)
     ! --- END DIAG ---
     
     ! adjustment of pseudo-viscosity
     IF ((l_iter_nw>n_iter_max_nw).OR.(conv_issue)) THEN
        ! handling of unsuccessful back-tracking Newton iteration
        IF (l_attempt>n_attempt_max) THEN
           conv_issue_max_ati=.TRUE.
           n_conv_fp=.FALSE.
        END IF
        !
        n_conv_nw=.TRUE.
        l_iter_nw=ZERO_int
        ! eta=eta_save*(ONE_r+ONE_r/(TWO_r**(l_attempt)))
        eta=ZERO_p_SIX_r*eta_save
        eta_save=eta
        mem_mult_fact(TWO_int:n_memory)=mem_mult_fact(ONE_int:(n_memory-ONE_int))
        mem_mult_fact(ONE_int)=ZERO_p_SIX_r
        !
        delta_lambda_ip1=delta_lambda_i
     ELSE
        ! successful inner Newton iteration -> evaluate Alg.2 lines 8-13
        !
        ! Eq.(37): KKT violation measure Phi_tilde_+, guard Dt*eta>0
        phi_tilde_plus=MAX(Phi_alpha_all,&
             &-delta_lambda_ip1/MAX(dt*eta,lb_invert))
        norm_Phi_tilde_inf=MAXVAL(ABS(phi_tilde_plus))
        !
        ! reference stress (tau0) to non-dimensionalise stress-scale checks
        Phi_ref=MAX(ABS(mat_param_inelast(ONE_int)),lb_invert)
        !
        ! Alg.2 line 8: KKT check (stress conditions scaled by Phi_ref)
        kkt_ok=ALL(Phi_alpha_all<=tol_kkt*Phi_ref).AND.&
             &ALL(delta_lambda_ip1>=-tol_kkt).AND.&
             &ALL(ABS(Phi_alpha_all*delta_lambda_ip1)<=tol_kkt*Phi_ref)
        !
        IF (kkt_ok) THEN
           ! Alg.2 line 9: converged
           ! compute final consistent tangent at converged state
           d_delta_lambda_ip1_d_C=comp_d_delta_lambda_ip1_d_C(&
                &active_slip_system,jacobian_nw,d_delta_lambda_i_d_C,&
                &F_p_inv_iter,unimod_C,Je,C_inv,elast_const,m_SH,symm_class,&
                &orth_basis,orth_proj,eta,dt,algo,deriv_visco_plast_potent,&
                &env,NDI,NSHR)
           n_conv_fp=.FALSE.
           ! --- DIAG: converged increment (CONVRG) ---
           WRITE(16,'(A,1X,I4,1X,ES16.8,1X,I3)') 'CONVRG',l_iter_fp,&
                &SQRT(DOT_PRODUCT(residual_nw,residual_nw)),COUNT(active_slip_system)
           ! --- END DIAG ---
        ELSE
           ! Alg.2 line 11: implicit tangent update Eq.(21)
           d_delta_lambda_ip1_d_C=comp_d_delta_lambda_ip1_d_C(&
                &active_slip_system,jacobian_nw,d_delta_lambda_i_d_C,&
                &F_p_inv_iter,unimod_C,Je,C_inv,elast_const,m_SH,symm_class,&
                &orth_basis,orth_proj,eta,dt,algo,deriv_visco_plast_potent,&
                &env,NDI,NSHR)
           !
           ! Alg.2 line 12: adaptive Newton tolerance update Eq.(39)
           tol_nw=eps_min_nw   ! C4: adaptive update disabled (fixed tol)
           !
           ! Alg.2 line 13: penalty update Eq.(36) with beta=10, gamma=0.25
           eta_save=eta
           IF (norm_Phi_tilde_inf>gamma_eta*norm_Phi_tilde_prev) THEN
              eta=beta_eta*eta
           END IF
           norm_Phi_tilde_prev=norm_Phi_tilde_inf
           !
           ! advance fixed-point iterate i -> i+1
           delta_lambda_i=delta_lambda_ip1
           d_delta_lambda_i_d_C=d_delta_lambda_ip1_d_C
           l_iter_fp=l_iter_fp+ONE_int
        END IF
        !
        ! check maximum number of allowed fixed-point iterations
        IF (l_iter_fp>n_iter_max_fp) THEN
           conv_issue_max_ati=.TRUE.
           n_conv_fp=.FALSE.
        END IF
        !
     END IF
     !
  END DO ! l_iter_fp
  !
  ! variable assignment in case of divergence (number of fixed-point
  ! iterations or number of attempts exceeded)
  IF (conv_issue_max_ati) THEN
     F_p_inv_iter=F_p_inv_n
     plastic_slip_iter=plastic_slip_n
     d_delta_lambda_ip1_d_C=ZERO_r
     active_slip_system=.FALSE.
  END IF
  !
END SUBROUTINE determine_incremental_lagrange_multiplier_augmented_lagrange
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! determine the incremental lagrange multipliers based on rate-dependent
! viscoplastic formulation (newton-raphson scheme with cubic back-tracking
! and adaptively adjustment of viscoplastic parameters,
! type of rate-dependence proposed by Cuitino, Ortiz, Modelling Simul.
! Mater. Sci. Eng.,1 (1992), pp. 225-263 and Miehe, Schroeder, IJNME,
! 50 (2001), pp. 273-298
! **********************************************************************
SUBROUTINE determine_incremental_lagrange_multiplier_miehe_viscoplastic(&
     &F_p_inv_iter,plastic_slip_iter,d_delta_lambda_ip1_d_C,&
     &active_slip_system,unimod_C,F_p_inv_n,C_inv,Je,plastic_slip_n,&
     &algo,dt,elast_const,mat_param_inelast,m_SH,symm_class,orth_basis,&
     &orth_proj,NDI,NSHR,env,hard_law,theta_KK,aver_adjust_per_step,&
     &conv_issue)
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:compute_trace_second_order_tensor
  USE constitutive_routines_interface, only:comp_d_delta_lambda_ip1_d_C,&
       & extract_viscoplastic_parameters_crystal_plast,&
       & update_visco_in_inelast_param_crystal_plast, &
       & adjust_visco_parameter_miehe_vp_log,&
       & newton_back_track_solver_crystal_visco_plasticity, &
       & residual_tangent_computation_nw_crystal_plasticity, &
       & update_internal_variable_nw_iter
  !
  ! variable declaration
  ! F_p_inv_()..........inverse of plastic part of deformation gradient
  ! plastic_slip_().....accumulated plastic slip on each system
  ! d_delta_lambda_ip1_d_C..deriv. of incr. lagrange multiplier w.r.t.
  !                         right Cauchy-Green tensor
  ! active_slip_system..list of active slip systems
  ! unimod_C............right Cauchy-Green (unimodular in case of vol/
  !                     isochoric split elasticity)
  ! C_inv...............inverse of right Cauchy-Green tensor
  ! Je..................determinant of elastic part of deformation
  !                     gradient
  ! algo................id for algorithm employed
  !                     0...augmented-lagrangian
  !                     1...NCP-function (Kanzow-Kleinmichel)
  !                     2...Perzyna viscoplastic formulation
  !                     3...Ortiz/Miehe viscoplastic formulation
  ! dt..................time increment
  ! elast_const.........material parameter for elastic behavior
  ! mat_param_inelast...material parameter assoc. with inelastic behavior
  ! m_SH................Seth-Hill parameter
  ! symm_class..........symmetry included in elastic law
  ! orth_basis..........all second order, orthogonal basis N(i), rotated
  !                     in case of cubic elasticity
  ! orth_proj...........symmetric dyadic products of orthogonal basis
  !                     i.e. N(i) dyadic N(i)
  ! env.................environment from which material is called
  ! hard_law............id of hardening law
  !                     0....nonlin. Taylor-type (Schmidt-Baldassari)
  !                     1....nonlin. interaction matrix (Cailletaud,Forest)
  ! theta_KK............algorithmic parameter for Kanzow-Kleinmichel
  !                     NCP-function
  ! conv_issue..........flag indicating convergence issue in adaptively
  !                     ramping of viscoplastic material constants
  !
  IMPLICIT NONE
  REAL(KIND=rk),DIMENSION(24),INTENT(OUT)::plastic_slip_iter
  REAL(KIND=rk),DIMENSION(3,3),INTENT(OUT)::F_p_inv_iter
  REAL(KIND=rk),DIMENSION(24,6),INTENT(OUT)::d_delta_lambda_ip1_d_C
  REAL(KIND=rk),INTENT(IN)::Je,dt,m_SH,theta_KK
  REAL(KIND=rk),INTENT(INOUT)::aver_adjust_per_step
  REAL(KIND=rk),DIMENSION(3),INTENT(IN)::elast_const
  REAL(KIND=rk),DIMENSION(13),INTENT(IN)::mat_param_inelast
  REAL(KIND=rk),DIMENSION(24),INTENT(IN)::plastic_slip_n
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::unimod_C,F_p_inv_n,C_inv
  REAL(KIND=rk),DIMENSION(6,3,3),INTENT(IN)::orth_basis
  REAL(KIND=rk),DIMENSION(6,6,6),INTENT(IN)::orth_proj
  REAL(KIND=rk)::eta,tol_nw,p_exp_vp_target,p_exp_vp_t,tau_relax_vp_t
  REAL(KIND=rk)::tau_relax_vp_target
  REAL(KIND=rk),DIMENSION(3)::mat_param_visco
  REAL(KIND=rk),DIMENSION(13)::mat_param_inelast_tp
  REAL(KIND=rk),DIMENSION(24)::delta_lambda_ip1,delta_lambda_i
  REAL(KIND=rk),DIMENSION(24)::delta_lambda_ip1_save,dummy_var
  REAL(KIND=rk),DIMENSION(24)::deriv_visco_plast_potent,Phi_alpha_dum
  REAL(KIND=rk),DIMENSION(24,6)::d_delta_lambda_i_d_C
  REAL(KIND=rk),DIMENSION(24,24)::jacobian_nw
  INTEGER(KIND=ik),INTENT(IN)::algo,symm_class,NDI,NSHR,hard_law
  INTEGER(KIND=ik)::adp_i,aver_adj_ps,i_dil,n_adjust
  INTEGER(KIND=ik),PARAMETER::n_iter_max_nw=20,n_max_bt=10,adp_i_max=200
  LOGICAL,INTENT(OUT)::conv_issue
  LOGICAL,DIMENSION(24),INTENT(OUT)::active_slip_system
  LOGICAL::adp_n_conv,conv_issue_nbt
  CHARACTER(4),INTENT(IN)::env
  !
  ! initialize variables for iterative procedure
  F_p_inv_iter=F_p_inv_n
  plastic_slip_iter=plastic_slip_n
  delta_lambda_ip1=ZERO_r
  delta_lambda_i=ZERO_r
  !
  ! variable for tangent computation
  d_delta_lambda_i_d_C=ZERO_r
  !
  ! parameter for augmented lagrangian approach
  eta=ONE_HALF_r
  !
  ! parameters for newton iteration and adaptive parameter adjustment
  tol_nw=tol_nw_red
  adp_i=ZERO_int; n_adjust=ZERO_int
  adp_n_conv=.TRUE.
  conv_issue=.FALSE.
  !
  delta_lambda_ip1_save=delta_lambda_ip1
  !
  ! extract viscoplastic material parameters and adjust their values
  ! according to previous convergence history
  mat_param_visco=extract_viscoplastic_parameters_crystal_plast(&
       &mat_param_inelast)
  p_exp_vp_target=mat_param_visco(1); tau_relax_vp_target=mat_param_visco(2)
  p_exp_vp_t=p_exp_vp_target; tau_relax_vp_t=tau_relax_vp_target;
  aver_adj_ps=MAX(ZERO_int,FLOOR((aver_adjust_per_step-ONE_r),ik))
  DO i_dil=ONE_int,aver_adj_ps
     CALL adjust_visco_parameter_miehe_vp_log (p_exp_vp_t,tau_relax_vp_t,&
          & 'decr_exp_incr_tau',p_exp_vp_target,tau_relax_vp_target)
  END DO
  mat_param_visco(1)=p_exp_vp_t; mat_param_visco(2)=tau_relax_vp_t
  mat_param_inelast_tp=update_visco_in_inelast_param_crystal_plast(&
       &mat_param_visco,mat_param_inelast)
  !
  ! iterative loop to adaptively ramp the viscoplastic parameters
  DO WHILE (adp_n_conv)
     adp_i=adp_i+ONE_int
     delta_lambda_ip1=delta_lambda_ip1_save
     CALL newton_back_track_solver_crystal_visco_plasticity(&
          &delta_lambda_ip1,conv_issue_nbt,unimod_C,F_p_inv_n,Je,&
          &plastic_slip_n,delta_lambda_i,algo,dt,elast_const,&
          &mat_param_inelast_tp,m_SH,symm_class,orth_basis,orth_proj,&
          &NDI,NSHR,env,hard_law,theta_KK,eta,tol_nw,n_iter_max_nw,&
          &n_max_bt)
     !
     IF (.NOT.(conv_issue_nbt)) THEN
        ! check if visco parameters ramped to target values
        IF((ABS(p_exp_vp_t-p_exp_vp_target)<zero_distinct).AND.&
             &(ABS(tau_relax_vp_t-tau_relax_vp_target)<zero_distinct)&
             &)THEN
           adp_n_conv=.FALSE.
        ELSE
           CALL adjust_visco_parameter_miehe_vp_log (p_exp_vp_t,&
                &tau_relax_vp_t,'incr_exp_decr_tau',p_exp_vp_target,&
                &tau_relax_vp_target)
           mat_param_visco(1)=p_exp_vp_t
           mat_param_visco(2)=tau_relax_vp_t
           mat_param_inelast_tp=update_visco_in_inelast_param_crystal_plast(&
                &mat_param_visco,mat_param_inelast)
           delta_lambda_ip1_save=delta_lambda_ip1
        END IF
     ELSE
        CALL adjust_visco_parameter_miehe_vp_log (p_exp_vp_t,&
             &tau_relax_vp_t,'decr_exp_incr_tau',p_exp_vp_target,&
             &tau_relax_vp_target)
        mat_param_visco(1)=p_exp_vp_t; mat_param_visco(2)=tau_relax_vp_t
        mat_param_inelast_tp=update_visco_in_inelast_param_crystal_plast(&
             &mat_param_visco,mat_param_inelast)
        n_adjust=n_adjust+ONE_int
     END IF
     !
     IF (adp_i>adp_i_max) THEN
        adp_n_conv=.FALSE.
        conv_issue=.TRUE.
     END IF
  END DO ! adp_n_conv
  !
  ! update number of adaptive steps and internal variables / 
  ! compute tangent and active set
  IF (.NOT.(conv_issue)) THEN
     aver_adjust_per_step=aver_adjust_per_step+REAL(n_adjust)
     dummy_var=ZERO_r
     CALL update_internal_variable_nw_iter(delta_lambda_ip1,&
          &plastic_slip_iter,F_p_inv_iter,dummy_var,plastic_slip_n,&
          &F_p_inv_n)
     CALL residual_tangent_computation_nw_crystal_plasticity(dummy_var,&
          &jacobian_nw,active_slip_system,deriv_visco_plast_potent,&
          &unimod_C,F_p_inv_iter,F_p_inv_n,Je,plastic_slip_iter,&
          &delta_lambda_i,delta_lambda_ip1,algo,dt,m_SH,elast_const,&
          &symm_class,orth_basis,orth_proj,mat_param_inelast,NDI,NSHR,&
          &env,eta,hard_law,theta_KK,Phi_alpha_dum)
  ELSE
     delta_lambda_ip1=ZERO_r
     FORALL (i_dil=1:24) jacobian_nw(i_dil, i_dil)=tau_relax_vp_target/dt
     active_slip_system=.FALSE.
     deriv_visco_plast_potent=ZERO_r
  END IF
  !
  ! compute the tangent contribution (derivative of incremental lagr.
  ! multiplier / hardening variables w.r.t. right Cauchy-Green tensor
  d_delta_lambda_ip1_d_C=comp_d_delta_lambda_ip1_d_C(active_slip_system,&
       &jacobian_nw,d_delta_lambda_i_d_C,F_p_inv_iter,unimod_C,Je,C_inv,&
       &elast_const,m_SH,symm_class,orth_basis,orth_proj,eta,dt,algo,&
       &deriv_visco_plast_potent,env,NDI,NSHR)
  !
END SUBROUTINE determine_incremental_lagrange_multiplier_miehe_viscoplastic
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! Compute the diagonal functions, their derivatives and projectors
! **********************************************************************
SUBROUTINE diagonal_funct_deriv_projectors(e_vect,d_vect,f_vect,Proj_tens,&
     &lambda,index_distinct,n_distinct,Ce,m_SH,tol)
  !
  ! based on: Miehe, Lambrecht, Algorithms of computation of stress and
  !             elasticity moduli in terms of Seth-Hill family of
  !             generalized strain tensors, CINME, 2001
  !
  USE Abaqus_Interface
  USE Constants
  USE constitutive_routines_interface, only:aux_prod_projector
  !
  ! declaration of variables
  !
  ! lambda...........vector of eigenvalues
  ! tol..............tolerance to equality of eigenvalues
  ! index_distinct...index of the distinct eigenvalues in lambda
  ! n_distinct.......number of distinct eigenvalues
  ! Ce...............elastic right Cauchy-Green tensor
  ! m_SH.............Seth-Hill parameter
  ! ()_vect..........diagonal functions and their 1st and 2nd derivative
  ! Proj_tens........projectors to transform stress and tangent
  !                  (stress dual to Seth-Hill strain to 2nd Piola
  !                  Kirchhoff and material tangent
  !
  IMPLICIT NONE
  !
  INTEGER(KIND=ik),INTENT(IN)::n_distinct
  INTEGER(KIND=ik),DIMENSION(3),INTENT(IN)::index_distinct
  INTEGER(KIND=ik)::i_count,ind_t,count,j_count
  !
  REAL(KIND=rk),DIMENSION(3),INTENT(IN)::lambda
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::Ce
  REAL(KIND=rk),INTENT(IN)::tol,m_SH
  REAL(KIND=rk),DIMENSION(3),INTENT(OUT)::e_vect,d_vect,f_vect
  REAL(KIND=rk),DIMENSION(3,3,3),INTENT(OUT)::Proj_tens
  REAL(KIND=rk),DIMENSION(3,3)::Proj_t,temp_var
  !
  e_vect=ZERO_r; d_vect=ZERO_r; f_vect=ZERO_r; Proj_tens=ZERO_r
  !
  DO i_count=ONE_int,n_distinct       
     ind_t=index_distinct(i_count)
     ! compute diagonal functions and the 1st and 2nd derivative
     IF (ABS(m_SH) > tol) THEN
        e_vect(i_count)=ONE_r/m_SH*(lambda(ind_t)**(m_SH/TWO_r)-ONE_r)
        d_vect(i_count)=lambda(ind_t)**(m_SH/TWO_r-ONE_r)
        f_vect(i_count)=TWO_r*(m_SH/TWO_r-ONE_r)*lambda(ind_t)**&
             &(m_SH/TWO_r-TWO_r)
     ELSE
        e_vect(i_count)=ONE_HALF_r*log(lambda(ind_t))
        d_vect(i_count)=ONE_r/lambda(ind_t)
        f_vect(i_count)=-TWO_r/lambda(ind_t)**TWO_int;   
     END IF
     ! compute projectors
     count=ONE_int
     Proj_t=ZERO_r
     DO j_count=1,n_distinct
        IF (j_count /= i_count) THEN
           IF (count==ONE_int) THEN
              Proj_t = aux_prod_projector(Ce,lambda,index_distinct,&
                   &i_count,j_count)

           ELSE
              temp_var=aux_prod_projector(Ce,lambda,&
                   &index_distinct,i_count,j_count)
              Proj_t = MATMUL(Proj_t,temp_var)

           END IF
           count=count+ONE_int
        END IF
     END DO ! j_count
     IF (n_distinct==ONE_int) THEN
        Proj_t=Proj_t+IDENT
     END IF

     Proj_tens(i_count,:,:)=Proj_t
  END DO ! i_count
  !
END SUBROUTINE diagonal_funct_deriv_projectors
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! evaluate the yield function of the current slip system and its 
! derivative w.r.t. increments in lagrange multipliers
! **********************************************************************
SUBROUTINE eval_yield_function_linearization(Phi_alpha,&
     &dPhi_alpha_d_delta_lambda_ip1_beta,tau_rss,&
     &d_tau_alpha_d_delta_lambda_ip1_beta,plastic_slip,i_slip,algo,&
     &mat_param_inelast,hard_law)
  !
  USE Abaqus_Interface
  USE Constants
  USE constitutive_routines_interface, only:extract_hardening_parameters_crystal_plast,&
       & hardening_function_cailletaud_forest, &
       & hardening_function_schmidt_baldassari
  !
  ! variable declaration
  ! Phi_alpha..........value of yield function of current slip system
  ! dPhi_alpha_d_delta_lambda_ip1_beta...deriv. yield function w.r.t.
  !                    all incremental lagrange multiplier beta
  ! tau_rss............resolved shear stress on current slip system
  ! d_tau_alpha_d_delta_lambda_ip1_beta..deriv. of resolved shear stress
  !                    w.r.t. incremental lagrange multiplier
  ! plastic_slip.......accumulated plastic slip on each slip system
  ! i_slip.............number of current slip system
  ! algo...............id for algorithm employed
  ! mat_param_inelast..material parameter assoc. with inelastic behavior
  ! hard_law...........id of hardening law
  !                    0....nonlin. Taylor-type (Schmidt-Baldassari)
  !                    1....nonlin. interaction matrix (Cailletaud,Forest)
  IMPLICIT NONE
  REAL(KIND=rk),INTENT(OUT)::Phi_alpha
  REAL(KIND=rk),DIMENSION(24),INTENT(OUT)::dPhi_alpha_d_delta_lambda_ip1_beta
  REAL(KIND=rk),INTENT(IN)::tau_rss
  REAL(KIND=rk),DIMENSION(13),INTENT(IN)::mat_param_inelast
  REAL(KIND=rk),DIMENSION(24),INTENT(IN)::d_tau_alpha_d_delta_lambda_ip1_beta
  REAL(KIND=rk),DIMENSION(24),INTENT(IN)::plastic_slip
  REAL(KIND=rk)::Y_alpha,Y
  REAL(KIND=rk),DIMENSION(3)::hard_param_sb
  REAL(KIND=rk),DIMENSION(10)::hard_param_cf
  REAL(KIND=rk),DIMENSION(24)::dY_alpha_d_delta_lambda_beta
  INTEGER(KIND=ik),INTENT(IN)::i_slip,algo,hard_law
  !
  ! current hardening (change in yield stress) and its derivative w.r.t.
  ! increments in lagrange multipliers
  SELECT CASE (hard_law)
  CASE (ZERO_int) ! nonlinear Taylor-type hardening
     hard_param_sb=extract_hardening_parameters_crystal_plast(&
          &mat_param_inelast,hard_law,THREE_int)
     CALL hardening_function_schmidt_baldassari(Y_alpha,&
          &dY_alpha_d_delta_lambda_beta,plastic_slip,hard_param_sb)
     Y=hard_param_sb(1)
     !
  CASE (ONE_int)  ! nonlinear hardening based on interaction matrix
     hard_param_cf=extract_hardening_parameters_crystal_plast(&
          &mat_param_inelast,hard_law,TEN_int)
     CALL hardening_function_cailletaud_forest(Y_alpha,&
          &dY_alpha_d_delta_lambda_beta,plastic_slip,i_slip,&
          &hard_param_cf)
     Y=hard_param_cf(1)
     !
  CASE DEFAULT
     WRITE(16,'(A)')'Unknown value of hard_law in: "eval_yield_function_linearization"'
     WRITE(16,'(A)')'STOPPING COMPUTATION'
     STOP
  END SELECT
  !
  ! yield function and corresponding derivatives
  Phi_alpha=tau_rss-(Y+Y_alpha)
  dPhi_alpha_d_delta_lambda_ip1_beta=d_tau_alpha_d_delta_lambda_ip1_beta-&
       &dY_alpha_d_delta_lambda_beta
  !
END SUBROUTINE eval_yield_function_linearization
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! hardening function for single crystal plasticity, incorporationg
! interaction matrix and exponential hardening
! based on: Busso, Cailletaud, IJP, 21 (2005), p. 2212-2231
!           Forest, Rubin, EJMA/Solids, 55 (2016), p. 278-288
! **********************************************************************
SUBROUTINE hardening_function_cailletaud_forest(Y_alpha,&
     &dY_alpha_d_delta_lambda_beta,plastic_slip,alpha_slip,&
     &hard_param_cf)
  !
  USE Abaqus_Interface
  USE Constants
  USE constitutive_routines_interface, only:interaction_matrix
  !
  ! variable declaration
  ! Y_alpha...........hardening on slip system alpha
  ! dY_alpha_d_delta_lambda_beta....deriv. of hardening w.r.t. increments
  !                   in lagrange multiplier
  ! plastic_slip......accumulated slip on each slip system
  ! hard_param_sb.....hardening parameters
  ! alpha_slip........current slip system
  !
  IMPLICIT NONE
  REAL(KIND=rk),INTENT(OUT)::Y_alpha
  REAL(KIND=rk),DIMENSION(24),INTENT(OUT)::dY_alpha_d_delta_lambda_beta
  REAL(KIND=rk),DIMENSION(10),INTENT(IN)::hard_param_cf
  REAL(KIND=rk),DIMENSION(24),INTENT(IN)::plastic_slip
  REAL(KIND=rk)::h,DeltaY,h_lin
  REAL(KIND=rk),DIMENSION(24)::aux_funct_vec
  REAL(KIND=rk),DIMENSION(24,24)::interact_matrix
  INTEGER(KIND=ik),INTENT(IN)::alpha_slip
  INTEGER(KIND=ik)::i_hfc
  !
  ! assign material parameters
  interact_matrix=interaction_matrix(hard_param_cf(4:9))
  DeltaY=hard_param_cf(2)
  h=hard_param_cf(3)
  h_lin=hard_param_cf(10)
  !
  DO i_hfc=ONE_int,TWENTY_FOUR_int
     aux_funct_vec(i_hfc)=EXP(-h*plastic_slip(i_hfc))
  END DO
  !
  Y_alpha=DOT_PRODUCT(interact_matrix(alpha_slip,:),(DeltaY*(ONE_r-&
       &aux_funct_vec)+h_lin*plastic_slip))
  dY_alpha_d_delta_lambda_beta=interact_matrix(alpha_slip,:)*&
       &(DeltaY*h*aux_funct_vec+h_lin)
  !
END SUBROUTINE hardening_function_cailletaud_forest
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! computes the hardening (increase in yield stress) on all slip systems
! based on the current estimate of plastic slip and derivative of
! hardening w.r.t. increment in lagrange multipliers
! Schmidt-Baldassari, CMAME 192 (2003), pp.1261-1280
! **********************************************************************
SUBROUTINE hardening_function_schmidt_baldassari(Y_alpha,&
     &dY_alpha_d_delta_lambda_beta,plastic_slip,hard_param_sb)
  !
  USE Abaqus_Interface
  USE Constants
  !
  ! variable declaration
  ! Y_alpha...........hardening on slip system alpha
  ! dY_alpha_d_delta_lambda_beta....deriv. of hardening w.r.t. increments
  !                   in lagrange multiplier
  ! plastic_slip......accumulated slip on each slip system
  ! hard_param_sb.....hardening parameters
  !
  IMPLICIT NONE
  REAL(KIND=rk),INTENT(OUT)::Y_alpha
  REAL(KIND=rk),DIMENSION(24),INTENT(OUT)::dY_alpha_d_delta_lambda_beta
  REAL(KIND=rk),DIMENSION(3),INTENT(IN)::hard_param_sb
  REAL(KIND=rk),DIMENSION(24),INTENT(IN)::plastic_slip
  REAL(KIND=rk)::aux_var,hardening_var_A,h,DeltaY
  !
  ! assign material parameters
  DeltaY=hard_param_sb(2)
  h=hard_param_sb(3)
  !
  ! hardening variable and auxilliary function
  hardening_var_A=SUM(plastic_slip)
  aux_var=EXP(-h*hardening_var_A)
  !
  ! hardening and corresponding derivative
  Y_alpha=DeltaY*(ONE_r-aux_var)
  dY_alpha_d_delta_lambda_beta=DeltaY*h*aux_var
  !
END SUBROUTINE hardening_function_schmidt_baldassari
! ......................................................................
! **********************************************************************
!
! ********************************************************************
! compute Mandel stress tensor (hydrostatic and deviatoric part) model
! employs decoupled volumetric / isochoric split of an isotropic
! neo-hookean type free energy function
! ********************************************************************
SUBROUTINE isotr_elasticity_mandel_stress_volumetric_isochoric_decoupled(&
     &M_hyd,M_dev,Je,unimod_Ce,kappa,mu)
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:compute_trace_second_order_tensor
  !
  ! variable declaration
  ! M_hyd..........hydrostatic Mandel stress 1/3 tr(M)
  ! M_dev..........deviatoric Mandel stress tensor
  ! Je.............determinant of elastic part of deformation gradient
  ! unimod_Ce......unimodular part of elastic Cauchy-Green tensor
  ! kappa, mu......elastic constants (bulk-, shear modulus)
  !
  IMPLICIT NONE
  REAL(KIND=rk),INTENT(OUT)::M_hyd
  REAL(KIND=rk),DIMENSION(3,3),INTENT(OUT)::M_dev
  REAL(KIND=rk),INTENT(IN)::Je,kappa,mu
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::unimod_Ce
  !
  M_hyd=kappa*LOG(Je)
  M_dev=mu*unimod_Ce-ONE_by_THREE_r*mu*&
       &compute_trace_second_order_tensor(unimod_Ce)*IDENT
  !
END SUBROUTINE isotr_elasticity_mandel_stress_volumetric_isochoric_decoupled
! ....................................................................
! ********************************************************************
!
! ********************************************************************
! linear elasticity (isotropic/cubic) based on subspace projection
! ********************************************************************
SUBROUTINE linear_elasticity_subspace_projection(stress,tangent,strain,&
     &elast_const,symm_class,orth_basis,orth_proj)
  !
  ! linear elasticity for isotropy and cubic anisotropy based on
  ! subspace projections
  ! Mahnken, Anisotropy in geometrically non-linear elasticity with
  ! generalized Seth-Hill strain tensors projected to invariant
  ! subspaces, CNME, 21, p.405, 2005
  ! in contrast to the original paper the isotropic elasticity is
  ! defined in terms of kappa and mu instead of lambda and mu
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:double_contraction_second_order_tensors
  !
  ! declaration of variables
  !
  ! strain........Seth-Hill strain tensor
  ! elast_const...elastic const. (isotr.: kappa, mu; cubic:C11,C12,C44)
  ! symm_class....symmetry class of crystal (1=isotropic, 2=cubic)
  ! orth_basis....all second order, orthogonal basis N(i), rotated
  !               in case of cubic elasticity
  ! orth_proj.....symmetric dyadic products of orthogonal basis
  !               i.e. N(i) dyadic N(i)
  ! stress........stress tensor conjugate to Seth-Hill strain
  ! tangent.......tangent consistent with above stress tensor
  !
  IMPLICIT NONE
  !
  INTEGER(KIND=ik),INTENT(IN)::symm_class
  INTEGER(KIND=ik)::n_proj,i_count,count,j_count
  INTEGER(KIND=ik),DIMENSION(3)::n_basis_list
  !
  REAL(KIND=rk),DIMENSION(3),INTENT(IN)::elast_const
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::strain
  REAL(KIND=rk),DIMENSION(6,3,3),INTENT(IN)::orth_basis
  REAL(KIND=rk),DIMENSION(6,6,6),INTENT(IN)::orth_proj
  !
  REAL(KIND=rk),DIMENSION(3,3),INTENT(OUT)::stress
  REAL(KIND=rk),DIMENSION(6,6),INTENT(OUT)::tangent
  REAL(KIND=rk),DIMENSION(3)::alpha_vec
  REAL(KIND=rk),DIMENSION(3,3)::orth_basis_t,proj_gen_strain_t
  REAL(KIND=rk),DIMENSION(6,6)::orth_proj_t
  !
  ! set number of projectors and basis for projector based on symmetry
  ! class
  SELECT CASE (symm_class)
  CASE (ONE_int) ! isotropic elasticity based on Seth-Hill strain
     n_proj=TWO_int
     n_basis_list=(/ONE_int,FIVE_int,FIVE_int/)
     alpha_vec=(/THREE_r*elast_const(1),TWO_r*elast_const(2),TWO_r*&
          &elast_const(2)/)
  CASE (TWO_int) ! cubic elasticity based on Seth-Hill strain
     n_proj=THREE_int
     n_basis_list=(/ONE_int,TWO_int,THREE_int/)
     alpha_vec=(/elast_const(1)+TWO_r*elast_const(2),elast_const(1)-&
          &elast_const(2),TWO_r*elast_const(3)/)
  CASE DEFAULT
     WRITE(16,'(A)')'Unknown symmetry class in "linear_elasticity_subspace_projection"'
     WRITE(16,'(A)')'STOPPING COMPUTATION'
     STOP
  END SELECT
  !
  ! compute stress and tangent based on strain, material parameters
  ! and projectors
  count=ZERO_int
  tangent=ZERO_r; stress=ZERO_r
  DO i_count=ONE_int,n_proj
     count=count+ONE_int
     orth_basis_t=orth_basis(count,:,:)
     orth_proj_t=orth_proj(count,:,:)
     proj_gen_strain_t=double_contraction_second_order_tensors(orth_basis_t,&
          &strain)*orth_basis_t
     !
     DO j_count=TWO_int,n_basis_list(i_count)
        count=count+ONE_int
        orth_basis_t=orth_basis(count,:,:)
        orth_proj_t=orth_proj_t+orth_proj(count,:,:)
        proj_gen_strain_t=proj_gen_strain_t + &
             & double_contraction_second_order_tensors(orth_basis_t,&
             &strain)*orth_basis_t
     END DO
     !
     tangent=tangent+orth_proj_t*alpha_vec(i_count)
     stress=stress+proj_gen_strain_t*alpha_vec(i_count)
     !
  END DO
END SUBROUTINE linear_elasticity_subspace_projection
! ....................................................................
! ********************************************************************
!
! ********************************************************************
! newton type solver with cubic back-tracking applied to crystal (visco)
! plasticity
! **********************************************************************
SUBROUTINE newton_back_track_solver_crystal_visco_plasticity(&
     &delta_lambda_ip1,conv_issue,unimod_C,F_p_inv_n,Je,plastic_slip_n,&
     &delta_lambda_i,algo,dt,elast_const,mat_param_inelast,m_SH,&
     &symm_class,orth_basis,orth_proj,NDI,NSHR,env,hard_law,theta_KK,&
     &eta,tol_nw,n_iter_max_nw,n_max_bt)
  !
  ! Based on: Press et al., Numerical recipes in Fortran, Cambridge
  ! University Press, 1996
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:SOLVE_ILL_COND_SYSTEM_OF_EQUATION,&
       & back_tracking_step_length_scaling_parameter
  USE constitutive_routines_interface, only:residual_tangent_nw_crystal_plasticity_generic_solver
  !
  ! declaration of variables
  ! delta_lambda_i/ip1..vector containing incremental lagrange
  !                     multipliers/ incremental hardening variables
  ! conv_issue..........error flag to handle non-convergence (number
  !                     nw iteration and/or back-tracking iteration
  !                     exceeds corresponding upper bound)
  ! unimod_C............right Cauchy-Green (unimodular in case of vol/
  !                     isochoric split elasticity)
  ! F_p_inv_()..........inverse of plastic part of deformation gradient
  ! Je..................determinant of elastic part of deformation
  !                     gradient
  ! plastic_slip_().....accumulated plastic slip on each system
  ! algo................id for algorithm employed
  !                     0...augmented-lagrangian
  !                     1...NCP-function (Kanzow-Kleinmichel)
  !                     2...Perzyna viscoplastic formulation
  !                     3...Ortiz/Miehe viscoplastic formulation
  ! dt..................time increment
  ! elast_const.........material parameter for elastic behavior
  ! mat_param_inelast...material parameter assoc. with inelastic behavior
  ! m_SH................Seth-Hill parameter
  ! symm_class..........symmetry included in elastic law
  ! orth_basis..........all second order, orthogonal basis N(i), rotated
  !                     in case of cubic elasticity
  ! orth_proj...........symmetric dyadic products of orthogonal basis
  !                     i.e. N(i) dyadic N(i)
  ! env.................environment from which material is called
  ! hard_law............id of hardening law
  !                     0....nonlin. Taylor-type (Schmidt-Baldassari)
  !                     1....nonlin. interaction matrix (Cailletaud,Forest)
  ! theta_KK............algorithmic parameter for Kanzow-Kleinmichel
  !                     NCP-function
  ! eta.................pseudo-viscosity (employed in augmented-lagr.)
  ! tol_nw..............tolerance for iterative improvement in newton-
  !                     scheme
  ! n_iter_max_nw.......maximum number of newton iterations,
  ! n_max_bt............maximum number of back-tracking iterations
  !                          
  IMPLICIT NONE
  REAL(KIND=rk),INTENT(IN)::Je,dt,m_SH,theta_KK,eta,tol_nw
  REAL(KIND=rk),DIMENSION(3),INTENT(IN)::elast_const
  REAL(KIND=rk),DIMENSION(13),INTENT(IN)::mat_param_inelast
  REAL(KIND=rk),DIMENSION(24),INTENT(IN)::delta_lambda_i
  REAL(KIND=rk),DIMENSION(24),INTENT(IN)::plastic_slip_n
  REAL(KIND=rk),DIMENSION(24),INTENT(INOUT)::delta_lambda_ip1
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::unimod_C,F_p_inv_n
  REAL(KIND=rk),DIMENSION(6,3,3),INTENT(IN)::orth_basis
  REAL(KIND=rk),DIMENSION(6,6,6),INTENT(IN)::orth_proj
  REAL(KIND=rk)::lambda,lambda_two,obj_f_two,obj_f,obj_f_save
  REAL(KIND=rk),DIMENSION(24)::residual_nw,grad_f,grad_f_save
  REAL(KIND=rk),DIMENSION(24)::delta_lambda_ip1_save_bt
  REAL(KIND=rk),DIMENSION(24)::dd_delta_lambda_ip1
  REAL(KIND=rk),DIMENSION(24,1)::dummy_var
  REAL(KIND=rk),DIMENSION(24,24)::jacobian_nw
  INTEGER(KIND=ik),INTENT(IN)::algo,symm_class,NDI,NSHR,hard_law
  INTEGER(KIND=ik),INTENT(IN)::n_iter_max_nw,n_max_bt
  INTEGER(KIND=ik)::l_iter_nw,back_track_count
  LOGICAL::n_conv_nw,back_tracking,back_track_off
  LOGICAL,INTENT(OUT)::conv_issue
  CHARACTER(4),INTENT(IN)::env
  !
  ! initialization newton iteration
  l_iter_nw=ZERO_int
  n_conv_nw=.TRUE.
  !
  ! initialization related to back-tracking
  lambda=ONE_r
  lambda_two=lambda
  back_tracking=.FALSE.
  back_track_count=ZERO_int
  obj_f_two=ZERO_r
  back_track_off=.FALSE.
  ! parameter for backtracking
  conv_issue=.FALSE.
  !
  ! newton iteration (l_iter_nw)
  DO WHILE (n_conv_nw)
     !
     ! evaluate residual and jacobian
     CALL residual_tangent_nw_crystal_plasticity_generic_solver(&
          &residual_nw,jacobian_nw,unimod_C,F_p_inv_n,Je,plastic_slip_n,&
          &delta_lambda_i,delta_lambda_ip1,algo,dt,m_SH,elast_const,&
          &symm_class,orth_basis,orth_proj,mat_param_inelast,NDI,NSHR,&
          &env,eta,hard_law,theta_KK)
     !
     ! compute quantities for backtracking
     obj_f=ONE_HALF_r*dot_product(residual_nw,residual_nw)
     grad_f=MATMUL(TRANSPOSE(jacobian_nw),residual_nw)
     !
     ! perform back-tracking if necessary
     IF ((l_iter_nw .GT. ZERO_int).AND.(.NOT.(back_track_off))) THEN
        CALL back_tracking_step_length_scaling_parameter(back_tracking,&
             &back_track_count,lambda,lambda_two,delta_lambda_ip1,&
             &delta_lambda_ip1_save_bt,dd_delta_lambda_ip1,obj_f,&
             &obj_f_save,obj_f_two,grad_f_save,TWENTY_FOUR_int)
     END IF
     !
     ! safe guard to avoid infinite loop
     IF (back_track_count>n_max_bt) THEN
        n_conv_nw=.FALSE.
        conv_issue=.TRUE.
     END IF
     !
     ! regular Newton step (switched off while back-tracking)
     IF ((.NOT.(back_tracking)).AND.(.NOT.(conv_issue))) THEN
        obj_f_save=obj_f
        grad_f_save=grad_f
        delta_lambda_ip1_save_bt=delta_lambda_ip1
        !
        ! compute iterative improvement
        dummy_var(:,1)=residual_nw
        dummy_var=SOLVE_ILL_COND_SYSTEM_OF_EQUATION(jacobian_nw,&
             &-dummy_var,ONE_int,TWENTY_FOUR_int,-eps_iter)
        dd_delta_lambda_ip1=dummy_var(:,1)
        !
        ! check convergence
        IF ((SQRT(DOT_PRODUCT(dd_delta_lambda_ip1,dd_delta_lambda_ip1))&
             &<tol_nw).OR.(SQRT(DOT_PRODUCT(residual_nw,residual_nw))&
             &<tol_nw_red))THEN
           n_conv_nw=.FALSE.
           !
        ELSE
           ! update and increase counter
           delta_lambda_ip1=delta_lambda_ip1+dd_delta_lambda_ip1
           l_iter_nw=l_iter_nw+ONE_int
        END IF
        !
        ! catch non-convergence or back-track failure
        IF (l_iter_nw>n_iter_max_nw) THEN
           n_conv_nw=.FALSE.
           conv_issue=.TRUE.
        END IF
     END IF
     !
  END DO ! l_iter_nw
  !
END SUBROUTINE newton_back_track_solver_crystal_visco_plasticity
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! computation of orthonormal basis for different subspaces
! **********************************************************************
SUBROUTINE orthonormal_basis_subspace(orth_basis,orth_proj,symm_class,&
     &Rot_tens)
  !
  ! compute the orthonormal bases for each subspace dependent on
  ! symmetry class chosen
  ! Based on: Mahnken, "Anisotropy in geometrically non-linear
  ! elasticity with generalized Seth-Hill strain tensors projected to
  ! invariant subspaces, CNME, 21, p.405, 2005
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:DYADIC_PRODUCT_SECOND_ORDER_TENSORS,&
       & arrange_condense_4_order_tensor
  !
  ! declaration of variables
  !
  ! symm_class......symmetry class of crystal (1=isotropic, 2=cubic)
  ! Rot_tens........rotation tensor
  ! N(i)............orthogonal basis (i=1...6)
  ! orth_basis......all second order, orthogonal basis N(i), rotated
  !                 in case of cubic elasticity
  ! orth_proj.......symmetric dyadic products of orthogonal basis
  !                 i.e. N(i) dyadic N(i)
  !
  IMPLICIT NONE
  !
  INTEGER(KIND=ik),INTENT(IN)::symm_class
  INTEGER(KIND=ik)::i_count
  !
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::Rot_tens
  REAL(KIND=rk),DIMENSION(3,3)::N1,N2,N3,N4,N5,N6
  REAL(KIND=rk),DIMENSION(9,9)::proj_9b9
  REAL(KIND=rk),DIMENSION(6,6)::proj_6b6
  REAL(KIND=rk),DIMENSION(6,3,3),INTENT(OUT)::orth_basis
  REAL(KIND=rk),DIMENSION(3,3,3,3)::proj_t
  REAL(KIND=rk),DIMENSION(6,6,6),INTENT(OUT)::orth_proj
  !
  ! setting the orthogonal basis
  N1=Ident/SQRT(THREE_r)
  N2=RESHAPE((/ZERO_r,ZERO_r,ZERO_r,&
       &ZERO_r,ONE_r,ZERO_r,&
       &ZERO_r,ZERO_r,-ONE_r/),(/3,3/))/SQRT(TWO_r)
  N3=RESHAPE((/-TWO_r,ZERO_r,ZERO_r,&
       &ZERO_r,ONE_r,ZERO_r,&
       &ZERO_r,ZERO_r,ONE_r/),(/3,3/))/SQRT(SIX_r)
  N4=RESHAPE((/ZERO_r,ONE_r,ZERO_r,&
       &ONE_r,ZERO_r,ZERO_r,&
       &ZERO_r,ZERO_r,ZERO_r/),(/3,3/))/SQRT(TWO_r)
  N5=RESHAPE((/ZERO_r,ZERO_r,ZERO_r,&
       &ZERO_r,ZERO_r,ONE_r,&
       &ZERO_r,ONE_r,ZERO_r/),(/3,3/))/SQRT(TWO_r)
  N6=RESHAPE((/ZERO_r,ZERO_r,ONE_r,&
       &ZERO_r,ZERO_r,ZERO_r,&
       &ONE_r,ZERO_r,ZERO_r/),(/3,3/))/SQRT(TWO_r)
  !
  SELECT CASE (symm_class)
  CASE (ZERO_int)! isotropic elasticity (volumetric-isochoric split)
  CASE (ONE_int) ! isotropic elasticity based on Seth-Hill strain
  CASE (TWO_int) ! cubic elasticity based on Seth-Hill strain
     ! N1=MATMUL(Rot_tens,MATMUL(N1,TRANSPOSE(Rot_tens)))
     ! N2=MATMUL(Rot_tens,MATMUL(N2,TRANSPOSE(Rot_tens)))
     ! N3=MATMUL(Rot_tens,MATMUL(N3,TRANSPOSE(Rot_tens)))
     ! N4=MATMUL(Rot_tens,MATMUL(N4,TRANSPOSE(Rot_tens)))
     ! N5=MATMUL(Rot_tens,MATMUL(N5,TRANSPOSE(Rot_tens)))
     ! N6=MATMUL(Rot_tens,MATMUL(N6,TRANSPOSE(Rot_tens)))
  CASE DEFAULT
     WRITE(16,'(A)')'Unknown symmetry class in "orthonormal_basis_subspace"'
     WRITE(16,'(A)')'STOPPING COMPUTATION'
     STOP
  END SELECT
  !
  orth_basis(1,:,:)=N1; orth_basis(2,:,:)=N2
  orth_basis(3,:,:)=N3; orth_basis(4,:,:)=N4
  orth_basis(5,:,:)=N5; orth_basis(6,:,:)=N6
  !
  DO i_count=ONE_int,SIX_int
     proj_t=DYADIC_PRODUCT_SECOND_ORDER_TENSORS(&
          &orth_basis(i_count,:,:),orth_basis(i_count,:,:))
     CALL arrange_condense_4_order_tensor(proj_t,proj_9b9,proj_6b6,&
          &'6b6')
     orth_proj(i_count,:,:)=proj_6b6
  END DO
  !
END SUBROUTINE orthonormal_basis_subspace
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! compute projection tensors for transformation of stress (dual stress
! to 2nd Piola Kirchhoff stress) and tangent moduli (elasticity moduli
! to material tangent moduli)
! **********************************************************************
SUBROUTINE projection_tensors_stress_tangent_transformation(P_tens,&
     &T_dc_L,Proj_tens,d_vect,f_vect,eta,xi,theta,n_distinct,T_stress,&
     &env)
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:symm_tensor_product_second_order_tensors
  !
  ! declaration of variables
  !
  ! P_tens......projection tensor for stress and tangent
  ! T_dc_L......projection tensor for tangent
  !
  IMPLICIT NONE
  !
  INTEGER(KIND=ik),INTENT(IN)::n_distinct
  INTEGER(KIND=ik)::p_pts,q_pts,r_pts
  !
  REAL(KIND=rk),DIMENSION(6,6),INTENT(OUT)::P_tens,T_dc_L
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::T_stress
  REAL(KIND=rk),DIMENSION(3,3)::Proj_stress_pp,Proj_stress_qq
  REAL(KIND=rk),DIMENSION(3,3)::Proj_stress_pq,Proj_stress_qp
  REAL(KIND=rk),DIMENSION(6,6)::D_tens_ppp,D_tens_ppp_t,K_tens_t
  REAL(KIND=rk),DIMENSION(6,6)::D_tens_p_qq,D_tens_qq_p
  REAL(KIND=rk),DIMENSION(6,6)::D_tens_q_pq,D_tens_pq_q
  REAL(KIND=rk),DIMENSION(6,6)::D_tens_q_qp,D_tens_qp_q
  REAL(KIND=rk),DIMENSION(6,6)::D_tens_r_pq,D_tens_pq_r
  REAL(KIND=rk),DIMENSION(3,3,3),INTENT(IN)::Proj_tens
  REAL(KIND=rk),DIMENSION(3),INTENT(IN)::d_vect,f_vect
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::xi,theta
  REAL(KIND=rk),INTENT(IN)::eta
  !
  CHARACTER(4),INTENT(IN)::env
  !
  P_tens=ZERO_r; T_dc_L=ZERO_r
  DO p_pts=ONE_int,n_distinct
     Proj_stress_pp=MATMUL(Proj_tens(p_pts,:,:),MATMUL(T_stress,&
          &Proj_tens(p_pts,:,:)))
     D_tens_ppp=symm_tensor_product_second_order_tensors(&
          &Proj_tens(p_pts,:,:),Proj_stress_pp,env)
     D_tens_ppp_t=symm_tensor_product_second_order_tensors(&
          &Proj_stress_pp,Proj_tens(p_pts,:,:),env)
     T_dc_L=T_dc_L+ONE_FOURTH_r*f_vect(p_pts)*(D_tens_ppp+D_tens_ppp_t)
     !
     K_tens_t=symm_tensor_product_second_order_tensors(&
          &Proj_tens(p_pts,:,:),Proj_tens(p_pts,:,:),env)
     P_tens=P_tens+ONE_HALF_r*d_vect(p_pts)*K_tens_t
     !
     DO q_pts=ONE_int,n_distinct
        IF (p_pts /= q_pts) THEN
           Proj_stress_qq=MATMUL(Proj_tens(q_pts,:,:),MATMUL(T_stress,&
                &Proj_tens(q_pts,:,:)))
           Proj_stress_pq=MATMUL(Proj_tens(p_pts,:,:),MATMUL(T_stress,&
                &Proj_tens(q_pts,:,:)))
           Proj_stress_qp=MATMUL(Proj_tens(q_pts,:,:),MATMUL(T_stress,&
                &Proj_tens(p_pts,:,:)))
           D_tens_p_qq=symm_tensor_product_second_order_tensors(&
                &Proj_tens(p_pts,:,:),Proj_stress_qq,env)
           D_tens_qq_p=symm_tensor_product_second_order_tensors(&
                &Proj_stress_qq,Proj_tens(p_pts,:,:),env)
           D_tens_q_pq=symm_tensor_product_second_order_tensors(&
                &Proj_tens(q_pts,:,:),Proj_stress_pq,env)
           D_tens_pq_q=symm_tensor_product_second_order_tensors(&
                &Proj_stress_pq,Proj_tens(q_pts,:,:),env)
           D_tens_q_qp=symm_tensor_product_second_order_tensors(&
                &Proj_tens(q_pts,:,:),Proj_stress_qp,env)
           D_tens_qp_q=symm_tensor_product_second_order_tensors(&
                &Proj_stress_qp,Proj_tens(q_pts,:,:),env)
           T_dc_L=T_dc_L+TWO_r*xi(p_pts,q_pts)*(D_tens_p_qq+&
                &D_tens_qq_p+D_tens_q_pq+D_tens_pq_q+D_tens_q_qp+&
                &D_tens_qp_q)
           !
           K_tens_t=symm_tensor_product_second_order_tensors(&
                &Proj_tens(p_pts,:,:),Proj_tens(q_pts,:,:),env)
           P_tens=P_tens+theta(p_pts,q_pts)*K_tens_t;
           DO r_pts=ONE_int,n_distinct
              IF ((r_pts/=p_pts).AND.(r_pts/=q_pts)) THEN
                 D_tens_r_pq=symm_tensor_product_second_order_tensors(&
                      &Proj_tens(r_pts,:,:),Proj_stress_pq,env)
                 D_tens_pq_r=symm_tensor_product_second_order_tensors(&
                      &Proj_stress_pq,Proj_tens(r_pts,:,:),env)
                 T_dc_L=T_dc_L+TWO_r*eta*(D_tens_r_pq+D_tens_pq_r)
              END IF
           END DO
        END IF
     END DO ! q
  END DO ! p
  !
END SUBROUTINE projection_tensors_stress_tangent_transformation
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! compute residual and jacobian for Newton-type method to determine
! incremental lagrange multipliers / incremental plastic slips
! **********************************************************************
SUBROUTINE residual_tangent_computation_nw_crystal_plasticity(&
     &residual_nw,jacobian_nw,active_slip_system,&
     &deriv_visco_plast_potent,unimod_C,F_p_inv_iter,F_p_inv_n,Je,&
     &plastic_slip_iter,delta_lambda_i,delta_lambda_ip1,algo,dt,m_SH,&
     &elast_const,symm_class,orth_basis,orth_proj,mat_param_inelast,&
     &NDI,NSHR,env,eta,hard_law,theta_KK,Phi_alpha_all)
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:comp_sec_ord_tens_rank_one, &
       & transform_symm_tens_vec_form
  USE constitutive_routines_interface, only:linearization_mandel_stress,&
       & linearization_mandel_stress_elasticity_SH_strain, &
       & assembly_residual_jacobian_newton, &
       & eval_yield_function_linearization, &
       & isotr_elasticity_mandel_stress_volumetric_isochoric_decoupled, &
       & resolved_shear_stress_linearization, slip_system_fcc, &
       & S2PK_IC_cubic_elasticity_seth_hill
  !
  ! variable declaration
  ! residual_nw...............residual vector
  ! jacobian_nw...............jacobian matrix
  ! deriv_visco_plast_potent..derivative of the viscoplast potential
  !                           aux. var in rate-independent case
  ! active_slip_system........list of active slip systems
  ! unimod_C..................unimodular right Cauchy-Green (only
  !                           unimodular for isotr-decoupled elasticity)
  ! F_p_inv_iter..............inverse of plastic deformation gradient
  ! F_p_inv_n.................inv. of plastic deform. grad. t_{n}
  ! Je........................determinant of elastic deform. grad.
  ! plastic_slip_iter.........accumulated plastic slip on slip systems
  ! delta_lambda_i............incr. lagr. multipl. (fix-point iter i)
  ! delta_lambda_ip1..........incr. lagr. multipl. (fix-point iter i+1)
  ! algo......................id for algorithm employed
  !                           0...augmented-lagrangian
  !                           1...NCP-function (Kanzow-Kleinmichel)
  !                           2...Perzyna viscoplastic formulation
  !                           3...Ortiz/Miehe viscoplastic formulation
  ! dt........................time increment
  ! m_SH......................Seth-Hill parameter (defining the strain measure)
  ! elast_const...............elastic constants (isotr.: kappa,mu;
  !                           cubic: C11,C12,C44)
  ! symm_class................symmetry class (1=isotropic, 2=cubic)
  ! orth_basis................container of orthogonal basis (second order
  !                           tensor)
  ! orth_proj.................container of dyadic product of orthogonal
  !                           basis
  ! mat_param_inelast.........mat. param.r assoc. with inelastic behavior
  ! env.......................environment in which routine is called
  ! eta.......................pseudo-viscosity in augmented lagrangian
  ! hard_law..................id of hardening law
  !                        0....nonlin. Taylor-type (Schmidt-Baldassari)
  !                        1....nonlin. interaction matrix (Cailletaud,Forest)
  ! theta_KK..................algorithmic parameter for Kanzow-Kleinmichel
  !                           NCP-function
  !
  IMPLICIT NONE
  REAL(KIND=rk),DIMENSION(24),INTENT(OUT)::residual_nw
  REAL(KIND=rk),DIMENSION(24),INTENT(OUT)::deriv_visco_plast_potent
  REAL(KIND=rk),DIMENSION(24),INTENT(OUT)::Phi_alpha_all
  REAL(KIND=rk),DIMENSION(24,24),INTENT(OUT)::jacobian_nw
  REAL(KIND=rk),INTENT(IN)::Je,dt,eta,m_SH,theta_KK
  REAL(KIND=rk),DIMENSION(3),INTENT(IN)::elast_const
  REAL(KIND=rk),DIMENSION(13),INTENT(IN)::mat_param_inelast
  REAL(KIND=rk),DIMENSION(24),INTENT(IN)::plastic_slip_iter
  REAL(KIND=rk),DIMENSION(24),INTENT(IN)::delta_lambda_i
  REAL(KIND=rk),DIMENSION(24),INTENT(IN)::delta_lambda_ip1
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::unimod_C,F_p_inv_iter
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::F_p_inv_n
  REAL(KIND=rk),DIMENSION(6,3,3),INTENT(IN)::orth_basis
  REAL(KIND=rk),DIMENSION(6,6,6),INTENT(IN)::orth_proj
  REAL(KIND=rk)::kappa,mu,M_hyd,tau_alpha_rss,Phi_alpha
  REAL(KIND=rk),DIMENSION(3)::normal_dir,slip_dir
  REAL(KIND=rk),DIMENSION(6)::S_2PK_vec
  REAL(KIND=rk),DIMENSION(24)::d_tau_alpha_d_delta_lambda_ip1_beta
  REAL(KIND=rk),DIMENSION(24)::dPhi_alpha_d_delta_lambda_ip1_beta
  REAL(KIND=rk),DIMENSION(3,3)::unimod_Ce,unimod_C_F_p_inv_iter
  REAL(KIND=rk),DIMENSION(3,3)::M_dev,S_2PK_IC,M_IC,N_alpha
  REAL(KIND=rk),DIMENSION(6,6)::C_mat_tang
  REAL(KIND=rk),DIMENSION(9,24)::d_Mandel_stress_d_delta_lambda_ip1_beta
  INTEGER(KIND=ik),INTENT(IN)::algo,symm_class,NDI,NSHR,hard_law
  INTEGER(KIND=ik)::i_rtc,i_slip
  !
  LOGICAL,DIMENSION(24),INTENT(OUT)::active_slip_system
  CHARACTER(4),INTENT(IN)::env
  !
  ! unimodular part of elastic right Cauchy-Green stretch (Ce is only
  ! unimodular if the isotr_decoup elastic law is employed)
  unimod_Ce=MATMUL(TRANSPOSE(F_p_inv_iter),MATMUL(unimod_C,F_p_inv_iter))
  !
  ! single contraction of unimodular C and inverse of F_p
  unimod_C_F_p_inv_iter=MATMUL(unimod_C,F_p_inv_iter)
  !
  SELECT CASE (symm_class)
     ! ---------------------------------------------------------------
  CASE (ZERO_int)! volumetric/isochoric decoupled, isotropic elasticity
     ! ---------------------------------------------------------------
     kappa=elast_const(1); mu=elast_const(2)
     CALL isotr_elasticity_mandel_stress_volumetric_isochoric_decoupled(&
          &M_hyd,M_dev,Je,unimod_Ce,kappa,mu)
     !
     ! derivative of Mandel stress w.r.t. increment in Lagrange multiplier
     d_Mandel_stress_d_delta_lambda_ip1_beta=linearization_mandel_stress(&
          &F_p_inv_n,unimod_C_F_p_inv_iter,mu,env)
     ! ---------------------------------------------------------------
  CASE (ONE_int) ! isotropic elasticity based on Seth-Hill strain
     ! ---------------------------------------------------------------
     CALL S2PK_IC_cubic_elasticity_seth_hill(S_2PK_vec,C_mat_tang,&
          &unimod_Ce,m_SH,elast_const,symm_class,orth_basis,orth_proj,&
          &NDI,NSHR,env)
     CALL transform_symm_tens_vec_form(S_2PK_IC,S_2PK_vec,ONE_int,env)
     !
     ! compute Mandel stress (Ce is stored in unimod_Ce)
     M_IC=MATMUL(unimod_Ce,S_2PK_IC);    M_dev=M_IC
     !
     ! derivative of Mandel stress w.r.t. increment in Lagrange multiplier
     d_Mandel_stress_d_delta_lambda_ip1_beta=&
          &linearization_mandel_stress_elasticity_SH_strain(F_p_inv_n,&
          &unimod_C_F_p_inv_iter,unimod_Ce,S_2PK_IC,C_mat_tang,env)
     !
     ! ---------------------------------------------------------------
  CASE (TWO_int) ! cubic elasticity based on Seth-Hill strain
     ! ---------------------------------------------------------------
     CALL S2PK_IC_cubic_elasticity_seth_hill(S_2PK_vec,C_mat_tang,&
          &unimod_Ce,m_SH,elast_const,symm_class,orth_basis,orth_proj,&
          &NDI,NSHR,env)
     CALL transform_symm_tens_vec_form(S_2PK_IC,S_2PK_vec,ONE_int,env)
     !
     ! compute Mandel stress (Ce is stored in unimod_Ce)
     M_IC=MATMUL(unimod_Ce,S_2PK_IC);    M_dev=M_IC
     !
     ! derivative of Mandel stress w.r.t. increment in Lagrange multiplier
     d_Mandel_stress_d_delta_lambda_ip1_beta=&
          &linearization_mandel_stress_elasticity_SH_strain(F_p_inv_n,&
          &unimod_C_F_p_inv_iter,unimod_Ce,S_2PK_IC,C_mat_tang,env)
     !
  CASE DEFAULT
     WRITE(16,'(A)')'Unknown symmetry class in "residual_tangent_computation_nw_crystal_plasticity"'
     WRITE(16,'(A)')'STOPPING COMPUTATION'
     STOP
  END SELECT
  !
  ! initialize residual and jacobian
  residual_nw=ZERO_r
  jacobian_nw=ZERO_r
  SELECT CASE (algo)
  CASE (ZERO_int) ! augmented lagrangian approach
     FORALL (i_rtc=1:24)  jacobian_nw(i_rtc, i_rtc) = ONE_r ! identity
  CASE (ONE_int)  ! NCP-function (Kanzow-Kleinmichel)
  CASE (TWO_int)  ! Perzyna viscoplastic formulation
     FORALL (i_rtc=1:24)  jacobian_nw(i_rtc, i_rtc) = ONE_r ! identity
  CASE (THREE_int)! Ortiz/Miehe viscoplastic formulation
  CASE DEFAULT
     WRITE(16,'(A)')'Unknown algorithm id used in "residual_tangent_computation_nw_crystal_plasticity"'
  END SELECT
  !
  active_slip_system=.FALSE.
  deriv_visco_plast_potent=ZERO_r
  Phi_alpha_all=ZERO_r
  !
  DO i_slip=ONE_int,TWENTY_FOUR_int
     CALL slip_system_fcc(normal_dir,slip_dir,i_slip)
     N_alpha=comp_sec_ord_tens_rank_one(slip_dir,normal_dir)
     !
     CALL resolved_shear_stress_linearization(tau_alpha_rss,&
          &d_tau_alpha_d_delta_lambda_ip1_beta,M_dev,N_alpha,&
          &d_Mandel_stress_d_delta_lambda_ip1_beta,env)
     !
     CALL eval_yield_function_linearization (Phi_alpha,&
          &dPhi_alpha_d_delta_lambda_ip1_beta,tau_alpha_rss,&
          &d_tau_alpha_d_delta_lambda_ip1_beta,plastic_slip_iter,&
          &i_slip,algo,mat_param_inelast,hard_law)
     Phi_alpha_all(i_slip)=Phi_alpha
     !
     ! assemble the contribution to the residual and the tangent
     CALL assembly_residual_jacobian_newton (residual_nw,jacobian_nw,&
          & active_slip_system,deriv_visco_plast_potent,Phi_alpha,&
          & dPhi_alpha_d_delta_lambda_ip1_beta,delta_lambda_i,&
          & delta_lambda_ip1,plastic_slip_iter,algo,dt,eta,i_slip,&
          & mat_param_inelast,theta_KK,hard_law) 
     !
  END DO
  !
END SUBROUTINE residual_tangent_computation_nw_crystal_plasticity
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! combination of update procedure and residual / jacobian computation in
! order to use it with generic solvers for nonlinear equations
! **********************************************************************
SUBROUTINE residual_tangent_nw_crystal_plasticity_generic_solver(&
     &residual_nw,jacobian_nw,unimod_C,F_p_inv_n,Je,plastic_slip_n,&
     &delta_lambda_i,delta_lambda_ip1,algo,dt,m_SH,elast_const,&
     &symm_class,orth_basis,orth_proj,mat_param_inelast,NDI,NSHR,env,&
     &eta,hard_law,theta_KK)
  !
  USE Abaqus_Interface
  USE Constants
  USE constitutive_routines_interface, only:update_plastic_deformation_gradient,&
       & residual_tangent_computation_nw_crystal_plasticity
  !
  ! variable declaration
  ! residual_nw...............residual vector
  ! jacobian_nw...............jacobian matrix
  ! unimod_C..................unimodular right Cauchy-Green (only
  !                           unimodular for isotr-decoupled elasticity)
  ! F_p_inv_n.................inverse of plastic deform. grad. t_{n}
  ! Je........................determinant of elastic deform. grad.
  ! plastic_slip_n............accumulated plastic slip on slip systems
  ! delta_lambda_i............incr. lagr. multipl. (fix-point iter i)
  ! delta_lambda_ip1..........incr. lagr. multipl. (fix-point iter i+1)
  ! algo......................id for algorithm employed
  !                           0...augmented-lagrangian
  !                           1...NCP-function (Kanzow-Kleinmichel)
  !                           2...Perzyna viscoplastic formulation
  !                           3...Ortiz/Miehe viscoplastic formulation
  ! dt........................time increment
  ! m_SH......................Seth-Hill parameter (defining the strain measure)
  ! elast_const...............elastic constants (isotr.: kappa,mu;
  !                           cubic: C11,C12,C44)
  ! symm_class................symmetry class (1=isotropic, 2=cubic)
  ! orth_basis................container of orthogonal basis (second order
  !                           tensor)
  ! orth_proj.................container of dyadic product of orthogonal
  !                           basis
  ! mat_param_inelast.........mat. param.r assoc. with inelastic behavior
  ! env.......................environment in which routine is called
  ! eta.......................pseudo-viscosity in augmented lagrangian
  ! hard_law..................id of hardening law
  !                        0....nonlin. Taylor-type (Schmidt-Baldassari)
  !                        1....nonlin. interaction matrix (Cailletaud,Forest)
  ! theta_KK..................algorithmic parameter for Kanzow-Kleinmichel
  !                           NCP-function
  IMPLICIT NONE
  REAL(KIND=rk),DIMENSION(24),INTENT(OUT)::residual_nw
  REAL(KIND=rk),DIMENSION(24,24),INTENT(OUT)::jacobian_nw
  REAL(KIND=rk),INTENT(IN)::Je,dt,eta,m_SH,theta_KK
  REAL(KIND=rk),DIMENSION(3),INTENT(IN)::elast_const
  REAL(KIND=rk),DIMENSION(13),INTENT(IN)::mat_param_inelast
  REAL(KIND=rk),DIMENSION(24),INTENT(IN)::plastic_slip_n
  REAL(KIND=rk),DIMENSION(24),INTENT(IN)::delta_lambda_i
  REAL(KIND=rk),DIMENSION(24),INTENT(IN)::delta_lambda_ip1
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::unimod_C
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::F_p_inv_n
  REAL(KIND=rk),DIMENSION(6,3,3),INTENT(IN)::orth_basis
  REAL(KIND=rk),DIMENSION(6,6,6),INTENT(IN)::orth_proj
  REAL(KIND=rk),DIMENSION(24)::deriv_visco_plast_potent,Phi_alpha_dum
  REAL(KIND=rk),DIMENSION(24)::plastic_slip_iter
  REAL(KIND=rk),DIMENSION(3,3)::F_p_inv_iter
  INTEGER(KIND=ik),INTENT(IN)::algo,symm_class,NDI,NSHR,hard_law
  !
  LOGICAL,DIMENSION(24)::active_slip_system
  CHARACTER(4),INTENT(IN)::env
  !
  ! compute hardening variables on all slip systems
  plastic_slip_iter=plastic_slip_n+delta_lambda_ip1
  !
  ! update plastic part of the deformation gradient
  F_p_inv_iter=update_plastic_deformation_gradient(F_p_inv_n,&
       &delta_lambda_ip1)
  !
  CALL residual_tangent_computation_nw_crystal_plasticity(residual_nw,&
       &jacobian_nw,active_slip_system,deriv_visco_plast_potent,&
       &unimod_C,F_p_inv_iter,F_p_inv_n,Je,plastic_slip_iter,delta_lambda_i,&
       &delta_lambda_ip1,algo,dt,m_SH,elast_const,symm_class,orth_basis,&
       &orth_proj,mat_param_inelast,NDI,NSHR,env,eta,hard_law,theta_KK,Phi_alpha_dum)
  !
END SUBROUTINE residual_tangent_nw_crystal_plasticity_generic_solver
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! compute the resolved shear stress on current slip system and the
! corresponding derivatives w.r.t. increment of lagrange multipliers
! **********************************************************************
SUBROUTINE resolved_shear_stress_linearization(tau_alpha_rss,&
     &d_tau_alpha_d_delta_lambda_ip1_beta,M_dev,N_alpha,&
     &d_Mandel_stress_d_delta_lambda_ip1_beta,env)
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:double_contraction_second_order_tensors,&
       & transform_tens_vec_form
  !
  ! variable declaration
  ! tau_alpha_rss........resolved shear stress on slip system alpha
  ! d_tau_alpha_d_delta_lambda_ip1_beta..deriv. of resolved shear
  !                                      stress on system alpha w.r.t.
  !                                      incr. lagrange multipl. beta
  ! M_dev................Mandel stress tensor (deviatoric if volumetric-
  !                      isochoric split isotropic elasticity is used)
  ! N_alpha..............rank one tensor defining the slip system
  !                      (slip direction dyad slip plane normal)
  ! d_Mandel_stress_d_delta_lambda_ip1_beta..deriv of Mandel stress
  !                      tensor w.r.t. incremental lagrange multipl.
  !
  IMPLICIT NONE
  REAL(KIND=rk),INTENT(OUT)::tau_alpha_rss
  REAL(KIND=rk),DIMENSION(24),INTENT(OUT)::d_tau_alpha_d_delta_lambda_ip1_beta
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::M_dev,N_alpha
  REAL(KIND=rk),DIMENSION(9,24),INTENT(IN)::d_Mandel_stress_d_delta_lambda_ip1_beta
  REAL(KIND=rk),DIMENSION(3,3)::d_Mandel_stress_d_delta_lambda_beta
  INTEGER(KIND=ik)::i_slip
  !
  CHARACTER(4),INTENT(IN)::env
  !
  ! resolved shear stress
  tau_alpha_rss=double_contraction_second_order_tensors(M_dev,N_alpha)
  !
  ! derivative of resolved shear stress w.r.t. incr. lagrange multiplier
  d_tau_alpha_d_delta_lambda_ip1_beta=ZERO_r
  !
  DO i_slip=ONE_int,TWENTY_FOUR_int
     CALL transform_tens_vec_form(d_Mandel_stress_d_delta_lambda_beta,&
          &d_Mandel_stress_d_delta_lambda_ip1_beta(:,i_slip),ONE_int,env)
     d_tau_alpha_d_delta_lambda_ip1_beta(i_slip) = &
          & double_contraction_second_order_tensors(&
          &d_Mandel_stress_d_delta_lambda_beta,N_alpha)
  END DO
  !
END SUBROUTINE resolved_shear_stress_linearization
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! provide the slip normal and slip direction of the given index
! **********************************************************************
SUBROUTINE slip_system_fcc(normal_dir,slip_dir,index)
  !
  USE Abaqus_Interface
  USE Constants
  !
  IMPLICIT NONE
  INTEGER(KIND=ik),INTENT(IN)::index
  REAL(KIND=rk),DIMENSION(3),INTENT(OUT)::normal_dir,slip_dir
  !
  SELECT CASE (index)
  CASE (1)
     normal_dir=(/ONE_by_sqrt_three_r,ONE_by_sqrt_three_r,ONE_by_sqrt_three_r/)
     slip_dir=(/ONE_by_sqrt_two_r,-ONE_by_sqrt_two_r,ZERO_r/)
  CASE (2)
     normal_dir=(/ONE_by_sqrt_three_r,ONE_by_sqrt_three_r,ONE_by_sqrt_three_r/)
     slip_dir=(/ONE_by_sqrt_two_r,ZERO_r,-ONE_by_sqrt_two_r/)
  CASE (3)
     normal_dir=(/ONE_by_sqrt_three_r,ONE_by_sqrt_three_r,ONE_by_sqrt_three_r/)
     slip_dir=(/ZERO_r,ONE_by_sqrt_two_r,-ONE_by_sqrt_two_r/)
     ! --------------------------------------------------------------------
  CASE (4)
     normal_dir=(/ONE_by_sqrt_three_r,ONE_by_sqrt_three_r,-ONE_by_sqrt_three_r/)
     slip_dir=(/ONE_by_sqrt_two_r,-ONE_by_sqrt_two_r,ZERO_r/)
  CASE (5)
     normal_dir=(/ONE_by_sqrt_three_r,ONE_by_sqrt_three_r,-ONE_by_sqrt_three_r/)
     slip_dir=(/ONE_by_sqrt_two_r,ZERO_r,ONE_by_sqrt_two_r/)
  CASE (6)
     normal_dir=(/ONE_by_sqrt_three_r,ONE_by_sqrt_three_r,-ONE_by_sqrt_three_r/)
     slip_dir=(/ZERO_r,ONE_by_sqrt_two_r,ONE_by_sqrt_two_r/)
     ! --------------------------------------------------------------------
  CASE (7)
     normal_dir=(/ONE_by_sqrt_three_r,-ONE_by_sqrt_three_r,ONE_by_sqrt_three_r/)
     slip_dir=(/ONE_by_sqrt_two_r,ONE_by_sqrt_two_r,ZERO_r/)
  CASE (8)
     normal_dir=(/ONE_by_sqrt_three_r,-ONE_by_sqrt_three_r,ONE_by_sqrt_three_r/)
     slip_dir=(/ONE_by_sqrt_two_r,ZERO_r,-ONE_by_sqrt_two_r/)
  CASE (9)
     normal_dir=(/ONE_by_sqrt_three_r,-ONE_by_sqrt_three_r,ONE_by_sqrt_three_r/)
     slip_dir=(/ZERO_r,ONE_by_sqrt_two_r,ONE_by_sqrt_two_r/)
     ! --------------------------------------------------------------------
  CASE (10)
     normal_dir=(/-ONE_by_sqrt_three_r,ONE_by_sqrt_three_r,ONE_by_sqrt_three_r/)
     slip_dir=(/ONE_by_sqrt_two_r,ONE_by_sqrt_two_r,ZERO_r/)
  CASE (11)
     normal_dir=(/-ONE_by_sqrt_three_r,ONE_by_sqrt_three_r,ONE_by_sqrt_three_r/)
     slip_dir=(/ONE_by_sqrt_two_r,ZERO_r,ONE_by_sqrt_two_r/)
  CASE (12)
     normal_dir=(/-ONE_by_sqrt_three_r,ONE_by_sqrt_three_r,ONE_by_sqrt_three_r/)
     slip_dir=(/ZERO_r,ONE_by_sqrt_two_r,-ONE_by_sqrt_two_r/)
     ! --------------------------------------------------------------------
  CASE (13)
     normal_dir=(/ONE_by_sqrt_three_r,ONE_by_sqrt_three_r,ONE_by_sqrt_three_r/)
     slip_dir=-(/ONE_by_sqrt_two_r,-ONE_by_sqrt_two_r,ZERO_r/)
  CASE (14)
     normal_dir=(/ONE_by_sqrt_three_r,ONE_by_sqrt_three_r,ONE_by_sqrt_three_r/)
     slip_dir=-(/ONE_by_sqrt_two_r,ZERO_r,-ONE_by_sqrt_two_r/)
  CASE (15)
     normal_dir=(/ONE_by_sqrt_three_r,ONE_by_sqrt_three_r,ONE_by_sqrt_three_r/)
     slip_dir=-(/ZERO_r,ONE_by_sqrt_two_r,-ONE_by_sqrt_two_r/)
     ! --------------------------------------------------------------------
  CASE (16)
     normal_dir=(/ONE_by_sqrt_three_r,ONE_by_sqrt_three_r,-ONE_by_sqrt_three_r/)
     slip_dir=-(/ONE_by_sqrt_two_r,-ONE_by_sqrt_two_r,ZERO_r/)
  CASE (17)
     normal_dir=(/ONE_by_sqrt_three_r,ONE_by_sqrt_three_r,-ONE_by_sqrt_three_r/)
     slip_dir=-(/ONE_by_sqrt_two_r,ZERO_r,ONE_by_sqrt_two_r/)
  CASE (18)
     normal_dir=(/ONE_by_sqrt_three_r,ONE_by_sqrt_three_r,-ONE_by_sqrt_three_r/)
     slip_dir=-(/ZERO_r,ONE_by_sqrt_two_r,ONE_by_sqrt_two_r/)
     ! ---------------------------------------------------------------------------
  CASE(19)
     normal_dir=(/ONE_by_sqrt_three_r,-ONE_by_sqrt_three_r,ONE_by_sqrt_three_r/)
     slip_dir=-(/ONE_by_sqrt_two_r,ONE_by_sqrt_two_r,ZERO_r/)
  CASE(20)
     normal_dir=(/ONE_by_sqrt_three_r,-ONE_by_sqrt_three_r,ONE_by_sqrt_three_r/)
     slip_dir=-(/ONE_by_sqrt_two_r,ZERO_r,-ONE_by_sqrt_two_r/)
  CASE(21)
     normal_dir=(/ONE_by_sqrt_three_r,-ONE_by_sqrt_three_r,ONE_by_sqrt_three_r/)
     slip_dir=-(/ZERO_r,ONE_by_sqrt_two_r,ONE_by_sqrt_two_r/)
     ! ---------------------------------------------------------------------------
  CASE(22)
     normal_dir=(/-ONE_by_sqrt_three_r,ONE_by_sqrt_three_r,ONE_by_sqrt_three_r/)
     slip_dir=-(/ONE_by_sqrt_two_r,ONE_by_sqrt_two_r,ZERO_r/)
  CASE(23)
     normal_dir=(/-ONE_by_sqrt_three_r,ONE_by_sqrt_three_r,ONE_by_sqrt_three_r/)
     slip_dir=-(/ONE_by_sqrt_two_r,ZERO_r,ONE_by_sqrt_two_r/)
  CASE(24)
     normal_dir=(/-ONE_by_sqrt_three_r,ONE_by_sqrt_three_r,ONE_by_sqrt_three_r/)
     slip_dir=-(/ZERO_r,ONE_by_sqrt_two_r,-ONE_by_sqrt_two_r/)
     ! ---------------------------------------------------------------------------
  CASE DEFAULT
     WRITE(16,'(A)')'Unknown number of slip system in: slip_system_fcc'
     WRITE(16,'(A)')'STOPPING COMPUTATION'
     STOP
  END SELECT
  !
END SUBROUTINE slip_system_fcc
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! performs the return mapping onto the set of active slip systems, if
! the yield criteria are violated and computes the stress tensors based
! on the updated crystallograpic slip / plastic deformation gradient
! **********************************************************************
SUBROUTINE stress_computation_single_crystal_plasticity(S2PK,S2PK_IC,&
     &F_p_inv_n1,plastic_slip_iter,d_delta_lambda_ip1_d_C,&
     &active_slip_system,Ce_n1,Ce_n1_inv,Je_n1,C,det_Fp_inv,&
     &F_p_inv_iter,F,dt,F_p_inv_n,plastic_slip_n,elast_const,&
     &mat_param_inelast,m_SH,symm_class,orth_basis,orth_proj,NDI,NSHR,&
     &algo,env,hard_law,theta_KK,aver_adj_stp,conv_issue_dil)
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:comp_det_Matrix, M33INV,&
       & transform_symm_tens_vec_form
  USE constitutive_routines_interface, only:S2PK_IC_isotr_elasticity,&
       & determine_incremental_lagrange_multiplier_augmented_lagrange, &
       & determine_incremental_lagrange_multiplier_miehe_viscoplastic, &
       & S2PK_IC_cubic_elasticity_seth_hill
  !
  ! variable declaration
  ! S2PK................second Piola-Kirchhoff stress (ref. conf.)
  ! S2PK_IC.............second Piola-Kirchhoff stress (intermed. conf.)
  ! F_p_inv_n1..........inverse of plastic part of deformation gradient
  !                     at t_{n+1}
  ! plastic_slip_iter...accumulated plastic slip on slip systems
  ! d_delta_lambda_ip1_d_C..deriv. of incremental lagrange multipliers
  !                         w.r.t. right Cauchy-Green tensor
  ! active_slip_system..list of active slip systems
  ! Ce_n1...............elastic right Cauchy-Green tensor
  ! Ce_n1_inv...........inverse of elastic right Cauchy-Green tensor
  ! Je_n1...............determinant of elastic part of deformation gradient
  ! C...................right Cauchy-Green tensor
  ! det_Fp_inv..........determinant of inverse of plastic part of
  !                     deformation gradient
  ! F_p_inv_iter........inverse of plastic part of deformation gradient
  ! F...................deformation gradient (t_{n+1})
  ! dt..................time increment
  ! F_p_inv_n...........inverse of plastic part of deformation gradient
  !                     at t_n
  ! elast_const.........material parameter for elastic behavior
  ! mat_param_inelast...material parameter assoc. with inelastic behavior  
  ! m_SH................Seth-Hill parameter
  ! symm_class..........symmetry included in elastic law
  ! orth_basis..........all second order, orthogonal basis N(i), rotated
  !                     in case of cubic elasticity
  ! orth_proj...........symmetric dyadic products of orthogonal basis i.e.
  !                     N(i) dyadic N(i)
  ! algo................id for algorithm employed
  !                     0...augmented-lagrangian
  !                     1...NCP-function (Kanzow-Kleinmichel)
  !                     2...Perzyna viscoplastic formulation
  !                     3...Ortiz/Miehe viscoplastic formulation
  ! env.................environment in which routine is called
  ! hard_law............id of hardening law
  !                     0....nonlin. Taylor-type (Schmidt-Baldassari)
  !                     1....nonlin. interaction matrix (Cailletaud,Forest)
  ! theta_KK............algorithmic parameter for Kanzow-Kleinmichel
  !                     NCP-function
  ! aver_adj_stp........average number of steps employed to adaptively
  !                     adjust the visco parameter
  ! conv_issue_dil......error flag indicating non-convergence in
  !                     determination of incremental lagr. multiplier
  !
  IMPLICIT NONE
  REAL(KIND=rk),INTENT(OUT)::Je_n1,det_Fp_inv
  REAL(KIND=rk),DIMENSION(24),INTENT(OUT)::plastic_slip_iter
  REAL(KIND=rk),DIMENSION(3,3),INTENT(OUT)::S2PK,S2PK_IC,F_p_inv_n1
  REAL(KIND=rk),DIMENSION(3,3),INTENT(OUT)::Ce_n1,Ce_n1_inv,C
  REAL(KIND=rk),DIMENSION(3,3),INTENT(OUT)::F_p_inv_iter
  REAL(KIND=rk),DIMENSION(24,6),INTENT(OUT)::d_delta_lambda_ip1_d_C
  REAL(KIND=rk),INTENT(IN)::dt,m_SH,theta_KK
  REAL(KIND=rk),INTENT(INOUT)::aver_adj_stp
  REAL(KIND=rk),DIMENSION(3),INTENT(IN)::elast_const
  REAL(KIND=rk),DIMENSION(13),INTENT(IN)::mat_param_inelast
  REAL(KIND=rk),DIMENSION(24),INTENT(IN)::plastic_slip_n
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::F,F_p_inv_n
  REAL(KIND=rk),DIMENSION(6,3,3),INTENT(IN)::orth_basis
  REAL(KIND=rk),DIMENSION(6,6,6),INTENT(IN)::orth_proj
  REAL(KIND=rk)::kappa,mu
  REAL(KIND=rk),DIMENSION(6)::S2PK_IC_vec
  REAL(KIND=rk),DIMENSION(3,3)::unimod_C,C_inv,F_e_n1,F_e_inv_n1
  REAL(KIND=rk),DIMENSION(6,6)::C_mat_tang
  INTEGER(KIND=ik),INTENT(IN)::algo,symm_class,NDI,NSHR,hard_law
  !
  LOGICAL::OK_FLAG
  LOGICAL,INTENT(OUT)::conv_issue_dil
  LOGICAL,DIMENSION(24),INTENT(OUT)::active_slip_system
  CHARACTER(4),INTENT(IN)::env
  !
  ! compute right Cauchy-Green tensor, its inverse and determinant of F
  ! plastic deformation is isochoric
  Je_n1=comp_det_Matrix(F) 
  C=MATMUL(TRANSPOSE(F),F)
  CALL M33INV (C,C_inv, OK_FLAG)
  IF (.NOT.(OK_FLAG)) THEN
     WRITE(16,'(A)')'ERROR in: "stress_computation_single_crystal_plasticity"'
     WRITE(16,'(A)')'M33INV encountered singular matrix'
     STOP
  END IF
  !
  ! assign unimodular right Cauchy-Green tensor (only for decoupled
  ! isotropic elasticity)
  SELECT CASE (symm_class)
  CASE (ZERO_int)! volumetric/isochoric decoupled, isotropic elasticity
     unimod_C=C*Je_n1**(-TWO_by_THREE_r)
  CASE (ONE_int) ! isotropic elasticity based on Seth-Hill strain
     unimod_C=C
  CASE (TWO_int) ! cubic elasticity based on Seth-Hill strain
     unimod_C=C
  CASE DEFAULT
     WRITE(16,'(2A)')'Unknown symm_class in: "stress_computation_single_crystal_plasticity"'
     WRITE(16,'(A)')'STOPPING COMPUTATION'
     STOP
  END SELECT
  !
  ! determine the incremental lagrange multipliers
  SELECT CASE (algo)
     ! ---------------------------------------------------------------
  CASE (ZERO_int) ! augmented lagrangian approach
     ! --------------------------------------------------------------
     CALL determine_incremental_lagrange_multiplier_augmented_lagrange(&
          &F_p_inv_iter,plastic_slip_iter,d_delta_lambda_ip1_d_C,&
          &active_slip_system,unimod_C,F_p_inv_n,C_inv,Je_n1,plastic_slip_n,&
          &algo,dt,elast_const,mat_param_inelast,m_SH,symm_class,&
          &orth_basis,orth_proj,NDI,NSHR,env,hard_law,theta_KK,&
          &conv_issue_dil)
     ! ---------------------------------------------------------------
  CASE (ONE_int)  ! NCP-function (Kanzow-Kleinmichel)
     ! ---------------------------------------------------------------
     WRITE(16,'(A)')'formulation based on NCP-function is not yet implemented in "stress_computation_single_crystal_plasticity"'
     WRITE(16,'(A)')'STOPPING COMPUTATION'
     STOP
     ! ---------------------------------------------------------------
  CASE (TWO_int)  ! Perzyna viscoplastic formulation
     ! ---------------------------------------------------------------
     WRITE(16,'(A)')'viscoplastic formulation based on Perzyna overstress function is not yet implemented in "stress_computation_single_crystal_plasticity"'
     WRITE(16,'(A)')'STOPPING COMPUTATION'
     STOP
     ! ---------------------------------------------------------------
  CASE (THREE_int)! Ortiz/Miehe viscoplastic formulation
     ! ---------------------------------------------------------------
     CALL determine_incremental_lagrange_multiplier_miehe_viscoplastic(&
          &F_p_inv_iter,plastic_slip_iter,d_delta_lambda_ip1_d_C,&
          &active_slip_system,unimod_C,F_p_inv_n,C_inv,Je_n1,plastic_slip_n,&
          &algo,dt,elast_const,mat_param_inelast,m_SH,symm_class,&
          &orth_basis,orth_proj,NDI,NSHR,env,hard_law,theta_KK,&
          &aver_adj_stp,conv_issue_dil)
     ! WRITE(16,'(A)')'viscoplastic formulation based on Ortiz/Miehe overstress function is not yet implemented in "stress_computation_single_crystal_plasticity"'
     ! WRITE(16,'(A)')'STOPPING COMPUTATION'
     ! STOP
  CASE DEFAULT
     WRITE(16,'(A)')'Unknown algorithm id used in: "stress_computation_single_crystal_plasticity"'
  END SELECT
  !
  ! ensure plastic incompressibility by projection technique
  det_Fp_inv=comp_det_Matrix(F_p_inv_iter)
  F_p_inv_n1=F_p_inv_iter*det_Fp_inv**(-ONE_by_THREE_r)
  !
  ! compute elastic deformation gradient and its inverse
  F_e_n1=MATMUL(F,F_p_inv_n1)
  CALL M33INV (F_e_n1,F_e_inv_n1,OK_FLAG)
  IF (.NOT.(OK_FLAG)) THEN
     WRITE(16,'(A)')'ERROR in: "stress_computation_single_crystal_plasticity (2)"'
     WRITE(16,'(A)')'M33INV encountered singular matrix'
     STOP
  END IF
  !
  ! compute elastic right Cauchy-Green tensor Ce and its inverse for
  ! tangent computation
  Ce_n1_inv=MATMUL(F_e_inv_n1,TRANSPOSE(F_e_inv_n1))
  Ce_n1=MATMUL(TRANSPOSE(F_e_n1),F_e_n1)
  Je_n1=comp_det_Matrix(F_e_n1)
  !
  SELECT CASE (symm_class)
     ! ---------------------------------------------------------------
  CASE (ZERO_int)! volumetric/isochoric decoupled, isotropic elasticity
     ! ---------------------------------------------------------------
     ! assign material parameters
     kappa=elast_const(1); mu=elast_const(2)
     S2PK_IC=S2PK_IC_isotr_elasticity(Ce_n1,Ce_n1_inv,Je_n1,kappa,mu)
     ! ---------------------------------------------------------------
  CASE (ONE_int) ! isotropic elasticity based on Seth-Hill strain
     ! ---------------------------------------------------------------
     CALL S2PK_IC_cubic_elasticity_seth_hill(S2PK_IC_vec,C_mat_tang,Ce_n1,&
          &m_SH,elast_const,symm_class,orth_basis,orth_proj,NDI,NSHR,env)
     !
     CALL transform_symm_tens_vec_form(S2PK_IC,S2PK_IC_vec,ONE_int,&
          &env)
     !
     ! ---------------------------------------------------------------
  CASE (TWO_int) ! cubic elasticity based on Seth-Hill strain
     ! ---------------------------------------------------------------
     CALL S2PK_IC_cubic_elasticity_seth_hill(S2PK_IC_vec,C_mat_tang,Ce_n1,&
          &m_SH,elast_const,symm_class,orth_basis,orth_proj,NDI,NSHR,env)
     CALL transform_symm_tens_vec_form(S2PK_IC,S2PK_IC_vec,ONE_int,&
          &env)
     !
  CASE DEFAULT
     WRITE(16,'(A)')'Unknown symmetry class in "stress_computation_single_crystal_plasticity (2)"'
     WRITE(16,'(A)')'STOPPING COMPUTATION'
     STOP
  END SELECT
  !
  S2PK=MATMUL(F_p_inv_n1,MATMUL(S2PK_IC,transpose(F_p_inv_n1)))
  !
END SUBROUTINE stress_computation_single_crystal_plasticity
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! Compute the stress and tangent modulus for isotropic/cubic elasticity
! quadratic strain energy function employing Seth-Hill strain measures
! **********************************************************************
SUBROUTINE S2PK_IC_cubic_elasticity_seth_hill(S_2PK,C_mat_tang,Ce,m_SH,&
     &elast_const,symm_class,orth_basis,orth_proj,NDI,NSHR,env)
  !
  ! based on: Miehe, Lambrecht, Algorithms of computation of stress and
  !             elasticity moduli in terms of Seth-Hill family of
  !             generalized strain tensors, CINME, 2001
  !           Mahnken, Anisotropy in geometrically non-linear elasticity
  !             with generalized Seth–Hill strain tensors projected to
  !             invariant subspaces, CINME, 2005
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:double_contraction_fourth_order_second_order_tensor_matmul,&
       & double_contraction_fourth_order_tensor_matmul, &
       & transform_symm_tens_vec_form
  USE constitutive_routines_interface, only:check_eigval_distinct,&
       & coefficients_for_projection, diagonal_funct_deriv_projectors, &
       & linear_elasticity_subspace_projection, &
       & projection_tensors_stress_tangent_transformation
  !
  ! declaration of variables
  !
  ! Ce..........elastic right Cauchy-Green tensor
  ! m_SH........Seth-Hill parameter (defining the strain measure)
  ! elast_const.elastic constants (isotr.: kappa,mu;cubic: C11,C12,C44)
  ! symm_class..symmetry class (1=isotropic, 2=cubic)
  ! orth_basis..container of orthogonal basis (second order tensor)
  ! orth_proj...container of dyadic product of orthogonal basis
  ! Ce_vec......vector form of Ce
  ! Ce_princ....principle/eigen values of Ce
  ! E_SH........Seth-Hill strain tensor
  ! ()_vect.....diagonal (generating) functions (stored in vector)
  ! env.........environment in which routine is called
  ! S_2PK.......2nd Piola-Kirchhoff stress (stored in vector)
  ! C_mat_tang..material tangent (stored in matrix format)
  !
  IMPLICIT NONE
  !
  INTEGER(KIND=ik),INTENT(IN)::symm_class,NDI,NSHR
  INTEGER(KIND=ik)::n_distinct,i_count
  INTEGER(KIND=ik),DIMENSION(3)::index_distinct
  !
  REAL(KIND=rk),INTENT(IN)::m_SH
  REAL(KIND=rk),DIMENSION(3),INTENT(IN)::elast_const
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::Ce
  REAL(KIND=rk),DIMENSION(6,3,3),INTENT(IN)::orth_basis
  REAL(KIND=rk),DIMENSION(6,6,6),INTENT(IN)::orth_proj
  REAL(KIND=rk),DIMENSION(6),INTENT(OUT)::S_2PK
  REAL(KIND=rk),DIMENSION(6,6),INTENT(OUT)::C_mat_tang
  REAL(KIND=rk),DIMENSION(3)::Ce_princ
  REAL(KIND=rk),DIMENSION(6)::Ce_vec,T_stress_vec
  REAL(KIND=rk),DIMENSION(3,3)::E_SH,T_stress
  REAL(KIND=rk),DIMENSION(6,6)::tangent,P_tens,T_dc_L
  REAL(KIND=rk),DIMENSION(3)::e_vect,d_vect,f_vect
  REAL(KIND=rk),DIMENSION(3,3,3)::Proj_tens
  REAL(KIND=rk)::eta
  REAL(KIND=rk),DIMENSION(3,3)::theta,xi
  !
  CHARACTER(4),INTENT(IN)::env
  !
  ! transform elastic right Cauchy-Green to vector format and compute
  ! the eigenvalues
  CALL transform_symm_tens_vec_form(Ce,Ce_vec,ZERO_int,env)
  CALL SPRINC(Ce_vec,Ce_princ,ONE_int,NDI,NSHR)
  !
  ! determine the number of distinct eigenvalues
  CALL check_eigval_distinct(n_distinct,index_distinct,Ce_princ,&
       &zero_distinct)
  !
  ! obtain diagonal functions and their derivatives and projectors
  CALL diagonal_funct_deriv_projectors(e_vect,d_vect,f_vect,Proj_tens,&
       &Ce_princ,index_distinct,n_distinct,Ce,m_SH,zero_distinct)
  !
  ! compute Seth-Hill strain tensor
  E_SH=ZERO_r
  DO i_count=ONE_int,n_distinct
     E_SH = E_SH + e_vect(i_count)*Proj_tens(i_count,:,:)
  END DO
  !
  ! evaluate the elastic law (dual stress and elasticity moduli)
  CALL linear_elasticity_subspace_projection(T_stress,tangent,E_SH,&
       &elast_const,symm_class,orth_basis,orth_proj)
  !
  CALL transform_symm_tens_vec_form(T_stress,T_stress_vec,ZERO_int,&
       &env)
  !
  ! compute coefficients for projection to Lagrangian stress / moduli
  CALL coefficients_for_projection(theta,xi,eta,e_vect,d_vect,Ce_princ,&
       &index_distinct,n_distinct)
  !
  ! compute projections
  CALL projection_tensors_stress_tangent_transformation(P_tens,T_dc_L,&
       &Proj_tens,d_vect,f_vect,eta,xi,theta,n_distinct,T_stress,env)
  !
  ! compute 2nd Piola Kirchhoff stress and corresponding moduli
  S_2PK=double_contraction_fourth_order_second_order_tensor_matmul(&
       &P_tens,T_stress_vec)
  !
  C_mat_tang=double_contraction_fourth_order_tensor_matmul(P_tens,&
       &double_contraction_fourth_order_tensor_matmul(tangent,P_tens))+&
       &T_dc_L
  !
END SUBROUTINE S2PK_IC_cubic_elasticity_seth_hill
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! update the internal variables in iterative procedure
! **********************************************************************
SUBROUTINE update_internal_variable_nw_iter(delta_lambda_ip1,&
     &plastic_slip_iter,F_p_inv_iter,dd_delta_lambda_ip1,&
     &plastic_slip_n,F_p_inv_n)
  !
  USE Abaqus_Interface
  USE Constants
  USE constitutive_routines_interface, only:update_plastic_deformation_gradient
  !
  ! variable declaration
  ! F_p_inv_()..........inverse of plastic part of deformation gradient
  ! delta_lambda_ip1....incremental lagrange multiplier
  ! plastic_slip_().....accumulated plastic slip on each slip system
  ! dd_delta_lambda_ip1.iterative improvement of incr. lagr. multipl.
  ! 
  IMPLICIT NONE
  REAL(KIND=rk),DIMENSION(24),INTENT(INOUT)::delta_lambda_ip1
  REAL(KIND=rk),DIMENSION(24),INTENT(OUT)::plastic_slip_iter
  REAL(KIND=rk),DIMENSION(3,3),INTENT(OUT)::F_p_inv_iter
  REAL(KIND=rk),DIMENSION(24),INTENT(IN)::dd_delta_lambda_ip1
  REAL(KIND=rk),DIMENSION(24),INTENT(IN)::plastic_slip_n
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::F_p_inv_n
  !
  ! update increment in Lagrange multiplier
  delta_lambda_ip1=delta_lambda_ip1+dd_delta_lambda_ip1
  !
  ! compute hardening variables on all slip systems
  plastic_slip_iter=plastic_slip_n+delta_lambda_ip1
  !
  ! update plastic part of the deformation gradient
  F_p_inv_iter=update_plastic_deformation_gradient(F_p_inv_n,&
       &delta_lambda_ip1)
  !
END SUBROUTINE update_internal_variable_nw_iter
! ......................................................................
! **********************************************************************
!
!#######################################################################
!
! **********************************************************************
! auxiliary function to compute product of projectors (required for
! elastic law based on Seth-Hill strain measure)
! **********************************************************************
FUNCTION aux_prod_projector(A_tens,lambda,index_distinct,p_app,q_app)
  !
  USE Abaqus_Interface
  USE Constants
  !
  ! declaration of variables
  !
  ! lambda...........vector of eigenvalues
  ! index_distinct...index of the distinct eigenvalues in lambda
  ! A_tens...........symmetric, second order tensor whose eigenvalues
  !                  are stored in lambda
  !
  IMPLICIT NONE
  !
  INTEGER(KIND=ik),INTENT(IN)::p_app,q_app
  INTEGER(KIND=ik),DIMENSION(3),INTENT(IN)::index_distinct
  !
  REAL(KIND=rk),DIMENSION(3),INTENT(IN)::lambda
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::A_tens
  REAL(KIND=rk),DIMENSION(3,3)::aux_prod_projector
  !
  aux_prod_projector=(A_tens-lambda(index_distinct(q_app))*IDENT)/ &
       & (lambda(index_distinct(p_app))-lambda(index_distinct(q_app)))
  !
END FUNCTION aux_prod_projector
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! compute additional contribution to consistent Jacobian stemming from
! the transformation of the Lie-derivative to the Jaumann-rate (i.e. a
! constitutive model, which employs Lie-derivative in its rate form is
! transformed to be consistent with Abaqus' Jaumann-rate of the Kirch-
! hoff stress)
! **********************************************************************
FUNCTION comp_add_stress_contrib_tang(STRESS)
  !
  USE Abaqus_Interface
  USE Constants
  !
  ! NOTE: storage order in Abaqus: 11,22,33,12,13,23
  !
  ! declaration of variables
  ! STRESS........Cauchy-stress tensor (vector form)
  !
  IMPLICIT NONE
  REAL(KIND=rk),DIMENSION(6),INTENT(IN)::STRESS
  REAL(KIND=rk),DIMENSION(6,6)::add_stress
  REAL(KIND=rk),DIMENSION(6,6)::comp_add_stress_contrib_tang
  !
  add_stress=RESHAPE((/TWO_r*STRESS(1),ZERO_r,ZERO_r,STRESS(4),&
       &STRESS(5),ZERO_r,&
       &ZERO_r,TWO_r*STRESS(2),ZERO_r,STRESS(4),ZERO_r,STRESS(6),&
       &ZERO_r,ZERO_r,TWO_r*STRESS(3),ZERO_r,STRESS(5),STRESS(6),&
       &STRESS(4),STRESS(4),ZERO_r,ONE_HALF_r*(STRESS(1)+STRESS(2)),&
       &ONE_HALF_r*STRESS(6),ONE_HALF_r*STRESS(5),&
       &STRESS(5),ZERO_r,STRESS(5),ONE_HALF_r*STRESS(6),ONE_HALF_r*&
       &(STRESS(1)+STRESS(3)),ONE_HALF_r*STRESS(4),&
       &ZERO_r,STRESS(6),STRESS(6),ONE_HALF_r*STRESS(5),ONE_HALF_r*&
       &STRESS(4),ONE_HALF_r*(STRESS(2)+STRESS(3))/),(/6,6/))
  !
  comp_add_stress_contrib_tang=add_stress
  !
END FUNCTION comp_add_stress_contrib_tang
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! computes the derivative of the incremental lagrange multipliers w.r.t.
! right Cauchy-Green tensor C
! **********************************************************************
FUNCTION comp_d_delta_lambda_ip1_d_C (active_slip_system,jacobian_nw,&
     &d_delta_lambda_i_d_C,F_p_inv_iter,unimod_C,Je,C_inv,elast_const,&
     &m_SH,symm_class,orth_basis,orth_proj,eta,dt,algo,&
     &deriv_visco_plast_potent,env,NDI,NSHR)
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:comp_sec_ord_tens_rank_one, &
       & double_contraction_fourth_order_tensor, &
       & DYADIC_PRODUCT_SECOND_ORDER_TENSORS, &
       & expand_fourth_order_tensor_matrix2tensor,&
       & left_double_contraction_fourth_order_tensor_second_order_tensor,&
       & SOLVE_SYSTEM_OF_EQUATION, &
       & SOLVE_ILL_COND_SYSTEM_OF_EQUATION, &
       & sq_prod_unsym_sec_ord_tens_dble_contr_fourth_order_sym_Ident,&
       & transform_symm_tens_vec_form
  USE constitutive_routines_interface, only:derivative_mandel_stress_ce_sh_elasticity,&
       & slip_system_fcc, S2PK_IC_cubic_elasticity_seth_hill
  !
  ! declaration of variables
  ! active_slip_system..........list of active slip systems
  ! jacobian_nw.................jacobian matrix
  ! d_delta_lambda_i_d_C........derivative of incremental Lagrange 
  !                             multipliers w.r.t. right Cauchy-Green 
  !                             tensor (at fixed-point iteration i, augmented 
  !                             lagrangian approach)
  ! F_p_inv_iter................inverse of plastic deformation gradient
  ! unimod_C....................unimodular right Cauchy-Green (truely 
  !                             unimodular in case of isochoric-volumetric
  !                             elastic model
  ! Je..........................determinant of elastic deformation grad
  ! C_inv.......................inverse of right Cauchy-Green tensor
  ! elast_const.................elastic constants
  ! m_SH........................Seth-Hill parameter
  ! symm_class..................symmetry class of crystal (0=isotropic 
  !                             isochoric/volumetric split, 1=isotropic, 
  !                             2=cubic)
  ! orth_basis..................all second order, orthogonal basis N(i),
  !                             rotated in case of cubic elasticity
  ! orth_proj...................symmetric dyadic products of orthogonal 
  !                             basis i.e. N(i) dyadic N(i)
  ! eta.........................pseudo-viscosity in augmented lagrangian
  ! dt..........................time increment
  ! algo........................id for algorithm employed
  !                             0...augmented-lagrangian
  !                             1...NCP-function (Kanzow-Kleinmichel)
  !                             2...Perzyna viscoplastic formulation
  !                             3...Ortiz/Miehe viscoplastic formulation
  ! deriv_visco_plast_potent....derivative of the viscoplast potential
  !                             aux. var in rate-independent case
  ! d_tau_alpha_d_C.............derivative of projected Mandel stress 
  !                             w.r.t. Ce and Ce w.r.t. C
  !
  IMPLICIT NONE
  REAL(KIND=rk),INTENT(IN)::Je,m_SH,eta,dt
  REAL(KIND=rk),DIMENSION(3),INTENT(IN)::elast_const
  REAL(KIND=rk),DIMENSION(24),INTENT(IN)::deriv_visco_plast_potent
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::F_p_inv_iter,unimod_C
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::C_inv
  REAL(KIND=rk),DIMENSION(24,6),INTENT(IN)::d_delta_lambda_i_d_C
  REAL(KIND=rk),DIMENSION(24,24),INTENT(IN)::jacobian_nw
  REAL(KIND=rk),DIMENSION(6,3,3),INTENT(IN)::orth_basis
  REAL(KIND=rk),DIMENSION(6,6,6),INTENT(IN)::orth_proj
  REAL(KIND=rk),DIMENSION(24,6)::comp_d_delta_lambda_ip1_d_C
  INTEGER(KIND=ik),INTENT(IN)::algo,symm_class,NDI,NSHR
  REAL(KIND=rk)::mu
  REAL(KIND=rk),DIMENSION(3)::normal_dir,slip_dir
  REAL(KIND=rk),DIMENSION(6)::S_2PK_IC_vec,d_tau_alpha_d_C_vec
  REAL(KIND=rk),DIMENSION(3,3)::unimod_Ce_iter,N_alpha
  REAL(KIND=rk),DIMENSION(3,3)::d_tau_alpha_d_C,S_2PK_IC
  REAL(KIND=rk),DIMENSION(6,6)::C_mat_tang_IC
  REAL(KIND=rk),DIMENSION(24,6)::rhs_d_delta_lambda_ip1_d_C
  REAL(KIND=rk),DIMENSION(24,6)::d_delta_lambda_ip1_d_C
  REAL(KIND=rk),DIMENSION(9,9)::sq_prod_F_p_inv_tr
  REAL(KIND=rk),DIMENSION(3,3,3,3)::sq_prod_F_p_inv_tr_tens
  REAL(KIND=rk),DIMENSION(3,3,3,3)::dyad_prod_unimodCe_Cinv
  REAL(KIND=rk),DIMENSION(3,3,3,3)::d_unimod_Ce_d_C,d_dev_Mandel_d_C
  REAL(KIND=rk),DIMENSION(3,3,3,3)::d_Mandel_stress_d_Ce
  INTEGER(KIND=ik)::i_slip
  !
  CHARACTER(4),INTENT(IN)::env
  LOGICAL,DIMENSION(24),INTENT(IN)::active_slip_system
  !
  ! elastic right Cauchy-Green tensor (unimodular in case of vol/isochoric
  ! decoupled isotropic elasticity)
  unimod_Ce_iter=MATMUL(TRANSPOSE(F_p_inv_iter),MATMUL(unimod_C,F_p_inv_iter))
  rhs_d_delta_lambda_ip1_d_C=ZERO_r
  !
  ! derivative of unimodular elastic right Cauchy-Green tensor Ce w.r.t.
  ! right Cauchy-Green tensor C
  sq_prod_F_p_inv_tr = &
       & sq_prod_unsym_sec_ord_tens_dble_contr_fourth_order_sym_Ident(&
       & TRANSPOSE(F_p_inv_iter),TRANSPOSE(F_p_inv_iter),env)
  sq_prod_F_p_inv_tr_tens = expand_fourth_order_tensor_matrix2tensor(&
       &sq_prod_F_p_inv_tr,env)
  !
  SELECT CASE (symm_class)
     ! ---------------------------------------------------------------
  CASE (ZERO_int)! volumetric/isochoric decoupled, isotropic elasticity
     ! ---------------------------------------------------------------
     dyad_prod_unimodCe_Cinv = DYADIC_PRODUCT_SECOND_ORDER_TENSORS(&
          &unimod_Ce_iter,C_inv)
     d_unimod_Ce_d_C=-ONE_by_THREE_r*dyad_prod_unimodCe_Cinv+&
          &sq_prod_F_p_inv_tr_tens*Je**(-TWO_by_THREE_r)
     !
     ! derivative of Mandel stress w.r.t. right Cauchy-Green tensor C
     mu=elast_const(TWO_int)
     d_dev_Mandel_d_C=mu*d_unimod_Ce_d_C
     !
     ! ---------------------------------------------------------------
  CASE (ONE_int) ! isotropic elasticity based on Seth-Hill strain
     ! ---------------------------------------------------------------
     CALL S2PK_IC_cubic_elasticity_seth_hill(S_2PK_IC_vec,C_mat_tang_IC,&
          &unimod_Ce_iter,m_SH,elast_const,symm_class,orth_basis,&
          &orth_proj,NDI,NSHR,env)
     CALL transform_symm_tens_vec_form(S_2PK_IC,S_2PK_IC_vec,ONE_int,&
          &env)
     !
     ! compute the derivative of the Mandel stress w.r.t. elastic right
     ! Cauchy-Green tensor Ce
     d_Mandel_stress_d_Ce = derivative_mandel_stress_ce_sh_elasticity(&
          &S_2PK_IC,C_mat_tang_IC,unimod_Ce_iter,env)
     d_dev_Mandel_d_C=double_contraction_fourth_order_tensor(&
          &d_Mandel_stress_d_Ce,sq_prod_F_p_inv_tr_tens)
     !
     ! ---------------------------------------------------------------
  CASE (TWO_int) ! cubic elasticity based on Seth-Hill strain
     ! ---------------------------------------------------------------
     CALL S2PK_IC_cubic_elasticity_seth_hill(S_2PK_IC_vec,C_mat_tang_IC,&
          &unimod_Ce_iter,m_SH,elast_const,symm_class,orth_basis,&
          &orth_proj,NDI,NSHR,env)
     CALL transform_symm_tens_vec_form(S_2PK_IC,S_2PK_IC_vec,ONE_int,&
          &env)
     !
     ! compute the derivative of the Mandel stress w.r.t. elastic right
     ! Cauchy-Green tensor Ce
     d_Mandel_stress_d_Ce = derivative_mandel_stress_ce_sh_elasticity(&
          &S_2PK_IC,C_mat_tang_IC,unimod_Ce_iter,env)
     d_dev_Mandel_d_C=double_contraction_fourth_order_tensor(&
          &d_Mandel_stress_d_Ce,sq_prod_F_p_inv_tr_tens)
     !
  CASE DEFAULT
     WRITE(16,'(2A)')'Unknown symm_class in: comp_d_delta_lambda_ip1_d_C'
     WRITE(16,'(A)')'STOPPING COMPUTATION'
     STOP
  END SELECT
  ! compute the derivative of the resolved shear stress w.r.t. right 
  ! Cauchy-Green tensor C
  DO i_slip=ONE_int,TWENTY_FOUR_int
     IF (active_slip_system(i_slip)) THEN
        ! tensor defining the slip system
        CALL slip_system_fcc(normal_dir,slip_dir,i_slip)
        N_alpha=comp_sec_ord_tens_rank_one(slip_dir,normal_dir)
        !
        ! derivative of resolved shear stress w.r.t. C
        d_tau_alpha_d_C=left_double_contraction_fourth_order_tensor_second_order_tensor(&
             &d_dev_Mandel_d_C,N_alpha)
        CALL transform_symm_tens_vec_form(d_tau_alpha_d_C,&
             &d_tau_alpha_d_C_vec,ZERO_int,env)
        !
        ! the factor 2.0 in front of dPhi_alpha_dC is due to the fact 
        ! that material tangent is contracted with 0.5*\Delta C
        SELECT CASE (algo)
        CASE (ZERO_int) ! augmented lagrangian approach
           rhs_d_delta_lambda_ip1_d_C(i_slip,:) = d_delta_lambda_i_d_C(i_slip,:)+&
                &TWO_r*eta*dt*d_tau_alpha_d_C_vec
        CASE (ONE_int)  ! NCP-function (Kanzow-Kleinmichel)
           rhs_d_delta_lambda_ip1_d_C(i_slip,:) =-TWO_r*&
                &deriv_visco_plast_potent(i_slip)*d_tau_alpha_d_C_vec
        CASE (TWO_int)  ! Perzyna viscoplastic formulation
           rhs_d_delta_lambda_ip1_d_C(i_slip,:) = TWO_r*&
                &deriv_visco_plast_potent(i_slip)*d_tau_alpha_d_C_vec
        CASE (THREE_int)! Ortiz/Miehe viscoplastic formulation
           rhs_d_delta_lambda_ip1_d_C(i_slip,:) = TWO_r*&
                &deriv_visco_plast_potent(i_slip)*d_tau_alpha_d_C_vec
        CASE DEFAULT
           WRITE(16,'(A)')'Unknown algorithm id used in: comp_d_delta_lambda_ip1_d_C'
        END SELECT
     END IF
  END DO
  !
  SELECT CASE (algo)
  CASE (ZERO_int) ! augmented lagrangian approach
     d_delta_lambda_ip1_d_C=SOLVE_SYSTEM_OF_EQUATION(jacobian_nw,&
          &rhs_d_delta_lambda_ip1_d_C,SIX_int,TWENTY_FOUR_int)
  CASE (ONE_int)  ! NCP-function (Kanzow-Kleinmichel)
     d_delta_lambda_ip1_d_C=SOLVE_ILL_COND_SYSTEM_OF_EQUATION(jacobian_nw,&
          &rhs_d_delta_lambda_ip1_d_C,SIX_int,TWENTY_FOUR_int,-eps_iter)
  CASE (TWO_int)  ! Perzyna viscoplastic formulation
     d_delta_lambda_ip1_d_C=SOLVE_ILL_COND_SYSTEM_OF_EQUATION(jacobian_nw,&
          &rhs_d_delta_lambda_ip1_d_C,SIX_int,TWENTY_FOUR_int,-eps_iter)
  CASE (THREE_int)! Ortiz/Miehe viscoplastic formulation
     d_delta_lambda_ip1_d_C=SOLVE_ILL_COND_SYSTEM_OF_EQUATION(jacobian_nw,&
          &rhs_d_delta_lambda_ip1_d_C,SIX_int,TWENTY_FOUR_int,-eps_iter)
  CASE DEFAULT
     WRITE(16,'(A)')'Unknown algorithm id used in: comp_d_delta_lambda_ip1_d_C (sec)'
  END SELECT
  comp_d_delta_lambda_ip1_d_C = d_delta_lambda_ip1_d_C
END FUNCTION comp_d_delta_lambda_ip1_d_C
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! compute the derivative of the elastic right Cauchy-Green tensor w.r.t.
! the inverse of the plastic deformation gradient
! Ce=F^(p-T)*C*F^(p-1) -> d_Ce_d_Fp_inv
! **********************************************************************
FUNCTION comp_deriv_d_C_e_d_F_p_inv(F_p_inv,C,env)
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:expand_fourth_order_tensor_matrix2tensor,&
       & spec_tensor_product_unsym_second_order_tensors, &
       & square_product_unsym_second_order_tensors
  !
  ! variable declaration
  ! F_p_inv.........inverse of plastic part of deformation gradient
  ! C...............right Cauchy-Green tensor
  !
  IMPLICIT NONE
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::F_p_inv,C
  REAL(KIND=rk),DIMENSION(3,3)::temp_tens
  REAL(KIND=rk),DIMENSION(3,3,3,3)::d_Ce_d_Fp_inv
  REAL(KIND=rk),DIMENSION(3,3,3,3)::comp_deriv_d_C_e_d_F_p_inv
  !
  CHARACTER(4),INTENT(IN)::env
  !
  temp_tens=MATMUL(TRANSPOSE(F_p_inv),C)
  d_Ce_d_Fp_inv=expand_fourth_order_tensor_matrix2tensor(&
       &spec_tensor_product_unsym_second_order_tensors(IDENT,temp_tens,&
       &env),env)
  d_Ce_d_Fp_inv=d_Ce_d_Fp_inv+expand_fourth_order_tensor_matrix2tensor(&
       &square_product_unsym_second_order_tensors(temp_tens,IDENT,env),&
       &env)
  !
  comp_deriv_d_C_e_d_F_p_inv=d_Ce_d_Fp_inv
  !
END FUNCTION comp_deriv_d_C_e_d_F_p_inv
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! compute the derivative of the second Piola-Kirchhoff stress in the 
! intermediate configuration w.r.t. elastic right Cauchy-Green tensor 
! C^e
! (ATTENTION: isotropic elasticity with volumetric-isochoric split is 
! assumed)
! tangent computation due to Miehe, Lambrecht, Algorithms for computation
! of stresses and elasticity moduli in terms of Seth–Hill’s family of 
! generalized strain tensors, CINME, 2001, p.337-353
! **********************************************************************
FUNCTION comp_derivative_d_S2PK_IC_d_Ce_IC (Je,Ce,Ce_inv,kappa,mu,env)
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:compute_trace_second_order_tensor,&
       & DYADIC_PRODUCT_SECOND_ORDER_TENSORS, &
       & expand_minor_symm_fourth_order_tensor_matrix2tensor, &
       & symm_tensor_product_second_order_tensors
  !
  ! variable declaration
  ! Je.........determinant of elastic part of deformation gradient
  ! Ce.........elastic right Cauchy-Green tensor
  ! Ce_inv.....inverse of elastic right Cauchy-Green tensor
  ! kappa......bulk modulus
  ! mu.........shear modulus
  !
  IMPLICIT NONE
  REAL(KIND=rk),INTENT(IN)::Je,kappa,mu
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::Ce,Ce_inv
  REAL(KIND=rk)::tr_Ce
  REAL(KIND=rk),DIMENSION(3,3)::Dev_S2PK_IC_bar
  REAL(KIND=rk),DIMENSION(3,3,3,3)::Ce_inv_dyad_Ce_inv,C_iso_4ot
  REAL(KIND=rk),DIMENSION(3,3,3,3)::Ident_Ce_inv_tens,Proj_Ce_inv
  REAL(KIND=rk),DIMENSION(3,3,3,3)::d_S2PK_IC_d_Ce
  REAL(KIND=rk),DIMENSION(3,3,3,3)::comp_derivative_d_S2PK_IC_d_Ce_IC
  REAL(KIND=rk),DIMENSION(6,6)::Ident_Ce_inv
  !
  CHARACTER(4),INTENT(IN)::env
  !
  Ce_inv_dyad_Ce_inv=DYADIC_PRODUCT_SECOND_ORDER_TENSORS(Ce_inv,Ce_inv)
  Ident_Ce_inv=symm_tensor_product_second_order_tensors(Ce_inv,Ce_inv,&
       &env)
  Ident_Ce_inv=ONE_HALF_r*Ident_Ce_inv
  Ident_Ce_inv_tens=expand_minor_symm_fourth_order_tensor_matrix2tensor(&
       &Ident_Ce_inv,env)
  Proj_Ce_inv=Je**(-TWO_by_THREE_r)*(Ident_Ce_inv_tens-ONE_by_THREE_r*&
       &Ce_inv_dyad_Ce_inv)
  tr_Ce=compute_trace_second_order_tensor(Ce)
  Dev_S2PK_IC_bar=Je**(-TWO_by_THREE_r)*mu*(IDENT-ONE_by_THREE_r*tr_Ce*&
       &Ce_inv)
  C_iso_4ot=TWO_by_THREE_r*mu*tr_Ce*Proj_Ce_inv-TWO_by_THREE_r*&
       &(dyadic_product_second_order_tensors(Dev_S2PK_IC_bar,Ce_inv)+&
       &dyadic_product_second_order_tensors(Ce_inv,Dev_S2PK_IC_bar))
  d_S2PK_IC_d_Ce=kappa*Ce_inv_dyad_Ce_inv-TWO_r*kappa*log(Je)*&
       &Ident_Ce_inv_tens+C_iso_4ot
  comp_derivative_d_S2PK_IC_d_Ce_IC=d_S2PK_IC_d_Ce
  !
END FUNCTION comp_derivative_d_S2PK_IC_d_Ce_IC
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! compute the consistent algorithmic tangent for single crystal plasticity
! model (e.g. Schmidt-Baldassari, CMAME, 192 (2003), pp.1261-1280)
! for both the rate-independent formulation based on augmented Lagrangian
! formulation, NCP-function based formulation and the rate-dependent
! formulation based on Perzyna-type or Miehe-Ortiz-typeviscoplastic
! formulation (can be interpreted as penalty formulation of constraint
! optimization problem)
! **********************************************************************
FUNCTION consistent_algorithmic_tangent_augmented_lagrangian_penalty(C,&
     &F_p_inv_n,active_slip_system,d_delta_lambda_ip1_d_C,det_Fp_inv,&
     &F_p_inv_iter,Je_n1,Ce_n1,Ce_n1_inv,F_p_inv_n1,S2PK_IC,elast_const,&
     &m_SH,symm_class,orth_basis,orth_proj,NDI,NSHR,env)
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:double_contraction_fourth_order_tensor,&
       & expand_fourth_order_tensor_matrix2tensor, &
       & expand_minor_symm_fourth_order_tensor_matrix2tensor,&
       & extract_symm_fourth_ord_tens_minor_symm, &
       & map_fourth_order_tensor_2_matrix_format, &
       & sq_prod_unsym_sec_ord_tens_dble_contr_fourth_order_sym_Ident, &
       & square_product_unsym_second_order_tensors
  
  USE constitutive_routines_interface, only:comp_deriv_d_C_e_d_F_p_inv, &
       & comp_derivative_d_S2PK_IC_d_Ce_IC, &
       & deriv_d_S_d_Fp_inv, &
       & derivative_inverse_plastic_def_grad_right_cauchy_green_tensor_c, &
       & mod_d_F_p_inv_d_C_isochoric_proj_constraint,&
       & S2PK_IC_cubic_elasticity_seth_hill
  !
  ! variable declaration
  ! C.......................right Cauchy-Green tensor
  ! F_p_inv_n...............inverse of plastic part of deformation
  !                         gradient at t_n
  ! active_slip_system......list of active slip systems
  ! d_delta_lambda_ip1_d_C..deriv. of incremental lagrange multipliers
  !                         w.r.t. right Cauchy-Green tensor
  ! det_Fp_inv..............determinant of inverse of plastic part of
  !                         deformation gradient
  ! F_p_inv_iter............inverse of plastic part of deform. grad. (iter)
  ! Je_n1...................determinant of elastic part of deform. grad.
  ! Ce_n1...................elastic right Cauchy-Green tensor
  ! Ce_n1_inv...............inverse of elastic right Cauchy-Green tensor
  ! F_p_inv_n1..............inverse of plastic part of deform. grad.
  !                         at t_{n+1}
  ! S2PK_IC.................second Piola-Kirchhoff stress in intermed. conf.
  ! env.....................environment in which routine is called
  !
  IMPLICIT NONE
  REAL(KIND=rk),INTENT(IN)::det_Fp_inv,Je_n1,m_SH
  REAL(KIND=rk),DIMENSION(3),INTENT(IN)::elast_const
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::C,F_p_inv_n,F_p_inv_iter
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::Ce_n1,Ce_n1_inv
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::S2PK_IC,F_p_inv_n1
  REAL(KIND=rk),DIMENSION(24,6),INTENT(IN)::d_delta_lambda_ip1_d_C
  REAL(KIND=rk),DIMENSION(6,3,3),INTENT(IN)::orth_basis
  REAL(KIND=rk),DIMENSION(6,6,6),INTENT(IN)::orth_proj
  REAL(KIND=rk)::kappa,mu
  REAL(KIND=rk),DIMENSION(6)::S_2PK_vec
  REAL(KIND=rk),DIMENSION(6,6)::C_mat_tang,d_S_d_C
  REAL(KIND=rk),DIMENSION(6,6)::d_S_d_Fp_inv_dc_d_Fp_inv_d_C
  REAL(KIND=rk),DIMENSION(6,6)::consistent_algorithmic_tangent_augmented_lagrangian_penalty
  REAL(KIND=rk),DIMENSION(9,9)::d_Ce_d_C,sq_prod_F_p_inv
  REAL(KIND=rk),DIMENSION(3,3,3,3)::d_F_p_inv_n1_d_C,d_S2PK_IC_d_Ce
  REAL(KIND=rk),DIMENSION(3,3,3,3)::d_Ce_d_C_tens,sq_prod_F_p_inv_tens
  REAL(KIND=rk),DIMENSION(3,3,3,3)::temp_tens,d_S_d_C_tens
  REAL(KIND=rk),DIMENSION(3,3,3,3)::d_S_d_Fp_inv_tens
  REAL(KIND=rk),DIMENSION(3,3,3,3)::d_Ce_d_Fp_inv,plast_tang_contrib
  INTEGER(KIND=ik),INTENT(IN)::symm_class,NDI,NSHR
  LOGICAL,DIMENSION(24),INTENT(IN)::active_slip_system
  !
  CHARACTER(4),INTENT(IN)::env
  !
  ! derivative of inverse of plastic deformation gradient w.r.t. right
  ! Cauchy-Green tensor C
  d_F_p_inv_n1_d_C=&
       &derivative_inverse_plastic_def_grad_right_cauchy_green_tensor_c(&
       &F_p_inv_n,active_slip_system,d_delta_lambda_ip1_d_C,env)
  !
  ! modification of the above derivative d_F_p_inv_n1_d_C due to
  ! linearization of isochoric projection technique
  d_F_p_inv_n1_d_C=mod_d_F_p_inv_d_C_isochoric_proj_constraint(&
       &d_F_p_inv_n1_d_C,det_Fp_inv,F_p_inv_iter)
  !
  SELECT CASE (symm_class)
  CASE (ZERO_int) ! isotropic elasticity (volumetric-isochoric split)
     kappa=elast_const(1); mu=elast_const(2)
     ! derivative of the second Piola-Kirchhoff in intermediate
     ! configuration w.r.t. elastic right Cauchy-Green tensor C^e
     d_S2PK_IC_d_Ce=comp_derivative_d_S2PK_IC_d_Ce_IC(Je_n1,Ce_n1,&
          &Ce_n1_inv,kappa,mu,env)
     !
  CASE (ONE_int)  ! isotropic elasticity based on Seth-Hill strain
     ! derivative of the second Piola-Kirchhoff in intermediate
     ! configuration w.r.t. elastic right Cauchy-Green tensor C^e
     CALL S2PK_IC_cubic_elasticity_seth_hill(S_2PK_vec,C_mat_tang,Ce_n1,&
          &m_SH,elast_const,symm_class,orth_basis,orth_proj,NDI,NSHR,&
          &env)
     d_S2PK_IC_d_Ce = &
          &expand_minor_symm_fourth_order_tensor_matrix2tensor(&
          &C_mat_tang,env)
     !
  CASE (TWO_int)  ! cubic elasticity based on Seth-Hill strain
     ! derivative of the second Piola-Kirchhoff in intermediate
     ! configuration w.r.t. elastic right Cauchy-Green tensor C^e
     CALL S2PK_IC_cubic_elasticity_seth_hill(S_2PK_vec,C_mat_tang,Ce_n1,&
          &m_SH,elast_const,symm_class,orth_basis,orth_proj,NDI,NSHR,&
          &env)
     d_S2PK_IC_d_Ce = &
          &expand_minor_symm_fourth_order_tensor_matrix2tensor(&
          &C_mat_tang,env)
     !
  CASE DEFAULT
     WRITE(16,'(A)')'Unknown symmetry class in "consistent_algorithmic_tangent_augmented_lagrangian_penalty"'
     WRITE(16,'(A)')'STOPPING COMPUTATION'
     STOP
  END SELECT
  !
  ! derivative of elastic right Cauchy-Green tensor C^e w.r.t. right
  ! Cauchy-Green tensor C
  d_Ce_d_C=sq_prod_unsym_sec_ord_tens_dble_contr_fourth_order_sym_Ident(&
       &TRANSPOSE(F_p_inv_n1),TRANSPOSE(F_p_inv_n1),env)
  d_Ce_d_C_tens=expand_fourth_order_tensor_matrix2tensor(d_Ce_d_C,env)
  !
  sq_prod_F_p_inv=square_product_unsym_second_order_tensors(F_p_inv_n1,&
       &F_p_inv_n1,env)
  sq_prod_F_p_inv_tens=expand_fourth_order_tensor_matrix2tensor(&
       &sq_prod_F_p_inv,env)
  !
  ! derivative of second Piola-Kirchhoff tensor w.r.t. right Cauchy-Green
  ! tensor
  temp_tens=double_contraction_fourth_order_tensor(d_S2PK_IC_d_Ce,&
       &d_Ce_d_C_tens)
  d_S_d_C_tens=&
       &double_contraction_fourth_order_tensor(sq_prod_F_p_inv_tens,&
       &temp_tens)
  d_S_d_C=extract_symm_fourth_ord_tens_minor_symm(&
       &map_fourth_order_tensor_2_matrix_format(d_S_d_C_tens,env))
  !
  ! derivative of the elastic right Cauchy-Green tensor C^e w.r.t. the
  ! inverse of the plastic deformation gradient
  d_Ce_d_Fp_inv=comp_deriv_d_C_e_d_F_p_inv(F_p_inv_n1,C,env)
  !
  ! derivative of second Piola-Kirchhoff tensor w.r.t. inverse of plastic
  ! deformation gradient (the factor 0.5 is due the fact that the elastic
  ! tangent stored in
  ! d_S2PK_IC_d_Ce = 2.0 * \frac{\partial \hat{S}}{\partial C^{e}})
  temp_tens=double_contraction_fourth_order_tensor(&
       &ONE_HALF_r*d_S2PK_IC_d_Ce,d_Ce_d_Fp_inv)
  d_S_d_Fp_inv_tens=double_contraction_fourth_order_tensor(&
       &sq_prod_F_p_inv_tens,temp_tens)
  d_S_d_Fp_inv_tens=d_S_d_Fp_inv_tens + deriv_d_S_d_Fp_inv(S2PK_IC,&
       &F_p_inv_n1,env)
  !
  ! contraction of derivative of second Piola-Kirchhoff w.r.t. inverse
  ! of deformation gradient with derivative of inverse of deformation
  ! gradient w.r.t. right Cauchy-Green tensor C
  plast_tang_contrib=double_contraction_fourth_order_tensor(&
       &d_S_d_Fp_inv_tens,d_F_p_inv_n1_d_C)
  d_S_d_Fp_inv_dc_d_Fp_inv_d_C=extract_symm_fourth_ord_tens_minor_symm(&
       &map_fourth_order_tensor_2_matrix_format(plast_tang_contrib,env))
  !
  ! add contributions of tangent
  consistent_algorithmic_tangent_augmented_lagrangian_penalty=d_S_d_C + &
       &d_S_d_Fp_inv_dc_d_Fp_inv_d_C
  !
END FUNCTION consistent_algorithmic_tangent_augmented_lagrangian_penalty
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! compute the material tangent, i.e. derivative of 2nd Piola-Kirchhoff
! stress tensor S w.r.t. right Cauchy-Green tensor C for single crystal
! (visco)-plasticity material model
! **********************************************************************
FUNCTION CSDA_material_tangent_single_crystal_plast(S2PK,F,dt,F_p_inv_n,&
     &plastic_slip_n,elast_const,mat_param_inelast,m_SH,symm_class,&
     &orth_basis,orth_proj,numtang,algo,NDI,NSHR,env,hard_law,theta_KK,&
     &aver_adj_stp)
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:approp_index,&
       & pertub_F_numtang_material_tang, M33INV, &
       & transform_symm_tens_vec_form
  USE constitutive_routines_interface, only:stress_computation_single_crystal_plasticity
  !
  ! variable declaration
  IMPLICIT NONE
  REAL(KIND=rk),INTENT(IN)::dt,m_SH,theta_KK,aver_adj_stp
  REAL(KIND=rk),DIMENSION(3),INTENT(IN)::elast_const
  REAL(KIND=rk),DIMENSION(13),INTENT(IN)::mat_param_inelast
  REAL(KIND=rk),DIMENSION(24),INTENT(IN)::plastic_slip_n
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::S2PK,F,F_p_inv_n
  REAL(KIND=rk),DIMENSION(6,3,3),INTENT(IN)::orth_basis
  REAL(KIND=rk),DIMENSION(6,6,6),INTENT(IN)::orth_proj
  REAL(KIND=rk)::Je_n1_ptb,det_Fp_inv_ptb,aver_adj_stp_ptb
  REAL(KIND=rk),DIMENSION(6)::S2PK_vec,S2PK_ptb_vec
  REAL(KIND=rk),DIMENSION(24)::plastic_slip_ptb
  REAL(KIND=rk),DIMENSION(3,3)::F_inv,pertub,F_pertub,S2PK_ptb
  REAL(KIND=rk),DIMENSION(3,3)::S2PK_IC_ptb,F_p_inv_n1_ptb,C_ptb
  REAL(KIND=rk),DIMENSION(3,3)::Ce_n1_ptb,Ce_n1_inv_ptb
  REAL(KIND=rk),DIMENSION(3,3)::F_p_inv_iter_ptb
  REAL(KIND=rk),DIMENSION(6,6)::C_tang
  REAL(KIND=rk),DIMENSION(6,6)::CSDA_material_tangent_single_crystal_plast
  REAL(KIND=rk),DIMENSION(24,6)::d_delta_lambda_ip1_d_C_ptb
  INTEGER(KIND=ik),INTENT(IN)::symm_class,algo,numtang,NDI,NSHR
  INTEGER(KIND=ik),INTENT(IN)::hard_law
  INTEGER(KIND=ik)::j_cmt
  INTEGER(KIND=ik),DIMENSION(6,2)::index_list
  !
  LOGICAL::OK_FLAG,conv_issue_dil_ptb,run_ptb_loop
  LOGICAL,DIMENSION(24)::active_slip_system_ptb
  CHARACTER(4),INTENT(IN)::env
  !
  ! parameter for numerical tangent computation
  index_list=approp_index(SIX_int,env,'CSDA_material_tangent_single_crystal_plast')
  CALL transform_symm_tens_vec_form(S2PK,S2PK_vec,ZERO_int,env)
  CALL M33INV (F, F_inv, OK_FLAG)
  IF (.NOT.(OK_FLAG)) THEN
     WRITE(16,'(A)')'ERROR in: C"SDA_material_tangent_single_crystal_plast"'
     WRITE(16,'(A)')'M33INV encountered singular matrix'
     STOP
  END IF
  !
  j_cmt=ONE_int; run_ptb_loop=.TRUE.
  !
  ! pertubation loop
  DO WHILE (run_ptb_loop)
     pertub=pertub_F_numtang_material_tang(F_inv,index_list(j_cmt,:))
     SELECT CASE (numtang)
     CASE (ONE_int) ! CSDA tangent
        WRITE(16,'(A)')'CSDA tangent computation is not implemented in "CSDA_material_tangent_single_crystal_plast"'
        WRITE(16,'(A)')'declaration of complex variables required'
        WRITE(16,'(A)')'STOPPING COMPUTATION'
        STOP
     CASE (TWO_int) ! FD tangent
        F_pertub = F + pertub*rel_step_size/TWO_r
     CASE DEFAULT
        WRITE(16,'(A)')'Unknown numtang ID in "CSDA_material_tangent_single_crystal_plast"'
        WRITE(16,'(A)')'STOPPING COMPUTATION'
        STOP
     END SELECT
     !
     aver_adj_stp_ptb=aver_adj_stp
     CALL stress_computation_single_crystal_plasticity(S2PK_ptb,&
          &S2PK_IC_ptb,F_p_inv_n1_ptb,plastic_slip_ptb,&
          &d_delta_lambda_ip1_d_C_ptb,active_slip_system_ptb,Ce_n1_ptb,&
          &Ce_n1_inv_ptb,Je_n1_ptb,C_ptb,det_Fp_inv_ptb,&
          &F_p_inv_iter_ptb,F_pertub,dt,F_p_inv_n,plastic_slip_n,&
          &elast_const,mat_param_inelast,m_SH,symm_class,orth_basis,&
          &orth_proj,NDI,NSHR,algo,env,hard_law,theta_KK,&
          &aver_adj_stp_ptb,conv_issue_dil_ptb)
     !
     CALL transform_symm_tens_vec_form(S2PK_ptb,S2PK_ptb_vec,ZERO_int,&
          &env)
     SELECT CASE (numtang)
     CASE (ONE_int) ! CSDA tangent
        WRITE(16,'(A)')'CSDA tangent computation is not implemented in "CSDA_material_tangent_single_crystal_plast (2)"'
        WRITE(16,'(A)')'declaration of complex variables required'
        WRITE(16,'(A)')'STOPPING COMPUTATION'
        STOP
     CASE (TWO_int) ! FD tangent
        C_tang(:,j_cmt)=(S2PK_ptb_vec-S2PK_vec)/rel_step_size
     CASE DEFAULT
        WRITE(16,'(A)')'Unknown numtang ID in "CSDA_material_tangent_single_crystal_plast (2)"'
        WRITE(16,'(A)')'STOPPING COMPUTATION'
        STOP
     END SELECT
     !
     ! increment counter of pertubations and exit loop
     j_cmt=j_cmt+ONE_int
     IF (j_cmt > SIX_int) THEN
        run_ptb_loop=.FALSE.
     END IF
     !
     ! emergency exit in case of unsuccessful stress computation
     IF (conv_issue_dil_ptb) THEN
        run_ptb_loop=.FALSE.
     END IF
     !
  END DO
  !
  CSDA_material_tangent_single_crystal_plast=C_tang
  !
END FUNCTION CSDA_material_tangent_single_crystal_plast
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! compute parts of the derivative of the second Piola-Kirchhoff tensor
! w.r.t. the inverse of the plastic deformation gradient
! **********************************************************************
FUNCTION deriv_d_S_d_Fp_inv(S2PK_IC,F_p_inv,env)
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:expand_fourth_order_tensor_matrix2tensor,&
       & spec_tensor_product_unsym_second_order_tensors, &
       & square_product_unsym_second_order_tensors
  !
  ! variable declaration
  ! S2PK_IC............second Piola-Kirchhoff stress in intermed. conf.
  ! F_p_inv............inverse of plastic part of deformation gradient
  IMPLICIT NONE
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::S2PK_IC,F_p_inv
  REAL(KIND=rk),DIMENSION(3,3)::prod_Fp_i_S2PK_IC
  REAL(KIND=rk),DIMENSION(3,3,3,3)::d_S_d_Fp_inv
  REAL(KIND=rk),DIMENSION(3,3,3,3)::deriv_d_S_d_Fp_inv
  !
  CHARACTER(4),INTENT(IN)::env
  !
  prod_Fp_i_S2PK_IC=MATMUL(F_p_inv,S2PK_IC)
  d_S_d_Fp_inv=expand_fourth_order_tensor_matrix2tensor(&
       &spec_tensor_product_unsym_second_order_tensors(prod_Fp_i_S2PK_IC,&
       &IDENT,env),env)
  d_S_d_Fp_inv=d_S_d_Fp_inv+expand_fourth_order_tensor_matrix2tensor(&
       &square_product_unsym_second_order_tensors(IDENT,prod_Fp_i_S2PK_IC,&
       &env),env)
  !
  deriv_d_S_d_Fp_inv=d_S_d_Fp_inv
  !
END FUNCTION deriv_d_S_d_Fp_inv
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! compute the derivative of the inverse of the plastic deformation
! gradient w.r.t. right Cauchy-Green tensor C
! **********************************************************************
FUNCTION derivative_inverse_plastic_def_grad_right_cauchy_green_tensor_c(&
     &F_p_inv_n,active_slip_system,d_delta_lambda_ip1_d_C,env)
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:comp_sec_ord_tens_rank_one, &
       & DYADIC_PRODUCT_SECOND_ORDER_TENSORS, &
       & transform_symm_tens_vec_form
  USE constitutive_routines_interface, only:slip_system_fcc
  !
  ! variable declaration
  ! F_p_inv_n...............inverse of plastic part of deformation
  ! active_slip_system......list of active slip systems
  ! d_delta_lambda_ip1_d_C..deriv. of incremental lagrange multipliers
  !                         w.r.t. right Cauchy-Green tensor
  IMPLICIT NONE
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::F_p_inv_n
  REAL(KIND=rk),DIMENSION(24,6),INTENT(IN)::d_delta_lambda_ip1_d_C
  REAL(KIND=rk),DIMENSION(3)::normal_dir,slip_dir
  REAL(KIND=rk),DIMENSION(3,3)::N_alpha,d_delta_lambda_alpha_ip1_d_C
  REAL(KIND=rk),DIMENSION(3,3)::F_p_inv_n_sing_contr_N_alpha
  REAL(KIND=rk),DIMENSION(3,3,3,3)::d_F_p_inv_n1_d_C
  REAL(KIND=rk),DIMENSION(3,3,3,3)::derivative_inverse_plastic_def_grad_right_cauchy_green_tensor_c
  INTEGER(KIND=ik)::i_slip
  LOGICAL,DIMENSION(24),INTENT(IN)::active_slip_system
  CHARACTER(4),INTENT(IN)::env
  !
  d_F_p_inv_n1_d_C=ZERO_r
  DO i_slip=ONE_int,TWENTY_FOUR_int
     IF (active_slip_system(i_slip)) THEN
        ! tensor defining the slip system
        CALL slip_system_fcc(normal_dir,slip_dir,i_slip)
        N_alpha=comp_sec_ord_tens_rank_one(slip_dir,normal_dir)
        CALL transform_symm_tens_vec_form(d_delta_lambda_alpha_ip1_d_C,&
             &d_delta_lambda_ip1_d_C(i_slip,:),ONE_int,env)
        F_p_inv_n_sing_contr_N_alpha=MATMUL(F_p_inv_n,N_alpha)
        d_F_p_inv_n1_d_C = d_F_p_inv_n1_d_C - &
             &dyadic_product_second_order_tensors(&
             &F_p_inv_n_sing_contr_N_alpha,&
             &d_delta_lambda_alpha_ip1_d_C)
     END IF
  END DO
  !
  derivative_inverse_plastic_def_grad_right_cauchy_green_tensor_c=&
       &d_F_p_inv_n1_d_C
  !
END FUNCTION derivative_inverse_plastic_def_grad_right_cauchy_green_tensor_c
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! compute the derivative of the Mandel stress w.r.t. elastic right
! Cauchy-Green tensor, employing quadratic elastic energy based on
! Seth-Hill generalized strain measure
! **********************************************************************
FUNCTION derivative_mandel_stress_ce_sh_elasticity(S_2PK_IC,C_mat_elast,&
     &Ce,env)
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:expand_fourth_order_tensor_matrix2tensor, &
       & expand_minor_symm_fourth_order_tensor_matrix2tensor, &
       & spec_tensor_product_unsym_second_order_tensors, &
       & square_product_unsym_second_order_tensors
  !
  ! declaration of variables
  ! S_2PK_IC...........second Piola-Kirchhoff stress in intermed. config
  ! C_mat_elast........material tangent in intermed. config., i.e. deriv
  !                    of S_2PK w.r.t. elastic right Cauchy-Green tensor
  ! Ce.................elastic right Cauchy-Green tensor
  !
  IMPLICIT NONE
  REAL(KIND=rk)::temp_var
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::S_2PK_IC,Ce
  REAL(KIND=rk),DIMENSION(6,6),INTENT(IN)::C_mat_elast
  REAL(KIND=rk),DIMENSION(3,3,3,3)::d_Mandel_stress_d_Ce
  REAL(KIND=rk),DIMENSION(3,3,3,3)::derivative_mandel_stress_ce_sh_elasticity
  REAL(KIND=rk),DIMENSION(3,3,3,3)::C_mat_elast_tens,Ce_sc_C_mat_elast
  INTEGER(KIND=ik)::m_dmsc,i_dmsc,j_dmsc,k_dmsc,l_dmsc
  !
  CHARACTER(4),INTENT(IN)::env
  !
  d_Mandel_stress_d_Ce=expand_fourth_order_tensor_matrix2tensor(&
       &ONE_HALF_r*(spec_tensor_product_unsym_second_order_tensors(IDENT,&
       &S_2PK_IC,env) + &
       &square_product_unsym_second_order_tensors (IDENT,&
       &S_2PK_IC,env)),env)
  C_mat_elast_tens=expand_minor_symm_fourth_order_tensor_matrix2tensor(&
       &C_mat_elast,env)
  DO i_dmsc=ONE_int,THREE_int
     DO j_dmsc=ONE_int,THREE_int
        DO k_dmsc=ONE_int,THREE_int
           DO l_dmsc=k_dmsc,THREE_int
              temp_var=ZERO_r
              DO m_dmsc=ONE_int,THREE_int
                 temp_var=temp_var+Ce(i_dmsc,m_dmsc)*&
                      &C_mat_elast_tens(m_dmsc,j_dmsc,k_dmsc,l_dmsc)
              END DO ! m
              Ce_sc_C_mat_elast(i_dmsc,j_dmsc,k_dmsc,l_dmsc)=temp_var
              IF (l_dmsc>k_dmsc) THEN
                 Ce_sc_C_mat_elast(i_dmsc,j_dmsc,l_dmsc,k_dmsc)=temp_var
              END IF
           END DO ! l
        END DO ! k
     END DO ! j
  END DO ! i
  !
  d_Mandel_stress_d_Ce=d_Mandel_stress_d_Ce+ONE_HALF_r*Ce_sc_C_mat_elast
  derivative_mandel_stress_ce_sh_elasticity=d_Mandel_stress_d_Ce
  !
END FUNCTION derivative_mandel_stress_ce_sh_elasticity
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! extract the hardening parameters from the list of material parameters
! associated with inelastic material behavior
! **********************************************************************
FUNCTION extract_hardening_parameters_crystal_plast(mat_param_inelast,&
     &hard_law,dim_hard)
  !
  USE Abaqus_Interface
  USE Constants
  !
  ! variable declaration
  ! mat_param_inelast..material parameter assoc. with inelastic behavior
  ! hard_law...........id of hardening law
  !                    0....nonlin. Taylor-type (Schmidt-Baldassari)
  !                    1....nonlin. interaction matrix (Cailletaud,Forest)
  !
  IMPLICIT NONE
  !
  INTEGER(KIND=ik),INTENT(IN)::hard_law,dim_hard
  REAL(KIND=rk),DIMENSION(13),INTENT(IN)::mat_param_inelast
  REAL(KIND=rk),DIMENSION(dim_hard)::extract_hardening_parameters_crystal_plast
  !
  SELECT CASE (hard_law)
  CASE (ZERO_int) ! nonlinear Taylor-type hardening
     extract_hardening_parameters_crystal_plast=&
          &mat_param_inelast(ONE_int:dim_hard)
  CASE (ONE_int)  ! nonlinear hardening based on interaction matrix
     extract_hardening_parameters_crystal_plast=&
          &mat_param_inelast(ONE_int:dim_hard)
  CASE DEFAULT
     WRITE(16,'(A)')'Unknown value of hard_law in: "extract_hardening_parameters_crystal_plast"'
     WRITE(16,'(A)')'STOPPING COMPUTATION'
     STOP
  END SELECT
  !
END FUNCTION extract_hardening_parameters_crystal_plast
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! extract the viscoplastic parameters from the list of material
! parameters associated with inelastic material behavior
! **********************************************************************
FUNCTION extract_viscoplastic_parameters_crystal_plast(&
     &mat_param_inelast)
  !
  USE Abaqus_Interface
  USE Constants
  !
  ! variable declaration
  ! mat_param_inelast..material parameter assoc. with inelastic behavior
  !
  IMPLICIT NONE
  !
  REAL(KIND=rk),DIMENSION(13),INTENT(IN)::mat_param_inelast
  REAL(KIND=rk),DIMENSION(3)::extract_viscoplastic_parameters_crystal_plast
  !
  extract_viscoplastic_parameters_crystal_plast=mat_param_inelast(11:13)
  !
END FUNCTION extract_viscoplastic_parameters_crystal_plast
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! interaction matrix for single crystal (visco) plasticity model,
! accounting for interaction of different slip systems
! **********************************************************************
FUNCTION interaction_matrix(hard_param_im)
  !
  USE Abaqus_Interface
  USE Constants
  !
  ! explanation taken from Khadyko et al., IJP 76(2016), p. 51-74
  ! Franciosi and Zaoui, Acta Metall. 30 (1982), p. 1627-1637
  ! h_0.......self-hardening
  ! h_1.......interaction of coplanar dislocations
  ! h_2.......Hirth lock on system pair with normal slip directions
  ! h_3.......interaction of slip systems with colinear slip directions
  ! h_4.......glissile lock, interaction of non-coplanar slip systems
  ! h_5.......sessile / Lomer-Cottrell lock, interaction of non-coplanar
  !           slip systems
  ! variable declaration
  ! hard_param_im......material parameter of hardening matrix
  IMPLICIT NONE
  REAL(KIND=rk),DIMENSION(6),INTENT(IN)::hard_param_im
  REAL(KIND=rk)::h_0,h_1,h_2,h_3,h_4,h_5
  REAL(KIND=rk),DIMENSION(12,12)::interact_mat_part
  REAL(KIND=rk),DIMENSION(24,24)::interaction_matrix
  !
  ! assign material parameters
  h_0=hard_param_im(1); h_1=hard_param_im(2); h_2=hard_param_im(3)
  h_3=hard_param_im(4); h_4=hard_param_im(5); h_5=hard_param_im(6)
  !
  interact_mat_part=RESHAPE(&
       &(/h_0, h_1, h_1, h_3, h_4, h_4, h_2, h_4, h_5, h_2, h_5, h_4, &
       &  h_1, h_0, h_1, h_4, h_2, h_5, h_4, h_3, h_4, h_5, h_2, h_4, &
       &  h_1, h_1, h_0, h_4, h_5, h_2, h_5, h_4, h_2, h_5, h_4, h_3, &
       &  h_3, h_4, h_4, h_0, h_1, h_1, h_2, h_5, h_4, h_2, h_4, h_5, &
       &  h_4, h_2, h_5, h_1, h_0, h_1, h_5, h_2, h_4, h_4, h_3, h_4, &
       &  h_4, h_5, h_2, h_1, h_1, h_0, h_4, h_4, h_3, h_5, h_4, h_2, &
       &  h_2, h_4, h_5, h_2, h_5, h_4, h_0, h_1, h_1, h_3, h_4, h_4, &
       &  h_4, h_3, h_4, h_5, h_2, h_4, h_1, h_0, h_1, h_4, h_2, h_5, &
       &  h_5, h_4, h_2, h_4, h_4, h_3, h_1, h_1, h_0, h_4, h_5, h_2, &
       &  h_2, h_5, h_5, h_2, h_4, h_5, h_3, h_4, h_4, h_0, h_1, h_1, &
       &  h_5, h_2, h_4, h_4, h_3, h_4, h_4, h_2, h_5, h_1, h_0, h_1, &
       &  h_4, h_4, h_3, h_5, h_4, h_2, h_4, h_5, h_2, h_1, h_1, h_0/),(/12,12/))
  !
  interaction_matrix(1:12,1:12)=interact_mat_part
  interaction_matrix(1:12,13:24)=interact_mat_part
  interaction_matrix(13:24,1:12)=interact_mat_part
  interaction_matrix(13:24,13:24)=interact_mat_part
  !
END FUNCTION interaction_matrix
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! compute the derivative of the (deviatoric) Mandel stress w.r.t.
! increment in lagrange multipliers
! ATTENTION: linearization relies on deviatoric character of Schmid-tensor 
! (i.e. dyadic product of slip direction and slip plane normal)
! **********************************************************************
FUNCTION linearization_mandel_stress(F_p_inv_n,unimod_C_F_p_inv_iter,mu,&
     &env)
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:comp_sec_ord_tens_rank_one,&
       & transform_tens_vec_form
  USE constitutive_routines_interface, only:slip_system_fcc
  !
  ! variable declaration
  ! F_p_inv_n.............inverse of plastic part of deformation gradient
  !                       at t_{n}
  ! unimod_C_F_p_inv_iter.single contraction of unimodular right Cauchy-
  !                       Green tensor and inverse of plastic def. grad.
  ! mu....................shear modulus
  IMPLICIT NONE
  REAL(KIND=rk),INTENT(IN)::mu
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::F_p_inv_n
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::unimod_C_F_p_inv_iter
  REAL(KIND=rk),DIMENSION(3)::normal_dir,slip_dir
  REAL(KIND=rk),DIMENSION(9)::d_Mandel_stress_d_delta_lambda_temp
  REAL(KIND=rk),DIMENSION(3,3)::N_beta,d_F_p_inv_d_delta_lambda_beta
  REAL(KIND=rk),DIMENSION(3,3)::d_unimod_Ce_d_delta_lambda_ip1_beta
  REAL(KIND=rk),DIMENSION(9,24)::linearization_mandel_stress
  INTEGER(KIND=ik)::i_slip_sys
  !
  CHARACTER(4),INTENT(IN)::env
  !
  DO i_slip_sys=ONE_int,TWENTY_FOUR_int
     CALL slip_system_fcc(normal_dir,slip_dir,i_slip_sys)
     N_beta=comp_sec_ord_tens_rank_one(slip_dir,normal_dir)
     d_F_p_inv_d_delta_lambda_beta=-MATMUL(F_p_inv_n,N_beta)
     !
     d_unimod_Ce_d_delta_lambda_ip1_beta=MATMUL(&
          & TRANSPOSE(d_F_p_inv_d_delta_lambda_beta),unimod_C_F_p_inv_iter)+&
          & MATMUL(TRANSPOSE(unimod_C_F_p_inv_iter),&
          &d_F_p_inv_d_delta_lambda_beta)
     !
     CALL transform_tens_vec_form(d_unimod_Ce_d_delta_lambda_ip1_beta,&
          &d_Mandel_stress_d_delta_lambda_temp,ZERO_int,env)
     linearization_mandel_stress(:,i_slip_sys)=mu*&
          &d_Mandel_stress_d_delta_lambda_temp
  END DO
  !
END FUNCTION linearization_mandel_stress
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! compute the derivative of the Mandel stress w.r.t. increment in 
! lagrange multipliers
! **********************************************************************
FUNCTION linearization_mandel_stress_elasticity_SH_strain(F_p_inv_n,&
     &unimod_C_F_p_inv_iter,Ce,S_2PK_IC,C_mat_elast,env)
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:comp_sec_ord_tens_rank_one,&
       & double_contraction_fourth_order_tensor_second_order_tensor,&
       & transform_tens_vec_form
  USE constitutive_routines_interface, only:derivative_mandel_stress_ce_sh_elasticity,&
       & slip_system_fcc
  !
  ! variable declaration
  ! F_p_inv_n.............inverse of plastic part of deformation gradient
  !                       at t_{n}
  ! unimod_C_F_p_inv_iter.single contraction of unimodular right Cauchy-
  !                       Green tensor and inverse of plastic def. grad.
  ! Ce....................elastic right Cauchy-Green tensor
  ! S_2PK_IC..............second Piola-Kirchhoff stress tensor in intermed.
  !                       conf.
  ! C_mat_elast...........elastic stiffness (deriv. S2PK w.r.t. Ce)
  !
  IMPLICIT NONE
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::F_p_inv_n,Ce,S_2PK_IC
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::unimod_C_F_p_inv_iter
  REAL(KIND=rk),DIMENSION(6,6),INTENT(IN)::C_mat_elast
  REAL(KIND=rk),DIMENSION(3)::normal_dir,slip_dir
  REAL(KIND=rk),DIMENSION(9)::d_Mandel_stress_d_delta_lambda_temp
  REAL(KIND=rk),DIMENSION(3,3)::N_beta,d_F_p_inv_d_delta_lambda_beta
  REAL(KIND=rk),DIMENSION(3,3)::d_unimod_Ce_d_delta_lambda_ip1_beta
  REAL(KIND=rk),DIMENSION(3,3)::d_Mandel_stress_d_Ce_dc_d_unimod_Ce_d_delta_lambda_ip1_beta
  REAL(KIND=rk),DIMENSION(9,24)::linearization_mandel_stress_elasticity_SH_strain
  REAL(KIND=rk),DIMENSION(3,3,3,3)::d_Mandel_stress_d_Ce
  INTEGER(KIND=ik)::i_slip_sys
  !
  CHARACTER(4),INTENT(IN)::env
  !
  ! compute derivative of Mandel stress w.r.t elastic right Cauchy-Green
  ! tensor
  d_Mandel_stress_d_Ce=derivative_mandel_stress_ce_sh_elasticity(&
       &S_2PK_IC,C_mat_elast,Ce,env)
  !
  DO i_slip_sys=ONE_int,TWENTY_FOUR_int
     CALL slip_system_fcc(normal_dir,slip_dir,i_slip_sys)
     N_beta=comp_sec_ord_tens_rank_one(slip_dir,normal_dir)
     d_F_p_inv_d_delta_lambda_beta=-MATMUL(F_p_inv_n,N_beta)
     !
     d_unimod_Ce_d_delta_lambda_ip1_beta=MATMUL(&
          & TRANSPOSE(d_F_p_inv_d_delta_lambda_beta),unimod_C_F_p_inv_iter)+&
          & MATMUL(TRANSPOSE(unimod_C_F_p_inv_iter),&
          &d_F_p_inv_d_delta_lambda_beta)
     !
     d_Mandel_stress_d_Ce_dc_d_unimod_Ce_d_delta_lambda_ip1_beta=&
          &double_contraction_fourth_order_tensor_second_order_tensor(&
          &d_Mandel_stress_d_Ce,d_unimod_Ce_d_delta_lambda_ip1_beta)
     CALL transform_tens_vec_form(&
          &d_Mandel_stress_d_Ce_dc_d_unimod_Ce_d_delta_lambda_ip1_beta,&
          &d_Mandel_stress_d_delta_lambda_temp,ZERO_int,env)
     !
     linearization_mandel_stress_elasticity_SH_strain(:,i_slip_sys)=&
          &d_Mandel_stress_d_delta_lambda_temp
  END DO
  !
END FUNCTION linearization_mandel_stress_elasticity_SH_strain
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! modification of the derivative of the inverse of the plastic
! deformation gradient w.r.t. right Cauchy-Green tensor due to the
! isochoric projection algorithm that enforces isochoric constraint of
! F_p_inv
! **********************************************************************
FUNCTION mod_d_F_p_inv_d_C_isochoric_proj_constraint(d_F_p_inv_i1_d_C,&
     &det_Fp_inv,F_p_inv_iter)
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:DYADIC_PRODUCT_SECOND_ORDER_TENSORS,&
       & left_double_contraction_fourth_order_tensor_second_order_tensor, &
       & M33INV
  !
  ! variable declaration
  ! d_F_p_inv_i1_d_C.......deriv. of inverse of plastic part of deform.
  !                        gradient w.r.t. right Cauchy-Green tensor
  ! det_Fp_inv.............determinant of inverse of plastic part of
  !                        deformation gradient
  ! F_p_inv_iter...........inverse of plastic part of deform. grad.
  !
  IMPLICIT NONE
  REAL(KIND=rk),INTENT(IN)::det_Fp_inv
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::F_p_inv_iter
  REAL(KIND=rk),DIMENSION(3,3,3,3),INTENT(IN)::d_F_p_inv_i1_d_C
  REAL(KIND=rk)::det_F_p_inv_p_m1_3
  REAL(KIND=rk),DIMENSION(3,3)::F_p_i1,d_det_F_p_inv_d_C
  REAL(KIND=rk),DIMENSION(3,3,3,3)::d_F_p_inv_n1_d_C
  REAL(KIND=rk),DIMENSION(3,3,3,3)::mod_d_F_p_inv_d_C_isochoric_proj_constraint
  LOGICAL::OK_FLAG
  !
  CALL M33INV (F_p_inv_iter, F_p_i1, OK_FLAG)
  IF (.NOT.(OK_FLAG)) THEN
     WRITE(16,'(A)')'ERROR in: "mod_d_F_p_inv_d_C_isochoric_proj_constraint"'
     WRITE(16,'(A)')'M33INV encountered singular matrix'
     STOP
  END IF
  d_det_F_p_inv_d_C=&
       &left_double_contraction_fourth_order_tensor_second_order_tensor(&
       &d_F_p_inv_i1_d_C,TRANSPOSE(F_p_i1))
  det_F_p_inv_p_m1_3=det_Fp_inv**(-ONE_by_THREE_r)
  d_det_F_p_inv_d_C=-ONE_by_THREE_r*det_F_p_inv_p_m1_3*d_det_F_p_inv_d_C
  !
  d_F_p_inv_n1_d_C=det_F_p_inv_p_m1_3*d_F_p_inv_i1_d_C+&
       & dyadic_product_second_order_tensors(F_p_inv_iter,&
       & d_det_F_p_inv_d_C)
  !
  mod_d_F_p_inv_d_C_isochoric_proj_constraint=d_F_p_inv_n1_d_C
  !
END FUNCTION mod_d_F_p_inv_d_C_isochoric_proj_constraint
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! compute the second Piola-Kirchhoff stress for isotropic volumetric-
! isochoric decoupled elasticity (neo-hookean type)
! employed in Schmidt-Baldassari, CMAME, 192 (2003), pp.1261-1280
! **********************************************************************
FUNCTION S2PK_IC_isotr_elasticity(Ce,Ce_inv,Je,kappa,mu)
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:compute_trace_second_order_tensor
  !
  ! declaration of variables
  ! Ce...............elastic right Cauchy-Green tensor
  ! Ce_inv...........inverse of elastic right Cauchy-Green tensor
  ! Je...............determinant of elastic part of deformation grad.
  ! kappa,mu.........isotropic elastic constants (bulk-, shear modulus)
  ! S2PK_IC_.........second Piola-Kirchhoff stress in intermed. conf.
  !
  IMPLICIT NONE
  REAL(KIND=rk),INTENT(IN)::Je,kappa,mu
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::Ce,Ce_inv
  REAL(KIND=rk)::tr_Ce
  REAL(KIND=rk),DIMENSION(3,3)::S2PK_IC_isotr_elasticity
  !
  tr_Ce=compute_trace_second_order_tensor(Ce)
  S2PK_IC_isotr_elasticity=kappa*log(Je)*Ce_inv+Je**(-TWO_by_THREE_r)*&
       &mu*(IDENT-ONE_by_THREE_r*tr_Ce*Ce_inv)
  !
END FUNCTION S2PK_IC_isotr_elasticity
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! update plastic part of the deformation gradient
! (íntegration based on euler-implicit, additional projection required
! to enforce isochoric nature of F_p)
! **********************************************************************
FUNCTION update_plastic_deformation_gradient(F_p_inv_n,delta_lambda_ip1)
  !
  USE Abaqus_Interface
  USE Constants
  USE auxiliary_routines_interface, only:comp_sec_ord_tens_rank_one
  USE constitutive_routines_interface, only:slip_system_fcc
  !
  ! variable declaration
  ! F_p_inv_n..........inverse of plastic part of deformation gradient
  ! delta_lambda_ip1...incremental lagrange multiplier
  !
  IMPLICIT NONE
  REAL(KIND=rk),DIMENSION(3)::normal_dir,slip_dir
  REAL(KIND=rk),DIMENSION(24),INTENT(IN)::delta_lambda_ip1
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::F_p_inv_n
  REAL(KIND=rk),DIMENSION(3,3)::update_plastic_deformation_gradient
  REAL(KIND=rk),DIMENSION(3,3)::relative_F_p_inv,N_alpha
  INTEGER(KIND=ik)::i_slip
  !
  relative_F_p_inv=IDENT
  !
  DO i_slip=ONE_int,TWENTY_FOUR_int
     CALL slip_system_fcc(normal_dir,slip_dir,i_slip)
     N_alpha=comp_sec_ord_tens_rank_one(slip_dir,normal_dir)
     relative_F_p_inv=relative_F_p_inv-delta_lambda_ip1(i_slip)*N_alpha
  END DO
  !
  ! F_p_inv_iter at t_{n+1}
  update_plastic_deformation_gradient=MATMUL(F_p_inv_n,relative_F_p_inv)
  !
END FUNCTION update_plastic_deformation_gradient
! ......................................................................
! **********************************************************************
!
! **********************************************************************
! update the vector of material parameters associated with inelastic
! material behavior due to the adjustment of viscoplastic parameters
! **********************************************************************
FUNCTION update_visco_in_inelast_param_crystal_plast(mat_param_visco_t,&
     &mat_param_inelast)
  !
  USE Abaqus_Interface
  USE Constants
  !
  ! declaration of variables
  ! mat_param_visco_t....material constants describing rate-dependent
  !                      behavior
  ! mat_param_inelast_t..material constants describing inelastic behavior
  !
  IMPLICIT NONE
  REAL(KIND=rk),DIMENSION(3),INTENT(IN)::mat_param_visco_t
  REAL(KIND=rk),DIMENSION(13),INTENT(IN)::mat_param_inelast
  REAL(KIND=rk),DIMENSION(13)::update_visco_in_inelast_param_crystal_plast
  INTEGER(KIND=ik)::i_uvi
  !
  update_visco_in_inelast_param_crystal_plast=mat_param_inelast
  DO i_uvi=ONE_int,TWO_int
     update_visco_in_inelast_param_crystal_plast(TEN_int+i_uvi)=&
          &mat_param_visco_t(i_uvi)
  END DO
  !
END FUNCTION update_visco_in_inelast_param_crystal_plast
! ......................................................................
! **********************************************************************
! END                  Constitutive_routines
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! UMAT for single crystal plasticity, employing isotropic and cubic
! elasticity based on quadratic form of the strain energy function in
! terms of general Seth-Hill strain measure and anisotropic hardening
! based on interaction matrix (up to 6 parameter), rate-dependent and
! rate-independent formulation with different algorithms, finite
! deformations framework
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,DRPLDE,&
     & DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,&
     & NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,&
     & DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
  !
  ! include definition of parameters and other modules
  USE Abaqus_Interface
  USE Abaqus_Interface_extended
  USE Constants
  USE auxiliary_routines_interface, only:comp_det_Matrix, &
       & ROTATION_MATRIX_BUNGE, pushrc,transform_symm_tens_vec_form,&
       & transform_tens_vec_form
  USE constitutive_routines_interface, only:comp_add_stress_contrib_tang,&
       & crystal_plasticity_augm_lagr, orthonormal_basis_subspace
  !
  ! ********************************************************************
  ! usage: abaqus job=<job-name> user=<path-to-umat> cpus=<# of cpus>
  !        computation on multiple cpus requires parallel mpi set up in
  !        abaqus_v6.env by enforcing fully parallel mpi
  !
  ! input file:
  ! ** MATERIALS
  ! *Material, name=<material-name>
  ! *Depvar
  !   38,
  ! *User Material, constants=25
  !  c1,  c2,  c3,  c4,  c5,  c6,  c7,  c8
  !  c9, c10, c11, c12, c13, c14, c15, c16
  ! c17, c18, c19, c20, c21, c22, c23, c24
  ! c25
  ! **
  ! c1............elastic law
  !               0...isotropic neo-hookean elasticity with
  !                   isochoric-volumetric split
  !               1...isotropic elasticity based on Seth-Hill strain
  !                   measure
  !               2...cubic elasticity based on Seth-Hill strain measure
  ! c2............Seth-Hill parameter (only used for c1=1 or c1=2)
  !               -2 <= c2 <= 2
  ! c3............formulation of inelastic mechanism
  !               0...augmented-lagrangian (rate-independent)
  !               1...NCP-function (Kanzow-Kleinmichel) (rate-independent)
  !               2...Perzyna-type viscoplastic formulation
  !               3...Ortiz/Miehe viscoplastic formulation
  ! c4............algorithmic parameter for NCP-function (0 <= c4 <= 4)
  ! c5............algorithmic tangent computation
  !               0...analytical tangent
  !               1...CSDA (implementation pending)
  !               2...FD tangent computation
  ! c6............hardening law
  !               0...nonlinear exponential-type Taylor hardening
  !               1...nonlinear exponential-type hardening based on
  !                   interaction matrix
  ! c7,c8,c9......orientation given in Euler angles due to Bunge in degree
  ! c10,c11,c12...elastic constants
  !               if c1=0 or c1=1: c10=bulk modulus, c11=shear modulus, c12=arbitr.
  !               if c1=2:         c10=C11, c11=C12, c12=C44
  !                                (three elastic constants for cubic elastic)
  ! c13...........initial yield
  ! c14...........asymptotic increase in initial yield
  ! c15...........shape parameter for expontial law
  !               controls how fast the asymptotic yield stress is reached
  ! c16-c21.......six parameters of interaction matrix
  ! c22...........hardening modulus for linear hardening
  ! c23...........exponent in rate-dependent formulation
  ! c24...........time-like parameter or 1/reference strain rate (viscoplastic)
  ! c25...........drag stress (only required for c3=2 Perzyna-type viscoplastic)
  !
  ! ********************************************************************
  ! declaration of variables
  !
  ! -------- Abaqus variables ------------
  ! STRESS..........Cauchy stress tensor (input and update)
  ! STRAN...........strain tensor at t_n
  !                 (changes necessary to account for thermal expansion)
  ! DSTRAN..........strain increment (incremental symm. part of velocity
  !                 gradient, changes necessary to account for thermal
  !                 expansion)
  ! DDSDDE..........algorithmic tangent at t_n+1
  ! PREDEF..........predefined field variables at t_n
  !                 (interpolated from nodal values)
  ! DPRED...........increment in predefined field variables
  ! TIME............time at beginning of the increment
  !                 (1)-total time of the current step
  !                 (2)-total time
  ! DTIME...........time increment
  ! DDSDDT..........derivative of Cauchy stress w.r.t. temperature
  ! RPL.............heat source density per time due to mechanical work
  !                 at t_n+1
  ! DRPLDE..........derivative of RPL w.r.t. strain increment
  ! DRPLDT..........derivative of RPL w.r.t. temperature
  ! DROT............incremental rotation matrix (due to Jaumann rate)
  ! DFGRD0..........deformation gradient at t_n
  ! DFGRD1..........deformation gradient at t_n+1
  ! COORDS..........coordinates of the material point (integration point)
  ! PROPS...........material constants employed in UMAT
  ! NPROPS..........# of material constants
  ! NTENS...........# of stress and strain components
  ! NSTATV..........# of internal state variables
  ! STATEV..........internal state variables in UMAT (input and update)
  ! SSE.............elastic stored energy (input and update)
  ! SPD.............plastic dissipation (input and update)
  ! SCD.............dissipation due to creep (input and update)
  ! TEMP............temperature at t_n
  ! DTEMP...........temperature increment
  ! CMNAME..........name of a UMAT (allows the use of several user
  !                 materials in a single UMAT subroutine/Abaqus simul.)
  ! NDI.............# of normal stress components at integration point
  ! NSHR............# of shear stress components at integration point
  ! PNEWDT..........fraction of proposed new time increment and current
  !                 time increment
  ! CELENT..........characteristic element length
  ! NOEL............element number
  ! NPT.............integration point number
  ! LAYER...........# of layers (Composite Shell oder Solid)
  ! KSPT............# of current layer
  ! JSTEP...........step number
  ! KINC............# of increment
  ! 
  ! -------- own variables ------------
  ! symm_class.....symmetry class of crystal (0=isotropic, vol/isochoric split,
  !                1=isotropic, 2=cubic)
  ! m_SH...........parameter generalized Seth-Hill strain
  ! algo...........formulation of inelastic mechanism
  !                0...augmented-lagrangian (rate-independent)
  !                1...NCP-function (Kanzow-Kleinmichel) (rate-independent)
  !                2...Perzyna-type viscoplastic formulation
  !                3...Ortiz/Miehe viscoplastic formulation
  ! theta_KK.......algorithmic parameter for NCP-function (0 <= c4 <= 4)
  ! numtang........algorithmic tangent computation
  ! hard_law.......hardening law
  !                0...nonlinear exponential-type Taylor hardening
  !                1...nonlinear exponential-type hardening based on
  !                    interaction matrix
  ! orient.........orientation given in Euler angles due to Bunge in degree
  ! elast_const....elastic constants
  !                if c1=0 or c1=1: c10=bulk modulus, c11=shear modulus, c12=arbitr.
  !                if c1=2:         c10=C11, c11=C12, c12=C44
  !                                 (three elastic constants for cubic elastic)
  ! mat_param_inelast..material constants describing inelastic behavior
  !                    c13....c25
  ! orth_basis.....container of all second order, orthogonal basis N(i)
  ! orth_proj......container of symmetric dyadic products of orthogonal
  !                basis, i.e. N(i) dyadic N(i)
  !
  ! explicit variable declaration
  IMPLICIT NONE
  !
  ! -------- Abaqus variables ------------
  INTEGER(KIND=ik),INTENT(IN)::NPROPS,NTENS,NSTATV,NDI,NSHR,NOEL
  INTEGER(KIND=ik),INTENT(IN)::NPT,LAYER,KSPT,KINC!,KSTEP
  !
  INTEGER(KIND=ik),DIMENSION(4),INTENT(IN)::JSTEP
  !
  REAL(KIND=rk),INTENT(IN)::PREDEF,DPRED,DTIME,SSE,SPD,SCD
  REAL(KIND=rk),INTENT(IN)::TEMP,DTEMP,CELENT
  REAL(KIND=rk),INTENT(INOUT)::PNEWDT
  REAL(KIND=rk),INTENT(OUT)::RPL,DRPLDT
  !
  REAL(KIND=rk),DIMENSION(2),INTENT(IN)::TIME
  REAL(KIND=rk),DIMENSION(3),INTENT(IN)::COORDS
  REAL(KIND=rk),DIMENSION(NTENS),INTENT(IN)::STRAN,DSTRAN
  REAL(KIND=rk),DIMENSION(NTENS),INTENT(INOUT)::STRESS
  REAL(KIND=rk),DIMENSION(NTENS),INTENT(OUT)::DDSDDT,DRPLDE
  REAL(KIND=rk),DIMENSION(NPROPS),INTENT(IN)::PROPS
  REAL(KIND=rk),DIMENSION(NSTATV),INTENT(INOUT)::STATEV
  REAL(KIND=rk),DIMENSION(NTENS,NTENS),INTENT(OUT)::DDSDDE
  REAL(KIND=rk),DIMENSION(3,3),INTENT(IN)::DROT,DFGRD0,DFGRD1
  !
  ! ------ own variables ---------
  INTEGER(KIND=ik)::symm_class,algo,numtang,hard_law
  ! scalars
  REAL(KIND=rk)::m_SH,J_n1,theta_KK,DTIME_tp
  ! fields
  REAL(KIND=rk),DIMENSION(3)::orient,elast_const
  REAL(KIND=rk),DIMENSION(6)::S_2PK_vec
  REAL(KIND=rk),DIMENSION(13)::mat_param_inelast
  REAL(KIND=rk),DIMENSION(3,3)::R_init,S2PK
  REAL(KIND=rk),DIMENSION(6,6)::C_mat_tang
  REAL(KIND=rk),DIMENSION(6,3,3)::orth_basis
  REAL(KIND=rk),DIMENSION(6,6,6)::orth_proj
  !
  ! ------ other variables ------
  CHARACTER(4),PARAMETER::env_used='ABAQ'
  CHARACTER(80)::CMNAME
  LOGICAL::conv_issue
  !
  INTRINSIC MATMUL,TRANSPOSE
  !
  ! ********************************************************************
  !                           UMAT source code
  ! ********************************************************************
  !
  ! initialize output
  STRESS=ZERO_r;  DDSDDE=ZERO_r; RPL=ZERO_r
  DDSDDT=ZERO_r;  DRPLDE=ZERO_r; DRPLDT=ZERO_r
  !
  ! assign material parameters
  symm_class=int(PROPS(1))        ! c1
  m_SH=PROPS(2)                   ! c2
  algo=int(PROPS(3))              ! c3
  theta_KK=PROPS(4)               ! c4
  numtang=int(PROPS(5))           ! c5
  hard_law=int(PROPS(6))          ! c6
  orient=PROPS(7:9)               ! c7...c9
  elast_const=PROPS(10:12)        ! c10...c12
  mat_param_inelast=PROPS(13:25)  ! c13...c25
  !
  IF (algo.EQ.ZERO_int) THEN
     DTIME_tp=ONE_HUNDREDTH_r
  ELSE
     DTIME_tp=DTIME
  END IF
  ! compute rotation matrix and orthonormal subspace basis / projections
  SELECT CASE (symm_class)
  CASE (ZERO_int)
     R_init=IDENT
  CASE (ONE_int)
     R_init=IDENT
  CASE (TWO_int)
     R_init=ROTATION_MATRIX_BUNGE(orient)
  CASE DEFAULT
     WRITE(16,'(A)')'Unknown symmetry class in UMAT'
     WRITE(16,'(A)')'STOPPING COMPUTATION'
     STOP
  END SELECT
  !
  ! status flag
  IF (STATEV(NSTATV)<-eps_iter) THEN
     R_init=ROTATION_MATRIX_BUNGE(orient)
     CALL transform_tens_vec_form(TRANSPOSE(R_init),STATEV(1:9),&
          &ZERO_int,env_used)
     STATEV(NSTATV)=ONE_r
  END IF
  !
  ! compute determinant of deformation gradient at t_{n+1}
  J_n1 = comp_det_Matrix(DFGRD1)
  !
  ! compute orthonormal basis and corresponding subspace projectors
  CALL orthonormal_basis_subspace(orth_basis,orth_proj,symm_class,&
       &R_init)  
  !
  ! compute 2nd Piola-Kirchhoff stress and material tangent based on
  ! phenomenological single crystal plasticity model
  CALL crystal_plasticity_augm_lagr(S2PK,C_mat_tang,STATEV,DFGRD1,&
       &DTIME_tp,elast_const,mat_param_inelast,m_SH,symm_class,&
       &orth_basis,orth_proj,algo,numtang,NDI,NSHR,NSTATV,env_used,&
       &hard_law,theta_KK,conv_issue)
  !
  ! request a time increment change due to convergence issue at GP
  IF (conv_issue) THEN
     PNEWDT=ONE_FOURTH_r
  END IF
  !
  CALL transform_symm_tens_vec_form(S2PK,S_2PK_vec,ZERO_int,env_used)
  !
  ! post-processing of stress and moduli
  ! (push 2nd Piola-Kirchhoff to Kirchhoff (STRESS)
  ! and material to spatial tangent (DDSDDE),Lie derivative of Kirchhoff)
  CALL pushrc(STRESS,DDSDDE,S_2PK_vec,C_mat_tang,DFGRD1,TWO_int,env_used)
  !  
  SELECT CASE (env_used)
  CASE ('ABAQ')
     ! convert Kirchhoff stress -> Cauchy stress
     STRESS = STRESS / J_n1
     !
     ! convert spatial tangent (Lie derivative of Kirchhoff stress) ->
     ! Jaumann-rate of Kirchhoff stress (Nguyen, Waas, ZAMP, 2016(67),35)
     DDSDDE = DDSDDE / J_n1
     DDSDDE = DDSDDE + comp_add_stress_contrib_tang(STRESS)
     !
  CASE('FEAP')
     WRITE(*,'(A)')'Correct definition of stress and tangent to be handed'
     WRITE(*,'(A)')'over to FEAP needs to be elaborated'
     WRITE(*,'(A)')'STOPPING COMPUTATION'
     STOP
  CASE DEFAULT
     WRITE(16,'(A)')'Unknown value of env_used in UMAT'
     WRITE(16,'(A)')'STOPPING COMPUTATION'
     STOP
  END SELECT
  !
END SUBROUTINE UMAT
! #########################################################################
! ###################### Abaqus utility routines ##########################
! #########################################################################
! Schreiben von Ergebnissen in externe Datenbank Auskommentieren fuer DEBUG
SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
!
  !INCLUDE 'ABA_PARAM.INC'
  DIMENSION TIME(2),VNLEPSEQ_OUT(8)
  CHARACTER*256 JOBNAME,OUTDIR
  
  CALL GETJOBNAME(JOBNAME,LENJOBNAME)
  CALL GETOUTDIR(OUTDIR,LENOUTDIR)
! 
  IF(LOP.EQ.0) THEN
     ! OPEN(15,FILE=TRIM(OUTDIR)//'/'//TRIM(JOBNAME)//'.esv',&
     !      & STATUS='REPLACE')
     OPEN(16,FILE=TRIM(OUTDIR)//'/'//TRIM(JOBNAME)//'.feh',&
          & STATUS='REPLACE')
     ! OPEN(17,FILE=TRIM(OUTDIR)//'/'//TRIM(JOBNAME)//'.foo',&
     !      & STATUS='REPLACE')
     ! OPEN(18,FILE=TRIM(OUTDIR)//'/'//TRIM(JOBNAME)//'.bar',&
     !      & STATUS='REPLACE')
     ! OPEN(101,FILE=TRIM(OUTDIR)//'/'//TRIM(JOBNAME)//'.fuh',&
     !      & STATUS='REPLACE')
     ! OPEN(102,FILE=TRIM(OUTDIR)//'/'//TRIM(JOBNAME)//'.fet',&
     !      & STATUS='REPLACE')
     ! OPEN(103,FILE=TRIM(OUTDIR)//'/'//TRIM(JOBNAME)//'.fep',&
     !      & STATUS='REPLACE')
     ! OPEN(104,FILE=TRIM(OUTDIR)//'/'//TRIM(JOBNAME)//'.fol',&
     !      & STATUS='REPLACE')
     ! OPEN(105,FILE=TRIM(OUTDIR)//'/'//TRIM(JOBNAME)//'.fop',&
     !      & STATUS='REPLACE')
     ! OPEN(106,FILE=TRIM(OUTDIR)//'/'//TRIM(JOBNAME)//'.fos',&
     !      & STATUS='REPLACE')
  ELSEIF(LOP.EQ.3) THEN
     ! CLOSE(15)
     CLOSE(16)
     ! CLOSE(17)
     ! CLOSE(18)
     ! CLOSE(101)
     ! CLOSE(102)
     ! CLOSE(103)
     ! CLOSE(104)
     ! CLOSE(105)
     ! CLOSE(106)
  END IF

  RETURN
END SUBROUTINE UEXTERNALDB
! ###########################################################################
! Initialization of SDVs (solution-dependent variables)
! ###########################################################################
SUBROUTINE SDVINI(STATEV,COORDS,NSTATV,NCRDS,NOEL,NPT,LAYER,KSPT)
  !
  !
  ! insert " *INITIAL CONDITIONS,TYPE=SOLUTION,USER " in input file to initialize
  ! SDV's
  !
  USE Abaqus_Interface
  !
  IMPLICIT NONE
  INTEGER(KIND=ik)::NSTATV,NCRDS,NOEL,NPT,LAYER,KSPT
  REAL(KIND=rk),PARAMETER::ZERO_r=0.0,ONE_r=1.0
  REAL(KIND=rk),DIMENSION(NSTATV)::STATEV
  REAL(KIND=rk),DIMENSION(NCRDS)::COORDS
  !
  ! initialize regular sdv
  STATEV(1:NSTATV-1) = ZERO_r
  !
  ! initialize multi purpose flag
  STATEV(NSTATV)=-ONE_r
  RETURN
END SUBROUTINE SDVINI
! ###########################################################################
!
! #########################################################################
! SUBROUTINE TO DEFINE THE HEAT SOURCE AND MATERIAL TANGENT ENTRY
! #########################################################################
SUBROUTINE HETVAL(CMNAME,TEMP,TIME,DTIME,STATEV,FLUX,PREDEF,DPRED)
  !
  ! insert "*Heat Generation" in material definition in input file to
  ! activate this routine
  !
  USE Abaqus_Interface
  !
  IMPLICIT NONE
  !
  ! CHARACTERS (MATERIALNAME)
  CHARACTER(80)::CMNAME
  REAL(KIND=rk),DIMENSION(2),INTENT(IN)::TEMP
  REAL(KIND=rk),DIMENSION(2),INTENT(IN)::TIME
  REAL(KIND=rk),INTENT(IN)::DTIME
  REAL(KIND=rk),DIMENSION(2),INTENT(OUT)::FLUX
  ! CHANGE THE DIMENSION OF STATEV TO THE CONSIDERED CASE
  REAL(KIND=rk),dimension(38)::STATEV
  REAL(KIND=rk),dimension(:), allocatable::PREDEF
  REAL(KIND=rk),dimension(:), allocatable::DPRED
  !
  ! PLUG IN THE STATEV VARIABLE NUMBER OF YOUR STORED DERIVATIVE OF SOURCE
  ! TERM W.R.T. TEMPERATURE (NONLOCAL VARIABLE)
  !
  FLUX(1) 	= 0.0_rk
  !
  FLUX(2) 	= 0.0_rk
  RETURN
END SUBROUTINE HETVAL
! #########################################################################
