    module DarkEnergyBarotropic
    ! Scalar field dark energy with generalized linear barotropic equation of state
    !   P_de = ca2*rho_de + C,
    ! with constant adiabatic sound speed ca2 = dP/drho and constant C, as in arXiv:1008.1943 (quintessence)
    ! and arXiv:1206.5194 (phantom). The model is parameterized by the present-day equation of state w0 = w(a=1)
    ! and ca2, giving (Eqs. 14, 15 of arXiv:1008.1943)
    !   w(a) = (1+ca2)(1+w0)/(1 + w0 - (w0 - ca2) a^(3(1+ca2))) - 1,
    !   rho_de(a)/rho_de(1) = ((1+w0) a^(-3(1+ca2)) + ca2 - w0)/(1+ca2),
    ! so that dw/dlna = 3(1+w)(w - ca2) (Eq. 12). For -1 < w0 <= 0, w evolves monotonically from ca2 at early times
    ! to w0 today. For w0 < -1 (phantom), ca2 <= w0 is required for the density to be positive at all a < 1; w then
    ! evolves from -1 at early times to w0 today (tending to ca2 in the future). w never crosses -1. ca2 = w0 gives
    ! constant w and w0 = -1 a cosmological constant.
    ! Perturbations are those of a fluid with constant rest-frame sound speed cs2_lam (=1 for canonical and phantom
    ! scalar fields), as in TDarkEnergyFluid, but with the exact dw/dlna for this equation of state.
    use DarkEnergyInterface
    use classes
    implicit none

    private

    type, extends(TDarkEnergyModel) :: TDarkEnergyBarotropic
        real(dl) :: w_lam = -1._dl ! w0 = w(a=1)
        real(dl) :: ca2 = -1._dl ! constant adiabatic sound speed squared (asymptotic value of w)
        real(dl) :: cs2_lam = 1._dl ! constant rest-frame sound speed squared
        real(dl), private :: pow = 0._dl ! 3(1+ca2), cached by Init
        real(dl), private :: one_plus_w0 = 0._dl ! 1 + w0, cached by Init
    contains
    procedure :: ReadParams => TDarkEnergyBarotropic_ReadParams
    procedure, nopass :: PythonClass => TDarkEnergyBarotropic_PythonClass
    procedure, nopass :: SelfPointer => TDarkEnergyBarotropic_SelfPointer
    procedure :: Init => TDarkEnergyBarotropic_Init
    procedure :: PrintFeedback => TDarkEnergyBarotropic_PrintFeedback
    procedure :: w_de => TDarkEnergyBarotropic_w_de
    procedure :: grho_de => TDarkEnergyBarotropic_grho_de
    procedure :: Effective_w_wa => TDarkEnergyBarotropic_Effective_w_wa
    procedure :: assume_scale_indep_lowz_growth => TDarkEnergyBarotropic_assume_scale_indep_lowz_growth
    procedure :: PerturbedStressEnergy => TDarkEnergyBarotropic_PerturbedStressEnergy
    procedure :: PerturbationEvolve => TDarkEnergyBarotropic_PerturbationEvolve
    end type TDarkEnergyBarotropic

    public TDarkEnergyBarotropic

    contains

    function expm1_over_x(x) result(res)
    ! (exp(x) - 1)/x, accurate also for small x (so results are smooth through ca2 = -1)
    real(dl), intent(in) :: x
    real(dl) :: res, u

    u = exp(x)
    if (u == 1._dl) then
        res = 1._dl
    else if (u == 0._dl) then
        res = -1._dl/x
    else
        res = (u - 1._dl)/log(u)
    end if

    end function expm1_over_x

    subroutine TDarkEnergyBarotropic_ReadParams(this, Ini)
    use IniObjects
    class(TDarkEnergyBarotropic) :: this
    class(TIniFile), intent(in) :: Ini

    call this%TDarkEnergyModel%ReadParams(Ini)
    if (Ini%HasKey('wa')) error stop 'Barotropic dark energy model uses ca2 rather than wa'
    this%w_lam = Ini%Read_Double('w', -1.d0)
    this%ca2 = Ini%Read_Double('ca2')
    this%cs2_lam = Ini%Read_Double('cs2_lam', 1.d0)

    end subroutine TDarkEnergyBarotropic_ReadParams

    function TDarkEnergyBarotropic_PythonClass()
    character(len=:), allocatable :: TDarkEnergyBarotropic_PythonClass

    TDarkEnergyBarotropic_PythonClass = 'DarkEnergyBarotropic'

    end function TDarkEnergyBarotropic_PythonClass

    subroutine TDarkEnergyBarotropic_SelfPointer(cptr, P)
    use iso_c_binding
    type(c_ptr) :: cptr
    type(TDarkEnergyBarotropic), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TDarkEnergyBarotropic_SelfPointer

    subroutine TDarkEnergyBarotropic_Init(this, State)
    use classes
    class(TDarkEnergyBarotropic), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State

    this%one_plus_w0 = 1._dl + this%w_lam
    this%pow = 3._dl*(1._dl + this%ca2)
    ! w0 = -1 gives w = -1 at all times for any ca2
    this%is_cosmological_constant = abs(this%one_plus_w0) < 1.e-6_dl
    if (this%is_cosmological_constant) then
        this%num_perturb_equations = 0
    else
        if (this%one_plus_w0 > 0) then
            ! trap dark energy with w>0 today or at early times
            if (max(this%w_lam, this%ca2) > 0) &
                error stop 'Barotropic dark energy: need w <= 0 and ca2 <= 0, otherwise w > 0 today or at early times'
        else if (this%ca2 > this%w_lam) then
            ! otherwise the density changes sign at some a < 1
            error stop 'Barotropic dark energy: phantom w < -1 requires ca2 <= w'
        end if
        this%num_perturb_equations = 2
    end if

    end subroutine TDarkEnergyBarotropic_Init

    subroutine TDarkEnergyBarotropic_PrintFeedback(this, FeedbackLevel)
    class(TDarkEnergyBarotropic) :: this
    integer, intent(in) :: FeedbackLevel

    if (FeedbackLevel > 0) write(*, '("(w0, ca2) = (", f8.5,", ", f8.5, ")")') &
        this%w_lam, this%ca2

    end subroutine TDarkEnergyBarotropic_PrintFeedback

    function TDarkEnergyBarotropic_w_de(this, a)
    class(TDarkEnergyBarotropic) :: this
    real(dl) :: TDarkEnergyBarotropic_w_de
    real(dl), intent(in) :: a
    real(dl) :: loga, y

    ! 1 + w = (1+w0)/(a^pow + 3(1+w0)(1 - a^pow)/pow), written in forms that do not overflow and are smooth
    ! through pow = 0 (ca2 = -1), where 1 + w = (1+w0)/(1 - 3(1+w0) ln a)
    loga = log(a)
    y = this%pow*loga
    if (y <= 0) then
        TDarkEnergyBarotropic_w_de = this%one_plus_w0/(exp(y) - 3*this%one_plus_w0*loga*expm1_over_x(y)) - 1
    else
        TDarkEnergyBarotropic_w_de = this%one_plus_w0*exp(-y)/(1 - 3*this%one_plus_w0*loga*expm1_over_x(-y)) - 1
    end if

    end function TDarkEnergyBarotropic_w_de

    function TDarkEnergyBarotropic_grho_de(this, a) ! relative density (8 pi G a^4 rho_de /grhov)
    class(TDarkEnergyBarotropic) :: this
    real(dl) :: TDarkEnergyBarotropic_grho_de
    real(dl), intent(in) :: a
    real(dl) :: loga

    if (a == 0._dl) then
        TDarkEnergyBarotropic_grho_de = 0._dl ! a^4 rho_de -> 0 for ca2 < 1/3
    else
        ! rho_de(a)/rho_de(1) = 1 + 3(1+w0)(a^(-pow) - 1)/pow = 1 - 3(1+w0) ln(a) (exp(x)-1)/x, x = -pow*ln(a)
        loga = log(a)
        TDarkEnergyBarotropic_grho_de = a**4*(1 - 3*this%one_plus_w0*loga*expm1_over_x(-this%pow*loga))
    end if

    end function TDarkEnergyBarotropic_grho_de

    subroutine TDarkEnergyBarotropic_Effective_w_wa(this, w, wa)
    ! w(a) ~ w + (1-a) wa near a=1, used e.g. by halofit fits
    class(TDarkEnergyBarotropic), intent(inout) :: this
    real(dl), intent(out) :: w, wa

    w = this%w_lam
    wa = -3*(1 + this%w_lam)*(this%w_lam - this%ca2)

    end subroutine TDarkEnergyBarotropic_Effective_w_wa

    function TDarkEnergyBarotropic_assume_scale_indep_lowz_growth(this)
    ! As TDarkEnergyEqnOfState: with cs2=1 dark energy does not cluster on sub-horizon scales
    class(TDarkEnergyBarotropic) :: this
    logical :: TDarkEnergyBarotropic_assume_scale_indep_lowz_growth

    TDarkEnergyBarotropic_assume_scale_indep_lowz_growth = this%cs2_lam > 0.99_dl

    end function TDarkEnergyBarotropic_assume_scale_indep_lowz_growth

    subroutine TDarkEnergyBarotropic_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TDarkEnergyBarotropic), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) :: a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix

    dgrhoe = ay(w_ix)*grhov_t
    dgqe = ay(w_ix + 1)*grhov_t*(1 + w)

    end subroutine TDarkEnergyBarotropic_PerturbedStressEnergy

    subroutine TDarkEnergyBarotropic_PerturbationEvolve(this, ayprime, w, w_ix, &
        a, adotoa, k, z, y)
    ! Same fluid equations as TDarkEnergyFluid (variables delta and v), with dw/dlna = 3(1+w)(w-ca2)
    class(TDarkEnergyBarotropic), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a, adotoa, w, k, z, y(:)
    integer, intent(in) :: w_ix
    real(dl) Hv3_over_k, one_plus_w

    one_plus_w = 1 + w
    Hv3_over_k = 3*adotoa*y(w_ix + 1)/k
    ! density perturbation
    ayprime(w_ix) = -3*adotoa*(this%cs2_lam - w)*(y(w_ix) + one_plus_w*Hv3_over_k) &
        - one_plus_w*k*(y(w_ix + 1) + z) &
        - 3*adotoa*one_plus_w*(w - this%ca2)*Hv3_over_k
    ! velocity
    if (abs(one_plus_w) > 1e-6) then
        ayprime(w_ix + 1) = -adotoa*(1 - 3*this%cs2_lam)*y(w_ix + 1) + &
            k*this%cs2_lam*y(w_ix)/one_plus_w
    else
        ayprime(w_ix + 1) = 0
    end if

    end subroutine TDarkEnergyBarotropic_PerturbationEvolve

    end module DarkEnergyBarotropic
