PROGRAM tabulation_example
    
    USE hoppet_v1
    
    IMPLICIT NONE
    
    ! Declare variables
    
    INTEGER                :: order, nloop, ix
    REAL(dp)               :: quark_masses(4:6), pdf_at_xQ(-6:6)
    REAL(dp)               :: dy, ymax, Q0, Q
    REAL(dp)               :: heralhc_xvals(200)
    REAL(dp), POINTER      :: pdf0(:, :)
    
    TYPE(grid_def)         :: grid, gdarray(4)
    TYPE(dglap_holder)     :: dh
    TYPE(pdf_table)        :: table
    TYPE(running_coupling) :: coupling
    
    ! Read in the x values that the evolved PDF will be returned at
    ! heralhc_xvals.txt is created automatically from the R code
    
    OPEN (20, FILE = "heralhc_xvals.txt", STATUS = 'old', ACCESS = 'sequential', FORM = 'formatted', ACTION = 'read')
    READ (20, *) heralhc_xvals
    CLOSE(20)
    
    ! Initialise the variables declared above
    
    order             = -6
    nloop             = 3
    quark_masses(4:6) = (/1.414213563_dp, 4.5_dp, 175.0_dp/)
    ymax              = 13.9_dp
    dy                = 0.1_dp
    
    ! Read in the initial energy scale and the scale that is being evolved to
    ! mu.txt is created automatically from the R code
    ! mu in the R code is called Q in this code
    
    OPEN (21, FILE = "mu.txt", STATUS = 'old', ACCESS = 'sequential', FORM = 'formatted', ACTION = 'read')
    READ (21, *) Q0, Q
    CLOSE(21)
    
    ! Create the grid of x values that will be used to perform the evolution
    
    CALL InitGridDef(gdarray(4), dy/27.0_dp, 0.2_dp, order = order)
    CALL InitGridDef(gdarray(3), dy/9.0_dp,  0.5_dp, order = order)
    CALL InitGridDef(gdarray(2), dy/3.0_dp,  2.0_dp, order = order)
    CALL InitGridDef(gdarray(1), dy,         ymax,   order = order)
    CALL InitGridDef(grid, gdarray(1:4), locked = .TRUE.)
    
    ! InitDglapHolder contains information about convolution using splitting functions
    ! nflo and nfhi are limits for the n_f value
    ! For gluons only, set both to 0
    
    CALL InitDglapHolder(grid, dh, factscheme = factscheme_MSbar, nloop = nloop, nflo = 0, nfhi = 0)
    
    ! Ask for the evolution to be done
    
    CALL AllocPDF(grid, pdf0)
    
    pdf0 = unpolarized_dummy_pdf(xValues(grid))
    
    CALL InitRunningCoupling(coupling, alfas = 0.35_dp, Q = Q0, nloop = nloop, quark_masses = quark_masses)
    CALL AllocPdfTable(grid, table, Qmin = 1.0_dp, Qmax = 10000.0_dp, dlnlnQ = dy/4.0_dp, freeze_at_Qmin = .TRUE.)
    CALL AddNfInfoToPdfTable(table, coupling)
    CALL EvolvePdfTable(table, Q0, pdf0, dh, coupling, nloop = nloop)
    
    ! Print the results line-by-line for each of the x values defined in heralhc_xvals
    
    DO ix = 1, SIZE(heralhc_xvals)
        
        CALL EvalPdfTable_xQ(table, heralhc_xvals(ix), Q, pdf_at_xQ)
        WRITE(6, '(ES11.5)') pdf_at_xQ(0)
        
    END DO
    
    ! Clean the grid
    
    CALL Delete(table)
    CALL Delete(pdf0)
    CALL Delete(dh)
    CALL Delete(coupling)
    CALL Delete(grid)
    
CONTAINS
    
    ! The function that calculates the evoltion
    
    FUNCTION unpolarized_dummy_pdf(xvals) RESULT(pdf)
        
        ! Declare the variables for use in this function
        
        REAL(dp), INTENT(IN) :: xvals(:)
        REAL(dp)             :: pdf(SIZE(xvals), ncompmin:ncompmax)
        REAL(dp)             :: uv(SIZE(xvals)), dv(SIZE(xvals))
        REAL(dp)             :: ubar(SIZE(xvals)), dbar(SIZE(xvals))
        REAL(dp)             :: A_g, lambda_g
        
        pdf = zero
        
        CALL LabelPdfAsHuman(pdf)
        
        ! Read in the parameters for the initial PDF
        
        OPEN (22, FILE = "parameters.txt", STATUS = 'old', ACCESS = 'sequential', FORM = 'formatted', ACTION = 'read')
        READ (22, *) A_g, lambda_g
        CLOSE(22)
        
        ! This initial form comes from Amir, R. et al. (2012) (arXiv:1212.2974v1)
        ! The parameters A_g and lambda_g are read in from the parameters.txt file above
        
        pdf(:, iflv_g) = A_g * xvals**(-lambda_g) * (1 - xvals)**5.6
        
        ! Set initial quarks to 0
        
        pdf(:,-iflv_s) = 0
        pdf(:, iflv_s) = 0
        pdf(:, iflv_u) = 0
        pdf(:,-iflv_u) = 0
        pdf(:, iflv_d) = 0
        pdf(:,-iflv_d) = 0
        
    END FUNCTION unpolarized_dummy_pdf
    
END PROGRAM tabulation_example
