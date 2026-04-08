!===============================================================================
! CODE ACCESSIBILITY - eNeuro Submission
! 
! Title: Macroscopic Description for Networks of Spiking Neurons
! Description: 
!   This program simulates the ASYMMETRIC paradigm of the neural mass model and
!   spiking network. It includes:
!   1. Asymmetric Short-Term Plasticity (STP applies to I-neurons and E-to-I 
!      connections, but NOT to recurrent E-to-E connections).
!   2. Top-down Theta (theta) rhythm forcing to generate Phase-Amplitude Coupling.
!   3. A specific constant current offset (0.03*sqrt(K)) unique to this regime.
!===============================================================================

MODULE PARAMETERS
  IMPLICIT NONE
  
  ! --- Network Size ---
  INTEGER(8), PARAMETER :: n = 30000       ! Total number of neurons
  INTEGER(8), PARAMETER :: ne = 25000      ! Number of Excitatory neurons
  INTEGER(8), PARAMETER :: ni = 5000      ! Number of Inhibitory neurons

  ! --- Time and Mathematical Constants ---
  REAL(8), PARAMETER :: pi = 3.14159265358979303846d0
  REAL(8), PARAMETER :: t_total = 1000.0d0  ! Total simulation time
  REAL(8), PARAMETER :: dt = 0.0001d0       ! Integration time step
  
  ! --- Neuron Membrane Parameters (Quadratic Integrate-and-Fire) ---
  REAL(8), PARAMETER :: v_p = 100.0d0       ! Peak membrane potential (Spike threshold)
  REAL(8), PARAMETER :: v_r = -100.0d0      ! Reset membrane potential
  REAL(8), PARAMETER :: t_refractory = 2.0d0 / v_p ! Refractory period
  REAL(8), PARAMETER :: t_syn_delay = 1.0d0 / v_p  ! Synaptic transmission delay

  ! --- Synaptic Coupling Strengths ---
  REAL(8), PARAMETER :: f0 = 1.0d0
  REAL(8), PARAMETER :: big_jee = 0.27d0 * f0
  REAL(8), PARAMETER :: big_jii = 0.953939d0 * f0
  REAL(8), PARAMETER :: big_jie = 0.3d0 * f0
  REAL(8), PARAMETER :: big_jei = 0.96286d0 * f0
  
  ! --- Short-Term Plasticity (STP) Parameters ---
  REAL(8), PARAMETER :: U_stp = 0.5d0       ! Utilization parameter (user's 'a=0.5')
  REAL(8), PARAMETER :: tau_d = 2.0d0       ! Recovery time constant for synaptic depression

  ! --- Top-Down Theta Forcing Parameters ---
  REAL(8), PARAMETER :: f_theta = 0.05d0    ! Theta rhythm frequency (5 Hz)
  REAL(8), PARAMETER :: t_stim_onset = 200.0d0 ! Time when theta forcing begins
  
  ! --- Data Measurement Windows ---
  REAL(8), PARAMETER :: delta_t1 = 0.02d0   ! Sliding window for macroscopic rate calculation
  REAL(8), PARAMETER :: t_infinity = 1.0d6  ! Large number to prevent continuous spiking
  REAL(8), PARAMETER :: t_record_start = MAX(0.0d0, t_total - 3000.0d0) ! Recording window

  ! --- Global Macroscopic (Mean-Field) Variables ---
  REAL(8) :: th_re, th_ve, th_ri, th_vi, th_xe, th_xi
  REAL(8) :: delta_0ii, delta_0ee, k0
  
  ! --- Global Microscopic Variables ---
  REAL(8), DIMENSION(n) :: v, spike_count, t_last_spike, t_syn_release, isi_sum, isi_sq_sum, num_isi, mean_x
  INTEGER(8), DIMENSION(:), ALLOCATABLE   :: ge2, gi2, gei2, gie2
  
  REAL(8) :: t, t_window, r_all, r_alle, r_alli, v_all, v_alle, v_alli
  REAL(8) :: delta_kee, delta_kii, j0ii, j0ee, j0ie, j0ei, big_jie0, big_jei0
  INTEGER(8) :: ge1(ne), gi1(ni), sumge1(ne+1), sumgi1(ni+1), in_degree(n), in_degree_cross(n)
  INTEGER(8) :: gie1(ne), gei1(ni), sumgie1(ne+1), sumgei1(ni+1)

END MODULE PARAMETERS

!===============================================================================
! MAIN PROGRAM
!===============================================================================
PROGRAM main_simulation_asymmetric
  USE PARAMETERS
  IMPLICIT NONE
  
  ! --- Local Variables ---
  INTEGER :: file_idx
  INTEGER(8) :: i, j, p, idx_target, l
  INTEGER(8) :: count_1, count_2
  REAL(8) :: rn, rn0, rand_val, l0, m0
  REAL(8) :: v_mean_temp, v2_mean_temp, v_i_mean(n), v_i2_mean(n)
  REAL(8) :: var_mean_v, mean_var_v(n), sync_rho, count_samples
  REAL(8) :: d_re, d_ve, d_ri, d_vi, start_time, stop_time, d_xe, d_xi
  REAL(8) :: count_active, count_active_e, count_active_i
  REAL(8) :: cv_mean, cv_e, cv_i, count_valid_cv, count_cve, count_cvi
  REAL(8) :: current_isi, local_std, local_cv
  
  ! Dynamic current variables
  REAL(8) :: i_0e_base, i_0i_base, ie_current, ii_current, theta_forcing_signal
  
  CHARACTER(len=20) :: cite, file_dk, file_time, file_theory, file_v0e, file_v0i, file_isse, file_issi, file_cv
  INTEGER(8), DIMENSION(:,:), ALLOCATABLE :: ge, gi, gei, gie

  CALL CPU_TIME(start_time)

  ! --- Parameter Initialization ---
  file_idx = 500
  k0 = 1000.0d0
  delta_0ee = 2.0d0 / f0
  delta_0ii = 0.3d0 / f0

  ! --- File Operations ---
  WRITE(cite, "(i8)") file_idx
  file_dk = 'adk' // TRIM(ADJUSTL(cite)) // '.dat'
  file_time = 'atime' // TRIM(ADJUSTL(cite)) // '.dat'
  file_theory = 'atheory' // TRIM(ADJUSTL(cite)) // '.dat'
  file_v0e = 'av0e' // TRIM(ADJUSTL(cite)) // '.dat'
  file_v0i = 'av0i' // TRIM(ADJUSTL(cite)) // '.dat'
  file_isse = 'aisse' // TRIM(ADJUSTL(cite)) // '.dat'
  file_issi = 'aissi' // TRIM(ADJUSTL(cite)) // '.dat'
  file_cv = 'acv' // TRIM(ADJUSTL(cite)) // '.dat'

  OPEN(100, file=file_dk, ACTION='WRITE')
  OPEN(200, file=file_time, ACTION='WRITE')
  OPEN(210, file=file_theory, ACTION='WRITE')
  OPEN(300, file=file_v0e, ACTION='WRITE')
  OPEN(310, file=file_v0i, ACTION='WRITE')
  OPEN(400, file=file_isse, ACTION='WRITE')
  OPEN(700, file=file_cv, ACTION='WRITE')
  OPEN(800, file=file_issi, ACTION='WRITE')

  ! --- Scaling currents and couplings by sqrt(K) ---
  delta_kee = delta_0ee * SQRT(k0)
  delta_kii = delta_0ii * SQRT(k0)
  big_jie0 = big_jie * k0 / DBLE(k0)
  big_jei0 = big_jei * k0 / DBLE(k0)
  j0ee = big_jee * SQRT(k0)
  j0ii = -big_jii * SQRT(k0)
  j0ie = big_jie0 * SQRT(k0)
  j0ei = -big_jei0 * SQRT(k0)

  ALLOCATE (ge(ne, ne), gi(ni, ni), gie(ni, ne), gei(ne, ni))
  ge = 0; gi = 0; gei = 0; gie = 0

  !=============================================================================
  ! NETWORK TOPOLOGY GENERATION 
  !=============================================================================
  CALL RANDOM_SEED()
  
  ! Build Excitatory-to-Excitatory (E-E) Network
  DO i = 1, ne
     CALL RANDOM_NUMBER(rand_val)
     l0 = k0 + delta_kee * TAN(0.5d0 * pi * 2.0d0 * (0.5d0 - rand_val))
     l = CEILING(l0)
     DO WHILE(l < 1 .OR. l > ne)
        CALL RANDOM_NUMBER(rand_val)
        l0 = k0 + delta_kee * TAN(0.5d0 * pi * 2.0d0 * (0.5d0 - rand_val))
        l = CEILING(l0)
     ENDDO
     IF(l <= ne) THEN
        DO j = 1, l
           idx_target = -1
           DO WHILE(idx_target < 0)
              DO p = 1, 10
                 CALL RANDOM_NUMBER(rn0)
                 idx_target = CEILING(ne * rn0)
                 IF(ge(i, idx_target) == 0 .AND. idx_target /= i) EXIT
              ENDDO
           ENDDO
           ge(i, idx_target) = 1
        ENDDO
     ENDIF
  ENDDO

  ! Build Inhibitory-to-Inhibitory (I-I) Network
  DO i = 1, ni
     CALL RANDOM_NUMBER(rand_val)
     l0 = k0 + delta_kii * TAN(0.5d0 * pi * 2.0d0 * (0.5d0 - rand_val))
     l = CEILING(l0)
     DO WHILE(l < 1 .OR. l > ni)
        CALL RANDOM_NUMBER(rand_val)
        l0 = k0 + delta_kii * TAN(0.5d0 * pi * 2.0d0 * (0.5d0 - rand_val))
        l = CEILING(l0)
     ENDDO
     IF(l <= ni) THEN
        DO j = 1, l
           idx_target = -1
           DO WHILE(idx_target < 0)
              DO p = 1, 10
                 CALL RANDOM_NUMBER(rn0)
                 idx_target = CEILING(ni * rn0)
                 IF(gi(i, idx_target) == 0 .AND. idx_target /= i) EXIT
              ENDDO
           ENDDO
           gi(i, idx_target) = 1
        ENDDO
     ENDIF
  ENDDO

  ! Build Cross-Connections (Uniform Random)
  DO i = 1, ni
     DO j = 1, ne
        CALL RANDOM_NUMBER(rand_val)
        IF(rand_val < k0 / ni .AND. gie(i, j) == 0) gie(i, j) = 1
     ENDDO
  ENDDO
  
  DO i = 1, ne
     DO j = 1, ni
       CALL RANDOM_NUMBER(rand_val)
       IF(rand_val < k0 / ne .AND. gei(i, j) == 0) gei(i, j) = 1
     ENDDO
  ENDDO

  ! Convert to Adjacency Lists
  ALLOCATE (ge2(SUM(ge)), gi2(SUM(gi)))
  sumge1 = 0; sumgi1 = 0; count_2 = 0
  DO i = 1, ne
     count_1 = 0
     DO j = 1, ne
        IF(ge(j, i) == 1) THEN
           count_1 = count_1 + 1
           count_2 = count_2 + 1
           ge2(count_2) = j
        ENDIF
     ENDDO
     ge1(i) = count_1
     sumge1(i+1) = count_2
  ENDDO

  count_2 = 0
  DO i = 1, ni
     count_1 = 0
     DO j = 1, ni
        IF(gi(j, i) == 1) THEN
           count_1 = count_1 + 1
           count_2 = count_2 + 1
           gi2(count_2) = j
        ENDIF
     ENDDO
     gi1(i) = count_1
     sumgi1(i+1) = count_2
  ENDDO

  ALLOCATE (gie2(SUM(gie)), gei2(SUM(gei)))
  sumgie1 = 0; sumgei1 = 0; count_2 = 0
  DO i = 1, ne
     count_1 = 0
     DO j = 1, ni
        IF(gie(j, i) == 1) THEN
           count_1 = count_1 + 1
           count_2 = count_2 + 1
           gie2(count_2) = j
        ENDIF
     ENDDO
     gie1(i) = count_1
     sumgie1(i+1) = count_2
  ENDDO

  count_2 = 0
  DO i = 1, ni
     count_1 = 0
     DO j = 1, ne
        IF(gei(j, i) == 1) THEN
           count_1 = count_1 + 1
           count_2 = count_2 + 1
           gei2(count_2) = j
        ENDIF
     ENDDO
     gei1(i) = count_1
     sumgei1(i+1) = count_2
  ENDDO

  ! Output in-degrees
  in_degree = 0; in_degree_cross = 0
  DO i = 1, ne
     in_degree(i) = SUM(ge(i, :))
     in_degree_cross(i) = SUM(gei(i, :))
  ENDDO
  DO i = 1, ni
     in_degree(i+ne) = SUM(gi(i, :))
     in_degree_cross(i+ne) = SUM(gie(i, :))
  ENDDO
  DO i = 1, n
     WRITE(100, *) i, in_degree(i), in_degree_cross(i)
  ENDDO

  DEALLOCATE (ge, gi, gei, gie)

  !=============================================================================
  ! INITIAL CONDITIONS
  !=============================================================================
  CALL RANDOM_SEED()
  DO i = 1, n
     CALL RANDOM_NUMBER(rn)
     v(i) = 80.0d0 * (0.5d0 - rn)
     CALL RANDOM_NUMBER(rand_val)
     mean_x(i) = 0.02d0 * rand_val
  ENDDO

  th_re = 0.10d0; th_ve = -1.0d0; th_xe = 1.0d0
  th_ri = 0.10d0; th_vi = -1.0d0; th_xi = 1.0d0
  
  spike_count = 0.0d0; t_last_spike = 0.0d0; t_syn_release = t_infinity
  t = 0.0d0; t_window = 0.0d0; r_all = 0.0d0
  
  isi_sum = 0.0d0; isi_sq_sum = 0.0d0; num_isi = 0.0d0
  v_mean_temp = 0.0d0; v2_mean_temp = 0.0d0; v_i_mean = 0.0d0; v_i2_mean = 0.0d0
  count_samples = 0.0d0; var_mean_v = 0.0d0; mean_var_v = 0.0d0
  
  sync_rho = 0.0d0; count_valid_cv = 0.0d0; cv_mean = 0.0d0
  count_cve = 0.0d0; count_cvi = 0.0d0; cv_e = 0.0d0; cv_i = 0.0d0

  !=============================================================================
  ! TIME EVOLUTION (Euler Integration)
  !=============================================================================
  DO WHILE(t <= t_total)
     t = t + dt
     
     ! --- External Current Logic (Theta Rhythm Forcing) ---
     IF (t < t_stim_onset) THEN
        i_0e_base = 0.1d0
        i_0i_base = i_0e_base / 1.02d0
        ie_current = i_0e_base * SQRT(k0)
        ii_current = i_0i_base * SQRT(k0)
        theta_forcing_signal = i_0e_base * SQRT(k0) * (1.0d0 - COS(2.0d0 * pi * f_theta * t)) ! Just for output log
     ELSE
        i_0e_base = 0.07d0
        i_0i_base = i_0e_base / 1.02d0
        ie_current = (i_0e_base * SQRT(k0)) * (1.0d0 - COS(2.0d0 * pi * f_theta * t))
        ii_current = (i_0i_base * SQRT(k0)) * (1.0d0 - COS(2.0d0 * pi * f_theta * t))
        theta_forcing_signal = ie_current
     ENDIF

     d_re = th_re; d_ve = th_ve; d_xe = th_xe
     d_ri = th_ri; d_vi = th_vi; d_xi = th_xi

     ! --- 1. Macroscopic Mean-Field Dynamics (ASYMMETRIC Paradigm) ---
     ! BIOLOGICAL ASYMMETRY: The E-to-E connections are NOT modulated by d_xe.
     th_re = th_re + (2.0d0 * d_re * d_ve + delta_0ee * big_jee * d_re / pi) * dt
     th_ve = th_ve + (d_ve**2 + ie_current + 0.03d0 * SQRT(k0) - (pi * d_re)**2 + j0ee * d_re - big_jei * d_xi * d_ri * SQRT(k0)) * dt
     th_ri = th_ri + (2.0d0 * d_ri * d_vi + delta_0ii * big_jii * d_xi * d_ri / pi) * dt
     th_vi = th_vi + (d_vi**2 + ii_current + (0.03d0 / 1.02d0) * SQRT(k0) - (pi * d_ri)**2 + j0ii * d_xi * d_ri + big_jie * d_xe * d_re * SQRT(k0)) * dt
     
     th_xe = th_xe + ((1.0d0 - d_xe) / tau_d - 0.8d0 * U_stp * d_xe * d_re) * dt
     th_xi = th_xi + ((1.0d0 - d_xi) / tau_d - U_stp * d_xi * d_ri) * dt

     ! Calculate mean voltage excluding refractory neurons
     count_active = 0.0d0; count_active_e = 0.0d0; count_active_i = 0.0d0
     v_all = 0.0d0; v_alle = 0.0d0; v_alli = 0.0d0
     
     DO i = 1, n
        IF(t - t_last_spike(i) > t_refractory) THEN
           v_all = v_all + v(i)
           count_active = count_active + 1.0d0
           IF(i <= ne) THEN
              v_alle = v_alle + v(i)
              count_active_e = count_active_e + 1.0d0
           ELSE
              v_alli = v_alli + v(i)
              count_active_i = count_active_i + 1.0d0
           ENDIF
        ENDIF
     ENDDO
     v_all  = v_all / count_active
     v_alle = v_alle / count_active_e
     v_alli = v_alli / count_active_i

     ! --- Sliding window data recording ---
     IF(t > t_window + delta_t1) THEN
        t_window = t
        r_all  = SUM(spike_count) / (n * delta_t1)
        r_alle = SUM(spike_count(1:ne)) / (ne * delta_t1)
        r_alli = SUM(spike_count(ne+1:n)) / (ni * delta_t1)
        
        IF(t > t_record_start) THEN
           ! Added an 8th column to log the exact theta forcing amplitude
           WRITE(200, "(8f16.8)") t, r_all, r_alle, r_alli, v_all, v_alle, v_alli, theta_forcing_signal
           WRITE(210, "(7f16.8)") t, th_re, th_ri, th_ve, th_vi, th_xe, th_xi
           
           count_samples = count_samples + 1.0d0
           v_mean_temp = v_mean_temp + v_all
           v2_mean_temp = v2_mean_temp + v_all * v_all
           DO i = 1, n
              v_i_mean(i) = v_i_mean(i) + v(i)
              v_i2_mean(i) = v_i2_mean(i) + v(i) * v(i)
           ENDDO
        ENDIF
        spike_count = 0.0d0
     ENDIF

     ! --- 2. Microscopic Spiking Dynamics (QIF Neurons) ---
     DO i = 1, ne
        IF(t - t_last_spike(i) >= t_refractory) THEN
           v(i) = v(i) + (v(i)**2 + ie_current + 0.03d0 * SQRT(k0)) * dt
        ENDIF
        mean_x(i) = mean_x(i) + (1.0d0 - mean_x(i)) / tau_d * dt
     ENDDO
     
     DO i = ne+1, n
        IF(t - t_last_spike(i) >= t_refractory) THEN
           v(i) = v(i) + (v(i)**2 + ii_current + (0.03d0 / 1.02d0) * SQRT(k0)) * dt
        ENDIF
        mean_x(i) = mean_x(i) + (1.0d0 - mean_x(i)) / tau_d * dt
     ENDDO

     ! --- 3. Spike Detection and Reset ---
     DO i = 1, n
        IF(v(i) >= v_p) THEN
           IF(t > t_record_start) THEN
              IF(i <= ne) WRITE(300, "(f16.8,i16)") t, i
              IF(i >  ne) WRITE(310, "(f16.8,i16)") t, i
              
              current_isi = t - t_last_spike(i)
              isi_sum(i) = isi_sum(i) + current_isi
              isi_sq_sum(i) = isi_sq_sum(i) + current_isi**2
              num_isi(i) = num_isi(i) + 1.0d0
           ENDIF
           
           t_last_spike(i) = t
           t_syn_release(i) = t
           spike_count(i) = spike_count(i) + 1.0d0
           v(i) = v_r
        ENDIF
     ENDDO

     ! --- 4. Synaptic Transmission (ASYMMETRIC STP) ---
     ! BIOLOGICAL ASYMMETRY: The E-to-E transmission is NOT modulated by mean_x.
     DO i = 1, ne
        IF((t - t_syn_release(i)) > t_syn_delay) THEN
           t_syn_release(i) = t_infinity
           DO j = sumge1(i)+1, sumge1(i) + ge1(i)
              v(ge2(j)) = v(ge2(j)) + j0ee / k0           ! <-- NO STP (mean_x) applied here
           ENDDO
           DO j = sumgie1(i)+1, sumgie1(i) + gie1(i)
              v(gie2(j)+ne) = v(gie2(j)+ne) + j0ie * mean_x(i) / k0  ! <-- STP applied here
           ENDDO
           mean_x(i) = mean_x(i) - 0.8d0 * U_stp * mean_x(i)
        ENDIF
     ENDDO
     
     ! Inhibitory transmission (STP applied to all targets)
     DO i = ne+1, n
        IF((t - t_syn_release(i)) > t_syn_delay) THEN
           t_syn_release(i) = t_infinity
           DO j = sumgi1(i-ne)+1, sumgi1(i-ne) + gi1(i-ne)
              v(gi2(j)+ne) = v(gi2(j)+ne) + j0ii * mean_x(i) / k0
           ENDDO
           DO j = sumgei1(i-ne)+1, sumgei1(i-ne) + gei1(i-ne)
              v(gei2(j)) = v(gei2(j)) + j0ei * mean_x(i) / k0
           ENDDO
           mean_x(i) = mean_x(i) - U_stp * mean_x(i)
        ENDIF
     ENDDO

  ENDDO

  !=============================================================================
  ! FINAL STATISTICAL CALCULATIONS (Synchronization & CV)
  !=============================================================================
  IF (count_samples > 0.0d0) THEN
     v_mean_temp = (v_mean_temp / count_samples)**2
     v2_mean_temp = v2_mean_temp / count_samples
     var_mean_v = v2_mean_temp - v_mean_temp
     
     DO i = 1, n
        v_i_mean(i) = (v_i_mean(i) / count_samples)**2
        v_i2_mean(i) = v_i2_mean(i) / count_samples
        mean_var_v(i) = v_i2_mean(i) - v_i_mean(i)
     ENDDO

     sync_rho = SQRT(var_mean_v / (SUM(mean_var_v(:)) / DBLE(n)))
  ENDIF
  
  DO i = 1, n
     IF(num_isi(i) > 0.0d0) THEN
        local_std = SQRT(ABS(isi_sq_sum(i)/num_isi(i) - (isi_sum(i)/num_isi(i))**2))
        local_cv = local_std / (isi_sum(i)/num_isi(i))
        
        count_valid_cv = count_valid_cv + 1.0d0
        cv_mean = cv_mean + local_cv
        
        IF(i <= ne) THEN
           count_cve = count_cve + 1.0d0
           cv_e = cv_e + local_cv
           WRITE(400, "(i8,f16.8)") in_degree(i), local_cv
        ELSE
           count_cvi = count_cvi + 1.0d0
           cv_i = cv_i + local_cv
           WRITE(800, "(i8,f16.8)") in_degree(i), local_cv
        ENDIF
     ENDIF
  ENDDO
  
  IF (count_valid_cv > 0) cv_mean = cv_mean / count_valid_cv
  IF (count_cve > 0) cv_e = cv_e / count_cve
  IF (count_cvi > 0) cv_i = cv_i / count_cvi

  WRITE(700, "(6f32.16)") cv_mean, cv_e, cv_i, sync_rho, var_mean_v, SUM(mean_var_v(:))
  
  DEALLOCATE (ge2, gi2, gei2, gie2)

  CALL CPU_TIME(stop_time)
  WRITE(*,*) 'Simulation completed in: ', stop_time - start_time, ' seconds.'

END PROGRAM main_simulation_asymmetric