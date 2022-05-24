
module filter_mod

contains

!----------------------------------------------------------------
subroutine filter_main()


! A bunch of initialization.

! Set up stages to write : input, preassim, postassim, output
call parse_stages_to_write(stages_to_write)

! Count and set up State copy numbers
num_state_ens_copies = count_state_ens_copies(ens_size, prior_inflate, post_inflate)
num_extras           = num_state_ens_copies - ens_size

! Initialize the obs_sequence; every pe gets a copy for now
call filter_setup_obs_sequence(seq, in_obs_copy, obs_val_index, input_qc_index, DART_qc_index, compute_posterior)

! Allocate model size storage and ens_size storage for metadata for outputting ensembles
model_size = get_model_size()

if(distributed_state) then
   call init_ensemble_manager(state_ens_handle, num_state_ens_copies, model_size)
   msgstring = 'running with distributed state; model states stay distributed across all tasks for the entire run'
else
   call init_ensemble_manager(state_ens_handle, num_state_ens_copies, model_size, transpose_type_in = 2)
   msgstring = 'running without distributed state; model states are gathered by ensemble for forward operators'
endif

call set_num_extra_copies(state_ens_handle, num_extras)

call initialize_file_information(num_state_ens_copies ,                     &
                                 file_info_input      , file_info_mean_sd,  &
                                 file_info_forecast   , file_info_preassim, &
                                 file_info_postassim  , file_info_analysis, &
                                 file_info_output)

! Read the state (ensemble of model forcasts)
call read_state(state_ens_handle, file_info_input, read_time_from_file, time1, &
                prior_inflate, post_inflate, perturb_from_single_instance)

! Compute mean and spread
call compute_copy_mean_sd(state_ens_handle, 1, ens_size, ENS_MEAN_COPY, ENS_SD_COPY)

time_step_number = -1


! Inifinte loop. 
! Run until you run out of observations
!  or
! Exit if you can not advance the model (run a forcast) inside filter
AdvanceTime : do

   time_step_number = time_step_number + 1

   !------------------------
   ! Check the time on the state.  
   call move_ahead(state_ens_handle, ens_size, seq, last_key_used, window_time, &
      key_bounds, num_obs_in_set, curr_ens_time, next_ens_time)
   !------------------------

   
   !------------------------
   ! Exit if there are no more observations in the sequence
   if(key_bounds(1) < 0) then
      call trace_message('No more obs to assimilate, exiting main loop', 'filter:', -1)
      exit AdvanceTime
   endif
   !------------------------

   ! if model state data not at required time, advance model
   ! Big models - filter exits if the model is not at the correct time
   ! Small models - 
   if (curr_ens_time /= next_ens_time) then
      ! Must be done before the model runs and updates the data.

      call allocate_vars(state_ens_handle)

      call all_copies_to_all_vars(state_ens_handle)

      call advance_state(state_ens_handle, ens_size, next_ens_time, async, &
              adv_ens_command, tasks_per_model_advance, file_info_output, file_info_input)

      call all_vars_to_all_copies(state_ens_handle)

      ! updated mean and spread after the model advance
      call compute_copy_mean_sd(state_ens_handle, 1, ens_size, ENS_MEAN_COPY, ENS_SD_COPY)

      ! update so curr time is accurate.
      curr_ens_time = next_ens_time
      state_ens_handle%current_time = curr_ens_time
      call set_time_on_extra_copies(state_ens_handle)

      ! only need to sync here since we want to wait for the
      ! slowest task to finish before outputting the time.
      call timestamp_message('After  running model', sync=.true.)
   else
      call trace_message('Model does not need to run; data already at required time', 'filter:', -1)
   endif

   call trace_message('Before setup for next group of observations')


   !------------------------
   ! Make space for the results of the forward operator
   ! Create an ensemble for the observations from this time plus
   ! obs_error_variance, observed value, key from sequence, global qc,
   ! then mean for each group, then variance for each group
   call init_ensemble_manager(obs_fwd_op_ens_handle, TOTAL_OBS_COPIES, &
                              int(num_obs_in_set,i8), 1, transpose_type_in = 2)

   ! Also need a qc field for copy of each observation
   call init_ensemble_manager(qc_ens_handle, ens_size, &
                              int(num_obs_in_set,i8), 1, transpose_type_in = 2)
   !------------------------

   !------------------------
   ! Inflation 
   !------------------------

   !------------------------
   ! Forward operator
   call get_obs_ens_distrib_state(state_ens_handle, obs_fwd_op_ens_handle, &
           qc_ens_handle, seq, keys, obs_val_index, input_qc_index, &
           OBS_ERR_VAR_COPY, OBS_VAL_COPY, OBS_KEY_COPY, OBS_GLOBAL_QC_COPY, &
           OBS_EXTRA_QC_COPY, OBS_MEAN_START, OBS_VAR_START, &
           isprior=.true., prior_qc_copy=prior_qc_copy)
   !------------------------
 
   !------------------------
   ! Obs space diagnostics
   !------------------------

   !------------------------
   ! Assimilation
   call filter_assim(state_ens_handle, obs_fwd_op_ens_handle, seq, keys, &
      ens_size, num_groups, obs_val_index, prior_inflate, &
      ENS_MEAN_COPY, ENS_SD_COPY, &
      PRIOR_INF_COPY, PRIOR_INF_SD_COPY, OBS_KEY_COPY, OBS_GLOBAL_QC_COPY, &
      OBS_MEAN_START, OBS_MEAN_END, OBS_VAR_START, &
      OBS_VAR_END, inflate_only = .false.)
   !------------------------ 

   
   !------------------------
   ! Posterior inflation
   !-------------------------


   !------------------------
   ! Posterior forward operator
   if (compute_posterior) then
   
       call get_obs_ens_distrib_state(state_ens_handle, obs_fwd_op_ens_handle, &
                qc_ens_handle, seq, keys, obs_val_index, input_qc_index, &
                OBS_ERR_VAR_COPY, OBS_VAL_COPY, OBS_KEY_COPY, OBS_GLOBAL_QC_COPY, &
                OBS_EXTRA_QC_COPY, OBS_MEAN_START, OBS_VAR_START, &
                isprior=.false., prior_qc_copy=prior_qc_copy)
   
      !------------------------
      ! Write posterior observation space diagnostics
      !------------------------
   endif
   !------------------------ 


   !------------------------
   ! adaptive state space posterior inflation
   ! this block computes the adaptive state space posterior inflation
   ! (it was applied earlier, this is computing the updated values for
   ! the next cycle.)
   !------------------------

   !------------------------
   ! Free up the obs ensemble space
   call end_ensemble_manager(obs_fwd_op_ens_handle)
   call end_ensemble_manager(qc_ens_handle)
   !------------------------

   call trace_message('Bottom of main advance time loop')

end do AdvanceTime


!------------------------
! Output obs_sequence and state restart files
! Only pe 0 outputs the observation space diagnostic file
if(my_task_id() == 0) call write_obs_seq(seq, obs_sequence_out_name)

! Output all restart files if requested
file_info_all = combine_file_info( &
                 (/file_info_input, file_info_mean_sd, file_info_forecast, &
                   file_info_preassim, file_info_postassim, file_info_analysis, &
                   file_info_output/) )

call write_state(state_ens_handle, file_info_all)
!------------------------


!------------------------
! Clean up
call end_assim_model()

call end_ensemble_manager(state_ens_handle)

! Free up the obs sequence
call destroy_obs_sequence(seq)
!------------------------

call     trace_message('Filter done')
call timestamp_message('Filter done')

end subroutine filter_main
