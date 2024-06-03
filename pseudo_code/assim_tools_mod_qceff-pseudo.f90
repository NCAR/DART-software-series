! Pseudo code for filter_assim

!  Arguements:
! ens_handle  - this is the ensemble (e.g. 80 model states)
!        Each processor has every ensemble member, but only for part of the state.
!        Anything with my_  indicates the section of the state processor owns.
!        The assimilation updates the ensemble.
!
! obs_ens_handle - this contains the ensemble of forward operators (what each ensemble
!        member thinks the observation should be)
!        Each processor has every ensemble member, but only for some of the observations.
!        The assimilation updates the ensemble of forward operators.
!
! obs_seq - this is the actual obervations.  e.g 80C in Denver, CO, and can
!  contain other data.
!
! keys - this gives the index into the obs_seq for the observations.
!
! ens_size - the number of ensemble members
!
! obs_val_index - tells you which index in the obs_sequece is the observation value.
!
! The CAPITAL_WORDS are parameters:
!
!   ENS_MEAN_COPY, ENS_SD_COPY,
!   ENS_INF_COPY, ENS_INF_SD_COPY, OBS_KEY_COPY, OBS_GLOBAL_QC_COPY,
!   OBS_PRIOR_MEAN_START, OBS_PRIOR_MEAN_END, OBS_PRIOR_VAR_START,
!   OBS_PRIOR_VAR_END
!
!  The parameters are really integer values. They make the code more human readable.
!  For example, for an 80 member ensemble:
!    ENS_MEAN_COPY=81
!    ENS_SD_COPY=82
!
!  MISSING_R8 is a value that indicates DART was unable to calulate something: -88888.
!
! Ignoring:
! when you look at the full code for filter assim you will see the following:
!  * groups - groups chops up the ensemble into groups
!  * vertical conversion - the observation coordinate and the localization unit
!    can be different. e.g. some observations Z is in pressure, but the localiztion
!    might be in meters.
!  * inflate_only - this is on option for not assimilating
!  * some of the more compilacated inflation options
!  * caching of get_close
! I've striped out these for the pseudo-code.

subroutine filter_assim(ens_handle, obs_ens_handle, obs_seq, keys,           &
   ens_size,  obs_val_index, inflate, ENS_MEAN_COPY, ENS_SD_COPY, &
   ENS_INF_COPY, ENS_INF_SD_COPY, OBS_KEY_COPY, OBS_GLOBAL_QC_COPY,          &
   OBS_PRIOR_MEAN_START, OBS_PRIOR_MEAN_END, OBS_PRIOR_VAR_START,            &
   OBS_PRIOR_VAR_END)


!-- Section 1 --------------------------------------------------
! This section is where each processor gets the infomation on
! the observations and the parts of the state it owns.
! DART treats the state as a vector, so even though the model maybe 3D: X,Y,Z, variable
! in DART this is one long vector(:). This section is getting the X,Y,Z, and variable
! (meta data) for each state element.

! Get info on my number and indices for obs
my_num_obs = get_my_num_vars(obs_ens_handle)
call get_my_vars(obs_ens_handle, my_obs_indx)

! Get the locations for all of my observations
call get_my_obs_loc(obs_ens_handle, obs_seq,   keys, my_obs_loc, my_obs_kind, my_obs_type, obs_time)

! Get info on my number and indices for state
my_num_state = get_my_num_vars(ens_handle)
call get_my_vars(ens_handle, my_state_indx)

! Get the quantity and location all my state variables
! e.g. pressure at 89deg log, 51deg lat.
do i = 1, ens_handle%my_num_vars
   call get_state_meta_data(my_state_indx(i), my_state_loc(i), my_state_kind(i))

   ! Convert all my state variables to appropriate probit space.  RESPECT THE BOUNDS!
   call transform_to_probit(ens_size, ens_handle%copies(1:ens_size, i), dist_for_state, &
      state_dist_params(i), probit_ens, .false., &
      bounded_below, bounded_above, lower_bound, upper_bound)
   ens_handle%copies(1:ens_size, i) = probit_ens

end do

! 
do i = 1, my_num_obs
   obs_qc = obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, i)
   ! Only do conversion of qc if forward operator is good
   if(nint(obs_qc) == 0) then
      ! Need to specify what kind of prior to use for each
      call probit_dist_info(my_obs_kind(i), .false., .false., dist_for_obs, &
         bounded_below, bounded_above, lower_bound, upper_bound)
   
      ! Convert all my fwd operators (extended state) variables to appropriate probit space. RESPECT THE BOUNDS!
      call transform_to_probit(ens_size, obs_ens_handle%copies(1:ens_size, i), dist_for_obs, &
         obs_dist_params(i), probit_ens, .false., &
         bounded_below, bounded_above, lower_bound, upper_bound)
      obs_ens_handle%copies(1:ens_size, i) = probit_ens
   endif
end do


!-- Section 2 --------------------------------------------------
! This section is setup for localilzation calculations
! Any time you see get_close, this is for localization.

! Initialize the method for getting state variables close to a given ob on my process
call get_close_init(gc_state, my_num_state, 2.0_r8*cutoff, my_state_loc)

! Initialize the method for getting obs close to a given ob on my process
call get_close_init(gc_obs, my_num_obs, 2.0_r8*cutoff, my_obs_loc)


!-- Section 3 --------------------------------------------------
! This is the guts of filter_assim. It is where the assimilation
! happens. Note the loop is around all obervations

! Loop through all the (global) observations sequentially
SEQUENTIAL_OBS: do i = 1, obs_ens_handle%num_vars

!-- Section 3a --------------------------------------------------
! Every processor finds out what the observation is for this step of the loop,
!   kind: e.g. Temperature
!   value: what the actual value of the observation is
!   location: lat, lon, height
!   owner: which processor owns the observations (the owner has done the forward
!          operator caluclation for the observation)
!
   ! Every pe has information about the global obs sequence
   call get_obs_from_key(obs_seq, keys(i), observation)
   call get_obs_def(observation, obs_def)
   base_obs_loc = get_obs_def_location(obs_def)
   obs_err_var = get_obs_def_error_variance(obs_def)
   base_obs_type = get_obs_def_type_of_obs(obs_def)
   base_obs_kind = get_quantity_for_type_of_obs(base_obs_type)

   ! Get the value of the observation
   call get_obs_values(observation, obs, obs_val_index)

   ! Find out who has this observation and where it is
   call get_var_owner_index(ens_handle, int(i,i8), owner, owners_index)

!-- Section 3b --------------------------------------------------
!  There is a split here.  The owner of the obs calculates the obs_prior
!  and broadcasts this out to all other processors.

   ! Following block is done only by the owner of this observation
   !-----------------------------------------------------------------------
   if(ens_handle%my_pe == owner) then
      ! each task has its own subset of all obs.
 
      obs_qc = obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, owners_index)
      ! Only value of 0 for DART QC field should be assimilated
      IF_QC_IS_OKAY: if(nint(obs_qc) ==0) then

         obs_prior = obs_ens_handle%copies(1:ens_size, owners_index)

         ! Compute the prior mean and variance for this observation
         orig_obs_prior_mean = obs_ens_handle%copies(OBS_PRIOR_MEAN_START:OBS_PRIOR_MEAN_END, owners_index)
         orig_obs_prior_var  = obs_ens_handle%copies(OBS_PRIOR_VAR_START:OBS_PRIOR_VAR_END, owners_index)

        ! convert this observation ensemble from probit to regular space
         call transform_from_probit(ens_size, obs_ens_handle%copies(1:ens_size, owners_index) , &
            obs_dist_params(owners_index), obs_ens_handle%copies(1:ens_size, owners_index))

      endif IF_QC_IS_OKAY

      !Broadcast the info from this obs to all other processes
      call broadcast_send(map_pe_to_task(ens_handle, owner), obs_prior, obs_inc, &
           net_a, scalar1=obs_qc, &
           scalar2=vertvalue_obs_in_localization_coord, scalar3=whichvert_real)

   ! Next block is done by processes that do NOT own this observation
   !-----------------------------------------------------------------------
   else
      ! I don't store this obs; receive the obs prior from broadcast
      ! Also get qc and inflation information if needed
      ! also a converted vertical coordinate if needed
      call broadcast_recv(map_pe_to_task(ens_handle, owner), obs_prior, obs_inc, &
           net_a, scalar1=obs_qc, &
           scalar2=vertvalue_obs_in_localization_coord, scalar3=whichvert_real)

   endif
   !-----------------------------------------------------------------------

!-- Section 3c --------------------------------------------------
!  If the observation is bad, for example the forward operator could
!  not be calculated we skip to the next observation

   ! Everybody is doing this section, cycle if qc is bad
   if(nint(obs_qc) /= 0) cycle SEQUENTIAL_OBS
  
   ! Need to specify what kind of prior to use for obs being assimilated
      call probit_dist_info(base_obs_kind, .false., .false., dist_for_obs, &
         bounded_below, bounded_above, lower_bound, upper_bound)

   ! Convert the prior and posterior for this observation to probit space
   call transform_to_probit(grp_size, obs_prior(grp_bot:grp_top), dist_for_obs, &
      temp_dist_params, probit_obs_prior(grp_bot:grp_top), .false., &
      bounded_below, bounded_above, lower_bound, upper_bound)
   call transform_to_probit(grp_size, obs_post(grp_bot:grp_top), dist_for_obs, &
      temp_dist_params, probit_obs_post(grp_bot:grp_top), .true., &
      bounded_below, bounded_above, lower_bound, upper_bound)

   ! Compute observation space increments
   call obs_increment(obs_prior(grp_bot:grp_top), grp_size, obs(1), &
         obs_err_var, obs_inc(grp_bot:grp_top), inflate, my_inflate,   &
         my_inflate_sd, net_a(group))

   ! compute prior mean and variance of obs
   obs_prior_mean = sum(obs_prior) / ens_size
   obs_prior_var = sum((obs_prior - obs_prior_mean)**2) / ens_size
  
!-- Section 3d --------------------------------------------------
!  For localization, we need to find which oberservations and which bits
!  of the state are within the localization radius


   ! Find observations on my proceess that are close to observation being assimilated
   call get_close_obs(gc_obs, base_obs_loc, base_obs_type, &
                      my_obs_loc, my_obs_kind, my_obs_type, &
                      num_close_obs, close_obs_ind, close_obs_dist, ens_handle)


   
   ! Find state variables on my process that are close to observation being assimilated
   call get_close_state(gc_state, base_obs_loc, base_obs_type, &
                        my_state_loc, my_state_kind, my_state_indx, &
                        num_close_states, close_state_ind, close_state_dist, ens_handle)

!-- Section 3e --------------------------------------------------
! Update the state

   ! Now everybody updates their close states
   ! Loop through to update each of my state variables that is potentially close
   STATE_UPDATE: do j = 1, num_close_states
      state_index = close_state_ind(j)

      ! Compute the covariance localization and adjust_obs_impact factors
      final_factor = cov_and_impact_factors(base_obs_loc, base_obs_type, my_state_loc(state_index), &
      my_state_kind(state_index), close_state_dist(j), cutoff_rev)

      if(final_factor <= 0.0_r8) cycle STATE_UPDATE

      ! Update the state variable ensemble members
      call obs_updates_ens(ens_size, num_groups, ens_handle%copies(1:ens_size, state_index), &
            my_state_loc(state_index), my_state_kind(state_index), obs_prior, obs_inc, &
            obs_prior_mean, obs_prior_var, base_obs_loc, base_obs_type, obs_time, &
            net_a, grp_size, grp_beg, grp_end, i, &
            my_state_indx(state_index), final_factor, correl, local_varying_ss_inflate, inflate_only)

   end do STATE_UPDATE

   !------------------------------------------------------

!-- Section 3f --------------------------------------------------
! Update the forward operators
! Note, only updating for observations that have not already been
! assimilated.

   ! Now everybody updates their obs priors (only ones after this one)
   OBS_UPDATE: do j = 1, num_close_obs
      obs_index = close_obs_ind(j)

      ! Only have to update obs that have not yet been used
      if(my_obs_indx(obs_index) > i) then

         ! If the forward observation operator failed, no need to
         ! update the unassimilated observations
         if (any(obs_ens_handle%copies(1:ens_size, obs_index) == MISSING_R8)) cycle OBS_UPDATE

         ! Compute the covariance localization and adjust_obs_impact factors (module storage)
         final_factor = cov_and_impact_factors(base_obs_loc, base_obs_type, my_obs_loc(obs_index), &
         my_obs_kind(obs_index), close_obs_dist(j), cutoff_rev)

         if(final_factor <= 0.0_r8) cycle OBS_UPDATE

         ! Do the update
         call obs_updates_ens(ens_size, num_groups, obs_ens_handle%copies(1:ens_size, obs_index), &
            my_obs_loc(obs_index), my_obs_kind(obs_index), obs_prior, obs_inc, &
            obs_prior_mean, obs_prior_var, base_obs_loc, base_obs_type, obs_time, &
            net_a, grp_size, grp_beg, grp_end, i, &
            -1*my_obs_indx(obs_index), final_factor, correl, .false., inflate_only)
 
      endif
   end do OBS_UPDATE


end do SEQUENTIAL_OBS


! Do the inverse probit transform for state variables
call transform_all_from_probit(ens_size, ens_handle%my_num_vars, ens_handle%copies, &
   state_dist_params, ens_handle%copies)


!-------------------------------------------------------------

! Assure user we have done something
if (print_trace_details >= 0) then
write(msgstring, '(A,I8,A)') &
   'Processed', obs_ens_handle%num_vars, ' total observations'
   call error_handler(E_MSG,'filter_assim:',msgstring)
endif

end subroutine filter_assim
