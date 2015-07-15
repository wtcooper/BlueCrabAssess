//###################################################################
//###################################################################

// Adapted from Chesapeake Bay 2010 blue crab assessment 
//By: M. Wilberg 2/28/2011
//
//Adapted by W Cooper 2012 
//Major changes:
//	1) Pulled out sex-specific due to limited biostatistical sampling of landings in Gulf 
//	2) Added stage-specific mortality, including in ref points
//	3) Added environmental influences on S-R and stage-specific M 
// 4) Added retrospective analyses and projections

//###################################################################
//###################################################################

TOP_OF_MAIN_SECTION
//increase number of estimated parameters
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(2000);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(2000400);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(10000000);
  arrmblsize = 10000000;


//###################################################################
//###################################################################
DATA_SECTION
  init_adstring DataFile;
  !! ad_comm::change_datafile_name(DataFile);

  init_int testing	//toggle to turn on/off console output for testing
  init_int fyear   //first year of the model run
  init_int lyear   //last year of the model run
  init_int retroYears
  init_int projYears

  int timeSteps  //number of time steps in a year
  int stages	//number of stages/ages to model
  !! timeSteps=1;
  !! stages=2;

  int retroSteps
  int projSteps
  !! retroSteps=retroYears*timeSteps;
  !! projSteps=projYears*timeSteps;

  //Because coded to be multiple time steps per year, need to define an indexing scheme that isn't year-based
  int mTimeSteps	// total model time steps for population dynamics
  int startIndex	//starting index for the model
  int endIndex		//ending index for the model
  !!mTimeSteps = (lyear-fyear+1)*timeSteps; //calculate the time steps in the model
  !!startIndex=1000;		//set a value, can be anything > the lag of environment time series 

  //substract off the retrospective period
  !!endIndex=startIndex+mTimeSteps-1-retroSteps;  //substract 1 since startIndex is first step



  //#################################
  //#Catch Data
  //#################################
  init_int ftcyear  //first year of total catch
  init_int ltcyear  //last year of total catch
  int cStartIndex
  int cTimeSteps  //catch time steps based on catch years
  !!cTimeSteps = (ltcyear-ftcyear+1)*timeSteps;
  !!cStartIndex = (ftcyear-fyear)*timeSteps+startIndex;
  init_vector com_TC_obs(cStartIndex,cStartIndex+cTimeSteps-1)  //total catch
  init_number C_cv					//catch data standard deviation
  init_int effFlag		//flag to use effort time series
  init_vector com_Eff_obs(cStartIndex,cStartIndex+cTimeSteps-1)  //total catch

  //#################################
  //#Survey data
  //#################################
  //Adults
  init_int numAdSurv
  init_int fsayear  //first year of adult surveys
  init_int lsayear  //last year of adult surveys
  int adTimeSteps
  int adStartIndex
  !!adTimeSteps=(lsayear-fsayear+1)*timeSteps;
  !!adStartIndex = (fsayear-fyear)*timeSteps+startIndex;
  init_matrix ad_survey_obs(1,numAdSurv,adStartIndex,adStartIndex+adTimeSteps-1)  //Adult survey CPUE
  init_matrix ad_survey_cv(1,numAdSurv,adStartIndex,adStartIndex+adTimeSteps-1)  	//Survey CVs for adults
  init_vector sa_time(1,numAdSurv)		//survey time

  //Recruits
  init_int numRecSurv
  init_int fsryear
  init_int lsryear
  int recTimeSteps
  int recStartIndex
  !!recTimeSteps=(lsryear-fsryear+1)*timeSteps;
  !!recStartIndex = (fsryear-fyear)*timeSteps+startIndex;
  init_matrix re_survey_obs(1,numRecSurv,recStartIndex,recStartIndex+recTimeSteps-1)  //Recruit survey CPUE
  init_matrix re_survey_cv(1,numRecSurv,recStartIndex,recStartIndex+recTimeSteps-1)  	//survey SDs for recruits
  init_vector sr_time(1,numRecSurv)

  //#################################
  //#Fishery params
  //#################################
  init_number p_rec  //Proportion of recreational harvest per region
  init_number p_under  //Proportion of harvest underreporting per region
  init_number maxF	//Max F for F_pen calculation
  init_number maxM	//Max M for F_pen calculation


  //#################################
  //#Adult Z estimates as prior
  //#################################
  init_number aveZ
  init_number Z_cv


  //#################################
  //#Life History params
  //#################################
  init_number sratio					//Sex ratio
  init_vector M(1,stages)			//mortality at age for each stage (e.g., for CS: recruits, post-recruits)
  vector Myr(1,stages)			//M rate on per year basis (for ref points)
  init_vector pSpawn(1,stages)  	//proportion of females spawning in each season for differential spawning throughout the year
  init_number sp_time  //proportion of the time step before spawning occurs
  init_int SRSwitch				//switch for recruit function formulation; 1=bev holt, 2=ricker

 LOCAL_CALCS
//convert M to the appropriate time frame from per year basis
  Myr=M;
  M=M/timeSteps;
 END_CALCS


  //#################################
  //#Environmental time series params/data
  //#################################
  init_int numEnvTS
  init_int feyear   //first year of the model
  init_int leyear   //last year of the model
  int eTimeSteps
  int eStartIndex		//this should be less than the startIndex
  !!eTimeSteps = (leyear-feyear+1)*timeSteps;
  !!eStartIndex = (feyear-fyear)*timeSteps+startIndex;
  init_matrix envObs(1,numEnvTS,eStartIndex,eStartIndex+eTimeSteps-1)	//environmental time series (regions, season,timesteps)
  init_vector env_cv(1,numEnvTS)	//environmental time series (regions, season,timesteps)
  init_int envRecTS					//time series # that influences recruitment
  init_int envRecLag				//lag of recruitment influence
  init_vector envMTS(1,stages)					//time series # that influences mortality
  init_vector envMLag(1,stages)						//lag in mortality influence
  matrix env(1,numEnvTS,eStartIndex,endIndex+projSteps+1) //add on one for the forward stepping recruitment


  //#################################
  //#Projections
  //#################################
  init_matrix envProj(1,numEnvTS,endIndex+1,endIndex+projSteps+1)
  init_vector effProj(endIndex+1,endIndex+projSteps+1) //START with terminal year since non-estimable F


  //#################################
  //#Parameters (initial val, min, max, phase)
  //#################################
  
  init_number init_NIni
  init_vector init_NParams(1,3)
  init_number init_RIni
  init_vector init_RParams(1,3)

  //F params
  init_number F_qIni
  init_vector F_qParams(1,3)
  init_number F_devIni
  init_vector F_devParams(1,3)
  init_number eff_cvIni
  init_vector eff_cvParams(1,3)

  //recruitment  params
  init_number rec_devIni
  init_vector rec_devParams(1,3)
  init_number rec_cvIni
  init_vector rec_cvParams(1,3)

  init_number S0Ini
  init_vector S0Params(1,3)
  init_number steepIni
  init_vector steepParams(1,3)

  init_number sr_beta_envIni
  init_vector sr_beta_envParams(1,3)
  init_vector M_beta_envIni(1,stages)
  init_vector M_beta_envParams(1,3)
  init_number M_cvIni
  init_vector M_cvParams(1,3)
  init_vector selIni(1,stages)
  init_vector selParams(1,3)


  //#################################
  //#Likelihood Weights
  //#################################
  
  init_number com_lambda		//survey weight
  init_vector sa_lambda(1,numAdSurv)		//survey weight
  init_vector sr_lambda(1,numRecSurv)
  init_number recDev_lambda		//survey weight
  init_number effort_lambda		//survey weight
  init_number aveZ_lambda		//survey weight


  //#################################
  //#Additional param control flags not addressed in data section
  //#################################
  init_number biasAdj					//Adjustment multiplier for bias correction factor


  //#################################
  //#Reference point calcs
  //#################################
  //number for reference point explorations
  init_number Fval_init  							//lowest value of F used in SPR calcs
  init_number Fval_max   						//highest value of F used in SPR calcs
  init_number Fval_inc   							//increment for F
  int Fval_num
  !!Fval_num=(Fval_max-Fval_init)/Fval_inc+1;
  init_int nspr          //number of values F-%SPR will be calculated at
  init_vector SPR_targ(1,nspr)  //Values of SPR for Fval reference point calculations


  //#################################
  //#EOF Test
  //#################################
  init_int test			//check that data read in appropriately



  //#################################
  //#Additional Variables
  //#################################
  
  //Total harvest including recreational
  vector TC_obs(cStartIndex,cStartIndex+cTimeSteps-1)  //total catch

  //Variances for data sets
  number C_var  //variance of catch
  number Z_var  //variance of catch
  matrix ad_survey_var(1,numAdSurv,adStartIndex,adStartIndex+adTimeSteps-1)  //variances for adult surveys
  matrix re_survey_var(1,numRecSurv,recStartIndex,recStartIndex+recTimeSteps-1)  //variances for recruitment surveys


  //Define index varaibles
  int y  //index variable for time step
  int s  //index variable for season
  int r  //index variable for region
  int i  //index variable 
  number year  //for report section
  int ispr
  int iter
  int iterMCMC
  !!iterMCMC=0;
  int index

 LOCAL_CALCS

  if (SRSwitch==2) steepParams(2)=5.0;	//make sure to bound steepness appropriately for Ricker

  C_var=log(C_cv*C_cv+1);  //variance of catch
  Z_var=log(Z_cv*Z_cv+1);  //variance of Z

  //commercial + recreational catch+prop. underreported 
  TC_obs=com_TC_obs*(1.+p_rec+p_under);

  //Calculate variances from CVs
  for (i=1; i<=numAdSurv; i++){
    for (y=adStartIndex; y<=adStartIndex+adTimeSteps-1; y++){
      ad_survey_var(i,y)=log(ad_survey_cv(i,y)*ad_survey_cv(i,y)+1);  //variances for adult surveys
    }
  }
  for (i=1; i<=numRecSurv; i++) {
    for (y=recStartIndex; y<=recStartIndex+recTimeSteps-1; y++){
      re_survey_var(i,y)=log(re_survey_cv(i,y)*re_survey_cv(i,y)+1);  //variances for recruitment surveys
    }
  }

  for (i=1; i<=numEnvTS; i++){
    for (y=eStartIndex; y<=endIndex; y++){
      env(i,y)=envObs(i,y);
    }
    for (y=endIndex+1; y<=endIndex+projSteps; y++){
      env(i,y)=envProj(i,y);
    }
  }

  if(test!=12345)  //check to make sure end of file number is correct
  {
    //if not correct, output the data and exit.
    cout << "Data not reading properly" << endl;
    cout << "Commercial\t" << com_TC_obs << endl;
    cout << "adults\t" <<ad_survey_obs << endl;
    cout << "recruits\t" <<re_survey_obs << endl;
    cout << "environment\t" <<env << endl;
    cout << "max F\t" <<maxF << endl;
    cout << "max M\t" <<maxM << endl;
    cout << "envMTS\t" <<envMTS << endl;
    cout << "env Lag\t" <<envMLag << endl;
    cout << "ini N\t" <<init_NIni << endl;
    cout << "steepIni\t" <<steepIni << endl;
    cout << "sel params\t" <<selParams<< endl;
    cout << "EOF test: " << test << endl;
    exit(1);
  }

 END_CALCS



//###################################################################
//###################################################################
PARAMETER_SECTION

//Copy the Parameter estimates to double vals so can put as bounds,phases below
 LOCAL_CALCS
//template: double xxxMin=log(xxxParams(1)); double xxxMax=log(xxxParams(2)); double xxxPhase=log(xxxParams(2));
  double log_init_NMin=log(init_NParams(1)); 			double log_init_NMax=log(init_NParams(2)); 			double init_NPhase=init_NParams(3);
  double log_init_RMin=log(init_RParams(1)); 			double log_init_RMax=log(init_RParams(2)); 			double init_RPhase=init_RParams(3);
  double log_F_qMin=log(F_qParams(1)); 			double log_F_qMax=log(F_qParams(2)); 		double F_qPhase=F_qParams(3);
  double log_F_devMin=log(F_devParams(1)); 			double log_F_devMax=log(F_devParams(2)); 			double F_devPhase=F_devParams(3);
  double log_eff_cvMin=log(eff_cvParams(1)); 		double log_eff_cvMax=log(eff_cvParams(2)); 			double eff_cvPhase=eff_cvParams(3);
  double log_rec_devMin=log(rec_devParams(1)); 	double log_rec_devMax=log(rec_devParams(2));		double rec_devPhase=rec_devParams(3);
  double log_rec_cvMin=log(rec_cvParams(1)); 		double log_rec_cvMax=log(rec_cvParams(2)); 			double rec_cvPhase=rec_cvParams(3);
  double log_S0Min=log(S0Params(1)); 					double log_S0Max=log(S0Params(2)); 						double S0Phase=S0Params(3);
  double log_steepMin=log(steepParams(1)); 			double log_steepMax=log(steepParams(2)); 				double steepPhase=steepParams(3);
  double sr_beta_envMin=sr_beta_envParams(1); 		double sr_beta_envMax=sr_beta_envParams(2); 		double sr_beta_envPhase=sr_beta_envParams(3);
  double M_beta_envMin=M_beta_envParams(1); 		double M_beta_envMax=M_beta_envParams(2); 		double M_beta_envPhase=M_beta_envParams(3);
  double log_M_cvMin=log(M_cvParams(1)); 		double log_M_cvMax=log(M_cvParams(2)); 			double M_cvPhase=M_cvParams(3);
  double log_selMin=log(selParams(1)); 					double log_selMax=log(selParams(2)); 					double selPhase=selParams(3);

 END_CALCS


  //initial R and N
  //template: (Min,Max,Phase)
  init_bounded_number log_init_N(log_init_NMin,log_init_NMax,init_NPhase)
  init_bounded_number log_init_R(log_init_RMin,log_init_RMax,init_RPhase)

  //Fishing mortality for each year
  init_bounded_number log_F_q(log_F_qMin,log_F_qMax,F_qPhase)
  init_bounded_vector log_F_dev(startIndex+1,endIndex-1,log_F_devMin,log_F_devMax,F_devPhase) //don't bother with terminal year deviation
  init_bounded_number log_eff_cv(log_eff_cvMin,log_eff_cvMax,eff_cvPhase)

  //Recruitment params
  init_bounded_dev_vector log_rec_dev(startIndex+1,endIndex,log_rec_devMin,log_rec_devMax,rec_devPhase)
  init_bounded_number log_rec_cv(log_rec_cvMin,log_rec_cvMax,rec_cvPhase)

  //Stock-recruitment parameters
  init_bounded_number log_S0(log_S0Min,log_S0Max,S0Phase)
  init_bounded_number log_steep(log_steepMin,log_steepMax,steepPhase)

  //S-R environmental link parameter
  init_bounded_number sr_beta_env(sr_beta_envMin,sr_beta_envMax,sr_beta_envPhase)

  //Adult environmental link parameter
  init_bounded_vector M_beta_env(1,stages,M_beta_envMin,M_beta_envMax,M_beta_envPhase)
  init_bounded_number log_M_cv(log_M_cvMin,log_M_cvMax,M_cvPhase)

  //Vulnerability at each stage
  init_bounded_vector log_sel(1,stages,log_selMin,log_selMax,selPhase)


  //############### Derived parameters ###############//
  
  //sdreport_matrix  for some of these
  sdreport_vector N(startIndex,endIndex+1+projSteps)   //abundance
  sdreport_vector R(startIndex,endIndex+1+projSteps)   //recruitment
  vector SP(startIndex,endIndex+projSteps) //number of spawners
  vector TC(startIndex,endIndex+projSteps)      //total catch
  vector u(startIndex,endIndex+projSteps)      //total exploitation rate
  sdreport_vector F(startIndex,endIndex+projSteps)   //fishing mortality rate
  vector effort(startIndex,endIndex+projSteps)   //fishing mortality rate
  vector sel(1,stages)				//selectivity (partial recruitment) of recruits to the fishery
  matrix Mt(1,stages,startIndex,endIndex+projSteps)		//M at time t; don't do for years to account for seasonal M
  vector Z(startIndex,endIndex+projSteps)		//total Z


  matrix ad_survey_est(1,numAdSurv,startIndex,endIndex)  //estimated adult survey indices
  matrix re_survey_est(1,numRecSurv,startIndex,endIndex)  //estimated recruitment survey indices
  vector qa(1,numAdSurv)   //catchability for adult surveys 
  vector qr(1,numRecSurv)   //catchability for recruitment surveys


  number S0
  number steep
  number R0
  number A0
  number alpha			//Alpha of the S-R relationship
  number beta			//Beta of the S-R relationship
  number rec_var		//variance for recruitment deviations
  number eff_var			//variance for effort residuals 
  number M_var			//variance for F deviations


  //Derived variables for reference point calculations 
  number Fval							//F for SPR calculations
  vector SPR(1,Fval_num)				//spawners per recruit (NOT spawning potential ratio)
  number SPR0							//virgin SPR
  number SPR1							//spawner per recruit per year
  vector tSPR(startIndex,endIndex+projSteps) //transitional SPR
  vector NPR(1,Fval_num)			//numbers per reruit
  vector YPR(1,Fval_num)				//yield per recruit
  vector SPRatio(1,Fval_num)		//spawning potential ratio
  vector N_eq(1,Fval_num)			//equilibruim numbers
  vector R_eq(1,Fval_num)			//equilibrium recruitment
  vector C_eq(1,Fval_num)			//equilibrium catch
  vector C_eqSort(1,Fval_num)			//sorted equilibrium catch

  //Reference points  
  vector u0_eq(1,Fval_num)			//equilibrium explotation rate for age0+
  vector u1_eq(1,Fval_num)			//equilibrium explotation rate for age1+
  vector uAll_eq(1,Fval_num)			//equilibrium explotation rate for age0+
  vector FSPR_ref(1,nspr)				//F%SPR reference points
  vector SPRDiff(1,nspr)				//temporary array to check if at proper F for SPR calcs
  vector Fvec(1,Fval_num)			//equilibrium explotation rate for age0+

  likeprof_number MSY				//MSY estimate
  number OFL				//Overfishing Limit (Ncurrent*uMSY)
  number u0MSY			//exploitation rate at MSY for age 0
  number u1MSY			//exploitation rate at MSY for age 1+
  number uMSY			//exploitation rate at MSY for age 0+
  number FMSY			//F rate at MSY
  number RMSY		//equilibrium recruitment at msy
  number NMSY			//Number at MSY
  number FLim			//F Limit (target)
  number NLim			//N Limit (target)
  number cLim			//c used in calculation of targets
  number FCurr			//Current F (geometric mean of last 3 years of model run, not including terminal year)
  vector MCurr(1,stages)			//Current M (geometric mean of last 3 years of model run)
  number NCurr			//Current N (geometric mean of last 3 years of model run)
  number SPRCurr			//Current SPR using FCurr and MCurr
  number FMSYRatio	//Fcurr/FMSY
  number NMSYRatio	//Ncurr/NMSY
  number UMSYRatio	//Ncurr/NMSY
  number FFLimRatio	//Fcurr/FLimit
  number NNLimRatio	//Ncurr/NLim
  number termToMSY  //for F calcs 


  //variables for likelihood function
  vector Lsr(1,numRecSurv)  //likelihood components for recruit surveys
  vector Lsa(1,numAdSurv)  //likelihood components for adult surveys
  number Lc        //likelihood components for catch time series 
  number Lz        //likelihood components for adult Z estimate (anchors F/N0/R0)  
  number Lrdev           //likelihood for recruitment deviations
  number Leff           //likelihood for effort residuals

  number F_pen           //penalty for F above the max F
  number M_pen           //penalty for M above the max M
  objective_function_value negLL 	//negative log-liklihood


  //############### Starting parameter values ###############//
 LOCAL_CALCS

  log_init_N=log(init_NIni);
  log_init_R=log(init_RIni);
  log_F_q=log(F_qIni);
  log_F_dev=log(F_devIni);
  log_M_cv=log(M_cvIni);
  log_rec_cv=log(rec_cvIni);
  log_rec_dev=log(rec_devIni);
  log_S0=log(S0Ini);
  log_steep=log(steepIni);
  sr_beta_env=sr_beta_envIni;
  M_beta_env=M_beta_envIni;
  log_sel=log(selIni);

 END_CALCS



//###################################################################
//###################################################################
PROCEDURE_SECTION
  set_initial_conditions();
  if (testing==1)  cout << "End set_initial_conditions()" << endl;
  calculate_abundance_and_catch();
  if (testing==1)  cout << "End calculate_abundance_and_catch()" << endl;
  calculate_predicted_indices();
  if (testing==1)  cout << "End calculate_predicted_indices()" << endl;
  calculate_objective_function();
  if (testing==1)  cout << "End calculate_objective_function()" << endl;
  mcmc();
  if (testing==1) cout << "End mcmc()" << endl;

  if (testing==1) {
    calculate_tSPR();
    obs_pred();
    MSY_estimates();
    HPD_estimates();
    general_report();

    cout << "Procedure section completed first cycle, now exiting"<< endl;
    exit(1); //exit if in testing phase -- runs model at initial parameter values
  }


//############### Main Functions ###############

FUNCTION set_initial_conditions

  //convert parameters from the log scale
  S0=exp(log_S0);
  steep=exp(log_steep);
  negLL=0.0;
  M_var=log(exp(log_M_cv)*exp(log_M_cv)+1);
  rec_var=log(exp(log_rec_cv)*exp(log_rec_cv)+1);
  eff_var=log(exp(log_eff_cv)*exp(log_eff_cv)+1);
  sel=exp(log_sel);
  F_pen=0;
  M_pen=0;


  //######## S-R Params ########  
  
  //Calculate virgin SPR, including proportion of recruits spawning  
  A0=sratio*pSpawn(2)*exp(-(M(1)+sp_time*M(2)))/(1.-exp(-(M(2)))) + sratio*pSpawn(1)*exp(-(sp_time*M(1)));
  R0=S0/A0;

  if (SRSwitch==1)	{ //Beverton-Holt
    alpha = S0*(1-steep)/(4*steep*R0);
    beta = (5*steep-1)/(4*steep*R0);
  }

  if (SRSwitch==2)	{ //Ricker
    beta = log(5*steep)/(0.8*S0);
    alpha =(exp((5.*log(5.*steep))/4.)*R0)/S0;
  }

  //calculate reference points after setting S-R params so can get FMSY for F projections
  calculate_reference_points();


  //######## M ########  
  
  //compute the yearly M accounting for seasonal differences and environmental differences
  //leave this here to deal with seasonality
  for(y=startIndex; y<=endIndex+projSteps; y++) {

    //only apply deviation + bias correction if active
    Mt(1,y)=M(1);
    if (active(M_beta_env)) Mt(1,y)=M(1)*exp(M_beta_env(1)*env(envMTS(1),y-envMLag(1)))*exp(-0.5*M_var);

    if (y<=endIndex) {
      posfun(maxM-(timeSteps*Mt(1,y)),.000001,M_pen);
      negLL+=100.*M_pen;
    }

    Mt(2,y)=M(2);
    if (active(M_beta_env)) Mt(2,y)=M(2)*exp(M_beta_env(2)*env(envMTS(2),y-envMLag(2)))*exp(-0.5*M_var);


    if (y<=endIndex) {
      posfun(maxM-(timeSteps*Mt(2,y)),.000001,M_pen);
      negLL+=100.*M_pen;
    }
  }


  //######## F ########  
  
  if (effFlag==0) effort=1.0;
  else {
    //set up effort to average and replace all missing data with average
    double avg_effort=mean(com_Eff_obs);

    //for any year prior to effort data, set to avg of all other years
    for (i=startIndex; i<=endIndex; i++) effort(i)=com_Eff_obs(i);
    for (i=endIndex+1; i<=endIndex+projSteps; i++) effort(i)=avg_effort; //effort(endIndex)+effProj(i)*termToMSY; //deviation off the last year
    effort/=avg_effort; //scale to observed years and not including projected years
  }

  //If want to estimate ave. q instead of 1st year q: change F_dev to bounded_dev_vector and adjust here
  F(startIndex)=exp(log_F_q+log(effort(startIndex)));
  posfun(maxF-(timeSteps*F(startIndex)),.000001,F_pen);
  negLL+=100.*F_pen;

  for(y=startIndex+1; y<endIndex; y++) { //don't include terminal year estimate
    //Computed as F=q*Eff*exp(dev)
    F(y)=exp(log_F_q+log(effort(y))+log_F_dev(y));
    posfun(maxF-(timeSteps*F(y)),.000001,F_pen);
    negLL+=100.*F_pen;
  }

  //for terminal year, use estimated deviation from previous year to keep scaled together
  F(endIndex)=exp(log_F_q+log(effort(endIndex))+log_F_dev(endIndex-1));

  termToMSY=FMSY-F(endIndex); //effort range from terminal year to MSY; negative if FMSY<termF 

  //no F deviations on terminal year or projected years
  for(y=endIndex+1; y<=endIndex+projSteps; y++) {
    F(y)=(F(endIndex)+effProj(y)*termToMSY);
    effort(y)=exp(log(F(y)))/exp(log_F_q);
  }

  F=F/timeSteps;


  //######## Adult Z ########  
  
  for(y=startIndex; y<=endIndex+projSteps; y++) Z(y)=F(y)+Mt(2,y);




FUNCTION calculate_abundance_and_catch

  N(startIndex)=exp(log_init_N);
  R(startIndex)=exp(log_init_R);

  for(y=startIndex; y<=endIndex+projSteps; y++) {

    //spawners also include some animals that were recruits in the beginning of the year
    SP(y)=sratio*(N(y)*exp(-sp_time*(Mt(2,y)+sel(2)*F(y))))*pSpawn(2)
      + sratio*(R(y)*exp(-sp_time*(Mt(1,y)+sel(1)*F(y))))*pSpawn(1);

    if (SRSwitch==1)	{ //Beverton-Holt
      // don't use recruit deviations for projection years
      if (y<endIndex) R(y+1)=(SP(y)/(SP(y)*beta+alpha))*exp(sr_beta_env*env(envRecTS,y+1-envRecLag))*exp(log_rec_dev(y+1)-biasAdj*0.5*rec_var);
      else R(y+1)=(SP(y)/(SP(y)*beta+alpha))*exp(sr_beta_env*env(envRecTS,y+1-envRecLag));
    }

    if (SRSwitch==2)	{ //Ricker
      // don't use recruit deviations for projection years
      if (y<endIndex) R(y+1)=(alpha*SP(y)*exp(-beta*SP(y)))*exp(sr_beta_env*env(envRecTS,y+1-envRecLag))*exp(log_rec_dev(y+1)-biasAdj*0.5*rec_var);
      else R(y+1)=(alpha*SP(y)*exp(-beta*SP(y)))*exp(sr_beta_env*env(envRecTS,y+1-envRecLag));
    }

    //abundance for the next year
    N(y+1)=R(y)*exp(-(Mt(1,y)+sel(1)*F(y)))+N(y)*exp(-(Mt(2,y)+sel(2)*F(y)));

    //Baranov catch equation
    TC(y)=N(y)*((sel(2)*F(y))/(sel(2)*F(y)+Mt(2,y))*(1.-exp(-(sel(2)*F(y)+Mt(2,y)))))
      + R(y)*((sel(1)*F(y))/(sel(1)*F(y)+Mt(1,y))*(1.-exp(-(sel(1)*F(y)+Mt(1,y)))));

    u(y)=TC(y)/(R(y)*((1-exp(-sel(1)*F(y)))/(1-exp(-F(y))))+N(y));


  }


  //Calculate year-dependent F/N Reference Point components (i.e., ratios)
  NCurr=mfexp((log(N(endIndex))+log(N(endIndex-1))+log(N(endIndex-2)))/3);
  FCurr=mfexp((log(F(endIndex-1))+log(F(endIndex-2)))/2);
  FMSYRatio=FCurr/FMSY;
  NMSYRatio=NCurr/NMSY;
  UMSYRatio=mfexp((log(u(endIndex))+log(u(endIndex-1))+log(u(endIndex-2)))/3)/uMSY;
  OFL=uMSY*mfexp((log(N(endIndex))+log(N(endIndex-1))+log(N(endIndex-2)))/3);  //N is by region, so need to take mean if more

  cLim=max(1-M(2),0.5);
  NLim=cLim*NMSY;
  FLim=FMSY;
  if (NCurr <= NLim) FLim=(FMSY*NCurr)/(cLim*NMSY);
  FFLimRatio=FCurr/FLim;
  NNLimRatio=NCurr/NLim;

  MCurr(1)=mfexp((log(Mt(1,endIndex))+log(Mt(1,endIndex-1))+log(Mt(1,endIndex-2)))/3);
  MCurr(2)=mfexp((log(Mt(2,endIndex))+log(Mt(2,endIndex-1))+log(Mt(2,endIndex-2)))/3);

  SPR0=sratio*pSpawn(2)*exp(-(MCurr(1)+sp_time*MCurr(2)))/(1.-exp(-(MCurr(2))))
    + sratio*pSpawn(1)*exp(-(sp_time*MCurr(1)));
  SPR1=sratio*pSpawn(2)*exp(-(MCurr(1)+sel(1)*FCurr+sp_time*(MCurr(2)+sel(2)*FCurr)))/(1.-exp(-(MCurr(2)+sel(2)*FCurr)))
    + sratio*pSpawn(1)*exp(-sp_time*(MCurr(1)+sel(1)*FCurr));
  SPRCurr=SPR1/SPR0;


FUNCTION calculate_predicted_indices


  //########## Recruits #########
  for (i=1; i<=numRecSurv; i++){

    qr(i)=0.0;
    double counter=0.0;

    for(y=startIndex; y<=endIndex; y++) {
      if (y<recStartIndex) continue;

      if(re_survey_obs(i,y)!=-999.) { //check to make sure year is not missing
        if(!last_phase()) {
          //small constant added to recruitment in earlier stages to 
          //increase numerical stability
          //NOTE: this formulation assumes survey occurs before vulnerable to harvest
          qr(i)+=log(re_survey_obs(i,y))-log(R(y)*exp(-(sr_time(i)*Mt(1,y)))+.000001);
        }
        else { //small constant not included in last estimation stage
          qr(i)+=log(re_survey_obs(i,y))-log(R(y)*exp(-(sr_time(i)*Mt(1,y))));
        }
        counter++;
      }
    }
    //calculate geometric mean
    qr(i)=exp(qr(i)/counter);
    //Calculate predicted index of abundance
    //NOTE: this formulation assumes survey occurs before vulnerable to harvest
    for(y=startIndex; y<=endIndex; y++) {
      re_survey_est(i,y)=qr(i)*(R(y)*exp(-(sr_time(i)*Mt(1,y))));
    }
  }


  //########## Adults #########
  for (i=1; i<=numAdSurv; i++){
    //calculate catchability for each sex-index combination
    double counter=0.0;
    qa(i)=0.0;

    for(y=startIndex; y<=endIndex; y++) {
      if (y<adStartIndex) continue;

      if(ad_survey_obs(i,y)!=-999.) { //check to make sure year is not missing

        qa(i)+=log(ad_survey_obs(i,y))-log(N(y)*exp(-sa_time(i)*(Mt(2,y)+sel(2)*F(y))));
        counter++;
      }
    }
    //calculate geometric mean
    qa(i)=exp(qa(i)/counter);

    //Calculate each predicted index of abundance
    for(y=startIndex; y<=endIndex; y++) {
      ad_survey_est(i,y)=qa(i)*(N(y)*exp(-sa_time(i)*(Mt(2,y)+sel(2)*F(y))));
    }
  }



FUNCTION calculate_objective_function

  double pi=3.141593;

  //calculate adult survey likelihood component
  for (i=1; i<=numAdSurv; i++){
    Lsa(i)=0.0;
    for(y=startIndex; y<=endIndex; y++) {
      if (y<adStartIndex) continue;
      if(ad_survey_obs(i,y)!=-999.) { //check to make sure year is not missing -- some holes
        Lsa(i)+=0.5*log(2.*pi)+0.5*log(ad_survey_var(i,y))+log(ad_survey_obs(i,y))+square(log(ad_survey_obs(i,y)+.000001)-log(ad_survey_est(i,y)+.000001))/(2*ad_survey_var(i,y));
      }
    }
    Lsa(i)=sa_lambda(i)*Lsa(i);
  }


  //calculate recruit survey likelihood component
  for (i=1; i<=numRecSurv; i++){
    Lsr(i)=0.0;
    for(y=startIndex; y<=endIndex; y++) {
      if (y<recStartIndex) continue;
      if(re_survey_obs(i,y)!=-999.) { //check to make sure year is not missing
        Lsr(i)+=0.5*log(2.*pi)+0.5*log(re_survey_var(i,y))+log(re_survey_obs(i,y))+square(log(re_survey_obs(i,y)+.000001)-log(re_survey_est(i,y)+.000001))/(2*re_survey_var(i,y));
      }
    }
    Lsr(i)=sr_lambda(i)*Lsr(i);
  }


  //calculate total catch likelihood component
  Lc=0.0;
  for(y=startIndex; y<=endIndex; y++) {
    if (y<cStartIndex) continue;
    if(TC_obs(y)!=-999.) { //check to make sure year is not missing
      Lc+=square(log(TC_obs(y)+.000001)-log(TC(y)+.000001));
    }
  }
  Lc=com_lambda*(0.5*log(2.*pi)*size_count(TC_obs(startIndex,endIndex))+0.5*log(C_var)*size_count(TC_obs(startIndex,endIndex))+sum(log(TC_obs(startIndex,endIndex)))+0.5*Lc/C_var);


  //calculate likelihood component for recruitment deviations
  Lrdev=recDev_lambda*(0.5*log(2.*pi)*size_count(log_rec_dev)+0.5*log(rec_var)*size_count(log_rec_dev)+sum(log_rec_dev)+0.5*norm2(log_rec_dev)/rec_var);


  //calculate likelihood component for effort residuals if effort time series is included
  Leff=0.0;
  if (effFlag==1) {
    for(y=startIndex; y<endIndex; y++) {
      Leff+=0.5*log(2.*pi)+0.5*log(eff_var)+log(effort(y))+0.5*square(log(effort(y))-(log(F(y))-log_F_q))/eff_var;
    }
    Leff=effort_lambda*Leff;
  }

  //calculate likelihood component for total Z of adults as prior, read from independent Z estimate
  Lz=aveZ_lambda*(0.5*log(2.*pi)+0.5*log(Z_var)+log(aveZ)+0.5*square(log(aveZ)-log(mean(Z(startIndex,endIndex-1))))/Z_var);

  negLL+=sum(Lsa)+sum(Lsr)+Lc+Lrdev+Leff+Lz;


FUNCTION calculate_reference_points

  //Reference point variables
  MSY=0.0;
  u1MSY=0.0;
  u0MSY=0.0;
  uMSY=0.0;
  i=0;
  OFL=0.0;
  SPRDiff=1e10;

  //With recruit spawners 
  SPR0=sratio*pSpawn(2)*(exp(-(Myr(1)+sp_time*Myr(2)))/(1.-exp(-(Myr(2)))))
    + sratio*pSpawn(1)*(exp(-((sp_time*Myr(1)))));

  Fval=Fval_init;
  for(i=1; i<=Fval_num; i++)
  {
    Fvec(i)=Fval; //record the F values

    SPR(i)=sratio*pSpawn(2)*(exp(-(Myr(1)+sel(1)*Fval+sp_time*(Myr(2)+sel(2)*Fval)))/(1.-exp(-(Myr(2)+sel(2)*Fval))))
      + sratio*pSpawn(1)*(exp(-(sp_time*(Myr(1)+sel(1)*Fval))));

    NPR(i)=exp(-(Myr(1)+sel(1)*Fval))/(1.-exp(-(Myr(2)+sel(2)*Fval)));

    YPR(i)=(sel(1)*Fval)/(sel(1)*Fval+Myr(1))*(1.-exp(-(sel(1)*Fval+Myr(1))))
      + ((sel(2)*Fval)/(sel(2)*Fval+Myr(2))*(1.-exp(-(sel(2)*Fval+Myr(2)))))*NPR(i);

    if (SRSwitch==1)	 R_eq(i)=(SPR(i)-alpha)/(SPR(i)*beta);
    if (SRSwitch==2)	R_eq(i)=(log(alpha)+log(SPR(i)))/(beta*SPR(i));

    N_eq(i) = NPR(i)*R_eq(i);
    C_eq(i)=YPR(i)*R_eq(i);

    //calculate exploitation rate
    //age 0+
    u0_eq(i)=(sel(1)*Fval)/(sel(1)*Fval+Myr(1))*(1.-exp(-(sel(1)*Fval+Myr(1))));
    //age 1+
    u1_eq(i)=(sel(2)*Fval)/(sel(2)*Fval+Myr(2))*(1.-exp(-(sel(2)*Fval+Myr(2))));
    //all ages
    if (i>1) uAll_eq(i)=C_eq(i)/(N_eq(i)+R_eq(i)*((1-exp(-sel(1)*Fval))/(1-exp(-Fval))));

    //MSY
    if (C_eq(i)>MSY) {
      MSY=C_eq(i);
      FMSY=Fval;
      NMSY=N_eq(i);
      RMSY=R_eq(i);
      u0MSY=u0_eq(i);
      u1MSY=u1_eq(i);
      uMSY=uAll_eq(i);
    }


    //loop through SPR targets and see if at the correct F for each target
    for (ispr=1; ispr<=nspr; ispr++){
      if (square(SPR(i)/SPR0-SPR_targ(ispr)) < SPRDiff(ispr)) {
        SPRDiff(ispr)=square(SPR(i)/SPR0-SPR_targ(ispr));
        FSPR_ref(ispr)=Fval;
      }
    }

    //increment the female F for the SPR
    Fval+=Fval_inc;
  }


//############### Reporting functions ###############//

FUNCTION mcmc
  //Code to write results of MCMC to file so we can access the chains
  if(mceval_phase()) {
    //Define output file stream for MCMC results
    if(iterMCMC==0) {
      ofstream mcmcout("cmsa_refs.mcmc");
      mcmcout <<"MSY\t"<<"FMSY\t"<< "NMSY\t"<<"uMSY\t"<<"FLim\t"<<"NLim\t"<<"FMSYRatio\t"<<"NMSYRatio\t"<<"UMSYRatio\t"<<"FFLimRatio\t"<<"NNLimRatio\t"<<"SPRCurrent"<<endl;

      //print out yearly F and N
      ofstream mcmcout2("cmsa_yearly.mcmc");

      year=fyear;
      for(y=startIndex;y<=endIndex;y++){
        mcmcout2 <<"N"<<year<<"\t";
        year=year+1.0/timeSteps;
      }

      year=fyear;
      for(y=startIndex;y<=endIndex;y++){
        mcmcout2 <<"R"<<year<<"\t";
        year=year+1.0/timeSteps;
      }

      year=fyear;
      for(y=startIndex;y<=endIndex;y++){
        if( y<endIndex) mcmcout2 <<"F"<<year<<"\t";
        else mcmcout2 <<"F"<<year << endl;
        year=year+1.0/timeSteps;
      }

      ofstream mcmcout3("cmsa_pars.mcmc");
      mcmcout3 <<"N0\t"<<"R0\t"<< "Fq\t"<<"S0\t"<<"h"<<endl;

      iterMCMC++;
    }


    ofstream mcmcout("cmsa_refs.mcmc",ios::app);
    mcmcout <<MSY<<"\t"<<FMSY<<"\t"<<NMSY<<"\t"<<uMSY<<"\t"<<FLim<<"\t"<<NLim<<"\t"<<FMSYRatio<<"\t"<<NMSYRatio<<"\t"<<UMSYRatio<<"\t"<<FFLimRatio<<"\t"<<NNLimRatio<<"\t"<<SPRCurr<<endl;

    //print out yearly F and N
    ofstream mcmcout2("cmsa_yearly.mcmc",ios::app);

    for(y=startIndex;y<=endIndex;y++){
      mcmcout2 <<N(y) << "\t";
    }

    for(y=startIndex;y<=endIndex;y++){
      mcmcout2 <<R(y) << "\t";
    }

    for(y=startIndex;y<=endIndex;y++){
      if (y<endIndex) mcmcout2 <<F(y) << "\t";
      else mcmcout2 <<F(y) <<endl;
    }

    //print out yearly F and N
    ofstream mcmcout3("cmsa_pars.mcmc",ios::app);
    mcmcout3 <<exp(log_init_N)<<"\t"<<exp(log_init_R) <<"\t"<<exp(log_F_q)<<"\t"<<exp(log_S0)<<"\t"<<exp(log_steep) << endl;

  }


FUNCTION calculate_tSPR

  ofstream ofs_tSPR("tSPR.dat");
  ofs_tSPR<< "year\ttSPR" << endl;
  year=fyear;

  for(y=startIndex;y<=endIndex+projSteps;y++){

    SPR0=sratio*pSpawn(2)*exp(-(Mt(1,y)+sp_time*Mt(2,y)))/(1.-exp(-(Mt(2,y))))
      + sratio*pSpawn(1)*exp(-(sp_time*Mt(1,y)));
    SPR1=sratio*pSpawn(2)*exp(-(Mt(1,y)+sel(1)*F(y)+sp_time*(Mt(2,y)+sel(2)*F(y))))/(1.-exp(-(Mt(2,y)+sel(2)*F(y))))
      + sratio*pSpawn(1)*exp(-sp_time*(Mt(1,y)+sel(1)*F(y)));

    tSPR(y)=SPR1/SPR0;
    ofs_tSPR << year << "\t" << tSPR(y) << endl;
    year=year+1.0/timeSteps;
  }



FUNCTION obs_pred

  ofstream ofs_op("obs_pred_results.dat");
  ofs_op << "survey year sex a_r s_c snum obs pred" << endl;
  year=fyear;

  for(y=startIndex; y<=endIndex; y++) {

    //total observed and predicted catch
    ofs_op << "0 "<< year << " t a c 0 " << TC_obs(y) << " " << TC(y) <<  endl;

    //adult surveys
    for (i=1; i<=numAdSurv; i++)
    ofs_op << i << " "<< year << " 0 a s 0 " << ad_survey_obs(i,y) << " " << ad_survey_est(i,y) << endl;

    //recruit surveys
    for (i=1; i<=numRecSurv; i++)
    ofs_op  << i << " "<< year << " 0 r s 0 " << re_survey_obs(i,y) << " " << re_survey_est(i,y) << endl;

    if (y==startIndex) ofs_op << "0 "<< year << " r r r 0 " <<  R(y) << " " << "NA"  <<  endl;
    else {
      if (SRSwitch==1) ofs_op << "0 "<< year << " r r r 0 " << R(y)<< " "   << SP(y-1)/(SP(y-1)*beta+alpha)  <<  endl;
      if (SRSwitch==2) ofs_op << "0 "<< year << " r r r 0 " << R(y)<< " "   << alpha*SP(y-1)*exp(-beta*(SP(y-1)))  <<  endl;
    } //recruitment deviations

    year=year+1.0/timeSteps;
  }


FUNCTION HPD_estimates
  ofstream ofs_hpd("HPD_results.dat");
  ofs_hpd << "year Adult Spawners Rec RecSurvey1 TC recM adM  F FMSYRatio NMSYRatio FFLimRatio NNLimRatio u0 u1 uAll SREnv MEnvRec MEnvAd" << endl;

  year=fyear; //for outputting the year if multiple time steps per year

  for(y=startIndex;y<=endIndex+projSteps;y++){
    ofs_hpd << year << " " << N(y) << " "<< SP(y)<< " " << R(y) << " " <<R(y)*exp(-sr_time(1)*Mt(1,y)) << " " << TC(y) << " " << Mt(1,y) << " " << Mt(2,y) << " " << F(y)<< " " << F(y)/FMSY<< " " << N(y)/NMSY<< " "<< F(y)/FLim<< " " << N(y)/NLim<< " "<< (sel(1)*F(y))/(sel(1)*F(y)+Mt(1,y))*(1.-exp(-(sel(1)*F(y)+Mt(1,y)))) << " "<< (sel(2)*F(y))/(sel(2)*F(y)+Mt(2,y))*(1.-exp(-(sel(2)*F(y)+Mt(2,y)))) << " "<< u(y) << " " << env(envRecTS,y-envRecLag) << " "<< env(envMTS(1),y-envMLag(1)) << " " << env(envMTS(2),y-envMLag(2)) << endl;
    year=year+1.0/timeSteps;
  }


FUNCTION MSY_estimates
  ofstream ofs_msy("MSY_results.dat");
  {
    //Column headings
    ofs_msy << "Fval\t" << "C_eq\t" << "N_eq\t" << "R_eq\t" << "YPR\t" << "SPR\t" << "SPRatio\t" << "u0_eq\t" << "u1_eq\t" << "uAll_eq\t" << endl;

    Fval=Fval_init;
    for(i=1; i<=Fval_num; i++) {
      ofs_msy << Fval << "\t" << C_eq(i) << "\t" << N_eq(i) << "\t" << R_eq(i) << "\t" <<YPR(i) << "\t" << SPR(i) << "\t" << SPR(i)/SPR0 << "\t" << u0_eq(i) << "\t" << u1_eq(i)<< "\t" << uAll_eq(i) << endl;
      Fval+=Fval_inc;
    }
  }

FUNCTION general_report
  ofstream ofs_gen("gen_results.dat");
  {
    ofs_gen << "Name Value" << endl;
    ofs_gen << "negLL " << negLL <<endl;
    ofs_gen << "Lsa "<< Lsa <<endl;
    ofs_gen << "Lsr "<< Lsr << endl;
    ofs_gen << "Lc "<<  Lc << endl;
    ofs_gen << "Lz "<<  Lz << endl;
    ofs_gen << "Lrdev "<< Lrdev << endl;
    ofs_gen << "Leff "<< Leff << endl;
    ofs_gen << "init_N "<< exp(log_init_N) << endl;
    ofs_gen << "init_R " << exp(log_init_R) << endl;
    for (i=1; i<=numAdSurv; i++) ofs_gen << "qa_i"<<i <<" " << qa(i) << endl;
    for (i=1; i<=numRecSurv; i++) ofs_gen << "qr_i"<<i <<" " << qr(i) << endl;
    ofs_gen << "F_q "<< exp(log_F_q) << endl;
    ofs_gen << "rec_cv "<< exp(log_rec_cv) << endl;
    ofs_gen << "p_rec "<< p_rec << endl;
    ofs_gen << "p_under " << p_under << endl;

    ofs_gen << "SRType " << SRSwitch << endl;
    ofs_gen << "S0 " << exp(log_S0) << endl;
    ofs_gen << "steepness " << exp(log_steep) << endl;
    ofs_gen << "alpha " << alpha << endl;
    ofs_gen << "beta " << beta << endl;
    ofs_gen << "sr_beta_env " << sr_beta_env << endl;
    ofs_gen << "M_beta_env_1 " << M_beta_env(1) << endl;
    ofs_gen << "M_beta_env_2 " << M_beta_env(2) << endl;
    ofs_gen << "sel_1 " << exp(log_sel(1)) << endl;
    ofs_gen << "sel_2 " << exp(log_sel(2)) << endl;
    ofs_gen << "Mr " << Myr(1) << endl;
    ofs_gen << "Ma " << Myr(2) << endl;
    ofs_gen << "Ma " << Myr(2) << endl;
    ofs_gen << "sp_time " << sp_time <<endl;
    ofs_gen  << "MSY " << MSY << endl;
    ofs_gen  << "FMSY " << FMSY << endl;
    ofs_gen  << "FMSYRatio " << FMSYRatio << endl;
    ofs_gen  << "NMSY " << NMSY << endl;
    ofs_gen  << "NMSYRatio " << NMSYRatio << endl;
    ofs_gen  << "RMSY " << RMSY << endl;
    ofs_gen  << "u0MSY " << u0MSY << endl;
    ofs_gen  << "u1MSY " << u1MSY<< endl;
    ofs_gen  << "uMSY " << uMSY<< endl;
    ofs_gen  << "UMSYRatio " << UMSYRatio << endl;
    ofs_gen  << "FLim " << FLim << endl;
    ofs_gen  << "FFLimRatio " << FFLimRatio << endl;
    ofs_gen  << "NLim " << NLim << endl;
    ofs_gen  << "NNLimRatio " << NNLimRatio << endl;
    ofs_gen  << "cLim " << cLim << endl;
    ofs_gen  << "OFL " << OFL<< endl;
    ofs_gen  << "projYears " << projYears<< endl;


    for(ispr=1; ispr<=nspr; ispr++) {
      ofs_gen  << "F"<<SPR_targ(ispr) << "% " << FSPR_ref(ispr)<< endl;
    }
  }


//###################################################################
//###################################################################
REPORT_SECTION
//Call reporting functions
  calculate_tSPR();
  obs_pred();
  MSY_estimates();
  HPD_estimates();
  general_report();

  report << "Likelihood Components" <<endl;
  report << "negLL\t" << negLL <<endl;
  report << "Lsa\t" << Lsa <<endl;
  report << "Lsr\t" << Lsr << endl;
  report << "Lc\t" <<  Lc << endl;
  report << "Leff\t" <<  Leff << endl;
  report << "Lz\t" <<  Lz << endl;
  report << "Lrdev\t" << Lrdev << endl;
  report << "F_pen\t" << F_pen << endl;
  report << "\nParameter Estimates (NOT log space unless marked)" <<endl;
  report << "init_N\t" << exp(log_init_N) << endl;
  report << "init_R\t" << exp(log_init_R) << endl;
  report << "F\t" << F << endl;
  report << "M_rec\t" << Mt(1) << endl;
  report << "M_ad\t" << Mt(2) << endl;
  report << "AveF\t" << mean(F(startIndex,endIndex-1)) << endl;
  report << "AveZ\t" << mean(Z(startIndex,endIndex-1)) << endl;
  report << "AveU\t" << mean(u(startIndex,endIndex)) << endl;
  report << "F_q "<< exp(log_F_q) << endl;
  report << "log_F_dev\t" << log_F_dev << endl;
  report << "mean(log_F_dev)\t" << mean(log_F_dev) << endl;
  report << "rec_dev\t" << exp(log_rec_dev) << endl;
  report << "mean(log_rec_dev)\t" << mean(log_rec_dev) << endl;
  report << "mean(rec_dev)\t" << mean(exp(log_rec_dev)) << endl;
  report << "rec_cv\t" << exp(log_rec_cv) << endl;
  for (i=1; i<=numAdSurv; i++) report << "qa_i"<<i <<"\t" << qa(i) << endl;
  for (i=1; i<=numRecSurv; i++) report << "qr_i"<<i <<"\t" << qr(i) << endl;
  if (SRSwitch==1)	 //Beverton-Holt
  report << "SR=Beverton-Holt\t" << endl;
  if (SRSwitch==2)	 //Ricker
  report << "SR=Ricker\t" << endl;
  report << "S0\t" << exp(log_S0) << endl;
  report << "steepness\t" << exp(log_steep) << endl;
  report << "alpha\t" << alpha << endl;
  report << "beta\t" << beta << endl;
  report << "sr_beta_env\t" << sr_beta_env << endl;
  report << "M_beta_env\t" << M_beta_env << endl;
  report << "M_beta_env_rec\t" << M_beta_env(1) << endl;
  report << "M_beta_env_ad\t" << M_beta_env(2) << endl;
  report << "sel\t" << exp(log_sel) << endl;
  report << "p_rec\t" << p_rec << endl;
  report << "p_under\t" << p_under << endl;
  report << "M\t" << Myr << endl;
  report << "rec_cv\t" << exp(log_rec_cv) << endl;
  report << "sp_time\t" << sp_time <<endl;

  report << "\nReference Point Calculations" <<endl;
  report << "negLL\t" << negLL <<endl;
  report  << "MSY\t" << MSY << endl;
  report  << "uMSY\t" << uMSY<< endl;
  report  << "NMSY " << NMSY << endl;
  report  << "UMSYRatio " << UMSYRatio << endl;
  report  << "NMSYRatio " << NMSYRatio << endl;
  report  << "FMSY\t" << FMSY << endl;
  report  << "FMSYRatio " << FMSYRatio << endl;
  report  << "RMSY\t" << RMSY << endl;
  report  << "u0MSY\t" << u0MSY << endl;
  report  << "u1MSY\t" << u1MSY<< endl;
  report  << "FLim\t" << FLim << endl;
  report  << "FFLimRatio\t" << FFLimRatio << endl;
  report  << "NLim\t" << NLim << endl;
  report  << "NNLimRatio\t" << NNLimRatio << endl;
  report  << "cLim\t" << cLim << endl;
  report  << "OFL\t" << OFL << endl;
  report  << "SPRCurr\t" << SPRCurr<< endl;

  report <<"\n"<< negLL <<endl;
  report  <<MSY << endl;
  report  <<FMSY<< endl;
  report  <<NMSY << endl;
  report  <<FFLimRatio << endl;
  report  <<NNLimRatio << "\n"<< endl;

  for(ispr=1; ispr<=nspr; ispr++) {
    report  << "F"<<SPR_targ(ispr) << "%\t" << FSPR_ref(ispr)<< endl;
  }

//###################################################################
//###################################################################
GLOBALS_SECTION

  #include "admodel.h"
  //define constant variable


//###################################################################
//###################################################################
RUNTIME_SECTION
  maximum_function_evaluations 5000,25000,20000,20000,20000,20000
  convergence_criteria 1.0e-8

  //Leave space below this line

