#### Author of this code: Eva-Maria Ahrer, adapted from Tiberius TransitGPPM model (author: J. Kirk)

import numpy as np

# from fitting_utils import lightcurve

import dynesty
import emcee
import copy

import sys
import pickle

from scipy import optimize,stats

import matplotlib.pyplot as plt
import multiprocessing as mp

# from fitting_utils import parametric_fitting_functions as pf
from fitting_utils import plotting_utils as pu
from fitting_utils import priors
from fitting_utils import joint_fitting as jf

from dynesty import plotting as dyplot
from dynesty.utils import resample_equal

class Sampling(object):
    def __init__(self, lightcurve_list, sampling_arguments, sampling_method, output_foldername):

        """


        Inputs:
        lightcurve_list    - list, light curve classes which includes the full model (transit, systematics, etc.)
        sampling_arguments - dict, parameters needed for dynesty / emcee; e.g. live points, precision criterion, nsteps, nwalkers
        sampling_method    - str, either dynesty, emcee, LM



        Can return:
        - dynesty result
        - emcee result
        """

        self.lightcurve_list = lightcurve_list
        self.nlightcurves    = len(lightcurve_list)
        self.sampling_method = sampling_method
        self.sampling_arguments = sampling_arguments
        self.output_foldername = output_foldername

        if sampling_method =='emcee' or sampling_method =='LM' or sampling_method == 'dynesty':

            self.joint_fitter = jf.JointFitter(lightcurve_list,verbose=True)

            self.param_dict      = self.joint_fitter.param_dict_global
            self.prior_dict      = self.joint_fitter.prior_dict_global
            self.param_list_free = self.joint_fitter.param_list_global

            if self.sampling_method == 'dynesty':
                self.nDims = len(self.param_list_free)

        else:
            self.lightcurve = lightcurve_list[0]
            self.param_dict      = self.lightcurve.param_dict
            self.param_list_free = self.lightcurve.param_list_free
            self.prior_dict      = self.lightcurve.prior_dict

    # -------------------- Dynesty methods -------------------- #
    
    def prior_setup(self, x):

        if self.sampling_method == 'dynesty':

            theta = [0] * self.nDims

            for i in range(self.nDims):
                if self.prior_dict['%s_prior'%self.param_list_free[i]] == 'N':
                    theta[i] = priors.GaussianPrior(self.prior_dict['%s_1'%self.param_list_free[i]], self.prior_dict['%s_2'%self.param_list_free[i]])(np.array(x[i]))
                elif self.prior_dict['%s_prior'%self.param_list_free[i]] == 'U':
                    theta[i] = priors.UniformPrior(self.prior_dict['%s_1'%self.param_list_free[i]], self.prior_dict['%s_2'%self.param_list_free[i]])(np.array(x[i]))
            return theta


    def loglikelihood(self,theta):

        logL_total = 0.0

        for ilc, lc in enumerate(self.lightcurve_list):
              
            if len(self.lightcurve_list)>1:
                lc_theta,lc_param_list = self.joint_fitter.map_theta(theta,ilc,self.param_list_free) # e.g. lc_param_list (t0,c1_lc0,c2_lc0)
            else:
                lc_theta = theta

            lc.update_model(lc_theta)
            noise = lc.return_flux_err()

            if lc.GP_used:
                model_calc = lc.calc(with_GP=False)
                logL = lc.GP_model.lnlike(model_calc,noise)
            else:
                residuals  = lc.calc_residuals()
                
                N = len(noise)
                logL = -N/2. * np.log(2*np.pi) - np.nansum(np.log(noise)) - np.nansum(residuals**2 / (2*noise**2))
            
            logL_total += logL
        
        return logL_total


    def run_dynesty(self, wavelength_bin=0, median_update=True):

        # For statistics
        self.best_fit_median = median_update

        namelist = self.param_list_free

        live_points = self.sampling_arguments['nlive_pdim']
        precision_criterion = self.sampling_arguments['precision_crit']
        sampler = dynesty.NestedSampler(self.loglikelihood, 
                                        self.prior_setup, 
                                        self.nDims,
                                        nlive=live_points*self.nDims, 
                                        bootstrap=0) #,sample='rslice')
        sampler.run_nested(dlogz=precision_criterion, print_progress=True)
        self.sampler_results = sampler.results

        # Get equal weight sampels
        self.samples = self.sampler_results.samples_equal()

        # Get sample with highest likelihood
        logl = self.sampler_results['logl']
        best_idx = np.argmax(logl)
        self.best_fit_values_logl = self.sampler_results.samples[best_idx]

        # Get median values
        _, _, med, _, _ = self.get_confidence_interval(wavelength_bin)
        self.best_fit_values_median = med

        # Decide which values to update the fitted lightcurve
        if median_update:
            best_val = self.best_fit_values_median
        else:
            best_val = self.best_fit_values_logl

        # Update fitted lightcurve
        for ilc,lc in enumerate(self.lightcurve_list):
            if len(self.lightcurve_list)>1:
                lc_theta, _ = self.joint_fitter.map_theta(best_val,ilc,namelist)
            else:
                lc_theta = best_val
            lc.update_model(lc_theta)

        return self.lightcurve_list

    # -------------------- EMCEE methods -------------------- #
    def logprior_emcee(self, theta):
        """Compute joint log prior based on global prior_dict."""
        lp = 0
        for i, pname in enumerate(self.param_list_free):
            prior_type = self.prior_dict[f'{pname}_prior']
            if prior_type == 'N':
                mu = self.prior_dict[f'{pname}_1']
                sigma = self.prior_dict[f'{pname}_2']
                lp += -0.5 * ((theta[i] - mu)/sigma)**2 - np.log(sigma*np.sqrt(2*np.pi))
            elif prior_type == 'U':
                lower = self.prior_dict[f'{pname}_1']
                upper = self.prior_dict[f'{pname}_2']
                if not (lower <= theta[i] <= upper):
                    return -np.inf
        return lp

    def logprobability_emcee(self, theta):
        """Full log-probability for emcee: lnprior + lnlike."""
        lp = self.logprior_emcee(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + self.loglikelihood(theta)

    def advance_chain(self,sampler,p0,nsteps,burn,save_chain,wavelength_bin):
        """The function that advances the emcee sampler chain with a progress bar

        Inputs:
        sampler - the emcee sampler, intitiated in run_emcee
        p0 - the array of (starting) parameter nvalues
        nsteps - the number of steps to advance the chain over
        burn - is this a burn chain? If so, don't save to file
        save_chain - True/False, do we want to save the chain to file?
        wavelength_bin - the number of the wavelength bin we're fitting, so that we can save the output correctly

        Returns:
        sampler - the inputted emcee sampler advanced by nsteps"""

        width = 100 # for progress bar
        highest_prob = 0
        print('Progress:') # for progress bar
        for i,(pos, prob, state) in enumerate(sampler.sample(p0,iterations=nsteps,store=True)):
            n = int((width+1) * float(i) / nsteps)
            sys.stdout.write("\r[{0}{1}]".format('#' * n, ' ' * (width - n))) # for progress bar
            
            if np.max(prob) > highest_prob:
                highest_prob_pars = pos[np.argmax(prob)]
                highest_prob = np.max(prob)

            if not burn and save_chain:
                f = open('prod_chain_wb%s.txt'%(str(wavelength_bin+1).zfill(4)),'a')
            for k in range(pos.shape[0]):
                # loop over all walkers and append to file
                thisPos = pos[k]
                thisProb = prob[k]
                if not burn and save_chain: # only save production chain to file and only file steps otherwise these files are huge!
                    if nsteps > 500 and i > nsteps/2.:
                        f.write("{0:4d} {1:s} {2:f}\n".format(k," ".join(map(str,thisPos)),thisProb ))
                    if nsteps < 500 and i > nsteps - 100:
                        f.write("{0:4d} {1:s} {2:f}\n".format(k," ".join(map(str,thisPos)),thisProb ))
            if not burn and save_chain:
                f.close()

        return sampler, highest_prob_pars, highest_prob

    def run_emcee(self, burn=False, save_chain=True, wavelength_bin=0):

        """Run emcee MCMC sampling."""
        nsteps   = self.sampling_arguments['nsteps']
        ndiscard = self.sampling_arguments['ndiscard']
        nwalkers = self.sampling_arguments['nwalkers']
        nthreads = self.sampling_arguments['nthreads']
        npars    = len(self.param_list_free)
        namelist = self.param_list_free
        nwalkers_total = nwalkers * npars
        self.burn = burn

        # Scatter walkers around starting parameters
        starting_values = np.array([self.param_dict[k].currVal for k in self.param_list_free])

        if burn:
            p0 = emcee.utils.sample_ball(starting_values, 1e-3*starting_values, size=nwalkers_total)
        else:
            p0 = [starting_values + 1e-8*np.random.randn(npars) for j in range(nwalkers_total)]

         # intiate emcee sampler object
        if npars > 1:
            sampler = emcee.EnsembleSampler(nwalkers_total,npars,self.logprobability_emcee,threads=nthreads)
        else: # from my own tests I find that for a single parameter, the acceptance fraction is too high. Increasing the stretch scale factor decreases the acceptance fraction to a more acceptable value. This is relevant for ingress/egress fitting for ingress/egress with just Rp/Rs
            sampler = emcee.EnsembleSampler(nwalkers_total,npars,self.logprobability_emcee,threads=nthreads,
                                            moves=emcee.moves.StretchMove(10))

        # run chains
        print('################')
        if burn:
            print("Running burn-in for bin %d..."%(wavelength_bin+1))
            # f = open('burn_chain_%d.txt'%(wavelength_bin+1),'w') # deciding to only save production chain

        else:
            print("Running production for bin %d..."%(wavelength_bin+1))
            if save_chain:
                f = open('prod_chain_wb%s.txt'%(str(wavelength_bin+1).zfill(4)),'w')
                f.close()

        if nsteps == "auto":
            not_converged = True
            chain_number = 1
            nsteps = 2000

            while not_converged:

                if chain_number == 1:
                    sampler, highest_prob_pars, highest_prob = self.advance_chain(sampler,p0,nsteps,burn,save_chain,wavelength_bin)
                else:
                    sampler, highest_prob_pars, highest_prob = self.advance_chain(sampler,sampler.get_last_sample(),nsteps,burn,save_chain,wavelength_bin)

                total_steps = chain_number*nsteps

                try:
                    auto_corr_time = np.round(total_steps/np.median(sampler.acor))

                    # ideal scenario, we're >= 50x the median autocorr time
                    if auto_corr_time >= 50: # taking DFM's estimate
                        not_converged = False
                        print("\n\nChains run for %d total steps"%(chain_number*nsteps))
                        nsteps = total_steps # updated nsteps for calculation of corner plots and parameter values later on

                    # not so good scenario but chains are getting long
                    elif auto_corr_time >= 20 and total_steps >= 10000:
                        print("\n\n After %d steps the number of steps is %dX the autocorrelation time, finishing chain"%(chain_number*nsteps,auto_corr_time))
                        nsteps = total_steps # updated nsteps for calculation of corner plots and parameter values later on
                        not_converged = False

                    # chains too long, probably won't converge now
                    elif total_steps >= 20000:
                        print("\n\n After %d steps the chains have not yet converged, exiting"%(chain_number*nsteps))
                        nsteps = total_steps
                        not_converged = False

                    else:
                        print("\n\n After %d steps the number of steps is %dX the autocorrelation time, running chain again"%(chain_number*nsteps,auto_corr_time))
                        chain_number += 1
                except:
                    if total_steps >= 20000:
                        print("\n\n After %d steps the chains have not yet converged, exiting"%(chain_number*nsteps))
                        nsteps = total_steps
                        not_converged = False
                    else:
                        print("\n\n After %d steps the chains have not yet converged, running chain again"%(chain_number*nsteps))
                        chain_number += 1
        else:
            sampler, highest_prob_pars, highest_prob = self.advance_chain(sampler,p0,nsteps,burn,save_chain,wavelength_bin)

        self.sampler_mcmc = sampler

        if nsteps >= 500:
            if burn:
                samples = sampler.chain[:, int(nsteps/2):, :].reshape((-1, npars))
            else:
                # samples = sampler.get_chain(discard=int(nsteps/4), thin=10, flat=True)
                samples = sampler.get_chain(discard=ndiscard, thin=10, flat=True)
        else:
            samples = sampler.chain[:, -100:, :].reshape((-1, npars))

        self.samples = samples

        print('\n')
        _, _, med, _, _ = self.get_confidence_interval(wavelength_bin, equal_prob=True)
        self.best_fit_values_median = med

        for ilc,lc in enumerate(self.lightcurve_list):
            if len(self.lightcurve_list)>1:
                lc_theta, _ = self.joint_fitter.map_theta(med,ilc,namelist)
            else:
                lc_theta = med
            lc.update_model(lc_theta)
            print(f'testing lc 1: {lc.systematic_model.poly_used}')

        if burn:
          print("...burn-in complete for bin %d"%(wavelength_bin+1))
        else:
          for ilc,lc in enumerate(self.lightcurve_list):
                pickle.dump(lc,open(f'fitted_lightcurve_model_lc{ilc}_wb{str(wavelength_bin+1).zfill(4)}.pickle','wb'))
          print("...production complete for bin %d"%(wavelength_bin+1))            

        sampler.reset()

        return self.lightcurve_list
    
    ### -------- Levenberg-Marquadt methods -------- ###
    def run_LM(self,wavelength_bin=0):
        """
        Run Levenberg-Marquardt optimization to minimize residuals
        """

        # Initial parameter vector
        theta0 = np.array([self.param_dict[p].currVal for p in self.param_list_free])

        # Define residual function
        def residuals(theta):
            # Check prior once before updating lightcurves
            prior_val = self.logprior_emcee(theta)
            if not np.isfinite(prior_val):
                return np.ones(sum(len(lc.flux_array) for lc in self.lightcurve_list))*np.inf

            # Update all lightcurves and collect residuals
            all_residuals = []
            for ilc, lc in enumerate(self.lightcurve_list):
                if theta is not None:
                    if len(self.lightcurve_list)>1:
                        lc_theta, lc_param_list = self.joint_fitter.map_theta(theta, ilc, self.param_list_free)
                    else:
                        lc_theta = theta
                    lc.update_model(lc_theta)

                lc_residuals = lc.calc_residuals() / lc.flux_err
                all_residuals.append(lc_residuals)

            return np.concatenate(all_residuals)

        # Run Levenberg-Marquardt fit
        result = optimize.least_squares(residuals, theta0, method='lm')
        self.sampler_results = result

        # Update model with best-fit parameters for all lightcurves
        for ilc, lc in enumerate(self.lightcurve_list):
            if len(self.lightcurve_list)>1:
                lc_theta, _ = self.joint_fitter.map_theta(result.x, ilc, self.param_list_free)
            else:
                lc_theta = result.x
            lc.update_model(lc_theta)

        # Estimate uncertainties from covariance matrix
        try:
            J = result.jac
            cov = np.linalg.inv(J.T.dot(J))*self.reducedChisq()
            uncertainties = np.sqrt(np.diag(cov))
        except:
            print("Unable to estimate uncertainties from Jacobian")
            uncertainties = np.zeros_like(result.x)
        
        return self.lightcurve_list


    def build_bounds(self, full_model=False):
        """
        Build parameter bounds using prior definitions and GP/kernel settings.
        Returns a list of (min, max) tuples in the order of self.param_list_free.
        """
        bnds = []

        for name in self.param_list_free:

            # Use prior definitions if available
            prior_type = self.prior_dict.get(f'{name}_prior', None)
            if prior_type == 'U':
                bnds.append((self.prior_dict[f'{name}_1'], self.prior_dict[f'{name}_2']))
            elif prior_type == 'N':
                # Optional: ±3σ around mean as bounds
                mu, sigma = self.prior_dict[f'{name}_1'], self.prior_dict[f'{name}_2']
                bnds.append((mu - 3*sigma, mu + 3*sigma))
            else:
                # Fallback: ±10% around current value
                curr = self.pars[name].currVal
                bnds.append((curr*0.9, curr*1.1))

        # GP bounds
        if self.lightcurve.GP_used:
            if self.lightcurve.wn_kernel:
                bnds.append((np.log((self.lightcurve.kernel_priors['min_WN_sigma'])**2),
                             np.log((self.lightcurve.kernel_priors['max_WN_sigma'])**2)))

            for i in range(self.gp_ndim):
                for key in [f'min_lniL_{i+1}', f'max_lniL_{i+1}']:
                    if key in self.lightcurve.kernel_priors:
                        bnds.append((self.lightcurve.kernel_priors[key], self.lightcurve.kernel_priors[key]))

            # GP amplitude
            if 'min_A' in self.lightcurve.kernel_priors and 'max_A' in self.lightcurve.kernel_priors:
                bnds.append((self.lightcurve.kernel_priors['min_A'], self.lightcurve.kernel_priors['max_A']))

        return bnds


    # -------------------- Statistical Evaluation Methods -------------------- #
    def chisq(self, theta=None):
        """
        Inputs:
        theta              - array, fitted parameters
        lc_idx (optional)  - int, lightcurve index
        Returns:
        chisq_dict         - dict, 'total' chisquare across all lightcurves, 'lc0' for lightcurve 0 etc..
        """
        
        if theta is not None:
            for ilc, lc in enumerate(self.lightcurve_list):
                if theta is not None:
                    if len(self.lightcurve_list)>1:
                        lc_theta,lc_param_list = self.joint_fitter.map_theta(theta,ilc,self.param_list_free) # e.g. lc_param_list (t0,c1_lc0,c2_lc0)
                    else:
                        lc_theta = theta
                    lc.update_model(lc_theta)

        chisq_dict = {}
        total_chisq = 0.0

        for ilc, lc in enumerate(self.lightcurve_list):
            resids = lc.calc_residuals() / lc.flux_err

            chisq_dict[f'lc{ilc}'] = np.sum(resids**2)
            total_chisq += np.sum(resids**2)

        chisq_dict['total'] = total_chisq
        return chisq_dict

    def reducedChisq(self, theta=None):
        """
        Returns:
        rchisq_dict    - dict, 'total' chisq of all lightcurves, 'lc0' chisq of lightcurve 0 etc..
        """
        chisq_dict  = self.chisq(theta)
        rchisq_dict = {}
        for ilc, lc in enumerate(self.lightcurve_list):
            chisq = chisq_dict[f'lc{ilc}']
            dof   = len(lc.flux_array) - len(lc.param_list_free)
            rchisq_dict[f'lc{ilc}'] = chisq/dof 
        
        tchisq = chisq_dict['total']

        tdof = sum(len(lc.flux_array) for lc in self.lightcurve_list) - len(self.param_list_free)
        rchisq_dict['total'] = tchisq / tdof

        return rchisq_dict

    def rms(self, theta=None):
        """
        Returns:
        rms_dict       - dict, 'lc0' RMS of lightcurve 0 etc.. 
        """
        rms_dict = {}
        if theta is not None:
            for ilc, lc in enumerate(self.lightcurve_list):
                if theta is not None:
                    if len(self.lightcurve_list)>1:
                        lc_theta,lc_param_list = self.joint_fitter.map_theta(theta,ilc,self.param_list_free) # e.g. lc_param_list (t0,c1_lc0,c2_lc0)
                    else:
                        lc_theta = theta 
                    lc.update_model(lc_theta)

        for ilc, lc in enumerate(self.lightcurve_list):
            resids = lc.calc_residuals()
            
            rms_dict[f'lc{ilc}'] = np.sqrt(np.mean(resids**2))

        return rms_dict

    def BIC(self, theta=None):
        # note we can use loglikelihood_emcee also for LM fit since the statistic is independent of the sampling method
        npars = len(self.param_list_free)
    
        if theta is None:
            if self.best_fit_median:
                theta = self.best_fit_values_median
            else:
                theta = self.best_fit_values_logl
          
        return npars * np.log(sum(len(lc.flux_array) for lc in self.lightcurve_list)) - 2 * self.loglikelihood(theta)

    def AIC(self, theta=None):
      
        npars = len(self.param_list_free)
        
        if theta is None:
            if self.best_fit_median:
                theta = self.best_fit_values_median
            else:
                theta = self.best_fit_values_logl
            
        return 2 * npars - 2 * self.loglikelihood(theta)    

    def red_noise_beta(self, theta=None):

        # Get the RMS of the residuals using the existing function

        if theta is None:
            if self.best_fit_median:
                theta = self.best_fit_values_median
            else:
                theta = self.best_fit_values_logl

        rms_dict = self.rms(theta)

        for ilc,lc in enumerate(self.lightcurve_list):
            rms_val   = rms_dict[f'lc{ilc}']

            time_diff = np.diff(lc.time_array).min() * 24 * 60  # in minutes
            max_points = int(np.round(30 / time_diff))
            bins = np.linspace(lc.time_array[0],
                            lc.time_array[-1],
                            int(len(lc.time_array) / max_points))

            # Rebin residuals
            _, binned_y, _ = pu.rebin(bins,
                                    lc.time_array,
                                    lc.flux_array - lc.calc(),
                                    lc.flux_err, weighted=False)

            max_rms = np.sqrt(np.nanmean(binned_y**2))

            gaussian_white_noise = np.array([1, 1/np.sqrt(max_points)])
            offset = np.max(gaussian_white_noise) / rms_val
            gaussian_white_noise /= offset
            beta_factor = max(max_rms / gaussian_white_noise)

        return beta_factor

    def save_results(self, wb, verbose):

        if self.sampling_method == 'dynesty':
            self.save_dynesty_results(wb, verbose)
        if self.sampling_method == 'emcee':  
            self.save_emcee_results(wb, verbose)

        if self.sampling_method == 'LM':
            # Estimate uncertainties from covariance matrix
            try:
                J = self.sampler_results.jac
                cov = np.linalg.inv(J.T.dot(J))*self.reducedChisq()
                uncertainties = np.sqrt(np.diag(cov))
            except:
                print("Unable to estimate uncertainties from Jacobian")
                uncertainties = np.zeros_like(self.sampler_results.x)
            self.save_LM_results(self.lightcurve_list,self.sampler_results.x,uncertainties,wb,verbose=verbose)
        return
    
    def save_emcee_results(self, wavelength_bin, verbose):

        output_foldername = self.output_foldername
        namelist = self.param_list_free
        ndiscard = self.sampling_arguments['ndiscard']
        npars    = len(namelist)

        self.arr_low2, self.arr_low1, _, self.arr_high1, self.arr_high2 = self.get_confidence_interval(bin_number=wavelength_bin+1,
                                                                                                                    equal_prob=True,
                                                                                                                    verbose=verbose,
                                                                                                                    save_results=True)

        # save plots of chains
        pu.plot_chains(self.sampler_mcmc, self.burn,wavelength_bin,
                       npars, namelist, ndiscard,
                       save_folder=output_foldername)

        # generate and save corner plot
        pu.make_corner_plot(self.samples,
                            bin_number=(wavelength_bin+1),
                            save_fig=True,
                            namelist=self.param_list_free,
                            save_folder=output_foldername)
        #                    parameter_modes=mode)
        
        return


    def save_dynesty_results(self, wb, verbose):

        output_foldername = self.output_foldername

        # Save results in pickle files
        pickle.dump(self.sampler_results, open(output_foldername + '/pickled_objects/' + 'dynesty_result_wb%s.pickle'%(str(wb+1).zfill(2)),'wb'))
        pickle.dump(self.sampler_results.samples, open(output_foldername + '/pickled_objects/' + 'dynesty_samples_wb%s.pickle'%(str(wb+1).zfill(2)),'wb'))
        pickle.dump(self.sampler_results.importance_weights(), open(output_foldername + '/pickled_objects/' + 'dynesty_weights_wb%s.pickle'%(str(wb+1).zfill(2)),'wb'))
        pickle.dump(self.samples, open(output_foldername + '/pickled_objects/' + 'dynesty_equal_weights_samples_wb%s.pickle'%(str(wb+1).zfill(2)),'wb'))
            
        # Save results for confidence interval
        self.arr_low2, self.arr_low1, _, self.arr_high1, self.arr_high2 = self.get_confidence_interval(bin_number=wb+1,
                                                                                                        best_fit_logl=self.best_fit_values_logl,
                                                                                                        verbose=verbose,
                                                                                                        save_results=True)
        
        # maybe move the corner plot somewhere else
        fig, ax = dyplot.cornerplot(self.sampler_results, color='dodgerblue',
                                truth_color='black', show_titles=True,
                                quantiles=None, max_n_ticks=3)
        fig.savefig(output_foldername + '/plots/' + 'dynesty_corner_plot_wb%s.png'%(str(wb+1).zfill(2)))

        fig, ax = dyplot.traceplot(self.sampler_results, truths=np.zeros(self.nDims),
                             truth_color='black', show_titles=True,
                             trace_cmap='viridis', connect=True,
                             connect_highlight=range(5))
        
        fig.savefig(output_foldername + '/plots/' + 'dynesty_trace_plot_wb%s.png'%(str(wb+1).zfill(2)))
        return


    def save_LM_results(self,param_medians,param_uncertainties,bin_number,verbose=True):
        """Function to save the results from an LM fit to a fitted_parameters.dat and LM_statistics.dat tables equivalent to emcee results.

        Inputs:
        fitted_lightcurve - the best fitting, resulting, Lightcurvemodel object
        param_medians - the best fitting parameters
        param_uncertainties - the 1 sigma uncertainties on the parameters
        bin_number - the bin number (correcting for Python indexing, i.e. adding 1)
        verbose - True/False - print the best-fitting results to terminal?

        Returns:
        Nothing but saving the results to fitted_parameters.dat and LM_statistics.txt"""

        ndim = len(param_medians)
        namelist = self.param_list_free
        output_foldername = self.output_foldername

        if bin_number == 0:
            new_tab = open(output_foldername + 'fitted_parameters.txt','w')
        else:
            new_tab = open(output_foldername + 'fitted_parameters.txt','a')

        print('\nSaving fitted parameters to table...\n')

        for i in range(ndim):
            # note, we repeat the uncertainties column twice here even though there is only one uncertainty value, this is so the other functions can better handle this table
            new_tab.write("%s_%d = %f + %f - %f \n"%(namelist[i].replace('$','').replace("\\",''),bin_number+1,param_medians[i],param_uncertainties[i],param_uncertainties[i]))

            if verbose:
                print("%s_%d = %f +/- %f"%(namelist[i].replace('$','').replace("\\",''),bin_number+1,param_medians[i],param_uncertainties[i]))

        new_tab.write('#------------------ \n')
        new_tab.close()

        return

    def write_fit_diagnostics(self,wavelength_bin,emcee_sampler=None,nsteps=None,burn=False):

        output_foldername = self.output_foldername

        if wavelength_bin == 0:
            read_mode = 'a'
        else:
            read_mode = 'w'

        if self.sampling_method == 'emcee':

            if burn:
                diagnostic_tab = open(output_foldername+ '/tables/' + 'burn_statistics.txt',read_mode)
            else:
                diagnostic_tab = open(output_foldername+ '/tables/' + 'prod_statistics.txt',read_mode)

        if self.sampling_method == 'LM':

            diagnostic_tab = open(output_foldername+ '/tables/' + 'LM_statistics.txt',read_mode)
        
        if self.sampling_method == 'dynesty':
            diagnostic_tab = open(output_foldername+ '/tables/' + 'dynesty_statistics.txt',read_mode)
            
        fitted_chi2 = self.chisq()
        fitted_reducedChi2 = self.reducedChisq()
        fitted_rms = self.rms()                     # Returns dictionary
        fitted_BIC = self.BIC()
        fitted_AIC = self.AIC()

        print('\nCalculating statistics for best fit...')
        if self.nlightcurves > 1:
            print("\n" + "="*40)
            print("Joint fit statistics")
            print("="*40)

            print("\nGlobal statistics:")
            print('Total chi2 = %.3f' % fitted_chi2[f'total'])
            print('Total reduced chi2 = %.3f' % fitted_reducedChi2[f'total'])

        print('BIC = %f' % fitted_BIC) # global (joint likelihood)
        print('AIC = %f' % fitted_AIC) # global
    
        print("\nIndividual light curve statistics:")
        for ilc,lc in enumerate(self.lightcurve_list):
            print(f'LC {ilc}')
            print('  Chi2 = %.3f' % fitted_chi2[f'lc{ilc}'])
            print('  Reduced chi2 = %.3f' % fitted_reducedChi2[f'lc{ilc}'])
            print('  Residual RMS (ppm) = %d' % (fitted_rms[f'lc{ilc}']*1e6))
        
        diagnostic_tab.write("\nBin %d \n" % (wavelength_bin))
        if self.nlightcurves > 1:
            diagnostic_tab.write("--- Joint fit statistics ---\n")
            print("\nIndividual light curve statistics:")
            for ilc,lc in enumerate(self.lightcurve_list):
                print(f'LC {ilc}')
                print('  Chi2 = %.3f' % fitted_chi2[f'lc{ilc}'])
                print('  Reduced chi2 = %.3f' % fitted_reducedChi2[f'lc{ilc}'])
                print('  Residual RMS (ppm) = %d' % (fitted_rms[f'lc{ilc}']*1e6))
                
        diagnostic_tab.write("\nBin %d \n" % (wavelength_bin))
        if self.nlightcurves > 1:
            diagnostic_tab.write("--- Joint fit statistics ---\n")

            diagnostic_tab.write("Global statistics:")
            diagnostic_tab.write('Total chi2 = %.3f \n' % fitted_chi2[f'total'])
            diagnostic_tab.write('Total reduced chi2 = %.3f \n' % fitted_reducedChi2[f'total'])

        diagnostic_tab.write('BIC = %f \n' % fitted_BIC) # global (joint likelihood)
        diagnostic_tab.write('AIC = %f \n' % fitted_AIC) # global

        if self.nlightcurves > 1:
            diagnostic_tab.write("Individual light curve statistics:\n")

        for ilc,lc in enumerate(self.lightcurve_list):
            diagnostic_tab.write('Residual RMS (ppm) = %d \n' % (fitted_rms[f'lc{ilc}']*1e6))

        if self.sampling_method == 'dynesty':
            print(self.sampler_results.summary())
            import sys
            o = sys.stdout
            # Redirect stdout to a file
            sys.stdout = diagnostic_tab
            print(self.sampler_results.summary())
            sys.stdout = o

        if emcee_sampler is not None:
            try:
                print('\nAutocorrelation time for each parameter = ',np.round(emcee_sampler.acor).astype(int))
                # Alternatively something like: emcee.autocorr.integrated_time(sampler.chain, low=10, high=None, step=1, c=5, full_output=True,axis=0, fast=False)
                diagnostic_tab.write('\nAutocorrelation time for each parameter = ')
                for ac in np.round(emcee_sampler.acor).astype(int):
                    diagnostic_tab.write('%d '%ac)
                diagnostic_tab.write('\n')

                print('nsamples/median(autocorrelation time) = %d'%np.round(nsteps/np.median(emcee_sampler.acor)))
                diagnostic_tab.write('nsamples/median(autocorrelation time) = %d \n'%(np.round(nsteps/np.median(emcee_sampler.acor))))
            except:
                print("\nAutocorrelation time can't be calculated - chains likely too short")
                diagnostic_tab.write("\nAutocorrelation time can't be calculated - chains likely too short \n")

            print('Acceptance fraction = %f'%(np.mean(emcee_sampler.acceptance_fraction)))

            diagnostic_tab.write('Acceptance fraction = %f \n'%(np.mean(emcee_sampler.acceptance_fraction)))
            diagnostic_tab.write('Total steps = %d \n'%(nsteps))

        diagnostic_tab.write('#------------------ \n')
        diagnostic_tab.close()

    def get_confidence_interval(self, bin_number, best_fit_logl=None, equal_prob=False, save_results=False, verbose=False):

        """ Get confidence interval from weighted CDF
        """

        # Call namelist
        namelist = self.param_list_free
        output_foldername = self.output_foldername

        # Defining the confidence intervals for 1, 2, and 3 sigma
        sig_1 = 0.5 + 0.6826/2.0
        sig_2 = 0.5 + 0.954/2.0

        if equal_prob:
            # Equal probability 
            samples = self.samples
            prob = np.empty(samples)
            prob.fill(1/samples)
        else:
            results = self.sampler_results
            samples = results.samples
            prob = results.importance_weights()

        if best_fit_logl is not None:
            print_best_fit = True
        else:
            print_best_fit = False

        npars = np.shape(samples)[1]

        arr_low2 = np.zeros(shape=(npars))
        arr_low1 = np.zeros(shape=(npars))
        arr_median = np.zeros(shape=(npars))
        arr_high1 = np.zeros(shape=(npars))
        arr_high2 = np.zeros(shape=(npars))

        if save_results:
        
            if bin_number == 1:
                new_tab = open(output_foldername + '/tables/' + 'fitted_parameters_median.txt','w')
                if print_best_fit:
                    new_tab_bf = open(output_foldername + '/tables/' + 'fitted_parameters_max_likelihood.txt','w')
            else:
                new_tab = open(output_foldername + '/tables/' + 'fitted_parameters_median.txt','a')
                if print_best_fit:
                    new_tab_bf = open(output_foldername + '/tables/' + 'fitted_parameters_max_likelihood.txt','a')

        # Loop to get confidence intervals for each parameter
        for i in range(npars):

            # Combine the probability and sample values into a single array for sorting
            arr_ordered = list(zip(prob[:], samples[:, i]))

            # Sort the array based on the sample values
            arr_ordered.sort(key=lambda x: x[1])
            arr_ordered = np.array(arr_ordered)

            # Making the CDF, sum of the prob for each row
            arr_ordered[:,0] = arr_ordered[:,0].cumsum()

            # Interpolate the CDF to find the sample values corresponding to the desired confidence intervals, x=CDF, y=sample values
            arr_ordered_interp = lambda x: np.interp(x, arr_ordered[:,0], arr_ordered[:,1],
                                                        left=arr_ordered[0,1], right=arr_ordered[-1,1])


            # And then you get the sample values corresponding to the desired confidence intervals
            arr_low2[i] = arr_ordered_interp(1-sig_2)
            arr_low1[i] = arr_ordered_interp(1-sig_1)
            arr_median[i] = arr_ordered_interp(0.5)
            arr_high1[i] = arr_ordered_interp(sig_1)
            arr_high2[i] = arr_ordered_interp(sig_2) 

        # Loop to print best fit parameters with median
        if save_results:

            print('\nFitted parameters with median...\n')

            for i in range(npars):
                key = namelist[i].replace('$','').replace("\\",'')
                mid_value = arr_median[i]

                write_parameters_med = f'{key}_{bin_number:d} = {mid_value:.6f} + {arr_high1[i] - mid_value:.6f} - {mid_value - arr_low1[i]:.6f}'
                new_tab.write(write_parameters_med + "\n")

                if verbose:
                    print(write_parameters_med)

            print('\nSaving to table...\n')
            new_tab.write('#------------------ \n')
            new_tab.close()

        # Loop to print best fit parameters with highest log likelihood
            if print_best_fit:

                new_tab_bf.write('# Uncertainties are still calculated from 68% confidence interval!\n')
                new_tab_bf.write('# If the uncertainties are lower than 0, values from the highest likelihood is outside of the confidence interval, which should NOT happen! \n')
                new_tab_bf.write('# If it does, please re-evaluate! \n')

                print('\nFitted parameters with highest log likelihood...\n')

                for i in range(npars):

                    key = namelist[i].replace('$','').replace("\\",'')
                    mid_value = best_fit_logl[i]
                    write_parameters_bf = f'{key}_{bin_number:d} = {mid_value:.6f} + {arr_high1[i] - mid_value:.6f} - {mid_value - arr_low1[i]:.6f}'

                    if (arr_high1[1] - mid_value < 0) or (mid_value - arr_low1[i] < 0):
                        print('Warning! Parameter from highest likelihood is outside the percentile! Re-evaluate!')

                    new_tab_bf.write(write_parameters_bf + "\n")

                    if verbose:
                        print(write_parameters_bf)

                print('\nSaving to table...\n')
                new_tab_bf.write('#------------------ \n')
                new_tab_bf.close()

        self.arr_low2 = arr_low2 
        self.arr_low1 = arr_low1
        self.arr_median = arr_median
        self.arr_high1 = arr_high1
        self.arr_high2 = arr_high2

        return arr_low2, arr_low1, arr_median, arr_high1, arr_high2

    # def get_arrays_for_sigma_plotting(self):

    #     low_2 = copy.deepcopy(self.lightcurve_list)
    #     low_1 = copy.deepcopy(self.lightcurve_list)
    #     high_1 = copy.deepcopy(self.lightcurve_list)
    #     high_2 = copy.deepcopy(self.lightcurve_list)

    #     list_lcs = [low_2, low_1, high_1, high_2]
    #     arrs = [self.arr_low2, self.arr_low1, self.arr_high1, self.arr_high2]

    #     for i, lc_array in enumerate(list_lcs):
    #         for ilc, lc in enumerate(lc_array):
    #             if len(list_lcs)>1:
    #                 lc_theta, _ = self.joint_fitter.map_theta(arrs[i],ilc,self.param_list_free)
    #             else:
    #                 lc_theta = arrs[i]
    #             lc.update_model(lc_theta)

    #     return [self.lightcurve_list, low_1, high_1, low_2, high_2]

    def recover_quartiles_single(self,samples,bin_number,verbose=True,save_result=False,burn=False):
        """
        Function that calculates the 16th, 50th and 84th percentiles from a numpy array / emcee chain and saves these to a table.

        Inputs:
        samples - the samples/chains from emcee
        namelist - the names of the parameters that were fit - needed for printing and saving to file
        bin_number - the number of the wavelength bin we're considering. Necessary for printing and saving to file.
        verbose - True/False: do we want to print the results to screen?
        save_result - True/False: do we want to save the results to a table?
        burn - True/False: is this a burn-in chain? If so, save to burn_parameters.dat, else save to fitted_parameters.dat

        Returns:
        (median, upper bound, lower bound) with shape (nparameters,3)
        """

        namelist = self.param_list_free
        lower = []
        median = []
        upper = []
        mode = []

        ndim = np.shape(samples)[1] # this is equal to the number of params

        length_nl = len(namelist)

        # generate dictionary of how many decimal places we want to round each parameter to before calculating the mode of the rounded distribution
        namelist_decimal_places = {"t0":6,"per":6,"rp":6,"a":2,"inc":2,"ecc":3,"w":2,\
                                "u1":2,"u2":2,"u3":2,"u4":2,"f":4,"s":3,"A":3,"step1":3,"step2":3,"breakpoint":0,"lniL":2,
                                "infl":3}

        # now pad the dictionray with the systematics coefficients which could be a large number of parameters (although much less than the 100 allowed for below)
        for i in range(100):
            namelist_decimal_places["r%s"%(i)] = 2
            namelist_decimal_places["c%s"%(i)] = 6
            namelist_decimal_places["lniL_%s"%(i)] = 2

        if save_result:
            if bin_number == 1:
                if burn:
                    new_tab = open('burn_parameters.txt','w')
                    # new_tab_2 = open('parameter_modes_burn.txt','w')
                else:
                    new_tab = open('fitted_parameters.txt','w')
                    new_tab_2 = open('parameter_modes_prod.txt','w')
            else:
                if burn:
                    new_tab = open('burn_parameters.txt','a')
                    # new_tab_2 = open('parameter_modes_burn.txt','a')
                else:
                    new_tab = open('fitted_parameters.txt','a')
                    new_tab_2 = open('parameter_modes_prod.txt','a')

        for i in range(ndim):
            par = samples[:,i]
            lolim,best,uplim = np.percentile(par,[16,50,84])
            lower.append(lolim)
            median.append(best)
            upper.append(uplim)

            # calculate mode of rounded sample array
            key = namelist[i].replace('$','').replace("\\",'')
            key = key.split("_")[0]
            rounded_par = np.round(par,namelist_decimal_places[key])
            mode_value, mode_count = stats.mode(rounded_par,keepdims=True)
            mode.append(mode_value[0])

            if save_result:
                new_tab.write("%s_%d = %f + %f - %f \n"%(namelist[i].replace('$','').replace("\\",''),bin_number,best,uplim-best,best-lolim))
                if not burn:
                    new_tab_2.write("%s_%d = %f (%d counts = %d%%) \n"%(namelist[i].replace('$','').replace("\\",''),bin_number,mode_value[0],mode_count[0],100*mode_count[0]/len(par)))

            if verbose:
                print("%s_%d = %f + %f - %f"%(namelist[i].replace('$','').replace("\\",''),bin_number,best,uplim-best,best-lolim))
                if not burn:
                    print("%s_%d (mode of posterior) = %f (%d counts = %.2f%%) \n"%(namelist[i].replace('$','').replace("\\",''),bin_number,mode_value[0],mode_count[0],100*mode_count[0]/len(par)))

        if save_result:
            print('\nSaving fitted parameters to table...\n')
            new_tab.write('#------------------ \n')
            new_tab.close()

            if not burn:
                new_tab_2.write('#------------------ \n')
                new_tab_2.close()

        return np.array(median),np.array(upper),np.array(lower),np.array(mode)


    def update_prior_file(self, input_prior_file, wb, ilc, best_fit_median=True, save_folder=None):
        """
        Create a new prior file where values are replaced by best-fit medians and uncertainties.

        Rules:
        - If prior_type == 'N': mean = median, sigma = avg_uncertainty
        - If prior_type == 'U' AND parameter is fixed: keep old bounds unchanged
        - If prior_type == 'U' AND parameter is free BUT "preserve U range": keep bounds unchanged
        - If prior_type == 'U' AND not preserving range:
                lower = median - err
                upper = median + err

        - Special: t0, a, inc, per are switched to "fixed"
        """

        output_foldername = self.output_foldername
        if best_fit_median:
            best_fit_file = "fitted_parameters_median.txt"
        else:
            best_fit_file = "fitted_parameters_max_likelihood.txt"

        # -------------------------------------------
        # Load best-fit parameters
        # -------------------------------------------
        bestfit = {}
        with open(output_foldername + "/tables/" + best_fit_file, "r") as bf:
            for line in bf:
                if "=" not in line:
                    continue

                parts = line.split()
                # Format: rp_1 = 0.092987 + 0.000074 - 0.000074
                name = parts[0].replace("_1", "")
                median = float(parts[2])
                plus = float(parts[4])
                minus = float(parts[6])
                avg_unc = 0.5*(plus + minus)

                bestfit[name] = (median, plus, minus, avg_unc)

        # -------------------------------------------
        # Output file
        # -------------------------------------------
        new_prior_file = input_prior_file.replace(".txt", "_fitted_lc%s_wb%s.txt"%(str(ilc), str(wb+1).zfill(4)))

        if save_folder is not None:
            new_prior_file = f'{save_folder}/{new_prior_file}'

        # Parameters that must be fixed
        force_fix = {"t0", "a", "inc", "per"}

        with open(input_prior_file, "r") as infile, open(new_prior_file, "w") as outfile:

            for line in infile:
                stripped = line.strip()

                # Keep blank lines or comments
                if stripped == "" or stripped.startswith("#"):
                    outfile.write(line)
                    continue

                parts = line.split()
                if len(parts) < 6:
                    outfile.write(line)
                    continue

                pname, fitflag, value, p1, p2, ptype = parts[:6]
                remainder = parts[6:]  # preserve trailing comments or columns

                # Not a parameter present in best-fit file
                if pname not in bestfit:
                    outfile.write(line)
                    continue

                median, plus, minus, avg_unc = bestfit[pname]

                # Determine if this param should become fixed
                is_forced_fixed = pname in force_fix

                # -------------------------------------------
                # If forced fixed
                # -------------------------------------------
                if is_forced_fixed:
                    new_fitflag = "fixed"
                    new_value = f"{median:.8f}"

                    # Preserve original priors because they are irrelevant once fixed
                    new_p1 = p1
                    new_p2 = p2

                # -------------------------------------------
                # If not forced fixed
                # -------------------------------------------
                else:

                    new_fitflag = fitflag  # usually "free"

                    if ptype == "N":
                        # Gaussian prior: mean = median ; sigma = avg unc
                        new_value = f"{median:.8f}"
                        new_p1 = f"{median:.8f}"
                        new_p2 = f"{avg_unc:.8f}"

                    elif ptype == "U":

                        # Preserve range for all uniform priors UNLESS you explicitly tighten
                        preserve_range = True

                        if preserve_range:
                            new_value = f"{median:.8f}"
                            new_p1 = p1
                            new_p2 = p2
                        else:
                            # Tighten the bounds to median ± error
                            new_value = f"{median:.8f}"
                            new_p1 = f"{(median - minus):.8f}"
                            new_p2 = f"{(median + plus):.8f}"

                    else:
                        # Unknown prior → keep original
                        new_value = f"{median:.8f}"
                        new_p1 = p1
                        new_p2 = p2

                # -------------------------------------------
                # Construct output line
                # -------------------------------------------
                newcols = [pname, new_fitflag, new_value, new_p1, new_p2, ptype]
                if remainder:
                    newcols += remainder

                outfile.write("    ".join(newcols) + "\n")

        return new_prior_file