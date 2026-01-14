import numpy as np 
import emcee 

from fitting_utils import sampling as s

class JointFitter(object):
    def __init__(self,lightcurve_list):


        """
        
        Inputs:
        lightcurve_list    - list of LightcurveModel objects
        sampling_list      - list of Sampling objects
        
        
        Constructs the global parameter vector.
        
        """

        self.lightcurve_list = lightcurve_list 
        #self.sampling_list   = sampling_list 
        #if len(sampling_list)!=len(lightcurve_list):
        #    raise ValueError("Sampling list is not the same length as lightcurve list.")
        self.nlightcurves    = len(lightcurve_list)

        # Collect all joint parameters (should be same across lightcurves)
        joint_param_names = set()
        for lc in self.lightcurve_list:
            joint_param_names.update(lc.param_list_joint)


        # Build a global parameter list: [joint_params, lc0_params, lc1_params, etc...]
        self.param_list_global = self.joint_param_names.copy()
        self.param_dict_global = {}
        self.prior_dict_global = {}

        # Add the joint parameters to dictionaries
        for param_name in joint_param_names:
            for ilc,lc in enumerate(self.lightcurve_list):
                if param_name not in lc.param_list_joint:
                    raise ValueError(f'Parameter {param_name} not defined as joint parameter in all lightcurves in this list. Check prior.txt files.')
                else:
                    if param_name not in self.param_dict_global.items():
                        # Add to the dictionary
                        self.param_dict_global[param_name]          = lc.param_dict[param_name]
                        self.prior_dict_global[param_name+'_1']     = lc.prior_dict[param_name+'_1']
                        self.prior_dict_global[param_name+'_2']     = lc.prior_dict[param_name+'_2']
                        self.prior_dict_global[param_name+'_prior'] = lc.prior_dict[param_name+'_prior']
                    #else: need to add a sanity check that the definitions are the same across datasets (in prior.txt)

        # Add the individual free parameters
        for ilc,lc in enumerate(self.lightcurve_list):
            for param_name in lc.param_list_free:
                self.param_list_global.append(f"{param_name}_lc{ilc}") # Update the list of free parameter names e.g. rp_lc0, rp_lc1
                self.param_dict_global[f"{param_name}_lc{ilc}"]       = lc.param_dict[param_name] # Update the global parameter dictionary (with values)
                self.prior_dict_global[f"{param_name}_lc{ilc}_1"]     = lc.prior_dict[param_name+'_1']
                self.prior_dict_global[f"{param_name}_lc{ilc}_2"]     = lc.prior_dict[param_name+'_2']
                self.prior_dict_global[f"{param_name}_lc{ilc}_prior"] = lc.prior_dict[param_name+'_prior']

        self.npars  = len(self.param_list_global)
        self.njoint = len(self.joint_param_names)

        print(f"Jointly fitting {self.njoint} across {self.nlightcurves} lightcurves.")
        print(f"Total free parameters: {self.npars}.")
        print(f"Joint parameters: {self.joint_param_names}.")

        

    def map_theta(self,theta):
        s

    def joint_logprobability(self,theta):

        """
        Inputs:
        theta         - global theta vector (joint_params, lc0_params, lc1_params, etc..)
        
        Returns:
        logprob_total - joint log probability (float)

        """
        log_prob_sum = 0.0

        for ilc,sampling_obj in enumerate(self.sampling_objects_list):
            # Extract parameters for this lightcurve 
            lc_theta = self.map_theta(theta)
            # And compute logprob
            log_prob = sampling_obj.logprobability_emcee(lc_theta)

            if not np.isfinite(log_prob):
                return -np.inf 
            
            log_prob_sum += log_prob
        
        return log_prob_sum