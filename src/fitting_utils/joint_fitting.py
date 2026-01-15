import numpy as np 
import emcee 
import re

from fitting_utils import sampling as s

class JointFitter(object):
    def __init__(self,lightcurve_list,verbose=False):


        """
        
        Inputs:
        lightcurve_list    - list of LightcurveModel objects
        sampling_list      - list of Sampling objects
        verbose            - bool, print statements
        
        
        Constructs the global parameter vector.
        
        """

        self.lightcurve_list = lightcurve_list
        self.nlightcurves    = len(lightcurve_list)

        # Collect all joint parameters (should be same across lightcurves)
        joint_param_names = set()
        for lc in self.lightcurve_list:
            joint_param_names.update(lc.param_list_joint)

        self.joint_param_names = list(joint_param_names)
        
        # Build a global parameter list: [joint_params, lc0_params, lc1_params, etc...]
        self.param_list_global = list(joint_param_names).copy()
        self.param_dict_global = {}
        self.prior_dict_global = {}

        # Add the joint parameters to dictionaries
        for param_name in joint_param_names:
            for ilc,lc in enumerate(self.lightcurve_list):
                if param_name not in lc.param_list_joint:
                    raise ValueError(f'Parameter {param_name} not defined as joint parameter in all lightcurves in this list. Check prior.txt files.')
                else:
                    if param_name not in self.param_dict_global:
                        # Add to the dictionary
                        self.param_dict_global[param_name]          = lc.param_dict[param_name]
                        self.prior_dict_global[param_name+'_1']     = lc.prior_dict[param_name+'_1']
                        self.prior_dict_global[param_name+'_2']     = lc.prior_dict[param_name+'_2']
                        self.prior_dict_global[param_name+'_prior'] = lc.prior_dict[param_name+'_prior']
                    #else: need to add a sanity check that the definitions are the same across datasets (in prior.txt)
                
        # Add the individual free parameters
        for ilc,lc in enumerate(self.lightcurve_list):
            for param_name in lc.param_list_free:
                if param_name not in joint_param_names:
                    self.param_list_global.append(f"{param_name}_lc{ilc}") # Update the list of free parameter names e.g. rp_lc0, rp_lc1
                    self.param_dict_global[f"{param_name}_lc{ilc}"]       = lc.param_dict[param_name] # Update the global parameter dictionary (with values)
                    self.prior_dict_global[f"{param_name}_lc{ilc}_1"]     = lc.prior_dict[param_name+'_1']
                    self.prior_dict_global[f"{param_name}_lc{ilc}_2"]     = lc.prior_dict[param_name+'_2']
                    self.prior_dict_global[f"{param_name}_lc{ilc}_prior"] = lc.prior_dict[param_name+'_prior']

        self.npars  = len(self.param_list_global)
        self.njoint = len(joint_param_names)

        if verbose:
            print(f"Jointly fitting {self.njoint} parameter across {self.nlightcurves} lightcurves.")
            print(f"Total free parameters: {self.npars}.")
            print(f"Joint parameters: {joint_param_names}.")
            print(f"Full global parameter list: {self.param_list_global}.")

            

    def map_theta(self,global_theta,lc_idx,global_param_list):

        """
        Take the global parameters and assign to each light curve.

        Inputs:
        global_theta      - global parameter vector [joint_params, lc0_params, lc1_params, etc...]
        lc_idx            - int, light curve index (in lightcurve_list)
        global_param_list - list, parameter names e.g. [t0, c1_lc0, c1_lc1]

        Returns:
        lc_theta      - lightcurve-specific parameter vector for input to Sampling.loglikelihood_emcee()

        """
        lc = self.lightcurve_list[lc_idx]

        lc_theta      = [] 
        lc_param_list = []

        for i, param_name in enumerate(global_param_list):
            # first the joint parameters 
            if param_name in self.joint_param_names:
                lc_theta.append(global_theta[i])
                lc_param_list.append(param_name)
            elif f'_lc{lc_idx}' in param_name:
                lc_theta.append(global_theta[i])
                lc_param_list.append(re.sub(f'_lc{lc_idx}','',param_name)) # remove lc label
        return lc_theta,lc_param_list