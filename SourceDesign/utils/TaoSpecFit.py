from math import exp
import copy

class MCShape:
    """
    MC Shape for a radioactive source
    """
    def __init__(self,source_name,energy,energy_leak_template):
        """
        args:
            - SourceName : Name of Source
            - GammaEnergy: Energy
            - energy_leak_template : energy leak template
        """
        self.energy = copy.deepcopy(energy)
        self.source_name = copy.deepcopy(source_name)
        self.e_leak_template = copy.deepcopy(energy_leak_template)

    def Gaus(self,x,par):
        """
        gaus
        """
        gaus = par[0]*exp(-0.5*pow((x[0]-par[1])/(par[1]*par[2]*0.01),2))
        return gaus
    
    def ELeakSpec(self,x,par):
        """
        Energy Leak spectrum
        args:
            x[0]  : energy
            par[0]: amplitude
            par[1]: peak value of full energy
        """
        escale = par[1]/self.energy
        return par[0]*self.e_leak_template.Eval(x[0]*1.0/escale)
    
    def __call__(self,x,par):
        gaus_par = [par[0],par[1],par[2]]
        gaus = self.Gaus(x,gaus_par)

        eleak_par = [par[0]*par[3],par[1]]
        eleak_value = self.ELeakSpec(x,eleak_par)

        return gaus + eleak_value

class TotalMCShape:
    """
    MC Shape for combine gamma source spectrum
    """
    def __init__(self,source_names,energies,energy_leak_templates):
        """
        args:
            - source_names : Names of Sources
            - energies: Energies
            - energy_leak_templates : energy leak template
        """
        self.energies = copy.deepcopy(energies)
        self.source_names = copy.deepcopy(source_names)
        self.e_leak_templates = copy.deepcopy(energy_leak_templates)
        self.mc_shapes = []
        for i in range(len(source_names)):
            mc_shape = MCShape(self.source_names[i],self.energies[i],self.e_leak_templates[i])
            self.mc_shapes.append(mc_shape)
    
    def __call__(self,x,pars):
        values = 0
        for i in range(len(self.source_names)):
            p_par = [pars[j] for j in range(i*4,i*4+4)]
            values += self.mc_shapes[i](x,p_par)
        return values

class LinearBkgShape:
    """
    Model for radioactive from n-H which have complicated shape.
    """
    def __init__(self,source_name,energy):
        """
        args:
            - SourceName : Name of Source
            - GammaEnergy: Energy
            - energy_leak_template : energy leak template
        """
        self.energy = copy.deepcopy(energy)
        self.source_name = copy.deepcopy(source_name)

    def Gaus(self,x,par):
        """
        gaus
        """
        gaus = par[0]*exp(-0.5*pow((x[0]-par[1])/(par[1]*par[2]*0.01),2))
        return gaus
    
    def LinearBkg(self,x,par):
        """
        Energy Leak spectrum
        args:
            x[0]  : energy
            par[0]: slop of linear bkg
            par[1]: bias of linear bkg
        """
        return par[0] * x[0] + par[1]
    
    def __call__(self,x,par):
        gaus_par = [par[0],par[1],par[2]]
        gaus = self.Gaus(x,gaus_par)

        eleak_par = [par[0]*par[3],par[0]*par[4]]
        eleak_value = self.LinearBkg(x,eleak_par)

        return gaus + eleak_value