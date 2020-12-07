classdef PhysicsConstants < handle
    properties (Constant)
        % CODATA
        PlanckConstant=6.62607015E-34;
        PlanckConstantReduced=6.62607015E-34/(2*pi);
        FineStructureConstant=7.2973525698E-3;
        ElectronMass=9.10938291E-31;
        GravitationalConstant=6.67384E-11;
        ProtonMass=1.672621777E-27;
        AtomicMassUnit=1.66053878283E-27;
        BohrRadius=0.52917721092E-10;
        BohrMagneton=927.400968E-26;
        BoltzmannConstant=1.380649E-23;
        StandardGravityAcceleration=9.80665;
        SpeedOfLight=299792458;
        StefanBoltzmannConstant=5.670373E-8;
        ElectronCharge=1.602176634E-19;
        VacuumPermeability=1.25663706212E-6;
        DielectricConstant=8.8541878128E-12;
        ElectronGyromagneticFactor=-2.00231930436153;
        AvogadroConstant=6.02214076E23;
        
        % Cs specific constants
        Cs133Mass=132.905429*1.66053878283E-27;
        Cs133NuclearSpin=7/2;
        CsD2Gamma = 2*pi*5.23413E6;
        CsD1Gamma = 2*pi*4.57512E6;
        CsD2lambda = 852.3472758227E-9;
        CsD1lambda = 894.5929598610E-9;
        CsHyperfineSplittingFreq = 9192631770; 
    end
    
    methods
        function pc = PhysicsConstants()
        end
    end
    
end
