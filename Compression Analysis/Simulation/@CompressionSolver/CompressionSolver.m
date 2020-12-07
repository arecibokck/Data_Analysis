classdef CompressionSolver < handle & matlab.mixin.Copyable
    
    properties (Access = public)
        timeStep;
        finalTime;
        timeSpan;
        numberOfAtoms;
        HDTWavelength;
        VDTWavelength;
        potentialType;
        compressionDirection;
        HDTBeamWaist;
        VDTBeamWaist;
        HDTPower;
        VDTPower;
        positionLimits;
        velocityLimits;
        potentialDepthVDT;
        potentialDepthHDT;
        FigurePosition;
        DebugMode;
        DoSave;
    end
    
    properties (SetAccess = private, GetAccess = public)
        InitialDistributionParameters
    end
    
    properties (Dependent, SetAccess = private)
       trappingFrequency; 
       relevantTFForTemperature; 
       latticeWavelength;
       beamWaist;
       potentialDepth;
       potentialDepthInFreq;
       potentialDepthInTemp;
       initialTemperature;
       recoilEnergy;
       initialKineticEnergy
    end
    
    properties (SetAccess=private, GetAccess=public)
        trapFreqVDTLongitudinal;
        trapFreqVDTTransverse;
        trapFreqHDTx
        trapFreqHDTy
        trapFreqHDTz
        physicalConstants;
        simulationResults;
        
        positionGrid;
        velocityGrid;
        
        initialPositions 
        initialVelocities
        initialKineticEnergyInuK
    end
    
    methods
        function cs  = CompressionSolver(varargin)
            
            p = inputParser;
            addParameter(p, 'TimeStep',       10e-06,...
                @(x) assert(isnumeric(x) && isscalar(x) && (x > 0)));  
            addParameter(p, 'FinalTime',       3e-03,...
                @(x) assert(isnumeric(x) && isscalar(x) && (x > 0)));  
            addParameter(p, 'VDTPower',       70e-3,...
                @(x) assert(isnumeric(x) && isscalar(x) && (x > 0)));  
            addParameter(p, 'HDTPower',       12e-3,...
                @(x) assert(isnumeric(x) && isscalar(x) && (x > 0)));  
            addParameter(p, 'PotentialType',   'gaussian',...
                @(x) assert(any(strcmpi(x,{'gaussian', 'harmonic', 'Incorrect entry! Check and try again!'}))));
            addParameter(p, 'CompressionDirection', 'vertical',...
                @(x) assert(any(strcmpi(x,{'vertical', 'horizontal', 'Incorrect entry! Check and try again!'}))));
            addParameter(p, 'PositionDistributionLimits',       [-1*20e-06 1*20e-06],...
                @(x) isnumeric(x));  
            addParameter(p, 'PositionSampling',       0.1e-6,...
                @(x) isnumeric(x));
            addParameter(p, 'VelocityDistributionLimits',       [-1*10e-03 1*10e-03],...
                @(x) isnumeric(x));
            addParameter(p, 'VelocitySampling',       0.05E-3,...
                @(x) isnumeric(x));
            
            
            addParameter(p, 'DebugMode',       false,...
                @islogical);
            addParameter(p, 'SaveData',        false,...
                @islogical);
            
            p.parse(varargin{:});
            
            cs.timeStep           = p.Results.TimeStep;
            cs.finalTime          = p.Results.FinalTime;  
            cs.VDTPower           = p.Results.VDTPower;
            cs.HDTPower           = p.Results.HDTPower;
            cs.HDTBeamWaist       = 25e-06;
            cs.VDTBeamWaist       = 50e-06;
            cs.HDTWavelength      = 865.9e-9;
            cs.VDTWavelength      = 1064e-9;
            cs.potentialType      = p.Results.PotentialType;
            cs.compressionDirection    = p.Results.CompressionDirection;
            posLimitsArray        = p.Results.PositionDistributionLimits;
            cs.positionLimits     = struct('Lower', posLimitsArray(1), 'Upper', posLimitsArray(2));
            velLimitsArray        = p.Results.VelocityDistributionLimits;
            cs.velocityLimits     = struct('Lower', velLimitsArray(1), 'Upper', velLimitsArray(2));
            cs.physicalConstants  = CompressionHelper.PhysicsConstants(); 
            
            set(0,'units','pixels');
            screensize = get(0,'ScreenSize');
            cs.FigurePosition     = [screensize(3)/3.5 screensize(4)/3.5];
            
            cs.DebugMode          = p.Results.DebugMode;
            cs.DoSave             = p.Results.SaveData;
            
            cs.reinitializeCompressionSolver(p.Results.PositionSampling,p.Results.VelocitySampling);
            
        end
        function U   = potential(this, x)
            switch this.potentialType
                case 'harmonic'
                    U = this.potentialDepth + 0.5 * this.physicalConstants.Cs133Mass * this.trappingFrequency^2 .* x.^2;
                case 'gaussian'
                    U = this.potentialDepth .* exp(-2 *((x/this.beamWaist).^2));
            end
        end
        function F   = force(this, x)
            % F = force(x, Trap) calculates the force [N] that the dipole-potential
            % exerts on a CS133-atom at distance x [m] from its center
            switch this.potentialType
                case 'harmonic'
                    F = - this.physicalConstants.Cs133Mass * this.trappingFrequency^2 * x;
                case 'gaussian'
                    F = this.potentialDepth * 4 * x/(this.beamWaist.^2) .* exp(-2 *((x/this.beamWaist).^2));
            end
        end % - Dipole force given a potential 
        function ret = odefcn(this, ~, x)
            % ret = odefcn(~, x, Trap) implements the ODE for the simulation
            ret = zeros(2,1);
            ret(1) = x(2); % Velocity
            ret(2) = this.force(x(1))/this.physicalConstants.Cs133Mass;
        end % - ODE to solve
    end % - lifecycle
    
    methods
        function set.timeStep(this, val)
            assert(val > 1e-06, 'Not time efficient to compute for time steps smaller than 1 microsecond!');
            this.timeStep = val;
        end
        function ret = get.timeStep(this)
            ret = this.timeStep;
        end
        function set.finalTime(this, val)
            assert(val <= 5e-03, 'Not time efficient to compute for time spans longer than 5 milliseconds!');
            this.finalTime = val;
        end
        function ret = get.finalTime(this)
            ret = this.finalTime;
        end
        function set.timeSpan(this, val)
            this.timeSpan = val;
        end
        function ret = get.timeSpan(this)
            ret = this.timeSpan;
        end
        function set.numberOfAtoms(this, val)
            assert(val <= 10000, 'Not time efficient to compute for atom numbers larger than 10,000!');
            this.numberOfAtoms = val;
        end
        function ret = get.numberOfAtoms(this)
            ret = this.numberOfAtoms;
        end
        function set.HDTWavelength(this, val)
            this.HDTWavelength = val;
        end
        function ret = get.HDTWavelength(this)
            ret = this.HDTWavelength;
        end
        function set.VDTWavelength(this, val)
            this.VDTWavelength = val;
        end
        function ret = get.VDTWavelength(this)
            ret = this.VDTWavelength;
        end
        function set.potentialType(this,type)
            %assert(any(strcmp(SimulationParameters.PotentialType, {'gaussian', 'harmonic'})),'Potential type has to be either Gaussian or Harmonic!');
            this.potentialType = validatestring(type,{'gaussian','harmonic'});
        end
        function ret = get.potentialType(this)
            ret = this.potentialType;
        end
        function set.compressionDirection(this,type)
            this.compressionDirection = validatestring(type,{'vertical','horizontal'});
        end
        function ret = get.compressionDirection(this)
            ret = this.compressionDirection;
        end
        function set.HDTBeamWaist(this, val)
            this.HDTBeamWaist = val;
        end
        function ret = get.HDTBeamWaist(this)
            ret = this.HDTBeamWaist;
        end
        function set.VDTBeamWaist(this, val)
            this.VDTBeamWaist = val;
        end
        function ret = get.VDTBeamWaist(this)
            ret = this.VDTBeamWaist;
        end
        function set.HDTPower(this,pow)
            this.HDTPower = pow;
        end
        function ret = get.HDTPower(this)
            ret = this.HDTPower;
        end
        function set.VDTPower(this,pow)
            this.VDTPower = pow;
        end
        function ret = get.VDTPower(this)
            ret = this.VDTPower;
        end
        function set.positionLimits(this, val)
            this.positionLimits = val;
        end
        function ret = get.positionLimits(this)
            ret = this.positionLimits;
        end
        function set.velocityLimits(this, val)
            this.velocityLimits = val;
        end
        function ret = get.velocityLimits(this)
            ret = this.velocityLimits;
        end
        function set.potentialDepthVDT(this, val)
            this.potentialDepthVDT = val;
        end
        function ret = get.potentialDepthVDT(this)
            ret = this.potentialDepthVDT;
        end
        function set.potentialDepthHDT(this, val)
            this.potentialDepthHDT = val;
        end
        function ret = get.potentialDepthHDT(this)
            ret = this.potentialDepthHDT;
        end
        function set.FigurePosition(this, val)
            this.FigurePosition = val;
        end
        function ret = get.FigurePosition(this)
            ret = this.FigurePosition;
        end
        function set.DebugMode(this, val)
            this.DebugMode = val;
        end
        function ret = get.DebugMode(this)
            ret = this.DebugMode;
        end
        function set.DoSave(this, val)
            this.DebugMode = val;
        end
        function ret = get.DoSave(this)
            ret = this.DoSave;
        end
    end % - setters and getters
    
    methods
        function ret       = get.trappingFrequency(this)
            switch this.compressionDirection
                case 'vertical'
                    ret = this.trapFreqHDTz;    
                case 'horizontal'
                    ret = this.trapFreqVDTTransverse;    
            end
        end
        function ret       = get.relevantTFForTemperature(this)
            switch this.compressionDirection
                case 'vertical'
                    ret = this.trapFreqVDTLongitudinal;
                case 'horizontal'
                    ret = this.trapFreqHDTx;
            end
        end
        function ret       = get.potentialDepth(this)
            switch this.compressionDirection
                case 'vertical'
                    ret = this.potentialDepthHDT;    
                case 'horizontal'
                    ret = this.potentialDepthVDT;    
            end
        end
        function ret       = get.latticeWavelength(this)
            switch this.compressionDirection
                case 'vertical'
                    ret = this.VDTWavelength;    
                case 'horizontal'
                    ret = this.HDTWavelength;    
            end
        end
        function ret       = get.beamWaist(this)
            switch this.compressionDirection
                case 'vertical'
                    ret = this.HDTBeamWaist;    
                case 'horizontal'
                    ret = this.VDTBeamWaist;    
            end
        end
        function ret = get.potentialDepthInFreq(this)
            % in uK
            ret = abs((this.potentialDepth/this.physicalConstants.PlanckConstant)) * 1e-3;     
        end
        function ret = get.potentialDepthInTemp(this)
            % in kHz
            ret = abs((this.potentialDepth/this.physicalConstants.BoltzmannConstant)) * 1e6;
        end
        function ret       = get.recoilEnergy(this)
           ret = (this.physicalConstants.PlanckConstantReduced * (2*pi/this.latticeWavelength))^2 / (2*this.physicalConstants.Cs133Mass); 
        end
        function ret = get.initialKineticEnergy(this)
            ret = 0.5 * this.physicalConstants.Cs133Mass * this.initialVelocities.^2;
        end
        function ret = get.initialKineticEnergyInuK(this)
            ret = this.initialKineticEnergy/this.physicalConstants.BoltzmannConstant*1E6;
        end
    end % - getters for dependent properties
    
    methods
        function plotAllPotentials(this)
           f_h = DQSIMhelper.getFigureByTag('potentials');
           a_h = get(f_h, 'CurrentAxes');
           if ~isempty(get(a_h, 'Children')) 
               clf(f_h);
           end
           this.plotPotential_private('gaussian'); 
           a_h = get(f_h, 'CurrentAxes');
           hold(a_h, 'on');
           this.plotPotential_private('harmonic'); 
        end
        function plotInitialConditions(this,varargin)
            % 
            
            % - inputHandling
            p = inputParser;
            p.addParameter('NumberOfBins',20,@(x) isnumeric(x) && isscalar(x) && mod(x,1)==0 && x>0 );
            p.parse(varargin{:});
            % - get Figure
            f_h = DQSIMhelper.getFigureByTag('InitialDistributionOfCloud');
            a_h = get(f_h, 'CurrentAxes');
            if ~isempty(get(a_h, 'Children')) 
            	clf(f_h);
            end
            f_h.Name = ['Initial conditions ' this.compressionDirection ' compression'];
            f_h.Units = 'pixel';
            f_h.Position = [this.FigurePosition 1200 600];
            set(groot,'CurrentFigure',f_h)
            ax1 = subplot(1,3,[1,1.8]);
            yyaxis left
            U = (this.potential(this.positionGrid)./this.physicalConstants.BoltzmannConstant) .*1e6;
            p_1 = plot(ax1,this.positionGrid*1e6, U, 'Color', [0, 0.4470, 0.7410], 'LineStyle', '-', 'LineWidth', 2); % plot the potential
            
            ylabel('Potential (uK)','FontSize', 14)
            xlim([min(this.positionGrid) max(this.positionGrid)]*1e6) 
            
            yyaxis right
            h1 = histogram(ax1,this.initialPositions*1e6,p.Results.NumberOfBins,'DisplayName','Position distribution','LineStyle', 'none');
            h1.Normalization = 'probability';
            ylabel('P(x)','Visible','off')
            xlabel([this.CapitalizeString(this.compressionDirection) ' position (um)'],'FontSize', 14)
            
            
            ax2 = subplot(1,3,3);
            h2 = histogram(this.initialKineticEnergyInuK,'Orientation','horizontal','DisplayName',['P_{GS} = ' num2str(this.InitialDistributionParameters.GroundStatePopulation)]); 
            h2.Normalization = 'probability';
            xlabel('Probability')
            ylabel('E_{kin} at release  [\mu K]')
            YLIM = get(ax2,'ylim');
            ylim([YLIM(1)  max(YLIM(2),2)])
            legend
            
        end
    end % - plotting
    
    methods(Access = private)
        plotPotential_private(this, varargin);
        function str = CapitalizeString(this,str)
            str(1) = upper(str(1));
        end
        
    end
    
    methods(Access = protected)
        function cp = copyElement(this)
            % Shallow copy object
            cp = copyElement@matlab.mixin.Copyable(this);
            
            % Forces the setter to redefine the function handles to the new copied object
            % cp.potentialType = this.potentialType;
            
            pl = properties(this);
            for k = 1:length(pl)
                sc = superclasses(this.(pl{k}));
                if any(contains(sc,{'matlab.mixin.Copyable'}))
                    cp.(pl{k}) = this.(pl{k}).copy();
                end
            end
        end
    end
    
end
