function [Archives,Time, IGDs, HVs, SPs]=SGEA(probID,N,runmax)
% <algorithm> 
% A Steady-State and Generational Evolutionary Algorithm for
% Dynamic Multi-objective Optimization (SGEA)
%
% INPUT:
%       prob_id   - problem id
%       N         - population size
%       runmax    - maximum number of runs
%
% Output:
%       Archives  - archived population over all changes
%       Metric    - a struct to keep performance metric values
%                   .Time: a 1 by runmax vector for computational time of each run
%                   .IGDs: a runmax by NC matrix to store the IGD values of NC environmental changes
%                   .HVs: a runmax by NC matrix to store the HV values
%                   .SPs: a runmax by NC matrix to store the SP values
%
% .. Author:
%       - Shouyong Jiang        13/07/2022
%       - Shouyong Jiang        01/01/2022    Add RM-MEDA as a production operator to estimate the changing PS and add more test problems
%
% .. Contact:             
%       - Dr. Shouyong Jiang at math4neu@gmail.com.

    %% Parameter setting
    sensor_ratio =0.1; % ratio of sensors to popsize
    eta = 0.1; % ratio of immigrants
    
    operator = @RMMEAD; % choose RMMEAD operator (default) for reproduction, another option is @EAreal.
    %operator = @EAreal;   
    
    %% Start of the Main Loop
    for run=1:runmax
    
        %%%% creat problem instance
        instance=problem(probID);

        dynPar=num2cell(instance.dynPar); % get frequency, severity of change
        [T0, taut, nt]=deal(dynPar{:});
        T=T0+10*nt*taut; %{T0 +10*nt*taut}

        objIdx=instance.varDim+1:instance.varDim+instance.objDim; % intermediate variable to simplify code
        
        Archive=[];
    
        tau=1;
        
        %%%%%%%%%%%%%%%%%% population initialisation
        Qop=init_population(instance,N);
        
        %%%%%%% evaluate initial population
        Pop=instance.func(Qop,0);
        
        
        [Pop, EP]=environmental_selection([],Pop,instance.objDim);
        
        
        changetime=0;
        C0=[];
    
        %%%%%%%%%%%%%%%%%%%%%Optimization
        tStart=tic;
        while tau<=T
    %         disp(['at generation',num2str(t)])

            t=instance.update_time(tau, taut, nt, T0); % update time instant
            
            checkvalue = false;
            responded = false;
            detections=0;
            
            Pop_=Pop(:,1:end-1);
            if isequal(operator,@RMMEAD)
                try % avoid boundary solutions
                    keep=filter_boundary(Pop(:, objIdx));
                    [Model,probability] = LocalPCA(Pop(keep,1:instance.varDim),instance.objDim,5);
                catch
                    [Model,probability] = LocalPCA(Pop(:,1:instance.varDim),instance.objDim,5);
                end
                      
            end

            % % Steady-state evolution
            for i=randperm(N)
                if ~responded && detections <floor(sensor_ratio*N)
                    checkvalue=detect_change(Pop(randi(N,1),1:end-1),instance, t); % change detection 
                    detections=detections +1;
                end
                if ~responded && checkvalue % change detected
                    % disp(['A change occurs at generation ',num2str(t)])
                    
                    changetime=changetime+1;% count the numbe of changes
                    
                    t_prev=instance.update_time(tau-1, taut, nt, T0); % previous time instant
                    [IGD,SP,HV]=Performance(Pop(:,objIdx),probID,t_prev);% calculate performance metrics
                    
                    IGDt(changetime)=IGD;
                    SPt(changetime)=SP;
                    HVt(changetime)=HV;
                    
                    Archive=[Archive;Pop(:,1:end-1)];% save populations for previous changes
                    %plot(Pop(:,end-2),Pop(:,end-1),'.');
                   
                    % respond to changes
                    [Pop, EP, C0]=respond_change(instance, EP, Pop, C0, eta, t);
                    Pop_=[];
                    
                    responded=true;
                end

                % reproduction
                if isequal(operator,@RMMEAD)
                    QWE=RMMEAD(Model,probability, 1, instance.varDim, instance.objDim);
                else
                    QWE=reproduction_EAreal(Pop, EP, instance.varDim, instance.bounds);
                end
                QWE=repair(QWE,instance.bounds); % repair solutions to be within bounds
                child=instance.func(QWE,t); % calculate objective values
                
                [Pop, EP]=update_pop(child, Pop, EP, instance.objDim); % update population and EP after generating a new individual
            end
    
            % generational environmental selection to primarily maintain diversity
            [Pop, EP]=environmental_selection(Pop_, Pop(:,1:end-1), instance.objDim);

            tau=tau+1; % increase generation counter by 1
        end
        tEnd=toc(tStart);
        
        % save results of current run
        Time(run)=tEnd;
        IGDs(run,:)=IGDt;
        SPs(run,:)=SPt;
        HVs(run,:)=HVt;
        Archives{run}=Archive;
    end
end





