classdef timereporter
    properties
        TimeStart
        Message
        Verbose
        Valid = false;
        CurrentProgress = 0;
        MaximumProgress = -1;
        ProgressPrinted = false;
        LastProgress = 0;
        ProgressStep = 1;
        ProgressDone = false;
    end

    methods
        function obj = timereporter(verbose, timeStart, progressStep)
            if(nargin < 2); timeStart = tic(); end
            obj.TimeStart = timeStart;
            obj.Verbose = verbose;
            if(nargin >= 3); obj.ProgressStep = progressStep; end
        end
        function [obj] = printRunning(obj, msg, verboseRequirement, nJobs, progressStep)
            if(verboseRequirement <= obj.Verbose)
                disp(['[Running] ', msg, '...']);
                obj.Message = msg;
                obj.Valid = true;
                obj.CurrentProgress = 0;
                obj.LastProgress = 0;
                obj.ProgressDone = false;
                if(nargin >= 4); obj.MaximumProgress = nJobs; end
                if(nargin >= 5); obj.ProgressStep = progressStep; end
            end
        end
        
        function [obj] = printProgress(obj, increment)
            if(~obj.Valid); return; end
            if(nargin < 2); increment = 1; end
            if(obj.MaximumProgress == -1)
                error(['The printProgress function cannot be used ', ...
                'without specifying the number of jobs.']);
            end
            obj.CurrentProgress = obj.CurrentProgress + increment;
            progressChange = obj.CurrentProgress - obj.LastProgress;
            if(progressChange < obj.ProgressStep); return; end
            obj.LastProgress = obj.CurrentProgress;
            if(obj.ProgressPrinted)
                fprintf(repmat('\b', 1, 11));
            else
                fprintf(repmat('\b', 1, 4));
            end
            progress = 100 * obj.CurrentProgress / obj.MaximumProgress;
            fprintf('%11s', [' (', sprintf('%4.1f',progress), '%)...']);
            obj.ProgressPrinted = true;
            if(progress >= 100)
                obj.ProgressDone = true;
                fprintf('\n');
            end
        end
        
        function [] = printDone(obj, msg, verboseRequirement)
            if(nargin >= 3 && verboseRequirement <= obj.Verbose)
                obj.Message = msg;
                obj.Valid = true;
            end
            if(obj.Valid)
                if(obj.ProgressPrinted && ~obj.ProgressDone)
                   fprintf('\n');
                end
                disp(['[Done] ', obj.Message, '.', ' TimePassed : ', ...
                    num2str(toc(obj.TimeStart)), ' s.']);
            end
        end
    end

end

