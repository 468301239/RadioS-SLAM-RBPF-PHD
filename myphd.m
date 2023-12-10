classdef myphd
    properties
        phd
        phd_old
        config
        reservedMeasurementsnext
        reservedlistsnext
        sumweight
        sumweight_old
    end
    
    methods
        % 构造函数
        function obj = myphd()
            % landmark motion model
            obj.config.stateTransitionFcn = @constpos;
            obj.config.stateTransitionJacobianFcn=@RngBrgMeasFcnVTjac;
            obj.config.hasAdditiveProcessNoise = true;
            obj.config.Q = 1e-17;

            % landmark measurement model
            obj.config.measurementFcn = @RngBrgMeasFcnVT;
            obj.config.measurementJacobianFcn = @RngBrgMeasFcnjac;
            obj.config.hasAdditiveMeasurementNoise = true;
            obj.config.Rngnoise=0.3;
            obj.config.Brgnoise=4*pi/180;
            obj.config.Rngnoiselos=0.5;
            obj.config.Brgnoiselos=2*pi/180;

            % birth strategy: sigma 2:0.9111, sigma 3: 0.9946
            obj.config.Rngthreshold=2.5;
            obj.config.Brgthreshold=40*pi/180;
            obj.config.wb=0.01;
            
            % detection FoV!!! Set it the same with data gen
            obj.config.fov=30;
            obj.config.fovscale=0.95;
            
            % detection, survival and clutter intensity
            obj.config.Pd=0.9;
            obj.config.Ps=1;
            obj.config.Kc=0.02/(pi*obj.config.fov^2);

            % phd start
            obj.config.m0=zeros(3,0);
            obj.config.P0=zeros(3,3,0);
            obj.config.w0=zeros(1,0);       
            
            obj.config.mergethreshold=3;
            obj.config.prunethreshold=0.0001;
            obj.config.prunethresholdoutfov=0.08;
            obj.config.extractthreshold=0.07;    
            
            obj.phd = gmphd(obj.config.m0,obj.config.P0,...
                    'Weights',obj.config.w0,...
                    'StateTransitionFcn',obj.config.stateTransitionFcn,...
                    'HasAdditiveProcessNoise',obj.config.hasAdditiveProcessNoise,...
                    'ProcessNoise',obj.config.Q,...
                    'MeasurementFcn',obj.config.measurementFcn,...
                    'HasAdditiveMeasurementNoise',obj.config.hasAdditiveMeasurementNoise);
            
            obj.phd_old = gmphd(obj.config.m0,obj.config.P0,...
                    'Weights',obj.config.w0,...
                    'StateTransitionFcn',obj.config.stateTransitionFcn,...
                    'HasAdditiveProcessNoise',obj.config.hasAdditiveProcessNoise,...
                    'ProcessNoise',obj.config.Q,...
                    'MeasurementFcn',obj.config.measurementFcn,...
                    'HasAdditiveMeasurementNoise',obj.config.hasAdditiveMeasurementNoise);
        
            obj.reservedMeasurementsnext=cell(0,1);
            obj.reservedlistsnext=zeros(1,0);
        end
        
        % 方法        
        function obj=prdupd(obj,measures,statexyb)
            predict(obj.phd,0.08);
            reservedMeasurements=obj.reservedMeasurementsnext;
            reservedlists=obj.reservedlistsnext;
            if ~isempty(reservedlists)
                birthPHD=getbirthVTPHD(reservedMeasurements,obj.config.wb,diag([obj.config.Rngnoise^2,obj.config.Brgnoise^2]));
                append(obj.phd, birthPHD);
            end

            if obj.phd.NumComponents==0
                birthPHD=getbirthVTPHD(measures,obj.config.wb,diag([obj.config.Rngnoise^2,obj.config.Brgnoise^2]));
                append(obj.phd, birthPHD);
                return;
            end
            
            mp = struct(OriginPosition = [statexyb(1),statexyb(2),statexyb(3)]);
            allmeas=obj.phd.MeasurementFcn(obj.phd.States,mp);
            if(~isempty(measures))
                obj.reservedMeasurementsnext=cell(0,1);
                obj.reservedlistsnext=zeros(1,0);
                for j=1:length(measures)
                    measdeviation=(allmeas-measures{j}.Measurement);
                    measdeviation(2,:)=normalizeAngles(measdeviation(2,:));
                    measabsdeviation=abs(measdeviation);
                    flagmatrix=[measabsdeviation(1,:)>obj.config.Rngthreshold;measabsdeviation(2,:)>obj.config.Brgthreshold];
                    if  min(sum(flagmatrix,1))>0
                        obj.reservedMeasurementsnext=[obj.reservedMeasurementsnext;measures(j)];
                        obj.reservedlistsnext=[j,obj.reservedlistsnext];
                    end
                end

                for j=1:length(obj.reservedlistsnext)
                    measures(obj.reservedlistsnext(j))=[];
                end

                obj.sumweight_old=sum(obj.phd.Weights);
                obj.phd_old=clone(obj.phd);
            end
            
            % exceeding the fov
            fovindex_expected=allmeas(1,:)<(obj.config.fov*obj.config.fovscale);
            outsidefov_expected=fovindex_expected<1;
            number_infov=sum(fovindex_expected);
            number_outsidefov=sum(outsidefov_expected);
            
            % divide the phd into inside fov and outside fov
            m0phd_infov=zeros(3,number_infov);
            P0phd_infov=zeros(3,3,number_infov);
            weights_infov=zeros(1,number_infov);
            m0phd_outsidefov=zeros(3,number_outsidefov);
            P0phd_outsidefov=zeros(3,3,number_outsidefov);
            weights_outsidefov=zeros(1,number_outsidefov);
            
            countinfov=0;
            countoutfov=0;
            for i=1:length(fovindex_expected)
                if(fovindex_expected(i)==1)                    
                    countinfov=countinfov+1;
                    m0phd_infov(:,countinfov)=obj.phd.States(:,i);
                    P0phd_infov(:,:,countinfov)=obj.phd.StateCovariances(:,:,i);
                    weights_infov(countinfov)=obj.phd.Weights(i);
                else
                    countoutfov=countoutfov+1;
                    m0phd_outsidefov(:,countoutfov)=obj.phd.States(:,i);
                    P0phd_outsidefov(:,:,countoutfov)=obj.phd.StateCovariances(:,:,i);
                    weights_outsidefov(countoutfov)=obj.phd.Weights(i);
                end
            end
            
            phd_infov = gmphd(m0phd_infov,P0phd_infov,...
                    'Weights',weights_infov,...
                    'StateTransitionFcn',obj.config.stateTransitionFcn,...
                    'HasAdditiveProcessNoise',obj.config.hasAdditiveProcessNoise,...
                    'ProcessNoise',obj.config.Q,...
                    'MeasurementFcn',obj.config.measurementFcn,...
                    'HasAdditiveMeasurementNoise',obj.config.hasAdditiveMeasurementNoise);
                
             phd_outsidefov = gmphd(m0phd_outsidefov,P0phd_outsidefov,...
                    'Weights',weights_outsidefov,...
                    'StateTransitionFcn',obj.config.stateTransitionFcn,...
                    'HasAdditiveProcessNoise',obj.config.hasAdditiveProcessNoise,...
                    'ProcessNoise',obj.config.Q,...
                    'MeasurementFcn',obj.config.measurementFcn,...
                    'HasAdditiveMeasurementNoise',obj.config.hasAdditiveMeasurementNoise);
                
            prune(phd_outsidefov, phd_outsidefov.Weights < obj.config.prunethresholdoutfov);  
          
            
            % start filtering process
            obj.phd=phd_infov;
            scale(obj.phd, obj.config.Ps);
            phdUndetected = clone(obj.phd);
            scale(phdUndetected, 1 - obj.config.Pd);
            if ~isempty(measures)
                obj.phd.Detections = measures;
                detectionGroups = logical(eye(length(measures))); 
                logqij = likelihood(obj.phd, detectionGroups);
                Lij = calculateScaling(logqij, obj.phd.Weights, obj.config.Pd, obj.config.Kc);
                correct(obj.phd, detectionGroups, Lij);
                append(obj.phd, phdUndetected);
            else
                obj.phd=phdUndetected;
            end    
            append(obj.phd, phd_outsidefov);
            
            obj.sumweight=sum(obj.phd.Weights);
            
            merge(obj.phd,obj.config.mergethreshold);
            prune(obj.phd, obj.phd.Weights < obj.config.prunethreshold);  
            return;
        end
        
        function obj2=clone(obj)
            obj2=myphd();
            obj2.phd=clone(obj.phd);
            obj2.phd_old=clone(obj.phd_old);
            obj2.config=obj.config;
            obj2.reservedMeasurementsnext=obj.reservedMeasurementsnext;
            obj2.reservedlistsnext=obj.reservedlistsnext;
            obj2.sumweight=obj.sumweight;
            obj2.sumweight_old=obj.sumweight_old;
            return;
        end
        
    end
end