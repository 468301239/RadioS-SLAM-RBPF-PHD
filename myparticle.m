classdef myparticle
    properties
        state
        weight
        map
        measures
        
        config
        A
        B
    end
    
    methods
        function obj=myparticle(stateinit,weightinit)
            obj.state=stateinit;
            obj.weight=weightinit;
            obj.map=myphd();
            
            obj.config.ax=0.5;
            obj.config.ay=0.5;
            obj.config.vb=0.005;
            config.dt=0.08;
            
            obj.A=eye(5);
            obj.A(1,3)=config.dt;
            obj.A(2,4)=config.dt;
            obj.B=[0.5*config.dt^2, 0,0;
                0, 0.5*config.dt^2,0;
                config.dt, 0,0;
                0, config.dt,0;
                0, 0,config.dt];
            return;
        end
        
        function obj=predict(obj)
            w=[normrnd(0,1,1)*obj.config.ax;normrnd(0,1,1)*obj.config.ay;normrnd(0,1,1)*obj.config.vb];
            obj.state=obj.A*obj.state+obj.B*w;
            return;
        end
        
        function obj=setmeas(obj,measures)
            obj.measures=measures;
            for i=1:length(measures)
                obj.measures{i}.MeasurementParameters.OriginPosition=[obj.state(1),obj.state(2),obj.state(5)];
            end
            return;
        end
        
        function obj=mapiter(obj)
            obj.map=obj.map.prdupd(obj.measures);
            return;
        end
        
        function obj=reweight(obj)
            x=obj.state(1);
            y=obj.state(2);
            b=obj.state(5);
            obj.weight=obj.weight*multifeaturereweight(obj.measures,obj.map,[x,y,b],100);
        end
        
        function obj=reweightlos(obj,losmeas)
            z=losmeas;
            mp = struct(OriginPosition = [obj.state(1),obj.state(2),obj.state(5)]);
            zexp=RngBrgMeasFcnVT([0;0;0],mp);
            zCov=diag([obj.map.config.Rngnoiselos^2,obj.map.config.Brgnoiselos^2]);
            error = z-zexp;
            error(2)=normalizeAngles(error(2));
            md2=error'/zCov*error;
            obj.weight=obj.weight*1/sqrt(det(2*pi*zCov))*exp(-1/2*md2);
        end
        
        function obj2=clone(obj)
            obj2=myparticle(0,0);
            obj2.state=obj.state;
            obj2.weight=obj.weight;
            obj2.measures=obj.measures;
            obj2.config=obj.config;
            obj2.A=obj.A;
            obj2.B=obj.B;
            obj2.map=clone(obj.map);
            return;
        end
    end
end