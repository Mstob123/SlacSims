
function BeamOut=YABP_FBPIC2YABP(BeamIn,RoundWeights)

    if ~exist('RoundWeights')
       RoundWeights=false; 
    end
    
    SI_constants;
    BeamOut=struct();
    
    
    
    NormWeight=round(BeamIn.weight/min(BeamIn.weight));
    if RoundWeights
        NormWeight=round(NormWeight/10)*10;
    end
        
    
    % figure(3)
    % histogram(NormWeight)
    qE = 1.6021766*10^-19;
    TotalBeamCharge=sum(BeamIn.weight)*qE;
    Beam_charge=ones(1,sum(NormWeight))*TotalBeamCharge/sum(NormWeight);

    % Beam_Gamma=[];

    % for n=1:length(NormWeight)
    %     MP_weight=NormWeight(n);
    %     
    %     MP_x=BeamIn.x(n);
    %    
    %     MP_z=BeamIn.z(n);
    %     
    %     MP_px=BeamIn.px(n);
    % 
    %     MP_pz=BeamIn.pz(n);
    %     
    % %     MP_gamma=BeamIn.pz(n)/(mE*c0);
    %     
    %     Beam_x=[Beam_x,ones(1,MP_weight)*MP_x];
    %     
    %     Beam_z=[Beam_z,ones(1,MP_weight)*MP_z];
    %     
    %     Beam_px=[Beam_px,ones(1,MP_weight)*MP_px];
    %     
    %     Beam_pz=[Beam_pz,ones(1,MP_weight)*MP_pz];
    % %     Beam_Gamma=[Beam_Gamma,ones(1,MP_weight)*MP_gamma];
    %     n
    % end
    Beam_x=repelem(BeamIn.x,NormWeight);
    Beam_y=repelem(BeamIn.y,NormWeight);
    Beam_z=repelem(BeamIn.z,NormWeight);
    Beam_px=repelem(BeamIn.px,NormWeight);
    Beam_py=repelem(BeamIn.py,NormWeight);
    Beam_pz=repelem(BeamIn.pz,NormWeight);

    Beam_thetaX=Beam_px./Beam_pz;
    Beam_thetaY=Beam_py./Beam_pz;
    
    Beam_Gamma=Beam_pz/(mE*c0);

    %%

    BeamOut.Bunch.x(1,:)=Beam_x;
    BeamOut.Bunch.x(2,:)=Beam_thetaX;
    BeamOut.Bunch.x(3,:)=Beam_y;
    BeamOut.Bunch.x(4,:)=Beam_thetaY;
    BeamOut.Bunch.x(5,:)=Beam_z;
    BeamOut.Bunch.x(6,:)=Beam_Gamma*0.511/1000;


    BeamOut.Bunch.Q=Beam_charge;
    BeamOut.Bunch.stop=zeros(size(Beam_charge));



    


       


        %%

    Witness2_Out=BeamOut;
end
