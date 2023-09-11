function Elec=FBPIC_ReadElectron(folder,file,species)
    Elec=struct();
%     folder='/Users/aknetschmac/Python/FBPIC/PWFA3/diags/hdf5/';

%     file='data00001250.h5';
    DatabaseInfo=h5info([folder,file]);

    DatabaseName=DatabaseInfo.Groups.Groups.Name;
    SlashPositions=strfind(DatabaseName,'/');
    TimeStep=str2num(DatabaseName([SlashPositions(end)+1:length(DatabaseName)]));
    Time=h5readatt([folder,file],[DatabaseName],'time');
    try
        Elec_Id=h5read([folder,file],[DatabaseName, '/particles/', species, '/id']);
    end
    Elec_Weight=h5read([folder,file],[DatabaseName,'/particles/', species, '/weighting']);

    Elec_Px=h5read([folder,file],[DatabaseName, '/particles/', species, '/momentum/x']);
    Elec_Pz=h5read([folder,file],[DatabaseName, '/particles/', species, '/momentum/z']);


    Elec_X=h5read([folder,file],[DatabaseName, '/particles/', species,  '/position/x']);
    Elec_Z=h5read([folder,file],[DatabaseName, '/particles/', species, '/position/z']);
    
    try
        Elec_Py=h5read([folder,file],[DatabaseName, '/particles/', species, '/momentum/y']);
        Elec_Y=h5read([folder,file],[DatabaseName, '/particles/', species,  '/position/y']);
        
        
    end
     
    Elec.x=Elec_X;
    try       
        Elec.y=Elec_Y;  
    end    
    Elec.z=Elec_Z;
    Elec.weight=Elec_Weight;
    Elec.px=Elec_Px;
    try
        Elec.py=Elec_Py;
    end
    Elec.pz=Elec_Pz;
    
    try
        Elec.id=double(Elec_Id);
    end
    
    
    Elec.Time=Time;    
    
end